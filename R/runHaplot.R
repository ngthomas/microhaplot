#' Run Haplotype Plot.
#'
#' Run Haplotype shiny
#' @path Path to shiny haPLOType app. Optional. If not specified, the path is default to local app path.
#' @export
#' @examples
#' runHaplotype()
runHaplotype <- function(path=system.file("shiny", "microhaplot", package = "microhaplot")) {
  if (path == "" || !file.exists(path)) {
    #stop("Could not find Shiny directory. Try re-installing `mypackage`.", call. = FALSE)
    stop("Could not find Shiny directory", call. = FALSE)
  }
  shiny::runApp(path, display.mode = "normal")
}

#' Redirect haPLOType app.
#'
#' makes a copy of shiny haPLOType to a different directory
#' @param path string. directory path. Required
#' @export
#' @examples
#' runHaplotype()
mvHaplotype <- function(path) {
  app.dir <- system.file("shiny", "microhaplot", package = "microhaplot")
  if (app.dir == "") {
    stop("Could not find shiny directory. Try re-installing `mypackage`.", call. = FALSE)
  }

  if(!file.exists(paste0(path))) dir.create(path)

  file.copy(paste0(app.dir, "/"),
            paste0(path, "/"),
            recursive =T)
}




#' Extract haplotype from alignment reads.
#'
#' The function \code{microhaplot} extracts haplotype from sequence alignment files through perl script \code{hapture} and returns a summary table of the read depth and read quality associate with haplotype.
#'
#' @param run.label character vector. Run label to be used to display in haPLOType. Required
#' @param sam.path string. Directory path folder containing all sequence alignment files (SAM). Required
#' @param label.path string. Label file path. This customized label file is a tab-separate file that contains entries of SAM file name, individual ID, and group label. Required
#' @param vcf.path string. VCF file path. Required
#' @param out.path string. Optional. If not specified, the intermediate files are created under \code{sam.path}, with the assumption that directory is granted for written permission.
#' @param add.filter boolean. Optional. If true, this removes any haplotype with unknown and deletion alignment characters i.e. "*" and "_", removes any locus with large number of haplotypes ( # > 40) , and remove any locus with fewer than half of the total individuals.
#' @param app.path string. Path to shiny haPLOType app. Optional. If not specified, the path is default to local app path.
#' @export
#' @examples
#' run.label<-"example 1"
#' sam.path<-"data/satro_sample"
#' label.path <- "data/satro_sample/sample_label.txt"
#' vcf.path <- "data/satro_sample/sebastes.vcf"
#' # runHaplot(run.label, sam.path, label.path, vcf.path)
runHaplot <- function(run.label, sam.path, label.path, vcf.path,
  out.path=sam.path,
  add.filter=FALSE,
  app.path=system.file("shiny", "microhaplot", package = "microhaplot")){

  run.label <- gsub(" +","_",run.label)
  haptureDir <- system.file("perl", "hapture", package = "microhaplot")

  # Need to check whether all path and files exist
  if(!file.exists(paste0(sam.path))) stop("the path for 'sam.path' - ", sam.path, " does not exist")
  if(!file.exists(paste0(label.path))) stop("the path for 'label.path' - ", label.path, " does not exist")
  if(!file.exists(paste0(vcf.path))) stop("the path for 'vcf.path' - ", vcf.path, " does not exist")
  if(!file.exists(paste0(out.path))) stop("the path for 'out.path' - ", out.path, " does not exist")

  # the perl script hapture should display any warning if the label field contains any missing or invalid elements


  system(paste0("rm -f ", out.path, "/runHapture.sh"))
  if(file.exists(paste0(out.path,"/intermed"))) system(paste0("rm -f ", out.path, "/intermed/", run.label, "_", "*.summary;"))
  if(!file.exists(paste0(out.path,"/intermed"))) dir.create(paste0(out.path,"/intermed"))

  # catch any problem in label file
  read.label <- tryCatch(read.table(label.path,sep="\t",stringsAsFactors = F), error = function(c) {
    c$message <- paste0(c$message, " (in ", label.path , ")")
    stop(c)
  })
  if (dim(read.label)[2]<3) stop(label.path, "contains less than 3 columns.")



  garb <- sapply(1:nrow(read.label), function(i) {

    line <- read.label[i,] %>% unlist
    if(!file.exists(paste0(sam.path,"/",line[1]))) stop("the SAM file, ",
                                                        sam.path,"/",line[1], ", does not exist")
    run.perl.script <- paste0("perl ", haptureDir,
      " -v ", vcf.path, " ",
      " -s ", sam.path, "/", line[1],
      " -i ", line[2],
      " -g ", line[3], " > ",
      out.path, "/intermed/", run.label, "_", line[2],"_",i,".summary &");

    write(run.perl.script,
      file=paste0(out.path, "/runHapture.sh"),
      append=T)
    if(i%%10==0) write("wait;",
                       file=paste0(out.path, "/runHapture.sh"),
                       append=T)
  })


  write(paste0("wait; exit 0;"),
    file=paste0(out.path, "/runHapture.sh"),
    append=T)

  cat("...running Hapture.pl to extract haplotype information (takes a while)...")
  system(paste0("bash ",out.path,"/runHapture.sh"))

  summary.tbl<-paste0(out.path,"/intermed/all.summary")

  concat.file <- paste0("cat ",out.path, "/intermed/", run.label, "_", "*.summary",">",summary.tbl)
  system(concat.file)

  haplo.sum <- read.table(summary.tbl, stringsAsFactors = FALSE, sep="\t") %>% dplyr::tbl_df()

  colnames(haplo.sum) <- c("group", "id", "locus", "haplo", "depth", "sum.Phred.C", "max.Phred.C")

  num.id <- length(unique(haplo.sum$id))

  cat(paste0("\n...Prepping feather file : ",out.path, "/",run.label,".feather\n"))

  if (add.filter) {

    haplo.cleanup <- haplo.sum %>%
      dplyr::filter(!grepl("[N]", haplo)) %>%
      dplyr::group_by(locus, id) %>%
      dplyr::mutate(n.haplo.per.indiv=n()) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(locus) %>%
      dplyr::mutate(n.indiv.per.locus = length(unique(id)), max.uniq.hapl=max(n.haplo.per.indiv)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(n.indiv.per.locus > num.id/2, max.uniq.hapl < 40)  %>%
      dplyr::select(group, id, locus, haplo, depth, sum.Phred.C, max.Phred.C) }
  else {
    haplo.cleanup <- haplo.sum %>% dplyr::select(group, id, locus, haplo, depth, sum.Phred.C, max.Phred.C)}

  haplo.add.balance <- haplo.cleanup %>%
    dplyr::arrange(desc(depth)) %>%
    dplyr::group_by(locus,id) %>%
    dplyr::mutate(allele.balance = depth/depth[1], rank=row_number() ) %>%
    dplyr::ungroup()

  vcf.pos.tbl <- read.table(vcf.path) %>%
    .[,1:2] %>% # grabbing locus name, and pos
    dplyr::group_by(V1) %>%
    dplyr::summarise(pos=paste0(V2, collapse=","))

  colnames(vcf.pos.tbl) <- c("locus","pos")

  feather::write_feather(haplo.add.balance, paste0(out.path, "/",run.label,".feather"))
  feather::write_feather(vcf.pos.tbl, paste0(out.path, "/",run.label,"_posinfo.feather"))

  #saveRDS(haplo.add.balance, paste0(out.path, "/",run.label,".rds"))
  #saveRDS(vcf.pos.tbl, paste0(out.path, "/",run.label,"_posinfo.rds"))

  cat(paste0("\n\nFeather file: copied into shiny directory: ",app.path, "/",run.label,"*.feather ",
             "\nRun runHaplotype() to open shiny app.\n\n"))

  system(paste0("cp ", out.path, "/",run.label,"*.feather ", app.path, "/.") )

  return(haplo.add.balance)
}
