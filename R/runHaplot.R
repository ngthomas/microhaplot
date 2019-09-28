#' Run shiny microhaplot
#'
#' Run shiny microhaplot app
#' @param path Path to shiny microhaplot app. Optional. If not specified, the path is default to local app path.
#' @return Runs shiny microhaplot application via \code{shiny::runApp} which typically doesn't return; interrupt R to stop the application (usually by pressing Ctrl+C or Esc).
#' @export
#' @examples
#' if(interactive()){
#' runShinyHaplot()
#' }
runShinyHaplot <- function(path = system.file("shiny", "microhaplot", package = "microhaplot")) {
  if (path == "" || !file.exists(path)) {
    #stop("Could not find Shiny directory. Try re-installing `mypackage`.", call. = FALSE)
    stop("Could not find Shiny directory", call. = FALSE)
  }
  shiny::runApp(path, display.mode = "normal")
}

#' Transfer a copy of microhaplot app.
#'
#' Moves shiny microhaplot app to a different directory
#' @param path string. directory path. Required
#' @return a logical value of whether \code{file.copy} is successfully transferred Shiny app to its new directory
#' @export
#' @examples
#'
#' mvShinyHaplot(tempdir())
#'
mvShinyHaplot <- function(path) {
  app.dir <- system.file("shiny", "microhaplot", package = "microhaplot")
  if (app.dir == "") {
    stop("Could not find shiny directory. Try re-installing `mypackage`.", call. = FALSE)
  }

  if (!file.exists(paste0(path))) dir.create(path)

  file.copy(app.dir, path, recursive = TRUE)
}




#' Extracts haplotype from alignment reads.
#'
#' The function \code{microhaplot} extracts haplotype from sequence alignment files through perl script \code{hapture} and returns a summary table of the read depth and read quality associate with haplotype.
#'
#' @param run.label character vector. Run label to be used to display in haPLOType. Required
#' @param sam.path string. Directory path folder containing all sequence alignment files (SAM). Required
#' @param label.path string. Label file path. This customized label file is a tab-separate file that contains entries of SAM file name, individual ID, and group label. Required
#' @param vcf.path string. VCF file path. Required
#' @param out.path string. Optional. If not specified, the intermediate files are created under \code{TEMPDIR}, with the assumption that directory is granted for written permission.
#' @param add.filter boolean. Optional. If true, this removes any haplotype with unknown and deletion alignment characters i.e. "*" and "_", removes any locus with large number of haplotypes ( # > 40) , and remove any locus with fewer than half of the total individuals.
#' @param app.path string. Path to shiny haPLOType app. Optional. If not specified, the path is default to \code{TEMPDIR}.
#' @param n.jobs positive integer. Number of SAM files to be parallel processed. Optional. This multithread is only available for non Window OS. Recommend two times the number of processors/core.
#' @return This function returns a dataframe of 9 columns i.e group, id, locus, haplotype, depth, sum of Phred score, max of Phred score, allele balance and haplotype rank from highest to lowest read depth. This dataframe will also be saved in \code{out.path}.
#' @export
#' @examples
#'
#' run.label <- "sebastes"
#'
#' sam.path <- tempdir()
#' untar(system.file("extdata",
#'                   "sebastes_sam.tar.gz",
#'                   package="microhaplot"),
#'       exdir = sam.path)
#'
#'
#' label.path <- file.path(sam.path, "label.txt")
#' vcf.path <- file.path(sam.path, "sebastes.vcf")
#'
#' mvShinyHaplot(tempdir())
#' app.path <- file.path(tempdir(), "microhaplot")
#'
#' # retrieve system Perl version number
#' perl.version <- as.numeric(system('perl -e "print $];"', intern=TRUE))
#'
#' if (perl.version >= 5.014) {
#' haplo.read.tbl <- prepHaplotFiles(run.label = run.label,
#'                             sam.path = sam.path,
#'                             out.path = tempdir(),
#'                             label.path = label.path,
#'                             vcf.path = vcf.path,
#'                             app.path = app.path)
#' }else {
#' message("Perl version is outdated. Must >= 5.014.")}
#'
prepHaplotFiles <- function(run.label, sam.path, label.path, vcf.path,
  out.path=tempdir(),
  add.filter=FALSE,
  app.path=tempdir(),
  n.jobs = 1){

  run.label <- gsub(" +","_",run.label)
  haptureDir <- system.file("perl", "hapture", package = "microhaplot")

  # Need to check whether all path and files exist
  if (!file.exists(sam.path)) stop("the path for 'sam.path' - ", sam.path, " does not exist")
  if (!file.exists(label.path)) stop("the path for 'label.path' - ", label.path, " does not exist")
  if (!file.exists(vcf.path)) stop("the path for 'vcf.path' - ", vcf.path, " does not exist")
  if (!file.exists(out.path)) stop("the path for 'out.path' - ", out.path, " does not exist")
  if (!file.exists(app.path)) stop("the path for 'app.path' - ", app.path, " does not exist; try to run mvShinyHaplot()")

  # check for numeric state

  if(!is.numeric(n.jobs)) stop("the n.jobs sepcificed is expected to be numeric")

  n.jobs <- round(n.jobs)
  if(n.jobs<=0) stop("the n.jobs is expected to be positive integer")

  # check whether perl is installed
  tryCatch({system("perl -v", intern=TRUE); message("Perl is found in system")},
           error= function(d) {
             if(.Platform$OS.type == "windows") {
               message("Perl is not found in system. Recommend installation from strawberryperl. Be sure to install version >=5.014")
             }else {
               message("Perl is not found in system. Recommend Installation from perl.org. Be sure to install version >=5.014")
             }
           })

  # ensure that the perl's version is at least 5.014
  perl.version <- system('perl -e "print $];"', intern=TRUE) %>% as.numeric
  if(perl.version < 5.014) stop ("The version Perl found in your current system is old-dated/incompatible. Microhaplot requires Perl v. >=5.014.")

  # the perl script hapture should display any warning if the label field contains any missing or invalid elements

  runHap.name <- ifelse(.Platform$OS.type == "windows", "runHapture.bat", "runHapture.sh")

  if (file.exists(file.path(out.path, runHap.name))) file.remove(file.path(out.path, runHap.name))

  if (file.exists(file.path(out.path,"intermed"))) file.remove(
    list.files(file.path(out.path, "intermed"),
               #pattern = paste0(run.label, "_*.summary"),
               full.names=TRUE))

  if (!file.exists(file.path(out.path,"intermed"))) dir.create(file.path(out.path,"intermed"))

  summary.path <- file.path(out.path, "intermed", "all.summary")

  if(file.exists(summary.path)) file.remove(summary.path)

  # catch any problem in label file
  read.label <- tryCatch(read.table(label.path,sep = "\t",stringsAsFactors = FALSE), error = function(c) {
    c$message <- paste0(c$message, " (in ", label.path , ")")
    stop(c)
  })
  if (dim(read.label)[2] < 3) stop(label.path, "contains less than 3 columns.")

  garb <- sapply(1:nrow(read.label), function(i) {

    line <- read.label[i,] %>% unlist
    if (!file.exists(file.path(sam.path,line[1]))) stop("the SAM file, ", file.path(sam.path,line[1]), ", does not exist")

    if(.Platform$OS.type == "windows") {
      run.perl.script <- paste0("perl ", haptureDir,
                                " -v ", vcf.path, " ",
                                " -s ", sam.path, "\\", line[1],
                                " -i ", line[2],
                                " -g ", line[3], " > ",
                                out.path, "\\intermed\\", run.label, "_", line[2],"_",i,".summary")
    } else {
      wait.ln <- ifelse(i %% n.jobs == 0," wait;"," ")
      run.perl.script <- paste0("perl ", haptureDir,
                                " -v ", vcf.path, " ",
                                " -s ", sam.path, "/", line[1],
                                " -i ", line[2],
                                " -g ", line[3], " > ",
                                out.path, "/intermed/", run.label, "_", line[2],"_",i,".summary &",
                                wait.ln)
    }

    write(run.perl.script,
      file = file.path(out.path, runHap.name),
      append = TRUE)
  })

  message("...running Hapture.pl to extract haplotype information (takes a while)...")

  if(.Platform$OS.type != "windows") {

    write(paste0("wait;"),
    file = file.path(out.path, runHap.name),
    append = TRUE)

    concat.cmd <- paste0("cat ",out.path, "/intermed/", run.label, "_", "*.summary"," > ",summary.path)

    # just in case if the user has loads of sam files and running out of buffer
    if (nrow(read.label) > 100) {
      concat.cmd <- paste0("find ",out.path, "/intermed -name ", run.label, "_", "*.summary",
                            "| while read F; do cat ${F} >>",summary.path, ";done")
    }

    write(concat.cmd,
          file = file.path(out.path, runHap.name),
          append = TRUE)

    system(paste0("bash ",out.path,"/runHapture.sh"))
  } else {
    concat.cmd <- paste0("type ",out.path, "\\intermed\\", run.label, "_", "*.summary"," > ",summary.path)

    write(concat.cmd,
          file = file.path(out.path, runHap.name),
          append = TRUE)


    system(file.path(out.path, runHap.name))
  }

  haplo.sum <- read.table(summary.path, stringsAsFactors = FALSE, sep = "\t") %>% dplyr::tbl_df()

  colnames(haplo.sum) <- c("group", "id", "locus", "haplo", "depth", "sum.Phred.C", "max.Phred.C")

  num.id <- length(unique(haplo.sum$id))

  message(paste0("\n...Prepping rds file : ",out.path, "/",run.label,".rds\n"))

  if (add.filter) {

    haplo.cleanup <- haplo.sum %>%
      dplyr::filter(!grepl("[N]", haplo)) %>%
      dplyr::group_by(locus, id) %>%
      dplyr::mutate(n.haplo.per.indiv = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(locus) %>%
      dplyr::mutate(n.indiv.per.locus = length(unique(id)), max.uniq.hapl = max(n.haplo.per.indiv)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(n.indiv.per.locus > num.id/2, max.uniq.hapl < 40)  %>%
      dplyr::select(group, id, locus, haplo, depth, sum.Phred.C, max.Phred.C)
    } else {
      haplo.cleanup <- haplo.sum %>% dplyr::select(group, id, locus, haplo, depth, sum.Phred.C, max.Phred.C)
    }

  haplo.add.balance <- haplo.cleanup %>%
    dplyr::arrange(dplyr::desc(depth)) %>%
    dplyr::group_by(locus,id) %>%
    dplyr::mutate(allele.balance = depth/depth[1], rank = dplyr::row_number() ) %>%
    dplyr::ungroup()

  vcf.pos.tbl <- read.table(vcf.path) %>%
    .[,1:2] %>% # grabbing locus name, and pos
    dplyr::group_by(V1) %>%
    dplyr::summarise(pos = paste0(V2, collapse = ","))

  colnames(vcf.pos.tbl) <- c("locus","pos")

  saveRDS(haplo.add.balance, paste0(out.path, "/",run.label,".rds"))
  saveRDS(vcf.pos.tbl, paste0(out.path, "/",run.label,"_posinfo.rds"))

  message(paste0("RDS file: copied into shiny directory: ",app.path, "/",run.label,"*.rds ", "\n\nRun the following to open shiny app: runShinyHaplot(\"",app.path,"\")"))

  if(file.exists(file.path(app.path, paste0(run.label,"*.rds")))) message("Overwritting previous version -  ", run.label, ".rds in ", app.path)


  if (out.path != app.path) system(paste0("cp ", out.path, "/",run.label,"*.rds ", app.path, "/.") )

  return(haplo.add.balance)
}
