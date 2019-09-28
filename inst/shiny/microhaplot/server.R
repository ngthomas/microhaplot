library("shiny")
library("shinyBS")
library("ggplot2")
library("plyr")
library("dplyr")
library("tidyr")
library("DT")
library("grid")
library("scales")
library("microhaplot")
library("ggiraph")

shinyServer(function(input, output, session) {


  output$about <- renderUI({
    HTML(
      paste("<strong>Microhaplot</strong> offers a streamline visual environment to assess quality of microhaplotype from short read alignment files.
          You can find most of the interactive features such as defining critera or adding locus comments at the top panel
while the bottom panel hosts a wide selection of tables and graphical summaries.",
            "",
            "<i>Summary Info</i>: microhaplotype summaries that are either grouped by group label, individual, or locus. These plots are great
          for gathering the big picture",
            "",
            "<i>Filter Analysis</i>: this section is useful to fine-tune your criteria at a single locus level",
            "",
            "<i>Inferential Analysis</i>: still in the works",
            "",
            "<i>Output</i>: You can view or download tables of raw or finalized microhaplotypes (in csv format)",
            "",
            "<b>Contact & Citation</b>",
            "If you have questions or suggest, feel free to contact me at <b>tngthomasng@gmail.com</b>",
            "",
            "<b>citation</b>: Ng, Thomas C. & Anderson Eric C. (2017, June 30). ngthomas/microhaplot: microhaplotype viewer. Zenodo. http://doi.org/10.5281/zenodo.821679",
            "<br/>",
            '<a href="https://doi.org/10.5281/zenodo.821679"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.821679.svg" alt="DOI"></a>',
            sep="<br/>"))
  })

  addTooltip(session, "downloadData",
             "Note: The download table is sensitive to the 'field selection' input")

  addPopover(session, "allReadDepth", "Histogram of total read depth", content=paste0("<b>double-click</b>: to set new minimum total read depth"), placement="top", trigger="hover")
  addPopover(session, "allAllelicRatio", "Histogram of allelic ratio", content=paste0("<b>double-click</b>: to set new minimum allelic ratio"), placement="top", trigger="hover")

  addPopover(session, "RDnARplot", "",
             content=paste0("<p><b>top</b>: # of unique qualified haplotype<br>
                            <b>bot</b>: # of indiv made the cutoff<br>
                            Report mean when 'all' loci are selected"),
             placement="top",
             trigger="hover")

  addPopover(session, "filterOpts","Options",
             content=paste0("<p>By default, all microhap. must pass the selected filter criteria on the left. ",
                            "The criteria can be overrided or further imposed by choosing these options:",
                            "</p>",
                            "<p></p>",
                            "<p><b>overrides only as min baseline</b>: this option set the current selected critera as the lower bound for every single locus</p>",
                            "<p></p>",
                            "<p><b>overrides all params</b>: this option applies the current selected critera to every single locus.</p> ",
                            "<p><b>relies on local locus param.</b>: this option disregards the current selection of min. read depth and allelic ratio. ",
                            "Instead, the filtering process uses parameters previously defined at a single-locus level</p>"),
             placement="bottom",
             trigger="hover")

  dirFiles <- list.files()

  rds.file <- grep(".rds", dirFiles)
  pos.rds.file <- grep("_posinfo.rds", dirFiles)
  annotate.rds.file <- grep("_annotate.rds", dirFiles)
  rds.file <- setdiff(rds.file, c(pos.rds.file,annotate.rds.file))

  file.name.ls <- rds.file

  haplo.sum <- NULL

  if (length(file.name.ls) > 0) {
    select.file.tem <- dirFiles[file.name.ls[1]]
    updateSelectInput(session,
                      "selectDB",
                      selected = select.file.tem,
                      choices = dirFiles[file.name.ls])
  }

  update.Haplo.file <- reactive({
    if (input$selectDB == "" ||
        is.null(input$selectDB) || !file.exists(input$selectDB))
      return()
    #message( "select DB_", input$selectDB, "_----\n")

    table.out <-readRDS(input$selectDB)  %>%
      ungroup() %>%
      mutate(id = as.character(id),
             locus =as.character(locus),
             group=as.character(group))

    if ("sum.Phred.C" %in% colnames(table.out)) {
      table.out <- table.out %>% select(-sum.Phred.C, -max.Phred.C)}

    table.out
  })

  extract.pos.file <- reactive({

    pos.file <- strsplit(input$selectDB, split=".rds") %>% unlist %>% paste0(.,"_posinfo.rds")
    if (!file.exists(pos.file)) return ()
    readRDS(pos.file)

  })

  # the format structure for the annotate file: contains locus ID (also contains ALL), min read depth, min allelic ratio,
  # keep status (1=keep), annotation

  extract.annotate.file <- reactive({

    annotate.file <- strsplit(input$selectDB, split=".rds") %>% unlist %>% paste0(.,"_annotate.rds")

    if (!file.exists(annotate.file)) {
      annotateTab$tbl <- data.frame(locus = as.character(panelParam$locus.label),
                                    ave.entropy = c("NA",rep(0, panelParam$n.locus)),
                                    n.alleles = rep(2, 0, panelParam$n.locus+1),
                                    min.rd = rep(0, panelParam$n.locus+1),
                                    min.ar = rep(0.5, panelParam$n.locus+1),
                                    max.ar.hm = rep(0.5, panelParam$n.locus+1),
                                    min.ar.hz = rep(0.5, panelParam$n.locus+1),
                                    status = c("NA",rep("Accept", panelParam$n.locus)),
                                    comment = rep("", panelParam$n.locus+1),
                                    stringsAsFactors = FALSE)

      return()
    }

    annotateTab$tbl <- readRDS(annotate.file)


    if((!"max.ar.hm" %in% colnames(annotateTab$tbl)) || (!"n.alleles" %in% colnames(annotateTab$tbl))) {
      #if("n.alleles" %in% colnames(annotateTab$tbl))
      #  annotateTab$tbl <- annotateTab$tbl %>% select(-n.alleles)

      annotateTab$tbl <- bind_cols(annotateTab$tbl,
                                   n.alleles = rep(2, 0, panelParam$n.locus+1),
                                   max.ar.hm = rep(0.5, panelParam$n.locus+1),
                                   min.ar.hz = rep(0.5, panelParam$n.locus+1)
      )

      saveRDS(annotateTab$tbl, annotate.file)

      annotateTab$tbl <- annotateTab$tbl %>%
        ungroup() %>%
        mutate(locus = as.character(locus),
               min.rd = as.numeric(min.rd),
               min.ar = as.numeric(min.ar),
               n.alleles = as.integer(n.alleles),
               max.ar.hm = as.numeric(max.ar.hm),
               min.ar.hz = as.numeric(min.ar.hz))


    }

    annotateTab$tbl <- annotateTab$tbl %>%
      group_by(locus) %>%
      mutate(max.ar.hm = ifelse(max.ar.hm > min.ar, min.ar, max.ar.hm),
             min.ar.hz = ifelse(min.ar.hz < min.ar, min.ar, min.ar.hz)) %>%
      ungroup()


  })

  observeEvent(input$max.ar.hm , {
    filterParam$max.ar.hm <- input$max.ar.hm
  })
  observeEvent(input$min.ar.hz, {
    filterParam$min.ar.hz <- input$min.ar.hz
  })

  update.annotate.field <- reactive({
    match.indx <- which(annotateTab$tbl$locus == input$selectLocus)

    updateTextInput(session,
                    "locusComment",
                    value=annotateTab$tbl$comment[match.indx])
    updateSelectInput(session,
                      "locusAccept",
                      selected=annotateTab$tbl$status[match.indx]
    )

    updateSliderTextInput(session, "n.alleles", selected=annotateTab$tbl$n.alleles[match.indx])
    updateSliderTextInput(session, "coverageMin", selected=annotateTab$tbl$min.rd[match.indx])
    updateSliderTextInput(session, "minAlleleRatio", selected=annotateTab$tbl$min.ar[match.indx])

    updateSliderInput(session, "max.ar.hm", max = annotateTab$tbl$min.ar[match.indx],
                      value = annotateTab$tbl$max.ar.hm[match.indx])

    updateSliderInput(session, "min.ar.hz", min = annotateTab$tbl$min.ar[match.indx],
                      value = annotateTab$tbl$min.ar.hz[match.indx])

    filterParam$minRD <- annotateTab$tbl$min.rd[match.indx]#input$coverageMin

    filterParam$minAR <- annotateTab$tbl$min.ar[match.indx]#input$minAlleleRatio
    filterParam$n.alleles <- annotateTab$tbl$n.alleles[match.indx]
    filterParam$max.ar.hm <- annotateTab$tbl$max.ar.hm[match.indx]
    filterParam$min.ar.hz <- annotateTab$tbl$min.ar.hz[match.indx]



    updateSliderTextInput(session, "max_read_depth", selected=1280)

    # addPopover(session, "filterSave","Saved Values (current)",
    #            content=paste0("<b>n.alleles</b>: ", annotateTab$tbl$n.alleles[match.indx], "<br>",
    #                           "<b>min.rd</b>: ", annotateTab$tbl$min.rd[match.indx], "<br>",
    #                           "<b>min.ar</b>: ", annotateTab$tbl$min.ar[match.indx], "<br><br>",
    #                           "<b>ar.hm</b>: ", annotateTab$tbl$max.ar.hm[match.indx], "<br>",
    #                           "<b>ar.hz</b>: ", annotateTab$tbl$min.ar.hz[match.indx]),
    #            trigger = "hover",  placement = "bottom")


    return()
  })

  Update.ave.entropy <- reactive ({
    haplo.filter <- Filter.haplo.by.RDnAR()
    if (is.null(haplo.filter))
      return ()

    entropy.tbl <- haplo.filter %>%
      filter(locus !="ALL", rank <= 2) %>%
      group_by(locus, group, haplo) %>% summarise(n=n()) %>%
      ungroup() %>%
      group_by(locus, group) %>%
      mutate(f = n/n()) %>%
      ungroup() %>% group_by(locus, haplo) %>%
      mutate(shan.entropy = f/sum(f)*log(sum(f)/f, base=2)) %>%
      ungroup() %>% group_by(locus) %>%
      summarise(ave.entropy.tem = round(mean(shan.entropy) ,3))

    annotateTab$tbl <- left_join(annotateTab$tbl, entropy.tbl, by="locus") %>%
      mutate(ave.entropy = ifelse(is.na(ave.entropy.tem), ave.entropy,
                                  ave.entropy.tem)) %>%
      select(-ave.entropy.tem)

  })

  # needs better label
  ranges <- reactiveValues(y = NULL, x = NULL, aip = NULL, alp = NULL)
  rangesH <- reactiveValues(y = NULL)

  locusPg <- reactiveValues(l = NULL, width = NULL, maxPg = 1)
  indivPg <- reactiveValues(i = NULL, width = NULL, maxPg = 1)
  groupPg <- reactiveValues(g = NULL, width = 1)
  hapPg <- reactiveValues(width = NULL, HW.exp= NULL, HW.obs= NULL, HW.tbl=NULL,
                          indiv.hap.grp=NULL)

  filterCriteriaPg <- reactiveValues(RD.detail.on =0,
                                     AR.detail.on =0)

  annotateTab <- reactiveValues(tbl=NULL)
  srhapPg <- reactiveValues(
    makePlot = FALSE,
    num.iter = NULL,
    frac.burn = NULL,
    random.seed = NULL,
    prior.model = NULL,
    locus.select = NULL,
    min.read.depth = NULL,
    data.table = NULL
  )

  filterParam <- reactiveValues(minRD = 1,
                                minAR = 0.5,
                                hover.minAR = 0, hover.minRD =0,
                                minRDhap = 1,
                                n.alleles=2,
                                max.ar.hm = 0.5,
                                min.ar.hz = 0.5,
                                opts = NULL)
  panelParam <- reactiveValues(
    n.locus = NULL,
    n.indiv = NULL,
    tot.indiv = NULL,
    locus.label.tbl = NULL,
    locus.label = NULL,
    locus.label.bare = NULL,
    indiv.label.tbl = NULL,
    indiv.label = NULL,
    indiv.label.bare = NULL,
    is.reject = NULL,
    n.group = 0,
    group.label = NULL,
    group.label.tbl = NULL,
    group.label.bare = NULL
  )


  observeEvent(input$selectDB, {

    message("select DB_", input$selectDB, "_", "----\n")
    if (input$selectDB == "" || is.null(input$selectDB))
      return()

    haplo.sum <- update.Haplo.file()

    panelParam$tot.indiv <- length(unique(haplo.sum$id))

    if (panelParam$tot.indiv ==0) return()
    if (!"group" %in% colnames(haplo.sum))
      haplo.sum <-
      cbind.data.frame("group" = "unlabel",
                       haplo.sum,
                       stringsAsFactors = FALSE) %>% tbl_df

    panelParam$n.locus <- length(unique(haplo.sum$locus))
    panelParam$n.indiv <- length(unique(haplo.sum$id))

    locus.sorted <- sort(unique(haplo.sum$locus))
    panelParam$locus.label.tbl <-
      data.frame(locus = locus.sorted, stringsAsFactors = FALSE) %>% tbl_df()
    panelParam$locus.label <- c("ALL", locus.sorted)
    panelParam$locus.label.bare <- locus.sorted

    indiv.sorted <- sort(unique(haplo.sum$id))
    panelParam$indiv.label.tbl <-
      data.frame(id = indiv.sorted, stringsAsFactors = FALSE) %>% tbl_df()
    panelParam$indiv.label <- c("ALL", indiv.sorted)
    panelParam$indiv.label.bare <- indiv.sorted

    group.sorted <- sort(unique(haplo.sum$group))
    panelParam$group.label.tbl <-
      data.frame(id = group.sorted, stringsAsFactors = FALSE) %>% tbl_df()
    panelParam$group.label <- c("ALL", group.sorted)
    panelParam$group.label.bare <- group.sorted
    panelParam$n.group <- length(group.sorted)

    updateSelectInput(session,
                      inputId="selectLocus",
                      selected = "ALL",
                      choices = panelParam$locus.label)


    updateSelectInput(session,
                      "selectIndiv",
                      selected = "ALL",
                      choices = panelParam$indiv.label)

    updateSelectInput(session,
                      "selectGroup",
                      selected = "ALL",
                      choices = panelParam$group.label)

    end.indx <- min(15, panelParam$n.locus)
    locusPg$l <- panelParam$locus.label.bare[1:end.indx]
    rangesH$y <- c(0, length(locusPg$l) + 1)

    end.indx <- min(15, panelParam$n.indiv)
    indivPg$i <- panelParam$indiv.label.bare[1:end.indx]
    ranges$y <- c(0, length(indivPg$i) + 1)

    session$sendCustomMessage(type = 'biplot_set', message = character(0))

    extract.annotate.file()
    update.annotate.field()

    Min.filter.haplo()

    #Filter.haplo.by.RDnAR()
    #Filter.haplo.sum()

    srhapPg$makePlot <- FALSE
    extract.pos.file()
    Filter.haplo.sum()

    message("Loading complete:" , input$selectDB, "_", "----\n" )
  }, priority = -3)


  # reacting to the locus & Indiv's previous and next button
  observeEvent(input$locusBack, {
    indx <- isolate(which(panelParam$locus.label == input$selectLocus))
    label <-
      ifelse(indx > 1, panelParam$locus.label[indx - 1], panelParam$locus.label[indx])
    updateSelectInput(session, "selectLocus", selected = label)
  })
  observeEvent(input$locusFor, {
    indx <- isolate(which(panelParam$locus.label == input$selectLocus))
    label <-
      ifelse(
        indx < length(panelParam$locus.label),
        panelParam$locus.label[indx + 1],
        panelParam$locus.label[indx]
      )
    updateSelectInput(session, "selectLocus", selected = label)
  })
  observeEvent(input$preL, {
    indx <- isolate(which(panelParam$locus.label == input$selectLocus))
    label <-
      ifelse(indx > 1, panelParam$locus.label[indx - 1], panelParam$locus.label[indx])
    updateSelectInput(session, "selectLocus", selected = label)
  })
  observeEvent(input$nextL, {
    indx <- isolate(which(panelParam$locus.label == input$selectLocus))
    label <-
      ifelse(
        indx < length(panelParam$locus.label),
        panelParam$locus.label[indx + 1],
        panelParam$locus.label[indx]
      )
    updateSelectInput(session, "selectLocus", selected = label)
  })
  observeEvent(input$indivBack, {
    indx <- isolate(which(panelParam$indiv.label == input$selectIndiv))
    label <-
      ifelse(indx > 1, panelParam$indiv.label[indx - 1], panelParam$indiv.label[indx])
    updateSelectInput(session, "selectIndiv", selected = label)
  })
  observeEvent(input$indivFor, {
    indx <- isolate(which(panelParam$indiv.label == input$selectIndiv))
    label <-
      ifelse(
        indx < length(panelParam$indiv.label),
        panelParam$indiv.label[indx + 1],
        panelParam$indiv.label[indx]
      )
    updateSelectInput(session, "selectIndiv", selected = label)
  })

  observeEvent(input$filterOpts, {
    if (1 %in% input$filterOpts || 2 %in% input$filterOpts) {

      sendSweetAlert(
        session = session,
        title = "Notification",
        text = "Under this option, parameters of loci might permeantly get overridden when 'apply' or 'save'!",
        type = "info"
      )
    }
  })

  # reacting to the filter update button
  observeEvent( input$updateFilter, {
    #input$updateFilter, {
    if(is.na(input$coverageMin) || input$coverageMin <0) {
      createAlert(session, "alert", "filterAlert", title = "Invalid input",
                  content = "Minimum read coverage must be a positive integer", append = FALSE)
      return()
    }

    closeAlert(session, "filterAlert")

    filterParam$minRD <- input$coverageMin
    filterParam$minAR <- input$minAlleleRatio
    filterParam$n.alleles <- input$n.alleles
    filterParam$opts <- input$filterOpts

    #filterParam$max.ar.hm <- ifelse(input$max.ar.hm > filterParam$minAR, filterParam$minAR, input$max.ar.hm)
    #filterParam$min.ar.hz <- ifelse(input$min.ar.hz < filterParam$minAR, filterParam$minAR, input$min.ar.hz)

    Filter.haplo.by.RDnAR()
    ranges$aip <- 0
    ranges$alp <- 0

    updateSliderInput(session, "max.ar.hm", max = filterParam$minAR, value= filterParam$max.ar.hm)
    updateSliderInput(session, "min.ar.hz", min = filterParam$minAR, value = filterParam$min.ar.hz)
  })

  observeEvent(input$updateRdMin,{
    #input$updateFilter, {
    if(is.na(input$rdMin) || input$rdMin <0) {
      return()
    }

    filterParam$minRDhap <- input$rdMin
  })

  observeEvent(input$selectGroup, {
    if (input$selectDB == "" || is.null(input$selectDB))
      return()

    indx <-
      isolate(which(panelParam$group.label.bare == input$selectGroup))
    haplo.sum <- update.Haplo.file()
    if (input$selectGroup != "ALL")
      haplo.sum <- haplo.sum %>% filter(group == input$selectGroup)

    panelParam$n.indiv <- length(unique(haplo.sum$id))

    indiv.sorted <- sort(unique(haplo.sum$id))
    panelParam$indiv.label.tbl <-
      data.frame(id = indiv.sorted, stringsAsFactors = FALSE) %>% tbl_df()

    panelParam$indiv.label <- c("ALL", indiv.sorted)
    panelParam$indiv.label.bare <- indiv.sorted
    updateSelectInput(session,
                      "selectIndiv",
                      selected = "ALL",
                      choices = panelParam$indiv.label)

    end.indx <- min(15, panelParam$n.indiv)
    indivPg$i <- panelParam$indiv.label.bare[1:end.indx]
    ranges$y <- c(0, length(indivPg$i) + 1)

    Min.filter.haplo()

    session$onFlushed(function(){
      session$sendCustomMessage(type = 'biplot_set', message = isolate(input$biplot_selected))
    },
    once=TRUE)

  }, priority = -2)

  observeEvent(input$max_read_depth,{
    session$onFlushed(function(){
      session$sendCustomMessage(type = 'biplot_set', message = isolate(input$biplot_selected))
    },
    once=TRUE)
  })



  observeEvent(input$selectLocus, {

    if (input$selectDB == "" || is.null(input$selectDB))
      return()

    indx <-isolate(which(panelParam$locus.label.bare == input$selectLocus))

    if(is.null(panelParam$n.locus)) return()

    output$locusSelect <- renderText({
      input$selectLocus
    })
    output$locusSelect1 <- renderText({
      input$selectLocus
    })

    output$DP1 <- renderText({
      ifelse(input$selectLocus == "ALL",
             "",
             "Distribution prior to RD & AR filter:")})
    output$DP2 <- renderText({ifelse(input$selectLocus == "ALL",
                                     "",
                                     "RD & AR Distribution post-filter:")})
    output$DP3 <- renderText({
      ifelse(input$selectLocus == "ALL",
             "",
             "within each rank, top/bottom = homoz/het ; Green = pass/call, Red = not pass/no call")})


    update.annotate.field()

    ## for each indiv locus
    if (input$selectLocus != "ALL") {
      output$maxlocusPage <- renderText({
        "1"
      })
      locusPg$maxPg <- 1
      filterParam$opts <- NULL
      updateNumericInput(session, "locusPage", value = 1, max = 1)

      updateCheckboxGroupInput(session,
                               "filterOpts",
                               choices=list())

      locusPg$l <- input$selectLocus
      rangesH$y <- c(0, 2)

      refineDiploidCall()
      # output$locusAcceptStatus <-
      #   renderText({
      #     ifelse(panelParam$is.reject[indx] == 0, "Accept", "Reject")
      #   })

      closeAlert(session,"hapLocusAlert")
      closeAlert(session,"cuthapLocusAlert")
      closeAlert(session, "ARLocusAlert")
    }
    else {
      updateCheckboxGroupInput(session,
                               "filterOpts",
                               choices=list("overrides only as min baseline"=1,
                                            "overrides all params"=2,
                                            "relies on local locus param."=3),
                               selected = c(0))
      filterParam$opts <- input$filterOpts

      # output$locusAcceptStatus <- renderText({
      #   "NA"
      # })
      output$maxlocusPage <-
        renderText({
          paste0(ceiling(
            as.numeric(panelParam$n.locus) /
              as.numeric(input$locusPerDisplay)
          ))
        })
      locusPg$maxPg <- ceiling(
        as.numeric(panelParam$n.locus) /
          as.numeric(input$locusPerDisplay)
      )
      updateNumericInput(
        session,
        "locusPage",
        value = 1,
        max = ceiling(panelParam$n.locus / 15)
      )
      updateSelectInput(session, "locusPerDisplay", selected = 15)
      end.indx <- min(15, panelParam$n.locus)
      locusPg$l <- panelParam$locus.label.bare[1:end.indx]
      rangesH$y <- c(0, end.indx + 1)

      createAlert(session, "hapAlert", "hapLocusAlert", title = "Single locus content",
                  content = "choose a locus to view content", append = FALSE)
      createAlert(session, "ARalert", "ARLocusAlert", title = "Single locus content",
                  content = "choose a locus to view content", append = FALSE)
      createAlert(session, "cutoffhapAlert", "cuthapLocusAlert", title = "choose single locus for more detail",
                  content = "choose a locus to view a more detailed breakdown of those filter criteria by microhaplotype", append = FALSE)
    }

    if(isolate(input$keepSel)){
      session$onFlushed(function(){
        session$sendCustomMessage(type = 'biplot_set', message = isolate(input$biplot_selected))
      },
      once=TRUE)
    } else {
      session$onFlushed(function(){
        session$sendCustomMessage(type = 'biplot_set', message = character(0))
      },
      once=TRUE)
    }
    # eventReactive(input$biplot,{
    #   #session$sendCustomMessage(type = 'biplot_set', message = character(0))
    #   session$onFlushed(function() {
    #     message("hidude")
    #     session$sendCustomMessage(type = 'biplot_set', message = isolate(input$biplot_selected))
    #     #session$sendCustomMessage(type = 'biplot_set', message = isolate(input$biplot_selected))
    #
    #   })

    Filter.haplo.by.RDnAR()
    #srhapPg$makePlot <- FALSE
  })

  observeEvent(input$locusPerDisplay, {
    if (is.null(panelParam$n.locus))
      return()

    if (input$selectLocus != "ALL") {
      output$maxlocusPage <- renderText({
        "1"
      })
      locusPg$maxPg <- 1
      updateNumericInput(session, "locusPage", value = 1, max = 1)
      locusPg$l <- input$selectLocus
      rangesH$y <- c(0, 2)
    }
    else {
      if (input$locusPerDisplay == 100) {
        output$maxlocusPage <- renderText({
          "1"
        })
        locusPg$maxPg <- 1
        updateNumericInput(session,
                           "locusPage",
                           value = 1,
                           max = 1)
        locusPg$l <- panelParam$locus.label.bare
        rangesH$y <- c(0, length(locusPg$l) + 1)
      }
      else{
        output$maxlocusPage <-
          renderText({
            paste0(ceiling(
              as.numeric(panelParam$n.locus) /
                as.numeric(input$locusPerDisplay)
            ))
          })
        locusPg$maxPg <- ceiling(
          as.numeric(panelParam$n.locus) /
            as.numeric(input$locusPerDisplay)
        )
        #message( "haha_", as.numeric(panelParam$n.locus)/as.numeric(input$locusPerDisplay), "_----\n")
        updateNumericInput(session, "locusPage", max = ceiling(
          as.numeric(panelParam$n.locus)
          / as.numeric(input$locusPerDisplay)
        ))

        end.indx <-
          min(as.numeric(panelParam$n.locus) ,
              as.numeric(input$locusPerDisplay))
        locusPg$l <- panelParam$locus.label.bare[1:end.indx]
        rangesH$y <- c(0, length(locusPg$l) + 1)
      }
    }

  })

  observeEvent(input$locusPage, {#updateLocusSizeDisplay, {
    if(is.null(input$locusPage) | is.na(input$locusPage) |!is.numeric(input$locusPage)) return()

    if(input$locusPage <=0) updateNumericInput(session, "locusPage", value = 1)

    if(input$locusPage > locusPg$maxPg) updateNumericInput(session, "locusPage", value = locusPg$maxPg)

    if (input$selectLocus == "ALL") {
      if (input$locusPerDisplay == 100) {
        locusPg$l <- panelParam$locus.label.bare
        rangesH$y <- c(0, length(locusPg$l) + 1)

      }
      else {
        pg <-
          min(as.numeric(input$locusPage),
              ceiling(
                as.numeric(panelParam$n.locus) /
                  as.numeric(input$locusPerDisplay)
              ))
        start.indx <-
          (as.numeric(input$locusPerDisplay) * (pg - 1)) + 1
        end.indx <-
          min(as.numeric(panelParam$n.locus) ,
              as.numeric(input$locusPerDisplay) * pg)
        locusPg$l <-
          panelParam$locus.label.bare[start.indx:end.indx]
        rangesH$y <- c(0, length(locusPg$l) + 1)

      }
    }
  })


  observeEvent(input$selectIndiv, {
    if(is.null(panelParam$n.locus)) return()

    output$indivSelect <- renderText({
      input$selectIndiv
    })
    if (input$selectIndiv != "ALL") {
      output$maxIndivPage <- renderText({
        "1"
      })
      indivPg$maxPg <- 1
      updateNumericInput(session, "indivPage", value = 1, max = 1)
      indivPg$i <- input$selectIndiv
      ranges$y <- c(0, 2)
    }
    else {
      output$maxIndivPage <-
        renderText({
          paste0(ceiling(
            as.numeric(panelParam$n.indiv) /
              as.numeric(input$indivPerDisplay)
          ))
        })
      indivPg$maxPg <- ceiling(
        as.numeric(panelParam$n.indiv) /
          as.numeric(input$indivPerDisplay))
      updateNumericInput(
        session,
        "indivPage",
        value = 1,
        max = ceiling(panelParam$n.indiv / 15)
      )
      updateSelectInput(session, "indivPerDisplay", selected = 15)
      end.indx <- min(15, panelParam$n.indiv)
      indivPg$i <- panelParam$indiv.label.bare[1:end.indx]
      ranges$y <- c(0, end.indx + 1)
    }
    Min.filter.haplo()

    session$onFlushed(function(){
      session$sendCustomMessage(type = 'biplot_set', message = isolate(input$biplot_selected))
    },
    once=TRUE)
  })

  observeEvent(input$indivPerDisplay, {
    if (is.null(panelParam$n.indiv))
      return()

    if (input$selectIndiv != "ALL") {
      output$maxIndivPage <- renderText({
        "1"
      })
      indivPg$maxPg <- 1
      updateNumericInput(session, "indivPage", value = 1, max = 1)
      indivPg$i <- input$selectIndiv
      ranges$y <- c(0, 2)
    }
    else {
      if (input$indivPerDisplay == 100) {
        output$maxIndivPage <- renderText({
          "1"
        })
        indivPg$maxPg <- 1
        updateNumericInput(session,
                           "indivPage",
                           value = 1,
                           max = 1)
        indivPg$i <- panelParam$indiv.label.bare
        ranges$y <- c(0, length(indivPg$i) + 1)
      }
      else{
        output$maxIndivPage <-
          renderText({
            paste0(ceiling(
              as.numeric(panelParam$n.indiv) /
                as.numeric(input$indivPerDisplay)
            ))
          })
        indivPg$maxPg <- ceiling(
          as.numeric(panelParam$n.indiv) /
            as.numeric(input$indivPerDisplay)
        )
        updateNumericInput(session, "indivPage", max = ceiling(
          as.numeric(panelParam$n.indiv)
          / as.numeric(input$indivPerDisplay)
        ))

        end.indx <-
          min(as.numeric(panelParam$n.indiv) ,
              as.numeric(input$indivPerDisplay))
        indivPg$i <- panelParam$indiv.label.bare[1:end.indx]
        ranges$y <- c(0, length(indivPg$i) + 1)
      }
    }

  })

  observeEvent(input$indivPage, { #updateIndivSizeDisplay
    if(is.null(input$indivPage) | is.na(input$indivPage) |!is.numeric(input$indivPage)) return()
    if(input$indivPage <=0){
      updateNumericInput(session, "indivPage", value = 1)
      #return()
    }
    if(input$indivPage > indivPg$maxPg) {
      updateNumericInput(session, "indivPage", value = indivPg$maxPg)
    }

    if (input$selectIndiv == "ALL") {
      if (input$indivPerDisplay == 100) {
        indivPg$i <- panelParam$indiv.label.bare
        ranges$y <- c(0, length(indivPg$i) + 1)

      }
      else {
        pg <-
          min(as.numeric(input$indivPage),
              ceiling(
                as.numeric(panelParam$n.indiv) /
                  as.numeric(input$indivPerDisplay)
              ))
        start.indx <-
          (as.numeric(input$indivPerDisplay) * (pg - 1)) + 1
        end.indx <-
          min(as.numeric(panelParam$n.indiv) ,
              as.numeric(input$indivPerDisplay) * pg)
        indivPg$i <-
          panelParam$indiv.label.bare[start.indx:end.indx]
        ranges$y <- c(0, length(indivPg$i) + 1)

      }
    }
  })

  observeEvent(input$submitSrMicroHap, {
    if (input$selectLocus != "ALL") {
      srhapPg$makePlot <- TRUE
      srhapPg$num.iter <- as.numeric(input$gibbIter)
      srhapPg$frac.burn <- input$fracBurn
      srhapPg$random.seed <- input$randomSeed
      srhapPg$prior.model <- input$selectPrior
      srhapPg$locus.select <- input$selectLocus
      srhapPg$min.read.depth <- input$coverageMin
      srhapPg$data.table <- update.Haplo.file() %>% filter(locus == input$selectLocus)

      Run.SrMicrohap()
    }
  })


  # observe events on locus annotation tab
  observeEvent(input$annotateSave, {
    if(is.null(annotateTab$tbl)) return()

    match.indx <- which(annotateTab$tbl$locus == input$selectLocus)

    # if(input$rewriteFilter[1]) {
    #   annotateTab$tbl$min.rd[match.indx] <-  filterParam$minRD
    #   annotateTab$tbl$min.ar[match.indx] <-  filterParam$minAR
    #   annotateTab$tbl$n.alleles[match.indx] <-  filterParam$n.alleles
    #   #filterParam$max.ar.hm <- ifelse(input$max.ar.hm > filterParam$minAR, filterParam$minAR, input$max.ar.hm)
    #   #filterParam$min.ar.hz <- ifelse(input$min.ar.hz < filterParam$minAR, filterParam$minAR, input$min.ar.hz)
    # }

    annotateTab$tbl$status[match.indx] <- input$locusAccept
    annotateTab$tbl$comment[match.indx] <- input$locusComment[1]

    annotate.file <- strsplit(input$selectDB, split=".rds") %>% unlist %>% paste0(.,"_annotate.rds")
    saveRDS(annotateTab$tbl, annotate.file)

  })

  # observe events on filter save tab
  observeEvent(input$filterSave, {
    if(is.null(annotateTab$tbl)) return()
    if(is.na(input$coverageMin) || input$coverageMin <0) {
      createAlert(session, "alert", "filterAlert", title = "Invalid input",
                  content = "Minimum read coverage must be a positive integer", append = FALSE)
      return()
    }

    closeAlert(session, "filterAlert")


    match.indx <- which(annotateTab$tbl$locus == input$selectLocus)

    if(filterParam$minRD == input$coverageMin && filterParam$minAR == input$minAlleleRatio &&
       filterParam$n.alleles == input$n.alleles && filterParam$minRD == annotateTab$tbl$min.rd[match.indx] &&
       filterParam$minAR == annotateTab$tbl$min.ar[match.indx] &&
       filterParam$n.alleles == annotateTab$tbl$n.alleles[match.indx]  &&
       filterParam$max.ar.hm == annotateTab$tbl$max.ar.hm[match.indx] &&
       filterParam$min.ar.hz ==  annotateTab$tbl$min.ar.hz[match.indx]) {
      return()
    }

    filterParam$minRD <- input$coverageMin
    filterParam$minAR <- input$minAlleleRatio
    filterParam$n.alleles <- input$n.alleles
    #filterParam$max.ar.hm <- ifelse(input$max.ar.hm > filterParam$minAR, filterParam$minAR, input$max.ar.hm)
    #filterParam$min.ar.hz <- ifelse(input$min.ar.hz < filterParam$minAR, filterParam$minAR, input$min.ar.hz)
    filterParam$opts <- input$filterOpts

    annotateTab$tbl$min.rd[match.indx] <-  filterParam$minRD
    annotateTab$tbl$min.ar[match.indx] <-  filterParam$minAR
    annotateTab$tbl$n.alleles[match.indx] <- filterParam$n.alleles
    annotateTab$tbl$max.ar.hm[match.indx] <- filterParam$max.ar.hm
    annotateTab$tbl$min.ar.hz[match.indx] <- filterParam$min.ar.hz

    updateSliderInput(session, "max.ar.hm", max = filterParam$minAR, value= filterParam$max.ar.hm)
    updateSliderInput(session, "min.ar.hz", min = filterParam$minAR, value = filterParam$min.ar.hz)

    #removePopover(session, "filterSave")
    # addPopover(session, "filterSave","Saved Values (current)",
    #            content=paste0("<b>n.alleles</b>: ", annotateTab$tbl$n.alleles[match.indx], "<br>",
    #                           "<b>min.rd</b>: ", annotateTab$tbl$min.rd[match.indx], "<br>",
    #                           "<b>min.ar</b>: ", annotateTab$tbl$min.ar[match.indx], "<br><br>",
    #                           "<b>ar.hm</b>: ", annotateTab$tbl$max.ar.hm[match.indx], "<br>",
    #                           "<b>ar.hz</b>: ", annotateTab$tbl$min.ar.hz[match.indx]),
    #            trigger = "hover",  placement = "bottom")


    annotate.file <- strsplit(input$selectDB, split=".rds") %>% unlist %>% paste0(.,"_annotate.rds")
    saveRDS(annotateTab$tbl, annotate.file)

    Filter.haplo.by.RDnAR()
  })


  haplo.summaryTbl <- reactive({
    haplo.sum <- Filter.haplo.sum() %>% ungroup() %>% filter(rank <=2)
    if (is.null(haplo.sum))
      return ()

    ### report only diploid-based summary!!!

    ar <- as.numeric((haplo.sum$rank>1)*haplo.sum$allele.balance)

    call.indx <- which(haplo.sum$allele.balance >= haplo.sum$min.ar.hz)

    toss.indiv.indx <- which( (ar >= haplo.sum$min.ar)*(ar < haplo.sum$min.ar.hz) +
                                (ar <= haplo.sum$min.ar)*(ar > haplo.sum$max.ar.hm) >0)


    toss.indiv <- haplo.sum[toss.indiv.indx,] %>% select(id, locus) %>% unique()

    haplo.all <- anti_join(haplo.sum[call.indx, ], toss.indiv, by=c("id", "locus")) %>% ungroup() %>%
      select(group, locus, id, haplo, depth, allele.balance, rank)

    haplo.1 <- haplo.all %>%
      rename("haplotype.1"=haplo, "read.depth.1" = depth, "ar" = allele.balance) %>%
      filter(rank == 1)

    if(nrow(haplo.1)== 0) return()
    haplo.1$ar <- 1

    haplo.joined <- full_join(haplo.1 %>% select(-rank),
                              haplo.all %>% rename("haplotype.2"=haplo, "read.depth.2" = depth, "ar.x" = allele.balance) %>% filter(rank == 2) %>% select(-rank),
                              by=c("group","locus", "id"))

    homoz.indx <- which(is.na(haplo.joined$haplotype.2))
    het.indx <- which(!is.na(haplo.joined$haplotype.2))

    haplo.joined$haplotype.2[homoz.indx] <- haplo.joined$haplotype.1[homoz.indx]
    haplo.joined$read.depth.2[homoz.indx] <- haplo.joined$read.depth.1[homoz.indx]

    haplo.joined$ar[het.indx] <- haplo.joined$ar.x[het.indx]

    haplo.joined %>% select(group, locus, id, haplotype.1, haplotype.2, read.depth.1, read.depth.2, ar)


  })

  haplo.freqTbl <- reactive({
    if (is.null(haplo.summaryTbl())) {
      return()
    }

    obs.freq.tbl <-  haplo.summaryTbl() %>%
      ungroup() %>%
      group_by(locus) %>%
      mutate(tot.haplo = n()) %>%
      group_by(locus, haplotype.1, haplotype.2) %>%
      summarise(obs.freq = n() / tot.haplo[1])

    expect.freq.tbl <-
      gather(obs.freq.tbl, whichHap,  hap1, 2:3) %>%
      group_by(locus, hap1) %>%
      summarise(n = sum(obs.freq / 2)) %>%
      mutate(hap2 = hap1, n1 = n) %>%
      expand(., nesting(hap1, n), nesting(hap2, n1)) %>%
      mutate(expected.freq = ifelse(hap1 == hap2, n * n1, 2 * n * n1)) %>%
      group_by(locus, hap1, hap2) %>%
      mutate(haplotype.1 = sort(c(hap1, hap2))[1],
             haplotype.2 = sort(c(hap1, hap2))[2]) %>%
      ungroup() %>%
      select(locus, haplotype.1, haplotype.2, expected.freq) %>%
      distinct()

    inner_join(obs.freq.tbl,
               expect.freq.tbl,
               by = c("locus", "haplotype.1", "haplotype.2"))
  })


  Min.filter.haplo <- reactive({
    haplo.filter <- update.Haplo.file()
    if (is.null(haplo.filter))
      return ()
    if (input$selectGroup != "ALL")
      haplo.filter <- haplo.filter %>% filter(group == input$selectGroup)
    if (input$selectLocus != "ALL")
      haplo.filter <- haplo.filter %>% filter(locus == input$selectLocus)
    if (input$selectIndiv != "ALL")
      haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv)

    hapPg$width <- haplo.filter %>% filter(rank <= 2) %>% select(haplo) %>% unique %>% unlist %>% length

    haplo.filter

  })

  Filter.haplo.by.RDnAR <- reactive({
    haplo.sum <- Min.filter.haplo()

    #message( "successfull pass on haplo.sum --- ", dim(haplo.sum),"=--\n")

    if (is.null(haplo.sum) || dim(haplo.sum)[1] ==0)
      return ()

    haplo.join.ar <- left_join(haplo.sum,
                               annotateTab$tbl,
                               by="locus") %>%
      group_by(locus)

    ## override if the existed value is below the minimal baseline
    if ("1" %in% filterParam$opts) {
      haplo.join.ar <- haplo.join.ar %>% mutate(min.rd = ifelse(filterParam$minRD>min.rd,
                                                                filterParam$minRD,
                                                                min.rd),
                                                min.ar= ifelse(filterParam$minAR>min.ar,
                                                               filterParam$minAR,
                                                               min.ar),
                                                n.alleles = ifelse(filterParam$n.alleles>n.alleles,
                                                                   filterParam$n.alleles,
                                                                   n.alleles))

      annotateTab$tbl <- annotateTab$tbl %>% mutate(min.rd = ifelse(filterParam$minRD>min.rd,
                                                                    filterParam$minRD,
                                                                    min.rd),
                                                    min.ar= ifelse(filterParam$minAR>min.ar,
                                                                   filterParam$minAR,
                                                                   min.ar),
                                                    n.alleles = ifelse(filterParam$n.alleles>n.alleles,
                                                                       filterParam$n.alleles,
                                                                       n.alleles))
    }
    # overriding all
    if ("2" %in% filterParam$opts) {
      haplo.join.ar <- haplo.join.ar %>% mutate(min.rd = filterParam$minRD,
                                                min.ar= filterParam$minAR,
                                                n.alleles = filterParam$n.alleles)

      annotateTab$tbl <- annotateTab$tbl %>% mutate(min.rd = filterParam$minRD,
                                                    min.ar= filterParam$minAR,
                                                    n.alleles = filterParam$n.alleles)
    }

    # using the general broad stroke
    if (! "3" %in% filterParam$opts)
      haplo.join.ar <- haplo.join.ar %>%
      mutate(min.rd = filterParam$minRD,
             min.ar = filterParam$minAR,
             n.alleles = filterParam$n.alleles)

    haplo.join.ar <- haplo.join.ar %>% ungroup() %>%
      filter(allele.balance >= min.ar,
             rank <= n.alleles) %>%
      group_by(group, id, locus) %>%
      mutate(tot.depth = sum(depth)) %>%
      ungroup() %>%
      filter(tot.depth >= min.rd) %>%
      select(-tot.depth) %>%
      ungroup()

    haplo.join.ar
  })

  # this reactive fn only used for Quality profiling opt (when min RD is replaced)
  Filter.haplo.by.RDhapnAR <- reactive({
    haplo.sum <- Min.filter.haplo()

    if (is.null(haplo.sum))
      return ()

    haplo.join.ar <- left_join(haplo.sum,
                               annotateTab$tbl,
                               by="locus") %>%
      group_by(locus)

    ## override if the existed value is below the minimal baseline
    if ("1" %in% filterParam$opts) {
      haplo.join.ar <- haplo.join.ar %>% mutate(min.rd = ifelse(filterParam$minRD>min.rd,
                                                                filterParam$minRD,
                                                                min.rd),
                                                min.ar= ifelse(filterParam$minAR>min.ar,
                                                               filterParam$minAR,
                                                               min.ar),
                                                n.alleles = ifelse(filterParam$n.alleles>n.alleles,
                                                                   filterParam$n.alleles,
                                                                   n.alleles))
    }
    # overriding all
    if ("2" %in% filterParam$opts) {
      haplo.join.ar <- haplo.join.ar %>% mutate(min.rd = filterParam$minRD,
                                                min.ar= filterParam$minAR,
                                                n.alleles = filterParam$n.alleles)
    }

    # using the general broad stroke
    if (! "3" %in% filterParam$opts)
      haplo.join.ar <- haplo.join.ar %>%
      mutate(min.rd = filterParam$minRD,
             min.ar = filterParam$minAR,
             n.alleles = filterParam$n.alleles)

    haplo.join.ar <- haplo.join.ar %>% ungroup() %>%
      filter(allele.balance >= min.ar,
             depth >= filterParam$minRDhap) %>%
      ungroup()

    haplo.join.ar
  })

  Filter.all <- reactive({
    haplo.sum <- update.Haplo.file()

    if (is.null(haplo.sum))
      return ()

    haplo.join.ar <- left_join(haplo.sum,
                               annotateTab$tbl,
                               by="locus") %>% ungroup() %>%
      filter(rank <= n.alleles) %>%
      arrange(group, id, locus, rank) %>%
      group_by(group, id, locus) %>%
      summarise(
        tot.depth = sum(depth),
        categ = ifelse(length(rank) > 1  &&
                         allele.balance[2] >= min.ar,
                       "Het","Homoz"),
        ar = ifelse(categ=="Het", allele.balance[2],
                    ifelse(length(rank)==1,0,allele.balance[2])),
        finalCall = ifelse(tot.depth < min.rd,
                           "NoCall",
                           ifelse(categ=="Het",
                                  ifelse(ar >= min.ar.hz,"call", "NoCall"),
                                  ifelse(ar <= max.ar.hm,"call", "NoCall"))))
  })


  Filter.haplo.sum <- reactive({

    haplo.filter <- Filter.haplo.by.RDnAR()

    if (is.null(haplo.filter) || nrow(haplo.filter) == 0)
      return ()

    Update.ave.entropy()
    groupPg$width <- dim(haplo.filter)[1]
    hapPg$width <- haplo.filter %>% filter(rank <= 2) %>% select(haplo) %>% unique %>% unlist %>% length

    #message( "successfull pass on filter.haplo.sum --- ", dim(haplo.filter),"=--\n")

    haplo.filter
  })

  Get.tbl.by.locus <- reactive({
    if (is.null(Filter.haplo.sum()))
      return()
    haplo.ct <- Filter.haplo.sum() %>%
      group_by(locus, id, group) %>%
      summarise(tot.hapl = n(), tot.depth = sum(depth))

  })

  Get.tbl.by.id <- reactive({
    if (is.null(Filter.haplo.sum()))
      return()
    haplo.ct <- Filter.haplo.sum() %>%
      group_by(id, locus) %>%
      summarise(tot.depth = sum(depth))
  })



  # BY LOCUS PANEL::

  output$haplDensityPlot <- renderPlot({
    if (is.null(input$selectLocus) ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l))
      return ()

    if (is.null(Get.tbl.by.locus()))
      return()
    if (dim(panelParam$locus.label.tbl)[1] == 0)
      return()

    haplo.tot.tbl <- Get.tbl.by.locus() %>%
      group_by(locus, tot.hapl) %>%
      summarise(ct = n()) %>%
      ungroup() %>%
      group_by(locus) %>%
      mutate(frac = ct / sum(ct))

    #message( "is it updating_", unlist(haplo.tot.tbl[1,]), "_----\n")



    uniqH.perI.tbl <-
      right_join(haplo.tot.tbl, panelParam$locus.label.tbl, by = "locus")
    if (is.null(uniqH.perI.tbl))
      return()
    uniqH.perI.tbl[is.na(uniqH.perI.tbl)] <- 0

    #message( "is it moving :_ _----\n")


    if (input$selectLocus != "ALL") {
      uniqH.perI.tbl <-
        uniqH.perI.tbl %>% filter(locus == input$selectLocus)
    }


    max.haplo <- max(uniqH.perI.tbl$tot.hapl)

    ggplot() +
      geom_point(data = uniqH.perI.tbl, aes(
        x = tot.hapl,
        y = locus,
        size = frac,
        color = frac
      )) +
      #scale_x_log10()+
      xlab("num of uniq. haplotypes \n(per indiv)") +
      ylab("Locus ID") +
      scale_color_continuous(guide = FALSE) + #"fraction")+
      scale_size_continuous(guide = FALSE) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0,colour = "white"),
        axis.line.y = element_line(color="grey", size = 0.5),
        plot.margin = unit(c(0, 3, 0, 0), "mm")
      ) +
      #scale_x_discrete(breaks= pretty_breaks())+
      ylim(locusPg$l) +
      scale_x_continuous(limits=c(0,max.haplo+1),breaks=round(seq(1,max.haplo, length.out=4)))+
      coord_cartesian(ylim = rangesH$y)
  }, height = function()
    ifelse(groupPg$width == 0, 0,
           as.numeric(input$lociHeight) * max(
             ifelse(
               input$selectLocus == "ALL",
               ifelse(
                 input$locusPerDisplay == 100,
                 10 * length(panelParam$locus.label),
                 10 * as.numeric(input$locusPerDisplay)
               ),
               1
             ), 400
           )))

  output$numHapPlot <- renderPlot({
    if (is.null(input$selectLocus) ||
        is.null(input$selectIndiv) || is.null(locusPg$l))
      return ()
    if (is.null(haplo.summaryTbl())) {
      return()
    }
    if (dim(panelParam$locus.label.tbl)[1] == 0)
      return()


    frac.calleable <-
      haplo.summaryTbl() %>% group_by(locus) %>% summarise(n = length(unique(c(
        haplotype.1, haplotype.2
      ))))

    frac.calleable <-
      right_join(frac.calleable, panelParam$locus.label.tbl, by = "locus")
    if (is.null(frac.calleable))
      return()

    frac.calleable[is.na(frac.calleable)] <- 0

    if (input$selectLocus != "ALL") {
      frac.calleable <-
        frac.calleable %>% filter(locus == input$selectLocus)
    }

    max.haplo <- max(frac.calleable$n)

    ggplot(frac.calleable, aes(x = n, y = locus)) +
      geom_point() +
      xlab("num of unique haplotypes\n(total indiv)") +
      ylab("") +
      scale_size_continuous(guide = FALSE) + #"fraction")+
      theme_bw() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0,colour = "white"),
        plot.margin = unit(c(0, 2, 0, 0), "mm")
      ) +
      ylim(locusPg$l) +
      scale_x_continuous(limits=c(0,max.haplo+1),breaks=round(seq(1,max.haplo, length.out=4)))+
      coord_cartesian(ylim = rangesH$y)
    #scale_x_discrete(limits=c(-1, max(frac.calleable$n)+1)) #breaks= pretty_breaks()
  }, height = function()
    ifelse(groupPg$width == 0, 0, as.numeric(input$lociHeight) * max(
      ifelse(
        input$selectLocus == "ALL",
        ifelse(
          input$locusPerDisplay == 100,
          10 *
            length(panelParam$locus.label),
          10 *
            as.numeric(input$locusPerDisplay)
        ),
        1
      ), 400
    )))

  output$fracIndivPlot <- renderPlot({
    if (is.null(input$selectLocus) ||
        is.null(input$selectIndiv) || is.null(locusPg$l))
      return ()
    if (is.null(haplo.summaryTbl())) {
      return()
    }
    if (dim(panelParam$locus.label.tbl)[1] == 0)
      return()

    nIndiv <-
      ifelse(input$selectIndiv == "ALL", panelParam$n.indiv, 1)

    frac.calleable <-
      haplo.summaryTbl() %>% group_by(locus) %>% summarise(f = n() / nIndiv)
    frac.calleable <-
      right_join(frac.calleable, panelParam$locus.label.tbl, by = "locus")
    if (is.null(frac.calleable))
      return()

    frac.calleable[is.na(frac.calleable)] <- 0

    if (input$selectLocus != "ALL") {
      frac.calleable <-
        frac.calleable %>% filter(locus == input$selectLocus)
    }

    frac.calleable <- frac.calleable %>% ungroup() %>% mutate(label.x = ifelse(f>0.5, f-0.3, f+0.3))

    #message( "checking", frac.calleable %>% filter(abs(f-0.5)<0.1) %>% select(f, label.x) %>% unlist(), "_----\n")

    ggplot(frac.calleable, aes(x = f, y = locus, color = f)) +
      geom_point() +
      geom_text(data=frac.calleable, aes(y=locus, x=label.x, label=round(f,2)), color="black")+
      xlab("fraction of indiv with\n calleable haplotype") +
      ylab("") +
      scale_color_continuous(guide = FALSE) + #"fraction")+
      theme_bw() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0,colour = "white"),
        plot.margin = unit(c(0, 3, 0, 0), "mm")
      ) +
      ylim(locusPg$l) +
      coord_cartesian(ylim = rangesH$y) +
      #xlim(c(0, 1))
      scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.2))

  }, height = function()
    ifelse(groupPg$width == 0, 0,
           as.numeric(input$lociHeight) * max(
             ifelse(
               input$selectLocus == "ALL",
               ifelse(
                 input$locusPerDisplay == 100,
                 10 * length(panelParam$locus.label),
                 10 * as.numeric(input$locusPerDisplay)
               ),
               1
             ), 400
           )))

  output$readDepthPerLocus <- renderPlot({
    if (is.null(input$selectLocus) ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l))
      return ()

    if (is.null(Get.tbl.by.locus()))
      return()
    if (dim(panelParam$locus.label.tbl)[1] == 0)
      return()

    readDepth.perI.tbl <-
      right_join(Get.tbl.by.locus(), panelParam$locus.label.tbl, by = "locus")
    if (is.null(readDepth.perI.tbl))
      return()
    readDepth.perI.tbl[is.na(readDepth.perI.tbl)] <- 0

    if (input$selectLocus != "ALL") {
      readDepth.perI.tbl <-
        readDepth.perI.tbl %>% filter(locus == input$selectLocus)
    }

    readDepth.perI.tbl <-
      readDepth.perI.tbl %>% group_by(locus) %>% mutate(mean.depth = mean(tot.depth))

    ggplot(readDepth.perI.tbl, aes(x = locus, y = tot.depth)) +
      xlab("") +
      ylab("tot read depth \n(per indiv)") +
      geom_violin() +
      geom_point(aes(x = locus, y = mean.depth),
                 cex = 3,
                 pch = 3) +
      #geom_point()+
      theme_bw() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0,colour = "white"),
        plot.margin = unit(c(0, 0, 0, 0), "mm")
      ) +
      scale_y_log10() +
      xlim(locusPg$l) +
      coord_flip(xlim = rangesH$y)
  }, height = function()
    ifelse(groupPg$width == 0, 0, as.numeric(input$lociHeight) * max(
      ifelse(
        input$selectLocus == "ALL",
        ifelse(
          input$locusPerDisplay == 100,
          10 *
            length(panelParam$locus.label),
          10 *
            as.numeric(input$locusPerDisplay)
        ),
        1
      ), 400
    )))



  ## BY INDIVIDUAL PANEL::

  output$AlleleRatioByIndiv <- renderPlot({
    if (is.null(input$selectLocus) ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l))
      return ()

    haplo.filter <- Filter.haplo.sum() %>%
      filter(allele.balance >= min.ar) %>% ungroup()


    if (is.null(haplo.filter) || dim(haplo.filter)[1] == 0)
      return ()
    if (dim(panelParam$indiv.label.tbl)[1] == 0)
      return()

    # haplo.filter <- haplo.filter %>%
    #   group_by(locus, id, group) %>%
    #   mutate(depth.ratio = ifelse(length(depth) == 1, 1, min(allele.balance)),
    #             depth.first = max(depth))

    #message( "should be safe \n")

    haplo.filter <-
      right_join(haplo.filter, panelParam$indiv.label.tbl, by = "id")
    if (is.null(haplo.filter))
      return()
    haplo.filter[is.na(haplo.filter)] <- 0

    #if (input$selectIndiv != "ALL") {
    #  haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv)
    #}


    ggplot(data = haplo.filter, aes(
      x = allele.balance, #depth.ratio,
      y = id,
      size = log(depth, 10), #depth.first
      color = group
    )) +
      geom_point(alpha = 0.4) +
      scale_size_continuous(guide = FALSE) + #"Read Depth of the most common haplotype (log 10)")+
      scale_color_discrete(
        guide = FALSE,
        drop = TRUE,
        limits = levels(panelParam$group.label.bare)
      ) +
      theme_bw() +
      ylab("individual ID") +
      xlab ("allelic balance ratio") +
      theme(
        legend.position = "bottom",
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0, colour = "white"),
        axis.line.y = element_line(color="grey", size = 0.5),
        plot.margin = unit(c(0, 2, 0, 0), "mm")
      ) +
      xlim(c(0, 1)) +
      ylim(indivPg$i) +
      coord_cartesian(ylim = ranges$y) +
      geom_vline(
        xintercept = filterParam$minAR,
        linetype = "dashed",
        color = "red"
      )
  }, height = function()
    ifelse(groupPg$width == 0, 0,
           as.numeric(input$indivHeight) * max(
             ifelse(
               input$selectIndiv == "ALL",
               ifelse(
                 input$indivPerDisplay == 100,
                 10 * length(panelParam$indiv.label),
                 10 * as.numeric(input$indivPerDisplay)
               ),
               1
             ), 400
           )))



  output$numUniqHapByIndiv <- renderPlot({
    if (is.null(input$selectLocus) ||
        is.null(input$selectIndiv) || is.null(indivPg$i))
      return ()
    if (dim(panelParam$locus.label.tbl)[1] == 0)
      return()

    if(is.null(Filter.haplo.sum())) return()

    filter.haplo <- Filter.haplo.sum() %>%
      filter(allele.balance >= min.ar) %>% ungroup() %>%
      group_by(id, locus) %>%
      summarise(n.hap.locus = n()) %>%
      ungroup() %>%
      group_by(id, n.hap.locus) %>%
      summarise(n.locus = n())


    tot.hap.per.indiv <-
      right_join(filter.haplo, panelParam$indiv.label.tbl, by = "id")
    if (is.null(tot.hap.per.indiv))
      return()

    tot.hap.per.indiv[is.na(tot.hap.per.indiv)] <- 0

    max.locus <- max(tot.hap.per.indiv$n.hap.locus)

    ggplot(tot.hap.per.indiv, aes(x = n.hap.locus, y = id, size = n.locus)) +
      geom_point() +
      xlab("num of haplotype\n(per locus)") +
      ylab("") +
      scale_size_continuous(guide = FALSE) + #"fraction")+
      theme_bw() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing = unit(c(0,0,0,0), 'mm'),
        panel.border = element_rect(size = 0,colour = "white"),
        plot.margin = unit(c(0, 3, 0, 0), "mm")
      ) +
      ylim(indivPg$i) +
      scale_x_continuous(limits=c(0,max.locus+1),breaks=round(seq(1,max.locus, length.out=4)))+
      coord_cartesian(ylim = ranges$y)
    #scale_x_discrete(limits=c(-1, max(frac.calleable$n)+1)) #breaks= pretty_breaks()
  }, height = function()
    ifelse(groupPg$width == 0, 0,
           as.numeric(input$indivHeight) * max(
             ifelse(
               input$selectIndiv == "ALL",
               ifelse(
                 input$indivPerDisplay == 100,
                 10 * length(panelParam$indiv.label),
                 10 * as.numeric(input$indivPerDisplay)
               ),
               1
             ), 400
           )))

  output$fracHaploPlot <- renderPlot({
    if (is.null(input$selectLocus) ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l))
      return ()
    if (is.null(haplo.summaryTbl())) {
      return()
    }

    nLocus <-
      ifelse(input$selectLocus == "ALL", panelParam$n.locus, 1)

    haplo.filter <-
      haplo.summaryTbl() %>% group_by(id) %>% summarise(f = n() / nLocus)
    if (dim(panelParam$indiv.label.tbl)[1] == 0)
      return()
    haplo.filter <-
      right_join(haplo.filter, panelParam$indiv.label.tbl, by = "id")
    if (is.null(haplo.filter))
      return()

    haplo.filter[is.na(haplo.filter)] <- 0

    #if (input$selectIndiv != "ALL") {
    #  haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv)
    #}
    haplo.filter <- haplo.filter %>% ungroup() %>% mutate(label.x = ifelse(f>0.5, f-0.3, f+0.3))


    ggplot(haplo.filter, aes(x = f, y = id, color = f)) +
      geom_point() +
      geom_text(data=haplo.filter, aes(y=id, x=label.x, label=round(f,2)), color="black")+
      xlab("fraction of calleable\n haplotypes") +
      ylab("") +
      scale_color_continuous(guide = FALSE) + #"fraction")+
      theme_bw() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0,colour = "white"),
        plot.margin = unit(c(0,4, 0, 0), "mm")
      ) +
      #plot.margin = unit(c(0, 0, 0, 0), "mm"))+
      ylim(indivPg$i) +
      coord_cartesian(ylim = ranges$y) +
      #xlim(c(0, 1))
      scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.2))

  }, height = function()
    ifelse(groupPg$width == 0, 0,
           as.numeric(input$indivHeight) * max(
             ifelse(
               input$selectIndiv == "ALL",
               ifelse(
                 input$indivPerDisplay == 100,
                 10 * length(panelParam$indiv.label),
                 10 * as.numeric(input$indivPerDisplay)
               ),
               1
             ), 400
           )))

  # output$meanReadDepthByIndiv <- renderPlot({
  #   if (is.null(input$selectLocus) || is.null(input$selectIndiv)|| input$selectDB == "" || is.null(input$selectDB) || is.null(locusPg$l) )
  #     return ()
  #
  #   if(is.null(Get.tbl.by.id())) return()
  #
  #   haplo.filter <- Get.tbl.by.id() %>%
  #     ungroup()%>%
  #     group_by(id) %>%
  #     summarise(mean.depth = mean(tot.depth))
  #
  #   if(dim(panelParam$indiv.label.tbl)[1]==0) return()
  #   haplo.filter <- right_join(haplo.filter, panelParam$indiv.label.tbl, by="id")
  #   if(is.null(haplo.filter)) return()
  #   haplo.filter[is.na(haplo.filter)]<- 0.0#0.0001 if turn on log 10 scale
  #
  #   #if (input$selectIndiv != "ALL") {
  #   #  haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv)
  #   #}
  #
  #   ggplot(haplo.filter, aes(y=id, x=mean.depth)) +
  #     geom_point()+
  #     ylab("")+
  #     xlab("mean locus read depth ")+
  #     theme_bw()+
  #     theme(axis.text.y=element_blank(),
  #           axis.ticks.y=element_blank(),
  #           panel.spacing = unit(0, 'mm'))+
  #     #plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  #     #scale_x_log10()+
  #     ylim(indivPg$i)+
  #     coord_cartesian(ylim=ranges$y)
  #
  # },height = function() ifelse(groupPg$width==0,0,
  #                              max(ifelse(input$selectIndiv=="ALL",
  #                                         ifelse(input$indivPerDisplay==100,
  #                                                9*length(panelParam$indiv.label),
  #                                                9*as.numeric(input$indivPerDisplay)),
  #                                         1),250)))

  output$readDepthByIndiv <- renderPlot({
    if (is.null(input$selectLocus) ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l))
      return ()
    if (is.null(Filter.haplo.sum()))
      return()
    if (dim(panelParam$indiv.label.tbl)[1] == 0)
      return()



    haplo.filter <-
      right_join(Filter.haplo.sum() %>%
                   filter(allele.balance >= min.ar) %>% ungroup(), panelParam$indiv.label.tbl, by = "id")
    if (is.null(haplo.filter))
      return()
    haplo.filter[is.na(haplo.filter)] <- 0

    haplo.filter <-
      haplo.filter %>% group_by(id) %>% mutate(mean.depth = mean(depth))

    #  if (input$selectIndiv != "ALL") {
    #    haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv)
    #  }

    ggplot(haplo.filter, aes(x = id, y = depth, group = id)) +
      xlab("") +
      ylab("haplotype read depth\n(per locus)") +
      geom_violin() +
      geom_point() +
      #geom_point(aes(x=id, y=mean.depth),pch=3, cex=3)+
      theme_bw() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0,colour = "white"),
        plot.margin = unit(c(0, 2, 0, 0), "mm")
      ) +
      #plot.margin = unit(c(0, 0, 0, 0), "mm"))+
      scale_y_log10()+
      xlim(indivPg$i) +
      coord_flip(xlim = ranges$y)

  }, height = function()
    ifelse(groupPg$width == 0, 0,
           as.numeric(input$indivHeight) * max(
             ifelse(
               input$selectIndiv == "ALL",
               ifelse(
                 input$indivPerDisplay == 100,
                 10 * length(panelParam$indiv.label),
                 10 * as.numeric(input$indivPerDisplay)
               ),
               1
             ), 400
           )))

  #   output$distPlot <- renderPlot({
  #     if (is.null(input$selectLocus) || input$selectLocus == "ALL" || is.null(input$selectIndiv))
  #       return ()
  #
  #     haplo.sample <- haplo.cutoff %>% filter(locus== input$selectLocus)
  #
  #     if (input$selectIndiv != "ALL")
  #       haplo.sample <- haplo.sample %>% filter(locus== input$selectLocus, id == input$selectIndiv)
  #     if (dim(haplo.sample)[1]==0)
  #       return ()
  #
  #     ggplot()+
  #       geom_segment(data=haplo.sample, aes(x = hapl.one.st, xend = hapl.one.end, y = id, yend = id, colour= "1"), size=2 )+
  #       geom_segment(data=haplo.sample, aes(x = hapl.three.pl.end, xend = hapl.one.st, y = id, yend = id, colour="2"), size=2 )+
  #       geom_segment(data=haplo.sample, aes(x = hapl.three.pl.st, xend = hapl.three.pl.end, y = id, yend = id, colour="3+"), size=1)+
  #       scale_x_log10()+
  #       theme_bw()+
  #       xlab("read coverage cutoff")+
  #       ylab("Individual ID")+
  #       scale_color_manual(name= "Haplotypes:", values=c("1"="light grey","2"= "#4BBA82", "3+"="#A48A82"))+
  #       theme(legend.position="bottom")+
  #       coord_cartesian(ylim=ranges$y)
  #   })
  #

  # observeEvent(input$ARplot_hover, {
  #   if (!is.null(input$ARplot_hover)) filterParam$hover.minAR = input$ARplot_hover$x
  # })


  observeEvent(input$hapRDClick,{
    if(!is.null(input$hapRDClick)){
      filterCriteriaPg$RD.detail.on = ifelse(filterCriteriaPg$RD.detail.on, 0, 1)
    }
  })

  observeEvent(input$hapARClick,{
    if(!is.null(input$hapARClick)){
      filterCriteriaPg$AR.detail.on = ifelse(filterCriteriaPg$AR.detail.on, 0, 1)
    }
  })


  observeEvent(input$ARplot_dblclick, {
    if (!is.null(input$ARplot_dblclick)){
      filterParam$minAR <- round(input$ARplot_dblclick$x,2)
      updateSliderTextInput(session, "minAlleleRatio",selected=filterParam$minAR)
      Filter.haplo.sum()
    }
  })

  # observeEvent(input$RDplot_hover, {
  #   if (!is.null(input$RDplot_hover)) filterParam$hover.minRD = input$RDplot_hover$x
  # })

  observeEvent(input$RDplot_dblclick, {
    if (!is.null(input$RDplot_dblclick)){
      filterParam$minRD <- ceiling(input$RDplot_dblclick$x)
      updateSliderTextInput(session, "coverageMin",selected=filterParam$minRD)
      Filter.haplo.sum()
    }
  })



  observeEvent(input$plot_dblclick, {
    brush <- input$plot_brush
    if (!is.null(brush)) {
      ranges$y <- c(brush$ymin, brush$ymax)
      ranges$x <- c(brush$xmin, brush$xmax)
    } else {
      ranges$y <- c(0, length(indivPg$i) + 1)
    }
  })

  #   observeEvent(input$plot1_dblclick, {
  #     brush <- input$plot1_brush
  #     if (!is.null(brush)) {
  #       ranges$y <- c(brush$ymin, brush$ymax)
  #       ranges$x <- c(brush$xmin, brush$xmax)
  #     } else {
  #       ranges$y <- NULL
  #       ranges$x <- NULL
  #     }
  #   })

  observeEvent(input$plotH_dblclick, {
    brush <- input$plotH_brush
    if (!is.null(brush)) {
      rangesH$y <- c(brush$ymin, brush$ymax)
    } else {
      rangesH$y <- c(0, length(locusPg$l) + 1)
    }
  })

  observeEvent(input$aip_dblclick, {
    brush <- input$aip_Brush
    if (!is.null(brush)) {
      ranges$aip <- c(brush$xmin, brush$xmax)
    } else {
      ranges$aip <- 0
    }
  })

  observeEvent(input$alp_dblclick, {
    brush <- input$alp_Brush
    if (!is.null(brush)) {
      ranges$alp <- c(brush$xmin, brush$xmax)
    } else {
      ranges$alp <- 0
    }
  })

  # by-group distribution panel

  # output$nIndivByGroupPlot <- renderPlot({
  #   if (is.null(input$selectLocus) ||
  #       is.null(input$selectIndiv) ||
  #       input$selectDB == "" || is.null(locusPg$l))
  #     return ()
  #
  #   filter.tbl <-
  #     Filter.haplo.sum() %>% group_by(group) %>% summarise(n.indiv = length(unique(id)))
  #
  #   ggplot(filter.tbl, aes(x = n.indiv, y = group, color = group)) +
  #     geom_point() +
  #     xlab("num of indiv") +
  #     ylab("") +
  #     scale_color_discrete(guide = FALSE) +
  #     theme_bw() +
  #     theme(
  #       #axis.text.y=element_blank(),
  #       axis.ticks.y = element_blank(),
  #       panel.spacing = unit(0, 'mm'),
  #       plot.margin = unit(c(0, 0, 0, 0), "mm")
  #     )
  # }, height = function() {
  #   ifelse(
  #     groupPg$width == 0,
  #     0,
  #     ifelse(input$selectGroup == "ALL", 100 * panelParam$n.group, 100)
  #   )
  # })

  output$fIndivByGroupPlot <- renderPlot({
    if (is.null(input$selectLocus) ||
        is.null(input$selectIndiv) || input$selectDB == "")
      return ()


    all.indiv <-
      update.Haplo.file() %>% group_by(group, locus) %>% summarise(nIndiv = ifelse(input$selectIndiv !=
                                                                                     "ALL", 1, length(unique(id))))


    filter.indiv <-
      Filter.haplo.sum() %>% group_by(group, locus) %>% summarise(fIndiv = length(unique(id)))
    frac.calleable <-
      left_join(filter.indiv, all.indiv, by = c("group", "locus")) %>% mutate(f =
                                                                                fIndiv / nIndiv)


    filter.indiv <- filter.indiv %>% group_by(group) %>% summarise(fIndiv = round(mean(fIndiv)))

    mean.f.tbl <- frac.calleable %>% group_by(group) %>%
      summarise(mean.f = mean(f, na.rm = TRUE))


    #if(is.null(frac.calleable)) return()
    #     frac.calleable[is.na(frac.calleable)]<- 0
    #
    #     if (input$selectLocus != "ALL") {
    #       frac.calleable <- frac.calleable %>% filter(locus == input$selectLocus)
    #     }


    ggplot()+
      geom_point(data=frac.calleable, aes(x = f, y = group, color = group), alpha = 0.5) +
      geom_point(
        data = mean.f.tbl,
        aes(x = mean.f, y = group),
        color = "black",
        pch = 3,
        cex = 3
      ) +
      geom_text(data=filter.indiv, aes(y=group, x=1.2, label=paste0(fIndiv," indiv")))+
      xlab("frac of indiv w/ calleable hap") +
      scale_x_continuous(limits=c(0,1.3), breaks=seq(0,1,0.2))+
      ylab("") +
      scale_color_discrete(guide = FALSE) + #"fraction")+
      theme_bw() +
      theme(
        #axis.text.y = element_blank(),
        panel.border = element_rect(size = 0,colour = "white"),
        axis.ticks.y = element_blank(),
        panel.spacing = unit(0, 'mm'),
        plot.margin = unit(c(0, 0, 0, 0), "mm")
      )
  }, height = function() {
    ifelse(
      groupPg$width == 0,
      0,
      ifelse(input$selectGroup == "ALL", 50 * panelParam$n.group, 100)
    )
  })

  # output$nLociByGroupPlot <- renderPlot({
  #   if (is.null(input$selectLocus) ||
  #       is.null(input$selectIndiv) ||
  #       input$selectDB == "" || is.null(locusPg$l))
  #     return ()
  #
  #   filter.tbl <-
  #     Filter.haplo.sum() %>% group_by(group) %>% summarise(n.locus = length(unique(locus)))
  #
  #   ggplot(filter.tbl, aes(x = n.locus, y = group, color = group)) +
  #     geom_point() +
  #     xlab("num of loci") +
  #     ylab("") +
  #     scale_color_discrete(guide = FALSE) +
  #     theme_bw() +
  #     theme(
  #       axis.text.y = element_blank(),
  #       axis.ticks.y = element_blank(),
  #       panel.spacing = unit(0, 'mm'),
  #       plot.margin = unit(c(0, 0, 0, 0), "mm")
  #     )
  # }, height = function() {
  #   ifelse(
  #     groupPg$width == 0,
  #     0,
  #     ifelse(input$selectGroup == "ALL", 100 * panelParam$n.group, 100)
  #   )
  # })

  output$fLociByGroupPlot <- renderPlot({
    if (is.null(input$selectLocus) ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" || is.null(locusPg$l))
      return ()

    all.locus <-
      update.Haplo.file() %>% group_by(group, id) %>% summarise(nLocus = ifelse(input$selectLocus !=
                                                                                  "ALL", 1, length(unique(locus))))
    filter.locus <-
      Filter.haplo.sum() %>% group_by(group, id) %>% summarise(fLocus = length(unique(locus)))
    frac.calleable <-
      left_join(filter.locus, all.locus, by = c("group", "id")) %>% mutate(f =
                                                                             fLocus / nLocus)
    filter.locus <- filter.locus %>% group_by(group) %>% summarise(fLocus = round(mean(fLocus)))


    mean.f.tbl <-
      frac.calleable %>% group_by(group) %>% summarise(mean.f = mean(f, na.rm =
                                                                       T))
    #     if(is.null(frac.calleable)) return()
    #     frac.calleable[is.na(frac.calleable)]<- 0
    #
    #     if (input$selectLocus != "ALL") {
    #       frac.calleable <- frac.calleable %>% filter(locus == input$selectLocus)
    #     }

    ggplot() +
      geom_point(data=frac.calleable, aes(x = f, y = group, color = group), alpha = 0.5) +
      geom_point(
        data = mean.f.tbl,
        aes(y = group, x = mean.f),
        color = "black",
        pch = 3,
        cex = 3
      ) +
      geom_text(data=filter.locus, aes(y=group, x=1.2, label=paste0(fLocus," indiv")))+
      xlab("frac of loci w/ calleable hap") +
      ylab("") +
      scale_color_discrete(guide = FALSE) + #"fraction")+
      scale_x_continuous(limits=c(0,1.3), breaks=seq(0,1,0.2))+
      theme_bw() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(size = 0,colour = "white"),
        panel.spacing = unit(0, 'mm'),
        plot.margin = unit(c(0, 0, 0, 0), "mm")
      )
  }, height = function() {
    ifelse(
      groupPg$width == 0,
      0,
      ifelse(input$selectGroup == "ALL", 50 * panelParam$n.group, 100)
    )
  })

  # panel that takes potential issues in calling haplotypes for individuals with more than two filtered haplotypes

  output$ambigIndivPlot <- renderPlot({

    if (is.null(input$selectLocus) ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l)) {
      return()
    }

    haplo.filter <- Filter.haplo.by.RDhapnAR()
    if(is.null(Filter.haplo.by.RDhapnAR)) {return()}

    haplo.rep <- haplo.filter %>%
      group_by(id, locus) %>%
      summarise(n.accept.haplo = 1*(max(rank)>filterParam$n.alleles)) %>%
      group_by(id) %>%
      summarise(n.loci.ambig = sum(n.accept.haplo)) %>%
      group_by(n.loci.ambig) %>%
      summarise(n.indiv = n(),
                indiv.label = ifelse(n.indiv > 10,
                                     paste0(c(id[1:10],"..."),collapse = "\n"),
                                     paste0(id,collapse = "\n")),
                label.position.y = ifelse(n.indiv%%2==1,0.8,0.77))

    max.x <- max(haplo.rep$n.loci.ambig)
    if(length(ranges$aip)==1 | is.null(ranges$aip[1])) ranges$aip <- c(0,max.x+1)

    ggplot(haplo.rep, aes(x=n.loci.ambig, y=1, size=n.indiv))+
      geom_point(aes(color=factor(n.loci.ambig)))+
      #geom_text(aes(label=n.indiv), color="white")+
      geom_text(aes(label=indiv.label, y=label.position.y), hjust="center", vjust=1, size=4)+
      xlab(paste0("total num of loci that have > ",isolate(filterParam$n.alleles)," haplotypes that meet min read depth and allelic ratio")) +
      ylab("")+
      theme_bw()+
      scale_color_discrete(guide=F)+
      scale_size_continuous(range = c(4, 18), guide = FALSE) +
      scale_x_continuous(limits=ranges$aip,
                         breaks=round(seq(ranges$aip[1],
                                          ranges$aip[2],
                                          length.out= ifelse(ranges$aip[2]-ranges$aip[1] >20,
                                                             20,
                                                             (ranges$aip[2]-ranges$aip[1])+1))),
                         minor_breaks = FALSE,
                         position = "top")+
      scale_y_continuous(limits=c(0,1.2),breaks=1,minor_breaks = NULL,labels = NULL)+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.x = element_blank(),
            axis.ticks.x = element_line(colour = "black"),
            panel.border=element_rect(size = 0,linetype = "blank"),
            plot.margin = unit(c(2, 0, 6, -5), "mm"))
  }, height = 400)

  output$indivProfileTbl <- DT::renderDataTable({
    if (is.null(input$selectLocus) ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l)) {
      return()
    }

    haplo.filter <- Filter.haplo.by.RDhapnAR()
    if(is.null(Filter.haplo.by.RDhapnAR)) {return()}

    haplo.fail <- haplo.filter %>%
      group_by(id, locus) %>%
      summarise(n.accept.haplo = 1*(max(rank)>filterParam$n.alleles)) %>%
      filter(n.accept.haplo==1)

    inner_join(haplo.filter, haplo.fail, by=c("id", "locus")) %>%
      mutate(ar = round(allele.balance, 3)) %>%
      select(id, locus, haplo, depth, ar, rank)
  })

  output$ambigLociPlot <- renderPlot({

    if (is.null(input$selectLocus) ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l)) {
      return()
    }

    haplo.filter <- Filter.haplo.by.RDhapnAR()
    if(is.null(Filter.haplo.by.RDhapnAR)) return()

    haplo.rep <- haplo.filter %>%
      #filter(depth >= filterParam$minRD,
      #       allele.balance >= filterParam$minAR) %>%
      group_by(id, locus) %>%
      summarise(n.accept.haplo = 1*(max(rank)>filterParam$n.alleles)) %>%
      group_by(locus) %>%
      summarise(n.indiv.ambig = sum(n.accept.haplo)) %>%
      group_by(n.indiv.ambig) %>%
      summarise(n.loci = n(),
                loci.label = ifelse(n.loci > 10,
                                    paste0(c(locus[1:10],"..."),collapse = "\n"),
                                    paste0(locus,collapse = "\n")),
                label.position.y = ifelse(n.loci%%2==1,0.8,0.77))

    max.x <- max(haplo.rep$n.indiv.ambig)
    if(length(ranges$alp)==1 | is.null(ranges$alp[1])) ranges$alp <- c(0,max.x+1)

    if(is.null(max.x)) return()
    ggplot(haplo.rep, aes(x=n.indiv.ambig, y=1, size=n.loci))+
      geom_point(aes(color=factor(n.indiv.ambig)))+
      #geom_text(aes(label=n.indiv), color="white")+
      geom_text(aes(label=loci.label, y=label.position.y), hjust="center", vjust=1, size=4)+
      xlab(paste0("total num of individuals that calls > ",isolate(filterParam$n.alleles)," haplotypes that meet min read depth and allelic ratio")) +
      ylab("")+
      theme_bw()+
      scale_color_discrete(guide=F)+
      scale_size_continuous(range = c(4, 18), guide = FALSE) +
      scale_x_continuous(limits=ranges$alp,
                         breaks=round(seq(ranges$alp[1],
                                          ranges$alp[2],
                                          length.out= ifelse(ranges$alp[2]-ranges$alp[1] >20,
                                                             20,
                                                             (ranges$alp[2]-ranges$alp[1])+1))),
                         minor_breaks = FALSE,
                         position = "top")+
      scale_y_continuous(limits=c(0,1.2),breaks=1,minor_breaks = NULL,labels = NULL)+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.x = element_blank(),
            axis.ticks.x = element_line(colour = "black"),
            panel.border=element_rect(size = 0,linetype = "blank"),
            plot.margin = unit(c(2, 0, 6, -5), "mm"))
  }, height = 400)

  output$uchaplabel <- renderPlot({

    if (is.null(input$selectLocus) ||
        input$selectLocus == "ALL" ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l)) {
      return()
    }

    haplo.filter <- refineDiploidCall()
    if(is.null(haplo.filter)) return()

    haplo.rep <- c(haplo.filter$hap.1, haplo.filter$hap.2) %>% unique %>% .[. != "-"] %>% unlist()

    haplo.rep.df <- data.frame(haplo=factor(haplo.rep, level=sort(haplo.rep,decreasing=T)),
                               x=0,
                               y=0)

    ggplot(haplo.rep.df, aes(x,y, fill=haplo))+
      geom_point(size=0, color="white")+
      #geom_density(adjust=0.5, color=NA)+
      facet_grid(haplo~., scales="free_y")+
      scale_x_continuous("",breaks=NULL,limits=c(0,0.1))+
      scale_fill_discrete(guide=FALSE,direction=-1)+
      scale_y_continuous("",breaks=NULL,limits=c(0,0.1))+ #breaks=1
      theme_bw()+
      theme(strip.text.y = element_text(angle =360, margin=margin(0,2,0,0)),
            strip.background  = element_rect(fill="white",size = 0,linetype = "blank"),
            panel.border=element_rect(fill="white",size = 0,linetype = "blank"),
            panel.spacing = unit(0, 'mm'),
            plot.margin = unit(c(2, 0, 6, -5), "mm"),
            aspect.ratio =1000)

  }, height = function() {
    ifelse(hapPg$width == 0 || input$selectLocus == "ALL",
           1,
           max(hapPg$width*30, 350))
  })

  output$uchapReadDepth <- renderPlot({

    if (is.null(input$selectLocus) ||
        input$selectLocus == "ALL" ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l)) {
      return()
    }

    haplo.filter <- refineDiploidCall()
    if(is.null(haplo.filter)) return()

    haplo.rep <- c(haplo.filter$hap.1, haplo.filter$hap.2) %>% unique %>% .[. != "-"] %>% unlist()

    haplo.match.rep.1 <- haplo.filter %>%
      filter(hap.1 %in% haplo.rep) %>%
      mutate(haplo=factor(hap.1, level=sort(haplo.rep,decreasing=T))) %>%
      group_by(haplo, id) %>%
      summarise(top.2 = 0.2 - 0.05*(categ[1]=="Het"),
                rank.mod = 1, depth = rd.1[1], geno.class = categ[1],
                is.incl = 1*(finalCall[1]=="call"))

    haplo.match.rep.2 <- haplo.filter %>% filter(hap.2 %in% haplo.rep) %>%
      mutate(haplo=factor(hap.2, level=sort(haplo.rep,decreasing=T))) %>%
      group_by(haplo, id) %>%
      summarise(top.2 = ifelse(rank.rel == 2, 0, -0.2) - 0.05*(categ[1]=="Het"),
                rank.mod = rank.rel,
                depth = rd.2[1], geno.class = categ[1],
                is.incl = 1*(finalCall[1]=="call" & categ == "Het"))

    haplo.match.rep <- bind_rows(haplo.match.rep.1, haplo.match.rep.2)
    if (nrow(haplo.match.rep)==0) return()

    rank.lab <- data.frame(rank.lab = c("1st ", "2nd " ,"3rd+"), pos.lab = c(0.175,-0.025,-0.215))

    ggplot()+
      geom_point(data=haplo.match.rep, aes(x=depth, y=top.2, color=factor(is.incl), pch=factor(geno.class)),size=1.8, alpha=0.6)+
      facet_grid(haplo~.)+#, scales="free_y")+
      scale_x_log10("distrib. of read depth")+
      scale_color_discrete(guide=FALSE)+
      scale_shape_discrete(guide=F)+
      #scale_y_continuous("",breaks=NULL)+ #breaks=1
      theme_bw()+
      geom_text(data= rank.lab, aes(x=1, y=pos.lab, label=rank.lab), size=3, hjust ="right")+
      # geom_vline(
      #   xintercept = filterParam$minRD,
      #   linetype = "dashed",
      #   color = "red"
      # )+
      scale_y_continuous("",limits=c(-0.3,0.3),breaks=0,minor_breaks = NULL,labels = NULL)+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            strip.text.y = element_text(angle =360,size=0, margin=margin(0,0,0,0)),
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0,colour = "white"),
            plot.margin = unit(c(2, 0, 2, 1), "mm"))

  }, height = function() {
    ifelse(hapPg$width == 0 || input$selectLocus == "ALL",
           1,
           max(hapPg$width*30, 350))
  })

  output$uchapAllelicRatio <- renderPlot({

    if (is.null(input$selectLocus) ||
        input$selectLocus == "ALL" ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l)) {
      return()
    }

    haplo.filter <- refineDiploidCall()
    if(is.null(haplo.filter)) return()

    haplo.rep <- c(haplo.filter$hap.1, haplo.filter$hap.2) %>% unique %>% .[. != "-"] %>% unlist()

    haplo.match.rep.1 <- haplo.filter %>%
      filter(hap.1 %in% haplo.rep) %>%
      mutate(haplo=factor(hap.1, level=sort(haplo.rep,decreasing=T))) %>%
      group_by(haplo, id) %>%
      summarise(top.2 = 0.2 - 0.05*(categ[1]=="Het"),
                rank.mod = 1, ar = 1, geno.class = categ[1],
                is.incl = 1*(finalCall[1]=="call"))

    haplo.match.rep.2 <- haplo.filter %>% filter(hap.2 %in% haplo.rep) %>%
      mutate(haplo=factor(hap.2, level=sort(haplo.rep,decreasing=T))) %>%
      group_by(haplo, id) %>%
      summarise(top.2 = ifelse(rank.rel == 2, 0, -0.2) - 0.05*(categ[1]=="Het"),
                rank.mod = rank.rel,
                ar = rd.2/rd.1, geno.class = categ[1],
                is.incl = 1*(finalCall[1]=="call" & categ == "Het"))

    haplo.match.rep <- bind_rows(haplo.match.rep.1, haplo.match.rep.2)
    if (nrow(haplo.match.rep)==0) return()

    ggplot()+
      geom_point(data=haplo.match.rep, aes(x=ar, y=top.2, color=factor(is.incl), pch=factor(geno.class)),size=1.8, alpha=0.6)+
      facet_grid(haplo~.)+#, scales="free_y")+
      scale_x_log10("distrib. of allelic ratio",breaks=c(0.01,0.1,0.2,0.5,1), limits=c(0.01,1))+
      scale_color_discrete(guide=FALSE)+
      scale_shape_discrete(guide=F)+
      #scale_y_continuous("",breaks=NULL)+ #breaks=1
      theme_bw()+
      geom_vline(
        xintercept = filterParam$minAR,
        linetype = "dashed",
        color = "red"
      )+
      scale_y_continuous("",limits=c(-0.3,0.3),breaks=0,minor_breaks = NULL,labels = NULL)+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            strip.text.y = element_text(angle =360,size=0, margin=margin(0,0,0,0)),
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0,colour = "white"),
            plot.margin = unit(c(2, 0, 2, 0), "mm"))

  }, height = function() {
    ifelse(hapPg$width == 0 || input$selectLocus == "ALL",
           1,
           max(hapPg$width*30, 350))
  })


  ## filter status:cutoff distribution panel

  refineDiploidCall <- reactive({

    if (is.null(input$selectLocus) ||
        input$selectLocus == "ALL" ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l) ||
        is.null(filterParam$n.alleles) || length(filterParam$n.alleles) == 0)
      return()

    haplo.filter <- Min.filter.haplo()
    if(is.null(Min.filter.haplo)) return()

    h <- haplo.filter %>%
      filter(rank <= filterParam$n.alleles) %>%
      arrange(group, id, locus, rank) %>%
      group_by(group, id, locus) %>%
      summarise(categ = ifelse(length(rank) >1 &&
                                 allele.balance[2] >= filterParam$minAR,
                               "Het","Homoz"),
                genotype = ifelse(categ=="Het",
                                  paste0(sort(c(haplo[1],haplo[2])),collapse = "/"),
                                  paste0(c(haplo[1],haplo[1]),collapse = "/") ),
                hap.1 = haplo[1],
                hap.2 =ifelse(length(rank)>1,haplo[2],"-"),
                ar = ifelse(categ=="Het", allele.balance[2],
                            ifelse(length(rank)==1,0,allele.balance[2])),
                rd.1 = depth[1],
                rd.2 = ifelse(length(rank)>1,depth[2],0),
                rank.rel = ifelse(length(rank)>1, 2, 1),
                finalCall = ifelse(sum(depth) < filterParam$minRD,
                                   "NoCall",
                                   ifelse(categ=="Het",
                                          ifelse(ar >= filterParam$min.ar.hz,"call", "NoCall"),
                                          ifelse(ar <= filterParam$max.ar.hm,"call", "NoCall"))),
                drop.haplo = ifelse(length(rank)>1 & categ == "Homoz",haplo[2], "-"))

    if (filterParam$n.alleles > 2) {

      h.2 <- left_join( haplo.filter %>%
                          filter(rank <= filterParam$n.alleles, rank > 2),
                        h, by = c("id", "locus", "group")) %>%
        arrange(group, id, locus, rank) %>%
        group_by(group, id, locus, rank) %>%
        summarise(categ = ifelse(allele.balance[1] >= filterParam$minAR,
                                 "Het","-"),
                  genotype = ifelse(categ=="Het",
                                    paste0(sort(c(hap.1[1],haplo[1])),collapse = "/"),
                                    paste0(c(hap.1[1],"-"),collapse = "/" ) ),
                  hap.1 = hap.1[1],
                  hap.2 = haplo[1],
                  ar = allele.balance[1],
                  rd.1 = rd.1[1],
                  rd.2 = depth[1],
                  rank.rel = rank[1],
                  finalCall = ifelse(finalCall[1] == "NoCall",
                                     "NoCall",
                                     ifelse(categ!="Het", "NoCall",
                                            ifelse(ar >= filterParam$min.ar.hz,"call", "NoCall"))),
                  drop.haplo = ifelse(categ!="Het",haplo[1], "-")) %>%
        ungroup() %>% select(-rank)

      h <- bind_rows(h, h.2)
    }
    h %>% select(-ar)

  })

  selected_pt <- reactive({
    if( is.null(input$biplot_selected)){
      character(0)
    } else input$biplot_selected
  })


  output$savedARstat <- renderUI({

    match.indx <- which(annotateTab$tbl$locus == input$selectLocus)

    return(HTML("<b>Saved AR for Homoz</b>: ", annotateTab$tbl$max.ar.hm[match.indx], "<br>",
                "<b>Saved AR for Het</b>: ", annotateTab$tbl$min.ar.hz[match.indx], "<br>"))
  })

  output$biPlotTbl <- DT::renderDataTable({

    if (is.null(input$selectLocus) ||
        input$selectLocus == "ALL" ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l)) {
      return()
    }

    dat <- refineDiploidCall()
    if(is.null(refineDiploidCall)) return()

    #message(paste0("Number of passeable entries: ", nrow(dat %>% filter(finalCall == "call"))))
    out <- dat[dat$id %in% selected_pt(), ]
    if( nrow(out) < 1  || length(out) == 0) return()
    row.names(out) <- NULL
    out %>% select(-drop.haplo, locus)

    # # With base graphics, we need to explicitly tell it which variables were
    # # used; with ggplot2, we don't.
    # res <- nearPoints(dat, input$biPlotClick,
    #                     threshold = 10,
    #                     maxpoints = 1,
    #                   xvar="second.depth", yvar="first.depth",
    #                     addDist = FALSE) %>% select(-locus) %>%
    #   mutate(ar= round(ar,3))

    #datatable(out %>% select(-locus) %>% mutate(ar= round(ar,3)))
  })

  output$pieCallChart <- renderPlot({
    if (is.null(input$selectLocus) ||
        input$selectLocus == "ALL" ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l)) return()


    gdepths <- refineDiploidCall()
    if(is.null(gdepths)) return()

    n.sample <- nrow(gdepths)
    pc.df <- gdepths %>% group_by(finalCall) %>%
      summarise(n=n(),
                frac=n/n.sample,
                perc = paste0(round(n*100/n.sample,1),"%"),
                call.label = ifelse(finalCall[1] == "call", "call", "no call"))

    locus.name <- isolate(input$selectLocus)
    if (nchar(locus.name) > 15) locus.name <- paste0(substring(locus.name, 0, 14),"...",collapse = "")

    ggplot(pc.df, aes(x="", y=n, fill=call.label))+
      geom_bar(width = 1, stat = "identity")+
      geom_text(aes(label = perc), position = position_stack(vjust = 0.65), size =4)+
      coord_polar("y", start=0, direction=-1)+
      scale_fill_brewer(palette="Blues", direction=-1,name=locus.name)+
      theme_minimal()+
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, face="bold")
      )+
      theme(axis.text.x=element_blank())

  },height=function(){ifelse(input$selectLocus == "ALL",1,150)})

  output$biplot <-renderggiraph({
    #output$biplot <-renderPlot({

    #renderPlot({
    if (is.null(input$selectLocus) ||
        input$selectLocus == "ALL" ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l)) {
      return()
    }

    gdepths <- refineDiploidCall()
    if(is.null(gdepths)) return()

    dataset <- gdepths %>%
      filter(rd.1 <= input$max_read_depth,
             rd.2 <= input$max_read_depth,
             finalCall != "NoCall")

    nc.dataset <- gdepths %>%
      filter(rd.1 <= input$max_read_depth,
             rd.2 <= input$max_read_depth,
             finalCall == "NoCall")



    g <- ggplot()
    if (nrow(dataset)>0) {
      g<-g+
        geom_point_interactive(data=dataset, aes(x = rd.2,
                                                 y = rd.1,
                                                 #color = factor(genotype),
                                                 fill = factor(genotype),
                                                 tooltip = paste0("id: ", factor(id), "\n",
                                                                  "categ: ",factor(categ), "\n",
                                                                  "geno: ", factor(genotype), "\n",
                                                                  "rank: ", rank.rel, "\n",
                                                                  "drop hap: ", drop.haplo),

                                                 data_id = id,
                                                 shape = factor(categ)),
                               size = 4, stroke=0.5, alpha=0.7)
    }
    # geom_point(data=dataset, aes(x = second.depth,
    #                              y = first.depth,
    #                              #color = factor(genotype),
    #                              fill = factor(genotype),
    #                              shape = factor(categ)),
    #            size = 4, stroke=0.7) +
    if (nrow(nc.dataset)>0) {
      g <- g + geom_text_interactive(data=nc.dataset, aes(x = rd.2,
                                                          y = rd.1,
                                                          tooltip = paste0("id: ", factor(id), "\n",
                                                                           "categ: Not Call\n",
                                                                           "geno: ", factor(genotype), "\n",
                                                                           "rank: ", rank.rel, "\n",
                                                                           "drop hap: ", drop.haplo),
                                                          data_id = id,
                                                          label ="x"),
                                     size = 5, alpha=0.7)
    }
    # geom_point(data=nc.dataset, aes(x = second.depth,
    #                              y = first.depth),
    #            size = 4, stroke=0.7, shape=4) +
    g <- g + geom_vline(xintercept = 0, colour = "black", size = 0.8) +
      geom_hline(yintercept = 0, colour = "black", size = 0.8) +
      scale_shape_manual("category",values = rep(c(23, 21)),"") +
      guides(fill=FALSE)+
      scale_fill_discrete("", labels=NULL)

    if(filterParam$max.ar.hm >0)
    { g <- g +
      geom_abline(intercept = 0, slope = 1/filterParam$max.ar.hm, colour = "#58cbeb",lwd=1.2)
    } else {
      g <- g +
        geom_vline(xintercept = 0, colour = "#58cbeb",lwd=1.2)
    }


    g <- g + geom_abline(intercept = 0, slope = 1/filterParam$min.ar.hz, colour = "#edd02d",lwd=1.2) +
      geom_abline(intercept = 0, slope = 1/filterParam$minAR, colour = "#f09d3e",lwd=1, lty=4) +
      geom_abline(intercept = filterParam$minRD, slope = -1, colour="#b157e8", lwd=1.2)+
      geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
      #coord_fixed(ratio = 0.4, expand = FALSE,clip = "off")+
      theme_bw() +
      theme(
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0,colour = "white"),
        axis.line.y = element_line(color="grey", size = 0.5),
        plot.margin = unit(c(0, 3, 0, 0), "mm")
      )+
      xlab("non-first read depth")+
      ylab("first read depth")

    #g

    girafe_options(girafe(code=print(g)),
                   opts_zoom(min = .5, max = 5),
                   opts_toolbar(position = "top"),
                   opts_selection(
                     type = "multiple", css = "fill:#FF0000;stroke:#FF0000;r:5;fill-opacity:1;"),
                   opts_hover(css = "fill:#DCECEE;stroke:black;cursor:pointer;fill-opacity:1"))

  }#, height = function() {
  # ifelse(hapPg$width == 0 || input$selectLocus == "ALL",
  #       1,600)}
  )


  output$allReadDepth <- renderPlot({

    if (is.null(input$selectLocus) ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l)) {
      return()
    }

    haplo.filter <- Min.filter.haplo() %>%
      filter(rank <= filterParam$n.alleles) %>%
      group_by(id, locus) %>%
      summarise(tot.depth = sum(depth))

    if(is.null(Min.filter.haplo)) return()

    ggplot(haplo.filter, aes(x=tot.depth))+
      #geom_density()+
      geom_histogram(binwidth = 0.1, boundary = -0.025)+
      scale_x_log10("distrib. of total read depth")+
      scale_y_continuous("", breaks=NULL)+
      theme_bw()+
      geom_vline(
        xintercept = filterParam$minRD,
        linetype = "dashed",
        color = "red"
      )+
      theme(strip.text.y = element_text(angle =360,size=0, margin=margin(0,0,0,0)),
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0,colour = "white"),
            #axis.line.x = element_line(color="grey", size = 0.5),
            plot.margin = unit(c(2, 0, 2, 0), "mm"))
    # geom_vline(
    #   xintercept = filterParam$hover.minRD,
    #   linetype = "dashed",
    #   color = "blue"
    # )

  }, height = 300)

  output$allAllelicRatio <- renderPlot({

    if (is.null(input$selectLocus) ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l)) {
      return()
    }

    haplo.filter <- Min.filter.haplo() %>% filter(rank <= filterParam$n.alleles)
    if(is.null(Min.filter.haplo)) return()

    ggplot(haplo.filter, aes(x=allele.balance))+
      #geom_density()+
      geom_histogram(binwidth = 0.05, boundary = -0.05)+
      scale_x_log10("distrib. of allelic ratio",#limits=c(min(haplo.filter$allele.balance),1),
                    breaks=c(0.01,0.1,0.2,0.5,1))+
      scale_y_continuous("",breaks=NULL)+
      theme_bw()+
      geom_vline(
        xintercept = filterParam$minAR,
        linetype = "dashed",
        color = "red"
      )+
      theme(strip.text.y = element_text(angle =360,size=0, margin=margin(0,0,0,0)),
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0,colour = "white"),
            plot.margin = unit(c(2, 0, 2, 1), "mm"))
    # geom_vline(
    #   xintercept = filterParam$hover.minAR,
    #   linetype = "dashed",
    #   color = "blue"
    # )
  }, height = 300)


  output$haplabel <- renderPlot({

    if (is.null(input$selectLocus) ||
        input$selectLocus == "ALL" ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l)) {
      return()
    }

    haplo.filter <- Min.filter.haplo() %>%
      filter(rank <= filterParam$n.alleles) %>%
      group_by(id, locus) %>%
      mutate(tot.depth = sum(depth)) %>%
      ungroup()

    if(is.null(Min.filter.haplo)) return()

    haplo.rep <- haplo.filter %>% filter(tot.depth >= filterParam$minRD,
                                         allele.balance >= filterParam$minAR) %>%
      select(haplo) %>% unique %>% unlist()

    if (length(haplo.rep)==0 | length(haplo.rep)>50) return()

    haplo.rep.df <- data.frame(haplo=factor(haplo.rep, level=sort(haplo.rep,decreasing=T)),
                               x=0,
                               y=0)

    ggplot(haplo.rep.df, aes(x,y, fill=haplo))+
      geom_point(size=0, color="white")+
      #geom_density(adjust=0.5, color=NA)+
      facet_grid(haplo~., scales="free_y")+
      scale_x_continuous("",breaks=NULL,limits=c(0,0.1))+
      scale_fill_discrete(guide=FALSE,direction=-1)+
      scale_y_continuous("",breaks=NULL,limits=c(0,0.1))+ #breaks=1
      theme_bw()+
      theme(strip.text.y = element_text(angle =360, margin=margin(0,2,0,0)),
            strip.background  = element_rect(fill="white",size = 0,linetype = "blank"),
            panel.border=element_rect(fill="white",size = 0,linetype = "blank"),
            panel.spacing = unit(0, 'mm'),
            plot.margin = unit(c(2, 0, 6, -5), "mm"),
            aspect.ratio =1000)

  }, height = function() {
    ifelse(hapPg$width == 0 || input$selectLocus == "ALL",
           1,
           max(hapPg$width*30, 300))
  })



  # figures: in reference to the two top most common haplotype
  output$hapReadDepth <- renderPlot({

    if (is.null(input$selectLocus) ||
        input$selectLocus == "ALL" ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l)) {
      return()
    }

    haplo.filter <- Min.filter.haplo() %>%
      filter(rank <= filterParam$n.alleles) %>%
      group_by(id, locus) %>%
      mutate(tot.depth = sum(depth)) %>%
      ungroup()

    if(is.null(Min.filter.haplo)) return()

    haplo.rep <- haplo.filter %>% filter(tot.depth >= filterParam$minRD,
                                         allele.balance >= filterParam$minAR) %>%
      select(haplo) %>% unique %>% unlist()

    if (length(haplo.rep)==0 | length(haplo.rep)>50) return()

    haplo.match.rep <- haplo.filter %>% filter(haplo %in% haplo.rep) %>%
      mutate(haplo=factor(haplo, level=sort(haplo.rep,decreasing=T))) %>%
      group_by(haplo) %>%
      mutate(n.accept.indiv = sum(tot.depth >=filterParam$minRD),
             n.reject.indiv = sum(tot.depth < filterParam$minRD)) %>%
      ungroup() %>%
      mutate(
        center.x.accept = mean(tot.depth[tot.depth>=filterParam$minRD]),
        center.x.reject = mean(tot.depth[tot.depth<filterParam$minRD]),
        center.y = 0.5+ min(n.accept.indiv/2, n.reject.indiv/2)) %>%
      group_by(haplo, id) %>%
      mutate(is.pass = 1*(tot.depth >=filterParam$minRD))


    if (nrow(haplo.match.rep)==0) return()
    # ggplot(haplo.match.rep, aes(x=depth, y=haplo, color=haplo))+
    #   #geom_density(adjust=1)+
    #   geom_point(alpha=0.5)+
    #   #stat_bin(binwidth = .1,geom="line")+
    #   scale_x_log10("distrib. of read depth")+
    #   scale_color_discrete(guide=FALSE, direction=-1)+
    #     theme_bw()+
    #     geom_vline(
    #     xintercept = filterParam$minRD,
    #     linetype = "dashed",
    #     color = "red"
    #   )+
    #   theme(panel.spacing = unit(0, 'mm'),
    #     plot.margin = unit(c(2, 0, 2, 0), "mm"))


    g <- ggplot()+
      geom_histogram(data=haplo.match.rep,
                     aes(x=depth, fill=factor(is.pass)), binwidth = 0.05, boundary = -0.025)+
      geom_point(data=haplo.match.rep, aes(x=depth, y=0), size=1.2, alpha=0.4)

    if (filterCriteriaPg$RD.detail.on)
      g <- g +
      geom_text(data=haplo.match.rep, aes(
        label=n.accept.indiv,
        x=center.x.accept,
        y=center.y), hjust=-1, vjust=-1)+
      geom_text(data=haplo.match.rep, aes(
        label=n.reject.indiv,
        x=center.x.reject,
        y=center.y), hjust=-1, vjust =-1)

    #geom_density(adjust=0.5, color=NA)+
    g + facet_grid(haplo~.)+#, scales="free_y")+
      scale_x_log10("distrib. of allelic read depth")+
      scale_fill_discrete(guide=FALSE,direction=-1)+
      scale_y_continuous("",breaks=NULL)+ #breaks=1
      theme_bw()+
      # geom_vline(
      #   xintercept = filterParam$minRD,
      #   linetype = "dashed",
      #   color = "red"
      # )+
      theme(strip.text.y = element_text(angle =360,size=0, margin=margin(0,0,0,0)),
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0,colour = "white"),
            #axis.line.y = element_line(color="grey", size = 1),
            plot.margin = unit(c(2, 0, 2, 1), "mm"))

  }, height = function() {
    ifelse(hapPg$width == 0 || input$selectLocus == "ALL",
           1,
           max(hapPg$width*30, 300))
  })

  output$hapAllelicRatio <- renderPlot({

    if (is.null(input$selectLocus) ||
        input$selectLocus == "ALL" ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l)) {
      return()
    }

    haplo.filter <- Min.filter.haplo() %>%
      filter(rank <= filterParam$n.alleles) %>%
      group_by(id, locus) %>%
      mutate(tot.depth = sum(depth)) %>%
      ungroup()

    if(is.null(Min.filter.haplo)) return()

    haplo.rep <- haplo.filter %>% filter(tot.depth >= filterParam$minRD,
                                         allele.balance >= filterParam$minAR) %>%
      select(haplo) %>% unique %>% unlist()

    if (length(haplo.rep)==0 | length(haplo.rep)>50) return()

    haplo.match.rep <- haplo.filter %>% filter(haplo %in% haplo.rep) %>%
      mutate(haplo=factor(haplo, level=sort(haplo.rep,decreasing=T))) %>%
      group_by(haplo) %>%
      mutate(n.accept.indiv = sum(allele.balance >=filterParam$minAR),
             n.reject.indiv = sum(allele.balance < filterParam$minAR)) %>%
      ungroup() %>%
      mutate(
        center.x.accept = min(mean(allele.balance[allele.balance>=filterParam$minAR]),
                              0.9),
        center.x.reject = min(mean(allele.balance[allele.balance<filterParam$minAR]),0.9),
        center.y = 0.5+ min(n.accept.indiv/2, n.reject.indiv/2))%>%
      group_by(haplo, id) %>%
      mutate(is.pass = 1*(allele.balance >=filterParam$minAR))

    if (nrow(haplo.match.rep)==0) return()

    g<- ggplot()+
      geom_histogram(data=haplo.match.rep, aes(x=allele.balance, fill=haplo), binwidth = 0.05,boundary = -0.05)+
      geom_point(data=haplo.match.rep, aes(x=allele.balance, y=0), size=1.2, alpha=0.4)
    #geom_density(adjust=0.1, color=NA)+

    if (filterCriteriaPg$AR.detail.on)
      g <- g +
      geom_text(data=haplo.match.rep, aes(
        label=n.accept.indiv,
        x=center.x.accept,
        y=center.y), hjust=2, vjust=-1)+
      geom_text(data=haplo.match.rep, aes(
        label=n.reject.indiv,
        x=center.x.reject,
        y=center.y), hjust=0, vjust =-1)


    g + facet_grid(haplo~.)+#, scales="free_y")+
      scale_x_log10("distrib. of allelic ratio",breaks=c(0.01,0.1,0.2,0.5,1), limits=c(0.01,1))+
      scale_fill_discrete(guide=FALSE, direction=-1)+
      scale_y_continuous("",breaks=NULL)+
      theme_bw()+
      geom_vline(
        xintercept = filterParam$minAR,
        linetype = "dashed",
        color = "red"
      )+
      theme(strip.text.y = element_text(angle =360,size=0, margin=margin(0,0,0,0)),
            panel.border = element_rect(size = 0,colour = "white"),
            panel.spacing = unit(0, 'mm'),
            plot.margin = unit(c(2, 0, 2, 0), "mm"))

  }, height = function() {
    ifelse(hapPg$width == 0 || input$selectLocus == "ALL",
           1,
           max(hapPg$width*30, 300))
  })

  ##ABOUT HAPLOTYPE distribution panel

  # panel display the variants connections in microhaplotype in relationship to position


  output$hapSeq <- renderPlot({
    if (is.null(input$selectLocus) ||
        input$selectLocus == "ALL" ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l))
      return ()

    haplo.filter <- Filter.haplo.sum()

    pos.file <- extract.pos.file()
    if(is.null(pos.file)) return()
    colnames(pos.file) <- c("locus", "pos")
    pos.str <- pos.file %>% filter(locus==input$selectLocus) %>% select(pos) %>% unlist
    position <- strsplit(pos.str, ",") %>% unlist %>% as.numeric

    if (is.null(haplo.filter) || nrow(haplo.filter) == 0) return()
    haplo.profile.frac <- haplo.filter %>%
      group_by(haplo) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      mutate(frac = n / sum(n))

    if (nrow(haplo.profile.frac) == 0)
      return()

    haplo.split.profile <-
      sapply(1:nrow(haplo.profile.frac), function(i) {
        char.split <- strsplit(haplo.profile.frac[i,]$haplo, "")
        #message( "character split_", unlist(char.split), "_----\n")
        n.char <- length(char.split[[1]])
        sapply(1:n.char, function(j)
          c(i, j, char.split[[1]][j], haplo.profile.frac[i,]$frac))
      }) %>%
      matrix(., ncol = 4, byrow = TRUE) %>%
      as.data.frame(stringsAsFactors = FALSE) %>%
      tbl_df()

    colnames(haplo.split.profile) <-
      c("group", "pos", "seq", "frac")
    haplo.split.profile <-
      haplo.split.profile %>% mutate(
        pos = as.numeric(position[as.numeric(pos)]),
        frac = as.numeric(frac),
        group = as.numeric(group)
      )


    g <-
      ggplot(haplo.split.profile,
             aes(
               x = pos,
               y = seq,
               group = group,
               size = frac,
               color = factor(group)
             )) +
      xlab("variant position") +
      ylab("sequence") +
      scale_size_continuous(range = c(3, 20), guide = FALSE) +
      scale_x_continuous(breaks = as.numeric(position),
                         labels = as.numeric(position))

    if (length(unique(haplo.split.profile$pos)) == 1) {
      g <- g + geom_point(alpha = 0.9)
    }
    else {
      g <- g + geom_path(alpha = 0.9)
    }


    g + scale_size_continuous(guide = FALSE) +
      scale_color_discrete(guide = FALSE) +
      theme_bw() +
      theme(legend.position = "bottom",panel.border = element_rect(size = 0,colour = "white"),
            axis.line.x = element_line(color="black", size = 0.5))
  }, height = function() {
    ifelse(groupPg$width == 0,
           0,
           ifelse(input$selectLocus == "ALL", 1, max(hapPg$width*20, 400)))
  })

  MakeHapFreqTable <- reactive({
    if (is.null(input$selectLocus) ||
        input$selectLocus == "ALL" ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l))
      return ()

    if (is.null(haplo.summaryTbl())) {
      return()
    }

    n.indiv <- nrow(haplo.summaryTbl())
    allelic.freq.tbl <- haplo.summaryTbl() %>%
      ungroup() %>%
      gather(., whichHap,  hap, 4:5) %>%
      group_by(locus, hap) %>%
      summarise(f= n()/(2*n.indiv),
                n.occur = n(),
                tot.rd =  sum((whichHap=="haplotype.1")*read.depth.1)+
                  sum((whichHap=="haplotype.2")*read.depth.2))

    if (nrow(allelic.freq.tbl) == 0)
      return()

    all.hap <- unique(allelic.freq.tbl$hap)
    allelic.freq.tbl <- allelic.freq.tbl %>%
      mutate(hap.label.w.num = paste0(hap,
                                      " (",
                                      as.numeric(factor(hap, levels=all.hap)),
                                      ")"))

  })

  # event when haplotype frequency plot is being interacted or clicked
  output$hapFreqClicked <- renderUI({

    default.msg <- "Haplotype frequency plot"
    allelic.freq.tbl <- MakeHapFreqTable()
    if (is.null(allelic.freq.tbl)) {
      return()
    }

    nearpoint <- nearPoints(allelic.freq.tbl, input$hapFreqPlotClick, xvar="f", yvar="hap.label.w.num",
                            threshold = 10,
                            maxpoints = 1)

    if((!is.null(nearpoint)) && nrow(nearpoint)>0) {
      return(HTML(paste0(
        "# of observed allele copies: <b>",nearpoint$n.occur,
        "</b>\t,\t",
        " # of read(s): <b>",nearpoint$tot.rd,"</b>")))
    }
    else if (!is.null(input$hapFreqPlotHover)) {
      return(HTML("Click on any haplotype point for stat"))
    }
    else{
      return(default.msg)
    }

  })


  output$hapFreq <- renderPlot({

    allelic.freq.tbl <- MakeHapFreqTable()
    if (is.null(allelic.freq.tbl)) {
      return()
    }

    ggplot(allelic.freq.tbl, aes(
      y = hap.label.w.num,
      x = f,
      color = factor(hap.label.w.num)
    )) +
      geom_point(size = 4) +
      scale_color_discrete(guide = FALSE) +
      #geom_text(data=haplo.tot.read.tbl, aes(y=hap1, x=1, label=paste0(tot.depth," reads")))+
      xlab("observed freq") +
      ylab("haplotype") +
      scale_x_continuous(limits=c(0,1.1), breaks=c(0,0.3,0.6,0.9))+
      #xlim(c(0,1.1))+
      theme_bw()+
      theme(plot.margin = unit(c(2, 0, 2, 0), "mm"),
            panel.border = element_rect(size = 0,colour = "white"))
  }, height = function() {
    ifelse(groupPg$width == 0,
           0,
           ifelse(input$selectLocus == "ALL", 1, max(hapPg$width*30, 400)))
  })


  # event when h-w plot is being interacted
  output$hwClicked <- renderUI({

    if (is.null(input$selectLocus) ||
        input$selectLocus == "ALL" ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l))
      return ()

    if (is.null(hapPg$HW.tbl)) return ()

    default.msg <- "Hardy-Weinberg plot"

    nearpoint <- nearPoints(hapPg$HW.tbl, input$HWplotClick, xvar="hap1", yvar="hap2", threshold = 10,
                            maxpoints = 1)

    if((!is.null(nearpoint)) && nrow(nearpoint)>0) {
      #hapPg$HW.exp <- nearpoint$size
      #hapPg$HW.obs <- nearpoint$n.x
      HTML(paste0("Individuals counts  - ",
                  " Observed # : <b>", ceiling(nearpoint$n.x),"</b>",
                  " , Expected # : <b>",ceiling(nearpoint$size),"</b>"))

    }
    else if (!is.null(input$HWPlotHover)) {
      return("Click on any microhaplotype pair for stat")
    }
    else{
      return(default.msg)
    }

  })



  output$PairWiseHap <- renderPlot({
    if (is.null(input$selectLocus) ||
        input$selectLocus == "ALL" ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l))
      return ()

    if (is.null(haplo.summaryTbl())) {
      return()
    }

    haplo.filter <- haplo.summaryTbl()


    haplo.filter <- haplo.filter %>%
      group_by(locus, id) %>%
      summarise(
        hap1 = sort(c(haplotype.1, haplotype.2))[1],
        hap2 = sort(c(haplotype.1, haplotype.2))[2]
      ) %>%
      ungroup() %>%
      group_by(locus, hap1, hap2) %>%
      summarise(n = n())

    n.hap <- 2 * sum(haplo.filter$n)
    freq.hap <- gather(haplo.filter, whichHap,  hap, 2:3) %>%
      group_by(locus, hap) %>%
      summarise(n = sum(n) / n.hap) %>%
      mutate(hap1 = hap, n1 = n) %>%
      expand(., nesting(hap, n), nesting(hap1, n1)) %>%
      mutate(freq = ifelse(hap == hap1, n * n1, 2 * n * n1)) %>%
      rename("hap1" = hap, "hap2" = hap1) %>%
      group_by(locus, hap1, hap2) %>%
      mutate(re.hap1 = sort(c(hap1, hap2))[1],
             re.hap2 = sort(c(hap1, hap2))[2])


    all.hap <- unique(c(freq.hap$re.hap1,freq.hap$re.hap2))
    haplo.filter <- haplo.filter %>% ungroup() %>% mutate(hap1 = factor(as.numeric(factor(hap1, levels=all.hap)),
                                                                        levels=1:length(all.hap)),
                                                          hap2 = factor(as.numeric(factor(hap2, levels=all.hap)),
                                                                        levels=1:length(all.hap)))

    freq.hap <- freq.hap %>% ungroup() %>% mutate(hap1 = factor(as.numeric(factor(re.hap1, levels=all.hap)),
                                                                levels=1:length(all.hap)),
                                                  hap2 = factor(as.numeric(factor(re.hap2, levels=all.hap)),
                                                                levels=1:length(all.hap)),
                                                  size = freq * n.hap / 2)


    hapPg$HW.tbl <- left_join(haplo.filter, freq.hap ,by=c("hap1", "hap2"))

    ggplot() +
      geom_point(data=haplo.filter,
                 aes(
                   x = hap1,
                   y = hap2,
                   size = n,
                   color = hap1 == hap2
                 )) +
      xlab("H-W plot: haplotype pair (x,y)") +
      ylab("") +
      geom_point(
        data = freq.hap,
        aes(
          x = hap1,
          y = hap2,
          size = freq * n.hap / 2
        ),
        shape = 21,
        fill = NA,
        color = "black"
      ) +
      scale_color_discrete(guide = FALSE) +
      scale_size_continuous(range = c(3, 20), guide = FALSE) +
      scale_x_discrete(limits=1:length(all.hap))+
      scale_y_discrete(limits=1:length(all.hap))+
      theme_bw()+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.spacing = unit(c(0,0,0,0), 'mm'),
            panel.border = element_rect(size = 0,colour = "white"),
            plot.margin = unit(c(2, 0, 2, 0), "mm"))
  }, height = function() {
    ifelse(groupPg$width == 0,
           0,
           ifelse(input$selectLocus == "ALL", 1, max(hapPg$width*30, 400)))
  })

  # event when haplotype by group
  output$hapByGroupPlotClicked <- renderUI({

    if (is.null(input$selectLocus) ||
        input$selectLocus == "ALL" ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l))
      return ()

    if (is.null(haplo.summaryTbl())) return ()

    n.indiv <- nrow(haplo.summaryTbl())
    allelic.freq.tbl <-  haplo.summaryTbl() %>%
      ungroup() %>%
      gather(., whichHap,  hap, 4:5) %>%
      group_by(group, locus, hap) %>%
      summarise(f= n()/(2*n.indiv),
                n.indiv = length(unique(id)))

    default.msg <- "Variance of haplotype by group"

    nearpoint <- nearPoints(allelic.freq.tbl, input$hapByGroupPlotClick, xvar="group", yvar="hap",
                            threshold = 10,
                            maxpoints = 1)

    if((!is.null(nearpoint)) && nrow(nearpoint)>0) {
      return(HTML(paste0("# of individuals carry at least one allele copy: <b>",nearpoint$n.indiv,"</b>")))
    }
    else if (!is.null(input$hapByGroupPlotHover)) {
      return("Click on any point for counts")
    }
    else{
      return(default.msg)
    }

  })

  output$hapByGroupPlot <- renderPlot({
    if (is.null(input$selectLocus) ||
        input$selectLocus == "ALL" ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l))
      return ()

    if (is.null(haplo.summaryTbl())) {
      return()
    }

    n.indiv <- nrow(haplo.summaryTbl())
    allelic.freq.tbl <-  haplo.summaryTbl() %>%
      ungroup() %>%
      gather(., whichHap,  hap, 4:5) %>%
      group_by(group, locus, hap) %>%
      summarise(f= n()/(2*n.indiv),
                n.indiv = length(unique(id)))

    ggplot(allelic.freq.tbl, aes(
      x = group,
      y = hap,
      color = hap,
      size = f
    )) +
      geom_point() +
      # geom_text(aes(label=n),
      #           colour="black", vjust = "bottom", hjust="left", size=3.5)+
      xlab("") +
      ylab("") +
      theme_bw() +
      scale_color_discrete(guide = FALSE) +
      scale_size_continuous(guide = FALSE)+
      theme(panel.border = element_rect(size = 0,colour = "white"))
    #axis.line.x = element_line(color="black", size = 0.5))

  }, height = function() {
    ifelse(groupPg$width == 0,
           0,
           ifelse(input$selectLocus == "ALL", 1, max(hapPg$width*20, 400)))
  })


  ## filter panel


  output$RDnARplot <- renderPlot({
    if (is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l) ||
        #input$selectLocus == "ALL" ||
        length(filterParam$minRD) ==0)
      return ()

    rd.lower.bound <- max(0, filterParam$minRD -1)
    rd.upper.bound <- filterParam$minRD+1
    ar.lower.bound <- max(0, filterParam$minAR - 0.05)
    ar.upper.bound <- min(1, filterParam$minAR +0.05)


    readDepthRange <- c(1, 5, 10, 30, 100, filterParam$minRD, rd.lower.bound, rd.upper.bound)
    allelicFracRange <- c(0, 0.1, 0.2, 0.5, 0.8,filterParam$minAR, ar.lower.bound, ar.upper.bound)

    if (input$selectLocus == "ALL") {
      readDepthRange <- c(1,filterParam$minRD,10,20,50)
      allelicFracRange <- c(0.1, filterParam$minAR, 0.5,0.8)
    }
    rd.af.grid <- expand.grid(unique(allelicFracRange), unique(readDepthRange))

    hap.sel <- Min.filter.haplo()

    if(is.null(hap.sel)) return()

    hap.top.n <- hap.sel %>%
      filter(rank <= filterParam$n.alleles) %>%
      group_by(locus, id) %>%
      mutate(tot.rd = sum(depth)) %>%
      ungroup()

    label.grid <- apply(rd.af.grid, 1, function(i){
      af<-i[1]
      rd<-i[2]

      hap.p <- hap.top.n %>%
        filter(allele.balance >= af,tot.rd >= rd)

      hap.n <- hap.p %>%
        group_by(locus) %>%
        summarise(num.hap = length(unique(haplo)), num.pass.indiv = length(unique(id))) %>%
        ungroup() %>%
        summarise(num.hap.pass = sum(num.hap),
                  num.indiv.pass = round(mean(num.pass.indiv,na.rm = TRUE)))

      hap.n %>%
        summarise(descr = paste0(num.hap.pass, " hap,\n", num.indiv.pass, " indiv"))
      # hap.all <- hap.sel %>% filter(allele.balance >= af) %>%
      #   group_by(locus, id) %>%
      #   mutate(max.depth = max(tot.depth)) %>% group_by(locus) %>%
      #   summarise(num.hap = length(unique(haplo)), num.ambig.indiv = sum(rank>filterParam$n.alleles)) %>%
      #   ungroup() %>%
      #   summarise(num.hap.pass.n = sum(num.hap),
      #             num.ambig.ind = round(mean(num.ambig.indiv,na.rm = TRUE)))
      #
      # cbind(hap.all, hap.n) %>%
      #   summarise(descr = ifelse(num.hap.pass == num.hap.pass.n,
      #                            paste0(num.hap.pass, " hap,\n",
      #                                   num.indiv.pass, " (",num.ambig.ind,") indiv"),
      #                            paste0(num.hap.pass," (",num.hap.pass.n,") hap,\n",
      #                                   num.indiv.pass, " (",num.ambig.ind,") indiv")
      #                            ))
    }) %>% bind_rows()

    rd.af.grid.tbl <- cbind(rd.af.grid, label.grid)
    colnames(rd.af.grid.tbl) <- c("af", "rd", "content")

    rd.af.grid.tbl <- rd.af.grid.tbl %>%
      group_by(af,rd)%>%
      mutate(color.grp = 1*(filterParam$minRD==rd) + 1*(filterParam$minAR==af))

    ggplot(rd.af.grid.tbl, aes(x=factor(rd), y=factor(af))) +
      geom_tile(color="black", aes(fill=factor(color.grp)), size=0.1, linetype="dashed")+
      geom_text(aes(label=content))+
      theme_bw()+
      theme(legend.position = "none",
            axis.ticks = element_line(size = 0))+
      scale_x_discrete(expand = c(0, 0))+
      scale_y_discrete(expand = c(0, 0))+
      xlab("min. total read depth")+
      ylab("min. allelic ratio")+
      scale_fill_manual(values=c("#FDE1E1","#FFAFAD","#FF8582"))



  }, height = function() {
    ifelse(input$selectLocus == "ALL",320, 400)})




  output$downloadData <- downloadHandler(
    filename = function() {
      if (input$selectTbl == "reported indiv haplotype (diploid)")
        return ("reported_diploid_haplotype.csv")
      if (input$selectTbl == "observed variants (filtered)")
        return ("observed_filtered_haplotype.csv")
      if (input$selectTbl == "observed variants (unfiltered)")
        return ("observed_unfiltered_haplotype.csv")
      if (input$selectTbl == "reported haplotype (flat)")
        return ("reported_full_haplotype.csv")
      if (input$selectTbl == "SNP report")
        return ("snp_report.csv")
      if (input$selectTbl == "locus annotation")
        return ("locus_annotation.csv")
    }
    ,
    content = function(file) {
      if (isolate(input$selectTbl) == "reported indiv haplotype (diploid)") {
        if (is.null(haplo.summaryTbl()))
          return()

        haplo.all <- haplo.summaryTbl() %>% rename("indiv.ID" = id) %>%
          mutate(ar = round(ar,3))

        write.csv(haplo.all, file)
      }
      if (isolate(input$selectTbl) ==  "observed variants (unfiltered)") {
        if (is.null(Min.filter.haplo())) return()
        haplo.all <- Min.filter.haplo() %>% rename("indiv.ID" = id)
        write.csv(haplo.all, file)
      }
      if (isolate(input$selectTbl) ==  "observed variants (filtered)") {
        if (is.null(Filter.haplo.sum())) return()
        haplo.all <- Filter.haplo.sum() %>% rename("indiv.ID" = id)
        write.csv(haplo.all, file)
      }

      if (isolate(input$selectTbl) == "SNP report") {
        haplo.summaryTable <- haplo.summaryTbl()
        n.base <- nchar(haplo.summaryTable$haplotype.1)
        haplo.all <-
          haplo.summaryTable[rep(seq(1, nrow(haplo.summaryTable)), n.base),] %>%
          group_by(locus, id, group) %>%
          mutate(snp.id = row_number(),
                 snp = paste0(
                   substr(haplotype.1, snp.id, snp.id),
                   "/",
                   substr(haplotype.2, snp.id, snp.id)
                 )) %>%
          select(-haplotype.1,-haplotype.2,-read.depth.1,-read.depth.2) %>%
          rename("indiv.ID" = id)

        write.csv(haplo.all, file)

      }
      if (isolate(input$selectTbl) ==  "locus annotation") {
        #        if (is.null(Filter.haplo.sum())) return()
        write.csv(annotateTab$tbl, file)
      }
      if (isolate(input$selectTbl) ==  "reported haplotype (flat)"){

        haplo.sum <- Filter.haplo.sum()
        if (is.null(haplo.sum))
          return ()

        call.indx <- which(haplo.sum$allele.balance >= haplo.sum$min.ar.hz)

        ar <- as.numeric((haplo.sum$rank>1)*haplo.sum$allele.balance)

        toss.indiv.indx <- which( (ar >= haplo.sum$min.ar)*(ar < haplo.sum$min.ar.hz) +
                                    (ar <= haplo.sum$min.ar)*(ar > haplo.sum$max.ar.hm) >0)

        toss.indiv <- haplo.sum[toss.indiv.indx,] %>% select(id, locus) %>% unique()

        haplo.all <- anti_join(haplo.sum[call.indx, ], toss.indiv, by=c("id", "locus")) %>%
          ungroup() %>%
          select(group, locus, id, haplo, allele.balance, depth) %>%
          mutate(allele.balance = round(allele.balance, 3)) %>%
          rename("indiv.ID" = id)

        write.csv(haplo.all, file)

      }

    }
  )

  output$haploTbl <- DT::renderDataTable({
    #observeEvent(input$updateTable, {
    #observeEvent(input$selectTbl, {

    if (input$selectTbl == "reported indiv haplotype (diploid)") {

      haplo.summary <- haplo.summaryTbl()
      if (is.null(haplo.summary))
        return()

      haplo.all <- haplo.summary %>% rename("indiv.ID" = id) %>%
        mutate(ar = round(ar,3))
    }

    if (input$selectTbl ==  "observed variants (unfiltered)") {
      haplo.filter <- Min.filter.haplo()
      if (is.null(haplo.filter)) return()

      haplo.all <- haplo.filter %>% rename("indiv.ID" = id) #%>%
      #select(-sum.Phred.C, -max.Phred.C)
    }


    if (input$selectTbl == "observed variants (filtered)") {
      if (is.null(Filter.haplo.sum())) return()

      haplo.all <-
        Filter.haplo.sum() %>% rename("indiv.ID" = id) #%>%
      #select(-sum.Phred.C, -max.Phred.C)

      #if ("mapq" %in% colnames(haplo.all)) haplo.all <- haplo.all %>% select(-mapq)

    }

    if (input$selectTbl == "SNP report") {
      if (is.null(haplo.summaryTbl())) return()


      haplo.summaryTable <- haplo.summaryTbl()
      n.base <- nchar(haplo.summaryTable$haplotype.1)
      haplo.all <-
        haplo.summaryTable[rep(seq(1, nrow(haplo.summaryTable)), n.base),] %>%
        group_by(locus, id, group) %>%
        mutate(snp.id = row_number(),
               snp = paste0(
                 substr(haplotype.1, snp.id, snp.id),
                 "/",
                 substr(haplotype.2, snp.id, snp.id)
               )) %>%
        select(-haplotype.1,-haplotype.2,-read.depth.1,-read.depth.2) %>%
        rename("indiv.ID" = id)
    }

    if (isolate(input$selectTbl) ==  "locus annotation") {
      haplo.all <- annotateTab$tbl
    }

    if (isolate(input$selectTbl) ==  "reported haplotype (flat)"){
      haplo.sum <- Filter.haplo.sum()
      if (is.null(haplo.sum))
        return ()

      call.indx <- which(haplo.sum$allele.balance >= haplo.sum$min.ar.hz)

      ar <- as.numeric((haplo.sum$rank>1)*haplo.sum$allele.balance)

      toss.indiv.indx <- which( (ar >= haplo.sum$min.ar)*(ar < haplo.sum$min.ar.hz) +
                                  (ar <= haplo.sum$min.ar)*(ar > haplo.sum$max.ar.hm) >0)

      toss.indiv <- haplo.sum[toss.indiv.indx,] %>% select(id, locus) %>% unique()

      haplo.all <- anti_join(haplo.sum[call.indx, ], toss.indiv, by=c("id", "locus")) %>%
        ungroup() %>%
        select(group, locus, id, haplo, allele.balance, depth) %>%
        mutate(allele.balance = round(allele.balance, 3)) %>%
        rename("indiv.ID" = id)

    }


    #output$haploTbl <- DT::renderDataTable({
    DT::datatable(haplo.all, options = list(lengthMenu = list(c(5, 15, -1), c(
      '5', '15', 'All'
    )),
    pageLength = 15))
    #})

  })

  #Run Senor Microhap

  Run.SrMicrohap <- reactive({
    if (!srhapPg$makePlot)
      return ()

    haplo.sum <- update.Haplo.file()
    if (is.null(haplo.sum))
      return ()

    #    message( input$gibbIter, "\t", input$randomSeed,"we got stuff----\n")

    progress <- shiny::Progress$new()
    progress$set(message = 'Initiate sampling', value = 0.1)
    on.exit(progress$close())
    RunSrMicrohap(
      haplo.sum,
      srhapPg$locus.select,
      srhapPg$num.iter,
      srhapPg$random.seed,
      srhapPg$prior.model
      #srhapPg$min.read.depth
    )

  })

  output$HapFreqByGroupPlot <- renderPlot({
    run.result <- Run.SrMicrohap()
    if (is.null(run.result))
      return()
    #if (run.result$n.haplo.pair==1) return()

    start.iter <-
      min(floor(run.result$n.sam * (srhapPg$frac.burn / 100)),
          run.result$n.sam - 2)


    freq.matrix <- array(
      unlist(run.result$save.pfreq),
      dim = c(run.result$n.group,
              run.result$n.haplo,
              run.result$n.sam)
    )[, , (start.iter:run.result$n.sam)]

    #message( dim(freq.matrix),"<- freq matrix----\n")
    if(run.result$n.group > 1) {
      hap.freq.stat <-
        data.frame(t(apply(expand.grid(1:(run.result$n.group), 1:(run.result$n.haplo)),
                           1,
                           function(i)
                             c(
                               group = run.result$group[i[1]],
                               hap = paste0(run.result$haplo[i[2]],
                                            " (",
                                            i[2],
                                            ")",
                                            collapse = ""),
                               mean = mean(freq.matrix[i[1], i[2], ]),
                               st.CI = quantile(freq.matrix[i[1], i[2], ], prob =
                                                  0.05),
                               end.CI = quantile(freq.matrix[i[1], i[2], ], prob =
                                                   0.95)
                             ))), stringsAsFactors = FALSE)
    }
    else{
      hap.freq.stat <-
        data.frame(t(apply(expand.grid(1:(run.result$n.group), 1:(run.result$n.haplo)),
                           1,
                           function(i)
                             c(
                               group = run.result$group[i[1]],
                               hap = paste0(run.result$haplo[i[2]],
                                            " (",
                                            i[2],
                                            ")",
                                            collapse = ""),
                               mean = mean(freq.matrix[i[2], ]),
                               st.CI = quantile(freq.matrix[i[2], ], prob =
                                                  0.05),
                               end.CI = quantile(freq.matrix[i[2], ], prob =
                                                   0.95)
                             ))), stringsAsFactors = FALSE)
    }

    ggplot(hap.freq.stat,
           aes(
             y = group,
             yend = group,
             x = as.numeric(st.CI.5.),
             xend = as.numeric(end.CI.95.)
           )) +
      geom_segment() +
      facet_grid(hap ~ ., space = "free") + #space = "free" #scales="free"
      theme(strip.text.y = element_text(angle = 0)) +
      geom_point(data = hap.freq.stat,
                 aes(x = as.numeric(mean), y = group),
                 pch = 3) +
      ylab("") +
      xlab("overall haplotype frequency (90% CI)")

  }, height = function() {
    ifelse(srhapPg$makePlot, 300, 0)
  })


  output$indivHapPosPlot <- renderPlot({
    run.result <- Run.SrMicrohap()

    progress <- shiny::Progress$new()
    progress$set(message = 'Generating figures', value = 0.5)
    on.exit(progress$close())

    if (is.null(run.result))
      return()
    if (run.result$n.haplo.pair == 1)
      return()


    start.iter <-
      min(floor(run.result$n.sam * (srhapPg$frac.burn / 100)),
          run.result$n.sam - 2)
    freq.matrix <- array(
      unlist(run.result$save.hap),
      dim = c(run.result$n.indiv,
              run.result$n.haplo,
              run.result$n.sam)
    )[, , (start.iter:run.result$n.sam)]
    length.iter <- run.result$n.sam - start.iter + 1

    ref.hap.matrix <-
      matrix(0, nrow = run.result$n.haplo.pair, ncol = run.result$n.haplo)
    ref.hap.matrix[cbind(1:run.result$n.haplo.pair, run.result$haplo.pair[, 1])] <-
      1
    ref.hap.matrix[cbind(1:run.result$n.haplo.pair, run.result$haplo.pair[, 2])] <-
      1 + ref.hap.matrix[cbind(1:run.result$n.haplo.pair, run.result$haplo.pair[, 2])]

    ref.hap.str <-
      apply(ref.hap.matrix, 1, function(x)
        paste(x, collapse = ""))

    ref.hap.label <-
      as.vector(apply(run.result$haplo.pair, 1, function(x)
        paste(x, collapse = "/")))
    names(ref.hap.label) <- ref.hap.str


    sample.hap.df <- lapply(1:run.result$n.indiv, function(i) {
      data.frame(indiv = i,
                 table(apply(freq.matrix[i, , ],
                             2,
                             function(x)
                               ref.hap.label[paste(x, collapse = "")])),
                 stringsAsFactors = FALSE)
    }) %>% bind_rows()

    sample.hap.df <- sample.hap.df %>%
      mutate(
        posterior = Freq / length.iter,
        indiv.id = run.result$indiv.id[indiv],
        group = run.result$group[run.result$grp.assoc.indiv[indiv]]
      )

    # gather empiricial result for comparison
    haplo.filter <- srhapPg$data.table %>%
      filter(depth >=filterParam$minRD,
             rank <= 2,
             allele.balance >= filterParam$minAR) %>%
      group_by(group, id) %>%
      arrange(-depth) %>%
      summarise(
        haplotype.pair = ifelse(length(depth) == 1,
                                paste0(sum(which(run.result$haplo==haplo[1])),
                                       "/",
                                       sum(which(run.result$haplo==haplo[1]))),
                                paste0(sort(c(sum(which(run.result$haplo==haplo[1])),
                                              sum(which(run.result$haplo==haplo[2])))),
                                       collapse = "/")))

    ggplot(sample.hap.df,
           aes(x = Var1, y = indiv.id, fill = posterior)) +
      geom_tile() +
      geom_tile(data=haplo.filter, aes(y=id, x=haplotype.pair), fill=NA, color="purple",lwd=1)+
      facet_grid(group ~ ., space = "free", scales = "free_y") +
      geom_text(aes(label = round(posterior, 2))) +
      theme_bw() +
      scale_fill_gradient(low = "white",
                          high = "light green",
                          guide = FALSE) +
      ylab("") +
      xlab("haplotype pair")

  }, height = function() {
    ifelse(srhapPg$makePlot,
           max(15 * panelParam$tot.indiv,
               450),
           0)
  })

})
