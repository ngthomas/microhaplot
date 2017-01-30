library(shiny)
library(shinyBS)
library("ggplot2")
library("plyr")
library("dplyr")
library("tidyr")
library("DT")
library("grid")
library("scales")
library("haplot")
library("reshape2")


shinyServer(function(input, output, session) {


  output$about <- renderUI({
    HTML(
      paste("<strong>Microhaplot</strong> offers a streamline visual environment to assess quality of microhaplotype extracted from short read alignment files.
          You can find most of the interactive features such as defining critera or adding locus comments at the top panel
while the bottom panel hosts a wide selection of tables and graphical summaries.",
            "<br/>",
            "<i>Summary Info</i>: microhaplotype summaries that are either grouped by group label, individual, or locus. These plots are great
          for gathering the big picture",
            "<i>Filter Analysis</i>: this section is useful to fine-tune your criteria at a single locus level",
            "<i>Inferential Analysis</i>: still in the works",
            "<i>Output</i>: You can view or download tables of raw or finalized microhaplotypes (in csv format)",
            "<br/>",
            "<b>Contact & Citation</b>",
            sep="<br/>"))
  })

  addTooltip(session, "downloadData",
             "Note: The download table is sensitive to the 'field selection' input")

  dirFiles <- list.files()
  rds.file <- grep(".rds", dirFiles)

  pos.rds.file <- grep("_posinfo.rds", dirFiles)
  annotate.rds.file <- grep("_annotate.rds", dirFiles)
  rds.file <- setdiff(rds.file, c(pos.rds.file,annotate.rds.file))
  haplo.sum <- NULL

  for (rdsFile in rds.file){
    colnames.file <- colnames(rdsFile)
    "group" %in% colnames.file

  }


  if (length(rds.file) > 0) {
    select.file.tem <- dirFiles[rds.file[1]]
    updateSelectInput(session,
                      "selectDB",
                      selected = select.file.tem,
                      choices = dirFiles[rds.file])
  }

  update.Haplo.file <- reactive({
    if (input$selectDB == "" ||
        is.null(input$selectDB) || !file.exists(input$selectDB))
      return()
    #cat(file=stderr(), "select DB_", input$selectDB, "_----\n")
    readRDS(input$selectDB)  %>% ungroup() %>% mutate(id = as.character(id))
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
      annotateTab$tbl <- data.frame(locus = panelParam$locus.label,
                                    ave.entropy = c("NA",rep(0, panelParam$n.locus)),
                                    min.rd = rep(0, panelParam$n.locus+1),
                                    min.ar = rep(0, panelParam$n.locus+1),
                                    status = c("NA",rep("Accept", panelParam$n.locus)),
                                    comment = rep("", panelParam$n.locus+1),
                                    stringsAsFactors = F)


      return()
    }


    annotateTab$tbl <- readRDS(annotate.file)
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
    updateNumericInput(session,
                       "minRD",
                       value=annotateTab$tbl$min.rd[match.indx])
    updateSliderInput(session,
                      "minAR",
                      value=annotateTab$tbl$min.ar[match.indx])

    filterParam$minRD <- annotateTab$tbl$min.rd[match.indx]#input$coverageMin
    filterParam$minAR <- annotateTab$tbl$min.ar[match.indx]#input$minAlleleRatio

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
  ranges <- reactiveValues(y = NULL, x = NULL)
  rangesH <- reactiveValues(y = NULL)

  locusPg <- reactiveValues(l = NULL, width = NULL)
  indivPg <- reactiveValues(i = NULL, width = NULL)
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

  filterParam <- reactiveValues(minRD = 0, minAR = 0, hover.minAR = 0, hover.minRD =0)
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
    cat(file = stderr(), "select DB_", input$selectDB, "_", "----\n")
    if (input$selectDB == "" || is.null(input$selectDB))
      return()

    haplo.sum <- update.Haplo.file()

    panelParam$tot.indiv <- length(unique(haplo.sum$id))

    if (!"group" %in% colnames(haplo.sum))
      haplo.sum <-
      cbind.data.frame("group" = "unlabel",
                       haplo.sum,
                       stringsAsFactors = F) %>% tbl_df
    cat(file = stderr(),
        "preview_",
        head(haplo.sum, 1) %>% unlist(),
        "_----\n")
    panelParam$n.locus <- length(unique(haplo.sum$locus))
    panelParam$n.indiv <- length(unique(haplo.sum$id))

    locus.sorted <- sort(unique(haplo.sum$locus))
    panelParam$locus.label.tbl <-
      data.frame(locus = locus.sorted, stringsAsFactors = F) %>% tbl_df()
    panelParam$locus.label <- c("ALL", locus.sorted)
    panelParam$locus.label.bare <- locus.sorted

    indiv.sorted <- sort(unique(haplo.sum$id))
    panelParam$indiv.label.tbl <-
      data.frame(id = indiv.sorted, stringsAsFactors = F) %>% tbl_df()
    panelParam$indiv.label <- c("ALL", indiv.sorted)
    panelParam$indiv.label.bare <- indiv.sorted

    group.sorted <- sort(unique(haplo.sum$group))
    panelParam$group.label.tbl <-
      data.frame(id = group.sorted, stringsAsFactors = F) %>% tbl_df()
    panelParam$group.label <- c("ALL", group.sorted)
    panelParam$group.label.bare <- group.sorted
    panelParam$n.group <- length(group.sorted)

    updateSelectInput(session,
                      "selectLocus",
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

    extract.annotate.file()
    update.annotate.field()

    Min.filter.haplo()
    #Filter.haplo.by.RDnAR()
    #Filter.haplo.sum()

    srhapPg$makePlot <- FALSE
    extract.pos.file()
    #cat(file=stderr(), "--", colnames(annotateTab$tbl) %>% unlist, "\n")

  }, priority = -3)

  ## updating Locus and individidual choice at the start of the session:
  #updateSelectInput(session, "selectLocus", selected="ALL", choices=locus.label)
  #updateSelectInput(session, "selectIndiv", selected="ALL", choices=indiv.label)

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

  # observeEvent(input$groupBack, {
  #   indx <- isolate(which(panelParam$group.label == input$selectGroup))
  #   label <-
  #     ifelse(indx > 1, panelParam$group.label[indx - 1], panelParam$group.label[indx])
  #   updateSelectInput(session, "selectGroup", selected = label)
  # })
  # observeEvent(input$groupFor, {
  #   indx <- isolate(which(panelParam$group.label == input$selectGroup))
  #   label <-
  #     ifelse(
  #       indx < length(panelParam$group.label),
  #       panelParam$group.label[indx + 1],
  #       panelParam$group.label[indx]
  #     )
  #   updateSelectInput(session, "selectGroup", selected = label)
  # })


  # reacting to the filter update button
  observeEvent(input$updateFilter, {
    if(is.na(input$coverageMin) || input$coverageMin <0) {
      createAlert(session, "alert", "filterAlert", title = "Invalid input",
                  content = "Minimum read coverage must be a positive integer", append = FALSE)
      return()
    }

    closeAlert(session, "filterAlert")

    if(input$coverageMin == filterParam$minRD && input$minAlleleRatio == filterParam$minAR)
      return ()

    filterParam$minRD <- input$coverageMin
    filterParam$minAR <- input$minAlleleRatio
    Filter.haplo.by.RDnAR()
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
      data.frame(id = indiv.sorted, stringsAsFactors = F) %>% tbl_df()

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
  }, priority = -2)



  observeEvent(input$selectLocus, {

    indx <-
      isolate(which(panelParam$locus.label.bare == input$selectLocus))

    output$locusSelect <- renderText({
      input$selectLocus
    })
    output$locusSelect1 <- renderText({
      input$selectLocus
    })

    update.annotate.field()

    if (input$selectLocus != "ALL") {
      output$maxlocusPage <- renderText({
        "1"
      })
      updateNumericInput(session, "locusPage", value = 1, max = 1)

      updateCheckboxGroupInput(session,
                               "filterOpts",
                               choices=list("keeps only top two haplotypes (per indiv)"=1))
      removePopover(session, "filterOpts")

      locusPg$l <- input$selectLocus
      rangesH$y <- c(0, 2)
      # output$locusAcceptStatus <-
      #   renderText({
      #     ifelse(panelParam$is.reject[indx] == 0, "Accept", "Reject")
      #   })

      closeAlert(session,"hapLocusAlert")
      closeAlert(session,"cuthapLocusAlert")
    }
    else {
      updateCheckboxGroupInput(session,
                               "filterOpts",
                               choices=list("keeps only top two haplotypes (per indiv)"=1,
                                            "relies only on locus-specific param."=2,
                                            "serves as the minimal baseline"=3))

      addPopover(session, "filterOpts","Options",
                 content=paste0("<p>By default, the observed microhaplotypes must meet the current selected filter criteria</p>",
                                "<p></p>",
                                "<p><b>relies only on locus-specific param.</b>: this option disregards the present selection of min. read depth and allelic ratio. ",
                                "Instead, the filtering process uses parameters defined at a single-locus level</p>",
                                "<p></p>",
                                "<p><b>serves as the minimal baseline</b>:all loci must pass the current selected filter values and
                                 locus-specific filter value</p> "),
                 placement="bottom",
                 trigger="hover")

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

      createAlert(session, "hapAlert", "hapLocusAlert", title = "No Locus selected",
                  content = "choose a locus to view content", append = FALSE)
      createAlert(session, "cutoffhapAlert", "cuthapLocusAlert", title = "choose single locus for more detail",
                  content = "choose a locus to view a more detailed breakdown of those filter criteria by microhaplotype", append = FALSE)
    }

    match.indx <- which(annotateTab$tbl$locus == input$selectLocus)
    updateNumericInput(session, "coverageMin", value=annotateTab$tbl$min.rd[match.indx])
    updateSliderInput(session, "minAlleleRatio", value=annotateTab$tbl$min.ar[match.indx])

    filterParam$minRD <- annotateTab$tbl$min.rd[match.indx]#input$coverageMin
    filterParam$minAR <- annotateTab$tbl$min.ar[match.indx]#input$minAlleleRatio

    #cat(file=stderr(), "min READ depth:", filterParam$minRD, "\n")
    #Min.filter.haplo()
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
      updateNumericInput(session, "locusPage", value = 1, max = 1)
      locusPg$l <- input$selectLocus
      rangesH$y <- c(0, 2)
    }
    else {
      if (input$locusPerDisplay == 100) {
        output$maxlocusPage <- renderText({
          "1"
        })
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
        #cat(file=stderr(), "haha_", as.numeric(panelParam$n.locus)/as.numeric(input$locusPerDisplay), "_----\n")
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
    if(input$locusPage <=0){
      updateNumericInput(session, "locusPage", value = 1)
      return()
    }

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
    output$indivSelect <- renderText({
      input$selectIndiv
    })
    if (input$selectIndiv != "ALL") {
      output$maxIndivPage <- renderText({
        "1"
      })
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
  })

  # observeEvent(input$acceptLocus, {
  #   if (input$selectLocus != "ALL") {
  #     indx <-
  #       isolate(which(panelParam$locus.label.bare == input$selectLocus))
  #     panelParam$is.reject[indx] <- 0
  #     output$locusAcceptStatus <-
  #       renderText({
  #         ifelse(panelParam$is.reject[indx] == 0, "Accept", "Reject")
  #       })
  #   }
  # })
  #
  # observeEvent(input$rejectLocus, {
  #   if (input$selectLocus != "ALL") {
  #     indx <-
  #       isolate(which(panelParam$locus.label.bare == input$selectLocus))
  #     panelParam$is.reject[indx] <- 1
  #     output$locusAcceptStatus <-
  #       renderText({
  #         ifelse(panelParam$is.reject[indx] == 0, "Accept", "Reject")
  #       })
  #   }
  # })

  observeEvent(input$indivPerDisplay, {
    if (is.null(panelParam$n.indiv))
      return()

    if (input$selectIndiv != "ALL") {
      output$maxIndivPage <- renderText({
        "1"
      })
      updateNumericInput(session, "indivPage", value = 1, max = 1)
      indivPg$i <- input$selectIndiv
      ranges$y <- c(0, 2)
    }
    else {
      if (input$indivPerDisplay == 100) {
        output$maxIndivPage <- renderText({
          "1"
        })
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
    if(input$indivPage <=0){
      updateNumericInput(session, "indivPage", value = 1)
      return()
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

    if(input$rewriteFilter[1]) {
      annotateTab$tbl$min.rd[match.indx] <-  filterParam$minRD
      annotateTab$tbl$min.ar[match.indx] <-  filterParam$minAR
    }

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
    filterParam$minRD <- input$coverageMin
    filterParam$minAR <- input$minAlleleRatio
    Filter.haplo.sum()

    match.indx <- which(annotateTab$tbl$locus == input$selectLocus)

    annotateTab$tbl$min.rd[match.indx] <-  filterParam$minRD
    annotateTab$tbl$min.ar[match.indx] <-  filterParam$minAR

    annotate.file <- strsplit(input$selectDB, split=".rds") %>% unlist %>% paste0(.,"_annotate.rds")
    saveRDS(annotateTab$tbl, annotate.file)

  })



  haplo.summaryTbl <- reactive({
    haplo.sum <- Filter.haplo.sum()
    if (is.null(haplo.sum))
      return ()
    haplo.filter <- haplo.sum %>% filter(rank <=2)

    haplo.filter <- haplo.filter %>%
      group_by(locus, id, group) %>%
      arrange(-depth) %>%
      summarise(
        haplotype.1 = ifelse(length(depth) == 1, haplo[1], sort(haplo)[1]),
        haplotype.2 = ifelse(length(depth) == 1, haplo[1], sort(haplo)[2]),
        read.depth.1 = depth[1],
        read.depth.2 = ifelse(length(depth) == 1, depth[1], depth[2])
      )
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

    if (is.null(haplo.sum))
      return ()

    haplo.join.ar <- left_join(haplo.sum,
                               annotateTab$tbl,
                               by="locus")

    if ("3" %in% input$filterOpts)
      haplo.join.ar <- haplo.join.ar %>% mutate(min.rd = ifelse(filterParam$minRD>min.rd,
                                                                filterParam$minRD,
                                                                min.rd),
                                                min.ar= ifelse(filterParam$minAR>min.ar,
                                                               filterParam$minAR,
                                                               min.ar))

    if (! "2" %in% input$filterOpts)
      haplo.join.ar <- haplo.join.ar %>%
      mutate(min.rd = filterParam$minRD,
             min.ar = filterParam$minAR)

    haplo.join.ar <- haplo.join.ar %>% filter(depth >= min.rd,
                                              allele.balance >= min.ar)

    if ("1" %in% input$filterOpts)
      haplo.join.ar <- haplo.join.ar %>% filter(rank <= 2)


    haplo.join.ar
  })


  Filter.haplo.sum <- reactive({

    haplo.filter <- Filter.haplo.by.RDnAR()

    if (is.null(haplo.filter))
      return ()

    Update.ave.entropy()
    groupPg$width <- dim(haplo.filter)[1]
    hapPg$width <- haplo.filter %>% filter(rank <= 2) %>% select(haplo) %>% unique %>% unlist %>% length

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

    #cat(file=stderr(), "is it updating_", unlist(haplo.tot.tbl[1,]), "_----\n")



    uniqH.perI.tbl <-
      right_join(haplo.tot.tbl, panelParam$locus.label.tbl, by = "locus")
    if (is.null(uniqH.perI.tbl))
      return()
    uniqH.perI.tbl[is.na(uniqH.perI.tbl)] <- 0

    #cat(file=stderr(), "is it moving :_ _----\n")


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
        panel.margin = unit(0, 'mm'),
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
           max(
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
        panel.margin = unit(0, 'mm'),
        panel.border = element_rect(size = 0,colour = "white"),
        plot.margin = unit(c(0, 2, 0, 0), "mm")
      ) +
      ylim(locusPg$l) +
      scale_x_continuous(limits=c(0,max.haplo+1),breaks=round(seq(1,max.haplo, length.out=4)))+
      coord_cartesian(ylim = rangesH$y)
    #scale_x_discrete(limits=c(-1, max(frac.calleable$n)+1)) #breaks= pretty_breaks()
  }, height = function()
    ifelse(groupPg$width == 0, 0, max(
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

    #cat(file=stderr(), "checking", frac.calleable %>% filter(abs(f-0.5)<0.1) %>% select(f, label.x) %>% unlist(), "_----\n")

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
        panel.margin = unit(0, 'mm'),
        panel.border = element_rect(size = 0,colour = "white"),
        plot.margin = unit(c(0, 3, 0, 0), "mm")
      ) +
      ylim(locusPg$l) +
      coord_cartesian(ylim = rangesH$y) +
      #xlim(c(0, 1))
      scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.2))

  }, height = function()
    ifelse(groupPg$width == 0, 0,
           max(
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
      ylab("read depth \n(per indiv)") +
      geom_violin() +
      geom_point(aes(x = locus, y = mean.depth),
                 cex = 3,
                 pch = 3) +
      #geom_point()+
      theme_bw() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.margin = unit(0, 'mm'),
        panel.border = element_rect(size = 0,colour = "white"),
        plot.margin = unit(c(0, 0, 0, 0), "mm")
      ) +
      scale_y_log10() +
      xlim(locusPg$l) +
      coord_flip(xlim = rangesH$y)
  }, height = function()
    ifelse(groupPg$width == 0, 0, max(
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

    haplo.filter <- Filter.haplo.sum()

    if (dim(haplo.filter)[1] == 0)
      return ()
    if (dim(panelParam$indiv.label.tbl)[1] == 0)
      return()

    haplo.filter <- haplo.filter %>%
      group_by(locus, id, group) %>%
      summarise(depth.ratio = ifelse(length(depth) == 1, 0, min(allele.balance)),
                depth.first = max(depth))

    haplo.filter <-
      right_join(haplo.filter, panelParam$indiv.label.tbl, by = "id")
    if (is.null(haplo.filter))
      return()
    haplo.filter[is.na(haplo.filter)] <- 0

    #if (input$selectIndiv != "ALL") {
    #  haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv)
    #}


    ggplot(data = haplo.filter, aes(
      x = depth.ratio,
      y = id,
      size = log(depth.first, 10),
      color = group
    )) +
      geom_point(alpha = 0.4) +
      scale_size_continuous(guide = FALSE) + #"Read Depth of the most common haplotype (log 10)")+
      scale_color_discrete(
        guide = FALSE,
        drop = T,
        limits = levels(panelParam$group.label.bare)
      ) +
      theme_bw() +
      ylab("individual ID") +
      xlab ("ratio of the 2nd : 1st common haplotype\n(allelic ratio)") +
      theme(
        legend.position = "bottom",
        panel.margin = unit(0, 'mm'),
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
           max(
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

    filter.haplo <- Filter.haplo.sum() %>%
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
        panel.margin = unit(c(0,0,0,0), 'mm'),
        panel.border = element_rect(size = 0,colour = "white"),
        plot.margin = unit(c(0, 3, 0, 0), "mm")
      ) +
      ylim(indivPg$i) +
      scale_x_continuous(limits=c(0,max.locus+1),breaks=round(seq(1,max.locus, length.out=4)))+
      coord_cartesian(ylim = ranges$y)
    #scale_x_discrete(limits=c(-1, max(frac.calleable$n)+1)) #breaks= pretty_breaks()
  }, height = function()
    ifelse(groupPg$width == 0, 0,
           max(
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
        panel.margin = unit(0, 'mm'),
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
           max(
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
  #           panel.margin = unit(0, 'mm'))+
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
      right_join(Filter.haplo.sum(), panelParam$indiv.label.tbl, by = "id")
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
        panel.margin = unit(0, 'mm'),
        panel.border = element_rect(size = 0,colour = "white"),
        plot.margin = unit(c(0, 2, 0, 0), "mm")
      ) +
      #plot.margin = unit(c(0, 0, 0, 0), "mm"))+
      scale_y_log10()+
      xlim(indivPg$i) +
      coord_flip(xlim = ranges$y)

  }, height = function()
    ifelse(groupPg$width == 0, 0,
           max(
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
      updateNumericInput(session, "minAlleleRatio",value=filterParam$minAR)
      Filter.haplo.sum()
    }
  })

  # observeEvent(input$RDplot_hover, {
  #   if (!is.null(input$RDplot_hover)) filterParam$hover.minRD = input$RDplot_hover$x
  # })

  observeEvent(input$RDplot_dblclick, {
    if (!is.null(input$RDplot_dblclick)){
      filterParam$minRD <- ceiling(input$RDplot_dblclick$x)
      updateNumericInput(session, "coverageMin",value=filterParam$minRD)
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
  #       panel.margin = unit(0, 'mm'),
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


    mean.f.tbl <-
      frac.calleable %>% group_by(group) %>% summarise(mean.f = mean(f, na.rm =
                                                                       T))
    #     if(is.null(frac.calleable)) return()
    #     frac.calleable[is.na(frac.calleable)]<- 0
    #
    #     if (input$selectLocus != "ALL") {
    #       frac.calleable <- frac.calleable %>% filter(locus == input$selectLocus)
    #     }

    ggplot(frac.calleable, aes(x = f, y = group, color = group)) +
      geom_point(alpha = 0.5) +
      geom_point(
        data = mean.f.tbl,
        aes(y = group, x = mean.f),
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
        panel.margin = unit(0, 'mm'),
        plot.margin = unit(c(0, 0, 0, 0), "mm")
      )
  }, height = function() {
    ifelse(
      groupPg$width == 0,
      0,
      ifelse(input$selectGroup == "ALL", 30 * panelParam$n.group, 100)
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
  #       panel.margin = unit(0, 'mm'),
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

    ggplot(frac.calleable, aes(x = f, y = group, color = group)) +
      geom_point(alpha = 0.5) +
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
        panel.margin = unit(0, 'mm'),
        plot.margin = unit(c(0, 0, 0, 0), "mm")
      )
  }, height = function() {
    ifelse(
      groupPg$width == 0,
      0,
      ifelse(input$selectGroup == "ALL", 30 * panelParam$n.group, 100)
    )
  })


  ## filter status:cutoff distribution panel

  output$allReadDepth <- renderPlot({

    if (is.null(input$selectLocus) ||
        is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l)) {
      return()
    }

    haplo.filter <- Min.filter.haplo()
    if(is.null(Min.filter.haplo)) return()

    ggplot(haplo.filter, aes(x=depth))+
      #geom_density()+
      geom_histogram(binwidth = 0.05, boundary = -0.025)+
      scale_x_log10("distrib. of read depth")+
      scale_y_continuous("", breaks=NULL)+
      theme_bw()+
      geom_vline(
        xintercept = filterParam$minRD,
        linetype = "dashed",
        color = "red"
      )+
      theme(strip.text.y = element_text(angle =360,size=0, margin=margin(0,0,0,0)),
            panel.margin = unit(0, 'mm'),
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

    haplo.filter <- Min.filter.haplo()
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
            panel.margin = unit(0, 'mm'),
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

    haplo.filter <- Min.filter.haplo()
    if(is.null(Min.filter.haplo)) return()

    haplo.rep <- haplo.filter %>% filter(rank <= 2,depth >= filterParam$minRD,
                                         allele.balance >= filterParam$minAR) %>%
      select(haplo) %>% unique %>% unlist()

    if (length(haplo.rep)==0) return()

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
      theme(strip.text.y = element_text(angle =360, margin=margin(0,0,0,0)),
            strip.background  = element_rect(fill="white",size = 0),
            panel.margin = unit(0, 'mm'),
            plot.margin = unit(c(2, 0, 6, 0), "mm"),
            aspect.ratio =1000)

  }, height = function() {
    ifelse(hapPg$width == 0,
           0,
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

    haplo.filter <- Min.filter.haplo()
    if(is.null(Min.filter.haplo)) return()

    haplo.rep <- haplo.filter %>% filter(rank <= 2,depth >= filterParam$minRD,
                                         allele.balance >= filterParam$minAR) %>%
      select(haplo) %>% unique %>% unlist()

    haplo.match.rep <- haplo.filter %>% filter(haplo %in% haplo.rep) %>%
      mutate(haplo=factor(haplo, level=sort(haplo.rep,decreasing=T))) %>%
      group_by(haplo) %>%
      mutate(n.accept.indiv = sum(depth >=filterParam$minRD),
             n.reject.indiv = sum(depth < filterParam$minRD)) %>%
      ungroup() %>%
      mutate(
             center.x.accept = mean(depth[depth>=filterParam$minRD]),
             center.x.reject = mean(depth[depth<filterParam$minRD])) %>%
      group_by(haplo, id) %>%
      mutate(is.pass = 1*(depth >=filterParam$minRD))


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
    #   theme(panel.margin = unit(0, 'mm'),
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
         y=10), hjust=-1, vjust=-1)+
      geom_text(data=haplo.match.rep, aes(
        label=n.reject.indiv,
        x=center.x.reject,
        y=10), hjust=-1, vjust =-1)

      #geom_density(adjust=0.5, color=NA)+
    g + facet_grid(haplo~.)+#, scales="free_y")+
      scale_x_log10("distrib. of read depth")+
      scale_fill_discrete(guide=FALSE,direction=-1)+
      scale_y_continuous("",breaks=NULL)+ #breaks=1
      theme_bw()+
      geom_vline(
        xintercept = filterParam$minRD,
        linetype = "dashed",
        color = "red"
      )+
      theme(strip.text.y = element_text(angle =360,size=0, margin=margin(0,0,0,0)),
            panel.margin = unit(0, 'mm'),
            panel.border = element_rect(size = 0,colour = "white"),
            #axis.line.y = element_line(color="grey", size = 1),
            plot.margin = unit(c(2, 0, 2, 1), "mm"))

  }, height = function() {
    ifelse(hapPg$width == 0,
           0,
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

    haplo.filter <- Min.filter.haplo()

    haplo.rep <- haplo.filter %>% filter(rank <= 2,depth >= filterParam$minRD,
                                         allele.balance >= filterParam$minAR) %>%
      select(haplo) %>% unique %>% unlist()

    haplo.match.rep <- haplo.filter %>% filter(haplo %in% haplo.rep) %>%
      mutate(haplo=factor(haplo, level=sort(haplo.rep,decreasing=T))) %>%
      group_by(haplo) %>%
      mutate(n.accept.indiv = sum(allele.balance >=filterParam$minAR),
             n.reject.indiv = sum(allele.balance < filterParam$minAR)) %>%
      ungroup() %>%
      mutate(
        center.x.accept = min(mean(allele.balance[allele.balance>=filterParam$minAR]),
                              0.9),
        center.x.reject = min(mean(allele.balance[allele.balance<filterParam$minAR]),0.9))%>%
      group_by(haplo, id) %>%
      mutate(is.pass = 1*(allele.balance >=filterParam$minAR))

    if (nrow(haplo.match.rep)==0) return()

    g<- ggplot(haplo.match.rep, aes(x=allele.balance, fill=haplo))+
      geom_histogram(binwidth = 0.05,boundary = -0.05)+
      geom_point(data=haplo.match.rep, aes(x=allele.balance, y=0), size=1.2, alpha=0.4)
      #geom_density(adjust=0.1, color=NA)+

    if (filterCriteriaPg$AR.detail.on)
      g <- g +
      geom_text(data=haplo.match.rep, aes(
        label=n.accept.indiv,
        x=center.x.accept,
        y=10), hjust=2, vjust=-1)+
      geom_text(data=haplo.match.rep, aes(
        label=n.reject.indiv,
        x=center.x.reject,
        y=10), hjust=0, vjust =-1)


      g + facet_grid(haplo~.)+#, scales="free_y")+
      scale_x_log10("distrib. of allelic ratio",breaks=c(0.01,0.1,0.2,0.5,1))+
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
            panel.margin = unit(0, 'mm'),
            plot.margin = unit(c(2, 0, 2, 0), "mm"))

  }, height = function() {
    ifelse(hapPg$width == 0,
           0,
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
        #cat(file=stderr(), "character split_", unlist(char.split), "_----\n")
        n.char <- length(char.split[[1]])
        sapply(1:n.char, function(j)
          c(i, j, char.split[[1]][j], haplo.profile.frac[i,]$frac))
      }) %>%
      matrix(., ncol = 4, byrow = T) %>%
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
           ifelse(input$selectLocus == "ALL", 0, max(hapPg$width*20, 400)))
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
    obs.freq.tbl <-  haplo.summaryTbl() %>%
      ungroup() %>%
      group_by(locus) %>%
      mutate(tot.haplo = n()) %>%
      group_by(locus, haplotype.1, haplotype.2) %>%
      summarise(obs.freq = n() / tot.haplo[1],
                n.occur = n(),
                is.homo = (haplotype.1[1] == haplotype.2[1]),
                tot.read.depth.1=sum(read.depth.1),
                tot.read.depth.2=sum(read.depth.2))

    allelic.freq.tbl <-
      gather(obs.freq.tbl, whichHap,  hap1, 2:3) %>%
      group_by(locus, hap1) %>%
      summarise(f = sum(obs.freq / 2),
                n = sum(ifelse(is.homo, n.occur/2, n.occur)),
                tot.read.depth = sum((whichHap=="haplotype.1")*tot.read.depth.1)+
                  sum((whichHap=="haplotype.2")*tot.read.depth.2)
      )

    if (nrow(allelic.freq.tbl) == 0)
      return()

    all.hap <- unique(allelic.freq.tbl$hap1)
    allelic.freq.tbl <- allelic.freq.tbl %>%
      mutate(hap.label.w.num = paste0(hap1,
                           " (",
                           as.numeric(factor(hap1, levels=all.hap)),
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
        "# of individual(s): <b>",nearpoint$n,
        "</b>\t,\t",
        " # of read(s): <b>",nearpoint$tot.read.depth,"</b>")))
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
           ifelse(input$selectLocus == "ALL", 0, max(hapPg$width*30, 400)))
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

    haplo.filter <- Filter.haplo.sum()

    if (is.null(haplo.filter))
      return()


    haplo.filter <- haplo.filter %>%
      filter(rank <= 2) %>%
      group_by(locus, id) %>%
      arrange(-depth) %>%
      summarise(
        hap1 = ifelse(length(depth) == 1, haplo[1], sort(haplo)[1]),
        hap2 = ifelse(length(depth) == 1, haplo[1], sort(haplo)[2])
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


    all.hap <- unique(c(freq.hap$re.hap1,freq.hap$re.hap1))
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

    ggplot(haplo.filter,
           aes(
             x = hap1,
             y = hap2,
             size = n,
             color = hap1 == hap2
           )) +
      geom_point() +
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
            panel.margin = unit(c(0,0,0,0), 'mm'),
            panel.border = element_rect(size = 0,colour = "white"),
            plot.margin = unit(c(2, 0, 2, 0), "mm"))
  }, height = function() {
    ifelse(groupPg$width == 0,
           0,
           ifelse(input$selectLocus == "ALL", 0, max(hapPg$width*30, 400)))
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

    obs.freq.tbl <-  haplo.summaryTbl() %>%
      ungroup() %>%
      group_by(group, locus) %>%
      mutate(tot.haplo = n()) %>%
      group_by(group, locus, haplotype.1, haplotype.2) %>%
      summarise(obs.freq = n() / tot.haplo[1],
                n.occur = n(),
                is.homo = (haplotype.1[1] == haplotype.2[1]))

    allelic.freq.tbl <-
      gather(obs.freq.tbl, whichHap,  hap1, 3:4) %>%
      group_by(group, locus, hap1) %>%
      summarise(f = sum(obs.freq / 2),
                n = sum(ifelse(is.homo, n.occur/2, n.occur)))

    default.msg <- "Variance of haplotype by group"

    nearpoint <- nearPoints(allelic.freq.tbl, input$hapByGroupPlotClick, xvar="group", yvar="hap1",
                            threshold = 10,
                            maxpoints = 1)

    if((!is.null(nearpoint)) && nrow(nearpoint)>0) {
      return(HTML(paste0("# of individual(s) : <b>",nearpoint$n,"</b>")))
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

    obs.freq.tbl <-  haplo.summaryTbl() %>%
      ungroup() %>%
      group_by(group, locus) %>%
      mutate(tot.haplo = n()) %>%
      group_by(group, locus, haplotype.1, haplotype.2) %>%
      summarise(obs.freq = n() / tot.haplo[1],
                n.occur = n(),
                is.homo = (haplotype.1[1] == haplotype.2[1]))

    allelic.freq.tbl <-
      gather(obs.freq.tbl, whichHap,  hap1, 3:4) %>%
      group_by(group, locus, hap1) %>%
      summarise(f = sum(obs.freq / 2),
                n = sum(ifelse(is.homo, n.occur/2, n.occur)))

    ggplot(allelic.freq.tbl, aes(
      x = group,
      y = hap1,
      color = hap1,
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
           ifelse(input$selectLocus == "ALL", 0, max(hapPg$width*20, 400)))
  })


  ## filter panel


  output$RDnARplot <- renderPlot({
    if (is.null(input$selectIndiv) ||
        input$selectDB == "" ||
        is.null(input$selectDB) || is.null(locusPg$l))
      return ()


    rd.lower.bound <- max(0, filterParam$minRD-1)
    rd.upper.bound <- filterParam$minRD+1
    ar.lower.bound <- max(0, filterParam$minAR - 0.05)
    ar.upper.bound <- min(1, filterParam$minAR +0.05)


    readDepthRange <- c(0, 1, 5, 10, 25, 50, 100, filterParam$minRD, rd.lower.bound, rd.upper.bound)
    allelicFracRange <- c(0, 0.1, 0.2, 0.4, 0.8, 1,filterParam$minAR, ar.lower.bound, ar.upper.bound)
    rd.af.grid <- expand.grid(unique(allelicFracRange), unique(readDepthRange))

    hap.sel <- Min.filter.haplo()
    if(is.null(hap.sel)) return()
    hap.sel <- hap.sel %>% filter(rank <=2)

    label.grid <- apply(rd.af.grid, 1, function(i){
      af<-i[1]
      rd<-i[2]
      hap.sel %>%
        filter(depth >= rd,
               allele.balance >= af) %>%
        group_by(locus) %>%
        summarise(num.hap = length(unique(haplo)), num.pass.indiv = length(unique(id))) %>%
        ungroup() %>%
        summarise(descr = paste0(sum(num.hap),
                                 " hap,\n",
                                 round(mean(num.pass.indiv,na.rm = T)),
                                 " indiv"
                                 #"(",
                                 #round(mean(num.pass.indiv)*100/panelParam$n.indiv,2),
                                 #"%)"
        ))
    }) %>% bind_rows()

    rd.af.grid.tbl <- cbind(rd.af.grid, label.grid)
    colnames(rd.af.grid.tbl) <- c("af", "rd", "content")

    rd.af.grid.tbl <- rd.af.grid.tbl %>%
      group_by(af,rd)%>%
      mutate(color.grp = 1*(filterParam$minRD==rd) + (1*filterParam$minAR==af))

    ggplot(rd.af.grid.tbl, aes(x=factor(rd), y=factor(af))) +
      geom_tile(color="black", aes(fill=factor(color.grp)), size=0.1, linetype="dashed")+
      geom_text(aes(label=content))+
      theme_bw()+
      theme(legend.position = "none",
            axis.ticks = element_line(size = 0))+
      scale_x_discrete(expand = c(0, 0))+
      scale_y_discrete(expand = c(0, 0))+
      xlab("min. read depth")+
      ylab("min. allelic ratio")+
      scale_fill_manual(values=c("#FDE1E1","#FFAFAD","#FF8582"))



  }, height =400)




  output$downloadData <- downloadHandler(
    filename = function() {
      if (input$selectTbl == "reported indiv haplotype")
        return ("reported_haplotype.csv")
      if (input$selectTbl == "observed variants (filtered)")
        return ("observed_filtered_haplotype.csv")
      if (input$selectTbl == "observed variants (unfiltered)")
        return ("observed_unfiltered_haplotype.csv")
      if (input$selectTbl == "SNP report")
        return ("snp_report.csv")
      if (input$selectTbl == "locus annotation")
        return ("locus_annotation.csv")
    }
    ,
    content = function(file) {
      if (isolate(input$selectTbl) == "reported indiv haplotype") {
        if (is.null(haplo.summaryTbl()))
          return()
        haplo.freq <-
          haplo.freqTbl() %>% mutate(
            obs.freq = round(obs.freq, 3),
            expected.freq = round(expected.freq, 3)
          )
        haplo.all.tbl <-
          haplo.summaryTbl() %>% rename("indiv.ID" = id)
        haplo.all <-
          left_join(haplo.all.tbl,
                    haplo.freq,
                    by = c("locus", "haplotype.1", "haplotype.2"))
        # haplo.isAccept <-
        #   data.frame(
        #     locus = panelParam$locus.label.bare,
        #     is.reject = panelParam$is.reject,
        #     stringsAsFactors = FALSE
        #   )
        # haplo.all <-
        #   left_join(haplo.all, haplo.isAccept, by = c("locus"))
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
          select(-haplotype.1,-haplotype.2,-read.depth.1,-read.depth.2)
        write.csv(haplo.all, file)

      }
      if (isolate(input$selectTbl) ==  "locus annotation") {
        #        if (is.null(Filter.haplo.sum())) return()
        write.csv(annotateTab$tbl, file)
      }


    }
  )

  output$haploTbl <- DT::renderDataTable({
    #observeEvent(input$updateTable, {
    #observeEvent(input$selectTbl, {

    if (input$selectTbl == "reported indiv haplotype") {
      if (is.null(haplo.freqTbl())) return()
      haplo.freq <-
        haplo.freqTbl() %>% mutate(
          obs.freq = round(obs.freq, 3),
          expected.freq = round(expected.freq, 3)
        )
      haplo.all.tbl <-
        haplo.summaryTbl() %>% rename("Individual ID" = id)
      haplo.all <-
        left_join(haplo.all.tbl,
                  haplo.freq,
                  by = c("locus", "haplotype.1", "haplotype.2"))
      # haplo.isAccept <-
      #   data.frame(
      #     locus = panelParam$locus.label.bare,
      #     is.reject = panelParam$is.reject,
      #     stringsAsFactors = FALSE
      #   )
      # haplo.all <-
      #   left_join(haplo.all, haplo.isAccept, by = c("locus"))
    }

    if (input$selectTbl ==  "observed variants (unfiltered)") {
      haplo.filter <- Min.filter.haplo()
      if (is.null(haplo.filter)) return()

      haplo.all <- haplo.filter %>% rename("indiv.ID" = id) %>%
        select(-sum.Phred.C, -max.Phred.C)

    }


    if (input$selectTbl == "observed variants (filtered)") {
      if (is.null(Filter.haplo.sum())) return()

      haplo.all <-
        Filter.haplo.sum() %>% rename("Individual ID" = id) %>%
        select(-sum.Phred.C, -max.Phred.C)

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
        select(-haplotype.1,-haplotype.2,-read.depth.1,-read.depth.2)
    }

    if (isolate(input$selectTbl) ==  "locus annotation") {
      haplo.all <- annotateTab$tbl
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

    #    cat(file=stderr(), input$gibbIter, "\t", input$randomSeed,"we got stuff----\n")

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

  # output$allHapFreqPlot <- renderPlot({
  #   run.result <- Run.SrMicrohap()
  #   if(is.null(run.result)) return()
  #
  # start.iter <- min(floor(run.result$n.sam * (srhapPg$frac.burn/100)),
  #                   run.result$n.sam-2)
  #   freq.matrix <- matrix(unlist(run.result$save.freq),
  #                         byrow = T,
  #                         ncol=run.result$n.haplo)[(start.iter:run.result$n.sam),]
  #
  #
  #   if(run.result$n.haplo ==1) {
  #     hap.freq.stat <- as.data.frame(t(c(quantile(freq.matrix, prob=c(0.05, 0.95) ),
  #                                                          mean=mean(freq.matrix),
  #                                                          hap=run.result$haplo)))
  #     colnames(hap.freq.stat) <- c("X5.","X95.", "mean","hap")
  #
  #     ggplot(hap.freq.stat,aes(x=mean, y=hap), pch=3)+
  #       geom_point()+
  #       theme_bw()+
  #       ylab("")+
  #       xlab("overall haplotype frequency (90% CI)")
  #
  #   } else {
  #   hap.freq.stat <- data.frame(t(apply(freq.matrix,2,
  #                                                   function(x) c(quantile(x, prob=c(0.05, 0.95) ),
  #                                                                 mean=mean(x)))),
  #                                    hap=run.result$haplo)
  #
  #   ggplot(hap.freq.stat,
  #          aes(y=hap, yend=hap, x=as.numeric(X5.), xend=as.numeric(X95.)))+
  #     geom_segment()+
  #     geom_point(data=hap.freq.stat, aes(x=mean, y=hap), pch=3)+
  #     theme_bw()+
  #     ylab("")+
  #     xlab("overall haplotype frequency (90% CI)")
  #
  #   }
  #
  #
  # }, height = function(){ifelse(srhapPg$makePlot,300,0)}
  # )

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

    #cat(file=stderr(), dim(freq.matrix),"<- freq matrix----\n")
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
                             ))), stringsAsFactors = F)
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
                             ))), stringsAsFactors = F)
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
                 stringsAsFactors = F)
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
                          guide = F) +
      ylab("") +
      xlab("haplotype pair")

  }, height = function() {
    ifelse(srhapPg$makePlot,
           max(15 * panelParam$tot.indiv,
               450),
           0)
  })

})
