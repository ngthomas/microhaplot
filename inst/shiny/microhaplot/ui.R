library(shiny)
library(shinyBS)
library(shinyWidgets)

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
    style = "padding-top: 200px;z-index: -1;", #180
    # Application title
    # titlePanel("Haplotype Viewer"),
    absolutePanel(
      style = "z-index:10",
      top = 0,
      left = 0,
      right = 0,
      fixed = TRUE,
      div(style = "padding: 10px 4px 0px 4px; border-bottom: 0px solid #CCC; margin-bottom: 0px; background: #F3FCFF",
          fluidRow(

            column(12,
                   tabsetPanel(
                     tabPanel("> field selection",
                              column(
                                12,
                                column(
                                  6,
                                  column(2,
                                         "Group:", style = "margin-top:21px;font-weight:bold; padding-left:0px",
                                         offset = 0),
                                  column(
                                    8,
                                    selectInput("selectGroup", label = "", "", selected = ""),
                                    style = "margin-top:-10px;margin-bottom:-10px;"
                                  ),
                                  style = "padding-left:0px"
                                ),
                                column(
                                  6,
                                  column(2, "Indiv:", style = "margin-top:21px;font-weight:bold; padding-left:0px"),
                                  column(
                                    8,
                                    selectInput("selectIndiv", label = "", "", selected = ""),
                                    style = "margin-top:-10px;margin-bottom:-10px;"
                                  ),
                                  column(1,
                                         actionButton(
                                           "indivBack", label = "<", width = "80%"
                                         ),
                                         style = "margin-top:10px; padding: 0 0% 0 0%"),
                                  column(1, actionButton(
                                    "indivFor", label = ">", width = "80%"
                                  ),
                                  style = "margin-top:10px;padding:0 0% 0 0%; margin-left: 0px;"),
                                  style = "padding-left:0px"
                                ),
                                column(10,
                                       column(2, "Locus:", style = "margin-top:22px;font-weight:bold; padding-left:0px"),
                                       column(
                                         6,
                                         selectInput("selectLocus", label = "", "", selected = "", width = "100%"),
                                         style = "margin-top:-10px;margin-bottom:-10px; margin-left: -5%;padding-left:0px"
                                       ),
                                       column(2,
                                              column(
                                                5,
                                                actionButton("locusBack", label = "<", width = "80%"),
                                                style = "margin-top:10px; padding: 0 0% 0 0%",
                                                offset = 0
                                              ),
                                              column(5, actionButton(
                                                "locusFor", label = ">", width = "80%"
                                              ),
                                              style = "margin-top:10px;padding:0 0% 0 0%; margin-left: 0px;"),
                                              style = "padding-left:0px"
                                       ),

                                       style = "padding-left:0px"

                                ),
                                style="background:#F7F7F7; padding-top:8px;padding-bottom:10px"
                              ),

                              style="padding-bottom:20px"
                     ),
                     tabPanel("+ read criteria",
                              column(
                                12,
                                column(
                                  7,
                                  column(4,
                                         chooseSliderSkin("Modern", color = "#3A75DC"),
                                         shinyWidgets::sliderTextInput("n.alleles",
                                                                       label = "Top n Alleles:",
                                                                       choices = c(0:20,Inf),
                                                                       selected = 2,
                                                                       grid = FALSE,
                                                                       width = 150),
                                         style = "margin-top:10px;text-align:center; margin-bottom:0px; padding-bottom:0px"),

                                  column(4,
                                         shinyWidgets::sliderTextInput("coverageMin",
                                                                       label = "Min Total Read Depth",
                                                                       choices = c(0:30,50,100,1000,Inf),
                                                                       selected = 2,
                                                                       grid = FALSE,
                                                                       width = 150),
                                         style = "margin-top:10px;text-align:center; margin-bottom:0px; margin-right:0px; padding-bottom:0px;padding-right:0px"),

                                  column(4,
                                         shinyWidgets::sliderTextInput("minAlleleRatio",
                                                                       label = "Min Allelic Ratio:",
                                                                       choices = seq(0,1,0.01),
                                                                       selected = 0,
                                                                       grid = FALSE,
                                                                       width = 150),
                                         style = "margin-top:10px;text-align:center; margin-bottom:0px; padding-bottom:0px"),
                                  style = "padding: 0 0 0 0; margin: -1% 0 0 0;"
                                ),
                                column(
                                  2,
                                  column(10, verticalLayout(actionButton("updateFilter", label = "apply",
                                                                         width="100%",
                                                                         style="text-align:center"),
                                                            actionButton(
                                                              "filterSave",
                                                              "save",
                                                              #icon("hdd-o"),
                                                              width="100%",
                                                              style = "margin-top:10px; text-align:center"
                                                            )), style =
                                           "margin-top:0px; padding-right: 0px;")
                                ),
                                column(
                                  3,
                                  checkboxGroupInput("filterOpts", label = "", choices=list()),

                                  style = "margin-top:-15px; margin-left:0px; padding: 0 0 0 0"
                                ),
                                style="background:#F7F7F7; padding-top:10px")
                     ),
                     tabPanel("+ locus annotation",
                              column(12,
                                     column(
                                       5,
                                       column(3, "Locus:", style ="margin-top:21px;font-weight:bold; padding-left:0px"),
                                       column(9, textOutput("locusSelect1"), style =
                                                "margin-left: 0px; color:grey; margin-top:21px;padding-left: 0px")
                                     ),
                                     column(5,
                                            column(2, "status:",style = "margin-top:23px;font-weight:normal; padding-left:0px"),
                                            column(4,
                                                   selectInput("locusAccept",
                                                               "",
                                                               choices=c("Accept","Reject","NA"),
                                                               selected="NA"),
                                                   style = "margin-top:-8px; margin-bottom:0px"
                                            )
                                     ),

                                     column(2,
                                            actionButton(
                                              "annotateSave",
                                              "save",
                                              icon("hdd-o"),
                                              style = "margin-top:10px; padding-right: 0px; text-align:center",
                                              width = "90%"
                                            ),
                                            offset=0
                                     ),

                                     column(12,
                                            column(1, "Note:", style ="margin-top:5px; padding-left:0px"),
                                            column(11,
                                                   textInput("locusComment", label = "", value = "", width="100%"),
                                                   style = "margin-top:-20px;padding-top:0px;padding-right: 0px;padding-left:0px")),
                                     style="background:#F7F7F7; padding-top:2px;padding-left:0px; margin-left:0px"
                              )
                     )
                   ))
          ))
    ),

    bsAlert("alert"),
    navbarPage(
      "",
      tabPanel("Data Set",
               fluidRow(column(
                 12,
                 column(
                   10,
                   selectInput(
                     "selectDB",
                     label = "Select Data Set:",
                     "",
                     selected = NULL,
                     width = "80%"
                   ),
                   style = "padding-left:0px;margin-top:10%;margin-bottom:-10px;",
                   offset = 1
                 )
               ))),

      # group label panel
      navbarMenu("Summary",
                 tabPanel(h5("By Group"),
                          fluidRow(
                            column(5, plotOutput("fIndivByGroupPlot", height =
                                                   "auto")),
                            column(5, plotOutput("fLociByGroupPlot", height =
                                                   "auto"))
                          ),
                          fluidRow(
                            div(style = "padding: 40px; border-bottom: 20px solid white; background: white")
                          )
                 ),

                 tabPanel(
                   h5("By Individual"),
                   fluidRow(
                     column(4, column(
                       12,
                       column(7, h6("Display (indiv/pg) : ")),
                       column(
                         5,
                         selectInput(
                           "indivPerDisplay",
                           label = NULL,
                           choices = list(
                             "15" = 15,
                             "30" =
                               30,
                             "60" =
                               60,
                             "ALL" =
                               100
                           ),
                           selected =
                             15
                         )
                       )
                     )),
                     column(4, column(
                       12,
                       column(1, h6("Page:")),
                       column(
                         5,
                         numericInput(
                           "indivPage",
                           label = NULL,
                           value = 1,
                           min = 1,
                           step = 1
                         ),
                         offset = 1
                       ),
                       column(1, h6("of ")),
                       column(1, h6(textOutput("maxIndivPage")))
                     )),
                     column(4,column(4,h6("Height:")),
                            column(8,
                                   selectInput(
                                     "indivHeight",
                                     label = NULL,
                                     choices = list("50%" = 0.5, "100%" = 1, "150%" = 1.5,
                                                    "200%" = 2,"500%" = 5),
                                     selected = 1)
                            )),
                     style = "border-bottom: 1px double #d9d9d9;  margin-bottom: 20px; padding-top:15px"
                   ),
                   fluidRow(
                     column(
                       5,
                       plotOutput(
                         "AlleleRatioByIndiv",
                         height = "auto",
                         dblclick = dblclickOpts(id = "plot_dblclick"),
                         brush = brushOpts(
                           id = "plot_brush",
                           direction = "y",
                           resetOnNew = TRUE
                         )
                       )
                     ),
                     column(2, plotOutput("numUniqHapByIndiv", height =
                                            "auto")),
                     column(2, plotOutput("fracHaploPlot", height = "auto")),
                     column(3, plotOutput("readDepthByIndiv", height =
                                            "auto"))
                   ),
                   column(12, h1("")),
                   br(),
                   fluidRow(
                     div(style = "padding: 10px; border-bottom: 8px solid white; background: white")
                   )
                 ),
                 tabPanel(
                   h5("By Locus"),
                   fluidRow(
                     column(4, column(
                       12,
                       column(7, h6("Display (loci/pg) : ")),
                       column(
                         5,
                         selectInput(
                           "locusPerDisplay",
                           label = NULL,
                           choices = list(
                             "15" = 15,
                             "30" =
                               30,
                             "60" =
                               60,
                             "ALL" =
                               100
                           ),
                           selected =
                             15
                         )
                       )
                     )),
                     column(4, column(
                       12,
                       column(1, h6("Page:")),
                       column(
                         5,
                         numericInput(
                           "locusPage",
                           label = NULL,
                           value = 1,
                           min = 1,
                           step = 1
                         ),
                         offset = 1
                       ),
                       column(1, h6(" of ")),
                       column(1, h6(textOutput("maxlocusPage")))
                     )),
                     column(4,column(4,h6("Height:")),
                            column(8,
                                   selectInput(
                                     "lociHeight",
                                     label = NULL,
                                     choices = list(
                                       "50%" = 0.5,
                                       "100%" =
                                         1,
                                       "150%" =
                                         1.5,
                                       "200%" =
                                         2,
                                       "500%" =
                                         5
                                     ),
                                     selected =
                                       1
                                   )
                            )),
                     style = "border-bottom: 1px double #d9d9d9;  margin-bottom: 20px; padding-top:15px"
                   ),
                   fluidRow(
                     column(
                       5,
                       plotOutput(
                         "haplDensityPlot",
                         height = "auto",
                         dblclick = dblclickOpts(id = "plotH_dblclick"),
                         brush = brushOpts(
                           id = "plotH_brush",
                           direction = "y",
                           resetOnNew = TRUE
                         )
                       )
                     ),
                     column(2, plotOutput("numHapPlot", height = "auto")),
                     column(2, plotOutput("fracIndivPlot", height = "auto")),
                     column(3, plotOutput("readDepthPerLocus", height =
                                            "auto"))
                   ),
                   fluidRow(
                     div(style = "padding: 40px; border-bottom: 20px solid white; background: white")
                   )
                 )
      ),
      navbarMenu("Criteria Cutoff",
                 # table panel
                 tabPanel(
                   h5("Global Scope"),
                   fluidRow(
                     column(5, plotOutput("allReadDepth", height = "auto",
                                          # hover = hoverOpts(
                                          #                 id = "RDplot_hover",
                                          #                 delay = 500,
                                          #                 delayType = "debounce"
                                          #               ),
                                          dblclick = dblclickOpts(id = "RDplot_dblclick")),
                            offset=2),
                     column(5, plotOutput("allAllelicRatio", height = "auto",
                                          dblclick = dblclickOpts(id = "ARplot_dblclick")
                     ))
                   ),
                   fluidRow(
                     column(12, h4(textOutput("DP1"))),
                     column(2, plotOutput("haplabel", height = "auto"),
                            style="padding-left:0px; padding-right:0px"
                     ),
                     column(5, plotOutput("hapReadDepth", height = "auto",
                                          click = "hapRDClick")),
                     column(5, plotOutput("hapAllelicRatio", height = "auto",
                                          click = "hapARClick"))
                   ),
                   fluidRow(
                     div(style = "padding: 10px; border-bottom: 2px solid white; background: white")
                   ),
                   fluidRow(
                     column(12, plotOutput("RDnARplot", height = "auto"))
                   ),
                   fluidRow(
                     div(style = "padding: 20px; border-bottom: 8px solid white; background: white")
                   )
                 ),

                 tabPanel(
                   h5("Quality Profiling"),
                   column(12, HTML("The following content disregards the <b>min total read depth</b>. Instead, it will relied on the following filter:")),
                   column(6,
                          shinyWidgets::sliderTextInput("rdMin",
                                                        label = "Min read depth (per haplotype)",
                                                        choices = c(0:30,50,100,1000,Inf),
                                                        selected = 2,
                                                        grid = FALSE),
                          style = "margin-top:10px;text-align:center; margin-bottom:0px; margin-right:0px; padding-bottom:0px;padding-right:0px", offset=3),
                   column(2, verticalLayout(actionButton("updateRdMin", label = "update",
                                                         width="100%",
                                                         style="text-align:center")),
                          style = "margin-top:25px;text-align:center; margin-bottom:0px; margin-right:0px; padding-bottom:0px;padding-right:0px;padding-left:0px",
                          offset=0),
                   column(6, h4("Individual list:")),
                   column(6, h4("Loci list:")),
                   # mod to add interactive table, ranks of ranks (based on rank of rd, calleable hap)
                   column(6, plotOutput("ambigIndivPlot", height = "auto",
                                        dblclick = dblclickOpts(id = "aip_dblclick"),
                                        brush = brushOpts(
                                          id = "aip_Brush",
                                          direction = "x",
                                          resetOnNew = TRUE
                                        ))),
                   column(6, plotOutput("ambigLociPlot", height = "auto",
                                        dblclick = dblclickOpts(id = "alp_dblclick"),
                                        brush = brushOpts(
                                          id = "alp_Brush",
                                          direction = "x",
                                          resetOnNew = TRUE
                                        ))),
                   column(6, h5("Of loci and individuals that have more than 'Top n' qualified alleles:")),
                   column(12, DT::dataTableOutput('indivProfileTbl')),
                   fluidRow(
                     div(style = "padding: 20px; border-bottom: 8px solid white; background: white")
                   )
                 )
      ),
      navbarMenu("Genotype Call",
                 # table panel
                 tabPanel(
                   h5("AR Refinement"),
                   fluidRow(
                     bsAlert("ARalert")
                   ),
                   fluidRow(
                     column(12,
                            column(5, checkboxInput("keepSel", label =  "keeps pt selection between loci"),
                                   style="margin-top:4px;"),
                            column(3, h6("Scale: max read depth displayed"),offset = 0,
                                   style="padding-top:0px; margin-top:8px"),
                            column(4, shinyWidgets::sliderTextInput("max_read_depth",
                                                                    label = "",
                                                                    choices = 10 * 2^(0:13),
                                                                    selected = 1280,
                                                                    grid = FALSE),
                                   style="padding-top:-20px; margin-top:-20px"),
                            style = "border-bottom: 1px double #d9d9d9;  margin-bottom: 20px; padding-top:-10px; margin-top:-10px"
                     ),
                     column(4,
                            column(12, plotOutput("pieCallChart",height="auto"),style="padding:0 0 0 0; margin:0 0 0 0"),
                            column(6, actionButton(
                              "preL",
                              "prev locus",
                              width="100%",
                              style = "margin-top:10px; text-align:center; margin-bottom:15px"
                            )),
                            column(6, actionButton(
                              "nextL",
                              "next locus",
                              width="100%",
                              style = "margin-top:10px; text-align:center; margin-bottom:15px"
                            )),
                            column(12, htmlOutput("savedARstat"), style="margin-bottom:4px"),
                            column(12, sliderInput("max.ar.hm",
                                                   label = "Max AR for Homoz (blue)",
                                                   min = 0,
                                                   max = 1,
                                                   value = 0.3,
                                                   ticks= FALSE,
                                                   step = 0.01)),
                            column(12, sliderInput("min.ar.hz",
                                                   label = "Min AR for Het (yellow)",
                                                   min = 0,
                                                   max = 1,
                                                   value = 0.4,
                                                   ticks=F,
                                                   step = 0.01))
                     ),
                     column(8,
                            ggiraph::ggiraphOutput("biplot", height = "450px"),
                            height="auto",
                            style= "width:'100%'"),
                     column(12,
                            DT::dataTableOutput("biPlotTbl")),
                     column(12, h4(textOutput("DP2")), textOutput("DP3")),
                     column(12, ""),
                     column(1, plotOutput("uchaplabel", height = "auto"),
                            style="padding-left:0px; padding-right:0px"
                     ),
                     column(6, plotOutput("uchapReadDepth", height = "auto")),
                     column(5, plotOutput("uchapAllelicRatio", height = "auto"))
                   ),
                   fluidRow(
                     div(style = "padding: 20px; border-bottom: 8px solid white; background: white")
                   )
                 ),
                 #locus assessement panel

                 tabPanel(
                   h5("Summaries"),
                   fluidRow(
                     bsAlert("hapAlert")
                   ),
                   fluidRow(
                     column(12,
                            column(6, htmlOutput("hapFreqClicked")),
                            column(6, htmlOutput("hwClicked"))),
                     column(6, plotOutput("hapFreq",
                                          height ="auto",
                                          click = "hapFreqPlotClick",
                                          hover=hoverOpts(id="hapFreqPlotHover",
                                                          delay=300,
                                                          delayType = "throttle")),
                            style="padding-right:0px"),
                     column(6, plotOutput("PairWiseHap", height =
                                            "auto",
                                          click = "HWplotClick",
                                          hover=hoverOpts(id="HWPlotHover",
                                                          delay=300,
                                                          delayType = "throttle")),
                            style="padding-left:0px"),
                     column(12, htmlOutput("hapByGroupPlotClicked")),
                     column(12, plotOutput("hapByGroupPlot",
                                           height ="auto",
                                           click = "hapByGroupPlotClick",
                                           hover=hoverOpts(id="hapByGroupPlotHover",
                                                           delay=300,
                                                           delayType = "throttle"))),
                     column(12, plotOutput("hapSeq", height = "auto"))
                   ),
                   fluidRow(
                     div(style = "padding: 20px; border-bottom: 8px solid white; background: white")
                   )
                 )
      ),
      # tabPanel(
      #   "Inferential Analysis",
      #   "SrMicroHap is sensitive only to the changes in \"Locus\" selector tab\n",
      #   wellPanel(fluidRow(
      #     column(
      #       12,
      #       column(
      #         2,
      #         numericInput(
      #           "gibbIter",
      #           "Num of Iter:",
      #           min = 1,
      #           value = 1000,
      #           step=1)
      #       ),
      #       column(
      #         2,
      #         numericInput(
      #           "fracBurn",
      #           "% for Burn-in",
      #           min = 0,
      #           max = 99,
      #           value = 0
      #         )
      #       ),
      #       column(
      #         2,
      #         numericInput(
      #           "randomSeed",
      #           label = "Random Seed",
      #           value = 34532,
      #           min = 1,
      #           step = 1
      #         ),
      #         offset=1
      #       ),
      #       column(
      #         2,
      #         selectInput(
      #           "selectPrior",
      #           label = "Prior Model",
      #           c("uniform", "empirical"),
      #           selected = "empirical"
      #         )
      #       ),
      #       column(
      #         2,
      #         actionButton(
      #           "submitSrMicroHap",
      #           "Submit",
      #           icon("random"),
      #           style = "margin-top:20px; padding-right: 0px;",
      #           width = "100%"
      #         ),
      #         offset=1
      #       )
      #     ),
      #     column(12, "Parameters define the characteristic of true haplotypes:"),
      #     column(12,
      #            column(3, numericInput(
      #              "minRDsr",
      #              label = "min read depth (total indiv)",
      #              value = 10,
      #              min = 1,
      #              step = 1
      #            )),
      #            column(3, numericInput(
      #              "minARsr",
      #              "min allelic ratio",
      #              min = 0,
      #              max = 1,
      #              value = 0.2,
      #              step=0.01
      #            ))
      #     ))),
      #   fluidRow(column(
      #     12, plotOutput("indivHapPosPlot", height = "auto")
      #   ),
      #   #column(12, plotOutput("allHapFreqPlot", height = "auto")),
      #   column(
      #     12, plotOutput("HapFreqByGroupPlot", height = "auto")
      #   )),
      #   fluidRow(
      #     div(style = "padding: 20px; border-bottom: 8px solid white; background: white")
      #   )
      # ),
      # table panel
      tabPanel(
        "Table",
        fixedPanel(
          #top="200px",
          style = "z-index:9;background-color:white; margin-top: -40px; padding-top: 15px; height:90px",
          width="100%",
          fluidRow(
            column(
              6,
              column(4, h5("Select Table:"), offset =
                       0),
              column(
                8,
                selectInput(
                  "selectTbl",
                  label = "",
                  c("observed variants (unfiltered)","observed variants (filtered)",
                    "reported indiv haplotype (diploid)", "reported haplotype (flat)", "SNP report", "locus annotation"),
                  selected = "observed variants (unfiltered)"
                ),
                style = "background-color:white; padding-right: 0px; margin-top:-20px;margin-bottom:-30px;padding-left:0%; padding-right: 0px;"
              )
            ),
            column(6, column(
              12,
              column(4, downloadButton('downloadData', 'Download'))
            )),
            style = "border-bottom: 1px double #d9d9d9;  margin-bottom: 40px; padding-top:15px; padding-bottom:20px")
        ),
        column(12, DT::dataTableOutput('haploTbl'),
               style = "padding-bottom: 40px; border-bottom: 8px solid white; background: white; margin-top: 100px;"),
        fluidRow(
          div(style = "padding: 20px; border-bottom: 8px solid white; background: white")
        )
      ),

      tabPanel(
        "About",
        column(12,
               column(8,
                      htmlOutput("about"),
                      offset=2
               ))
      ),

      position = "fixed-bottom"
    ),
    titlePanel("", windowTitle = "HapPLOType: a view to your haplotypes")

  )
)
