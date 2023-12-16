
######### ######### ######### ######### #########
cirmigen_ui = tabPanel(
  "circRNA-miRNA-mRNA interaction",
  #  hidden(
  div(
    id = "app-content",
    style = "font-size:80%;",
    sidebarPanel(
      width = 3,
      fluidRow(
        #div(style="margin-top:-1em"),
        #h4("Input circRNAs"),
        div("Input circRNAs", style = "text-align: center; background-color: grey; color:white; font-size:110%"),
        #div(style="margin-top:-1em"),
        textAreaInput(
          inputId = "circRNAinput",
          label = p("circRNAs are separated by a space, a line, or a comma", style = "font-size:8pt"),
          value = "5__102210394__102227535", #5__102210394__102227535, 14__99472961__99473012
          cols = 60,
          rows = 3
        ),
        
        ### input circRNA sequence
        div(style = "margin-top:-1em"),
        checkboxInput(inputId = "cbShowcircRNAseqlist", label = "Input cirRNA sequences"),
        div(style = "margin-top:-1em"),
        hidden(
          textAreaInput(
            inputId = "circRNAseqlist",
            label = p("circRNA sequences in fasta format", style = "font-size:8pt"),
            value = "",
            cols = 60,
            rows = 10
          )
        ),
        
        actionButton("goButton", "Execute", styleclass = "info")
        
        
      ), #end of fluidRow
      
      div(style = "margin-top:1em"),
      fluidRow(
        div("Settings for circRNA-miRNA", style = "text-align: center; background-color: grey; color:white; font-size:110%"),
        
        checkboxInput(
          inputId = "cbUseTargetScan",
          label = "Use TargetScan results",
          value = TRUE
        ),
        div(style = "margin-top:-1.5em"),
      ),
      
      fluidRow(column(width = 6,
                      hidden(
                        textInput(
                          inputId = "TargetScan_6mer",
                          label = "Min 6mer",
                          value = "0",
                          width = "100px"
                        )
                      )),
               column(width = 6,
                      hidden(
                        textInput(
                          inputId = "TargetScan_7mer1a",
                          label = "Min 7mer-1a",
                          value = "0",
                          width = "100px"
                        )
                      ))),
      
      fluidRow(column(width = 6,
                      hidden(
                        textInput(
                          inputId = "TargetScan_7merm8",
                          label = "Min 7mer-m8",
                          value = "0",
                          width = "100px"
                        )
                      )),
               column(width = 6,
                      hidden(
                        textInput(
                          inputId = "TargetScan_8mer1a",
                          label = "Min 8mer-1a",
                          value = "0",
                          width = "100px"
                        )
                      )),), 
  
      fluidRow(
        checkboxInput(
          inputId = "cbUseRNAhybrid",
          label = "Use RNAhybrid results",
          value = TRUE
        ),
        #    tags$div(title="Add the description of RNAhybrid here",p("Description of RNAhybrid settings",style = "font-size:8pt")),
        div(style = "margin-top:-1.5em"),
        
        fluidRow(column(width = 6,
                        hidden(
                          textInput(
                            inputId = "RNAhybrid_pvalThres",
                            label = "Max p-val",
                            value = "0.05",
                            width = "100px"
                          )
                        )),
                 #selectInput(inputId = "RNAhybrid_pvalThres", label = "Max p-value", choices =  c(seq(0.001, 0.009, 0.001),seq (0.01, 0.05,0.01),1), multiple=FALSE, selected = 0.05)),
                 column(width = 6,
                        hidden(
                          textInput(
                            inputId = "RNAhybrid_mfeThres",
                            label = "Max mfe",
                            value = "-18",
                            width = "100px"
                          )
                        )))
      ),
      #selectInput(inputId = "RNAhybrid_mfeThres", label = "Max mfe", choices = c(-seq(10, 50,1)), multiple=FALSE, selected = -18))
      #    ),
      fluidRow(
        checkboxInput(
          inputId = "cbUseMiRanda",
          label = "Use miRanda results",
          value = TRUE
        ),
        #    tags$div(title="Add the description of miRanda here",p("Description of the miRanda settings",style = "font-size:8pt")),
        div(style = "margin-top:-1.5em"),
        hidden(
          sliderInput(
            inputId = "miRandaScoreThres",
            label = "Min miRanda Score",
            min = 1,
            max = 300,
            value = 140
          )
        ),
        
        
      ), #end of fluidRow
      
      fluidRow(
        div("Settings for miRNA-mRNA", style = "text-align: center; background-color: grey; color:white; font-size:110%"),
        #    tags$div(title="miRNA-mRNA interactions are collected from targetscan version 72",em(h5("miRNA-mRNA interaction"))),
        #tags$div(title="Add the description of targetscan here",p("Description of the targetscan settings",style = "font-size:8pt")),
        # Select the limited number of genes of each miRNA to plot
        selectInput(
          inputId = "genenumlim",
          label = "Max number of genes",
          choices = c("All", as.character(c(
            1, 2, 3, 4, seq(5, 100, 5)
          ))),
          multiple = FALSE,
          selected = "40"
        ),
        sliderInput(
          inputId = "weightScoremigeneThres",
          label = "Max weighted-context score",
          min = -1,
          max = 0,
          value = 0
        ),
        sliderInput(
          inputId = "weightPercentilemigeneThres",
          label = "Min weighted-context score percentile",
          min = 0,
          max = 100,
          value = 99
        )
      )
      
    ),
    #end of sidebarPanel
    # Show a plot of the generated distribution
    mainPanel(
      width = 9,
      
      fluidRow(
        ### input interested miRNA
        checkboxInput(
          inputId = "cbShowcircbaseID_cirmi",
          label = "Display using circRNA IDs from circBase",
          value = TRUE
        ),
        ### input interested miRNA
        div(style = "margin-top:-2em"),
        checkboxInput(inputId = "cbShowmiRNAList", label = "Select miRNAs of interest"),
        div(style = "margin-top:-2em"),
        hidden(textAreaInput(
          inputId = "miRNAList",
          label = p("miRNAs are separated by a space, a line, or a comma", style = "font-size:8pt"),
          value = "",
          cols = 60,
          rows = 3
        )),
        
        ### input interested genes
        div(style = "margin-top:-2em"),
        checkboxInput(inputId = "cbShowGeneList", label = "Select genes of interest"),
        div(style = "margin-top:-1em"),
        hidden(textAreaInput(
          inputId = "geneList",
          label = p("Genes are separated by a space, a line, or a comma", style = "font-size:8pt"),
          value = "",
          #value="EGFR FAK IGF1R YY1 ROCK1 CDK6 IGF2BP3 STAT3 SphK1 CDK4 ELF3 MCL1 ADAM17 HPSE VEGF"
          cols = 60,
          rows = 3
        )),
        div(style = "margin-top:-2em"),
        checkboxInput(inputId = "cbdlhtml", label = "Save the network to html"),
        hidden(downloadButton('downloadNetwork', 'Download network')
               #        downloadButton('downloadCirmi', 'Download circRNA-miRNA interactions'),
               #        downloadButton('downloadmim', 'Download miRNA-mRNA interactions')
        ),
        #      div(style="margin-top:-2em"),
        #      checkboxInput("getSelectedNodes", "Print out selected miRNAs & genes"),
        visNetworkOutput("circRNAnetwork", width = "100%", height = "500px")
      ),
      
      # fluidRow(
      #   conditionalPanel(
      #     condition = "input.getSelectedNodes == 1",
      #     
      #     tags$style(type = 'text/css', '#txt_out {white-space: pre-wrap;}'),
      #     column(
      #       width = 6,
      #       h5("Selected miRNAs"),
      #       verbatimTextOutput("selectedMiRNAs"),
      #     ),
      #     column(
      #       width = 6,
      #       h5("Selected genes"),
      #       verbatimTextOutput("selectedGenes")
      #     )
      #     
      #   )
      # ),
      
      # fluidRow(conditionalPanel(
      #   condition = "input.cbShowMiRNA == 1 ",
      #   h5("List of miRNA"),
      #   verbatimTextOutput("showMiRNA")
      # )),
      # fluidRow(conditionalPanel(
      #   condition = " input.cbShowGene == 1 ",
      #   h5("List of genes"),
      #   verbatimTextOutput("showGene")
      # )),
      

      fluidRow(
        div(
          "CircRNA-miRNA interactions",
          style = "text-align: center;
                  background-color: #4472C4; color:white; font-size:100%"
        ),
        div(style = "margin-top:-1em"),
        checkboxInput(inputId = "cbShowCirmiTable", label = "Show table"),
        div(style = "margin-top:-1em"),
        conditionalPanel(
          condition = "input.cbShowCirmiTable == 1 ",
        column(12, align = "center",
               DT::dataTableOutput("cirmiTableOut"))
        ),
        div(style = "margin-top:-1em"),
        checkboxInput(inputId = "cbShowTargetScanPlots", label = "Plots of TargetScan"),
        div(style = "margin-top:-2em"),
        conditionalPanel(
          condition = "input.cbShowTargetScanPlots == 1 ",
          plotOutput("cirmiTargetScanPlots", width = "100%", height = "400px")
        ),
        
        div(style = "margin-top:-2em"),
        checkboxInput(inputId = "cbShowRNAhybridPlots", label = "Plots of RNAhybrid"),
        div(style = "margin-top:-2em"),
        conditionalPanel(
          condition = "input.cbShowRNAhybridPlots == 1 ",
          plotOutput("cirmiRNAhybridPlots", width = "100%", height = "400px")
        ),
        div(style = "margin-top:-2em"),
        checkboxInput(inputId = "cbShowMiRandaPlots", label = "Plots of miRanda"),
        div(style = "margin-top:-1em"),
        conditionalPanel(
          condition = "input.cbShowMiRandaPlots == 1 ",
          plotOutput("cirmiMiRandaPlots", width = "100%", height = "400px")
        )
        
      ),
      
      fluidRow(
        div(
          "miRNA-mRNA interactions",
          style = "text-align: center;
                    background-color: #4472C4; color:white; font-size:100%"
        ),
        div(style = "margin-top:-1em"),
        checkboxInput(inputId = "cbShowMimRNATable", label = "Show table"),
        div(style = "margin-top:-1em"),
        conditionalPanel(
          condition = "input.cbShowMimRNATable == 1 ",
        column(12, align = "center",
               DT::dataTableOutput("migenTableOut"))
        ),
        
        div(style = "margin-top:-1em"),
        checkboxInput(inputId = "cbShowMimRNAPlot", label = "Plots of TargetScan"),
        div(style = "margin-top:-1em"),
        conditionalPanel(
          condition = "input.cbShowMimRNAPlot == 1 ",
          plotOutput("mimRNAPlots", width = "100%", height = "400px")
         )
      )
      
      # fluidRow(
      #   div(
      #     "Statistics summary of circRNA-miRNA and miRNA-mRNA interactions",
      #     style = "text-align: center;
      #             background-color: #4472C4; color:white; font-size:100%"
      #   )
      # )
      
    )#end of mainPanel
    
    
  )#div
  # )#hidden
)#tabPanel


