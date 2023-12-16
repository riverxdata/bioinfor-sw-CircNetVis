GSEA_ui = tabPanel(
  "Pathway analysis",
  #  hidden(
  div(
    id = "app-content",
    style = "font-size:80%;",
    sidebarPanel(width = 3,
                 fluidRow(
                   #h4("Setting"),
                   div("Settings for GSEA", style = "text-align: center; background-color: grey; color:white; font-size:110%"),
                   #h5("Gene set enrichment analysis (GSEA)"),
                   #h6("By default, genes from the current network are used for GSEA. If you want to use all genes interact with  the miRNAs, use this option:"),
                   # Select the limited number of pathways to plot
                   #sliderInput(inputId = "topPW", label = "Top pathways to plot:", min=1, max=100, value=10),
                   #checkboxInput("cbgset", label = "Select all genes from TargetScan", value = FALSE),
                   radioButtons(
                     "rbgset",
                     label = "Select gene set",
                     c(
                       "Current genes in the miRNA-mRNA network" = "net",
                       "All genes interacting with the miRNAs" = "All"
                     )
                   ),
                   actionButton("goPWButton", "Go", styleclass = "info")
                 )
    ),
    #end of sidebarPanel
    
    # Show a plot of the generated distribution
    mainPanel(width = 9,
              div(
                "Gene set enrichment analysis",
                style = "text-align: center;
                    background-color: #4472C4; color:white; font-size:100%"
              ),
              fluidRow(DT::dataTableOutput("reactomeAPI")),
              div(
                "GO annotation",
                style = "text-align: center;
                    background-color: #4472C4; color:white; font-size:100%"
              ),
              fluidRow(DT::dataTableOutput("goAnno"))
              
              )#end of mainPanel
              
  )#div
  # )#hidden
)#tabPanel
