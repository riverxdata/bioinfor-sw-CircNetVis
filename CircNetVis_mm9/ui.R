library(shiny)
library(shinythemes)
#library(shinysky)
library(shinyjs)
library(DT)
library(visNetwork)
#library(shinybusy)

source("./ui/cirmigen_ui.R", local=environment())
source("ui/RBP_ui.R", local=environment())
source("ui/GSEA_ui.R", local=environment())
source("ui/help_ui.R", local=environment())

intro_ui = tabPanel("Getting Started",
      column(10, offset = 1,
      tags$div(class="jumbotron", 
               tags$strong(h4("CircNetVis: an interactive web application for visualizing interaction networks of circular RNAs")) 
               ),
      
      selectInput(
        inputId = "inputSpecies",
        label = "Select species",
        choices = c("Mouse - mm9",
                    "Human - hg19"),
        multiple = FALSE,
        selected = "Mouse - mm9"
      ),
      

      #includeHTML("introduction.html"),
      #tags$head(
      #  tags$link(rel = "stylesheet", type = "text/css", href = "ui.css")
      #)
      ),
      useShinyjs()
      #inlineCSS(appCSS)
)

ui = navbarPage("Menu",
                intro_ui,
                cirmigen_ui,
#                RBP_ui,
                GSEA_ui,
                help_ui,
                theme = shinytheme("flatly")
                #cerulean, cosmo, cyborg, darkly, flatly, journal, lumen, paper, readable, sandstone, simplex, slate, spacelab, superhero, united, yeti
  )

