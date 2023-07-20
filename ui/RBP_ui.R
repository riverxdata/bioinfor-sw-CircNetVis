RBP_ui = tabPanel(
  "circRNA-RBP interaction",
  div(id = "app-content",
      style = "font-size:80%;",
      sidebarPanel(
        fluidRow(
          #h4("circRNA-RBP interaction"),
          div("Settings for circRNA-RBP", style = "text-align: center; background-color: grey; color:white; font-size:110%"),
          
          ####Show only top sites
          checkboxInput(inputId = "cbShowTopSiteRBP", label = "Show only top site RBP"),
          ####Show site num on edge
          checkboxInput(inputId = "cbShowSiteNum", label = "Show number of sites on network"),
          
          ### input interested RBP
          checkboxInput(inputId = "cbShowRBPList", label = "Input RBP of interests"),
          
          hidden(
            textAreaInput(
              inputId = "RBPList",
              label = p("RBP are separated by a space, a line, or a comma", style = "font-size:8pt"),
              value = "",
              cols = 60,
              rows = 3
            )
          ),
          
          # Select the limited number of sites
          #sliderInput(inputId = "siteNumLim", label = "Minimum number of sites", min=0, max=2000, value=100),
          #selectInput(
          textInput(
            inputId = "siteNumLim",
            label = "Minimum number of sites",
            #choices =  as.character(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 50, 100, 2000)),
#            multiple = FALSE,
#            selected = "1"
            value = "1"
          ),
          
        ),
      ),
      #end of sidebarPanel
  
      # Show a plot of the generated distribution
      mainPanel(fluidRow(
        #checkboxInput(inputId = "cbdlRBPhtml", label = "Save the current network to an HTML file"),
        #hidden(
        #  downloadButton('downloadRBPNetwork', 'Download network')
        #),
        checkboxInput(inputId = "cbShowcircbaseID_cirRBP", label = "Display using circRNA IDs from circBase"),
        visNetworkOutput("circRBPnetwork", width = "100%", height = "600px")
      ),
      fluidRow(
        div(
          "circRNA-RBP interactions",
          style = "text-align: center;
                    background-color: #4472C4; color:white; font-size:100%"
        ),
        column(12, align = "center",
               DT::dataTableOutput("circRBPTableOut"))
      ),
      fluidRow(
        div(
          "Distribution of the number of sites in the interactions",
          style = "text-align: center;
                  background-color: #4472C4; color:white; font-size:100%"
        ),
        plotOutput("circRBPStatsPlots", width = "100%", height = "400px")
      ))#end of mainPanel
      
      
  )#div
  # )#hidden
)#tabPanel
