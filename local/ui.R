####################################################
# ui.R File: Layout for PSHACalculatoR             #
#                                                  #
# For details and program execution see: runIT.R   #
# Michael Haas                                     #
# 15.01.2014                                       #
# mhaas@gfz-potsdam.de                             #      
####################################################

library(shiny)

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("PSHACalculatoR"),
  
  #sidebar to select source type
  sidebarPanel(
    tabsetPanel(
      tabPanel("Sources",
         selectInput("sourceType","Source Type:",
                     c("Point" = "pt","Fault" = "ft","Area" = "ar")),
         #Interface for coordinate values
         numericInput("x1", "x1", 10 ),
         numericInput("y1", "y1", 10 ),
         numericInput("x2", "x2 (Fault/Area)",20 ),
         numericInput("y2", "y2 (Fault/Area)",20 ),
         numericInput("x3", "x3 (Area)", 30 ),
         numericInput("y3", "y3 (Area)", 20 ),
         numericInput("x4", "x4 (Area)", 20 ),
         numericInput("y4", "y4 (Area)", 10 )
         ),
      tabPanel("Gutenberg-Richter",
         # Gutenberg-Richter parameters 
         numericInput("aGR","a",3),        
         numericInput("bGR","b",1),
         numericInput("MminGR","Mmin",4),
         numericInput("MmaxGR","Mmax",7)
         ),
      tabPanel("GMPE",
               # GMPE parameters 
               numericInput("gmpeMw","Mw",5),
               sliderInput("gmpeSig", "Sigma", min=0.41, max=0.9, value=0.57, step = 0.01,
                           round = FALSE, format = "#,##0.#####", locale = "us",
                           ticks = TRUE, animate = FALSE)
              ),
      tabPanel("Exceedance Probability/Mean Annual Rate",
               numericInput("OP","Period of interest [y]",50))
         ) # end tabsetPanel
        ), # end sidebarPanel
  
  mainPanel(
    tabsetPanel(
      tabPanel("Source Plot",plotOutput("sourcePlot")),
      tabPanel("Distance probability",plotOutput("distPlot")),
      tabPanel("Gutenberg-Richter",plotOutput("grPlot")),
      tabPanel("GMPE",plotOutput("gmpePlot"),h3(textOutput("excP"))),
      tabPanel("Exceedance Probability/Mean Annual Rate",plotOutput("annRatePlot"),
               plotOutput("excProbPlot"))
     )#end tabsetPanel
    )#end mainPanel
))