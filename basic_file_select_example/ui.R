
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinyFiles)

shinyUI(fluidPage(

  # Application title
  titlePanel("CTDQC"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      shinyDirButton('directory', 'Folder select', 'Please select a folder'),
      selectInput('files', "Select File", choices = 'filelist'),
      sliderInput("bins",
                  "Number of bins:",
                  min = 1,
                  max = 50,
                  value = 30)
    ),

    # Show a plot of the generated distribution
    mainPanel(
      verbatimTextOutput('directorypath'),
      verbatimTextOutput('rawfilelist'),
      plotOutput("distPlot")
    )
  )
))
