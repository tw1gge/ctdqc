
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinyFiles)
library(leaflet)

shinyUI(fluidPage(

  # Application title
  titlePanel("CTDQC"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(width = 4,
                 fluidRow(
      shinyDirButton('directory', 'Select CNV Folder', 'Please select a folder'),
      textOutput('directory'),
      actionButton('read_files', "Read CNV files", icon=icon("book")),
      actionButton('trim', "Trim CTD", icon=icon("floppy-o")),
      actionButton('revert', "Revert", icon=icon("floppy-o")),
      actionButton('dump', "DUMP", icon=icon("floppy-o")),
      selectInput('select_profile', "Select profile", choices = NULL)
                 )
    ),

    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("Scan Plot", plotOutput("scan_plot")),
        tabPanel("Profile Plot", plotOutput("profile_plot")),
        tabPanel("Map", leafletOutput("map")),
        tabPanel("Summary", verbatimTextOutput("summary"))
      )
    )
  )
))
