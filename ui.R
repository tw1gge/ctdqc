
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinyFiles)
library(leaflet)

vchannels = c("v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7")

shinyUI(fluidPage(

  # Application title
  titlePanel("CTDQC"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      textOutput('directory'),
      shinyDirButton('directory', 'Select CNV Folder', 'Please select a folder'),
      actionButton('read_files', "Read CNV files", icon=icon("book")),
      actionButton('revert', "Revert", icon=icon("repeat")),
      actionButton('dump', "DUMP", icon=icon("floppy-o")),
      selectInput('select_profile', "Select profile", choices = NULL),

      h4("Trim"),
      actionButton('trim', "Trim", icon=icon("scissors")),
      actionButton('autotrim', "Autotrim", icon=icon("scissors")),

      h4("Decimate"),
      numericInput('bin_size', "Bin Size (m)", min = 0.5, max = 5, value = 0.5, step = 0.5),
      actionButton('decimate', "Decimate", icon=icon("delicious"))


    ),

    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("Scan Plot", plotOutput("scan_plot", brush = brushOpts("scan_brush", direction = "x"))),
        tabPanel("Profile Plot",
                 fluidRow(
                   column(4, selectInput("y", "Y axis", choices = NULL)),
                   column(4, selectInput("x1", "Blue X axis", choices = NULL)),
                   column(4, selectInput("x2", "Red X axis", choices = NULL))
                 ),
                 plotOutput("profile_plot")),
        tabPanel("TS Plot", plotOutput("TS_plot")),
        tabPanel("Map", leafletOutput("map")),
        tabPanel("Sensors",
                 inputPanel(
                   h4("Optode"),
                   selectInput("optode_foil", "Optode foil Batch #", choices = unique(optode_coefs$batch)),
                   selectInput('optode_T_channel', "Optode Temperature channel", choices = vchannels, selected = "v6"),
                   selectInput('optode_dphase_channel', "Optode dPhase channel", choices = vchannels, selected = "v7"),
                   actionButton('optode', "Process Optode", icon=icon("life-ring"))
                 ),
                 inputPanel(
                   h4("RINKO"),
                   selectInput('rinko_vchannels', "RINKO Temperature channel", choices = vchannels, selected = "v6"),
                   actionButton('rinko', "Process RINKO", icon=icon("times-circle-o"))
                 ),
                 inputPanel(
                   h4("Flurometer")
                 ),
                 verbatimTextOutput("xml")
                 ),
        tabPanel("Summary", verbatimTextOutput("summary"))
      ),
      verbatimTextOutput("debug")
    )
  )
))
