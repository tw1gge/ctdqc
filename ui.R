
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
      shinyDirButton('directory', 'Select Folder', 'Please select a folder'),
      actionButton('read_files', "Read .cnv files", icon=icon("book")),
      actionButton('read_bottle', "Read .bl files", icon=icon("flask")),
      actionButton('read_rdata', "Read Rdata file", icon=icon("registered")),
      actionButton('revert', "Revert", icon=icon("repeat")),
      selectInput('select_profile', "Select profile", choices = NULL),

      h4("Decimate"),
      numericInput('bin_size', "Bin Size (m)", min = 0.5, max = 5, value = 0.5, step = 0.5, width = "50%"),
      actionButton('decimate', "Decimate", icon=icon("delicious")),

      h4("Save"),
      actionButton('write_rdata', "Write Rdata", icon=icon("bookmark")),
      actionButton('write_csv', "Write csv's", icon=icon("table"))

    ),

    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("Scan Plot",
                 plotOutput("scan_plot", brush = brushOpts("scan_brush", direction = "x")),
                 h4("Trim"),
                 actionButton('pumped', "Subset pump", icon=icon("battery-1")),
                 actionButton('trim', "Trim", icon=icon("scissors")),
                 actionButton('autotrim', "Autotrim", icon=icon("scissors"))
                 ),
        tabPanel("Profile Plot",
                 fluidRow(
                   column(4, selectInput("y", "Y axis", choices = NULL)),
                   column(4, selectInput("x1", "Primary (Blue) X axis", choices = NULL)),
                   column(4, selectInput("x2", "Secondary (Red) X axis", choices = NULL))
                 ),
                 plotOutput("profile_plot"),
                 h6("Flags and factors are applied to the primary axis only"),
                 h4("Flag"),
                 actionButton('apply_flag', "Apply Flag", icon=icon("flag")),
                 h4("Factor / Offset"),
                 fluidRow(
                   column(6, numericInput('factor', label = "Factor", value = 1)),
                   column(6, numericInput('offset', label = "Offset", value = 0))
                 ),
                 actionButton('apply_factor', "Apply Factor + Offset", icon=icon("flag"))
                 ),
        tabPanel("TS Plot", plotOutput("TS_plot")),
        tabPanel("Table", dataTableOutput("datatable")),
        tabPanel("Map", leafletOutput("map")),
        tabPanel("Bottles", verbatimTextOutput("bottles")),
        tabPanel("Sensors",
                 inputPanel(
                   h4("Optode"),
                   fluidRow(
                     column(6, selectInput('optode_T_channel', "Optode Temperature channel", choices = vchannels, selected = "v6") ),
                     column(6, selectInput('optode_dphase_channel', "Optode dPhase channel", choices = vchannels, selected = "v7") )
                   ),
                   selectInput("optode_foil", "Optode foil Batch #", choices = NULL),
                   actionButton('optode', "Process Optode", icon=icon("life-ring"))
                 ),
                 inputPanel(
                   h4("RINKO"),
                   fluidRow(
                     column(6, selectInput('rinko_T_channel', "RINKO Temperature channel", choices = vchannels, selected = "v4") ),
                     column(6, selectInput('rinko_O_channel', "RINKO Oxygen channel", choices = vchannels, selected = "v5") )
                   ),
                   fluidRow(
                     column(6, numericInput('rinko_G', label = "G", value = 0)),
                     column(6, numericInput('rinko_H', label = "H", value = 1))
                   ),
                   actionButton('rinko', "Process RINKO", icon=icon("times-circle-o"))
                 ),
                 inputPanel(
                   h4("Licor PAR"),
                   fluidRow(
                     column(6, selectInput('par_channel', "PAR channel", choices = vchannels, selected = "v0") )
                   ),
                   numericInput('licor_factor', label = "Licor factor", value = 0.22345679),
                   numericInput('licor_offset', label = "Licor offset", value = 3.3737),
                   actionButton('licor', "Process Licor PAR", icon=icon("beer"))
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


