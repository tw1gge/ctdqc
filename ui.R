library(shiny, quietly=T)
library(shinyFiles, quietly=T)
library(leaflet, quietly=T)
library(rhandsontable, quietly=T)

shinyUI(fluidPage(

  # Application title
  titlePanel("CTDQC"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(width = 4,
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

      h4("Progress"),
      actionButton('mark_complete_QC2', "Mark QC2 done", icon=icon("check-square")),
      actionButton('mark_complete_all', "Mark All QC complete", icon=icon("check-square")),
      DT::dataTableOutput("progress")
      ),

    mainPanel(width = 8,
      tabsetPanel(
        tabPanel("Summary",
                 verbatimTextOutput("summary"),
                 verbatimTextOutput("xml")
                 ),
        tabPanel("Scan Plot",
                 plotOutput("scan_plot", brush = brushOpts("scan_brush", direction = "x"), height="800px"),
                 h4("Trim"),
                 actionButton('pumped', "Subset pump", icon=icon("battery-1")),
                 actionButton('trim', "Trim", icon=icon("scissors")),
                 actionButton('autotrim', "Autotrim", icon=icon("scissors")),
                 actionButton('remove_pressure_inversions', "Remove pressure inversions", icon=icon("angle-double-down")),
                 numericInput('decent_threshold', "minimum decent rate threshold (m/s)", value=0.15, step=0.05),
                 numericInput('inversion_window', "decent rate window (s)", value=2, step=0.5)
                 ),
        tabPanel("Profile Plot",
                 fluidRow(
                   column(4, selectInput("y", "Y axis", choices = NULL)),
                   column(4, selectInput("x1", "Primary (Blue) X axis", choices = NULL)),
                   column(4, selectInput("x2", "Secondary (Red) X axis", choices = NULL))
                   ),
                 plotOutput("profile_plot", brush = brushOpts("flag_brush", direction = "xy"), height="800px"),
                 h6("Flags and factors are applied to the primary axis only and for all dips"),
                 h4("Flag"),
                 actionButton('apply_flag', "Apply Flag", icon=icon("flag")),
                 inputPanel(
                   h4("Factor / Offset"),
                   fluidRow(
                     column(6, numericInput('factor', label="Factor", value=1.0, step=0.01)),
                     column(6, numericInput('offset', label="Offset", value=0.0, step=0.01))
                     ),
                   actionButton('apply_factor', "Apply Factor + Offset", icon=icon("flag"))
                   ),
                 h4("CTD / Niskin regressions"),
                 tableOutput("bottle_coef")
                 ),
        tabPanel("TS Plot", plotOutput("TS_plot")),
        tabPanel("Filter Plot",
                 fluidRow(
                   column(3, selectInput("filter_x1", "Parameter", choices = NULL)),
                   column(3, numericInput('filter_t', label="time constant (s)", value=0.15, step=0.01)),
                   column(6,
                     br(),
                     actionButton('prev_filter', "Preview filter", icon=icon("flag")),
                     actionButton('apply_filter', "Apply filter", icon=icon("flag"))
                     )
                   ),
                 plotOutput("filter_plot", height="600px", dblclick="filter_plot_dblclick", brush= brushOpts(id="filter_plot_brush", resetOnNew=T))
                 ),
        tabPanel("Table", dataTableOutput("datatable")),
        tabPanel("Map", leafletOutput("map")),
        tabPanel("Sensors",
          tableOutput("config"),
          uiOutput("sensor_ui")
          ),
        tabPanel("Bottles",
          rHandsontableOutput("bottles"),
          br(),
          fluidRow(
            h4("CTD / Niskin regressions"),
            selectInput("Plot_bottle_select", NULL, choices = c("Select parameter" = "", "Salinity", "Oxygen Optode", "Oxygen RINKO", "Chlorophyll")),
            plotOutput("bottle_plot")
            )
          ),
        tabPanel("Publish",
          # dynamically generated UI
          # selectInput('publish_variables', 'Variables to publish', sensor_metadata$parameter, multiple=T, selectize=T),
          actionButton('make_netcdf', "Publish NetCDF", icon=icon("object-group")),
          actionButton('write_csv', "Write csv's", icon=icon("table")),
          fluidRow(
            column(6,
              h4("Global metadata"),
              uiOutput("edit_metadata")
              ),
            column(6)
            ),
          verbatimTextOutput("metadata")
        )
      )
    )
  )
))


