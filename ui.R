library(shiny, quietly=T)
library(shinyFiles, quietly=T)
library(leaflet, quietly=T)
library(rhandsontable, quietly=T)

shinyUI(fluidPage(

  # Application title
  titlePanel("CTDQC"),

  # Sidebar with a slider input for number of bins
  # Sidebar ----
  sidebarLayout(
    sidebarPanel(width = 4,
      textOutput('directory'),
      shinyDirButton('directory', 'Select Folder', 'Please select a folder'),
      actionButton('read_files', "Read .cnv files", icon=icon("book")),
      actionButton('read_bottle', "Read .bl files", icon=icon("flask")),
      actionButton('read_rdata', "Read Rdata file", icon=icon("registered")),
      actionButton('revert', "Revert", icon=icon("repeat")),
      selectInput('select_profile', "Select profile", choices = NULL),

      h4("Save"),
      actionButton('write_rdata', "Write Rdata", icon=icon("bookmark")),

      h4("Progress"),
      actionButton('mark_complete_QC2', "Mark QC2 done", icon=icon("check-square")),
      actionButton('mark_complete_all', "Mark All QC complete", icon=icon("check-square")),
      DT::dataTableOutput("progress")
      ),

    mainPanel(width = 8,
      tabsetPanel(
        # summary ----
        tabPanel("Summary",
                 verbatimTextOutput("summary"),
                 verbatimTextOutput("xml")
                 ),
        # Scan plot ----
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
        # Profile plot ----
        tabPanel("Profile Plot",
                 fluidRow(
                   column(4, selectInput("y", "Y axis", choices = NULL)),
                   column(4, selectInput("x1", "Primary (Blue) X axis", choices = NULL)),
                   column(4, selectInput("x2", "Secondary (Red) X axis", choices = NULL))
                   ),
                 plotOutput("profile_plot", brush = brushOpts("flag_brush", direction = "xy"), height=800),
                 h6("Factors are applied to the primary axis only and for all dips"),
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
        # TS plot ----
        tabPanel("TS Plot", plotOutput("TS_plot", height=700, width=700)),
        # Drift plot ----
        tabPanel("Drift Plot", plotOutput("drift_plot", height=700, width=700)),
        # Hysteresis plot ----
        tabPanel("Hysteresis Plot",
                 fluidRow(
                   checkboxInput("CT_mode", label="CT alignment mode", value=F),
                   column(3, numericInput("lag", label="Lag (scans)", value=0, step=1)),
                   column(3, actionButton("apply_lag", "Apply lag", icon=icon("clock")))
                 ),
                 plotOutput("hyst_plot", height=800)),
        # Filter plot ----
        tabPanel("Filter Plot",
                 fluidRow(
                   column(3, selectInput("filter_x1", "Parameter", choices = NULL)),
                   column(3, numericInput('filter_t', label="time constant (s)", value=0.15, step=0.01)),
                   column(6,
                     br(),
                     actionButton('prev_filter', "Preview filter", icon=icon("eye")),
                     actionButton('apply_filter', "Apply filter", icon=icon("filter"))
                     )
                   ),
                 plotOutput("filter_plot", height=800, dblclick="filter_plot_dblclick", brush= brushOpts(id="filter_plot_brush", resetOnNew=T))
                 ),
        # Table ----
        tabPanel("Table", dataTableOutput("datatable")),
        # Map ----
        tabPanel("Map", leafletOutput("map")),
        # Sensors ----
        tabPanel("Sensors",
          checkboxInput("apply_global", "Apply Global", value = T),
          tableOutput("config"),
          wellPanel(
            fluidRow(
              column(5,
                h4("Remove variable"),
                selectInput("remove_variable_var", NULL, choices = NULL),
                actionButton("remove_variable", "Remove", icon=icon("minus-square"))
                )
            )
          ),
          wellPanel(
            fluidRow(
              column(5,
                h4("Calculate oxygen saturation"),
                  actionButton("oxygen_sat", "Calculate", icon=icon("skyatlas"))
                )
            )
          ),
          uiOutput("sensor_ui")
          ),
        # Bottles ----
        tabPanel("Niskin",
          actionButton('write_bottles', "Save Niskin .CSV", icon=icon("prescription-bottle")),
          rHandsontableOutput("bottles"),
          br(),
          fluidRow(
            h4("CTD / Niskin regressions"),
            selectInput("Plot_bottle_select", NULL, choices = c("Select parameter" = "", "Salinity", "Oxygen Optode", "Oxygen RINKO", "Chlorophyll")),
            plotOutput("bottle_plot")
            )
          ),
        # Publish ----
        tabPanel("Publish",
          # dynamically generated UI
          actionButton('decimate', "Decimate", icon=icon("delicious")),
          numericInput('bin_size', "Bin Size (m)", min = 0.5, max = 5, value = 0.5, step = 0.5),
          br(),
          actionButton('write_netcdf', "Publish NetCDF", icon=icon("object-group")),
          actionButton('write_csv', "Write csv's", icon=icon("table")),
          fluidRow(
            column(6,
              h4("Global metadata"),
              uiOutput("edit_metadata"),
              rHandsontableOutput("publish_param")
              ),
            column(6)
            ),
          verbatimTextOutput("metadata")
        )
      )
    )
  )
))


