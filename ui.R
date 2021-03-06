library(shiny, quietly=T)
library(shinyFiles, quietly=T)
library(leaflet, quietly=T)
library(rhandsontable, quietly=T)

vchannels = c("v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7")
temperature_serials = c("5558", "5823", "6267", "6268")
conductivity_serials = c("4499", "4523", "4724", "4725")
pressure_serials = c("1274", "1343")
par_serials = c("49")
altimeter_serials = c("68799", "73082")
turbidity_serials = c("11618", "14426")
fluorometer_serials = c("2315", "3817")
optode_serials = c("752")
rinko_serials = c("0263")

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
      tableOutput("progress")
      ),

    # Show a plot of the generated distribution
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
                 numericInput('decent_threshold', "decent rate threshold", value=0.0)
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
        tabPanel("Table", dataTableOutput("datatable")),
        tabPanel("Map", leafletOutput("map")),
        tabPanel("Sensors",
          h5("Check sensor serials against xml!"),
          hr(),
          fluidRow(
            column(2,
              h4("Optode"),
              selectInput("serial_optode", "Optode Serial", choices=optode_serials)
              ),
            column(5,
              selectInput('optode_T_channel', "Optode Temperature channel", choices = vchannels, selected = "v7"),
              selectInput("optode_foil", "Optode foil Batch #", choices = NULL)
              ),
            column(5,
              selectInput('optode_dphase_channel', "Optode dPhase channel", choices = vchannels, selected = "v6"),
              actionButton('optode', "Process Optode", icon=icon("life-ring"))
              )
            ),
          hr(),
          fluidRow(
            column(2,
              h4("RINKO"),
              selectInput("serial_rinko", "Rinko Serial", choices=rinko_serials)
              ),
            column(5,
              selectInput('rinko_T_channel', "RINKO Temperature channel", choices = vchannels, selected = "v5"),
              numericInput('rinko_G', label = "G Coefficent", value = 0)
              ),
            column(5,
              selectInput('rinko_O_channel', "RINKO Oxygen channel", choices = vchannels, selected = "v4"),
              numericInput('rinko_H', label = "H Coefficent", value = 1),
              actionButton('rinko', "Process RINKO", icon=icon("times-circle-o"))
              )
            ),
          hr(),
          fluidRow(
            column(2,
              h4("Licor PAR"),
              selectInput("serial_par", "PAR Serial", choices=par_serials)
              ),
            column(5,
              selectInput('par_channel', "PAR channel", choices = vchannels, selected = "v0"),
              actionButton('flag_par', "Flag all PAR for selected dip (Night)", icon=icon("moon-o"))
              ),
            column(5,
              numericInput('licor_factor', label = "Licor factor", value = 0.22345679),
              numericInput('licor_offset', label = "Licor offset", value = 3.3737),
              actionButton('licor', "Process Licor PAR", icon=icon("beer"))
              )
            ),
          hr(),
          fluidRow(
            column(2,
              h4("Fluorometer"),
              selectInput("serial_flu", "Fluorometer Serial", choices=fluorometer_serials)
              ),
            column(5,
              numericInput('par_flu_threshold', label = "Chlorophyll quenching PAR threshold", value = 1),
              actionButton('flag_flu', "Flag quenched chlorophyll fluorometry", icon=icon("ban")),
              tableOutput("chl_coef")
              ),
            column(5,
              numericInput('chl_factor', label="Chl Factor", value=1.0, step=0.01),
              numericInput('chl_offset', label="Chl Offset", value=0.0, step=0.01),
              actionButton('calc_flu', "derive Chlorophyll from flu regression", icon=icon("leaf"))
              )
            ),
          hr(),
          fluidRow(
            column(2,
              h4("Secondary CT"),
              selectInput("serial_temp", "Temperature Serial", choices=temperature_serials),
              selectInput("serial_cond", "Conductivity Serial", choices=conductivity_serials)
              ),
            column(5,
              actionButton('secondCT', "Overwrite Primary CT with secondary", icon=icon("reply-all"))
              )
            ),
          hr(),
          fluidRow(
            column(2,
              h4("Pressure"),
              selectInput("serial_prs", "Pressure Serial", choices=pressure_serials)
              )
            ),
          hr(),
          fluidRow(
            column(2,
              h4("Altimeter"),
              selectInput("serial_alt", "Altimeter Serial", choices=altimeter_serials)
              )
            )
          ),
        tabPanel("Bottles",
          rHandsontableOutput("bottles"),
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


