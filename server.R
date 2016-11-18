
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinyFiles)
library(oce)
library(data.table)
library(cefasMOS)
library(leaflet)

shinyServer(function(input, output, session) {

  volumes <- c("C:"="C:\\")
  shinyDirChoose(input, 'directory', roots=volumes, session=session, restrictions=system.file(package='base'))
  output$directory = renderText({parseDirPath(volumes, input$directory)})

    # make dynamic file list
  filelist = reactive({ list.files(parseDirPath(volumes, input$directory), full.names = F, pattern = "*.cnv") })

  profiles = reactiveValues(data = NULL)

    # read CNV files
  observeEvent(input$read_files, {
    dir = parseDirPath(volumes, input$directory)
    d = list()
    withProgress(message = 'loading files...', value = 0, {
      for(i in filelist()){
        # Increment the progress bar, and update the detail text.
        incProgress(1/length(filelist()), detail = paste("loading", i))
        d[[i]] = read.ctd.sbe(paste0(dir,"/",i))
        print(paste0(dir,i))
      }
    })
    profiles$data = d
    profiles$original = d
    profiles$positions = rbindlist(lapply( profiles$data , function(x) `@`( x , metadata)[c("startTime", "station", "longitude", "latitude")]))
  })

  observeEvent(input$trim,{
    profiles$data[[input$select_profile]] = ctdTrim(profiles$data[[input$select_profile]])
  })

  observeEvent(input$revert,{
    profiles$data[[input$select_profile]] = profiles$original[[input$select_profile]]
  })

  observeEvent(input$save,{
    save(profiles$data, file = "working.rdata")
  })

    # update select input when filelist changes
  observe({
    updateSelectInput(session, "select_profile", choices = names(profiles$data))
    })
  output$summary <- renderPrint( summary(profiles$data[[input$select_profile]]@metadata) )
  output$scan_plot = renderPlot({
    # workaround for plot function not liking null data
    if(!is.null(profiles$data[[input$select_profile]]))
    plotScan(profiles$data[[input$select_profile]])
    })
  output$profile_plot = renderPlot({
    # workaround for plot function not liking null data
    if(!is.null(profiles$data[[input$select_profile]]))
    plotProfile(profiles$data[[input$select_profile]])
    })
  output$map = renderLeaflet(
    leaflet(profiles$positions) %>%
      addTiles() %>%
      addCircles(~longitude, ~latitude)
  )
})
