
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
    profiles$positions = rbindlist(lapply( profiles$data , function(x) `@`( x , metadata)[c("filename", "startTime", "station", "longitude", "latitude")]))
  })

  observeEvent(input$trim,{
    profiles$data[[input$select_profile]] = ctdTrim(profiles$data[[input$select_profile]],
                                                    method = "scan", parameters = round(c(input$scan_brush$xmin, input$scan_brush$xmax)))
  })
  observeEvent(input$autotrim,{
    profiles$data[[input$select_profile]] = ctdTrim(profiles$data[[input$select_profile]], parameters = list(pmin=1))
  })

  observeEvent(input$decimate,{
    profiles$data[[input$select_profile]] = ctdDecimate(profiles$data[[input$select_profile]],
                                                        p = input$bin_size)
  })

  observeEvent(input$revert,{
    profiles$data[[input$select_profile]] = profiles$original[[input$select_profile]]
  })

  observeEvent(input$save,{
    save(profiles$data, file = "working.rdata")
  })

    # update select input when filelist changes
  observe({
    updateSelectInput(session, "select_profile", choices = names(profiles$original))
    })
  observe({
    # workaround for oce function not liking null data
      if(!is.null(profiles$data[[input$select_profile]])){
        choices = names(profiles$data[[input$select_profile]]@data)
        if("salinity" %in% choices){default_x1 = "salinity"}else{default_x1 = NULL}
        if("temperature" %in% choices){default_x2 = "temperature"}else{default_x2 = NULL}
        if("pressure" %in% choices){default_y = "pressure"}else{default_y = NULL}
        updateSelectInput(session, "x1", choices = choices, selected = default_x1)
        updateSelectInput(session, "x2", choices = choices, selected = default_x2)
        updateSelectInput(session, "y", choices = choices, selected = default_y)
      }
    })

  output$summary <- renderPrint(
    # workaround for oce function not liking null data
    if(!is.null(profiles$data[[input$select_profile]]))
    summary(profiles$data[[input$select_profile]])
    )
  output$scan_plot = renderPlot({
    # workaround for plot function not liking null data
    if(!is.null(profiles$data[[input$select_profile]]))
    plotScan(profiles$data[[input$select_profile]])
    abline(v = input$trim_scans)
    })
  output$profile_plot = renderPlot({
      if(!is.null(profiles$data[[input$select_profile]])){
        x1 = unlist(profiles$data[[input$select_profile]]@data[input$x1])
        x2 = unlist(profiles$data[[input$select_profile]]@data[input$x2])
        y = unlist(profiles$data[[input$select_profile]]@data[input$y])
        ylim = rev(range(y))
        plot(x = x1, y = y, type = "l", ylim = ylim, xlab = input$x1, ylab = input$y)
        par(new = T)
        plot(x = x2, y = y, type = "l", ylim = ylim, axes = F, xlab = NA, ylab = NA, col = "red")
        axis(side = 3)
        mtext(input$x2, side = 3, line = 3)
      }
    })
  output$TS_plot = renderPlot({
    # workaround for plot function not liking null data
    if(!is.null(profiles$data[[input$select_profile]]))
    plotTS(profiles$data[[input$select_profile]])
    })
  output$map = renderLeaflet(
    # workaround for oce function not liking null data
    if(!is.null(profiles$data[[input$select_profile]]))
    leaflet(profiles$positions) %>%
      addTiles() %>%
      addMarkers(~longitude, ~latitude, popup = ~paste("Station", station,"\r", "@", startTime)) %>%
      setView(lat = profiles$data[[input$select_profile]]@metadata$latitude,
              lng = profiles$data[[input$select_profile]]@metadata$longitude,
              zoom = 7)
  )
  output$debug = renderPrint({NULL})
})
