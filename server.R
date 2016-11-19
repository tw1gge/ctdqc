
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

rinko_coefs = list(serial = "#0263 ARO-CAV", A = -42.34162, B = +127.6475, C = -0.3677435, D = +0.01137, E = +0.0046, F = +7.57e-05)

shinyServer(function(input, output, session) {

  volumes = getVolumes()
  shinyDirChoose(input, 'directory', roots=volumes, session=session, restrictions=system.file(package='base'))
  output$directory = renderText({paste0(parseDirPath(volumes, input$directory), "/")})

  optode_coefs = read.csv("optode_coefs.csv")

  # make dynamic file list

  profiles = reactiveValues(data = NULL)

    # read CNV files
  observeEvent(input$read_files, {
    filelist = list.files(parseDirPath(volumes, input$directory), full.names = F, pattern = "*.cnv")
    dir = parseDirPath(volumes, input$directory)
    d = list()
    withProgress(message = 'loading files...', value = 0, {
      for(i in filelist){
        # Increment the progress bar, and update the detail text.
        incProgress(1/length(filelist), detail = paste("loading", i))
        d[[i]] = read.ctd.sbe(paste0(dir,"/",i))
        # check if bl file exists and if so load it
      }
    })
    profiles$data = d
    profiles$original = d
    profiles$positions = rbindlist(lapply( profiles$data , function(x) `@`( x , metadata)[c("filename", "startTime", "station", "longitude", "latitude")]))
  })

  observeEvent(input$read_bottle, {
    filelist = list.files(parseDirPath(volumes, input$directory), full.names = F, pattern = "*.cnv")
    bottlelist = list.files(parseDirPath(volumes, input$directory), full.names = F, pattern = "*.bl")
    dir = parseDirPath(volumes, input$directory)
    b = list()
    withProgress(message = 'loading files...', value = 0, {
      for(i in filelist){
        # Increment the progress bar, and update the detail text.
        incProgress(1/length(filelist), detail = paste("loading", i))
        # matching
        botfile = gsub(".cnv", ".bl", i)
        if(botfile %in% bottlelist){
          dat = read.csv(paste0(dir,"/",botfile), skip = 2, header = F)
          colnames(dat) = c("seq", "bottle", "dateTime", "start", "end")
        }else{dat = "No bottles"}
        b[[i]] = dat
      }
    })
    profiles$bottle_scans = b
  })

  observeEvent(input$read_rdata,{
    dir = parseDirPath(volumes, input$directory)
    load(paste0(dir, "/CTDQC.rdata"))
    profiles$data = session$data
    profiles$original = session$original
    profiles$positions = session$positions
  })

  observeEvent(input$pumped,{
    profiles$data[[input$select_profile]] = subset(profiles$data[[input$select_profile]], pumpStatus == 1)
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

  ## SENSORS

  observeEvent(input$optode, {
    optode_temperature = optode.analogtemp(
        profiles$data[[input$select_profile]]@data[[input$optode_T_channel]]
        )
    optode_Dphase = optode.analogDphase(
        profiles$data[[input$select_profile]]@data[[input$optode_dphase_channel]]
        )
    optode_oxygen = optode.phaseCalc(optode_Dphase, optode_temperature, subset(optode_coefs, batch == input$optode_foil))
    # TODO salinity correction
    profiles$data[[input$select_profile]] = ctdAddColumn(profiles$data[[input$select_profile]],
                                                         optode_temperature, "temperature_optode", label = "temperature",
                                                         unit = list(name=expression(degree*C), scale="Optode"))
    profiles$data[[input$select_profile]] = ctdAddColumn(profiles$data[[input$select_profile]],
                                                         optode_oxygen, "oxygen_optode",
                                                         unit = list(name = expression(mmol~m-3), scale="Optode"))
    processingLog(profiles$data[[input$select_profile]]) = paste("Optode processed with foil batch", input$optode_foil)

  })

  observeEvent(input$rinko, {
    rinko_temperature = rinko_temp(
        profiles$data[[input$select_profile]]@data[[input$rinko_T_channel]]
        )
    pressure =  unlist(profiles$data[[input$select_profile]]@data["pressure"])
    rinko_oxygen = rinko_o2(
        profiles$data[[input$select_profile]]@data[[input$rinko_O_channel]],
        rinko_temperature,
        S = profiles$data[[input$select_profile]]@data[["salinity"]],
        oC = rinko_coefs,
        G = input$rinko_G,
        H = input$rinko_H
        )
    profiles$data[[input$select_profile]] = ctdAddColumn(profiles$data[[input$select_profile]],
                                                         rinko_temperature, "temperature_RINKO", label = "temperature",
                                                         unit = list(name=expression(degree*C), scale="RINKO"))
    profiles$data[[input$select_profile]] = ctdAddColumn(profiles$data[[input$select_profile]],
                                                         rinko_oxygen, "oxygen_RINKO", label = "oxygen",
                                                         unit = list(name = expression(mmol~m-3), scale="RINKO"))
    processingLog(profiles$data[[input$select_profile]]) = paste("RINKO processed with coefs", rinko_coefs[["serial"]],
                                                                 ",G =", input$rinko_G,
                                                                 ",H =", input$rinko_H)

  })

  observeEvent(input$licor, {
    licor_par = par_from_voltage(
        profiles$data[[input$select_profile]]@data[[input$par_channel]],
        input$licor_factor, input$licor_offset)

    profiles$data[[input$select_profile]] = ctdAddColumn(profiles$data[[input$select_profile]],
                                                         licor_par, "par", label = "par",
                                                         unit = list(name = expression(uE), scale = "PAR/Irradiance, Cefas Licor"))
    processingLog(profiles$data[[input$select_profile]]) = paste("PAR processed with factor =", input$licor_factor, ",offset =", input$licor_offset)
  })

  observeEvent(input$write_rdata,{
    dir = parseDirPath(volumes, input$directory)
    session = profiles
    save(session, file = paste0(dir, "/CTDQC.rdata"))
  })

  observeEvent(input$write_csv,{
    dir = parseDirPath(volumes, input$directory)
    withProgress(message = 'writing files...', value = 0, {
      for(p in names(profiles$data)){
        incProgress(1/length(profiles$data), detail = paste("writing", p))
        write.ctd(profiles$data[[p]], file = paste0(p, ".csv"))
      }
    })
  })

  observeEvent(input$apply_flag,{
    warning("not implemented")
  })

  observeEvent(input$apply_factor,{
    raw = profiles$data[[input$select_profile]]@data[[input$x1]]
    mod = (raw * input$factor) + input$offset
    profiles$data[[input$select_profile]]@data[[input$x1]] = mod
    processingLog(profiles$data[[input$select_profile]]) = paste(input$x1,
                                                                 ",adjusted with factor", input$factor,
                                                                 ", offset", input$offset)
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
        updateSelectInput(session, "select_factor", choices = choices)
        updateSelectInput(session, "select_flag", choices = choices)
      }
    })
  observe({
    updateSelectInput(session, "optode_foil", choices = unique(optode_coefs$batch), selected = "1707")
    })

  output$summary <- renderPrint(
    # workaround for oce function not liking null data
    if(!is.null(profiles$data[[input$select_profile]]))
    summary(profiles$data[[input$select_profile]])
    )
  output$xml <- renderPrint(
    # workaround for oce function not liking null data
      if(!is.null(profiles$data[[input$select_profile]])){
        profiles$data[[input$select_profile]]@metadata["header"]
      }
    )
  output$scan_plot = renderPlot({
    # workaround for plot function not liking null data
    if(!is.null(profiles$data[[input$select_profile]]))
    plotScan(profiles$data[[input$select_profile]])
    abline(v = input$trim_scans)
    })
  output$profile_plot = renderPlot({
      if(!is.null(profiles$data[[input$select_profile]])){
        x1 = profiles$data[[input$select_profile]]@data[[input$x1]]
        x2 = profiles$data[[input$select_profile]]@data[[input$x2]]
        y = profiles$data[[input$select_profile]]@data[[input$y]]
        ylim = rev(range(y))
        plot(x = x2, y = y, type = "l", ylim = ylim, xlab = input$x1, ylab = input$y, col = "blue")
        par(new = T)
        plot(x = x1, y = y, type = "l", ylim = ylim, axes = F, xlab = NA, ylab = NA, col = "red")
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
  output$datatable = renderDataTable({
    data.frame(profiles$data[[input$select_profile]]@data)
  })
  output$bottles = renderPrint({
    scans = profiles$bottle_scans[[input$select_profile]]
    dt = data.table(data.frame(profiles$data[[input$select_profile]]@data))
    for(b in unique(scans$bottle)){
      dt[scan > scans$start[b] & scan < scans$end[b], c("bottle", "dateTime") := list(scans$bottle[b], scans$dateTime[b])]
    }
    potential = c("depth", "salinity", "salinity2", "fluorescence", "oxygen_optode")
    SDcols = names(dt)[names(dt) %in% potential]
    dt[,lapply(.SD, mean), .SDcols = SDcols, by = list(bottle, dateTime)]
  })
  output$debug = renderPrint({ NULL })
})

optode.analogtemp <- function(v){
  return( ((v * 45) / 5) - 5 )
}
optode.analogDphase <- function(v){
  return( 10 + (v / 5) * 60 )
}
optode.phaseCalc <- function(DPhase, Temp, coefs){
  # for mkl optodes 3830 & 3835
    with(coefs, {
      print(paste("using foil batch coefs", batch[1]))
      (C0[1]+C0[2]*Temp+C0[3]*Temp^2+C0[4]*Temp^3) +
      (C1[1]+C1[2]*Temp+C1[3]*Temp^2+C1[4]*Temp^3) *
      DPhase+(C2[1]+C2[2]*Temp+C2[3]*Temp^2+C2[4]*Temp^3) *
      DPhase^2+(C3[1]+C3[2]*Temp+C3[3]*Temp^2+C3[4]*Temp^3) *
      DPhase^3+(C4[1]+C4[2]*Temp+C4[3]*Temp^2+C4[4]*Temp^3) *
      DPhase^4
    })
}
par_from_voltage <- function(x, factor, offset){
    return(factor * exp(offset * x))
}

