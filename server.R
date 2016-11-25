
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinyFiles)
library(oce)
library(data.table)
library(leaflet)
library(rhandsontable)

source("functions.R", local = T)

shinyServer(function(input, output, session) {

  # find OS disk drives
  volumes = getVolumes()
  shinyDirChoose(input, 'directory', roots=volumes, session=session, restrictions=system.file(package='base'))
  output$directory = renderText({paste0(parseDirPath(volumes, input$directory), "/")})

  # make dynamic file list for storing the CTD objects, a list of S4 objects
  profiles = reactiveValues(data = NULL)

   ## read data
  observeEvent(input$read_files, {
    filelist = list.files(parseDirPath(volumes, input$directory), full.names = F, pattern = "*.cnv")
    dir = parseDirPath(volumes, input$directory)
    d = list()
    withProgress(message = 'loading files...', value = 0, {
      for(i in filelist){
        # Increment the progress bar, and update the detail text.
        incProgress(1/length(filelist), detail = paste("loading", i))
        d[[i]] = read.ctd.sbe(paste0(dir,"/",i))
      }
    })
      # insert data into data slot
    profiles$data = d
    profiles$untrimmed = d
      # make a backup for use by revert
    profiles$original = d
      # make a summary of the positions for the map
    profiles$positions = rbindlist(lapply( profiles$data , function(x) `@`( x , metadata)[c("filename", "startTime", "station", "longitude", "latitude")]))

  })

  observeEvent(input$read_bottle, {
      # make two file lists for comparison
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
          scans = cbind(i, read.csv(paste0(dir,"/",botfile), skip = 2, header = F))
          colnames(scans) = c("profile", "fire_seq", "niskin", "dateTime", "start_scan", "end_scan")
        }
        b[[i]] = scans
      }
    })
    scans = rbindlist(b)
      # for each ctd in the list, extract the data, then for each extract convert to data.frame, then combine
    dat = rbindlist(lapply(lapply(profiles$untrimmed , function(x) `@`( x , data)), data.frame), idcol = "profile", fill = T)
      # select only the columns we want, if they are available
    avail_names = colnames(dat)
    want_names = c("depth", "salinity", "salinity2", "fluorescence", "oxygen_optode", "oxygen_RINKO")
    names = na.omit(want_names[chmatch(avail_names, want_names)])
      # for each bottle, match to scans and profile
    dat = dat[, c("profile", "scan", names), with = F]
    dat[, scan0 := scan] # extra column needed for foverlaps
    setkey(scans, profile, start_scan, end_scan)
    dat = foverlaps(dat, scans, by.x=c("profile", "scan", "scan0"), nomatch = 0)
      # calc_mean
    dat = dat[,lapply(.SD, mean), by = list(profile, fire_seq, niskin, dateTime), .SDcols = names]
    dat = cbind(dat, data.frame("bottle_sal" = NA, "bottle_O2" = NA, "bottle_Chl" = NA))
    profiles$bottle_scans = dat
  })

  observeEvent(input$read_rdata,{
    dir = parseDirPath(volumes, input$directory)
    load(paste0(dir, "/CTDQC.rdata"))
      # copy data from loaded .rdata file to correct slots
    profiles$data = session$data
    profiles$original = session$original
    profiles$positions = session$positions
  })

  ## Processes

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
    profiles$data[[input$select_profile]] = ctdDecimate(profiles$data[[input$select_profile]], p = input$bin_size)
  })

  observeEvent(input$revert,{
    profiles$data[[input$select_profile]] = profiles$original[[input$select_profile]]
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


  ## SENSORS

  observeEvent(input$optode, {
    # for working data
    optode_temperature = optode.analogtemp(
        profiles$untrimmed[[input$select_profile]]@data[[input$optode_T_channel]]
        )
    optode_Dphase = optode.analogDphase(
        profiles$untrimmed[[input$select_profile]]@data[[input$optode_dphase_channel]]
        )
    salinity = profiles$untrimmed[[input$select_profile]]@data[["salinity"]]
    depth = profiles$untrimmed[[input$select_profile]]@data[["depth"]]
    optode_oxygen = optode.phaseCalc(optode_Dphase, optode_temperature, subset(optode_coefs, batch == input$optode_foil))
    optode_oxygen = optode.correction(optode_oxygen, optode_temperature, salinity, depth)
    # now add to untrimmed
    profiles$untrimmed[[input$select_profile]] = ctdAddColumn(profiles$untrimmed[[input$select_profile]],
                                                         optode_temperature, "temperature_optode", label = "temperature",
                                                         unit = list(name=expression(degree*C), scale="Optode"))
    profiles$untrimmed[[input$select_profile]] = ctdAddColumn(profiles$untrimmed[[input$select_profile]],
                                                         optode_oxygen, "oxygen_optode",
                                                         unit = list(name = expression(mmol~m-3), scale="Optode"))
    # now subset and apply to $data, don't write subset to log
    x = subset(profiles$untrimmed[[input$select_profile]],
               scan > min(profiles$data[[input$select_profile]][["scan"]] &
               scan < max(profiles$data[[input$select_profile]][["scan"]])))
    profiles$data[[input$select_profile]] = ctdAddColumn(profiles$data[[input$select_profile]],
                                                         x[["temperature_optode"]], "temperature_optode", label = "temperature",
                                                         unit = list(name=expression(degree*C), scale="Optode"))
    profiles$data[[input$select_profile]] = ctdAddColumn(profiles$data[[input$select_profile]],
                                                         x[["oxygen_optode"]], "oxygen_optode",
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
    profiles$untrimmed[[input$select_profile]] = ctdAddColumn(profiles$original[[input$select_profile]],
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

  ## Write out

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

  ## Ui and controls
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

  ## Output

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
    # workaround for plot function not liking null data
      if(!is.null(profiles$data[[input$select_profile]])){
        # extract data from nested S4 objects (CTD)
        x1 = profiles$data[[input$select_profile]]@data[[input$x1]]
        x2 = profiles$data[[input$select_profile]]@data[[input$x2]]
        y = profiles$data[[input$select_profile]]@data[[input$y]]
        ylim = rev(range(y))
        plot(x = x2, y = y, type = "l", ylim = ylim, xlab = input$x2, ylab = input$y, col = "red")
        par(new = T)
        plot(x = x1, y = y, type = "l", ylim = ylim, axes = F, xlab = NA, ylab = NA, col = "blue")
        axis(side = 3)
        mtext(input$x1, side = 3, line = 3)
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

  output$bottles = renderRHandsontable({
    # TODO make it work for all bottles
    # scans = profiles$bottle_scans[[input$select_profile]]
    # dt = data.table(data.frame(profiles$data[[input$select_profile]]@data))
    # for(b in unique(scans$bottle)){
    #   dt[scan > scans$start[b] & scan < scans$end[b], c("bottle", "dateTime") := list(scans$bottle[b], scans$dateTime[b])]
    # }
    # potential = c("depth", "salinity", "salinity2", "fluorescence", "oxygen_optode")
    # SDcols = names(dt)[names(dt) %in% potential]
    # dt[,lapply(.SD, mean), .SDcols = SDcols, by = list(bottle, dateTime)]
      # editable table
    rhandsontable(profiles$bottle_scans, readOnly = T) %>%
      hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE) %>%
      hot_col("bottle_sal", readOnly = F) %>%
      hot_col("bottle_O2", readOnly = F) %>%
      hot_col("bottle_Chl", readOnly = F)
  })
  # output$debug = renderText({
  #   print(profiles$untrimmed[[input$select_profile]])
  # })
})


