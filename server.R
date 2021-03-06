library(shiny,quietly=T)
library(shinyFiles,quietly=T)
library(oce, quietly=T)
library(data.table, quietly=T)
library(leaflet, quietly=T)
library(rhandsontable, quietly=T)
library(xml2, quietly=T)

source("functions.R", local = T)
CTDQC_version = "1.6"
editable_metadata = c("id", "title", "summary", "processing_level", "comment", "acknowledgment", "licence", "project", "creator", "creator_email")
sensor_metadata = fread("sensor_table.csv")
ctd_columns = list(
  PAR = list(name="par/sat/log", unit=list(expression(), scale="umol photon m-2"))
  )

shinyServer(function(input, output, session) {

  # find OS disk drives
  volumes = getVolumesFast()
  shinyDirChoose(input, 'directory', roots=volumes, session=session, restrictions=system.file(package='base'))
  output$directory = renderText({paste0(parseDirPath(volumes, input$directory), "/")})
  # make dynamic file list for storing the CTD objects, a list of S4 objects
  profiles = reactiveValues(data = NULL, bottles = NA)

  observeEvent(input$read_files, {
   ## read data
    if(is.null(input$directory)){ # stop crashing when you missclick
      showNotification("no folder selected!", type="error")
      return(NULL)
      }
    filelist = list.files(parseDirPath(volumes, input$directory), full.names = F, pattern = "*.cnv")
    if(length(filelist) == 0){ # don't crash if folder is empty
      showNotification("no .cnv files found in directory", type="error")
      return(NULL)
      }
    dir = parseDirPath(volumes, input$directory)
    d = list()
    m = list()
    withProgress(message = 'loading files...', value = 0, {
      for(i in filelist){
        # Increment the progress bar, and update the detail text.
        incProgress(1/length(filelist), detail = paste("loading", i))
        d[[i]] = read.ctd.sbe(paste0(dir,"/",i), columns=ctd_columns)
        m[[i]] = parse_sbe_xml(d[[i]])
      }
    })
      # check if filter has been applied
    headers = extract.oce.metadata(d, "header")
    filtered = stringr::str_count(headers, "filter_low_pass")
    if(filtered < length(d)){showNotification("pressure filter has not been applied for all profiles!", duration=NULL, type="warning") }

      # insert data into data slot
    profiles$data = d
    profiles$metadata = m
    profiles$untrimmed = d
      # make a backup for use by revert
    profiles$original = d
      # make a summary of the positions for the map
    profiles$positions = extract.oce.metadata(profiles$data, c("filename", "startTime", "station", "longitude", "latitude", "cruise"))
    profiles$global_metadata = netcdf.metadata(profiles$data, profiles$positions)
    profiles$global_metadata_default = profiles$global_metadata
    if(length(unique(profiles$positions$cruise)) > 1){showNotification("Cruise ID differ between cnv files!", duration=NULL, type="warning")}
    if(length(unique(m)) != 1){ showNotification("xml header differs between files", type="warning", duration=NULL) }
  })

  observeEvent(input$read_bottle, {
      # make two file lists for comparison
    if(is.null(input$directory) | is.null(profiles$data)){ return(NULL) } # stop crashing when you missclick
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
    scans = data.table(rbindlist(b))
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
    dat = foverlaps(dat, na.omit(scans), by.x=c("profile", "scan", "scan0"), nomatch = 0)
      # calc_mean
    dat = dat[,lapply(.SD, mean), by = list(profile, fire_seq, niskin, dateTime), .SDcols = names]
    dat = cbind(dat, data.frame("bottle_sal" = 0, "bottle_O2" = 0, "bottle_Chl" = 0))
    profiles$bottles = dat
  })

  observeEvent(input$read_rdata,{
    dir = parseDirPath(volumes, input$directory)
    load(paste0(dir, "/CTDQC.rdata"))
      # copy data from loaded .rdata file to correct slots
    profiles$data = session$data
    profiles$metadata = session$metadata
    profiles$original = session$original
    profiles$untrimmed = session$untrimmed
    profiles$positions = session$positions
    profiles$bottles = session$bottles
    profiles$global_metadata = session$global_metadata
    profiles$global_metadata_default = profiles$global_metadata
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

  observeEvent(input$remove_pressure_inversions,{
    # calculates smoothed decent rate from pressure (as per SBE data processing)
    prs = profiles$data[[input$select_profile]]@data[["pressure"]]
    decent_diff = c(0, diff(prs))
    decent_calc = zoo::rollmean(decent_diff, 48, fill=NA, align="right")*24 # for 9plus 2 second window
    # decent_calc = zoo::rollmean(decent_diff, 8, fill=NA, align="right")*4 # for 19plus 2 second window

    profiles$data[[input$select_profile]]@data = lapply(profiles$data[[input$select_profile]]@data, function(x) {
      x[decent_calc < input$decent_threshold] = NA
      return(x)
    })
    processingLog(profiles$data[[input$select_profile]]) = paste("Pressure inversions removed, minimum speed = ", input$decent_threshold,"m/s")
  })

  observeEvent(input$decimate,{
    profiles$data = lapply(profiles$data, function(x){
        # find data which is all NA
      for(param in names(x@data)){
        if(all(is.na(x[[param]]))){
          showNotification(paste("no data for", param), type="warning")
          x[[param]] = NULL
        }
      }
      return(ctdDecimate(x, p=input$bin_size))
    })
  })

  observeEvent(input$revert,{
    profiles$data[[input$select_profile]] = profiles$original[[input$select_profile]]
  })

  observeEvent(input$apply_flag,{
      # find the range of y to flag x on, based on plot brush
    y_select = round(c(input$flag_brush$ymin, input$flag_brush$ymax), 3)
    x_select = round(c(input$flag_brush$xmin, input$flag_brush$xmax), 3)
      # extract the y variable
    y = profiles$data[[input$select_profile]]@data[[input$y]]
    x = profiles$data[[input$select_profile]]@data[[input$x1]]
      # replace with NA sections of x which match the selected y
    profiles$data[[input$select_profile]]@data[[input$x1]][y %between% y_select & x %between% x_select] = NA
      # write details to ctd log
    log = paste("flag applied to", input$x1, "for",
                input$y, "between", y_select[1], "-", y_select[2],
                "and", input$x1, "between", x_select[1], "-", x_select[2])
    processingLog(profiles$data[[input$select_profile]]) = log
  })

  observeEvent(input$apply_factor,{
    # lapply though and apply to all dips
    profiles$data = lapply(profiles$data, function(x) {
      raw = x@data[[input$x1]]
      mod = (raw * input$factor) + input$offset
      x@data[[input$x1]] = mod
      log = paste(input$x1, ",adjusted with factor", input$factor, ", offset", input$offset)
      processingLog(x) = log
      return(x)
    })
  })

  observeEvent(input$calc_flu,{
    # lapply though and apply to all dips
    profiles$data = lapply(profiles$data, function(x) {
      raw = x@data[["fluorescence"]]
      mod = (raw * input$chl_factor) + input$chl_offset
      x@data[["chlorophyll"]] = mod
      log = paste("Chlorophyll derived with factor", input$chl_factor, ", offset", input$chl_offset)
      processingLog(x) = log
      return(x)
    })
  })

  observeEvent(input$mark_complete_QC2,{
    processingLog(profiles$data[[input$select_profile]]) = paste("Manual (QC2) Complete")
  })

  observeEvent(input$mark_complete_all,{
    processingLog(profiles$data[[input$select_profile]]) = paste("All QC Complete")
  })

  ## SENSORS

  observeEvent(input$optode, {
    # first calculate for all data (using untrimmed)
    for(i in names(profiles$data)){
      optode_temperature = optode.analogtemp(profiles$untrimmed[[i]]@data[[input$optode_T_channel]])
      optode_Dphase = optode.analogDphase(profiles$untrimmed[[i]]@data[[input$optode_dphase_channel]])
      salinity = profiles$untrimmed[[i]]@data[["salinity"]]
      depth = profiles$untrimmed[[i]]@data[["depth"]]
      optode_oxygen = optode.phaseCalc(optode_Dphase, optode_temperature, subset(optode_coefs, batch == input$optode_foil))
      optode_oxygen = optode.correction(optode_oxygen, optode_temperature, salinity, depth)
      # add to untrimmed

      # profiles$untrimmed[[i]] = oceSetData(profiles$untrimmed[[i]], "par", licor_par,
      #                                      unit = list(unit=expression(umol~s-1~m-2), scale = "PAR/Irradiance, Cefas Licor PAR"))


      profiles$untrimmed[[i]] = oceSetData(profiles$untrimmed[[i]], "temperature_optode", optode_temperature,
                                             unit = list(unit=expression(degree*C), scale="Optode"))
      profiles$untrimmed[[i]] = oceSetData(profiles$untrimmed[[i]], "oxygen_optode", optode_oxygen,
                                             unit = list(unit=expression(mmol~m-3), scale="Optode"))
      # now subset and apply to $data, don't write subset to log
      x = subset(profiles$untrimmed[[i]],
                 scan >= min(profiles$data[[i]][["scan"]], na.rm=T) &
                 scan <= max(profiles$data[[i]][["scan"]], na.rm=T))
      profiles$data[[i]] = oceSetData(profiles$data[[i]], "temperature_optode", x[["temperature_optode"]],
                                        unit = list(unit=expression(degree*C), scale="Optode"))
      profiles$data[[i]] = oceSetData(profiles$data[[i]], "oxygen_optode", x[["oxygen_optode"]],
                                        unit = list(unit=expression(mmol~m-3), scale="Optode"))
      processingLog(profiles$data[[i]]) = paste("Optode processed with foil batch", input$optode_foil)
    }
  })

  observeEvent(input$rinko, {
    for(i in names(profiles$data)){
      rinko_temperature = rinko_temp(profiles$untrimmed[[i]]@data[[input$rinko_T_channel]])
      pressure =  unlist(profiles$untrimmed[[i]]@data["pressure"])
      rinko_oxygen = rinko_o2(profiles$untrimmed[[i]]@data[[input$rinko_O_channel]],
                              rinko_temperature,
                              S = profiles$untrimmed[[i]]@data[["salinity"]],
                              oC = rinko_coefs,
                              G = input$rinko_G, H = input$rinko_H)
      # add to untrimmed
      profiles$untrimmed[[i]] = oceSetData(profiles$untrimmed[[i]], "temperature_RINKO", rinko_temperature,
                                             unit = list(unit=expression(degree*C), scale="RINKO"))
      profiles$untrimmed[[i]] = oceSetData(profiles$untrimmed[[i]], "oxygen_RINKO", rinko_oxygen,
                                             unit = list(unit=expression(mmol~m-3), scale="RINKO"))
      # now subset and apply to $data, don't write subset to log
      x = subset(profiles$untrimmed[[i]],
                 scan >= min(profiles$data[[i]][["scan"]], na.rm=T) &
                 scan <= max(profiles$data[[i]][["scan"]], na.rm=T))

      profiles$data[[i]] = oceSetData(profiles$data[[i]], "temperature_RINKO", x[["temperature_RINKO"]],
                                        unit = list(unit=expression(degree*C), scale="RINKO"))
      profiles$data[[i]] = oceSetData(profiles$data[[i]], "oxygen_RINKO", x[["oxygen_RINKO"]],
                                        unit = list(unit=expression(mmol~m-3), scale="RINKO"))
      processingLog(profiles$data[[i]]) = paste("RINKO processed with coefs",
                                                rinko_coefs[["serial"]],
                                                ",G =", input$rinko_G,
                                                ",H =", input$rinko_H)
    }
  })

  observeEvent(input$licor, {
      # loop and apply to all dips
    for(i in names(profiles$data)){
      # calculate using untrimmed data
      licor_par = par_from_voltage(profiles$untrimmed[[i]]@data[[input$par_channel]], input$licor_factor, input$licor_offset)
      profiles$untrimmed[[i]] = oceSetData(profiles$untrimmed[[i]], "par", licor_par,
                                           unit = list(unit=expression(umol~s-1~m-2), scale = "PAR/Irradiance, Cefas Licor PAR"))
      # now subset and add to data
      x = subset(profiles$untrimmed[[i]],
                 scan >= min(profiles$data[[i]][["scan"]], na.rm=T) &
                 scan <= max(profiles$data[[i]][["scan"]], na.rm=T))

      log = paste("PAR processed with factor =", input$licor_factor, ",offset =", input$licor_offset)
      profiles$data[[i]] = oceSetData(profiles$data[[i]], "par", x[["par"]],
                                      unit = list(unit=expression(umol~s-1~m-2), scale = "PAR/Irradiance, Cefas Licor PAR"),
                                      note = log)
    }
  })

  observeEvent(input$flag_par, {
    if("par" %in% names(profiles$data[[input$select_profile]]@data)){
      profiles$data[[input$select_profile]][["par"]] = NA
      log = paste("All PAR removed (Night)")
      processingLog(profiles$data[[input$select_profile]]) = log
    }
    else{
      showNotification("PAR has not been calculated, PAR not flagged", type="error", duration=NULL)
    }
  })


  observeEvent(input$flag_flu, {
    if("par" %in% names(profiles$data[[input$select_profile]]@data)){
      profiles$data[[input$select_profile]][["fluorescence"]] [profiles$data[[input$select_profile]][["par"]] > input$par_flu_threshold] = NA
      profiles$untrimmed[[input$select_profile]][["fluorescence"]] [profiles$untrimmed[[input$select_profile]][["par"]] > input$par_flu_threshold] = NA
      log = paste("fluorometry values where PAR >", input$par_flu_threshold, " have been flagged")
      processingLog(profiles$data[[input$select_profile]]) = log
    }
    else{
      showNotification("PAR has not been calculated, Fluorometry not flagged", type="error", duration=NULL)
    }
  })

  observeEvent(input$secondCT, {
    try({
      profiles$data[[input$select_profile]][["temperature"]] = profiles$data[[input$select_profile]][["temperature2"]]
      profiles$data[[input$select_profile]][["conductivity"]] = profiles$data[[input$select_profile]][["conductivity2"]]
      profiles$data[[input$select_profile]][["salinity"]] = profiles$data[[input$select_profile]][["salinity2"]]
      profiles$data[[input$select_profile]][["salinityDifference"]] = 0
      processingLog(profiles$data[[input$select_profile]]) = paste("Primary CT data replaced with that from secondary")
      showNotification("Primary CT data replaced with that from secondary")
    })
  })

  ## Write out

  observeEvent(input$write_rdata,{
    dir = parseDirPath(volumes, input$directory)
    session = list()
    session$data = profiles$data
    session$untrimmed = profiles$untrimmed
    session$original = profiles$original
    session$positions = profiles$positions
    session$global_metadata = profiles$global_metadata
    session$metadata = profiles$metadata
    if(!is.null(profiles$bottles)){
      session$bottles = profiles$bottles
    }
    if(!is.null(input$bottles)){
      session$bottles = hot_to_r(input$bottles)
    }
    save(session, file = paste0(dir, "/CTDQC.rdata"))
  })

  observeEvent(input$write_csv,{
    dir = parseDirPath(volumes, input$directory)
    withProgress(message = 'writing files...', value = 0, {
      data = lapply(profiles$data , function(x) `@`( x , data))
      metadata = lapply(profiles$data , function(x) `@`( x , metadata))
      log = lapply(profiles$data , function(x) `@`( x , processingLog))
      for(p in names(profiles$data)){
        incProgress(1/length(profiles$data), detail = paste("writing", p))
        d = as.data.table(data[[p]])
        if(!"latitude" %in% colnames(d)){
          d$latitude = metadata[[p]]$latitude
          d$longitude = metadata[[p]]$longitude
        }
        d$startTime = metadata[[p]]$startTime
        d$station = metadata[[p]]$station
        d[, dateTime := startTime + time]
        if(any(grepl("ctdDecimate", log[[p]]$value))){
          d = d[,-c("scan", "time", "dateTime"), with=F]
        }
        write.csv(d, file = paste0(dir, "/", p, ".csv"))
        write.csv(as.data.table(log[[p]]), file = paste0(dir, "/", p, ".log"))
      }
    })
  })

  observeEvent(input$make_netcdf, {
    sensor_metadata[variable == "temperature", serial := input$serial_temp]
    sensor_metadata[variable == "conductivity", serial := input$serial_cond]
    sensor_metadata[variable == "salinity", serial := input$serial_prs]
    sensor_metadata[variable == "pressure", serial := input$serial_prs]
    sensor_metadata[variable == "altimeter", serial := input$serial_alt]
    sensor_metadata[variable == "turbidity", serial := input$serial_ftu]
    sensor_metadata[variable == "fluorescence", serial := input$serial_flu]
    sensor_metadata[variable == "oxygen_RINKO", serial := input$serial_rinko]
    sensor_metadata[variable == "oxygen_optode", serial := input$serial_optode]
    sensor_metadata[variable == "par", serial := input$serial_par]
    write.ctd.netcdf(profiles, sensor_metadata)
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

  output$summary <- renderPrint({
    # workaround for oce function not liking null data
    if(!is.null(profiles$data[[input$select_profile]]))
    summary(profiles$data[[input$select_profile]])
    })
  output$xml <- renderText({
    # workaround for oce function not liking null data
    if(!is.null(profiles$data[[input$select_profile]])){
      paste(profiles$data[[input$select_profile]]@metadata$header, collapse="\n")
    }
    })
  output$scan_plot = renderPlot({
      # check if there is data, give warning if not
    validate(need(!is.null(profiles$data[[input$select_profile]]), "Data not loaded"))
    validate(
      need(
        any(grepl("filter_low_pass", profiles$data[[input$select_profile]]@metadata$header)),
           "Pressure filter has not been applied for this profile")
      )
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
        ylim = rev(range(y, na.rm=T))
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
  output$map = renderLeaflet({
    # workaround for oce function not liking null data
    if(!is.null(profiles$data[[input$select_profile]]))
    leaflet(profiles$positions) %>%
      addTiles() %>%
      addMarkers(~longitude, ~latitude, popup = ~paste("Station", station,"\r", "@", startTime)) %>%
      setView(lat = profiles$data[[input$select_profile]]@metadata$latitude,
              lng = profiles$data[[input$select_profile]]@metadata$longitude,
              zoom = 7)
  })
  output$datatable = renderDataTable({
    validate(need(profiles$untrimmed, "data not loaded"))
    data.frame(profiles$untrimmed[[input$select_profile]]@data)
  })
  output$bottles = renderRHandsontable({
      # editable table

    validate(need(profiles$bottles, "bottle file not loaded"))

    avail_sal_cols = colnames(profiles$bottles)[chmatch(c("salinity", "salinity2", "bottle_sal"), colnames(profiles$bottles))]
    rhandsontable(profiles$bottles, readOnly = T, digits = 6, highlightRow = T) %>%
      hot_col(c("bottle_sal", "bottle_O2", "bottle_Chl"), readOnly = F) %>%
      hot_col(na.omit(avail_sal_cols), format = "0.0000") %>%
      hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE)
  })

  output$bottle_plot = renderPlot({
    if(input$Plot_bottle_select == "Salinity"){
      dat = hot_to_r(input$bottles)[bottle_sal != 0]
      m = lm(data = dat, bottle_sal ~ salinity)
      par(mfrow = c(1, 2))
      plot(dat$salinity, dat$bottle_sal,
           col = "blue", ylab = "Salinity bottle", xlab = "CT", main = "CT vs Niskin, Primary CT",
           sub = paste0("y = ", round(coef(m)[1], 3), " + ", round(coef(m)[2], 3),"x"))
      abline(m)
      hist(m$residuals, main = "Residuals (Primary CT)")
      profiles$bottle_coef[["salinity"]] = list(var = "salinity", slope = coef(m)[2], intercept = coef(m)[1])
    }
    if(input$Plot_bottle_select == "Oxygen Optode"){
      dat = hot_to_r(input$bottles)[bottle_O2 != 0]
      m = lm(data = dat, bottle_O2 ~ oxygen_optode)
      par(mfrow = c(1, 2))
      plot(dat$oxygen_optode, dat$bottle_O2,
           col = "green", ylab = "Winkler", xlab = "Optode", main = "Optode vs Winkler",
           sub = paste0("y = ", round(coef(m)[1], 3), " + ", round(coef(m)[2], 3),"x"))
      abline(m)
      hist(m$residuals, main = "Residuals")
      profiles$bottle_coef[["oxygen_optode"]] = list(var = "oxygen_optode", slope = coef(m)[2], intercept = coef(m)[1])
    }
    if(input$Plot_bottle_select == "Oxygen RINKO"){
      dat = hot_to_r(input$bottles)[bottle_O2 != 0]
      m = lm(data = dat, bottle_O2 ~ oxygen_RINKO)
      par(mfrow = c(1, 2))
      plot(dat$oxygen_RINKO, dat$bottle_O2,
           col = "green", ylab = "Winkler", xlab = "RINKO", main = "RINKO vs Winkler",
           sub = paste0("y = ", round(coef(m)[1], 3), " + ", round(coef(m)[2], 3),"x"))
      abline(m)
      hist(m$residuals, main = "Residuals")
      profiles$bottle_coef[["oxygen_RINKO"]] = list(var = "oxygen_RINKO", slope = coef(m)[2], intercept = coef(m)[1])
    }
    if(input$Plot_bottle_select == "Chlorophyll"){
      dat = hot_to_r(input$bottles)[bottle_Chl != 0]
      m = lm(data = dat, bottle_Chl ~ fluorescence)
      par(mfrow = c(1, 2))
      plot(dat$fluorescence, dat$bottle_Chl,
           col = "green", ylab = "Chlorophyll (mg/l)", xlab = "Fluorescence", main = "Seapoint Vs Chlorophyll",
           sub = paste0("y = ", round(coef(m)[1], 3), " + ", round(coef(m)[2], 3),"x"))
      abline(m)
      hist(m$residuals, main = "Residuals")
      profiles$bottle_coef[["fluorescence"]] = list(var = "fluorescence", slope = coef(m)[2], intercept = coef(m)[1])
    }
  })

  output$bottle_coef = renderTable({
    print(rbindlist(profiles$bottle_coef))
    rbindlist(profiles$bottle_coef)
  })
  output$chl_coef = renderTable({
    data.frame(profiles$bottle_coef[["fluorescence"]])
  })
  output$progress = renderTable({
    # fetch processing log then grep for string to check progress
    log = lapply(profiles$data , function(x) `@`( x , processingLog))
    sensor = grepl("processed", log, ignore.case=T)
    trim = grepl("ctdTrim", log, ignore.case=T)
    QC2 = grepl("QC2", log, ignore.case=T)
    done = grepl("QC Complete", log, ignore.case=T)
    data.frame(
      "dip" = names(profiles$data),
      "trim" = trim,
      "sensor" = sensor,
      "QC2" = QC2,
      "done" = done
      )
  })
  output$edit_metadata = renderUI({
    validate(need(profiles$global_metadata_default, "data not loaded"))

    lapply(editable_metadata, function(i){
      # generate UI dynamically
      id = paste0("netcdf-", i)
      default_value = profiles$global_metadata_default[[i]] # use copy to avoid instant-writeback
      textInput(id, i, value=default_value)
    })
  })
  output$metadata = renderText({
    validate(need(profiles$global_metadata, ""))
    for(i in editable_metadata){
      id = paste0("netcdf-", i)
      profiles$global_metadata[[i]] = input[[id]]
    }
    paste(names(profiles$global_metadata), "; ", profiles$global_metadata, collapse="\n")
  })

})


