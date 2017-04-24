library(shiny,quietly=T)
library(shinyFiles,quietly=T)
library(oce, quietly=T)
library(data.table, quietly=T)
library(leaflet, quietly=T)
library(rhandsontable, quietly=T)
library(xml2, quietly=T)

source("functions.R", local = T)
CTDQC_version = 1.2
editable_metadata = c("id", "title", "summary", "processing_level", "comment", "acknowledgment", "licence", "project", "creator", "creator_email")
sensor_metadata = fread("sensor_table.csv")

shinyServer(function(input, output, session) {

  # find OS disk drives
  volumes = getVolumesFast()
  shinyDirChoose(input, 'directory', roots=volumes, session=session, restrictions=system.file(package='base'))
  output$directory = renderText({paste0(parseDirPath(volumes, input$directory), "/")})

  # make dynamic file list for storing the CTD objects, a list of S4 objects
  profiles = reactiveValues(data = NULL, bottles = NA)

   ## read data
  observeEvent(input$read_files, {
    if(is.null(input$directory)){ return(NULL) } # stop crashing when you missclick
    filelist = list.files(parseDirPath(volumes, input$directory), full.names = F, pattern = "*.cnv")
    dir = parseDirPath(volumes, input$directory)
    d = list()
    m = list()
    withProgress(message = 'loading files...', value = 0, {
      for(i in filelist){
        # Increment the progress bar, and update the detail text.
        incProgress(1/length(filelist), detail = paste("loading", i))
        d[[i]] = read.ctd.sbe(paste0(dir,"/",i))
        m[[i]] = parse_sbe_xml(d[[i]])
      }
    })
      # check if filter has been applied
    headers = extract.metadata(d, "header")
    filtered = stringr::str_count(headers, "filter_low_pass_A_vars = prDM")



      # insert data into data slot
    profiles$data = d
    profiles$metadata = m
    profiles$untrimmed = d
      # make a backup for use by revert
    profiles$original = d
      # make a summary of the positions for the map
    profiles$positions = extract.metadata(profiles$data, c("filename", "startTime", "station", "longitude", "latitude", "cruise"))
    profiles$global_metadata = netcdf.metadata(profiles$data, profiles$positions)
    profiles$global_metadata_default = profiles$global_metadata
    if(length(unique(profiles$positions$cruise)) > 1){warning("WARNING - Cruise ID differ between cnv files!")}
    if(filtered < length(d)){warning("WARNING - pressure filter has not been applied for all profiles!") }
    if(length(unique(m)) != 1){ warning("WARNING - xml header differs between files") }
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
    dat = cbind(dat, data.frame("bottle_sal" = 0, "bottle_O2" = 0, "bottle_Chl" = 0))
    profiles$bottles = dat
  })

  observeEvent(input$read_rdata,{
    dir = parseDirPath(volumes, input$directory)
    load(paste0(dir, "/CTDQC.rdata"))
      # copy data from loaded .rdata file to correct slots
    profiles$data = session$data
    profiles$original = session$original
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

  observeEvent(input$decimate,{
    profiles$data = lapply(profiles$data, ctdDecimate, p=input$bin_size)
  })

  observeEvent(input$revert,{
    profiles$data[[input$select_profile]] = profiles$original[[input$select_profile]]
  })

  observeEvent(input$apply_flag,{
      # find the range of y to flag x on, based on plot brush
    y_select = round(c(input$flag_brush$ymin, input$flag_brush$ymax), 3)
      # extract the y variable
    y = profiles$data[[input$select_profile]]@data[[input$y]]
      # replace with NA sections of x which match the selected y
    profiles$data[[input$select_profile]]@data[[input$x1]][y %between% y_select] = NA
      # write details to ctd log
    log = paste("flag applied to", input$x1, "for", input$y, "between", y_select[1], "and", y_select[2])
    processingLog(profiles$data[[input$select_profile]]) = log
  })

  observeEvent(input$apply_factor,{
    raw = profiles$data[[input$select_profile]]@data[[input$x1]]
    mod = (raw * input$factor) + input$offset
    profiles$data[[input$select_profile]]@data[[input$x1]] = mod
    processingLog(profiles$data[[input$select_profile]]) = paste(input$x1,
                                                                 ",adjusted with factor", input$factor,
                                                                 ", offset", input$offset)
  })

  observeEvent(input$mark_complete,{
    processingLog(profiles$data[[input$select_profile]]) = paste("QC Complete")
  })

  ## SENSORS

  observeEvent(input$optode, {
    # first calculate for all data (using untrimmed)
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
    # add to untrimmed
    profiles$untrimmed[[input$select_profile]] = ctdAddColumn(profiles$untrimmed[[input$select_profile]],
                                                         optode_temperature, "temperature_optode", label = "temperature",
                                                         unit = list(name=expression(degree*C), scale="Optode"))
    profiles$untrimmed[[input$select_profile]] = ctdAddColumn(profiles$untrimmed[[input$select_profile]],
                                                         optode_oxygen, "oxygen_optode",
                                                         unit = list(name = expression(mmol~m-3), scale="Optode"))
    # now subset and apply to $data, don't write subset to log
    x = subset(profiles$untrimmed[[input$select_profile]],
               scan >= min(profiles$data[[input$select_profile]][["scan"]]) &
               scan <= max(profiles$data[[input$select_profile]][["scan"]]))
    profiles$data[[input$select_profile]] = ctdAddColumn(profiles$data[[input$select_profile]],
                                                         x[["temperature_optode"]], "temperature_optode", label = "temperature",
                                                         unit = list(name=expression(degree*C), scale="Optode"))
    profiles$data[[input$select_profile]] = ctdAddColumn(profiles$data[[input$select_profile]],
                                                         x[["oxygen_optode"]], "oxygen_optode",
                                                         unit = list(name = expression(mmol~m-3), scale="Optode"))
    processingLog(profiles$data[[input$select_profile]]) = paste("Optode processed with foil batch", input$optode_foil)
  })

  observeEvent(input$rinko, {
    # first calculate for all data (using untrimmed)
    rinko_temperature = rinko_temp(
        profiles$untrimmed[[input$select_profile]]@data[[input$rinko_T_channel]]
        )
    pressure =  unlist(profiles$untrimmed[[input$select_profile]]@data["pressure"])
    rinko_oxygen = rinko_o2(
        profiles$untrimmed[[input$select_profile]]@data[[input$rinko_O_channel]],
        rinko_temperature,
        S = profiles$untrimmed[[input$select_profile]]@data[["salinity"]],
        oC = rinko_coefs,
        G = input$rinko_G,
        H = input$rinko_H
        )
    # add to untrimmed
    profiles$untrimmed[[input$select_profile]] = ctdAddColumn(profiles$untrimmed[[input$select_profile]],
                                                         rinko_temperature, "temperature_RINKO", label = "temperature",
                                                         unit = list(name=expression(degree*C), scale="RINKO"))
    profiles$untrimmed[[input$select_profile]] = ctdAddColumn(profiles$untrimmed[[input$select_profile]],
                                                         rinko_oxygen, "oxygen_RINKO", label = "oxygen",
                                                         unit = list(name = expression(mmol~m-3), scale="RINKO"))
    # now subset and apply to $data, don't write subset to log
    x = subset(profiles$untrimmed[[input$select_profile]],
               scan >= min(profiles$data[[input$select_profile]][["scan"]]) &
               scan <= max(profiles$data[[input$select_profile]][["scan"]]))
    profiles$data[[input$select_profile]] = ctdAddColumn(profiles$data[[input$select_profile]],
                                                         x[["temperature_RINKO"]], "temperature_RINKO", label = "temperature",
                                                         unit = list(name=expression(degree*C), scale="RINKO"))
    profiles$data[[input$select_profile]] = ctdAddColumn(profiles$data[[input$select_profile]],
                                                         x[["oxygen_RINKO"]], "oxygen_RINKO", label = "oxygen",
                                                         unit = list(name = expression(mmol~m-3), scale="RINKO"))
    processingLog(profiles$data[[input$select_profile]]) = paste("RINKO processed with coefs", rinko_coefs[["serial"]],
                                                                 ",G =", input$rinko_G,
                                                                 ",H =", input$rinko_H)
  })

  observeEvent(input$licor, {
    licor_par = par_from_voltage(
        profiles$data[[input$select_profile]]@data[[input$par_channel]],
        input$licor_factor, input$licor_offset)

    # add to untrimmed
    profiles$untrimmed[[input$select_profile]] = oceSetData(profiles$untrimmed[[input$select_profile]],
                                                         "par", licor_par,
                                                         units = list(unit=expression(umol~s-1~m-2), scale = "PAR/Irradiance, Cefas Licor PAR"))
    profiles$data[[input$select_profile]] = oceSetData(profiles$data[[input$select_profile]],
                                                         "par", licor_par,
                                                         units = list(unit=expression(umol~s-1~m-2), scale = "PAR/Irradiance, Cefas Licor PAR"))
    processingLog(profiles$data[[input$select_profile]]) = paste("PAR processed with factor =", input$licor_factor, ",offset =", input$licor_offset)
  })

  observeEvent(input$flag_flu, {
    if("par" %in% names(profiles$data[[input$select_profile]]@data)){
      profiles$data[[input$select_profile]][["fluorescence"]] [profiles$data[[input$select_profile]][["par"]] > input$par_flu_threshold] = NA
      profiles$untrimmed[[input$select_profile]][["fluorescence"]] [profiles$untrimmed[[input$select_profile]][["par"]] > input$par_flu_threshold] = NA
      log = paste("fluorometry values where PAR >", input$par_flu_threshold, " have been flagged")
      processingLog(profiles$data[[input$select_profile]]) = log
    }
    else{
      warning("PAR has not been calculated")
    }
  })

  observeEvent(input$secondCT, {
    profiles$data[[input$select_profile]][["temperature"]] = profiles$data[[input$select_profile]][["temperature2"]]
    profiles$data[[input$select_profile]][["conductivity"]] = profiles$data[[input$select_profile]][["conductivity2"]]
    profiles$data[[input$select_profile]][["salinity"]] = profiles$data[[input$select_profile]][["salinity2"]]
    profiles$data[[input$select_profile]][["salinityDifference"]] = 0
    processingLog(profiles$data[[input$select_profile]]) = paste("Primary CT data replaced with that from secondary")
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
      for(p in names(profiles$data)){
        incProgress(1/length(profiles$data), detail = paste("writing", p))
        write.ctd(profiles$data[[p]], file = paste0(p, ".csv"))
      }
    })
  })

  observeEvent(input$make_netcdf, {
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
      need(any(grepl("filter_low_pass_A_vars = prDM",
                     profiles$data[[input$select_profile]]@metadata$header)),
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
    data.frame(profiles$data[[input$select_profile]]@data)
  })
  output$bottles = renderRHandsontable({
      # editable table
    validate(need(profiles$bottles, "bottle file not loaded"))
    rhandsontable(profiles$bottles, readOnly = T, digits = 6, highlightRow = T) %>%
      hot_col(c("bottle_sal", "bottle_O2", "bottle_Chl"), readOnly = F) %>%
      hot_col(c("salinity", "salinity2", "bottle_sal"), format = "0.0000") %>%
      hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE)
  })
  output$bottle_plot = renderPlot({
    if(input$Plot_bottle_select == "Salinity"){
      dat = hot_to_r(input$bottles)[bottle_sal != 0]
      plot(dat$bottle_sal, dat$salinity2,
           col = "green", xlab = "Niskin", ylab = "CT", main = "CT vs Niskin, Blue = Primary CT")
      points(dat$bottle_sal, dat$salinity, col = "blue")
      m = lm(data = dat, salinity ~ bottle_sal)
      abline(m)
    }
  })
  output$bottle_plot = renderPlot({
    if(input$Plot_bottle_select == "Salinity"){
      dat = hot_to_r(input$bottles)[bottle_sal != 0]
      m = lm(data = dat, salinity ~ bottle_sal)
      par(mfrow = c(1, 2))
      plot(dat$bottle_sal, dat$salinity,
           col = "blue", xlab = "Salinity bottle", ylab = "CT", main = "CT vs Niskin, Primary CT",
           sub = paste0("y = ", round(coef(m)[1], 3), " + ", round(coef(m)[2], 3),"x"))
      abline(m)
      hist(m$residuals, main = "Residuals (Primary CT)")
      profiles$bottle_coef[["salinity"]] = list(var = "salinity", slope = coef(m)[2], intercept = coef(m)[1])
    }
    if(input$Plot_bottle_select == "Oxygen Optode"){
      dat = hot_to_r(input$bottles)[bottle_O2 != 0]
      m = lm(data = dat, oxygen_optode ~ bottle_O2)
      par(mfrow = c(1, 2))
      plot(dat$bottle_O2, dat$oxygen_optode,
           col = "green", xlab = "Winkler", ylab = "Optode", main = "Optode vs Winkler",
           sub = paste0("y = ", round(coef(m)[1], 3), " + ", round(coef(m)[2], 3),"x"))
      abline(m)
      hist(m$residuals, main = "Residuals")
      profiles$bottle_coef[["oxygen_optode"]] = list(var = "oxygen_optode", slope = coef(m)[2], intercept = coef(m)[1])
    }
    if(input$Plot_bottle_select == "Oxygen RINKO"){
      dat = hot_to_r(input$bottles)[bottle_O2 != 0]
      m = lm(data = dat, oxygen_RINKO ~ bottle_O2)
      par(mfrow = c(1, 2))
      plot(dat$bottle_O2, dat$oxygen_RINKO,
           col = "green", xlab = "Winkler", ylab = "RINKO", main = "RINKO vs Winkler",
           sub = paste0("y = ", round(coef(m)[1], 3), " + ", round(coef(m)[2], 3),"x"))
      abline(m)
      hist(m$residuals, main = "Residuals")
      profiles$bottle_coef[["oxygen_RINKO"]] = list(var = "oxygen_RINKO", slope = coef(m)[2], intercept = coef(m)[1])
    }
  })
  output$bottle_coef = renderTable({
    print(rbindlist(profiles$bottle_coef))
    rbindlist(profiles$bottle_coef)
  })
  output$progress = renderTable({
    # fetch processing log then grep for string to check progress
    log = lapply(profiles$data , function(x) `@`( x , processingLog))
    done = grepl("complete", log, ignore.case=T)
    data.frame(
      "dip" = names(profiles$data),
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


