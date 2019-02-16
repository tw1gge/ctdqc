library(shiny,quietly=T)
library(shinyFiles,quietly=T)
library(oce, quietly=T)
library(data.table, quietly=T)
library(leaflet, quietly=T)
library(rhandsontable, quietly=T)
library(xml2, quietly=T)
library(stringr, quietly=T)
library(zoo, quietly=T)
library(lubridate, quietly=T)
library(uuid, quietly=T)
library(ggplot2, quietly=T)
library(DT, quietly = T)
library(signal, quietly = T)
library(geosphere, quietly=T)

source("functions.R", local = T)
CTDQC_version = "2.0"
editable_metadata = c("id", "title", "summary", "processing_level", "comment", "acknowledgment", "licence", "project", "creator", "creator_email")
sensor_metadata = fread("sensor_table.csv")
ctd_columns = list(
  PAR = list(name="par/sat/log", unit=list(expression(), scale="umol photon m-2"))
  )
vchannels = c("v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7")

# runUrl("https://bitbucket.org/betascoo8/ctdqc/get/dev.zip")

shinyServer(function(input, output, session) {

  # find OS disk drives
  volumes = getVolumesFast()
  # volumes = c("A:" = "C:/Users/th05/Dropbox (CEFAS)/CTD", volumes)
  volumes = c("A:" = "H:/Dropbox (CEFAS)/CTD", volumes)
  # if(dir.exists("\\\\lowfilecds\\Function\\SmartBuoyData\\CTD - SBE")){
    # volumes = c("SmartBuoyData:" = "\\\\lowfilecds\\Function\\SmartBuoyData\\CTD - SBE", "C:" = "C:")
  # }
  shinyDirChoose(input, 'directory', roots=volumes, session=session, restrictions=system.file(package='base'), updateFreq=500)
  output$directory = renderText({paste0(parseDirPath(volumes, input$directory), "/")})
  # make dynamic file list for storing the CTD objects, a list of S4 objects
  profiles = reactiveValues(data = NULL, bottles = NA)

  observeEvent(input$read_files, {
   ## read data
    if(is.null(input$directory)){ # stop crashing when you missclick
      showNotification("no folder selected!", type="error")
      return(NULL)
      }
    filelist = list.files(parseDirPath(volumes, input$directory), recursive=T, full.names = F, pattern = "*.cnv")
    if(length(filelist) == 0){ # don't crash if folder is empty
      showNotification("no .cnv files found in directory", type="error")
      return(NULL)
      }
    dir = parseDirPath(volumes, input$directory)
    d = list() # data
    m = list() # metadata
    h = list() # header (the xml)
    config = list()
    withProgress(message = 'loading files...', value = 0, {
      for(i in filelist){
        # Increment the progress bar, and update the detail text.
        incProgress(1/length(filelist), detail = paste("loading", i))
        print(i)
        d[[i]] = read.ctd.sbe(paste0(dir,"/",i), columns = ctd_columns) # oce data
        # d[[i]] = calc_descent_rate(d[[i]])
        if(!"pumpStatus" %in% names(d[[i]]@data)){
          print("pumpStatus not found, assuming pump is on")
          d[[i]]@data[["pumpStatus"]] = rep(1, length(d[[i]]@data$scan))
        }
        d[[i]] = flag_extra_pump(d[[i]], 100)
        h[[i]] = parse_sbe_xml(d[[i]])
        config[[i]] = extract.xml_channel_config(h[[i]])
      }
    })
      # check if processing steps have been completed
    headers = paste(extract.oce.metadata(d, "header"), collapse = ",")
    if(stringr::str_count(headers, "filter_low_pass_[AB]_vars = pr[Dd]M") < length(d)){
      showNotification("pressure filter has not been applied for all profiles!", duration=NULL, type="warning")
      }
    if(stringr::str_count(headers, "celltm_date") < length(d)){
      showNotification("Cell thermal mass correction has not been applied for all profiles!", duration=NULL, type="warning")
      }
    if(stringr::str_count(headers, "Derive_date") < length(d)){
      showNotification("Derived variables hve not been calculated for all profiles!", duration=NULL, type="warning")
    }

      # insert data into data slot
    profiles$data = d
    profiles$header = h
    profiles$config = config
    profiles$untrimmed = d
    profiles$original = d # make a backup for use by revert
      # make a summary of the positions for the map
    profiles$positions = extract.oce.metadata(profiles$data, c("filename", "startTime", "station", "longitude", "latitude", "cruise"))
    profiles$global_metadata = netcdf.metadata(profiles$data, profiles$positions)
    profiles$global_metadata_default = profiles$global_metadata
    profiles$global_config = unique(rbindlist(profiles$config))[order(channel)]

    if(length(unique(profiles$positions$cruise)) > 1){showNotification("Cruise ID differ between cnv files!", duration=NULL, type="warning")}
    if(length(unique(config)) != 1){
      showNotification("channel configuration differs between files", type="warning", duration=NULL)
      updateCheckboxInput(session, inputId = "apply_global", value = F)
      }
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
          # check for no-bottles-fired
          botfile = paste0(dir,"/",botfile)
          if(length(readLines(botfile)) > 2){
            scans = cbind(i, read.csv(botfile, skip = 2, header = F))
            colnames(scans) = c("profile", "fire_seq", "niskin", "dateTime", "start_scan", "end_scan")
            b[[i]] = scans
          }
        }
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
      # check CTDQC.rdata is valid
    if(exists("CTDQC_version", where=session)){
      profiles$data = session$data
      profiles$header = session$header
      profiles$config = session$config
      profiles$untrimmed = session$untrimmed
      profiles$original = session$original
      profiles$positions = session$positions
      profiles$global_config = session$global_config
      profiles$global_metadata = session$global_metadata
      profiles$global_metadata_default = session$global_metadata}
      if(!is.null(session$bottles)){
        profiles$bottles = session$bottles
      }
    else{
        showNotification("CTDQC error", type="error")
        warning("CTDQC error")
    }
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
    rate = 1 / profiles$data[[input$select_profile]]@metadata$sampleInterval # 2second window
    print(rate)
    decent_calc = zoo::rollmean(decent_diff, rate * input$inversion_window, fill=NA, align="right") * rate
    profiles$data[[input$select_profile]]@data = lapply(profiles$data[[input$select_profile]]@data, function(x) {
      # remove data for all parameters!
      x[decent_calc < input$decent_threshold] = NA
      return(x)
    })
    processingLog(profiles$data[[input$select_profile]]) = paste("Pressure inversions removed, minimum speed = ", input$decent_threshold,"m/s")
  })

  observeEvent(input$decimate,{
    # TODO move all decimate function to csv writer
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
    profiles$untrimmed = lapply(profiles$untrimmed, function(x) {
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
    if(input$apply_global){
      ilst = names(profiles$data)
    }else{
      ilst = input$select_profile
    }
    for(i in ilst){
      # first calculate for all data (using untrimmed)
      optode_temperature = optode.analogtemp(profiles$untrimmed[[i]]@data[[input$optode_T_channel]])
      optode_Dphase = optode.analogDphase(profiles$untrimmed[[i]]@data[[input$optode_dphase_channel]])
      salinity = profiles$untrimmed[[i]]@data[["salinity"]]
      depth = profiles$untrimmed[[i]]@data[["depth"]]
      optode_oxygen = optode.phaseCalc(optode_Dphase, optode_temperature, subset(optode_coefs, batch == input$optode_foil))
      optode_oxygen = optode.correction(optode_oxygen, optode_temperature, salinity, depth)
      # add to untrimmed
      profiles$untrimmed[[i]] = oceSetData(profiles$untrimmed[[i]], "temperature_optode", optode_temperature,
                                             unit = list(unit=expression(degree*C), scale="Optode"))
      profiles$untrimmed[[i]] = oceSetData(profiles$untrimmed[[i]], "oxygen_optode", optode_oxygen,
                                             unit = list(unit=expression(mmol~m-3), scale="Optode"))
      # now subset and apply to $data, don't write subset to log
      scans = profiles$data[[i]][["scan"]]
      x = subset(profiles$untrimmed[[i]],
                 scan >= scans[1] &
                 scan <= scans[length(scans)])
      profiles$data[[i]] = oceSetData(profiles$data[[i]], "temperature_optode", x[["temperature_optode"]],
                                        unit = list(unit=expression(degree*C), scale="Optode"))
      profiles$data[[i]] = oceSetData(profiles$data[[i]], "oxygen_optode", x[["oxygen_optode"]],
                                        unit = list(unit=expression(mmol~m-3), scale="Optode"))
      processingLog(profiles$data[[i]]) = paste("Optode processed with foil batch", input$optode_foil)
    }
  })

  observeEvent(input$rinko, {
    if(input$apply_global){
      ilst = names(profiles$data)
    }else{
      ilst = input$select_profile
    }
    for(i in ilst){
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
      scans = profiles$data[[i]][["scan"]]
      x = subset(profiles$untrimmed[[i]],
                 scan >= scans[1] &
                 scan <= scans[length(scans)])

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
    if(input$apply_global){
      ilst = names(profiles$data)
    }else{
      ilst = input$select_profile
    }
    for(i in ilst){
      # calculate using untrimmed data
      licor_par = par_from_voltage(profiles$untrimmed[[i]]@data[[input$par_channel]], input$licor_factor, input$licor_offset)
      profiles$untrimmed[[i]] = oceSetData(profiles$untrimmed[[i]], "par", licor_par,
                                           unit = list(unit=expression(umol~s-1~m-2), scale = "PAR/Irradiance, Cefas Licor PAR"))
      # now subset and add to data
      scans = profiles$data[[i]][["scan"]]
      x = subset(profiles$untrimmed[[i]],
                 scan >= scans[1] &
                 scan <= scans[length(scans)])

      log = paste("PAR processed with factor =", input$licor_factor, ",offset =", input$licor_offset)
      profiles$data[[i]] = oceSetData(profiles$data[[i]], "par", x[["par"]],
                                      unit = list(unit=expression(umol~s-1~m-2), scale = "PAR/Irradiance, Cefas Licor PAR"),
                                      note = log)
    }
  })

  observeEvent(input$flag_par, {
    # TODO make work with global? and based on lat/lon/time
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
    # TODO make work with global?
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
    # TODO make work with global?
    try({
      profiles$data[[input$select_profile]][["temperature"]] = profiles$data[[input$select_profile]][["temperature2"]]
      profiles$data[[input$select_profile]][["conductivity"]] = profiles$data[[input$select_profile]][["conductivity2"]]
      profiles$data[[input$select_profile]][["salinity"]] = profiles$data[[input$select_profile]][["salinity2"]]
      profiles$data[[input$select_profile]][["salinityDifference"]] = 0
      processingLog(profiles$data[[input$select_profile]]) = paste("Primary CT data replaced with that from secondary")
      profiles$untrimmed[[input$select_profile]][["temperature"]] = profiles$data[[input$select_profile]][["temperature2"]]
      profiles$untrimmed[[input$select_profile]][["conductivity"]] = profiles$data[[input$select_profile]][["conductivity2"]]
      profiles$untrimmed[[input$select_profile]][["salinity"]] = profiles$data[[input$select_profile]][["salinity2"]]
      profiles$untrimmed[[input$select_profile]][["salinityDifference"]] = 0
      showNotification("Primary CT data replaced with that from secondary")
    })
  })

  observeEvent(input$remove_variable,{
    #
    try({
      profiles$data[[input$select_profile]][[input$remove_variable_var]] = NULL
      profiles$untrimmed[[input$select_profile]][[input$remove_variable_var]] = NULL
      processingLog(profiles$data[[input$select_profile]]) = paste("removed", input$remove_variable_var, "from data")
    })
  })

  ## Write out
  observeEvent(input$write_rdata,{
    dir = parseDirPath(volumes, input$directory)
    session = list()
    # 7 + 1 items
    session$data = profiles$data
    session$header = profiles$header
    session$config = profiles$config
    session$untrimmed = profiles$untrimmed
    session$original = profiles$original
    session$positions = profiles$positions
    session$global_config = profiles$global_config
    session$global_metadata = profiles$global_metadata
    session$global_metadata_default = profiles$global_metadata
    session$CTDQC_version = CTDQC_version
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
        p = substr(p, 1, nchar(p)-4) # drop the .cnv from filename
        write.csv(d, file = paste0(dir, "/", p, ".csv"))
        write.csv(as.data.table(log[[p]]), file = paste0(dir, "/", p, ".log"))
      }
    })
  })

  observeEvent(input$write_netcdf, {
    publish_param = data.table(hot_to_r(input$publish_param))[publish == T]
    fn = write.ctd.netcdf(profiles, sensor_metadata, publish_param, input$decimate, parseDirPath(volumes, input$directory))
    showNotification(paste("NetCDF written", fn))
  })

  observeEvent(input$prev_filter, {
    time_constant = input$filter_t
    sample_rate = 1 / profiles$data[[input$select_profile]]@metadata$sampleInterval # 2second window
    Wn = (1 / time_constant) / (sample_rate * 2)
    flt = signal::butter(2, Wn, "low")
    var = profiles$untrimmed[[input$select_profile]]@data[[input$filter_x1]]
    var = zoo::na.approx(var, na.rm=F, rule=2)
    profiles$prev_filter = signal::filtfilt(flt, var)
  })

  observeEvent(input$apply_filter, {
    time_constant = input$filter_t
    sample_rate = 1 / profiles$data[[input$select_profile]]@metadata$sampleInterval
    Wn = (1 / time_constant) / (sample_rate * 2)
    flt = signal::butter(2, Wn, "low")
    untrimmed = profiles$untrimmed[[input$select_profile]]@data[[input$filter_x1]]
    untrimmed = zoo::na.approx(untrimmed, na.rm=F, rule=2)
    untrimmed = signal::filtfilt(flt, untrimmed)
    profiles$untrimmed[[input$select_profile]]@data[[input$filter_x1]] = untrimmed
    # pull out the scans from untrimmed that we need, don't rerun filter.
    scan_range = range(profiles$data[[input$select_profile]]@data$scan, na.rm=T)
    data_ = profiles$untrimmed[[input$select_profile]]
    data_ = subset(data_, scan %between% scan_range)
    profiles$data[[input$select_profile]]@data[[input$filter_x1]] = data_@data[[input$filter_x1]]
    processingLog(profiles$data[[input$select_profile]]) = paste(input$filter_x1, "low pass filterd with ", input$filter_t, "second time constant")
  })


  ## Ui and controls
    # update select input when filelist changes
  observe({
    updateSelectInput(session, "select_profile", choices = names(profiles$original))
    })
  observe({
      validate(need(!is.null(profiles$data[[input$select_profile]]), "Data not loaded"))
      choices = names(profiles$data[[input$select_profile]]@data)

      input_x1_prev = isolate(input$x1)
      input_x2_prev = isolate(input$x2)
      input_y_prev = isolate(input$y)

      if("salinity" %in% choices){default_x1 = "salinity"}else{default_x1 = NULL}
      if("temperature" %in% choices){default_x2 = "temperature"}else{default_x2 = NULL}
      if("pressure" %in% choices){default_y = "pressure"}else{default_y = NULL}

      if(input_x1_prev %in% choices){default_x1 = input_x1_prev}
      if(input_x2_prev %in% choices){default_x2 = input_x2_prev}
      if(input_y_prev %in% choices){default_y = input_y_prev}

      updateSelectInput(session, "x1", choices = choices, selected = default_x1)
      updateSelectInput(session, "filter_x1", choices = choices, selected = default_x1)
      updateSelectInput(session, "x2", choices = choices, selected = default_x2)
      updateSelectInput(session, "y", choices = choices, selected = default_y)
      updateSelectInput(session, "select_factor", choices = choices)
      updateSelectInput(session, "select_flag", choices = choices)
      updateSelectInput(session, "remove_variable_var", choices = choices)
    })

  ## Output

  output$summary <- renderPrint({
    # workaround for oce function not liking null data
    validate(need(!is.null(profiles$data[[input$select_profile]]), "Data not loaded"))
    summary(profiles$data[[input$select_profile]])
    })
  output$xml <- renderText({
    validate(need(!is.null(profiles$data[[input$select_profile]]), "Data not loaded"))
    paste(profiles$data[[input$select_profile]]@metadata$header, collapse="\n")
    })
  output$config <- renderTable({
    validate(need(!is.null(profiles$data[[input$select_profile]]), "Data not loaded"))
    profiles$config[[input$select_profile]]

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
      validate(need(!is.null(profiles$data[[input$select_profile]]), "Data not loaded"))
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
    })
  filter_plot_ranges = reactiveValues(x=NULL, y=NULL)
  observeEvent(input$filter_plot_dblclick, {
    if(!is.null(input$filter_plot_brush)){
      filter_plot_ranges$x = c(input$filter_plot_brush$xmin, input$filter_plot_brush$xmax)
      filter_plot_ranges$y = c(input$filter_plot_brush$ymin, input$filter_plot_brush$ymax)

    }else{
      filter_plot_ranges$x = NULL
      filter_plot_ranges$y = NULL
    }
  })
  output$filter_plot = renderPlot({
    validate(need(!is.null(profiles$untrimmed[[input$select_profile]]), "Data not loaded"))
    pd = as.data.table(profiles$untrimmed[[input$select_profile]]@data)
    if(!is.null(profiles$prev_filter)){
      pd[, filter := profiles$prev_filter]
    }else{
      pd[, filter := NaN] # so ggplot won't plot it
    }
    if(exists("pumpStatus", where=pd)){
      pd = pd[pumpStatus == 1 & pressure > 3]
    }else{
      pd = pd[pressure > 3]
    }
    ggplot(pd) +
      geom_path(aes(time, filter), color="red") +
      geom_path(aes(time, get(input$filter_x1)), alpha=0.5) +
      theme_bw() +
      labs(x="time (seconds)", y=input$filter_x1) +
      coord_cartesian(xlim = filter_plot_ranges$x, ylim = filter_plot_ranges$y, expand=F)
  })
  output$TS_plot = renderPlot({
    validate(need(!is.null(profiles$data[[input$select_profile]]), "Data not loaded"))
    plotTS(profiles$data[[input$select_profile]])
  })
  output$hyst_plot = renderPlot({
    validate(need(!is.null(profiles$untrimmed[[input$select_profile]]), "Data not loaded"))
    sample_rate = 1 / profiles$data[[input$select_profile]]@metadata$sampleInterval
    if(input$lag < 0){
      lag_type = "lead"
      lag_mod = -1
    }else{
      lag_type = "lag"
      lag_mod = 1
    }
    if(input$CT_mode){
      pd = as.data.table(profiles$data[[input$select_profile]]@data)
      pd[, lagged_C := shift(pd$conductivity, type=lag_type, input$lag * lag_mod)]
      pd[, recalc_S := oce::swSCTp(lagged_C, temperature, pressure, conductivityUnit="S/m")]
      ggplot(pd) +
        geom_path(aes(salinity, pressure), color="black", alpha=0.2) +
        geom_path(aes(recalc_S, pressure), color="red",  alpha=0.5) +
        theme_bw() + scale_y_reverse() + scale_color_discrete("") +
        labs(x = "salinity", y="pressure", title=paste("sample rate", sample_rate))
    }else{
      pd = as.data.table(profiles$untrimmed[[input$select_profile]]@data)
      scan_max_pressure = min(pd[pressure == max(pressure)]$scan)
      pd[, dir := "down"]
      pd[scan > scan_max_pressure, dir := "up"]
      if(exists("pumpStatus", where=pd)){
        pd = pd[pumpStatus == 1 & pressure > 3]
      }else{
        pd = pd[pressure > 3]
      }
      lagged = shift(pd[[input$filter_x1]], type=lag_type, input$lag * lag_mod)
      ggplot(pd) +
        geom_path(aes(get(input$filter_x1), pressure, color=dir), alpha=0.3) +
        geom_path(aes(lagged, pressure, color=dir)) +
        theme_bw() + scale_y_reverse() + scale_color_discrete("") +
        labs(x = input$filter_x1, y="pressure", title=paste("sample rate", sample_rate)) +
        theme(legend.position="bottom")
    }
  })
  observeEvent(input$apply_lag,{
    if(input$lag < 0){
      lag_type = "lead"
      lag_mod = -1
    }else{
      lag_type = "lag"
      lag_mod = 1
    }
    lagged = profiles$data[[input$select_profile]]@data[[input$filter_x1]]
    lagged = shift(lagged, type=lag_type, input$lag * lag_mod)
    profiles$data[[input$select_profile]]@data[[input$filter_x1]] = lagged
    processingLog(profiles$data[[input$select_profile]]) = paste(input$filter_x1, "lagged by", input$lag, "scans")

    lagged = profiles$untrimmed[[input$select_profile]]@data[[input$filter_x1]]
    lagged = shift(lagged, type=lag_type, input$lag * lag_mod)
    profiles$untrimmed[[input$select_profile]]@data[[input$filter_x1]] = lagged
  })
  output$map = renderLeaflet({
    validate(need(!is.null(profiles$data[[input$select_profile]]), "Data not loaded"))
    leaflet(profiles$positions) %>%
      addTiles() %>%
      addMarkers(~longitude, ~latitude, popup = ~paste("Station", station,"\r", "@", startTime)) %>%
      setView(lat = profiles$data[[input$select_profile]]@metadata$latitude,
              lng = profiles$data[[input$select_profile]]@metadata$longitude,
              zoom = 7)
  })
  output$drift_plot = renderPlot({
    validate(need(!is.null(profiles$untrimmed[[input$select_profile]]), "Data not loaded"))
    pd = as.data.table(profiles$untrimmed[[input$select_profile]]@data)
    if(exists("pumpStatus", where=pd)){
      pd = pd[pumpStatus == 1 & pressure > 1]
    }
    start_pos = pd[1,.(longitude, latitude)]
    end_pos = pd[nrow(pd),.(longitude, latitude)]
    dist_covered = round(geosphere::distGeo(start_pos, end_pos))
    ggplot(pd) +
      geom_path(aes(longitude, latitude, color=time)) +
      viridis::scale_color_viridis() + theme_bw() +
      labs(title=paste("distance covered =", dist_covered, "(m)"), x="", y="") +
      coord_map()
  })
  output$datatable = renderDataTable({
    validate(need(profiles$untrimmed, "data not loaded"))
    data.frame(profiles$untrimmed[[input$select_profile]]@data)
  })
  output$bottles = renderRHandsontable({
    validate(need(profiles$bottles, "bottle file not loaded"))
        # identify the salinity columns so we can add more decimal places later
    avail_sal_cols = colnames(profiles$bottles)[chmatch(c("salinity", "salinity2", "bottle_sal"), colnames(profiles$bottles))]
      # editable table
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
  }, digits = 4)
  output$chl_coef = renderTable({
    data.frame(profiles$bottle_coef[["fluorescence"]])
  })
  output$progress = DT::renderDataTable({
    # fetch processing log then grep for string to check progress
    validate(need(profiles$data, ""))
    log = lapply(profiles$data , function(x) `@`( x , processingLog))
    sensor = grepl("processed", log, ignore.case=T)
    trim = grepl("ctdTrim", log, ignore.case=T)
    QC2 = grepl("QC2", log, ignore.case=T)
    done = grepl("QC Complete", log, ignore.case=T)
    data.frame("dip" = names(profiles$data), "trim" = trim, "sensor" = sensor, "QC2" = QC2, "done" = done )
  }, server=T)

  output$edit_metadata = renderUI({
    validate(need(profiles$global_metadata_default, "data not loaded"))
    lapply(editable_metadata, function(i){
      # generate UI dynamically
      id = paste0("netcdf-", i)
      default_value = profiles$global_metadata_default[[i]] # use copy to avoid instant-writeback
      textInput(id, i, value=default_value)
    })
  })

  output$publish_param = renderRHandsontable({
    validate(need(profiles$global_config, "data not loaded"))
    tbl = generate_parameter_table(profiles, sensor_metadata)
    tbl$publish = T
    rhandsontable(tbl, readOnly=F) %>%
      hot_col(c("parameter", "instid", "sbe_name"), readOnly=T)
  })

    # generate dynamic UI depending on which sensors are in the config
  output$sensor_ui = renderUI({
    validate(need(profiles$config, "data not loaded"))
    ui = list()

    if(any(grepl("Conductivity, 2", profiles$global_config$name))){
      # add ui to ui list
      ui <- list(ui, wellPanel(
        h4("Secondary CT"),
        actionButton('secondCT', "Overwrite Primary CT with secondary", icon=icon("reply-all"))
      ))
    }

    if(any(grepl("optode", profiles$global_config$comment, ignore.case=T))){
      ui = list(ui, wellPanel(
        h4("Optode"),
        selectInput("optode_foil", "Optode foil Batch #", choices = unique(optode_coefs$batch), selected="1707", width="200px"),
        selectInput('optode_T_channel', "Optode Temperature channel", choices = vchannels, selected = "v7", width="200px"),
        selectInput('optode_dphase_channel', "Optode dPhase channel", choices = vchannels, selected = "v6", width="200px"),
        actionButton('optode', "Process Optode", icon=icon("life-ring"))
      ))
    }

    if(any(grepl("rinko", profiles$global_config$comment, ignore.case=T))){
      ui = list(ui, wellPanel(
        h4("RINKO"),
        selectInput('rinko_T_channel', "RINKO Temperature channel", choices = vchannels, selected = "v5", width="200px"),
        selectInput('rinko_O_channel', "RINKO Oxygen channel", choices = vchannels, selected = "v4", width="200px"),
        numericInput('rinko_G', label = "G Coefficent", value = 0, width="200px"),
        numericInput('rinko_H', label = "H Coefficent", value = 1, width="200px"),
        actionButton('rinko', "Process RINKO", icon=icon("times-circle-o"))
      ))
    }

    if(any(grepl("licor", c(profiles$global_config$name, profiles$global_config$comment), ignore.case=T))){
      ui = list(ui, wellPanel(
        h4("LiCor PAR"),
        selectInput('par_channel', "PAR channel", choices = vchannels, selected = "v0", width="200px"),
        numericInput('licor_factor', label = "Licor factor", value = 0.22345679, width="200px"),
        numericInput('licor_offset', label = "Licor offset", value = 3.3737, width="200px"),
        actionButton('licor', "Process Licor PAR", icon=icon("beer"))
      ))
    }
    if(any(grepl("par", profiles$global_config$name, ignore.case=T))){
      ui = list(ui, wellPanel(
        h4("PAR"),
        actionButton('flag_par', "Flag all PAR for selected dip (Night)", icon=icon("moon-o"))
      ))
    }

    if(any(grepl("fluoro", profiles$global_config$name, ignore.case=T))){
      ui = list(ui, wellPanel(
        h4("Fluorometer"),
        numericInput('par_flu_threshold', label = "Chlorophyll quenching PAR threshold", value = 1, width="200px"),
        actionButton('flag_flu', "Flag quenched chlorophyll fluorometry", icon=icon("ban")),
        tableOutput("chl_coef"),
        numericInput('chl_factor', label="Chl Factor", value=1.0, step=0.01, width="200px"),
        numericInput('chl_offset', label="Chl Offset", value=0.0, step=0.01, width="200px"),
        actionButton('calc_flu', "derive Chlorophyll from flu regression", icon=icon("leaf"))
      ))
    }

    ui
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


