
  #### Sensor coefficents

rinko_coefs = list(serial = "#0263 ARO-CAV", A = -42.34162, B = +127.6475, C = -0.3677435, D = +0.01137, E = +0.0046, F = +7.57e-05)
optode_coefs = read.csv("optode_coefs.csv")

  #### Sensor functions

optode.analogtemp <- function(v, TempLimit0=-5, TempLimit1=35){
  TempLimit0 + (v / 5) * diff(c(TempLimit0, TempLimit1))
}

optode.analogDphase <- function(v, PhaseLimit0=10, PhaseLimit1=70){
  PhaseLimit0 + (v / 5) * diff(c(PhaseLimit0, PhaseLimit1))
}

optode.phaseCalc <- function(DPhase, Temp, coefs){
    with(coefs, {
      if(coef[1] == "SVU"){
        # For 4831 multipoint calibrated optodes
        print(paste("using SVU foil batch coefs", batch[1]))
        Ksv = C0 + C1*Temp + C2*Temp^2
        P0 = C3 + C4*Temp
        Pc = C5 + C6*DPhase # actually calphase
        ((P0/Pc)-1) / Ksv
      }else{
        # for mkl optodes 3830 & 3835
        print(paste("using foil batch coefs", batch[1]))
        (C0[1]+C0[2]*Temp+C0[3]*Temp^2+C0[4]*Temp^3) +
        (C1[1]+C1[2]*Temp+C1[3]*Temp^2+C1[4]*Temp^3) *
        DPhase+(C2[1]+C2[2]*Temp+C2[3]*Temp^2+C2[4]*Temp^3) *
        DPhase^2+(C3[1]+C3[2]*Temp+C3[3]*Temp^2+C3[4]*Temp^3) *
        DPhase^3+(C4[1]+C4[2]*Temp+C4[3]*Temp^2+C4[4]*Temp^3) *
        DPhase^4
      }
    })
}

optode.correction <- function(O2, t, S, depth = 0, optode_salinity = 0){
  # corrects optode measurements for salinity and depth
    # oxygen units returned same as input
  # Solubility and salinity comp based on Garcia and Gordon, 1992. Limno Ocean
  pCoef = 0.032 # empricial derived pressure compensation coef from Uchida et al, 2008. J Atmos Ocean Tech
  B0 = -6.24097E-03
  B1 = -6.93498E-03
  B2 = -6.90358E-03
  B3 = -4.29155E-03
  C0 = -3.11680E-07

  Ts = log((298.15-t)/(273.15+t)) # scaled temperature

  O2c = O2 * exp(S * (B0 + B1 * Ts + B2 * Ts^2 + B3 * Ts^3) + C0 * S^2) /
      exp(optode_salinity* (B0 + B1 * Ts + B2 * Ts^2 + B3 * Ts^3) + C0^2 * optode_salinity^2) *
      (1 + (pCoef * abs(depth)) / 1000)
  # sal_factor = exp((S - optode_salinity_setting) * (B0 + B1*Ts + B2*Ts^2 + B3*Ts^3) + C0 * (S^2 - optode_salinity_setting^2))
  # prs_factor = (((abs(depth))/1000)*pCoef) + 1
  return(O2c)
}

par_from_voltage <- function(x, factor, offset){
  # for LiCor Li-192 PAR with Cefas low light amplifier
    return(factor * exp(offset * x))
}

rinko_temp <- function(V, tC = list(A = -5.326887e+00, B = +1.663288e+01, C = -2.123968e+00, D = +4.543014e-01)){
    # RINKO III defaults to #0263 ARO-CAV

  temp = tC$A + tC$B * V + tC$C * V^2 + tC$D * V^3
  return(temp)
}

rinko_o2 <- function(V, t, S, oC = list(A = -4.234162e+01,
                                     B = +1.276475e+02,
                                     C = -3.677435e-01,
                                     D = +1.137000e-02,
                                     E = +4.600000e-03,
                                     F = +7.570000e-05),
                     p = 0.1, G = 0, H = 1){

  # V = output voltage
  # t = tempeture from rinko_temp
  # p = in-situ pressure in decibar
  # G & H = RINKO calibration coefs (alpha and beta)

    # RINKO III #0263 ARO-CAV

  P1 = oC$A / (1 + oC$D * (t - 25) + oC$F * (t - 25)^2)
  P2 = oC$B / (V * (1 + oC$D * (t - 25) + oC$F * (t - 25)^2) + oC$C)
  P = P1 + P2

    # G and H are calibration coefs
  DO = G + H * P
    # pressure correction
  d = p * 0.01 # convert from decibar to MPa
  DO = DO * (1 + oC$E * d)

  # from garcia and gordon
    A0 = 2.00856
    A1 = 3.224
    A2 = 3.99063
    A3 = 4.80299
    A4 = 0.978188
    A5 = 1.71069
    B0 = -0.00624097
    B1 = -0.00693498
    B2 = -0.00690358
    B3 = -0.00429155
    C0 = -3.1168E-07
    Ts = log((298.15-t)/(273.15+t))

    Cstar = exp(A0+(A1*Ts)+(A2*Ts^2)+
                    (A3*Ts^3)+(A4*Ts^4)+(A5*Ts^5)+
                    S*(B0+(B1*Ts)+(B2*Ts^2)+(B3*Ts^3))+
                    (C0*S^2))

    DO = ((Cstar * 44.614 * DO) / 100) * 0.0319988 # mg/l
    DO = (DO/31.9988) * 1000 # mg/l to mmol m-3

  return(DO)
}

getVolumesFast <- function (exclude){
  # faster than default version in shinyFiles
  osSystem <- Sys.info()["sysname"]
  if (osSystem == "Darwin") {
    volumes <- list.files("/Volumes/", full.names = T)
    names(volumes) <- basename(volumes)
  }
  else if (osSystem == "Linux") {
    volumes <- c(Computer = "/")
    media <- list.files("/media/", full.names = T)
    names(media) <- basename(media)
    volumes <- c(volumes, media)
  }
  else if (osSystem == "Windows") {
    volumes <- system("wmic logicaldisk get Caption", intern = T)
    volumes <- sub(" *\\r$", "", volumes)
    keep <- !tolower(volumes) %in% c("caption", "")
    volumes <- volumes[keep]
    names(volumes) <- volumes
  }
  else {
    stop("unsupported OS")
  }
  volumes
}

extract.oce.metadata <- function(oce, vars){
  # function for fetching metadata variables from oce object
    rbindlist(lapply(oce , function(x) `@`( x , metadata)[vars]))
}

parse_sbe_xml <- function(oce){
  hdr = oce@metadata$header # extract information from sbe xml header
  startLine = grep(hdr, pattern="<Sensors count")
  endLine = grep(hdr, pattern="</Sensors>")
  hdr = hdr[startLine:endLine] # subset each element of list
  hdr = substring(hdr, 2) # drop the '#'
  hdr = paste(hdr, collapse="\n") # unvector it
  hdr = xml2::read_xml(hdr)
  return(hdr)
}

extract.xml_channel_config <- function(sbe_xml){
  # takes output of parse_sbe_xml and returns table of channels
  channel = xml_find_all(sbe_xml, "//sensor/@Channel")
  config = list()
  for(i in xml_integer(channel)){
    name = xml_text(xml_find_all(sbe_xml, paste0("//sensor[@Channel=",i,"]/comment()")))
    sen = xml_find_first(sbe_xml, paste0("//sensor[@Channel=",i,"]/*"))
    sbename = xml_name(sen)
    if(!is.na(sbename)){
      serial = xml_text(xml_find_first(sen, "SerialNumber"))
      comment = xml_text(xml_find_first(sen, "SensorName"))
      config[[i]] = data.table(name = name, sbename = sbename, serial = serial, comment = comment)
    }else{
      config[[i]] = data.table(name = "empty", sbename = "empty", serial = NA, comment = NA)
    }
  }
  config = rbindlist(config, idcol="channel")
  return(config)
}


calc_descent_rate <- function(oce, window_size=2){
  if(!exists("descentRate", where=oce@data)){
    prs = oce@data[["pressure"]]
    d_prs = c(0, diff(prs))
    rate = 1 / oce@metadata$sampleInterval
    decent_calc = zoo::rollmean(decent_diff, rate * input$inversion_window, fill=NA, align="right") * rate
    oce@data[["descentRate"]] = descentRate
  }
  return(oce)
}

netcdf.metadata <- function(d, positions){
  # generates metadata from file
  metadata = list()
  metadata[["id"]] = paste0("SBE_CTD_", paste(unique(positions$cruise), collapse="&"))
  metadata[["title"]] = paste("CTD profiles from", paste(unique(positions$cruise), collapse=", "))
  metadata[["ncei_template_version"]] = "NCEI_NetCDF_Profile_Orthogonal_Template_v2.0"
  metadata[["featureType"]] = "profile"
  metadata[["summary"]] = ""
  metadata[["keywords"]] = ""
  metadata[["keywords_vocabluary"]] = "GCMD:GCMD Keywords"
  metadata[["Conventions"]] = "CF-1.6, ACDD-1.3"
  metadata[["naming_authority"]] = "uk.co.cefas"
  metadata[["history"]] = ""
  metadata[["source"]] = "CTD rossette system"
  metadata[["processing_level"]] = paste0("QC and processing done following Cefas SBE CTD QC SOP v", CTDQC_version)
  metadata[["comment"]] = ""
  metadata[["acknowledgement"]] = ""
  metadata[["licence"]] = "OGLv3.0 https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/"
  metadata[["standard_name_volcabulary"]] = "CF Standard Name Table v39"
  metadata[["date_created"]] = strftime(lubridate::now(tzone="UTC"), "%Y-%M-%dT%H:%M:%SZ")
  metadata[["creator"]] = "Tom Hull"
  metadata[["creator_email"]] = "tom.hull@cefas.co.uk"
  metadata[["creator_url"]] = "http://www.cefas.co.uk"
  metadata[["institution"]] = "Centre for Environment Fisheries and Aquaculture Science"
  metadata[["project"]] = ""
    # note geospatial ranges are calculated later
  metadata[["geospatial_lat_units"]]= "degrees_north"
  metadata[["geospatial_lon_units"]] = "degrees_east"
  metadata[["time_coverage_start"]] = strftime(min(positions$startTime), "%Y-%M-%dT%H:%M:%SZ")
  metadata[["time_coverage_end"]] = strftime(max(positions$startTime), "%Y-%M-%dT%H:%M:%SZ")
  metadata[["cdm_data_type"]] = "Station"
  metadata[["uuid"]] = uuid::UUIDgenerate(use.time=T)
  metadata[["references"]] = "doi::"
  return(metadata)
}

flag_extra_pump <- function(oce, scans=50){
  if(exists("pumpStatus", where=oce@data)){
    pumpStatus = oce@data[["pumpStatus"]]
    scan = oce@data[["scan"]]
    if(any(pumpStatus == 1)){
      pumpStatus[1:min(scan[pumpStatus == 1]) + scans] = 0
    }
    oce@data[["pumpStatus"]] = pumpStatus
  }
  return(oce)
}

write.ctd.netcdf <- function(session, sensor_metadata){
  require(RNetCDF)
  require(uuid)
  require(reshape2)
  # load("C:/CEND_22_16/CTDQC.rdata")
  # sensor_metadata = fread("sensor_table.csv")

  # validates that all QC steps are done
  log = rbindlist(lapply(session$data , function(x) as.data.frame(`@`( x , processingLog))), idcol=T)
  logsummary = log[,.(QC = any(grepl("QC Complete", value))), by=.id] # are any of the values...

  if(all(logsummary$QC) == F){
    warning("WARNING - QC not complete")
    return(NULL)
  }

  # reshape and process data
  pressures = lapply(session$data , function(x) as.data.frame(`@`( x , data)["pressure"]))
  pressures = rbindlist(pressures, idcol="id")
  pressures[, max := max(pressure), by=id]
  deepest_dip = pressures[max(max) == max]

  v_lat = lapply(session$data , function(x) as.data.frame(`@`( x , metadata)["latitude"]))
  v_lat = rbindlist(v_lat, idcol="id")
  v_lon = lapply(session$data , function(x) as.data.frame(`@`( x , metadata)["longitude"]))
  v_lon = rbindlist(v_lon, idcol="id")
  v_time = lapply(session$data , function(x) as.data.frame(`@`( x , metadata)["startTime"]))
  v_time = rbindlist(v_time, idcol="id")
  v_station = lapply(session$data , function(x) as.data.frame(`@`( x , metadata)["station"]))
  v_station = rbindlist(v_station, idcol="id")

    # check station numbers unique

  sb = rbindlist(lapply(session$data , function(x) as.data.frame(`@`( x , data))),
                 idcol="id", use.names=T, fill=T)
  sb = merge(v_station, sb, by="id", all.x = T)

  sensor_metadata = sensor_metadata[parameter %in% colnames(sb)]

  # setup netcdf
  nc = create.nc(paste0(session$global_metadata$id, ".nc"))

  # Dimensions
  dim.def.nc(nc, "pressure", length(deepest_dip$pressure))
  dim.def.nc(nc, "profile", length(session$data))

  var.def.nc(nc, "pressure", "NC_FLOAT", "pressure")
    att.put.nc(nc, "pressure", "standard_name", "NC_CHAR", "sea_water_pressure")
    att.put.nc(nc, "pressure", "units", "NC_CHAR", "dbar")
    att.put.nc(nc, "pressure", "axis", "NC_CHAR", "pressure")
    att.put.nc(nc, "pressure", "positive", "NC_CHAR", "down")
    att.put.nc(nc, "pressure", "comment", "NC_CHAR", "pressure from Digiquartz pressure transducer")
    var.put.nc(nc, "pressure", deepest_dip$pressure)

  var.def.nc(nc, "profile", "NC_INT", "profile")
    att.put.nc(nc, "profile", "long_name", "NC_CHAR", "Unique identifier for each feature instance")
    att.put.nc(nc, "profile", "cf_role", "NC_CHAR", "profile_id")
    var.put.nc(nc, "profile", as.integer(v_station$station))

  var.def.nc(nc, "time", "NC_DOUBLE", "profile")
    att.put.nc(nc, "time", "standard_name", "NC_CHAR", "time")
    att.put.nc(nc, "time", "units", "NC_CHAR", "seconds since 1970-01-01 00:00:00 0:00")
    att.put.nc(nc, "time", "calendar", "NC_CHAR", "gregorian") # note POSIXct assumes gregorian calendar
    att.put.nc(nc, "time", "axis", "NC_CHAR", "T")
    att.put.nc(nc, "time", "_FillValue", "NC_DOUBLE", -99999)
    att.put.nc(nc, "time", "comment", "NC_CHAR", "Time from ship GPS")
    var.put.nc(nc, "time", as.numeric(v_time$startTime))

  var.def.nc(nc, "lat", "NC_DOUBLE", "profile")
    att.put.nc(nc, "lat", "standard_name", "NC_CHAR", "latitude")
    att.put.nc(nc, "lat", "units", "NC_CHAR", "degrees_north")
    att.put.nc(nc, "lat", "axis", "NC_CHAR", "Y")
    att.put.nc(nc, "lat", "valid_min", "NC_FLOAT", -90)
    att.put.nc(nc, "lat", "valid_max", "NC_FLOAT", 90)
    att.put.nc(nc, "lat", "grid_mapping", "NC_CHAR", "crs")
    att.put.nc(nc, "lat", "comment", "NC_CHAR", "Position from ship GPS")
    var.put.nc(nc, "lat", v_lat$latitude)

  var.def.nc(nc, "lon", "NC_DOUBLE", "profile")
    att.put.nc(nc, "lon", "standard_name", "NC_CHAR", "longitude")
    att.put.nc(nc, "lon", "units", "NC_CHAR", "degrees_east")
    att.put.nc(nc, "lon", "axis", "NC_CHAR", "X")
    att.put.nc(nc, "lon", "valid_min", "NC_FLOAT", -180)
    att.put.nc(nc, "lon", "valid_max", "NC_FLOAT", 180)
    att.put.nc(nc, "lon", "grid_mapping", "NC_CHAR", "crs")
    att.put.nc(nc, "lon", "comment", "NC_CHAR", "Position from ship GPS")
    var.put.nc(nc, "lon", v_lon$longitude)

  var.def.nc(nc, "crs", "NC_INT", NA) # WGS84
    att.put.nc(nc, "crs", "grid_mapping_name", "NC_CHAR", "latitude_longitude")
    att.put.nc(nc, "crs", "espg_code", "NC_CHAR", "ESPG:4326")
    att.put.nc(nc, "crs", "longitude_of_prime_meridian", "NC_DOUBLE", 0.0)
    att.put.nc(nc, "crs", "semi_major_axis", "NC_DOUBLE", 6378137.0)
    att.put.nc(nc, "crs", "inverse_flattening", "NC_DOUBLE", 298.257223563)

  for(var in sensor_metadata$variable){
      # data
    metadata = sensor_metadata[variable == var]
      # covert to array for netcdf
    v_ = sb[,c("pressure", "station", metadata$parameter), with=F]
    v_ = reshape2::acast(v_, pressure ~ station)

    var.def.nc(nc, var, "NC_DOUBLE", c("pressure", "profile"))
      att.put.nc(nc, var, "long_name", "NC_CHAR", metadata$long_name)
      att.put.nc(nc, var, "standard_name", "NC_CHAR", metadata$standard_name)
      att.put.nc(nc, var, "units", "NC_CHAR", metadata$units)
      att.put.nc(nc, var, "_FillValue", "NC_DOUBLE", as.numeric(metadata$`_fillValue`))
      att.put.nc(nc, var, "valid_min", "NC_DOUBLE", as.numeric(metadata$valid_min))
      att.put.nc(nc, var, "valid_max", "NC_DOUBLE", as.numeric(metadata$valid_max))
      att.put.nc(nc, var, "coordinates", "NC_CHAR", "time lat lon pressure" )
      att.put.nc(nc, var, "instrument", "NC_CHAR", metadata$instid)
      att.put.nc(nc, var, "grid_mapping", "NC_CHAR", "crs" )
      att.put.nc(nc, var, "coverage_content_type", "NC_CHAR", "physicalMeasurement" )
      var.put.nc(nc, var, v_)

    var.def.nc(nc, metadata$instid, "NC_BYTE", "profile")
      att.put.nc(nc, metadata$instid, "long_name", "NC_CHAR", metadata$sensor_name)
      att.put.nc(nc, metadata$instid, "nodc_name", "NC_CHAR", metadata$sensor_nodc)
      att.put.nc(nc, metadata$instid, "make_model", "NC_CHAR", metadata$sensor_make)
      att.put.nc(nc, metadata$instid, "serial_number", "NC_CHAR", as.character(metadata$serial))
      att.put.nc(nc, metadata$instid, "precision", "NC_CHAR", as.character(metadata$precision))
  }

  for(g in names(session$global_metadata)){
    att.put.nc(nc, "NC_GLOBAL", g, "NC_CHAR", session$global_metadata[[g]])
  }
  att.put.nc(nc, "NC_GLOBAL", "geospatial_lat_min", "NC_FLOAT", min(v_lat$latitude))
  att.put.nc(nc, "NC_GLOBAL", "geospatial_lat_max", "NC_FLOAT", max(v_lat$latitude))
  att.put.nc(nc, "NC_GLOBAL", "geospatial_lon_min", "NC_FLOAT", min(v_lon$longitude))
  att.put.nc(nc, "NC_GLOBAL", "geospatial_lon_max", "NC_FLOAT", max(v_lon$longitude))

  sync.nc(nc)
  close.nc(nc)
}

read.sbe.ros <- function(file=NA, folder){
  if(is.na(file)){
    file_list = list.files(folder, pattern="*.ros")
    output = list()
    for(f in file_list){
      print(f)
      oce = oce::read.ctd.sbe(paste0(folder,"/", f))
      output[[f]] = as.data.table(oce@data)
      output[[f]][, dateTime := oce@metadata$startTime + time]
    }
    output = rbindlist(output, idcol="filename")
  }else{
      oce = oce::read.ctd.sbe(file)
      output = as.data.table(oce@data)
      output[, filename := file]
      output[, dateTime := oce@metadata$startTime + time]
  }
  output = output[,lapply(.SD, median), by=list(filename, bottlesFired)]
  output[, dateTime := as.POSIXct(dateTime, origin="1970-01-01", tz="UTC")]

  return(output)
}
