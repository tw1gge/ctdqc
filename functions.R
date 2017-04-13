
  #### Sensor coefficents

rinko_coefs = list(serial = "#0263 ARO-CAV", A = -42.34162, B = +127.6475, C = -0.3677435, D = +0.01137, E = +0.0046, F = +7.57e-05)
optode_coefs = read.csv("optode_coefs.csv")

  #### Sensor functions

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

extract.metadata <- function(oce, vars){
  # function for fetching metadata variables from oce object
    rbindlist(lapply(oce , function(x) `@`( x , metadata)[vars]))
}

parse_sbe_xml <- function(oce){
  # function to extract information from sbe xml header
  x = lapply(oce , function(x) `@`( x , metadata)[["header"]])
  startLines = lapply(x, grep, pattern="<Sensors count")
  endLines = lapply(x, grep, pattern="</Sensors>")
  x2 = lapply(1:length(x), function(i) x[[i]][startLines[[i]]:endLines[[i]]]) # subset each element of list
  if(length(unique(x2)) != 1){warning("xml header differs between files")} # check headers are all the same
  x = substring(x2[[1]], 2) # take one copy and remove first "#"
  x = xml2::read_xml(paste(x, collapse="\n")) # convert to xml
  sensors = melt(as_list(xml2::xml_find_all(x, "//sensor"))) # extract sensor nodes and convert to list
  sensors = subset(sensors, L3 %in% c("SensorName", "SerialNumber", "CalibrationDate", "GainSetting")) # subset just params we want
  sensors = dcast(sensors, L1 + L2 ~ L3) # cast to wide table
  return(sensors)
}

netcdf.metadata <- function(d, positions){
  # generates metadata from file
  metadata = list()
  metadata[["id"]] = paste0("SBE_CTD_", paste(unique(positions$cruise), collapse="&"))
  metadata[["title"]] = paste("CTD profiles from", paste(unique(positions$cruise), collapse=", "))
  metadata[["ncei_template_version"]] = "NCEI_NetCDF_TimeSeriesProfile_Orthogonal_Template_v2.0"
  metadata[["featureType"]] = "timeSeriesProfile"
  metadata[["summary"]] = ""
  metadata[["keywords"]] = ""
  metadata[["keywords_vocabluary"]] = "GCMD:GCMD Keywords"
  metadata[["Conventions"]] = "CF-1.6, ACDD-1.3"
  metadata[["naming_authority"]] = "uk.co.cefas"
  metadata[["history"]] = ""
  metadata[["source"]] = "CTD rossette system"
  metadata[["processing_level"]] = "QC and processing done following Cefas SBE CTD QC SOP v1.1"
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

check.qc.done <- function(session){
  # validates that all QC steps are done
  log = rbindlist(lapply(session$data , function(x) `@`( x , processingLog)), idcol=T)
  number_dips = length(session$data)
  log = as.data.frame(table(log$value))
  grep("QC Complete", log)
  return(T)
}

write.ctd.netcdf <- function(session){
  require(RNetCDF)
  require(uuid)
  require(reshape2)

  # reshape and process data
  pressures = rbindlist(lapply(session$data , function(x) `@`( x , data)["pressure"]), idcol=T)
  pressures[, max := max(pressure), by=.id]
  deepest_dip = pressures[max(max) == max]

  v_lat = rbindlist(lapply(session$data , function(x) `@`( x , metadata)["latitude"]), idcol=T)
  v_lon = rbindlist(lapply(session$data , function(x) `@`( x , metadata)["longitude"]), idcol=T)
  v_time = rbindlist(lapply(session$data , function(x) `@`( x , metadata)["startTime"]), idcol=T)
  v_station = rbindlist(lapply(session$data , function(x) `@`( x , metadata)["station"]), idcol=T)

  sb = rbindlist(lapply(session$data , function(x) `@`( x , data)), idcol=T, use.names=T, fill=T)
  sb = merge(v_station, sb, by=".id", all.x = T)

    # covert to array for netcdf
  v_temp = reshape2::acast(sb[,.(pressure, station, temperature)], pressure ~ station)

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
    att.put.nc(nc, "pressure", "instrument", "NC_CHAR", "inst_ctd")
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
    att.put.nc(nc, "time", "comment", "NC_CHAR", "Time from GPS NMEA")
  var.put.nc(nc, "time", as.numeric(v_time$startTime))

  var.def.nc(nc, "lat", "NC_DOUBLE", "profile")
    att.put.nc(nc, "lat", "standard_name", "NC_CHAR", "latitude")
    att.put.nc(nc, "lat", "units", "NC_CHAR", "degrees_north")
    att.put.nc(nc, "lat", "axis", "NC_CHAR", "Y")
    att.put.nc(nc, "lat", "valid_min", "NC_FLOAT", -90)
    att.put.nc(nc, "lat", "valid_max", "NC_FLOAT", 90)
    att.put.nc(nc, "lat", "grid_mapping", "NC_CHAR", "crs")
    att.put.nc(nc, "lat", "comment", "NC_CHAR", "Position from GPS NMEA")
  var.put.nc(nc, "lat", v_lat$latitude)

  var.def.nc(nc, "lon", "NC_DOUBLE", "profile")
    att.put.nc(nc, "lon", "standard_name", "NC_CHAR", "longitude")
    att.put.nc(nc, "lon", "units", "NC_CHAR", "degrees_east")
    att.put.nc(nc, "lon", "axis", "NC_CHAR", "X")
    att.put.nc(nc, "lon", "valid_min", "NC_FLOAT", -180)
    att.put.nc(nc, "lon", "valid_max", "NC_FLOAT", 180)
    att.put.nc(nc, "lon", "grid_mapping", "NC_CHAR", "crs")
    att.put.nc(nc, "lon", "comment", "NC_CHAR", "Position from GPS NMEA")
  var.put.nc(nc, "lon", v_lon$longitude)

  var.def.nc(nc, "crs", "NC_INT", NA) # WGS84
    att.put.nc(nc, "crs", "grid_mapping_name", "NC_CHAR", "latitude_longitude")
    att.put.nc(nc, "crs", "espg_code", "NC_CHAR", "ESPG:4326")
    att.put.nc(nc, "crs", "longitude_of_prime_meridian", "NC_DOUBLE", 0.0)
    att.put.nc(nc, "crs", "semi_major_axis", "NC_DOUBLE", 6378137.0)
    att.put.nc(nc, "crs", "inverse_flattening", "NC_DOUBLE", 298.257223563)

    # data
  var.def.nc(nc, "temperature", "NC_DOUBLE", c("pressure", "profile"))
    att.put.nc(nc, "temperature", "long_name", "NC_CHAR", "Water temperature, IPTS-90" )
    att.put.nc(nc, "temperature", "standard_name", "NC_CHAR", "sea_water_temperature" )
    att.put.nc(nc, "temperature", "units", "NC_CHAR", "degree_Celsius" )
    att.put.nc(nc, "temperature", "_FillValue", "NC_DOUBLE", -99999 )
    att.put.nc(nc, "temperature", "valid_min", "NC_DOUBLE", -5 )
    att.put.nc(nc, "temperature", "valid_max", "NC_DOUBLE", 35 )
    att.put.nc(nc, "temperature", "coordinates", "NC_CHAR", "time lat lon pressure" )
    att.put.nc(nc, "temperature", "instrument", "NC_CHAR", "inst_temp" )
    att.put.nc(nc, "temperature", "grid_mapping", "NC_CHAR", "crs" )
    att.put.nc(nc, "temperature", "coverage_content_type", "NC_CHAR", "physicalMeasurement" )
  var.put.nc(nc, "temperature", v_temp)

  var.def.nc(nc, "conductivity", "NC_DOUBLE", c("pressure", "profile"))
    att.put.nc(nc, "conductivity", "long_name", "NC_CHAR", "Water conductivity" )
    att.put.nc(nc, "conductivity", "standard_name", "NC_CHAR", "sea_water_electrical_conductivity" )
    att.put.nc(nc, "conductivity", "units", "NC_CHAR", "S m-1" )
    att.put.nc(nc, "conductivity", "_FillValue", "NC_DOUBLE", -99999 )
    att.put.nc(nc, "conductivity", "valid_min", "NC_DOUBLE", 0 )
    att.put.nc(nc, "conductivity", "valid_max", "NC_DOUBLE", 7 )
    att.put.nc(nc, "conductivity", "coordinates", "NC_CHAR", "time lat lon pressure" )
    att.put.nc(nc, "conductivity", "instrument", "NC_CHAR", "inst_cond" )
    att.put.nc(nc, "conductivity", "grid_mapping", "NC_CHAR", "crs" )
    att.put.nc(nc, "conductivity", "coverage_content_type", "NC_CHAR", "physicalMeasurement" )
  var.put.nc(nc, "conductivity", v_cond)

  var.def.nc(nc, "salinity", "NC_DOUBLE", c("pressure", "profile"))
    att.put.nc(nc, "salinity", "long_name", "NC_CHAR", "Practical Salinity" )
    att.put.nc(nc, "salinity", "standard_name", "NC_CHAR", "sea_water_practical_salinity" )
    att.put.nc(nc, "salinity", "units", "NC_CHAR", "1e-3" )
    att.put.nc(nc, "salinity", "_FillValue", "NC_DOUBLE", -99999 )
    att.put.nc(nc, "salinity", "valid_min", "NC_DOUBLE", 0 )
    att.put.nc(nc, "salinity", "valid_max", "NC_DOUBLE", 36 )
    att.put.nc(nc, "salinity", "coordinates", "NC_CHAR", "time lat lon pressure" )
    att.put.nc(nc, "salinity", "instrument", "NC_CHAR", "inst_ctd" )
    att.put.nc(nc, "salinity", "grid_mapping", "NC_CHAR", "crs" )
    att.put.nc(nc, "salinity", "coverage_content_type", "NC_CHAR", "physicalMeasurement" )
  var.put.nc(nc, "salinity", v_sal)

  var.def.nc(nc, "altimeter", "NC_DOUBLE", c("pressure", "profile"))
    att.put.nc(nc, "altimeter", "long_name", "NC_CHAR", "height above sea floor (altimeter)" )
    att.put.nc(nc, "altimeter", "standard_name", "NC_CHAR", "height_above_sea_floor" )
    att.put.nc(nc, "altimeter", "units", "NC_CHAR", "m" )
    att.put.nc(nc, "altimeter", "_FillValue", "NC_DOUBLE", -99999 )
    att.put.nc(nc, "altimeter", "valid_min", "NC_DOUBLE", 1)
    att.put.nc(nc, "altimeter", "valid_max", "NC_DOUBLE", 100)
    att.put.nc(nc, "altimeter", "coordinates", "NC_CHAR", "time lat lon pressure" )
    att.put.nc(nc, "altimeter", "instrument", "NC_CHAR", "inst_alt" )
    att.put.nc(nc, "altimeter", "grid_mapping", "NC_CHAR", "crs" )
    att.put.nc(nc, "altimeter", "coverage_content_type", "NC_CHAR", "physicalMeasurement" )
  var.put.nc(nc, "altimeter", v_alt)

  var.def.nc(nc, "turbidity", "NC_DOUBLE", c("pressure", "profile"))
    att.put.nc(nc, "turbidity", "long_name", "NC_CHAR", "Turbidity expressed as NTU/FTU (Nephelometric Turbidity Units)" )
    att.put.nc(nc, "turbidity", "standard_name", "NC_CHAR", "sea_water_turbidity" )
    att.put.nc(nc, "turbidity", "units", "NC_CHAR", "1" )
    att.put.nc(nc, "turbidity", "_FillValue", "NC_DOUBLE", -99999 )
    att.put.nc(nc, "turbidity", "valid_min", "NC_DOUBLE", 0.01)
    att.put.nc(nc, "turbidity", "valid_max", "NC_DOUBLE", 750)
    att.put.nc(nc, "turbidity", "coordinates", "NC_CHAR", "time lat lon pressure" )
    att.put.nc(nc, "turbidity", "instrument", "NC_CHAR", "inst_turb" )
    att.put.nc(nc, "turbidity", "grid_mapping", "NC_CHAR", "crs" )
    att.put.nc(nc, "turbidity", "coverage_content_type", "NC_CHAR", "physicalMeasurement" )
  var.put.nc(nc, "turbidity", v_turb)

  var.def.nc(nc, "par", "NC_DOUBLE", c("pressure", "profile"))
    att.put.nc(nc, "par", "long_name", "NC_CHAR", "Photosynthetically Active Radiation (PAR)" )
    att.put.nc(nc, "par", "standard_name", "NC_CHAR", "downwelling_photosynthetic_photon_flux_in_sea_water" )
    att.put.nc(nc, "par", "units", "NC_CHAR", "mol m-2 s-1" )
    att.put.nc(nc, "par", "_FillValue", "NC_DOUBLE", -99999 )
    att.put.nc(nc, "par", "valid_min", "NC_DOUBLE", 0.0001)
    att.put.nc(nc, "par", "valid_max", "NC_DOUBLE", 10)
    att.put.nc(nc, "par", "coordinates", "NC_CHAR", "time lat lon pressure" )
    att.put.nc(nc, "par", "instrument", "NC_CHAR", "inst_par" )
    att.put.nc(nc, "par", "grid_mapping", "NC_CHAR", "crs" )
    att.put.nc(nc, "par", "coverage_content_type", "NC_CHAR", "physicalMeasurement" )
  var.put.nc(nc, "par", v_par*1000) # convert to mol from umol

  # instruments
var.def.nc(nc, "inst_temp", "NC_BYTE", c("profile"))
  att.put.nc(nc, "inst_temp", "long_name", "NC_CHAR", "pumped temperature sensor")
  att.put.nc(nc, "inst_temp", "nodc_name", "NC_CHAR", "CTD")
  att.put.nc(nc, "inst_temp", "make_model", "NC_CHAR", "SBE 3plus" )
  att.put.nc(nc, "inst_temp", "precision", "NC_CHAR", "0.0002" )

var.def.nc(nc, "inst_cond", "NC_BYTE", c("profile"))
  att.put.nc(nc, "inst_cond", "long_name", "NC_CHAR", "pumped conductivity cell")
  att.put.nc(nc, "inst_cond", "nodc_name", "NC_CHAR", "CTD")
  att.put.nc(nc, "inst_cond", "make_model", "NC_CHAR", "SBE 4C" )
  att.put.nc(nc, "inst_cond", "precision", "NC_CHAR", "0.00004" )

var.def.nc(nc, "inst_ctd", "NC_BYTE", c("profile"))
  att.put.nc(nc, "inst_ctd", "long_name", "NC_CHAR", "pumped CTD")
  att.put.nc(nc, "inst_ctd", "nodc_name", "NC_CHAR", "CTD")
  att.put.nc(nc, "inst_ctd", "make_model", "NC_CHAR", "SBE 9plus")
  att.put.nc(nc, "inst_ctd", "precision", "NC_CHAR", "0.63 m" )

var.def.nc(nc, "inst_alt", "NC_BYTE", c("profile"))
  att.put.nc(nc, "inst_alt", "long_name", "NC_CHAR", "Sonar altimeter")
  att.put.nc(nc, "inst_alt", "nodc_name", "NC_CHAR", "altimeter")
  att.put.nc(nc, "inst_alt", "make_model", "NC_CHAR", "Teledyne/Benthos PSA-916")
  att.put.nc(nc, "inst_alt", "precision", "NC_CHAR", "0.1" )

var.def.nc(nc, "inst_turb", "NC_BYTE", c("profile"))
  att.put.nc(nc, "inst_turb", "long_name", "NC_CHAR", "Optical backscatter")
  att.put.nc(nc, "inst_turb", "nodc_name", "NC_CHAR", "turbidity meter")
  att.put.nc(nc, "inst_turb", "make_model", "NC_CHAR", "Seapoint STM")
  att.put.nc(nc, "inst_turb", "precision", "NC_CHAR", "0.1" )

var.def.nc(nc, "inst_par", "NC_BYTE", c("profile"))
  att.put.nc(nc, "inst_par", "long_name", "NC_CHAR", "Photosynthetically Active Radiation Sensor")
  att.put.nc(nc, "inst_par", "nodc_name", "NC_CHAR", "PAR sensor")
  att.put.nc(nc, "inst_par", "make_model", "NC_CHAR", "Licor Li-192 with Cefas low light amplifier")


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

