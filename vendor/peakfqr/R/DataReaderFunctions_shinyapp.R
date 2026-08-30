#############          Data reader functions   ##############

########   Function to format peak_flow data from RDB tab-delim or NWIS
#' Format Peak Flow Data
#' @description
#' Function to format peak flow data from RDB delimited file or NWIS query response.
#' @param pf_data_df A dataframe of peak flow data in RDB format
#'
#' @return pf_data   A re-formatted dataframe of peak flow data 
#' @keywords internal
#' 
#' @import dplyr
#' @importFrom dataRetrieval calcWaterYear
#' @importFrom stringr str_pad
#' @importFrom rlang .data
#'
format_pf_data <- function(pf_data_df){
  
  #If the first row is not data, drop (applies for RDB files with column definition row) - https://waterdata.usgs.gov/nwis/?tab_delimited_format_info
  if(pf_data_df[1,1] == "5s"){
      pf_data_df <- pf_data_df[-1,]
  }
  
  pf_data <- pf_data_df %>%
    select(.data$site_no, .data$peak_dt, .data$peak_va, .data$peak_cd, .data$gage_ht, .data$gage_ht_cd, .data$year_last_pk) %>% #retain columns as in Watstore file, plus 'year_last_pk' for code handling
    #Configure the data - Add broken out date components and water year
    mutate(peak_va = as.numeric(.data$peak_va), 
           gage_ht = as.numeric(.data$gage_ht), 
           year_last_pk = as.numeric(.data$year_last_pk),
           peak_year = as.numeric(substr(.data$peak_dt, 1, 4)), 
           peak_month = as.numeric(substr(.data$peak_dt, 6, 7)), 
           peak_day = as.numeric(substr(.data$peak_dt, 9, 10))) %>% 
    mutate(peak_month = ifelse(.data$peak_month == 0, 1, .data$peak_month), #provide a number for missing months - January will place it in the same water year as the calendar year. 
           peak_day = ifelse(.data$peak_day == 0, 1, .data$peak_day),  #provide a number for missing days - will not affect water year.
           peak_WY = calcWaterYear(as.character(paste0(.data$peak_year,"-", str_pad(.data$peak_month, width = 2, pad = 0),"-", str_pad(.data$peak_day, width = 2, pad = 0)), format = '%Y-%m-%d')))
  
  return(pf_data)
  
}

#########  Function to read delimited peak flow data
#' Read Peak Flow Delimited
#' @description 
#' Function to read RDB/delimited peak flow data. See the following link for information on 
#' peak flow data formats: 
#' <a href="https://help.waterdata.usgs.gov/output-formats#peak_rdb">https://help.waterdata.usgs.gov</a>
#' @param file_path A filepath to peak flow data in RDB/tab delimited format 
#'
#' @return A formatted dataframe of peak flow data 
#' @keywords internal
#'
read_pf_delim <- function(file_path){
  
  #Try to read in and format data from tab-delimited format
  z <- try(
          {read.delim(file_path, sep = "\t", comment.char = "#", na.strings = "", colClasses = "character") %>%
            format_pf_data()}, 
            silent = TRUE)
  
  #Error handling for bad file load 
  if(is(z, 'try-error') | is(z, 'error')) {
    stop("Please ensure that the selected file is peak flow data in tab-delimited format. Note that WATSTORE format files are not accepted.", call. = FALSE)
  } else { 
    #Store the data read in from the tab-delimited file. 
    pf_data <- z
  }
  
  return(pf_data)
  
}

########  Function to read NWIS data for sites
#' Read NWIS
#' @description
#' Read and format NWIS peak flow data for a specified list of sites. 
#' @param station_list A character string of NWIS site IDs, separated by 
#' commas. Site numbers should be between 8 and 15 digits long. 
#'
#' @return Returns a list containing the following: 
#'  
#' 1. pf_data: a formatted dataframe of peak flow data 
#' 
#' 2. inval_sites: a list of invalid site numbers that did not return data 
#' 
#' 3. nodata_sites: a list of site numbers for which no peak flow data is available
#' 
#' 4. dupyear_sites: a list of site numbers for which peak flow data include duplicate records for the same water year. These sites are not loaded with the data. 
#' @keywords internal
#'
#' 
#' @import dplyr
#' @importFrom dataRetrieval readNWISpeak
#' @importFrom dataRetrieval whatNWISsites
#' 
read_pf_nwis <- function(station_list = ""){
  
  e <- try(readNWISpeak(station_list, asDateTime = FALSE, convertType = FALSE), silent = TRUE)
  
  #Vector to save out invalid stations
  inval_sites <- c()
  #Vector to save out okay stations 
  retry_list <- c()
  #vector to save out stations with multiple peak flows for the same water year
  dupyear_sites <- c()
  
  #If an invalid site number generates an error, figure out which site number(s) contributed. If there are others, still run with remaining sites. 
  if(is(e, 'try-error') | is(e, 'error')){
    
    #Break down station list
    stations <- unlist(strsplit(station_list, ", "))
    
    #If there's an invalid site number, figure out which one(s)
    for(i in stations){
      
      y <- try(whatNWISsites(sites = i), silent = TRUE)
      
      if(is.null(y)){
        inval_sites <- inval_sites %>% append(i)
      }else{
        retry_list <- retry_list %>% append(i)
      }
    }
    
    #Give a warning of invalid sites numbers
    warning(paste0("No data could be retrieved for the following site(s) due to an invalid site number: ", toString(inval_sites)), 
            call. = FALSE)
    
    #Figure out if there are any valid sites. Rerun call if so.
    if(length(retry_list) > 0){
      retry_station_list <- toString(retry_list)
      nwis_response <- readNWISpeak(retry_station_list, asDateTime = FALSE, convertType = FALSE)
    }else{ # If there are no valid sites, throw error. 
      stop("Error in NWIS query. Check site number(s).", call. = FALSE)
    }
  }else{
    nwis_response <- e
  }
  
  #Handle if no sites returned peak flow data (ALL missing - individual sites in warning below)
  if(nrow(nwis_response) == 0){
    stop("No peak flow data for requested NWIS site(s).", call. = FALSE)
  }
  
  #Get site info from attributes
  pf_site_info <- attributes(nwis_response)$siteInfo
  
  if(nrow(pf_site_info) != 0){
    site_metadata <- pf_site_info[, c("site_no", "station_nm", "dec_lat_va", "dec_long_va")]
    colnames(site_metadata) <- c("site_no", "station_nm", "Latitude", "Longitude")
  } else {
    site_metadata <- matrix(ncol = 4) %>% data.frame()
    names(site_metadata) <- c("site_no", "station_nm", "Latitude", "Longitude") #blank dataframe with the right structure
  }
  
  #Process & reformat peak flow data returned by NWIS
  pf_data <- nwis_response %>%
             format_pf_data()
  
  #Remove any sites with multiple peak flows for the same WY. 
  #(Yes, this occurs in NWIS as of 2023)
  pf_data <- pf_data[!duplicated(pf_data[, c("site_no", "peak_WY")]),]
  
    #Safe fail if all sites removed
    if(nrow(pf_data) == 0){
      stop("All queried sites have been excluded because they report multiple peak flow values for the same water year.", call. = FALSE)
    }
    
    #Warn if only some sites removed.
    if(length(dupyear_sites) > 0){
      warning(paste0("Following site(s) excluded because they report multiple peak flow values for the same water year: ", toString(dupyear_sites)),
              call. = FALSE)
    }
    
  #Warnings if any of the entered sites did not return data
  request_sites <- trimws(unlist(strsplit(station_list, ",")))
  actual_sites <- pf_data$site_no %>% unique()
  missing_sites <- request_sites[request_sites %!in% actual_sites]
  nodata_sites <- missing_sites[missing_sites %!in% c(inval_sites, dupyear_sites)]
  
  if(length(nodata_sites) > 0){
    
    warning(paste0("No peak flow data are available for for the following site(s): ", toString(nodata_sites)), call. = FALSE)
    
  }
  
  #Add metadata to pf_data as attribute
  attributes(pf_data)$site_info <- site_metadata
  
  list(pf_data = pf_data, inval_sites = toString(inval_sites), nodata_sites = toString(nodata_sites), dupyear_sites = toString(dupyear_sites))
  
}

########  Function to input peak flow data as file path 
#' Read Peak Flow Data from File Path
#' @description
#' Read peak flow data from a file path poiting to an RDB or watstore file or a folder containing specifications and peak flow data.  
#' @param data_path Filepath to either a file or a folder containing peak flow data
#'
#' @return A formatted dataframe of peak flow data 
#' @keywords internal
#'
#'
pf_filepath_load <- function(data_path) {
  
  if(dir.exists(data_path)){ #If the upload is a folder, rather than a file: Folder with subfolder. Pull files from each subfolder
    
    # "Open/read multiple sites in a folder structure with one site per folder"
    pf_data <- data.frame(matrix(ncol = 10))
    names(pf_data) <- c("site_no","peak_dt","peak_va","peak_cd","gage_ht","gage_ht_cd","peak_year","peak_month","peak_day","peak_WY" )
    
    for (folder in dir(data_path)){
      sel_folder <- paste0(data_path, folder)
      sel_file <- paste0(sel_folder, "/", dir(sel_folder)) #Assuming the folder contains only the data file. Discuss need for further handling. 
      sel_pf_data <- read_pf_delim(sel_file)
      pf_data <- sel_pf_data %>% rbind(pf_data) #Add this data onto dataframe
    }
    
  }else{ #Handle file inputs - tab delimited or Watstore 
    
    loadError=F
    a=try({pf_data <- readWATSTOREAll(data_path)})
    loadError <- (is(a, 'try-error')|is(a,'error'))
    if(loadError==T){ 
      pf_data <- read_pf_delim(data_path)
    } else if(loadError==F) {
      pf_data <- data.frame(matrix(nrow = 0, ncol = 10))
      names(pf_data) <- c("site_no","peak_dt","peak_va","peak_cd","gage_ht","gage_ht_cd","peak_year","peak_month","peak_day","peak_WY" )
    }
  }
  return(pf_data)
}  

########  Function to read site info 
#' Read Site Info
#' @description
#' Reads in site information from a provided file path 
#' @param data_path A path to a tab delimited file containing site information 
#' including site number, station name, latitude, and longitude. 
#'
#' @return A dataframe of site number, station name, latitude, and longitude 
#' @keywords internal
#' 
#'
#'
read_site_info <- function(data_path) {
  
  site_data_all <- read.delim(data_path, comment.char = "#")

  #If the first row is not data, drop (applies for RDB files with column definition row) - https://waterdata.usgs.gov/nwis/?tab_delimited_format_info
  if(site_data_all[1,1] == "5s"){
    site_data_all <- site_data_all[-1,]
  }
  
  #Catch and fail if bad format.
  t <- try(site_data_all <- site_data_all[, c("site_no", "station_nm", "dec_lat_va", "dec_long_va")])
  
  if(is(t, "error") | is(t, "try-error")){
    stop("Site information could not be loaded. Please check the format of the selected file.", call. = FALSE)
  }
  
  colnames(t) <- c("site_no", "station_nm", "Latitude", "Longitude")

  return(t)
  
}

#######  Function to read psf linked file
#' Read PSF Linked File
#' @description
#' Read data file specified in a PSF file. 
#' Tries to read in and format data from WATSTORE or tab-delimited format.
#' 
#' @param pf_datapath Path to peak flow data from specification file.
#'
#' @return A formatted dataframe of peak flow data 
#' @keywords internal
#'
#' 
read_psf_linked_file <- function(pf_datapath) {
  #Try to read in and format data from WATSTORE or tab-delimited format
  z <- try({
    tryCatch(readWATSTOREAll(pf_datapath), 
             error = function(ex){
               read_pf_delim(pf_datapath)
             })},
    silent = TRUE)
  
  #Error handling for bad file load 
  if(is(z, 'try-error') | is(z, 'error')) {
    stop("The data file associated with this specification file could not be loaded.
           Please ensure that the data file still has the same name and filepath indicated by the specification file.", call. = FALSE)
  } 
  else if(nrow(z) == 0){
    stop("The data file associated with this specification file contains no valid data. 
           Please ensure that the data file is correctly formatted.", call. = FALSE)
    
  }
  
  else { 
    #Store the data read in from the tab-delimited file. 
    pf_data <- z
  }
  
  return(pf_data)
  
}   

######## Function to get site data 
#' Get Site Data
#' @description
#' Create a summary of peak flow data by site 
#' @param pf_data A formatted dataframe of peak flow data  
#'
#' @return Returns a dataframe with of the beginning year, end year, and historical record length of peak flow data for each site. 
#' 
#' @keywords internal
#' 
#' @import dplyr
#' @importFrom rlang .data
#'
#' 
get_site_data <- function(pf_data) {
  
  site_data <- pf_data %>% 
    group_by(.data$site_no) %>% 
    summarise(BegYear = min(.data$peak_WY), #Based on 8/21 input from  Seth - "populated from first... year of input data for the site"
              EndYear = 1, max(.data$peak_WY), #Based on 8/21 input from  Seth - "populated from... last year of input data for the site"
              HistPeriod = as.integer(sum(!is.na(.data$peak_va)))) %>% #Based on 8/21 input from  Seth - "populated from the record length computed as the number of rows of data for the site (with valid peak discharge values, some may only have peak gage height, which we don't want to count)." 
    select(.data$site_no, .data$HistPeriod, .data$BegYear, .data$EndYear)
  
  return(site_data)
  
}


######## Function to retrieve skew value and std err based on coordinates
#' Read Map Skew Value
#' @description
#' Retrieve map skew value and standard error based on latitude and longitude. 
#' References a shapefile of published skew coefficients and standard errors, and, where appropriate, rasters with skew coefficients.
#' 
#' @param latitude Latitude for site in decimal degrees
#' @param longitude Longitude for site in decimal degrees
#'
#' @return Dataframe with the following columns: 
#' GenSkew - the regional skew coefficient for the provided coordinates
#' SkewSE - the standard error for the provided coordinates
#' Citation - citation for the regional skew value
#' ShortName - short name for the citation, typically an SIR number
#' 
#' @details
#' If a regional skew coefficient is not available for the input coordinates, 
#' a dataframe of NA values will be returned. 
#' 
#' See the Regional Skew vignette for details on the areas for which 
#' regional skew coefficients are available:
#' \code{vignette("RegionalSkew", package = "peakfq")}
#' 
#' @examples
#' #Get regional skew for Rio Grande at Embudo, NM
#' latitude <- 	36.20555556
#' longitude <- -105.9639722
#' 
#' d <- read_skew_map_pt(36.20555556, -105.9639722)
#' 
#' regional_skew_coefficient <- d$GenSkew
#' regional_skew_Standard_Error <- d$SkewSE
#' regional_skew_source <- d$MapSkewSource
#' 
#' 
#' @export
#'
#' @import sf
#' @importFrom raster raster
#' @importFrom raster extract
#' 
read_skew_map_pt <- function(latitude, longitude){
  
  #Initialize return values to NA
  map_val_skew <- NA
  map_val_std_err <- NA
  map_val_cite <- NA
  map_val_name <- NA
  
  #Ensure a valid latitude and longitude have been provided. If not, alert and safe fail.
  if(!is.numeric(latitude) | !is.numeric(longitude)){
    stop(paste0("Invalid coordinates: ", longitude, ",", latitude))
  }
  
  #Reference map data
  regSkew_shp <- sf::read_sf(system.file("extdata", "BGLS_Skewmap/RegSkewPolys_PeakFQ_G100.shp", package = "peakfq"))

  
  #Look up skew values from shapefile
  input_proj <- "+proj=longlat +datum=WGS84 +no_defs" #input crs for coordinates
  
  pt_coords <- cbind(longitude, latitude) %>% 
    as.data.frame() %>% 
    sf::st_as_sf(coords = c("longitude", "latitude")) %>%
    sf::st_set_crs(input_proj) %>%
    sf::st_transform(sf::st_crs(regSkew_shp)) #ensure point layer in same crs as shapefile of skew values
  
  pt_val <- sf::st_intersection(regSkew_shp, pt_coords) %>% suppressWarnings()
  
  #Ensure an intersection was returned. If not, alert & safe fail.
  if(nrow(pt_val) == 0){
    warning("The coordinates provided fall outside the bounds of the skew map.")
    return(data.frame("GenSkew" = NA, "SkewSE" = NA, "Citation" = NA, "ShortName" = NA)) 
  }
  
  #Fail if map returns multiple skew values
  if(nrow(pt_val) > 1){
    warning(paste0("Multiple values in skew map for location ",
                   longitude, ",", latitude))
    return(data.frame("GenSkew" = NA, "SkewSE" = NA, "Citation" = NA, "ShortName" = NA)) 
  }
  
  map_val_cite <- pt_val$Citation
  map_val_name <- pt_val$Short_Name
  
  #Find is skew value is constant for polygon
  # 1 indicates constant skew for polygon
  # 0 indicates variable skew published in raster 
  # -99 indicates skew values are unavailable (study in progress or California EQ)
  map_val_const <- pt_val$Valu_Const
  
  if(map_val_const == 1){
    #Find corresponding skew value from shapefile
    map_val_skew <- pt_val$Skew_value 
    
    #Find corresponding std. error from shapefile
    map_val_std_err <- pt_val$Std_Err 
    
  }
  
  
  # If indicated by shapefile, reference .tif:
  # "If the Skew_val field is -99, check if the Tiff_file field is populated. 
  # If so, look up the skew value from the appropriate raster from the Tiff_file field 
  # and take the standard error from the Std_Err field. 
  else if(map_val_const == -99){
    
    #Find corresponding std. error from shapefile
    map_val_std_err <- pt_val$Std_Err
    
    #Get scale for the tiff skew raster from shapefile
    map_val_scale <- pt_val$Tiff_Scale
    
    tiff <- pt_val$Tiff_file
    
    tif_path <- system.file("extdata", paste0("skew_tif/", pt_val$Tiff_file),
                            package = "peakfq")
    
    #Fail if skew raster file doesn't exist
    if(!file.exists(tif_path)){
      warning(paste0("No raster file available for regional skew for location ",
                     longitude, ",", latitude))
      return(data.frame("GenSkew" = NA, "SkewSE" = NA, "Citation" = NA, "ShortName" = NA)) 
    }
    
    sel_tif <- raster::raster(tif_path)
    
    #Transform coordinates to raster coordinate system
    pt_coords <- pt_coords %>% sf::st_transform(sf::st_crs(sel_tif))
    
    #Extract raster value and divide by scale
    map_val_skew <- raster::extract(sel_tif, pt_coords)/map_val_scale
    

  }
  
  return(data.frame("GenSkew" = as.numeric(map_val_skew), "SkewSE" = as.numeric(map_val_std_err), "MapSkewSource" = as.character(map_val_cite), "MapSkewSourceText" = as.character(map_val_name))) 
}

######## Function to retrieve skew value and std err based on coordinates in station specifications dataframe
#' Read Map Skew Value
#' @description
#' Retrieve map skew value and standard error based on location coordinates in station specifications dataframe.
#' Utilizes map skew option and coordinates from station specifications dataframe and 
#' references USGS shapefile of skew values and standard errors, and, where appropriate, tiffs with skew values.
#' @param station_specs_df A dataframe of station specifications, created with function "generate_station_specifications"
#'
#' @return Updated station specifications dataframe containing skew and standard error values from map, where available. 
#' @keywords internal
#' 
#' @import dplyr
#'
#' 
get_map_skew <- function(station_specs_df){

  missing_lat_long <- c()
  missing_skew_val <- c()
  
  #Get sites to retrieve map skew for
  map_sites <- station_specs_df[!is.na(station_specs_df$MapSkew) & (station_specs_df$MapSkew == TRUE),]
  missing_lat_long <- map_sites$site_no[is.na(map_sites$Latitude) | is.na(map_sites$Longitude)]
  
  
  if(length(missing_lat_long) > 0){
    warning(paste0("A latitude and longitude must be provided to utilize the map skew. Please provide site information including the latitude and longitude or enter a skew value manually for site(s): ",
                   paste0(missing_lat_long, collapse=", ")))
  }
  
  #Remove sites with bad lat/longs
  map_sites <- map_sites[map_sites$site_no %!in% missing_lat_long,]
  
  if(nrow(map_sites) > 0){
    map_output <- mapply(read_skew_map_pt, map_sites$Latitude, map_sites$Longitude, SIMPLIFY = FALSE)
    map_output <- dplyr::bind_rows(map_output) #Convert to single dataframe
    map_output$site_no <- map_sites$site_no
    
    missing_skew_val <- map_output$site_no[is.na(map_output$GenSkew) | is.na(map_output$SkewSE)]
    
    #Make sure nothing went wrong matching site numbers and add to dataframe
    if(all(map_output$site_no == map_sites$site_no)){
      station_specs_df$GenSkew[station_specs_df$site_no %in% map_sites$site_no] <- map_output$GenSkew
      station_specs_df$SkewSE[station_specs_df$site_no %in% map_sites$site_no] <- map_output$SkewSE
      station_specs_df$MapSkewSource[station_specs_df$site_no %in% map_sites$site_no] <- map_output$MapSkewSource
      station_specs_df$MapSkewSourceText[station_specs_df$site_no %in% map_sites$site_no] <- map_output$MapSkewSourceText
      
    }
    else{
      stop("Non-matching site numbers when retrieving map skew")
    }
    
  }
  
  
  return(list(station_specs_df = station_specs_df, missing_lat_long = missing_lat_long, missing_skew_val = missing_skew_val))
}

