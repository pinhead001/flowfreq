#' Read WATSTORE
#'
#' Function to read peak discharge and gage height from a WATSTORE file into an R dataframe.
#' Other fields supported in the legacy WATSTORE format will be ignored.
#'
#' @param filepath character - path to WATSTORE file
#' @param convertType logical - whether or not to convert dates to R date type and peak discharge and gage height values to numeric type. If false, these columns are returned as text.
#'
#'
#' @return A dataframe containing peak flow and gage height information in the WATSTORE file.
#'
#'
#' @keywords internal
#'
#' @importFrom utils read.fwf
#' @importFrom stats complete.cases
readWATSTORE <- function(filepath, convertType=TRUE){

  #Modified from original code by Roy Sando

  if(!file.exists(filepath)){
    stop(paste("Specified WATSTORE file does not exist: ", filepath))
  }

  #Read in watstore file. Can include multiple gages.
  watstore <- utils::read.fwf(filepath, widths = c(16,8,7,12,8,4,4)) #Read based on the specification in the PeakFQ documentation. Fields on the right will be ignored later.

  colnames(watstore) <- c("site_no", "peak_dt", "peak_va", "peak_cd", "gage_ht", "gage_ht_cd", "year_last_pk")#Rename columns in watstore file
  #Remove non-peak-data rows
  H <- grep( "H", watstore[,"site_no"])
  Y <- grep( "Y", watstore[,"site_no"])
  N <- grep( "N", watstore[,"site_no"])
  Z <- grep( "Z", watstore[,"site_no"])

  watstore <- watstore[-c(H, Y, N, Z), ]

  #Remove leading 3 from site number
  watstore$site_no <- trimws(substring(watstore$site_no, 2))

  #Clean up codes and values
  watstore$peak_dt <- as.character(watstore$peak_dt)
  watstore$peak_va <- as.character(watstore[, "peak_va"])
  watstore$peak_cd <- trimws(as.character(watstore[, "peak_cd"]))
  watstore$gage_ht_cd <- trimws(as.character(watstore[, "gage_ht_cd"]))
  watstore$gage_ht <- as.character(watstore[, "gage_ht"])
  watstore$year_last_pk <- as.character(trimws(watstore[, "year_last_pk"]))


  #Add columns for month, year, and day
  watstore$peak_year <-  as.numeric(substring(watstore$peak_dt, 1, 4))
  watstore$peak_month <-  as.numeric(substring(watstore$peak_dt, 5, 6))
  watstore$peak_day <-  as.numeric(substring(watstore$peak_dt, 7, 8))

  if (NA %in% watstore$peak_month){
    message("One or more entries in WATSTORE file is missing a month")
  }
  if (NA %in% watstore$peak_day){
    message("One or more entries in WATSTORE file is missing a day")
  }

  if(convertType){
    #Convert to numeric if specified
    watstore$peak_dt <- as.Date(as.character(watstore$peak_dt), "%Y%m%d")
    watstore$peak_va <- as.numeric(watstore$peak_va)
    watstore$gage_ht <- as.numeric(watstore$gage_ht)
    watstore$year_last_pk <- as.integer(watstore$year_last_pk)
  }





  #Compute water year of peak
  #For peaks with a valid month before September, year in WATSTORE is the water year
  watstore[!is.na(watstore$peak_month) & watstore$peak_month < 10,"peak_WY"] <- watstore[!is.na(watstore$peak_month) & watstore$peak_month < 10,"peak_year"]
  #For peaks with a valid month after September, year in WATSTORE is the year prior to the water year
  watstore[!is.na(watstore$peak_month) & watstore$peak_month >= 10,"peak_WY"] <- watstore[!is.na(watstore$peak_month) & watstore$peak_month >= 10,"peak_year"] + 1
  watstore[is.na(watstore$peak_month),"peak_WY"] <- watstore[is.na(watstore$peak_month),"peak_year"] #For peaks without a month, the year in WATSTORE file is the water year

  #Remove blank rows
  watstore <- watstore[stats::complete.cases(watstore$site_no),]
  watstore <- watstore[trimws(watstore$site_no) != "",]

  #Separate codes with a comma
  watstore$peak_cd <- separateCodes(watstore$peak_cd)
  watstore$gage_ht_cd <- separateCodes(watstore$gage_ht_cd)

  watstore$peak_cd[watstore$peak_cd == ""] <- NA #Replace blanks with NA to match dataRetrieval package output
  watstore$gage_ht_cd[watstore$gage_ht_cd == ""] <- NA #Replace blanks with NA to match dataRetrieval package output


  if (nrow(watstore) != nrow(unique(watstore))){
    warning("WATSTORE file has duplicate rows")
    watstore <- unique(watstore)
  }


  return(watstore)

}

#' Read WATSTORE Info
#'
#' Function to read header information from a WATSTORE file into an R data frame
#'
#' @param filepath character - path to WATSTORE file.
#'
#'
#' @return A dataframe containing header information from the WATSTORE file. Reading aquifer type codes is not currently supported.
#'
#' @keywords internal
#'
#'
#' @importFrom utils read.fwf
#' @importFrom stats complete.cases
readWATSTOREinfo <- function(filepath){
  
  if(!file.exists(filepath)){
    stop(paste("Specified WATSTORE file does not exist: ", filepath))
  }
  
  #Read in watstore file. Can include multiple gages.
  watstoreNames <- utils::read.fwf(filepath, widths = c(16,48)) #Read correct widths for station name
  
  colnames(watstoreNames) <- c("site_no", "station_nm") #Rename columns in watstore file
  
  watstoreNames <- watstoreNames[grep("N", watstoreNames$site_no), ] #Select only rows with station names
  
  #Remove leading N from site number
  watstoreNames$site_no <- trimws(substring(watstoreNames$site_no, 2))
  
  #Clean up names
  watstoreNames$station_nm <- as.character(trimws(watstoreNames$station_nm))
  
  
  watstoreInfo <- utils::read.fwf(filepath, widths = c(16,6,7,2, 2, 2, 3, 2, 8, 7, 7, 8, 9)) #Read correct widths for station name
  
  colnames(watstoreInfo) <- c("site_no", "lat_va", "long_va", "seq_num", "state_cd", "district_cd", "county_cd", "site_tp_cd", "huc_cd", "drain_area_va", "contrib_drain_area_va", "alt_va", "well_depth_va") #Rename columns in watstore file
  
  watstoreInfo <- watstoreInfo[grep("H", watstoreInfo$site_no), ] #Select only rows with station headers
  
  watstoreInfo$site_no <- trimws(substring(watstoreInfo$site_no, 2)) #Remove leading N from site number
  
  watstoreInfo <- merge(watstoreInfo, watstoreNames, by="site_no")
  
  watstoreInfo <- watstoreInfo[,c("site_no", "station_nm", "site_tp_cd", "lat_va", "long_va", "district_cd", "state_cd", "county_cd", "alt_va", "huc_cd", "drain_area_va", "contrib_drain_area_va", "well_depth_va")]
  
  #Do a bunch of type conversions
  watstoreInfo$site_tp_cd <- as.character(watstoreInfo$site_tp_cd)
  watstoreInfo$lat_va <- as.numeric(as.character(watstoreInfo$lat_va))
  watstoreInfo$long_va <- as.numeric(as.character(watstoreInfo$long_va))
  watstoreInfo$district_cd <- as.character(watstoreInfo$district_cd)
  watstoreInfo$state_cd <- as.character(watstoreInfo$state_cd)
  watstoreInfo$county_cd <- as.character(watstoreInfo$county_cd)
  watstoreInfo$alt_va <- as.numeric(as.character(watstoreInfo$alt_va))
  watstoreInfo$drain_area_va <- as.numeric(as.character(watstoreInfo$drain_area_va))
  watstoreInfo$contrib_drain_area_va <- as.numeric(as.character(watstoreInfo$contrib_drain_area_va))
  watstoreInfo$well_depth_va  <- as.numeric(as.character(watstoreInfo$well_depth_va ))
  
  
  
  #Remove blank rows
  watstoreInfo <- watstoreInfo[stats::complete.cases(watstoreInfo$site_no),]
  watstoreInfo <- watstoreInfo[trimws(watstoreInfo$site_no) != "",]
  
  
  if (nrow(watstoreInfo) != nrow(watstoreInfo)){
    warning("WATSTORE file has duplicate rows")
    watstoreInfo <- unique(watstoreInfo)
  }
  
  
  return(watstoreInfo)
  
}

#' Read WATSTORE Peak Data and Site Information
#'
#' Function to read peak discharge and gage height from a WATSTORE file into an R dataframe, 
#' with site information returned as a dataframe attribute. 
#'
#' @param filepath character - path to WATSTORE file
#' @param convertType logical - whether or not to convert dates to R date type and peak discharge and gage height values to numeric type. If false, these columns are returned as text.
#'
#'
#' @return A dataframe containing peak flow and gage height information in the WATSTORE file. Site information from the WATSTORE file is returned as a dataframe attribute. 
#' 
#' 
#' @examples
#' \dontrun{
#' 
#' WSfile <- "./path/to/old_WATSTORE_file.txt" #Location of WATSTORE file
#' 
#' specs <- readWATSTOREAll(WSfile) #Read WATSTORE data to R
#' 
#' }
#'
#' @export
#'  
readWATSTOREAll <- function(filepath, convertType=TRUE){
  
  WSdata <- readWATSTORE(filepath, convertType)
  WSinfo <- readWATSTOREinfo(filepath)
  
  attributes(WSdata)$siteInfo <- WSinfo
  
  return(WSdata)
  
  
}


#' Read RDB peak data file
#'
#' Function to read peak streamflow and gage height data in RDB format from NWIS
#'
#' @param file_path character - path to data file
#'
#' @return dataframe containing peak flow data
#' 
#' @details
#' See [NWIS](https://help.waterdata.usgs.gov/faq/about-tab-delimited-output#:~:text=The%20principal%20machine-readable%20data%20format%20supported%20by%20the,header%20section%20containing%20zero%20or%20more%20comment%20lines)
#' documentation for description of the input file format. 
#' 
#' The returned dataframe will include all columns in the input tab-delimited
#' file as well as additional columns `peak_year`, `peak_month`, `peak_day`, and
#' `peak_WY` with the calendar year, month, day, and water year of the annual 
#' peak flow values respectively. 
#' 
#' @examples
#' \dontrun{
#' #Location of RDB data file
#' file_path <- system.file("extdata", "ExampleSites_nwis_peak.txt", package="peakfq")
#' 
#' d <- readRDB(file_path)
#' 
#' } 
#'  
#'
#' @export
#'
readRDB <- function(file_path){
  
  #Read in the data from the tab-delim file
  pf_data <- read.delim(file_path, sep = "\t", comment.char = "#", na.strings = "")
  pf_data <- pf_data[-1,] #drop first row after header - not data.
  
  if("peak_va" %!in% colnames(pf_data)){
    stop("Peak flow data file is missing required column 'peak_va'")
  }
  
  if("peak_dt" %!in% colnames(pf_data)){
    stop("Peak flow data file is missing required column 'peak_dt'")
  }
  
  #Configure the data - Add broken out date components and water year
  pf_data$peak_va <- as.numeric(pf_data$peak_va)
  pf_data$peak_year <- as.numeric(substr(pf_data$peak_dt, 1, 4))
  pf_data$peak_month <- as.numeric(substr(pf_data$peak_dt, 6, 7)) 
  pf_data$peak_day <- as.numeric(substr(pf_data$peak_dt, 9, 10)) 
  pf_data$peak_WY <- ifelse(pf_data$peak_month < 10 , pf_data$peak_year, pf_data$peak_year + 1) #Note, if month is zero (missing) year in data file is the water year
  
  return(pf_data)
  
}

#' Read RDB Site Names
#'
#' Function to read site names from NWIS RDB tab-delimited files
#'
#' @param inFile character - path to NWIS RDB tab-delimited peak flow data file
#'
#' @return Dataframe with the following columns:
#'  1. `site_no` - NWIS site numbers
#'  2. `station_nm` - NWIS site names
#'  
#' @details
#' Function only returns names for sites with an agency code of "USGS"
#' 
#'  
#' @examples
#' \dontrun{
#' #Location of RDB data file
#' file_path <- system.file("extdata", "ExampleSites_nwis_peak.txt", package="peakfq")
#' 
#' site_names <- RDB_station_nm(file_path)
#' 
#' } 
#'
#' @export
#'
RDB_station_nm <- function(inFile){
  
  d <- readLines(inFile)
  d <- d[grepl("# ", d)]
  
  #Get row number to start and data extraction
  startrow <- grep("# Sites in this file include:", d) + 1
  endrow <- grep("# Peak Streamflow-Qualification Codes\\(peak_cd\\):", d) - 1
  
  namerows <- d[startrow:endrow]
  
  #Now, to actually decode this...
  
  if(!all(grepl("#  USGS", namerows))){
    warning("Sites with agencies codes other than USGS were removed")
    namerows <- namerows[grepl("#  USGS", namerows)]
  }
  
  namerows <- trimws(sub("#  USGS", "", namerows))
  
  # We apply an AI generated regular expression to split the site name and number
  rexp <- "^(\\S+) (.*)$"
  n <- data.frame(site_no=sub(rexp,"\\1",namerows), 
                  station_nm=sub(rexp,"\\2",namerows))
  
  return(n)
  
}

#' Separate Peak Codes
#'
#' Function to separate peak flow qualifier codes
#'
#' @param codes character - vector of WATSTORE peak flow codes
#'
#' @return character - codes separated with a comma
#' 
#' @keywords internal
#' 
#' @importFrom stringr str_replace_all
#'
separateCodes <- function(codes){


  codes <- stringr::str_replace_all(codes, c("1" = "1,",
                                             "2" = "2,",
                                             "3" = "3,",
                                             "4" = "4,",
                                             "5" = "5,",
                                             "6" = "6,",
                                             "7" = "7,",
                                             "8" = "8,",
                                             "9" = "9,",
                                             "A" = "A,",
                                             "Bd" = "Bd,",
                                             "Bm" = "Bm,",
                                             "C" = "C,",
                                             "F" = "F,",
                                             "O" = "O,",
                                             "R" = "R,"
  ))

  codes <- sub(",$", "", codes) #Remove last comma

  return(codes)
}


#' Read peakfq Specifications File
#'
#' Reads a peakfq specifications (.psf) file into an R dataframe. 
#' 
#' 
#' 
#' @param PSFfile character -  path to peakfq specifications file
#'
#'
#' @return A list of 4 dataframes: 
#'   1. Dataframe with regional skew information, low-outlier test option,
#'      and urban/regulated peak usage options for each site in the analyis. 
#'   2. Dataframe with perception thresholds
#'   3. Dataframe with flow intervals
#'   4. Dataframe with user-specified peak values
#'
#' The following attributes of the first dataframe contain global analysis options: 
#' - `inFormat` - Format of the input peak flow data, either the USGS RDB format 
#' or legacy WATSTORE format
#' - `inFile` - path to input peak flow data file
#' - `siteinfoFile` - path to tab-delimited site information file
#' - `outFile` - path to the output log file. The name of this file will be used
#'             to determine the names for numeric output files. For example, if
#'             the outFile is specified as ./some/directory/my_project.txt, the
#'             numeric output files will be stored in ./some/directory and named
#'             my_project_xxx.csv, where xxx indicates the type of numeric
#'             output. 
#' - `ConfInterval` - Confidence interval for analysis
#' - `PlotPosition` - Plotting position computed as (m-a)/(N+1-2a) where m is order number, N is 
#' total number of peaks, and a is a parameter where:
#' a = 0.00 for Weibull
#' a = 0.30 for Median/Beard
#' a = 0.375 for Blom
#' a = 0.4 for Cunnane
#' a = 0.44 for Gringorten
#' 
#' - `PlotFormat` - Output image format. Either JPEG, PNG, SVG, or NONE. 
#' - `extended` - logical for use of extended analysis option
#' 
#' @section Specifications File Format: 
#' 
#' 
#' The header fields of the specification have the letters I or O in front. 
#' Each header field may be used only once in a specification file. 
#' Fields starting with I indicate input files. The three valid fields for
#' input files are:
#' - `I ASCI` - path to legacy WATSTORE data file
#' - `I RDB` - path to RDB tab-delimited data file
#' - `I INFO` - path to RDB tab-delimited site information file
#' 
#' Either the `I ASCI` or `I RDB` field is required. These specify the relative path
#' to the input peak-flow data. The `I INFO` field is an optional input for the 
#' relative path to a tab-delimited file with site information. This is only
#' used to input latitude and longitude values for the Shiny app. 
#' 
#' Fields starting with O indicate output options. All output fields are optional. 
#' The valid options for outputs are:
#' - `O Extended` - set to `YES` to run extended analysis, computing annual exceedance
#' probabilities as small as 10^-5. 
#' - `O File` - path to data output files. 
#' - `O ConfInterval` - confidence interval value as decimal (0.90 corresponds
#' to lower 5 percent and upper 95 percent confidence limits)
#' - `O Plot Format` - format of output plots. Valid values are `JPG`, `PNG`, `SVG`, 
#' or `NONE`. 
#' 
#' \itemize{
#' \item{- `O Plot Position` - Plotting position computed as (m-a)/(N+1-2a) where 
#' m is order number, N is  total number of peaks, and a is the plotting position parameter where:
#' \itemize{
#' \item{a = 0.00 for Weibull}
#' \item{a = 0.30 for Median/Beard}
#' \item{a = 0.375 for Blom}
#' \item{a = 0.4 for Cunnane}
#' \item{a = 0.44 for Gringorten}
#' }
#' }
#' 
#' 

#' }
#' 
#' 
#' The station specifications fields are used to set the analysis options for 
#' the one more sites being analyzed. The Station field is used to start
#' specifications for a new station. The identifier following the Station 
#' keyword must be unique and correspond to the identifier in the input 
#' peak-flow data file. The following fields are used to set analysis 
#' options for each station:
#' 
#' 
#' - `PCPT_Thresh` - perception threshold, set as
#'  `PCPT_Thresh StartYear EndYear Lower Upper`
#' - Interval - flow interval, set as
#'  `Interval Year Lower Upper`
#' - Peak - peak flow value, set as
#'  `Peak Year Value`. The `Peak` keyword is considered deprecated; using the
#'  `Interval` keyword with the same lower and upper intervals is preferred. 
#'
#' - `SkewOpt` - skew option. Set to `Station` to use only at-site skew, 
#' `Weighted` to use weighted skew, or `Generalized`  to use only the
#' generalized, or regional skew. 
#' - `GenSkew` - regional (generalized) skew coefficient (required if the skew
#' option is set to Weighted or Generalized)
#' - `SkewSE` - standard error for regional skew coefficient(required if the skew
#' option is set to Weighted or Generalized)
#' - `LOType` - low outlier test set to `MGBT` to use to multiple Grubbs-Beck test
#' or `FIXED` to set the low outlier threshold to a specific value
#' - `LoThresh` - low outlier threshold (required if the low outlier test is set
#' to `FIXED`)
#' - `Urb/Reg` - set to `Yes` to include urbanized and regulated peaks in analysis
#' - `WeightOpt` - skew weighting option. Selects one of three algorithms to use
#' for weighting the at-site and regional skew values. Valid options are: 
#' 1. INV - inverse mean-squared error weighting. This is the algorithm used in 
#' PeakFQ 7.4.1 and earlier
#' 2. ERL - effective record length weighting. This algorithm is used in HEC-SSP 2.3
#' 3. HWN - a generalization of the PeakFQ 7.4.1 algorithm using an optimized
#' adjustment factor when censored data are present. Results are identical to
#' PeakFQ 7.4.1 when no censored data are present.
#' 
#' The `PCPT_Thresh`, `Interval`, and `Peak` fields may be used as many times as
#' needed for each site. The other station specification fields may 
#' only be used once for each site. 
#' 
#' @section Example Specifications:
#' 
#' ```
#' I RDB 06666667_nwis_peak.txt
#' I INFO 06666667_nwis_site.txt
#' O File 06666667_log.txt
#' O Plot Format PNG
#' O Plot Position 0
#' O ConfInterval 0.9
#' O Extended YES
#' Station 06666667
#'     PCPT_Thresh 1974 2021 0 1E+20 Default
#'     PCPT_Thresh 2005 2005 0.5 1E+20 PEAK < STATED VALUE
#'     Interval 2005 0 0.5 PEAK < STATED VALUE
#'     SkewOpt Weighted
#'     GenSkew -0.3681
#'     SkewSE 0.64
#'     LOType FIXED
#'     LoThresh 10
#'     WeightOpt INV
#' 
#' ```
#' 
#' 
#' @examples
#' \dontrun{
#' #Location of analysis specifications file
#' psf_path <- system.file("extdata", "ExampleSpecifications.psf", package="peakfq")
#' 
#' specs <- readPSF(psf_path)
#' 
#' specs_analysis <- specs[[1]]
#' specs_perception_thresholds <- specs[[2]]
#' specs_flow_intervals <- specs[[3]]
#' specs_peaks <- specs[[4]]
#' }
#' 
#'
#' @export
readPSF <- function(PSFfile){


  #Check if input file exists
  if (!file.exists(PSFfile)) {
    stop(paste("PSF file doesn't exist:", PSFfile))
  }


  PSFfile <- file(PSFfile, "r") #Open connection to read file
  on.exit(close(PSFfile)) #Close .EXP connection on exit


  #Setup dataframe for site specifications
  PSFcolumns <- c("site_no", "SkewOpt", "GenSkew", "SkewSE", "LOType", "LoThresh", "Urb/Reg", "WeightOpt")
  PSFcolumns_old <- c("Analyze", "BegYear", "EndYear", "HistPeriod")
  PSFdata <- data.frame(matrix(data=NA, 0, length(PSFcolumns))) #Create empty dataframe to store .PSF data in
  colnames(PSFdata) <- PSFcolumns #Name the columns with fields in .PSF file
  
  #Create empty dataframe to store thresholds, intervals, and peaks
  PSFthresholds <- data.frame(matrix(nrow=0, ncol=7))
  colnames(PSFthresholds) <- c("site_no", "label", "start", "end","min", "max", "comment")
  
  PSFintervals <- data.frame(matrix(nrow=0, ncol=6))
  colnames(PSFintervals) <- c("site_no", "label", "peak_WY", "interval_low", "interval_up", "comment")
  
  PSFpeaks <- data.frame(matrix(nrow=0, ncol=5))
  colnames(PSFpeaks) <- c("site_no", "label", "peak_WY", "peak_va", "comment")


  #Set global specification values to NA initially
  inFormat <- NA
  inFile <- NA
  inInfo <- NA #Path to site information file
  outFile <- "" #Ok if this is empty string, so default to this
  ConfInterval <- NA
  PlotPosition <- NA
  PlotFormat <- NA
  extended <- FALSE

  #Read the first site
  reading <- TRUE #Create a variable for whether the file is being read

  while (reading == TRUE){

    nextLine <- readLines(PSFfile, n=1, warn = FALSE) #Read next line of file

    if(length(nextLine) > 0){

      if(is.na(strsplit(nextLine, " ")[[1]][1])){
        print("Blank line in .psf file")
      }

      else if(strsplit(nextLine, " ")[[1]][1] == "Station"){

        if(exists("newRow")){ #This should be modified to not use exists()
          PSFdata <- rbind(PSFdata, newRow) #If a row of data has been assembled, add it to the dataframe
        }

        siteNum <- strsplit(nextLine, " ")[[1]][2] #Get station number of the new station
        newRow <- data.frame(matrix(data=NA, 1, length(PSFcolumns)))
        colnames(newRow) <- PSFcolumns
        newRow$site_no <- siteNum

      }


      else if(substring(strsplit(trimws(nextLine), " ")[[1]][1],1,1) != "'"){
        #If the next line isn't a comment, read it

        fieldName <- strsplit(trimws(nextLine), " ")[[1]][1]

        if(fieldName == "PCPT_Thresh" | fieldName== "Interval" | fieldName == "Peak"){
          #Handle fields which can occur multiple times for a site

          fieldValue <- trimws(gsub("\\s+", " ", nextLine)) #Condense multiple spaces to one
          
          
          switch(fieldName, 
                 PCPT_Thresh = PSFthresholds <- rbind(PSFthresholds,
                                          data.frame(site_no = siteNum, 
                                          label = strsplit(fieldValue, " ")[[1]][1], 
                                          start = as.integer(strsplit(fieldValue, " ")[[1]][2]), 
                                          end = as.integer(strsplit(fieldValue, " ")[[1]][3]), 
                                          min = as.numeric(strsplit(fieldValue, " ")[[1]][4]), 
                                          max = as.numeric(strsplit(fieldValue, " ")[[1]][5]), 
                                          comment = trimws(gsub('^(?:[^ ]* ){5}', "", 
                                                                paste0(fieldValue, " "))) #use a regular expression to get the last part of the line, paste0() needed for case of no comment
                                          )),  
                 Interval = PSFintervals <- rbind(PSFintervals,
                                       data.frame(site_no = siteNum, 
                                       label = strsplit(fieldValue, " ")[[1]][1], 
                                       peak_WY = as.integer(strsplit(fieldValue, " ")[[1]][2]), 
                                       interval_low = as.numeric(strsplit(fieldValue, " ")[[1]][3]), 
                                       interval_up = as.numeric(strsplit(fieldValue, " ")[[1]][4]), 
                                       comment = trimws(gsub('^(?:[^ ]* ){4}', "",
                                                             paste0(fieldValue, " "))) #use a regular expression to get the last part of the line, paste0() needed for case of no comment
                                       )), 
                 Peak = PSFpeaks <- rbind(PSFpeaks, 
                                   data.frame(site_no = siteNum, 
                                   label = strsplit(fieldValue, " ")[[1]][1], 
                                   peak_WY = as.integer(strsplit(fieldValue, " ")[[1]][2]), 
                                   peak_va = as.numeric(strsplit(fieldValue, " ")[[1]][3]), 
                                   comment = trimws(gsub('^(?:[^ ]* ){3}', "",
                                                         paste0(fieldValue, " "))) #use a regular expression to get the last part of the line, paste0() needed for case of no comment
                 ))
                 )
          

        }
        
        else if(fieldName == "I"){
          
          #Read input filepath
          #This is slightly complicated due to change to 
          #taking tab delimited data and need to read
          #separate site information file
          
          inKey <- strsplit(trimws(nextLine), " ")[[1]][2]
          
          if(inKey == "ASCI"){
            inFormat <- strsplit(trimws(nextLine), " ")[[1]][2]
            inFile <- strsplit(trimws(nextLine), " ")[[1]][3]
            warning("Support for WATSTORE input data is deprecated")
          }
          else if(inKey == "RDB"){
            inFormat <- strsplit(trimws(nextLine), " ")[[1]][2]
            inFile <- strsplit(trimws(nextLine), " ")[[1]][3]
          }
          else if(inKey == "INFO"){
            inInfo <- strsplit(trimws(nextLine), " ")[[1]][3]
          }
          else{
            stop(paste("Unknown input file format", inFormat))
          }
          
          
          
        }
        
        else if(fieldName == "O"){
          
          #Read in output options
          
          fieldName <- strsplit(trimws(nextLine), " ")[[1]][2]
          fieldValue <- strsplit(trimws(nextLine), " ")[[1]][3]
          
          #Read in confidence interval if supplied
          if(fieldName == "ConfInterval"){
            ConfInterval <- as.numeric(fieldValue)
          }
          
          else if(fieldName == "File"){
            outFile <- fieldValue
          }
          
          #Handle plot controls
          else if(fieldName == "Plot"){
            fieldName <- strsplit(trimws(nextLine), " ")[[1]][3]
            fieldValue <- strsplit(trimws(nextLine), " ")[[1]][4]
            
            if(fieldName == "Position"){
              PlotPosition <- as.numeric(fieldValue)
            }
            else if(fieldName == "Format"){
              PlotFormat <- fieldValue
            }
            else{
              warning(paste("Invalid input keyword", fieldName))
            }
            
            
            
          }
          
          else if(fieldName == "EXTENDED"){
            if(toupper(fieldValue) == "YES"){
              extended <- TRUE
            }
            else if(toupper(fieldValue) != "NO"){
              extended <- FALSE
              warning(paste("Invalid input keyword in specifications:", fieldName, fieldValue))
            }
            #Don't need else case since variable is initialized FALSE
          }
          
          
          
        }



        else{

          fieldValue <- strsplit(trimws(nextLine), " ")[[1]][2]
          
          if(fieldName %in% PSFcolumns){
            newRow[1,fieldName] <- fieldValue
          }
          else if(fieldName %in% PSFcolumns_old){
            warning(paste("Deprecated PSF field name", fieldName, "will not be read"))
          }
          else{
            warning(paste("Unknown PSF field name:", fieldName))
          }
          
          
        }

      }
    }

    else{
      reading <- FALSE
      if(exists("newRow")){
        PSFdata <- rbind(PSFdata, newRow) #If a row of data has been assembled, add it to the dataframe
      }
    }

  }
  
  
  attributes(PSFdata)$peakFormat <-  inFormat
  attributes(PSFdata)$peakFile <- inFile
  attributes(PSFdata)$siteinfoFile <- inInfo
  attributes(PSFdata)$outFile <- outFile
  attributes(PSFdata)$ConfInterval <- ConfInterval
  attributes(PSFdata)$PlotFormat <- PlotFormat
  attributes(PSFdata)$PlotPosition <- PlotPosition
  attributes(PSFdata)$extended <- extended

  return(list(PSFdata, PSFthresholds, PSFintervals, PSFpeaks))

}


#' Build QT Dataframe
#'
#' Function to build dataframe of flow intervals and perception thresholds
#'
#' @param site_no character - unique identifier for site for which to build dataframe
#' @param inData dataframe containing peak flow data with the following columns: 
#'   1. peak_WY integer - water year of peak flow value
#'   2. peak_va numeric - annual peak flow value in cubic feet per second
#'   3. peak_cd character - NWIS peak flow code (\url{https://nwis.waterdata.usgs.gov/nwis/peak?help#flow_qual_cd})
#'   
#' @param psfIntervals dataframe containing flow interval data with the following columns:
#'   1. site_no character - unique identifier for site
#'   2. peak_WY integer - water year for the flow interval
#'   3. interval_low numeric - lower flow interval value
#'   4. interval_up numeric - upper flow interval value
#' @param psfThresholds  dataframe containing perception threshold data with the following columns:
#'   1. site_no character - unique identifier for site
#'   2. start integer - start water year for perception threshold
#'   3. end integer - end water year for perception threshold
#'   3. min numeric - minimum perceptible value
#'   4. max numeric - maximum perceptible value
#' @param psfPeaks dataframe containing user-specified peak flow data with the following columns:
#'   1. site_no character - unique identifier for site
#'   2. peak_WY integer - water year for peak flow value
#'   3. peak_va numeric - peak flow value
#'   
#' @param removeUrbReg logical - option indicating whether or not to remove data with NWIS peak flow code 6 or code C.
#' 
#' @param psfSpecs dataframe - optional dataframe with .psf specifications, will 
#' override removeUrbReg option if supplied
#' 
#' @param keepNoInfo logical - option indicating whether rows with no information
#' (lower perception threshold of infinity) should be kept
#'
#' @return Dataframe with the following columns:
#'   1. ql numeric - lower flow interval values
#'   2. qu numeric - upper flow interval values
#'   3. tl numeric - lower perception threshold values
#'   4. tu numeric - upper perception threshold values
#'   5. dtype integer - indicator if value is historic peak. 1 = historic peak, 0 = not historic peak. 
#'   
#'   The returned dataframe has the site number stored in the attribute `site_no`. 
#'   
#' @section Processing Rules: 
#' 
#' 
#' The following procedure is used when processing the input data. 
#' If the specifications include multiple thresholds, peaks, or intervals, for
#' the same year, then the last one specified is given priority. 
#' 
#' \enumerate{
#'   \item Start and end water years are taken as first and last years for which a
#' perception threshold is defined.
#'\itemize{
#'   \item All input data with water years outside of this range are thrown out. No warning is given
#'   \item If perception thresholds are not specified for one or more years between the
#' start and end years of an analysis, an error is thrown. 
#' }
#'   \item Flow intervals and perception thresholds for all years in the input data
#' are populated with initial values. The lower and upper flow interval are set 
#' equal to the data value, the lower perception threshold is set to zero and the
#' upper perception threshold is set to infinity.
#' 
#' \item Qualification codes are applied to the flow intervals
#'\itemize{
#'   \item Code 4 peaks have the lower flow interval set to zero
#'   \item Code 8 peaks have the upper flow interval set to infinity
#'   \item Code 3 and O peaks are removed from the analysis by setting the lower flow
#' interval to zero and the upper interval to infinity
#'  \item Code 6 and C peaks are removed from the analysis, unless the Urban/Regulated
#' option is set to TRUE, by setting the lower flow
#' interval to zero and the upper interval to infinity
#' }
#' 
#' \item `Peak` keyword from specifications is processed. Upper and lower flow
#' intervals for the specified years are set to the associated peak value in the specifications.
#' 
#' \item `Interval` keyword from specifications is processed. Upper and lower flow 
#' intervals for the specified years are set to the associated interval values
#' in the specifications. This will overwrite the `Peak` keyword if applied to
#' the same year. A warning is given if this happens.
#' 
#' \item Perception thresholds from the specifications are applied in the order
#' they appear in the specifications. If a perception threshold is not applied
#' for all years in the analysis period, a an error is thrown.
#' 
#' \item If a year has a perception threshold set, but no flow interval has been
#' processed from the input data or specifications, the lower flow interval is set to
#' zero and the upper flow interval is set to the lower perception threshold
#' associated with the interval.
#' 
#' \item The lower perception threshold is set to infinity for years with a flow
#' interval of zero to infinity
#' 
#' \item Years with a lower perception threshold of infinity (indicating no
#' information) are removed from the analysis, unless the `keepNoInfo` input is
#' set to `TRUE`.
#' 
#' }
#'
#'
#'
#' @keywords internal
#'
siteQT <- function(site_no, inData, psfIntervals, psfThresholds, psfPeaks, removeUrbReg = TRUE, psfSpecs = NA, keepNoInfo = FALSE){

  #Get data specific to site
  inData <- inData[inData$site_no == site_no,]
  psfIntervals <- psfIntervals[psfIntervals$site_no == site_no,]
  psfThresholds <- psfThresholds[psfThresholds$site_no == site_no,]
  psfPeaks <- psfPeaks[psfPeaks$site_no == site_no,]
  
  #Get urban/regulated option from psfSpecs if supplied
  if(is.data.frame(psfSpecs)){
    
    if(!is.na(psfSpecs$`Urb/Reg`[psfSpecs$site_no == site_no]) &
       toupper(psfSpecs$`Urb/Reg`[psfSpecs$site_no == site_no]) == "YES"){
      removeUrbReg = FALSE
    }
    
  }

  #Get year range of analysis
  startWY <- min(psfThresholds$start)
  endWY <- max(psfThresholds$end)
  
  if(!all(c("peak_WY", "peak_va", "peak_cd") %in% colnames(inData))){
    stop("Input inData is missing one of the required input columns")
  }
  
  QT <- inData[, c("peak_WY", "peak_va", "peak_cd")]

  QT <- QT[complete.cases(QT$peak_va),] #Remove GH only peaks
  
  if(any(QT$peak_va < 0)){
    stop(paste("Negative discharge in specifications for site", site_no))
  }
  
  
  if(nrow(QT) < 1) stop(paste("No peak flow data for site", site_no))

  if(!all(complete.cases(QT$peak_WY))) stop(paste("Some peak doesn't have associated year for site", site_no))

  #Remove peaks outside of PT range
  QT <- QT[QT$peak_WY >= startWY & QT$peak_WY <= endWY,]

  #Set default values for intervals and thresholds
  QT$ql <- QT$peak_va
  QT$qu <- QT$peak_va
  QT$tl <- NA
  QT$tu <- NA

  #Handle less than/greater than peaks
  QT$ql[grepl("4", QT$peak_cd)] <- Qmin #Less than peaks
  QT$qu[grepl("8", QT$peak_cd)] <- Qmax #Greater than peaks

  #Set dam breaks and opportunistic peaks as intervals from zero to infinity by default (can be overridden by later specifications)
  QT$ql[grepl("3", QT$peak_cd)] <- Qmin
  QT$qu[grepl("3", QT$peak_cd)] <- Qmax
  
  QT$ql[grepl("O", QT$peak_cd)] <- Qmin
  QT$qu[grepl("O", QT$peak_cd)] <- Qmax
  

  #Remove urban/regulated peaks if specified
  if(removeUrbReg){
    QT$ql[grepl("6", QT$peak_cd) | grepl("C", QT$peak_cd)] <- Qmin
    QT$qu[grepl("6", QT$peak_cd) | grepl("C", QT$peak_cd)] <- Qmax
    
    QT$tl[grepl("6", QT$peak_cd) | grepl("C", QT$peak_cd)] <- Qmax
    QT$tu[grepl("6", QT$peak_cd) | grepl("C", QT$peak_cd)] <- Qmax
  }

  QT <- merge(QT, data.frame(peak_WY = c(startWY:endWY)), all = T)
  
  #Check for duplicated years in input data
  if(any(duplicated(QT$peak_WY))){
    
    stop(paste("Duplicated water year", QT$peak_WY[duplicated(QT$peak_WY)],
               "in input for site", site_no))
    
  }

  #Start with specified discrete peaks in the .psf
  if(nrow(psfPeaks) > 0){

    warning(paste("Deprecated keyword 'Peak' used in .psf file for site", site_no))

    for(r in 1:nrow(psfPeaks)){

      WY <- psfPeaks$peak_WY[r]

      if(WY %!in% QT$peak_WY){
        stop(paste("Water year", WY, "is not in analysis range for site", site_no))
      }
      
      if(psfPeaks[r, c("peak_va")] >= 0){
        QT[QT$peak_WY == WY, c("ql", "qu")] <- psfPeaks[r, c("peak_va")]
      }
      else{
        
        if(psfPeaks[r, c("peak_va")] == -8888){
          stop(paste("Negative discharge", psfPeaks[r, c("peak_va")], 
                     "in specifications for water year", WY, "for site", site_no, "\n", 
                     "This error is usually caused by loading legacy specification files where 'Peak XXXX -8888' was used to exclude data from water year XXXX from the analysis. This can be fixed by replacing all occurrences of 'Peak XXXX -8888' in the specifications file with 'Interval XXXX 0 1E+20'"))
        }
        else{
          stop(paste("Negative discharge", psfPeaks[r, c("peak_va")], 
                     "in specifications for water year", WY, "for site", site_no))
        }
        
        
      }

      

    }
  }

  #Next do intervals
  if(nrow(psfIntervals) > 0){
    
    #Give warning if peak is overwritten by interval
    if(nrow(psfPeaks) > 0){
      if (any(psfPeaks$peak_WY %in% psfIntervals$peak_WY)){
        warning(paste("The following peaks in the specifications will be
                overwritten by intervals in the specification: ", 
                      paste(psfPeaks$peak_WY[psfPeaks$peak_WY %in% psfIntervals$peak_WY], collapse = ", " )))
      }
      
    }
    
    for(r in 1:nrow(psfIntervals)){

      WY <- psfIntervals$peak_WY[r]

      if(WY %!in% QT$peak_WY){
        stop(paste("Water year", WY, "is not in analysis range for site", site_no))
      }

      QT[QT$peak_WY == WY, c("ql", "qu")] <- psfIntervals[r, c("interval_low", "interval_up")]

    }
  }




  #Then do thresholds
  for(r in 1:nrow(psfThresholds)){

    Tstart <- psfThresholds$start[r]
    Tend <- psfThresholds$end[r]

    if(Tstart %!in% QT$peak_WY){
      stop(paste("Water year", Tstart, "is not in analysis range for site", site_no))
    }

    if(Tend %!in% QT$peak_WY){
      stop(paste("Water year", Tend, "is not in analysis range for site", site_no))
    }

    QT[QT$peak_WY >= Tstart & QT$peak_WY <= Tend, c("tl", "tu")] <- psfThresholds[r, c("min", "max")]

  }
  
  if(length(QT$peak_WY[!complete.cases(QT[, c("tl", "tu")])]) > 0 ){
    stop(paste("The following years in analysis period are missing perception thresholds:", 
               paste(QT$peak_WY[!complete.cases(QT[, c("tl", "tu")])], collapse = ",") , "for", site_no))
  }


  #Fill in years with intervals implied by PTs
  if(!all(QT$tl[is.na(QT$ql)] > Qmin)){
    warning(paste0("ZERO PERCEPTION THRESHOLD FOR MISSING YEARS FOR SITE ",
                  site_no, ". THESE YEARS WILL BE TREATED AS MISSING DATA."))
    
    #Set missing years with perception threshold of zero to missing data
    QT$ql[is.na(QT$ql) & QT$tl <= Qmin] <- Qmin
    QT$qu[is.na(QT$qu) & QT$tl <= Qmin] <- Qmax
    
  }
  QT$ql[is.na(QT$ql)] <- Qmin
  QT$qu[is.na(QT$qu)] <- QT$tl[is.na(QT$qu)]

  #Make sure all intervals and PTs are greater than Qmin
  QT$ql[QT$ql < Qmin] <- Qmin
  QT$qu[QT$qu < Qmin] <- Qmin
  QT$tl[QT$tl < Qmin] <- Qmin
  QT$tu[QT$tu < Qmin] <- Qmin #Not sure why this case would come up, but..
  
  
  QT$tl[QT$ql == Qmin & QT$qu == Qmax] <- Qmax #Set lower PT to Qmax for years with no information
  
  #Add column with check for historic peaks
  QT$dtype <- 0
  QT$dtype[grepl("7", QT$peak_cd)] <- 1
  
  attributes(QT)$site_no <- site_no
  
  if(!keepNoInfo){
    QT <- QT[QT$tl < Qmax, ] #Remove years with no information
  }
  
  if(nrow(QT) == 0){
    warning(paste0("No data available after applying analysis specifications for site ", site_no,
                  ". Make sure option to include urban/regulated peaks is set appropriately."))
  }
  
  else if(nrow(QT) < 10){
    warning(paste("Less than 10 years of data for site", site_no, 
            "Bulletin 17C recommends at least 10 years of data for analysis."))
  }

  return(QT)
}

