# C      + + + PURPOSE + + +
# C     This file contains the main function for
# C     processing for frequency analysis and output. Replaces FORTRAN subroutine
# C     J407 -- USGS-WRC FLOOD FREQUENCY ANALYSIS PER WRC BULL 17-B, 1981.
# C
# C
# C     + + + HISTORY + + +
# C     VER 76.00 BY WKIRBY, WRD-NR, MAY 1976. (BULL.17)
# C     VER 2.0 BY WKIRBY, WRD-NR, APRIL 77.  (BULL.17-A)
# C     VER 3.0 BY WKIRBY, WRD-NR, MAY 1979.
# C     VER 3.7P - PRIME REVISIONS - K.FLYNN 12/83.
# C     VER 3.8P - WK 12/86,  7/88.
# C     SET ARGUMENTS = 0 FOR NON-ANNIE/NON-WDM USE
# C     VER 3.9P - WK, AML 12/88
# C     VER 3.9A-P - WK, AML 2/89
# C     MODIFIED 8/9/89 AML (deleted BLOCKDATA)
# C     Modified 6/93 AML to coding convention and add requirements
# C                       for distribution by Texas, changed to an
# C                       80 char/record print file, made Z,H,N,Y
# C                       input records optional 
# C     Updated for batch version of PEAKFQ, 9/03
# C     Paul Hummel of AQUA TERRA Consultants
# C     VER 8.0 BY SETH SIEFKEN, WY-MT WSC 2024
# C             Converted to R package, removed all FORTRAN not used for EMA
# C             math routines. 
# C     



#Define global variables
Qmin <- 10^-20 #Minimum flow
Qmax <- 10^20 #Maximum flow


#Define a "not in" operator
'%!in%' <- function(x,y)!('%in%'(x,y)) 



#' Run peakfq Analysis
#'
#' Run peak-flow frequency analysis for all sites in a specifications file. 
#' 
#'
#' @param inPSF character - path to specifications file for peakfq
#' @param quietly logical - if `TRUE`, messages from the function are not printed
#'
#' @returns A list containing:
#' \itemize{
#'   \item lp3 -  Estimated moments and other information
#'   \item qnt - Computed frequency curve
#'   \item emp - Empirical plotting positions and EMA data representation
#'   \item mgb - MGBT  p values
#'   \item trd - Trend test results
#'   \item ffa_plots - Flood frequency curve plots
#'   \item inp_plots - Input data plots
#'   \item ema_plots - EMA data representation plots
#'   \item log - Log text with output summary and warnings
#' }
#'  
#' @details
#' 
#' Runs peak-flow frequency analysis to following the guidelines in 
#' [Bulletin 17C](https://doi.org/10.3133/tm4B5) (England and others, 2019).
#' 
#' See \code{vignette("peakfq_app", package = "peakfq")} for a description of 
#' the implementation of Bulletin 17C guidelines in the peakfq package. 
#' 
#' Numeric data from the analysis are returned in a list of dataframes and also
#' saved as .csv files if the `O File` field in the specifications file is set 
#' to a valid file. For example: `O File ./PATH/BASE_FILE_NAME.txt`
#' 
#' If output files are generated, they will be named as: 
#' \itemize{
#' \item{basename_lp3.csv - CSV with estimated moments, PILF information, and other information used to fit the log-Pearson Type 3 distribution}
#' \item{basename_qnt.csv - CSV with quantile estimates, variances, and confidence intervals }
#' \item{basename_emp.csv - CSV with input data representation and empirical plotting positions }
#' \item{basename_mgb.csv - CSV with MGBT p-values }
#' \item{basename_trd.csv - CSV with Mann-Kendall trend test results}
#' \item{basename.txt - Log file with errors or warnings from the computation code.}
#'}
#'
#'
#' *Please note that in the above text, "basename" refers to the `BASE_FILE_NAME` supplied in the `O File` field in the specifications file.* 
#'
#' Output data plots for each site will be saved if the `O File` field in the specifications file is set and
#' the `O Plot Format` field is set to one of the supported plot formats, such as `O Plot Format SVG`.
#' Output plots will be named as: 
#' 
#' \itemize{
#' \item{siteno_data.* - Plot of data input data }
#' \item{siteno_ema.* - Plot of data used in the expected moments algorithm after low outlier handling has been applied }
#' \item{siteno.* - Plot of flood frequency curve }
#'}
#'  
#' 
#' @section Numeric Outputs: 
#' 
#' The numeric data returned as dataframes and optionally saved as .csv files
#' are formatted as follows: 
#' 
#' \enumerate{
#'   \item lp3 -  Estimated moments and other information
#'   \itemize{
#'   \item{site_no - unique site identifier}
#'   \item{station_nm - site (station) name}
#'   \item{BegYear - Beginning year of the analysis}
#'   \item{EndYear - Ending year of the analysis}
#'   \item{RecordLength - Number of years of data in the analysis. Does not include missing data years for which no information is available (lower perception threshold of infinity)}
#'   \item{SkewOption - Analysis skew option - either "Station", "Weighted", or "Generalized"}
#'   \item{Mean - Mean (in log units) of data}
#'   \item{StandDev - Standard deviation (in log units) of data}
#'   \item{Skew - Skew of log-transformed data}
#'   \item{AtSiteMean - Mean (in log units) of data, without regional skew}
#'   \item{AtSiteStandDev - Standard deviation (in log units) of data, without regional skew}
#'   \item{AtSiteSkew - Skew of log-transformed data, without regional skew}
#'   \item{AtSiteMSEG - Mean Squared Error of the at-site skew coefficient}
#'   \item{AtSiteMSEG_GagedOnly - Mean Squared Error of the at-site skew coefficient without censoring}
#'   \item{RegSkew - Regional skew coefficient}
#'   \item{RegMSEG - Mean Squared Error of the regional skew coefficient}
#'   \item{HistPeaks - Number of historical peaks in the analysis}
#'   \item{MGBTPeaks - Number of peaks available for Multiple Grubbs-Beck low outlier test. Computed as number of exactly known (lower flow interval value equal to upper flow interval value) input peaks, excluding historic peaks, plus any not exactly known peaks with the upper flow interval value less than the smallest exactly known peak value.}
#'   \item{PILF_Method - Potentially Influential Low Flow threshold method; either "MGBT" or "FIXED"}
#'   \item{PILF_Thresh - Potentially Influential Low Flow threshold}
#'   \item{PILFs - Number of values below the PILF threshold}
#'   \item{PILF_0s - Number of zero values below the PILF threshold. Zero values are always counted as PILFs.}
#'   \item{weightOpt - Option for algorithm to computed weighted skew coefficient. Possible values are HWN (default; a generalization of the PeakFQ 7.4.1 algorithm using an optimized adjustment factor when censored data are present, results are identical to PeakFQ 7.4.1 when no censored data are present), INV (identical to PeakFQ 7.4.1 and earlier), and ERL (identical to HEC-SSP 2.3)}
#'   \item{WeightCo - Censored data adjustment factor for computing weighted skew coefficient}
#'   \item{RunTime - Time to run analysis, in seconds}
#' }
#'   \item qnt - Computed frequency curve
#' \itemize{
#'   \item{site_no - unique site identifier}
#'   \item{EXC_Prob - Annual exceedance probability}
#'   \item{Estimate - Estimated discharge, in cubic feet per second}
#'   \item{Variance - Variance of estimate, in log units}
#'   \item{Conf_Low - Lower confidence limit, in cubic feet per second}
#'   \item{Conf_Up - Upper confidence limit, in cubic feet per second}
#'   \item{K_Value - K value from Bulletin 17B Appendix 3}
#'   \item{Conf_Interval - confidence interval. Default is 0.90 (90-percent confidence interval)}
#' }   
#'   \item emp - Empirical plotting positions and EMA data representation
#' \itemize{
#'   \item{site_no - unique site identifier}
#'   \item{peak_WY - water year of peak flow}
#'   \item{peak_va - peak flow value}
#'   \item{plot_pos - plotting position of peak flow, computed using Hirsch-Stedinger plotting positions}
#'   \item{ql - lower flow interval value, in cubic feet per second}
#'   \item{qu - upper flow interval value, in cubic feet per second}
#'   \item{tl - lower perception threshold, in cubic feet per second}
#'   \item{tu - upper perception threshold, in cubic feet per second}
#'   \item{historic - true or false indicator of whether or not data value is a historic peak.}
#'   \item{ql_ema - lower flow interval value used in Expected Moments Algorithm, in cubic feet per second. Will be identical to `ql` if no PILF values are identified.}
#'   \item{qu_ema - upper flow interval value used in Expected Moments Algorithm, in cubic feet per second. Will be identical to `qu` if no PILF values are identified.}
#'   \item{tl_ema - lower perception threshold used in Expected Moments Algorithm, in cubic feet per second.}
#'   \item{tu_ema - upper perception threshold used in Expected Moments Algorithm, in cubic feet per second.}
#' }   
#'   \item mgb - MGBT  p values
#' \itemize{
#'   \item{site_no - unique site identifier}
#'   \item{peak_va - peak flow value identified as a low outlier by the Multiple Grubbs-Beck Test. NA values indicate no low outliers were identified for the site.}
#'   \item{p_value - p-value associated with peak flow value. NA values indicate no low outliers were identified for the site.}
#' }   
#'   \item trd - Trend test results
#' \itemize{
#'   \item{site_no - unique site identifier}
#'   \item{tau - Kendall's tau value}
#'   \item{p_value - p-value from Mann-Kendall test}
#'   \item{median_slope - Median slope from Mann-Kendall test}
#'   \item{n_peaks - number of peak values used in Mann-Kendall test}
#' }
#' }
#' 
#'  
#'  
#' @examples
#' \dontrun{
#' library(peakfq)
#' 
#' #Location of analysis specifications file
#' psf_path <- system.file("extdata", "ExampleSpecifications.psf", package="peakfq")
#' 
#' #Run analysis from the specifications file
#' results <- peakfq(psf_path)
#'  
#' #View Outputs
#'  
#' #Information on moments computed to fit the Log-Pearson Type III distribution
#' results[["lp3"]] 
#'  
#' #Annual exceedance probability estimates with variance and confidence intervals
#' results[["qnt"]]
#'  
#' #Empirical plotting positions and EMA data representation
#' results[["emp"]]
#'  
#' #MGBT p values
#' results[["mgb"]]
#'  
#' #Trend test results
#' results[["trd"]]
#'  
#' #Plots
#' #Note that plotting window must be correctly sized for plots to display properly;
#' #plots are optimized for a 9 inch wide by 6 inch high window.
#' results[["ffa_plots"]][[1]] #Frequency curve from first site
#' results[["inp_plots]][[2]] #Input data plot for second site
#' results[["ema_plots]][[3]] #EMA data representation plot for third site
#' }
#' 
#' @references England, J.F., Jr., Cohn, T.A., Faber, B.A., Stedinger, J.R., Thomas, W.O., Jr., Veilleux, A.G., Kiang, J.E., and Mason, R.R., Jr., 2019, Guidelines for determining flood flow frequency — Bulletin 17C (ver. 1.1, May 2019): U.S. Geological Survey Techniques and Methods, book 4, chap. B5, 148 p., \doi{doi:10.3133/tm4B5}.
#'  
#' @importFrom svglite svglite
#' @importFrom stringr str_replace  
#'
#' @export
peakfq <- function(inPSF, quietly=TRUE){
  self = environment()
  
  wText <- list() #List of warnings
  
  #Read input data and save warnings
  PSFlist <- withCallingHandlers(readPSF(inPSF), 
                     warning = function(w){
                                           self$wText <- c(self$wText, conditionMessage(w))
                                           }
                     )
  
  specs <- PSFlist[[1]]
  psfThresholds <- PSFlist[[2]]
  psfIntervals <- PSFlist[[3]]
  psfPeaks <- PSFlist[[4]]
  
  inDataPath <- file.path(dirname(inPSF), attributes(specs)$peakFile)

  if(!file.exists(inDataPath)){   #Check if input data file exists
    inDataPath <- attributes(specs)$peakFile # set inDataPath to reflect temp file location if running on server
    if(!file.exists(inDataPath)){ # Check if input data file exist
      stop(paste("Input data file",
               inDataPath,
               "not found"))
    }
  }
  
  #Read input data from file in .psf 
  if(attributes(specs)$peakFormat == "ASCI"){
    
    inData <- readWATSTOREAll(inDataPath)
    siteInfo <- attributes(inData)$siteInfo
  }
  else if(attributes(specs)$peakFormat == "RDB"){
    
    inData <- readRDB(inDataPath)
    siteInfo <- RDB_station_nm(inDataPath)
    
  }
  else{
    stop(paste("Invalue input data keyword",
               attributes(specs)$peakFormat,
               "in specifications"))
  }
  
  outDir <- file.path(dirname(inPSF), dirname(attributes(specs)$outFile))
  outFile <- strsplit(basename(attributes(specs)$outFile), "\\.")[[1]][1] 
  
  
  ConfInterval <- 0.90 #Default to 0.90 if not supplied in specifications
  
  if(!is.null(attr(specs, "ConfInterval"))){
    if(!is.na(attr(specs, "ConfInterval"))){
      
      if(attr(specs, "ConfInterval") <= 0 | attr(specs, "ConfInterval") >= 1){
        warning(paste("Invalid confidence interval option: ", attr(specs, "ConfInterval"), 
                      "Value must be greater than 0 and less than 1. Defaulting to", ConfInterval))
      }
      else{
        ConfInterval <- attributes(specs)$ConfInterval #Set to value is valid and supplied
      }
    }
  }

  
  
  
  if(is.numeric(attr(specs, "PlotPosition"))){
    
    PlotPosition_in <- as.numeric(attr(specs, "PlotPosition"))
    
    if(!is.na(PlotPosition_in) & PlotPosition_in >= 0 & PlotPosition_in <=0.5 ){
      PlotPosition <- PlotPosition_in
    }
    else{
      
      w <- paste("Input plotting position parameter", PlotPosition_in,
                       "is outside the valid range of [0, 0.5]. Default value of 0 will be used.")
      self$wText <- c(self$wText, w)
      warning(w)
      
      PlotPosition <- 0 #Default to for Weibull plotting if not supplied in specifications
    }
    
    
  }
  else{
    PlotPosition <- 0 #Default to for Weibull plotting if not supplied in specifications
  }
  
  
  
  
  sites <- specs$site_no
  
  #Process data for EMA and save warnings
  siteData <-  withCallingHandlers(lapply(sites, siteQT, inData, psfIntervals, psfThresholds, psfPeaks, TRUE, specs), 
                                 warning = function(w){
                                   self$wText <- c(self$wText, conditionMessage(w))
                                 }
  )
  
  
  
  names(siteData) <- sites
  
  
  #Remove any sites with less than 10 rows of data or 8 of exact data
  exactCount <- function(d){
    if(nrow(d[d$ql == d$qu,]) < 8){
    return(TRUE)
      }
    else{
      return(FALSE)
    }
  }
  
  sites_short <- sites[unlist(lapply(siteData, nrow)) < 10]
  sites_short <- c(sites_short, sites[unlist(lapply(siteData, exactCount))])
  sites_short <- unique(sites_short)
  
  if(length(sites_short) > 0){
    
    w <- paste("The following sites have less than 10 total peak values or less than 8 uncensored peak values used in the analysis and have been removed:", 
               paste(sites_short, collapse=", "))
    
    warning(w)
    
    wText <- c(wText, w) #Save warning for log file
    
    siteData <- siteData[names(siteData) %!in% sites_short]
    sites <- sites[sites %!in% sites_short]
    specs <- specs[specs$site_no %!in% sites_short,]
  }
  
  #Make sure we still have some valid sites. 
  if(length(sites) == 0){
    stop("No valid sites for analysis after removing sites with less than 10 total peak values or less than 8 uncensored peak values.")
  }
  
  
  
  validLOType <- c("MGBT", "FIXED", "NONE")
  
  if(!all(specs$LOType %in% validLOType)){
    
    stop(paste("Invalid low outlier option in specifications for sites:", 
               paste(specs$site_no[specs$LOType %!in% validLOType], collapse = ",")))
  }
  
  #Set up numerical inputs from specifications
  specs$LoThresh <- as.numeric(specs$LoThresh)
  specs$LoThresh[specs$LOType == "MGBT"] <- 0
  specs$LoThresh[specs$LOType == "NONE"] <- Qmin
  
  if(!all(!is.na(specs$LoThresh))){
    stop(paste("Invalid low outlier threshold in specifications for sites:", 
               paste(specs$site_no[is.na(specs$LoThresh)], collapse = ",")))
  }
  
  #Set up skew options
  specs$SkewOpt <- tolower(specs$SkewOpt)
  validSkewOpt <- c("station", "weighted", "regional")
  
  if(!all(specs$SkewOpt %in% validSkewOpt)){
    
    badSkewOpts <- specs$SkewOpt[specs$SkewOpt %!in% validSkewOpt]
    stop(paste("Invalid skew option(s) in specifications file:",badSkewOpts ))
  }
  
  specs$GenSkew <- as.numeric(specs$GenSkew)
  specs$SkewSE <- as.numeric(specs$SkewSE)
  
  specs$GenSkew[specs$SkewOpt == "station"] <- 0
  
  #Get specs for sites needing regional skew coefficients and SEs for checking
  rSkewSpecs <- specs[specs$SkewOpt == "weighted" | specs$SkewOpt == "regional",]
  
  if(!all(!is.na(rSkewSpecs$GenSkew))){
    badSites <- rSkewSpecs$site_no[is.na(rSkewSpecs$GenSkew)]
    stop(paste("Invalid regional skew coefficient at sites", 
               paste(badSites, collapse = ", ")))
  }
  
  if(!all(!is.na(rSkewSpecs$SkewSE))){
    badSites <- rSkewSpecs$site_no[is.na(rSkewSpecs$SkewSE)]
    stop(paste("Invalid regional skew standard error at sites", 
         paste(badSites, collapse = ", ")))
  }
  
  
  if(!all(rSkewSpecs$SkewSE > 0)){
    badSites <- rSkewSpecs$site_no[rSkewSpecs$SkewSE < 0]
    stop(paste("Input skew standard error is less than zero at sites", 
         paste(badSites, collapse = ", ")))
  }
  
  specs$SkewMSE <- -1e99 #Initialize all MSEs to -1e90
  
  if(length(specs$SkewMSE[specs$SkewOpt == "weighted"]) > 0){
    specs$SkewMSE[specs$SkewOpt == "weighted"] <- (specs$SkewSE[specs$SkewOpt == "weighted"])^2
  }
  
  if(length(specs$SkewMSE[specs$SkewOpt == "regional"]) > 0){
    #Regional skew option with non-zero skew MSE passes the regional skew mse as a negative values
    specs$SkewMSE[specs$SkewOpt == "regional"] <- -1*(specs$SkewSE[specs$SkewOpt == "regional"])^2
  }
  
  
  validWeightOpt <- c("INV", "ERL", "HWN")
  
  #Default to Halloween Method
  specs$WeightOpt[is.na(specs$WeightOpt)] <- "HWN"
  
  if(!all(specs$WeightOpt %in% validWeightOpt)){
    
    stop(paste("Invalid skew weighting option in specifications for sites:", 
               paste(specs$site_no[specs$WeightOpt %!in% validWeightOpt], collapse = ",")))
  }
  
  
  
  
  if(attributes(specs)$extended){
    AEPs <- rep(list(c(0.9950, 0.9900, 0.9800, 0.9750, 0.9600, 0.9500, 0.9000, 0.8000, 0.7000, 0.6667, 0.6000, 0.5704, 0.5000, 0.4292, 0.4000, 0.3000, 0.2000, 0.1000, 0.0500, 0.0400, 0.0250, 0.0200, 0.0100, 0.0050, 0.0020, 0.0010, 0.0001)), length(siteData))
  }
  else{
    AEPs <- rep(list(c(0.9950, 0.9900, 0.9800, 0.9750, 0.9600, 0.9500, 0.9000, 0.8000, 0.7000, 0.6667, 0.6000, 0.5704, 0.5000, 0.4292, 0.4000, 0.3000, 0.2000, 0.1000, 0.0500, 0.0400, 0.0250, 0.0200, 0.0100, 0.0050, 0.0020)), length(siteData))
      
  }
  
  siteResults <- mapply(emafit, siteData, specs$LoThresh, specs$GenSkew, specs$SkewMSE, ConfInterval, specs$WeightOpt, AEPs, sites, quietly, SIMPLIFY = FALSE)
  
  names(siteResults) <- sites
  
  #Compute empirical plotting positions
  empResults <- lapply(siteData, plotPosHS, PlotPosition)
  EMPoutputs <- EMPout(siteResults, empResults)
  
  
  emp <- do.call("rbind", EMPoutputs[[1]]) #Need do.call for multiple dataframes
  EMArepList <- EMPoutputs[[2]]
  #Trend test
  trendResults <- dplyr::bind_rows(lapply(siteData, trend_test))
  
  

  
  #Compile output dataframes
  EMAoutputs <- EMAout(siteResults)
  info <- EMAoutputs[[1]]
  qnt <- EMAoutputs[[2]]
  MGBTpvals  <- EMAoutputs[[3]]
  qntList <- EMAoutputs[[4]]
  
  #Fill in site names
  info$station_nm <- siteInfo[match(info$site_no, siteInfo$site_no), "station_nm"]
  info$station_nm[is.na(info$station_nm)] <- "" #Set missing names to empty string
  info <- info[, c(1, ncol(info), 2:(ncol(info)-1))] #rearrange columns
  
  #Add confidence level column
  qnt$Conf_Interval <- ConfInterval
  
  
  ffaPlots <- list()
  dataPlots <- list()
  emaPlots <- list()
  
  validPlotFormats <- c("png", "svg", "jpg", "jpeg", "none") #List of valid plot formats
  plotext <- tolower(attributes(specs)$PlotFormat)
  
  if(plotext %!in% validPlotFormats){
    warning(paste("Invalid plot format:", plotext))
    plotext <- "none"
  }
  else if(plotext != "none"){
    plotInfoText <- plotText(info) #Create text for ffa curve plots
	
	#Identify intentional missing data years
	nodatThresh <- psfThresholds[psfThresholds$min == Qmax, ] #Get thresholds indicating intentional missing data
	
	
	nodatYears <- mapply(seq, nodatThresh$start, nodatThresh$end, SIMPLIFY = FALSE)
	names(nodatYears) <- nodatThresh$site_no
	
	#Initialize to list with empty vector if no intentional no data years
	if(length(nodatYears) == 0) nodatYears <- list(integer(0))
	
	#Initialize empty list of missing years for each site
	nodatList <- replicate(length(info$site_no), integer(0), simplify = FALSE)
	names(nodatList) <- info$site_no

	#Combine no data years into a single vector for each site
	nodatYearsComb <- tapply(unlist(nodatYears, use.names = FALSE), rep(names(nodatYears), lengths(nodatYears)), FUN = c, simplify = FALSE)
	
	nodatList <- utils::modifyList(nodatList, nodatYearsComb)
    
    #Create file names for plots
    plotFiles <- paste0(outDir, "/", info$site_no, ".", plotext)
    dataplotFiles <- paste0(outDir, "/", info$site_no, "_inp.", plotext)
    emaplotFiles <- paste0(outDir, "/", info$site_no, "_ema.", plotext)
  
    #Make the plots
    ffaPlots <- mapply(ffaPlot, EMPoutputs[[1]], qntList, ConfInterval, info$site_no,
                       info$station_nm, info$PILF_Thresh, plotInfoText,
                       SIMPLIFY = FALSE)
    dataPlots <- mapply(pkPlot, empResults, info$site_no, info$station_nm, 
                        "identity", "User Input Peak-Discharge Data", 
                        nodatList,
                        SIMPLIFY = FALSE)
    emaPlots <- mapply(pkPlot, EMArepList, info$site_no, info$station_nm, 
                       "identity", "EMA Representation of Peak-Discharge Data", 
                       nodatList,
                       SIMPLIFY = FALSE)
  }
  
  
  
  #Output results
  outFiles <- c()
  
  if(!dir.exists(outDir)){
    warning(paste("Cannot save output. Directory does not exist:", outDir))
  }else if(is.na(outFile)){
    #No output file is specified, so no files will be saved. 
    message("Output files not saved")
  }
  else{
    
    #if(interactive()) {
      outFiles <- c(paste0(outDir, "/", outFile, "_lp3.csv"), 
                    paste0(outDir, "/", outFile, "_qnt.csv"), 
                    paste0(outDir, "/", outFile, "_emp.csv"), 
                    paste0(outDir, "/", outFile, "_mgb.csv"), 
                    paste0(outDir, "/", outFile, "_trd.csv"))
      
      
      withCallingHandlers(write.csv(info, outFiles[1], row.names = FALSE), 
                          error= function(e){
                            stop(paste("Cannot save output file", outFiles[1]))
                          }
      )
      
      withCallingHandlers(write.csv(qnt, outFiles[2], row.names = FALSE), 
                          error= function(e){
                            stop(paste("Cannot save output file", outFiles[1]))
                          }
      )
      withCallingHandlers(write.csv(emp, outFiles[3], row.names = FALSE), 
                          error= function(e){
                            stop(paste("Cannot save output file", outFiles[1]))
                          }
      )
      withCallingHandlers(write.csv(MGBTpvals, outFiles[4], row.names = FALSE), 
                          error= function(e){
                            stop(paste("Cannot save output file", outFiles[1]))
                          }
      )
      withCallingHandlers(write.csv(trendResults, outFiles[5], row.names = FALSE), 
                          error= function(e){
                            stop(paste("Cannot save output file", outFiles[1]))
                          }
      )
      
      
      
      
       
    #}
    
    
    
    if(plotext == "svg"){
      
      mapply(function(plotFile, plot){svglite::svglite(plotFile,
                                                       width = 9, 
                                                       height = 6,
                                                       fix_text_size = FALSE)
        print(plot)
        dev.off()},
        plotFiles,
        ffaPlots)
      
      mapply(function(plotFile, plot){svglite::svglite(plotFile,
                                                       width = 9, 
                                                       height = 6,
                                                       fix_text_size = FALSE)
        print(plot)
        dev.off()},
        dataplotFiles,
        dataPlots)
      
      mapply(function(plotFile, plot){svglite::svglite(plotFile,
                                                       width = 9, 
                                                       height = 6,
                                                       fix_text_size = FALSE)
        print(plot)
        dev.off()},
        emaplotFiles,
        emaPlots)
    }
    else if(plotext == "png"){
      mapply(function(plotFile, plot){png(plotFile,
                                          width = 9,
                                          height = 6,
                                          units = "in",
                                          res = 150, 
                                          bg = "white")
        print(plot)
        dev.off()},
        plotFiles,
        ffaPlots)
      
      mapply(function(plotFile, plot){png(plotFile,
                                          width = 9,
                                          height = 6,
                                          units = "in",
                                          res = 150, 
                                          bg = "white")
        print(plot)
        dev.off()},
        dataplotFiles,
        dataPlots)
      
      mapply(function(plotFile, plot){png(plotFile,
                                          width = 9,
                                          height = 6,
                                          units = "in",
                                          res = 150, 
                                          bg = "white")
        print(plot)
        dev.off()},
        emaplotFiles,
        emaPlots)
    }
    else if(plotext == "jpg" | plotext == "jpeg" ){
      mapply(function(plotFile, plot){jpeg(plotFile,
                                           width = 9,
                                           height = 6,
                                           units = "in",
                                           res = 150)
        print(plot)
        dev.off()},
        plotFiles,
        ffaPlots)
      
      mapply(function(plotFile, plot){jpeg(plotFile,
                                           width = 9,
                                           height = 6,
                                           units = "in",
                                           res = 150)
        print(plot)
        dev.off()},
        dataplotFiles,
        dataPlots)
      
      mapply(function(plotFile, plot){jpeg(plotFile,
                                           width = 9,
                                           height = 6,
                                           units = "in",
                                           res = 150)
        print(plot)
        dev.off()},
        emaplotFiles,
        emaPlots)
    }
    
  }
  
  #Create log file text
  logText <- logText(inPSF=inPSF, specs=specs, results_info=info, results_qnt = qnt, outFiles = outFiles, warnList = wText)
  
  if(!is.na(outFile) & dir.exists(outDir)){
    
    withCallingHandlers(write(logText, paste0(outDir, "/", outFile, "_log.txt")), 
                        error= function(e){
                          stop(paste("Cannot save output file", outFiles[1]))
                        }
    )
    
    
  }
  
  
  outList <- list(info, qnt, emp, MGBTpvals, trendResults, ffaPlots, dataPlots, emaPlots, logText)
  names(outList) <- c("lp3", "qnt", "emp", "mgb", "trd", "ffa_plots", "inp_plots", "ema_plots", "log")
  
  return(outList)
  
}





