
#' Build MOVE3 Frame
#'
#' Function that converts data from an RDB file to an R dataframe for MOVE3 analysis
#'
#' @param inRDB character - path to WATSTORE file.
#' @param excludeCodes character - vector peak discharge qualifier codes for which peaks should be excluded. See NWIS web help for list of qualifier codes: \cr
#' (https://nwis.waterdata.usgs.gov/nwis/peak?help#flow_qual_cd).
#'
#' @return Dataframe containing peak flow data for MOVE3 analysis.
#'
#' @importFrom stats complete.cases
#'
#' @keywords internal
#'
buildMOVE3frame_RDB <- function(inRDB, excludeCodes=c("3", "4", "8", "A", "O")){
  
  MOVE3frame <- data.frame(matrix(nrow=0, ncol=1))
  colnames(MOVE3frame) <- c("peak_WY")
  
  peakData <- peakfq::readRDB(inRDB) #Read input data from file
  
  peakData <- peakData[!grepl(paste(excludeCodes, collapse="|"), peakData$peak_cd), ] #Get peaks excluding those with certain qualifiers
  peakData <- peakData[stats::complete.cases(peakData$peak_va),] #Remove rows without peak flow value
  
  dupPeaks <- peakData[duplicated(peakData[, c("site_no", "peak_WY")]), c("site_no", "peak_WY")]
  
  if(nrow(dupPeaks) > 0){
    stop(paste("Multiple peak values for the same water year and site: \n"), 
         paste(capture.output(print(dupPeaks)), collapse = "\n"))
  }
  
  
  sites <- unique(peakData$site_no) #Get vector of site numbers
  
  #Loop over the sites and add them to the data frame
  for (site in sites){
    
    siteData <- peakData[peakData$site_no == site,] #Get peaks for the current site
    sitePeaks <- siteData[, c("peak_WY", "peak_va")] #Create dataframe with only site peak values and years
    sitePeaks$peak_va <- as.numeric(as.character(sitePeaks$peak_va)) #Convert peak value to numeric
    colnames(sitePeaks) <- c("peak_WY", as.character(site)) #Name column with site number
    
    MOVE3frame <- merge(sitePeaks, MOVE3frame, by="peak_WY", all=TRUE) #Add column from current station to MOVE3 dataframe
    
  }
  
  
  row.names(MOVE3frame) <- MOVE3frame[,1] #Set row names to first column
  
  MOVE3frame <- MOVE3frame[, 2:ncol(MOVE3frame), drop=FALSE] #Remove first column
  
  return(MOVE3frame)
  
}

#' MOVE3
#'
#' Function to run a MOVE3 analysis on a site
#'
#' @param peaks dataframe - peak values from the sites to be used in analysis, created by the buildMOVE3frame function. Row names are years of peak values, column names are gage numbers used in analysis.
#' @param site_no character - site number of the site for which an extended record is to be generated with MOVE3 (target site).
#' @param begyr integer - beginning year of MOVE3 analysis. Defaults to first year for which data are available.
#' @param endyr integer - end year of MOVE3 analysis. Defaults to last year for which data are available.
#' @param yrmin integer - minimum number of years of overlap required between target and index sites.
#' @param rhomin numeric - minimum correlation required between target and index sites.
#' @param Bulletin17C logical - whether or not Bulletin 17C's modified MOVE3 procedure should be used. Default is FALSE.
#' @param controlBulletin17Cskew logical - whether or not to attempt to control the skew of a Bulletin 17C MOVE3 analysis following the method of Siefken and McCarthy (2022). Default is TRUE.
#' @param fixneg logical - whether or not to increase the number of peak flow values estimated so that the numerator of Bulletin 17C equation 8-24 is not negative
#'
#' @return A list containing the following:
#' 1. First element is a matrix containing the results of MOVE3 analysis for the target site.
#' 2. Second element is a dataframe containing diagnostic information.
#' 3. Third element is a dataframe containing summary information for the target site and index sites used to synthesize peaks.
#' 4. Fourth element is a matrix of peaks from sites meeting the criteria for number of overlap years and correlation with the target site.
#'
#' @examples
#'
#' #Extend the record at NWIS site 06210000 with data from 06209500
#'
#' #Read peak-flow data from NWIS
#' require(dataRetrieval)
#' NWISdata <- readNWISpeak(c("06209500", "06210000"))
#'
#' #Vector peak discharge qualifier codes for which peaks should be excluded.
#' excludeCodes <- c("3", "4", "8", "A", "O")
#'
#' #Build correctly formatted dataframe
#' MOVE3frame <- data.frame(matrix(nrow=0, ncol=1))
#' colnames(MOVE3frame) <- c("peak_WY")
#'
#' #Get peaks excluding those with certain qualifiers
#' peakData <- NWISdata[!grepl(paste(excludeCodes, collapse="|"), NWISdata$peak_cd), ]
#'
#' #Remove rows without peak flow value
#' peakData <- peakData[stats::complete.cases(peakData$peak_va),]
#'
#' #Compute peak flow water year
#' peakData$peak_WY <-calcWaterYear(peakData$peak_dt)
#'
#' #Get vector of site numbers
#' sites <- unique(peakData$site_no)
#'
#' #Loop over the sites and add them to the data frame
#' for (site in sites){
#'
#'   #Get peaks for the current site
#'   siteData <- peakData[peakData$site_no == site,]
#'
#'   #Create dataframe with only site peak values and years
#'   sitePeaks <- siteData[, c("peak_WY", "peak_va")]
#'
#'   #Convert peak value to numeric
#'   sitePeaks$peak_va <- as.numeric(as.character(sitePeaks$peak_va))
#'
#'   #Name column with site number
#'   colnames(sitePeaks) <- c("peak_WY", as.character(site))
#'
#'   #Add column from current station to MOVE3 dataframe
#'   MOVE3frame <- merge(sitePeaks, MOVE3frame, by="peak_WY", all=TRUE)
#'
#' }
#'
#' #Set row names to first column
#' row.names(MOVE3frame) <- MOVE3frame[,1]
#'
#' #Remove first column
#' MOVE3frame <- MOVE3frame[, 2:ncol(MOVE3frame), drop=FALSE]
#'
#' #Run MOVE3
#' MOVE3(MOVE3frame, "06210000", endyr=2018)
#'
#' @export
#'
#' @importFrom stats cor sd var
#'
MOVE3 <- function(peaks, site_no, begyr="", endyr="", yrmin=10, rhomin=0.80, Bulletin17C=TRUE, controlBulletin17Cskew=TRUE, fixneg=FALSE){
  
  #controlBulletin17Cskew=FALSE #Disable skew control in version 1.0
  
  #Set year range
  if (!is.numeric(begyr)) {
    begyr <- min(as.integer(row.names(peaks))) #If no beginning year is set, start at first year in dataframe
  }
  if (!is.numeric(endyr)) {
    endyr <- max(as.integer(row.names(peaks),-1)) #If no end year is set, end at last year in dataframe
  }
  
  
  totalgages<-ncol(peaks) #Get number of gages in analysis
  
  if(totalgages < 2){
    message("Input data contains less than 2 sites.")
    return()
  }
  
  stations<-colnames(peaks) #Get vector of station numbers
  sta <- match(site_no, stations) #Get index of site to run analysis on
  
  if(is.na(sta)){
    message(paste("Specified target site number", site_no, "does not match input data site number."))
    return()
  }
  
  stanumber <- colnames(peaks)[sta]
  #print(paste("Target Station:", stanumber)) #Print target station number
  
  
  #Compute rows corresponding to set year range
  years <- rownames(peaks)
  
  begrow <- which(years == min(years[years >= begyr])) #Look up row number corresponding to beginning year
  endrow <- which(years == max(years[years <= endyr])) #Look up row number corresponding to end year
  
  
  
  #Prepare matrix of peak data
  Qmatrix2 <- as.matrix(peaks)
  
  if(any(Qmatrix2 <= 0, na.rm = TRUE)){
    message("All values in MOVE3 analysis must be positive.")
    return()
  }
  
  logQmatrix2 <- log10(Qmatrix2)
  
  imatrix <- logQmatrix2[begrow:endrow,, drop=FALSE] #Select years of interest
  imatrix <- imatrix[, c(sta, (1:ncol(imatrix))[-sta]), drop=FALSE] #Put current station in first column
  gageinfo <- data.frame(colnames(imatrix), row.names=colnames(imatrix)) #Create gage info frame with stations numbers in first column
  colnames(gageinfo) <- c("site")
  gageinfo$numPeaks <- colSums(!is.na(imatrix)) #Add number of peaks for each gage to gage info
  
  
  if(gageinfo$numPeaks[1] < yrmin){
    message("Target site has fewer years of record than selected minimum years of overlapping peaks. Unable to synthesize peaks.")
    return(list(NULL, gageinfo, NULL, NULL)) #Return the gageorder matrix and NULL for everything else, since no analysis can be performed
  }
  
  if(max(gageinfo$numPeaks[-1]) < yrmin){
    message("All index sites have fewer years of record than selected minimum years of overlapping peaks. Unable to synthesize peaks.")
    return(list(NULL, gageinfo, NULL, NULL)) #Return the gageorder matrix and NULL for everything else, since no analysis can be performed
  }
  
  
  #Compute correlation with target gage
  corMatrix2 <- stats::cor(imatrix,method="pearson",use="pairwise.complete.obs") #Compute correlation matrix of gages
  gageinfo$pearson <- corMatrix2[1:nrow(corMatrix2),1] #Add correlation with target site to gage info
  
  if(sum(is.na(imatrix[,1])) < 2){
    message("All index gages have less than 2 peak flow values to extend target gage record. Unable to synthesize peaks.")
    return(list(NULL, gageinfo, NULL, NULL)) #Return the gageorder matrix and NULL for everything else, since no analysis can be performed
  }
  
  gageinfo$concurrentPeaks <- colSums(!is.na(imatrix[!is.na(imatrix[,1]),])) #Count number of overlapping peaks with target gage for each index gage
  gageinfo$additionalPeaks <- colSums(!is.na(imatrix[is.na(imatrix[,1]),])) #Count number of peaks at the index gage that are not at the target gage
  
  gageinfo2 <- gageinfo[gageinfo$concurrentPeaks >= yrmin,] #Select gages with sufficient years of overlap
  gageinfo3 <- gageinfo2[gageinfo2$pearson >= rhomin,] #Select gages with sufficient correllation
  gageinfo3 <- gageinfo3[gageinfo3$additionalPeaks >= 2 | gageinfo3$site == gageinfo3$site[1],] #Select only the target gage and index gages with at least 2 additional peaks
  gageinfo3$rank <- rank(-gageinfo3[,"pearson"]) #Rank the correlation of the gages from greatest to least
  
  gageorder <- gageinfo3[order(gageinfo3$rank),] #Order based on correllation
  
  
  #Print error and return if no index gages meet criteria
  if (nrow(gageorder) < 2){
    message("No index gages meet correlation and length of overlap criteria. Unable to synthesize peaks.")
    return(list(NULL, gageinfo, NULL, NULL)) #Return the gageorder matrix and NULL for everything else, since no analysis can be performed
  }
  
  
  imatrix2 <- imatrix[,row.names(gageorder), drop=FALSE] #Create new matrix with index gages sorted by correlation
  
  
  
  #Need to select highest correllation index station which will fill AT LEAST 2 years of record
  imatrix3 <- imatrix2
  imatrix3[is.na(imatrix2[,1]), -1] <- NA #Set all peaks not concurrent with target peaks to NA
  sitesRemoved <- c() #Create an empty vector to store the numbers of any index gages removed because they don't fill at least two years
  for (col in 2:ncol(imatrix2)){
    
    if(length(imatrix2[rowSums(!is.na(imatrix2[,1:col-1, drop=FALSE])) == 0 & !is.na(imatrix2[,col]), col]) > 1){
      #SAS 2023-09-11 Fixed bug in logic above where sometimes sites with n < 2
      #were not thrown out, as the counting of additional peaks was incorrect.
      #Copying of highest correlated peaks below was correct, so bug only
      #resulted in occasional crashes if a site needed to be thrown out
      #and wasn't. Did not cause computation errors.
      imatrix3[rowSums(!is.na(imatrix3[,1:col-1, drop=FALSE])) == 0, col] <- imatrix2[rowSums(!is.na(imatrix3[,1:col-1, drop=FALSE])) == 0, col]
      
    }
    else{
      message(paste("Skipping index gage selected for less than two peaks-", colnames(imatrix2)[col]))
      sitesRemoved <- append(sitesRemoved, colnames(imatrix2)[col]) #Add to list of sites removed
    }
  }
  
  #Remove all index sites not selected because they couldn't fill 2 years
  imatrix3 <- imatrix3[, colnames(imatrix3) %!in% sitesRemoved]
  gageorder <- gageorder[rownames(gageorder) %!in% sitesRemoved,]
  
  
  
  conmatrix3 <- matrix(nrow=nrow(imatrix3), ncol=ncol(imatrix3), dimnames=list(row.names(imatrix3))) #Create matrix for index peaks concurrent with target peaks
  nonconmatrix3 <- matrix(nrow=nrow(imatrix3), ncol=ncol(imatrix3), dimnames=list(row.names(imatrix3), colnames(imatrix3))) #Create matrix for index peaks not concurrent with target peaks
  ymatrix3 <- matrix(nrow=nrow(imatrix3), ncol=ncol(imatrix3), dimnames=list(row.names(imatrix3), colnames(imatrix3))) #Create matrix of target peaks concurrent with index peaks
  conmatrix3[!is.na(imatrix3[,1]),] <- imatrix3[!is.na(imatrix3[,1]),] #Select index peaks concurrent with target peaks
  nonconmatrix3[is.na(imatrix3[,1]),] <- imatrix3[is.na(imatrix3[,1]),] #Select index peaks not concurrent with target peaks
  
  #Select target peaks concurrent with each index peak
  for (col in 1:ncol(ymatrix3)){
    ymatrix3[!is.na(conmatrix3[,col]), col] <- imatrix3[!is.na(conmatrix3[,col]), 1]
  }
  
  
  # Equation number references below refer to equations in Bulletin 17C Appendix 8
  gageorder$logmean_y1 <- colMeans(ymatrix3, na.rm = TRUE, dims = 1) # EQ 8-1
  gageorder$logmean_x1 <- colMeans(conmatrix3, na.rm = TRUE, dims = 1) # EQ 8-2
  gageorder$logmean_x2 <- colMeans(nonconmatrix3, na.rm = TRUE, dims = 1) # EQ 8-3
  gageorder$logvar_y1 <- colVars(ymatrix3) # EQ 8-4
  gageorder$logvar_x1 <- colVars(conmatrix3) # EQ 8-5
  gageorder$logvar_x2 <- colVars(nonconmatrix3) # EQ 8-6
  
  #Compute standard deviations from variances
  sdy <- sqrt(gageorder$logvar_y1[-1])
  sdx1 <- sqrt(gageorder$logvar_x1[-1])
  sdx2 <- sqrt(gageorder$logvar_x2[-1])
  
  betam <- matrix()
  
  for(i in 2:ncol(ymatrix3)){
    
    betam[i-1] <- sum((conmatrix3[,i]-gageorder$logmean_x1[i])*(conmatrix3[,1]-gageorder$logmean_y1[i]), na.rm=TRUE)/sum((conmatrix3[,i]-gageorder$logmean_x1[i])^2,na.rm=TRUE) # EQ 8-10
    
    gageorder$additionalPeaks[i] <- sum(!is.na(nonconmatrix3[,i])) # Replace number of nonconcurrent peaks with number of peaks filled
    
    
  }
  
  
  
  n1m <- gageorder$concurrentPeaks[-1]
  n2m <- gageorder$additionalPeaks[-1]
  ybar <- gageorder$logmean_y1[-1]
  x1bar <- gageorder$logmean_x1[-1]
  x2bar <- gageorder$logmean_x2[-1]
  vary <- sdy^2
  varx1 <- sdx1^2
  varx2 <- sdx2^2
  
  
  
  alpha2 <- (n2m*(n1m-4)*(n1m-1))/((n2m-1)*(n1m-3)*(n1m-2))  # EQ 8-11
  rhom <- betam * sdx1 / sdy # EQ 8-9
  sigma2 <- (1 / (n1m + n2m - 1)) * ((n1m - 1) * vary + (n2m - 1) * betam ^ 2 * varx2 + (n2m - 1) * alpha2 * (1 - rhom ^ 2) * vary + (n1m * n2m / (n1m + n2m)) * betam ^ 2 * (x2bar - x1bar) ^ 2) #EQ 8-8
  muy <- ybar + n2m / (n1m + n2m) * betam * (x2bar - x1bar) # EQ 8-7
  
  
  bigA <- ((n2m + 2) * (n1m - 6) * (n1m - 8)) / (n1m - 5) + (n1m - 4) * ((n1m * n2m * (n1m - 4)) / ((n1m - 3) * (n1m - 2)) - (2 * n2m * (n1m - 4) / (n1m - 3)) - 4) # EQ 8-14
  bigB <- (6 * (n2m + 2) * (n1m - 6)) / (n1m - 5) + 2 * (n1m ^ 2 - n1m - 14) + (n1m - 4) * ((2 * n2m * (n1m - 5)) / (n1m - 3) - 2 * (n1m + 3) - (2 * n1m * n2m * (n1m - 4)) / ((n1m - 3) * (n1m - 2))) # EQ 8-15
  bigC <-2 * (n1m + 1) + 3 * (n2m + 2) / (n1m - 5) - ((n1m + 1) * (2 * n1m + n2m - 2) * (n1m - 3)) / (n1m - 1) + (n1m - 4) * (2 * n2m / (n1m - 3) + 2 * (n1m + 1) + (n1m * n2m * (n1m - 4)) / ((n1m - 3) * (n1m - 2))) # EQ 8-16
  
  var_sigma2<- 2 * vary ^ 2 / (n1m - 1) + (n2m * vary ^ 2) / ((n1m + n2m - 1) ^ 2 * (n1m - 3)) * (bigA * rhom ^ 4 + bigB * rhom ^ 2 + bigC) #EQ 8-13
  
  nen1_26 <- 2 / ((2 / (n1m - 1) + n2m / ((n1m + n2m - 1) ^ 2 * (n1m - 3)) * (bigA * rhom ^ 4 + bigB * rhom ^ 2 + bigC))) + 1 #EQ 8-19
  
  ne <- nen1_26 - n1m #Use ne computed according to Bulletin 17C equation 8-19
  
  # Compute minimum years of record for B17C MOVE3 to work
  # Equations to be published in regional methods SIR
  momdiffs <- (n1m -1)*(sigma2 - vary) - n1m*(muy - ybar)^2
  nmin <- (-1*momdiffs + sqrt(momdiffs^2 + 4*n1m^2*sigma2*(muy - ybar)^2))/(2*sigma2)
  
  
  #Generate predictions
  pred3 <- 10^(as.matrix(conmatrix3[,1])) #First column is recorded peaks at target gage
  
  if(!Bulletin17C){
    ns <- n2m #Set number of peaks to synthesize to ne
    #Compute coefficients to extend record by n2 years
    aprime <- ((n1m + n2m) * muy - n1m * ybar) / n2m # EQ 8-23
    b <-(((n1m + n2m - 1) * sigma2 - (n1m - 1) * vary - n1m * (ybar - muy) ^ 2 - n2m * (aprime - muy) ^ 2) / ((n2m - 1) * varx2))^0.5 # EQ 8-24
    
    for(i in 2:ncol(nonconmatrix3)){
      x<-10^(aprime[i-1]+b[i-1]*(nonconmatrix3[,i]-x2bar[i-1])) #Compute synthesized peaks for each index gage
      pred3<-cbind(pred3,x)
    }
  }
  else{
    #Check if Bulletin 17C extension is mathematically feasible
    ne_round <- round(ne) #Round ne to the nearest integer
    ns <- ne_round #Set number of peaks to synthesize to ne
    
    
    k <- (n1m + ne_round - 1) * sigma2 - (n1m - 1) * vary - n1m*(n1m/ne_round + 1)*(ybar - muy)^2
    
    if(any(k < 0) & !fixneg){
      message("Numerator of Bulletin 17C EQ 8-24 is negative for at least one index site. Cannot run analysis.")
      message(paste("Minimum years of record: ", nmin))
      return(list(NULL, gageinfo, NULL, NULL)) #Return the gageorder matrix and NULL for everything else, since no analysis can be performed
    }
    else if(any(k < 0) & fixneg){
      
      #Set number of peaks to synthesize to the max of ne and the minimum number of peaks
      ns <- pmax(ceiling(nmin), ne_round)
      message(paste0("Numerator of Bulletin 17C EQ 8-24 is negative using ne = ",
                     paste(ne_round, collapse = ",")))
      message(paste0("Using minimum required synthesized peaks nmin = ",
                     paste(ns, collapse = ",")))
      
    }
    
    #Recompute k using ns, the actual number of peaks synthesized
    k <- (n1m + ns - 1) * sigma2 - (n1m - 1) * vary - n1m*(n1m/ns + 1)*(ybar - muy)^2
    
    aprime <- ((n1m + ns) * muy - n1m * ybar) / ns # EQ 8-23
    
    
    #Select ne index peaks to use for record extension
    indexPeaks <- nonconmatrix3
    xebar <- rep(NA, (ncol(nonconmatrix3)-1)) #Create vector to store mean of index peaks
    varxe <-  rep(NA, (ncol(nonconmatrix3)-1)) #Create vector to store variance of index peaks
    
    if(!controlBulletin17Cskew){
      for(i in 2:ncol(nonconmatrix3)){
        #Select which ne peaks at the index site to use record extension
        indexPeaks[!is.na(indexPeaks[,i]),][1:(n2m[i-1]-ns[i-1]), i] <- NA #Set all but newest ne peaks to NA
        #print(indexPeaks)
        xebar[i-1] <- mean(indexPeaks[,i], na.rm=TRUE)
        varxe[i-1] <- stats::var(indexPeaks[,i], na.rm=TRUE)
      }
    }
    
    #print(paste("varxe", varxe))
    
    else{
      
      
      
      
      for(i in 2:ncol(nonconmatrix3)){
        
        if(ns[i-1] < 3){
          message(paste("Number of synthesized peaks is less than 3: ", ns[i-1]))
          return(list(NULL, gageinfo, NULL, NULL)) #Return the gageorder matrix and NULL for everything else, since no analysis can be performed
        }
        
        
        #Testing on skewness of Bulletin 17C MOVE3
        ypeaks <- ymatrix3[,i] #Select the peaks at the target site that overlap with this index site
        ypeaks <- ypeaks[!is.na(ypeaks)]
        x2peaks <- nonconmatrix3[!is.na(nonconmatrix3[,i]),i]
        T1 <- sum((ypeaks - muy[i-1])^3) #Compute the third moment
        
        xe_best <- skewCheck(n1 = n1m[i-1], n2 = n2m[i-1], ne = ns[i-1], ybar1 = ybar[i-1], sy1=sdy[i-1] , muy=muy[i-1], sigmay2=sigma2[i-1], T1=T1, x2=x2peaks)
        indexPeaks[rownames(indexPeaks) %!in% names(xe_best), i] <- NA #Select set peaks not in the xe best peaks to NA
        
        
        xebar[i-1] <- mean(xe_best, na.rm=TRUE)
        varxe[i-1] <- stats::var(xe_best, na.rm=TRUE)
      }
    }
    
    
    
    #Compute coefficients to extend record by ne years (Bulletin 17C method)
    
    b <-(((n1m + ns - 1) * sigma2 - (n1m - 1) * vary - n1m * (ybar - muy) ^ 2 - ns * (aprime - muy) ^ 2) / ((ns - 1) * varxe))^0.5 # EQ 8-24
    
    for(i in 2:ncol(nonconmatrix3)){
      x <- 10^(aprime[i-1]+b[i-1]*(indexPeaks[,i]-xebar[i-1])) #Compute synthesized peaks for each index gage
      pred3 <- cbind(pred3,x)
    }
    
  }
  
  
  
  
  
  
  #Prediction error calculations
  sep<-matrix()
  for(i in 2:ncol(conmatrix3)){
    x<-conmatrix3[,1]-(aprime[i-1]+b[i-1]*(conmatrix3[,i]-x1bar[i-1])) #Computes difference between recorded peaks and synthesized peaks for concurrent years
    sep[i-1] <- stats::sd(x,na.rm=TRUE) #Compute standard error of prediction
    
  }
  sep<-as.matrix(sep)
  
  OLSse<-as.matrix(sum(sep[]*n2m,na.rm=TRUE)/sum(n2m,na.rm=TRUE))
  
  rwtd<-sum(rhom*n2m, na.rm=TRUE)/sum(n2m,na.rm=TRUE)
  
  m3se<-(OLSse)*(2/(1+rwtd))^0.5
  m3se_percent=100*((exp((m3se*2.3026)^2)-1)^0.5)
  
  
  
  sumne<-sum(ne,na.rm=TRUE)
  
  gageorder$SEP[2:nrow(gageorder)] <- sep
  gageorder$Eq_years[2:nrow(gageorder)] <- ne
  gageorder$N_synth[2:nrow(gageorder)] <- ns
  gageorder$mu_y[2:nrow(gageorder)] <- muy
  gageorder$sigma2_y[2:nrow(gageorder)] <- sigma2
  
  colnames(pred3)<-gageorder$site
  
  
  #Output results to dataframe
  pred4 <- data.frame(matrix(nrow=nrow(pred3), ncol=2), row.names=rownames(pred3))
  colnames(pred4) <- c("peak_va", "from_gage")
  for (col in 1:ncol(pred3)){
    pred4[!is.na(pred3[,col]), 1] <- pred3[!is.na(pred3[,col]),col] #Copy the peaks from the current gage
    pred4[!is.na(pred3[,col]), 2] <- colnames(pred3)[col] #Fill in the gage number
    
    #Identify years of nonconcurrent data from each site
    siteYears <- rownames(pred3[!is.na(nonconmatrix3[,col]),col, drop=FALSE])
    siteYearsText <- peakFrame2Text(peakRanges(as.numeric(siteYears)))
    gageorder$indexYears[col] <- siteYearsText
    
    #Identify years of synthesized data from each site
    siteSynthYears <- rownames(pred3[!is.na(pred3[,col]),col, drop=FALSE])
    siteSynthYearsText <- peakFrame2Text(peakRanges(as.numeric(siteSynthYears)))
    gageorder$synthYears[col] <- siteSynthYearsText
    
  }
  
  
  #Round final results according the USGS criteria
  pred4$peak_va <- roundPeak(pred4$peak_va)
  
  
  #Create summary table for MOVE3
  MOVE3sum <- gageorder[gageorder$additionalPeaks > 0, c("site", "additionalPeaks", "concurrentPeaks", "pearson", "Eq_years", "N_synth", "SEP",  "indexYears", "synthYears")]
  MOVE3sum$pearson <- round(MOVE3sum$pearson, 2)
  MOVE3sum$Eq_years <- round(MOVE3sum$Eq_years, 1)
  MOVE3sum$targetStation <- site_no
  MOVE3sum$targetPeaks <- gageorder[site_no, "concurrentPeaks"]
  #MOVE3sum$synthesizedPeaks <- sum(MOVE3sum$additionalPeaks)
  #MOVE3sum$percentSynthesized <- MOVE3sum$synthesizedPeaks/(MOVE3sum$synthesizedPeaks + MOVE3sum$targetPeaks) * 100 #Compute percentage of peaks synthesized
  MOVE3sum$targetWaterYears <- peakFrame2Text(peakRanges(as.numeric(rownames(peaks[!is.na(peaks[, site_no]),])))) #Get water years of peaks at target site
  MOVE3sum$synthesizedWaterYears <- peakFrame2Text(peakRanges(as.numeric(rownames(pred4[!is.na(pred4$from_gage) & pred4$from_gage != site_no,])))) #Get water years of synthesized peaks
  MOVE3sum$weightedPearson <- round(rwtd, 2) #Added weighted Pearson coefficient
  MOVE3sum$MOVE3StdErrorEst_percent <- round(as.numeric(m3se_percent), 1) #Add estimated MOVE3 standard error, in percent
  
  #print(MOVE3sum)
  
  return(list(pred4, gageorder, MOVE3sum, imatrix2))
  
}





#' Column Variance
#'
#' Compute variance of columns in a dataframe. NA values are automatically removed.
#'
#' @param x dataframe
#'
#' @return numeric - vector of variance values for each column. NA values are automatically removed.
#'
#' @keywords internal
colVars <- function(x) {
  
  n <- colSums(!is.na(x))
  
  variances <- (colMeans(x*x, na.rm=TRUE) - (colMeans(x, na.rm=TRUE))^2)* n/(n-1)
  
  return(variances)
  
}




## Function below is not used in version 1.0

#' SKew Check
#'
#' Function to evaluate technique for controlling skew of MOVE3 record extension
#'
#' @param n1 integer - as defined in Bulletin 17C appendix 8
#' @param n2 integer - as defined in Bulletin 17C appendix 8
#' @param ne numeric - as defined in Bulletin 17C appendix 8
#' @param ybar1 numeric - as defined in Bulletin 17C appendix 8
#' @param sy1 numeric - as defined in Bulletin 17C appendix 8
#' @param muy numeric - as defined in Bulletin 17C appendix 8
#' @param sigmay2 numeric - as defined in Bulletin 17C appendix 8
#' @param T1 numeric - sum from i = 1 to n1 of (yi - y1bar) where yi is the ith value at the target site in the overlap period
#' @param x2 numeric vector - x2 data series from the index site
#' @param maxcomb numeric - maximum number of combinations to evaluate
#'
#' @return numeric vector containing the xe values for which the skew of the BUlletin 17C MOVE3 record extension will most closely match the skew of a full MOVE3 record extension
#'
#' @importFrom e1071 skewness
#' @importFrom utils combn
#'
#' @keywords internal
skewCheck <- function(n1, n2, ne, ybar1, sy1, muy, sigmay2, T1, x2, maxcomb = 1e6){
  
  if(ne < 3){
    stop(paste("Number of synthesized peaks is less than 3: ", ne))
  }
  
  
  # Compute k using n2 (EQ XX Siefken and McCarthy)
  k_n2 <- (n1+n2-1)*sigmay2 - (n1-1)*sy1^2 - n1*(n1/n2 + 1)*(muy - ybar1)^2
  
  # Compute k using ne (EQ XX Siefken and McCarthy)
  k_ne <- (n1+ne-1)*sigmay2 - (n1-1)*sy1^2 - n1*(n1/ne + 1)*(muy - ybar1)^2
  
  skew_x2 <- e1071::skewness(x2, type=2) #Need to use type 2 to match Bulletin 17C definition
  
  # Compute skewness of MOVE3 extension using the full n2 years of additional record
  skew_full <- (n1+n2)/(sigmay2^1.5*(n1+n2-1)*(n1+n2-2))*(T1 + n2*(n1/n2*(muy-ybar1))^3 + 3*k_n2*n1/n2*(muy-ybar1) + (n2-1)*(n2-2)/n2*(k_n2/(n2-1))^1.5*skew_x2)
  
  
  # Compute the optimal skewness of the additional ne peaks (that required for the extended record to exactly match a full MOVE3 extension)
  skew_xe_opt <- ne/((ne-1)*(ne-2))*(k_ne/(ne-1))^(-1.5)*(skew_full*sigmay2^1.5*(n1+ne-1)*(n1+ne-2)/(n1+ne) - T1 - ne*(n1/ne*(muy-ybar1))^3 - 3*k_ne*(n1/ne)*(muy-ybar1))
  
  
  
  #Compute possible skews of xe additional peaks
  allComb <- choose(n2, ne) #compute n choose r with n = n2 and r = ne
  consecComb <- n2 - ne + 1 #compute combinations of consecutive peaks
  
  
  if (allComb < maxcomb){
    
    xeoptions <- utils::combn(x2, ne, simplify = FALSE)
    
  }
  
  else{
    #If too many combinations are possible, we'll restrict the options to a series of consecutive peaks
    xeoptions <- list(c())
    
    for (i in 1:(length(x2)-ne)){
      xeoptions[[i]] <- x2[i:(i+ne-1)]
    }
    
  }
  
  #Remove any combinations where all values are equal, since variance will be zero, resulting in division by zero
  xeoptions <- Filter(function(x) length(unique(x)) > 1, xeoptions)
  
  
  xeSkews <- sapply(xeoptions, e1071::skewness, type=2) #Compute skewness for each combination
  
  
  # print(paste("Minimum xe skew:", min(xeSkews)))
  # print(paste("Maximum xe skew:", max(xeSkews)))
  
  # Determine the best available skewness of the additional ne peaks
  xeSkew_best <- xeSkews[abs(xeSkews - skew_xe_opt) == min(abs(xeSkews - skew_xe_opt))]
  
  # print(paste("Best possible xe skew:", xeSkew_best))
  
  #Get the x2 values corresponding to the optimal skew
  xe_best <- xeoptions[abs(xeSkews - skew_xe_opt) == min(abs(xeSkews - skew_xe_opt))]
  
  
  
  # Compute the skewness of the extended recorded resulting from using the best combination of ne peaks
  skew_y_best <- (n1+ne)/(sigmay2^1.5*(n1+ne-1)*(n1+ne-2))*(T1 + n1^3/ne^2*(muy-ybar1)^3 + 3*k_ne*n1/ne*(muy-ybar1) + (ne-1)*(ne-2)/ne*(k_ne/(ne-1))^1.5*xeSkew_best)
  
  # print(paste("Best possible extended series skew:", skew_y_best))
  
  return(xe_best[[1]])
  
  
}



#' Peak Ranges
#'
#' Determines the water year ranges for which peak flows are available
#'
#' @param peakYears numeric - vector containing the water years for which peak flows are available.
#'
#'
#' @return Dataframe with the first and second columns containing the start water years and water end years respectively of periods for which peak flows are available.
#'
#' @examples
#'
#' peaks <- c(1982, 1983, 1984, 1990, 1991, 1992, 2000) #Years of peak-flow data
#' ranges <- peakfq:::peakRanges(peaks)
#' print(ranges)
#'
#' @keywords internal
peakRanges <- function(peakYears){
  
  peakYears <- sort(peakYears, decreasing=FALSE) #Make sure the years are sorted in ascending order
  
  breaks <- 0 #Counter for number of breaks in the data
  
  #Create Vectors for storing the start and end years of the ranges
  rangeStarts <- min(peakYears)
  rangeEnds <- max(peakYears)
  
  if (length(peakYears) > 1){
    
    #Loop over the years
    for (n in seq(2,length(peakYears))){
      
      if ((peakYears[n] - peakYears[n-1]) != 1) {
        
        breaks <- breaks + 1
        rangeEnds <- append(rangeEnds, peakYears[n-1], length(rangeEnds) - 1) #Append as second to last, since last is already included
        rangeStarts <- append(rangeStarts, peakYears[n])
        
      }
    }
  }
  
  ranges <- data.frame(startYears=rangeStarts, endYears=rangeEnds) #Convert vectors to data frame
  
  return(ranges)
  
}

#' Peak Years Text
#'
#' Converts data frame of year ranges for which peak flows are available to human-readable text i.e. "1923-1928, 1935, 1945-2004".
#'
#' @param peakFrame dataframe - consists of two columns, the first column containing the years for which a range of peak values starts, the second column containing the end years of the ranges.
#' This dataframe can be created from the `peakRanges()` function.
#'
#' @return Character vector containing the ranges of water years for which peak flows are available in human readable format.
#'
#' @examples
#' peaks <- c(1982, 1983, 1984, 1990, 1991, 1992, 2000) #Years of peak-flow data
#' ranges <- peakfq:::peakRanges(peaks)
#' peakfq:::peakFrame2Text(ranges)
#'
#' @keywords internal
peakFrame2Text <- function(peakFrame){
  
  rangeText <- c("") #Create empty string
  
  for (n in 1:nrow(peakFrame)){
    
    if (peakFrame[n,1] == peakFrame[n,2]) {
      #If raange is only one year, append only one year
      rangeText <- paste(rangeText, peakFrame[n,1], ", ", sep="")
    }
    
    else {
      rangeText <- paste(rangeText, peakFrame[n,1], "-", peakFrame[n,2], ", ", sep="")
    }
    
    
  }
  
  rangeText <- substr(rangeText, 1, nchar(rangeText) - 2) #Remove last comma and space
  
  return(rangeText)
  
}


#' Round Peak
#'
#' Rounds peak flow values in accordance with the rounding criteria in Rantz and others (1982)
#'
#' @param peaks numeric - vector of peak flow values to be rounded.
#'
#' @return Numeric vector containing the rounded peak flow values.
#'
#' @examples
#' peakFlows <- c(456.23, 12.6532, 1.0325, 14562, 9823) #Vector of peak-flow values
#' roundPeak(peakFlows) #Round the peak-flow values
#'
#' @references
#'
#' Rantz, S.E. and others, 1982, Measurement and Computation of Streamflow: Vol. 2 Computation of Discharge: U.S. Geological Survey Water Supply Paper 2175, https://doi.org/10.3133/wsp2175
#'
#' @export
roundPeak <- function(peaks){
  
  peaks[peaks<1 & !is.na(peaks)] <- round(peaks[peaks<1 & !is.na(peaks)], 2)#Flows less than 1 cfs are rounded to the nearest hundreth of a cfs
  peaks[peaks>=1 & peaks<10 & !is.na(peaks)] <- round(peaks[peaks>=1 & peaks<10 & !is.na(peaks)], 1)#Flows from 1 to 10 cfs are rounded to the tenth
  peaks[peaks>=10 & peaks<1000 & !is.na(peaks)] <- round(peaks[peaks>=10 & peaks<1000 & !is.na(peaks)], 0)#Flows from 10 to 1000 cfs are rounded to the nearest cfs
  peaks[peaks>= 1000 & !is.na(peaks)] <- signif(peaks[peaks>= 1000 & !is.na(peaks)], 3)#Flow of 1,000 cfs or greater are rounded to 3 significant figures
  
  return(peaks)
  
}













#' MOVE3 from RDB File
#'
#' Function to run a MOVE3 analysis for specified site in an RDB peak data file
#'
#' @param inRDB character - path to RDB file containing data to be used in the analysis.
#' @param outDir character - path to directory for output files.
#' @param targetGages character - vector of  target site numbers for which to perform MOVE3 record extension.
#' @param includeGages character - vector of site numbers to be included in the analysis as index sites. All target sites will also be considered as index sites. Default will use all gages in WATSTORE file.
#' @param begyr integer - beginning year of MOVE3 analysis. Defaults to first year for which data are available.
#' @param endyr integer - end year of MOVE3 analysis. Defaults to last year for which data are available.
#' @param yrmin integer - minimum number of years of overlap required between target and index sites.
#' @param rhomin numeric - minimum correlation required between target and index sites.
#' @param excludeCodes character - vector of peak discharge qualifier codes for which peaks should be excluded from MOVE3 analysis. See NWIS web help for list of qualifier codes (https://nwis.waterdata.usgs.gov/nwis/peak?help#flow_qual_cd).
#' @param Bulletin17C logical - whether or not Bulletin 17C's modified MOVE3 procedure should be used. Defaults to FALSE.
#' @param controlBulletin17Cskew logical - whether or not to attempt to control the skew of a Bulletin 17C MOVE3 analysis following the method of Siefken and McCarthy (2022).
#' @param fixneg logical - whether or not to increase the number of peak flow values estimated so that the numerator of Bulletin 17C equation 8-24 is not negative
#' @param siteInfoFile character - path to tab-delimited site information file. 
#'
#' @return A list containing the following:
#' 1. First element is a character vector of gage numbers for which  additional peaks were successfully synthesized with MOVE3.
#' 2. Second element is a character vector of gage numbers for which additional peaks could not be synthesized.
#'
#'
#' @examples
#'
#' dir.create(file.path(tempdir(), "MOVE3")) #Create directory for output files
#'
#' #Tab-delimited (RDB format) data file containing peak-flow data for target and index sites
#' inputData <- system.file("extdata",
#'                          "B17C_MOVE3example_nwis_peak.txt",
#'                          package="peakfq")
#'
#' #Directory for MOVE3 output files
#' outputFolder <- file.path(tempdir(), "MOVE3")
#'
#' #Site number of target site for which record will extended
#' targetSite <- c("02334885")
#'
#' #Minimum number of years of overlap required between target and index sites.
#' #Bulletin 17C recommends at least 10.
#' minOverlap <- 10
#'
#' #Minimum Pearson correlation coefficient between target and index sites.
#' minCorrel <- 0.8
#'
#' #Codes in WATSTORE file which will cause peaks to be excluded from MOVE3 analysis
#' excludeCodes <- c("3", "4", "8", "A")
#'
#' #Run MOVE3 analysis, setting Bulletin17C = FALSE will result in a full record extension
#' MOVE3_RDB(inputData,
#'             outputFolder,
#'             targetSite,
#'             yrmin = minOverlap,
#'             rhomin = minCorrel,
#'             excludeCodes = excludeCodes,
#'             Bulletin17C = TRUE)
#'
#'
#' @export
#'
#' @importFrom stats qqnorm qqline
#' @importFrom utils combn write.csv
#' @importFrom grDevices dev.off png
#' @importFrom graphics plot
#'
MOVE3_RDB <- function(inRDB, outDir, targetGages, includeGages=c(), begyr="", endyr="", yrmin=10, rhomin=0.80, excludeCodes=c("3", "4", "8", "A", "O"), Bulletin17C=TRUE, controlBulletin17Cskew=TRUE, fixneg=FALSE, siteInfoFile=""){
  
  
  MOVE3input <- buildMOVE3frame_RDB(inRDB, excludeCodes) #Build dataframe from WATSTORE

  #Read in site information if supplied
  siteInfo <- data.frame()
  if(file.exists(siteInfoFile)){
    
    siteInfo <- read.delim(siteInfoFile, comment.char = "#")
    
    #If the first row is not data, drop (applies for RDB files with column definition row) - https://waterdata.usgs.gov/nwis/?tab_delimited_format_info
    if(siteInfo[1,1] == "5s"){
      siteInfo <- siteInfo[-1,]
    }
  }
  
  #Consider revising behavior as follows:
  # 1. Either use ALL sites as index gages or only those listed in includeGages input
  # 2. Don't include all target sites as index sites (I don't remember why it did this originally)
  if(length(includeGages) > 0){
    #If includeGages is used, remove gages not included
    indexGages <- unique(c(targetGages, includeGages)) #Make list of all unique sites to include

    MOVE3input <-MOVE3input[,indexGages] #Select only sites included in analysis
  }
  
  #Create vectors to store gages with peaks successfully synthesized or not
  successGages <- c()
  failGages <- c()
  
  for (gage in targetGages){
    
    MOVE3out <- MOVE3(MOVE3input, gage, begyr=begyr, endyr=endyr, yrmin=yrmin, rhomin=rhomin, Bulletin17C=Bulletin17C, controlBulletin17Cskew=controlBulletin17Cskew, fixneg=fixneg) #Run MOVE3 analysis
    
    if(is.null(MOVE3out[[1]])){
      #If analysis failed to synthesize peaks, only output diagnostics
      message(paste("Failed to synthesize peaks for", gage))
      MOVE3diagnostics <- MOVE3out[[2]] #Get diagnostic info from MOVE3 output
      utils::write.csv(MOVE3diagnostics, paste(outDir,"/", gage, "_MOVE3diag.csv", sep=""))
      failGages <- append(failGages, gage)
    }
    else{
      #If analysis sucessfully synthesized peaks, prepare the output
      MOVE3peaks <- MOVE3out[[1]] #Get peaks from MOVE3 output
      
      MOVE3diagnostics <- MOVE3out[[2]] #Get diagnostic info from MOVE3 output
      MOVE3summary  <- MOVE3out[[3]]
      availablePeaks <- MOVE3out[[4]] #Dataframe of peaks meeting overlap and correlation requirements for MOVE3
      
      #Export MOVE3 generated peaks to text file 
      synthPeaks <- MOVE3peaks[!is.na(MOVE3peaks$from_gage) & MOVE3peaks$from_gage != gage,] #Select only synthesized peaks, need !is.na to clean up rows without synthesized peaks
      synthPeaks$comment <- paste("MOVE3 peak from site", synthPeaks$from_gage)
      
      #Need to identify years used in MOVE3, but not synthesized and output appropriate perception thresholds
      noIndexYears <- rownames(availablePeaks[apply(availablePeaks, 1, function(x) all(is.na(x))), ]) #Select rows with all NA values indicating no gage data available
      notsynthPeaks <- rownames(MOVE3peaks)[is.na(MOVE3peaks$from_gage)]
      notsynthPeaks <- notsynthPeaks[notsynthPeaks %!in% noIndexYears]
      
      if(length(notsynthPeaks) > 0){
        notsynthRanges <- peakRanges(as.integer(notsynthPeaks))
        threshText <- paste("PCPT_Thresh", notsynthRanges$startYears, notsynthRanges$endYears, "1E+20 1E+20 DATA USED IN MOVE3 BUT NOT SYNTHESIZED")
      } 
      else{
        threshText <- c() #Initial as empty vector
      }
      
      intervalText <- paste("Interval", rownames(synthPeaks), synthPeaks$peak_va, synthPeaks$peak_va, synthPeaks$comment)
      
      PSFtext <- c(threshText, intervalText)
      
      
      #Make nice output dataframe
      MOVE3summary$site <- as.character(MOVE3summary$site) #Convert factor to character
      MOVE3summary$targetName <- NA
      MOVE3summary$targetDrainageArea <- NA
      
      indexInfo <- data.frame()
      
      if(nrow(siteInfo) > 0){
        targetGageInfo <- siteInfo[siteInfo$site_no == gage,] #Get site info for current target gage
        
        if(nrow(targetGageInfo) == 1){
          MOVE3summary$targetName <- targetGageInfo$station_nm #Add target site name
          MOVE3summary$targetDrainageArea <- targetGageInfo$drain_area_va #Add target site drainage area
        }
        
        indexInfo <- siteInfo[siteInfo$site_no %in% MOVE3summary$site,] #Get site info for current index gages
       
        if(nrow(indexInfo) > 0){
          indexInfo <- indexInfo[,c("site_no", "station_nm", "drain_area_va")] #Select the relevant information
          colnames(indexInfo) <- c("site", "indexName", "indexDrainageArea") #rename columns
          MOVE3summary <- merge(MOVE3summary, indexInfo, by="site", all=TRUE)#Add index info to summary table
        }
        else{
          MOVE3summary$indexName <- NA
          MOVE3summary$indexDrainageArea <- NA
        }
        
        
      } 
      else{
        MOVE3summary$indexName <- NA
        MOVE3summary$indexDrainageArea <- NA
      }

      
      MOVE3summary <- MOVE3summary[,c("targetStation", "targetName", "targetDrainageArea", "targetPeaks", "targetWaterYears",
                                      "site", "indexName", "indexDrainageArea", "indexYears", "concurrentPeaks", "additionalPeaks",  "pearson", "weightedPearson", "MOVE3StdErrorEst_percent", "Eq_years", "N_synth",
                                       "synthYears")]
      
      
      
      write(PSFtext, paste(outDir,"/", gage, "_MOVE3peaks.txt", sep=""))
      utils::write.csv(MOVE3diagnostics, paste(outDir,"/", gage, "_MOVE3diag.csv", sep=""))
      utils::write.csv(MOVE3summary, paste(outDir,"/", gage, "_MOVE3summary.csv", sep=""))
      utils::write.csv(availablePeaks, paste(outDir,"/", gage, "_MOVE3availablePeaks.csv", sep=""))
      
      
      #Create Normal Probability plots for gages meeting correlation criteria
      for(colName in colnames(availablePeaks)){
        
        grDevices::png(paste0(outDir, "/", colName, "normalPlot.png"))
        stats::qqnorm(availablePeaks[, colName], main=paste(colName, "Normal Probability Plot"))
        stats::qqline(availablePeaks[, colName])
        grDevices::dev.off()
        
        # Prevents generating a scatter plot of a site against itself
        if(colName == gage) {next}
        
        #Plot concurrent peaks at target and index sites
        grDevices::png(paste0(outDir, "/", gage, colName, "-", "scatterPlot.png"))
        graphics::plot(availablePeaks[, 1], availablePeaks[, colName], xlab = gage, ylab = colName, main = "Log-Transformed Concurrent Peaks")
        graphics::abline(stats::lm(availablePeaks[, colName] ~ availablePeaks[, 1]))
        grDevices::dev.off()
        
      }
      
      
      
      
      
      
      successGages <- append(successGages, gage)
    }
    
    
  }
  
  
  return(list(successGages, failGages))
  
}


