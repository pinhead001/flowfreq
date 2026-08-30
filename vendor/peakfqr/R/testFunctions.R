#' Check Frequency Curve
#'
#' Function to check frequency curve against results from PeakFQ 7.5.1
#'
#' @param PSFfile character - path to PeakFQ 7.5.1 specifications file.
#' @param EXPfile character - path to .EXP file exported from PeakFQ 7.5.1 using the specifications in PSFfile.
#' @param excludeSites character - character vector of site numbers to exclude from test
#' @param quietly logical - if `TRUE`, messages from the peakfq are not printed
#'
#' @return None
#' 
#' @keywords internal
#'
#' @details
#' Runs the supplied specifications file in with the `peakfq()` function and
#' then runs tests comparing the output to supplied output from PeakFQ 7.5.1. 
#' The tests check the computed moments (mean, standard deviation, and
#' at-site skew), quantiles, confidence limits, variances, and k-values.  
#' 
#'
#'
checkFreqCurve <- function(PSFfile, EXPfile, excludeSites = c(""), quietly=TRUE){
  
  r <- suppressWarnings(peakfq(PSFfile, quietly = quietly)) 
  exp <- readEXP(EXPfile)
  
  r_lp3 <- r[[1]]
  r_qnt <- r[[2]]
  r_emp <- r[[3]]
  r_mgb <- r[[4]]
  
  exp_lp3 <- exp[[1]]
  exp_qnt <- exp[[2]]
  exp_lp3$site_no <- exp_lp3$Station
  
  for(n in 1:length(exp_qnt)){
    exp_qnt[[n]]$site_no <- names(exp_qnt)[n]
  }
  
  exp_qnt <- dplyr::bind_rows(exp_qnt)
  
  rownames(exp_qnt) <- 1:nrow(exp_qnt)
  
  exp_qnt$K_Value <- exp_qnt$`K-Value` #Copy column with desired name
  
  
  
  r_lp3 <- r_lp3[r_lp3$site_no %!in% excludeSites,]
  exp_lp3 <- exp_lp3[exp_lp3$site_no %!in% excludeSites,]
  
  r_qnt <- r_qnt[r_qnt$site_no %!in% excludeSites,]
  exp_qnt <- exp_qnt[exp_qnt$site_no %!in% excludeSites,]
  
  testthat::test_that("PILFs", {
    
    PILFTol <- 0
    
    testCols <- c("site_no", "Mean", "StandDev", "AtSiteSkew", "HistPeaks",  "PILFs", "PILF_0s" )
    
    
    qcheck <- merge(r_lp3[,testCols],
                    exp_lp3[,testCols], by = c("site_no"))
    qcheck$diff <- qcheck$PILFs.x - qcheck$PILFs.y
    
    failedRows <- qcheck[abs(qcheck$diff) > PILFTol,]
    
    if(nrow(failedRows) > 0){
      print(failedRows)
    }
    
    
    testthat::expect_lte(max(abs(qcheck$diff)), 0)
    
  })
  
  
  testthat::test_that("Skew", {
    
    skewTol <- 0.005
    
    testCols <- c("site_no", "Mean", "StandDev", "AtSiteSkew", "HistPeaks",  "PILFs", "PILF_0s" )
    
    # testthat::expect_equal(r_lp3[r_lp3$site_no %!in% excludeSites, testCols],
    #                        exp_lp3[exp_lp3$site_no %!in% excludeSites, testCols],
    #                        tolerance=0.01)
    
    qcheck <- merge(r_lp3[,testCols],
                    exp_lp3[,testCols], by = c("site_no"))
    qcheck$diff <- qcheck$AtSiteSkew.x - qcheck$AtSiteSkew.y
    
    failedRows <- qcheck[abs(qcheck$diff) >= skewTol,]
    
    if(nrow(failedRows) > 0){
      print(failedRows)
    }
    
    
    testthat::expect_lt(max(abs(qcheck$diff)), skewTol)
    
  })
  
  testthat::test_that("Quantiles", {
    
    qntTol <- 0.005
    
    testCols <- c("site_no", "EXC_Prob", "Estimate")
    
    qcheck <- merge(r_qnt[,testCols],
                    exp_qnt[,testCols], by = c("site_no", "EXC_Prob"))
    qcheck$diff <- qcheck$Estimate.x - qcheck$Estimate.y
    qcheck$reldiff <- qcheck$diff/qcheck$Estimate.y
    
    failedRows <- qcheck[abs(qcheck$reldiff) >= qntTol,]

    if(nrow(failedRows) > 0){
      print(failedRows)
    }
    
    testthat::expect_lt(max(abs(qcheck$reldiff)), qntTol)
    
  })
  
  testthat::test_that("ConfidenceLevels", {
    
    confTol <- 0.01
    
    testCols <- c("site_no", "EXC_Prob", "Conf_Low", "Conf_Up")
    
    qcheck <- merge(r_qnt[,testCols],
                    exp_qnt[,testCols], by = c("site_no", "EXC_Prob"))
    qcheck$low_diff <- qcheck$Conf_Low.x - qcheck$Conf_Low.y
    qcheck$up_diff <- qcheck$Conf_Up.x - qcheck$Conf_Up.y
    
    qcheck$low_reldiff <- qcheck$low_diff/qcheck$Conf_Low.y
    qcheck$up_reldiff <- qcheck$low_diff/qcheck$Conf_Up.y
    
    failedRows <- qcheck[abs(qcheck$low_reldiff) >= confTol | abs(qcheck$up_reldiff) >= confTol,]
    
    if(nrow(failedRows) > 0){
      print(failedRows)
    }
    
    testthat::expect_lt(max(abs(qcheck$low_reldiff)), confTol)
    testthat::expect_lt(max(abs(qcheck$up_reldiff)), confTol)
    
  })
  
  testthat::test_that("Variance", {
    
    varTol <- 0.0003
    
    testCols <- c("site_no", "EXC_Prob", "Variance")
    
    qcheck <- merge(r_qnt[,testCols],
                    exp_qnt[,testCols], by = c("site_no", "EXC_Prob"))
    qcheck$diff <- qcheck$Variance.x - qcheck$Variance.y
    
    failedRows <- qcheck[abs(qcheck$diff) >= varTol,]

    if(nrow(failedRows) > 0){
      print(failedRows)
    }
    
    #Don't check relative difference for variance, as can be very large due to rounding
    testthat::expect_lt(max(abs(qcheck$diff)), varTol)
    
  })
  
  testthat::test_that("K Values", {
    
    kTol <- 0.005
    
    testCols <- c("site_no", "EXC_Prob", "K_Value")
    
    qcheck <- merge(r_qnt[,testCols],
                    exp_qnt[,testCols], by = c("site_no", "EXC_Prob"))
    qcheck$diff <- qcheck$K_Value.x - qcheck$K_Value.y
    
    failedRows <- qcheck[abs(qcheck$diff) >= kTol,]
    
    if(nrow(failedRows) > 0){
      print(failedRows)
    }
    
    testthat::expect_lt(max(abs(qcheck$diff)), kTol)
    
  })
  
}

#Included only for internal use with tests
#' Read EXP File
#'
#' Function to read .EXP files from PeakFQ
#'
#' @param EXPfile character - path to .EXP file.
#'
#' @return A list containing:
#'  1. A dataframe containing computed sample moments and trend test results for each site in the .EXP file.
#'  2. Another list containing dataframes of peak-flow frequency analysis results and confidence intervals for each site in the .EXP file.
#'
#' @keywords internal
#'
#' @importFrom tidyr pivot_wider
#' @importFrom utils read.delim
readEXP <- function(EXPfile){
  #Function to read .EXP files from PeakFQ
  
  #EXPfile = path to .EXP file
  #Returns
  # 1. A dataframe containing data from the first portion of the .EXP file for each site in the .EXP file
  # 2. A list containing dataframes of EMA results and confidence intervals for each site in the .EXP file
  
  #Check if input file exists
  if (!file.exists(EXPfile)) {
    stop("EXP file doesn't exist")
  }
  
  
  #Count number of lines between site header and frequency curve table
  EXP <- file(EXPfile, "r") #Open connection to read file
  on.exit(close(EXP)) #Close .EXP connection on exit
  
  lines <- scan(EXP, what="character", sep="\n") #Get vector of EXP file lines
  infoLines <- which(grepl("EXC_Prob", lines))[1] - 2
  close(EXP)
  
  
  EXPfile <- file(EXPfile, "r") #Open connection to read file
  on.exit(close(EXPfile)) #Close .EXP connection on exit
  #Read the first site
  siteInfo <- readLines(EXPfile, n=infoLines + 1)
  info <- utils::read.delim(text=siteInfo, skip=1, sep="\t", header=FALSE) #Get information section of the site data
  infoCols <- t(gsub(" ", "", as.character(info[,1]), fixed=TRUE)) #Get info column names
  infoVals <- t(trimws(as.character(info[,2]))) #Get info values
  
  infoDF <- data.frame(matrix(nrow=1, ncol=length(infoCols))) #Create a dataframe to hold
  infoDF[1,] <- infoVals
  colnames(infoDF) <- trimws(infoCols) #Name info columns
  
  siteNum <- infoVals[1] #Get site number
  siteNums <- c(infoVals[1]) #Start vector of site numbers
  
  #Read the EMA results
  siteResults <- utils::read.delim(text=readLines(EXPfile, n=6), sep="\t", header=FALSE) #Read EMA results
  resultLabels <- as.character(siteResults[,1]) #Convert from factor to character
  siteResults <- data.frame(t(siteResults[,-1]))
  colnames(siteResults) <- trimws(resultLabels) #Label columns
  
  resultStack <- list(siteResults) #Add dataframe to list of results
  names(resultStack) <- siteNums #Label the result list with site numbers
  
  
  
  reading <- TRUE #Create a variable for whether the file is being read
  
  while (reading == TRUE){
    
    nextLine <- readLines(EXPfile, n=1) #Read next line of file
    
    
    if(length(nextLine) > 0){
      
      if(strsplit(nextLine, " - ")[[1]][1] == "Station"){
        
        siteInfo <- readLines(EXPfile, n=infoLines)
        info <- utils::read.delim(text=siteInfo, skip=0, sep="\t", header=FALSE) #Get information section of the site data
        
        infoSiteVals <- t(trimws(as.character(info[,2]))) #Get info values
        infoDF[nrow(infoDF) + 1,] <- infoSiteVals
        
        siteNums <- append(siteNums, infoSiteVals[1]) #Add site number to list of site numbers
        
        #Read the EMA results
        siteResults <- utils::read.delim(text=readLines(EXPfile, n=6), sep="\t", header=FALSE) #Read EMA results
        resultLabels <- as.character(siteResults[,1]) #Convert from factor to character
        siteResults <- data.frame(t(siteResults[,-1]))
        colnames(siteResults) <- trimws(resultLabels) #Label columns
        
        resultStack <- append(resultStack, list(siteResults)) #Add site to list of EMA results
        names(resultStack) <- siteNums #Label the result list with site numbers
        
      }
      
      else{
        warning("Unexpected line - EXP file may not have been read correctly")
        reading <- FALSE
      }
    }
    
    else{
      #print("Reached end of .EXP file")
      reading <- FALSE
    }
    
    
  }
  
  #Check if first row of column can be converted to numeric, if so convert all rows
  numCols <- !is.na(sapply(infoDF, as.numeric))[1,]
  numCols[1] <- FALSE #Don't convert first column with station number
  infoDF[,numCols] <- sapply(infoDF[,numCols], as.numeric) #Convert numeric columns to numeric type
  
  
  
  return(list(infoDF, resultStack))
  
}