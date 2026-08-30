#' Perform EMA fit
#'
#' Function to call FORTRAN EMA fitting routine
#'
#' @param QT dataframe - dataframe containing flow intervals, perception thresholds, and identifier for historic peaks in the columns:
#'   1. ql numeric - lower flow interval values
#'   2. qu numeric - upper flow interval values
#'   3. tl numeric - lower perception threshold values
#'   4. tu numeric - upper perception threshold values
#'   5. dtype integer - indicator if value is historic peak. 1 = historic peak, 0 = not historic peak. 
#'   6. peak_WY (optional) integer - water year
#'   
#'   
#' @param LOthresh numeric - low outlier test threshold:
#' - Values >= 1e-6 set manual low outlier value 
#' - Values < 1e-6, > 1e-99 disable low outlier test
#' - Values <= 1e-99 run MGBT
#' @param rG numeric - regional skew value
#' @param rGmse numeric - Mean square error (MSE) for regional skew. Set to -1e99 to use only at-site skew. 
#' @param eps numeric - confidence interval coverage. Set to 0.90 for 90% confidence intervals with lower 5% and upper 95 confidence limits. 
#' @param weightOpt character - skew weighting option. Valid options are: 
#' 1. INV - inverse mean-squared error weighting. This is the algorithm used in 
#' PeakFQ 7.4.1 and earlier
#' 2. ERL - effective record length weighting. This algorithm is used in HEC-SSP 2.3
#' 3. HWN - a generalization of the PeakFQ 7.4.1 algorithm using an optimized
#' adjustment factor when censored data are present. Results are identical to
#' PeakFQ 7.4.1 when no censored data are present.
#' 
#' @param AEPs numeric - vector of annual exceedance probabilities (AEPs) for which values are to be estimated
#' @param site_no character - site number, will be printed if `quietly` is set to `FALSE`. Not used otherwise. 
#' @param quietly logical - if `TRUE`, messages from the function are not printed
#' 
#' @returns A list containing:
#'  1. Moments and Summary Statistics - dataframe containing fitted moments and summary statistics for the analysis
#'  
#'  2. Frequency curve - dataframe containing the fitted frequency curve with estimates, confidence intervals, and variance of the estimate 
#'  for each AEP in the analysis
#'  
#'  3. EMA representation of the data
#'  
#'  4. Multiple Grubbs-Beck Test p-values for the analysis
#'  
#' @examples
#' 
#' 
#' #Initialize all flow intervals to discrete value
#' qu <- ql <- c(0.9, 1.4, 2.5, 3.0, 1.2, 1.4, 2.1, 2.6, 3.2, 2.7)
#'  
#' ql[1] <- 1e-20 #Censor first data value
#' # Set up data frame with flow intervals and perception thresholds
#' # Use same lower perception threshold for all observations
#' QT <- data.frame(ql = ql, 
#'                  qu = qu, 
#'                  tl = rep(0.9, length(qu)), 
#'                  tu = rep(1e20, length(qu)),
#'                  dtype = rep(0, length(qu)) #No historic peaks
#'                   )
#' 
#' r <- emafit(QT)
#'
#' @useDynLib peakfq, .registration = TRUE
#'
#' @export
#'
emafit <- function(QT, LOthresh = 0, rG = 0, rGmse = -1e99, eps = 0.90, weightOpt = "HWN", AEPs = c(0.9950, 0.9900, 0.9800, 0.9750, 0.9600, 0.9500, 0.9000, 0.8000, 0.7000, 0.6667, 0.6000, 0.5704, 0.5000, 0.4292, 0.4000, 0.3000, 0.2000, 0.1000, 0.0500, 0.0400, 0.0250, 0.0200, 0.0100, 0.0050, 0.0020), site_no = "", quietly = TRUE){
  
  if(!quietly & site_no != "") print(site_no) #Print site number for debugging if specified 
  
  if(!all(c("ql", "qu", "tl", "tu", "dtype") %in% colnames(QT))){
    stop("QT dataframe is missing on of the required columns")
  }
  
  if(nrow(QT) < 3){
    #Three is absolute minimum number of data points to compute skew
    stop("QT dataframe must have at least 3 rows of data")
  }
  
  if(any(QT[, c("ql", "qu", "tl", "tu")] <= 0)){
    stop("All flow intervals and perception thresholds in QT dataframe must be positive")
  }
  
  #Log-transform and convert data for FORTRAN
  QT_log <- data.frame(ql = as.double(log10(QT$ql)), 
                       qu = as.double(log10(QT$qu)), 
                       tl = as.double(log10(QT$tl)),
                       tu = as.double(log10(QT$tu)), 
                       dtype = as.integer(QT$dtype))
  
  if(any(is.na(QT_log))){
    stop("Invalid data in QT dataframe")
  }
  
  if(LOthresh >= 1e-6){
    
    
    if(any(QT$tu < LOthresh)){
      stop(paste0("One or more upper perception thresholds is less than the fixed low outlier value for site ", 
                  site_no, "."))
    }
    
    PILF_Method <- "FIXED"
    LOthresh <- log10(LOthresh)
  }
  else if(LOthresh > 1e-99){
    #For very small thresholds, set manual PILF to minimum discharge,
    #effectively disabling low outlier test
    PILF_Method <- "NONE"
    LOthresh <- log10(Qmin)
  }
  else{
    #Check that the smallest upper perception threshold is greater than the median peak to avoid numeric issues
    exactPeaks <- QT$ql[QT$ql == QT$qu & QT$dtype == 0] #Select 'exactly' known systematic peaks
    mPeak <- stats::median(exactPeaks)
    
    if(any(QT$tu < mPeak)){
      stop(paste0("One or more upper perception thresholds is less than the median (", mPeak, ") of peak values for site ", 
                  site_no, ". Upper perception thresholds less than the median are not allowed when using the Multiple Grubbs Beck Test as they can lead to numerical issues with the Expected Moments Algorithm if a PILF is detected greater than the upper perception threshold."))
    }
    
    #Check data requirement for number of nonzero peaks
    if(length(exactPeaks[exactPeaks > Qmin]) <= 5){
      stop(paste0("Need more than 5 non-zero peaks for multiple Grubbs-Beck test", 
                  site_no, "."))
    }
    
    LOthresh <- -99
    PILF_Method <- "MGBT"
  }
  
  
  #This should be modified in the FORTRAN source code to make it more sane
  if(rGmse <= 0 & rGmse >= -98){
    SkewOption = "Generalized"
  }
  else if(rGmse > 0){
    SkewOption = "Weighted"
  }
  else {
    SkewOption = "Station"
  }
  
  
  if(weightOpt %!in%  c("INV", "ERL", "HWN")){
    stop(paste("Invalid skew weighting option: ", weightOpt, 
               "Valid options are INV, ERL, or HWN"))
    
  }
  if(weightOpt == "INV"){
    weightOptInt = 3
  }
  else if(weightOpt == "ERL"){
    weightOptInt = 2
  }
  else{
    weightOptInt = 1 #Set skew weighting option to Halloween method as default
  }
  
  if(eps <= 0 | eps >= 1){
    stop(paste("Invalid confidence interval option: ", eps, 
               "Value must be greater than 0 and less than 1."))
  }
  
  startT <- Sys.time()
  
  
  n <- nrow(QT)
  
  cmoms <- matrix(rep(as.double(-99), 9), nrow=3, ncol=3) #matrix of central moments, yes it's a matrix not vector, tread carefully
  
  pq <- 1-AEPs
  nquantile <- length(AEPs) #Number of quantiles computed
  outVec <- as.double(rep(-99, nquantile))
  
  
  EMAout <- .Fortran("emafitpr", n,
                     QT_log$ql,
                     QT_log$qu,
                     QT_log$tl,
                     QT_log$tu,
                     QT_log$dtype,
                     as.double(0), #reg mean
                     as.double(-1e99), #reg mean MSE
                     as.double(1), #reg SD
                     as.double(1e99), #reg SD MSE
                     as.double(rG), #regional skew
                     as.double(rGmse), #regional skew MSE
                     as.double(LOthresh), #Fixed LO value
                     as.double(pq), #vector of quantiles to be estimated
                     as.integer(nquantile), #number of quantiles to be estimated
                     as.double(eps), #confidence interval coverage
                     as.integer(weightOptInt), #Skew weighting algorithm option
                     as.double(0), #gbcrit
                     as.integer(0), #number of flows used in MGBT (called "systematic" flows in FORTRAN, but definition is complicated)
                     as.integer(0), #number of zero  PILFs
                     as.integer(0), #number of PILFs
                     as.double(rep(-99, 20000)), #MGBT PILF p values. Length of 20000 is hard-coded into FORTRAN
                     as.double(rep(-99, 20000)), #peaks used in MGBT. Length of 20000 is hard-coded into FORTRAN
                     as.double(0), #as_G_mse
                     as.double(0), #as_G_mse_Syst
                     as.double(0), #as_G_ERL
                     cmoms,  #Estimated central moments
                     outVec, outVec, outVec, outVec, #vectors of estimates, CIs, and variances
                     as.double(0), #Wd
                     as.double(rep(-99, 25000)),
                     as.double(rep(-99, 25000)), 
                     as.double(rep(-99, 25000)),
                     as.double(rep(-99, 25000)), 
                     outVec
  )
  
  endT <- Sys.time()
  
#  print("FORTRAN DONE")
#  print(EMAout[[37]])
  
  
  cmoms <- EMAout[[27]]
  #cmoms is a matrix containing the first 3 central moments
  #estimated using:
  #regional skew with final MSE formula (column 1)
  #at-site skew (column 2)
  #regional skew with B17B MSE formula (column 3)
  MGBTcrit <- EMAout[[18]]
  PILFthresh <- 10^MGBTcrit
  nPILFs <- EMAout[[21]] 
  
  
  
  EXPout <- data.frame(
    BegYear = NA,
    EndYear = NA,
    RecordLength = n, 
    SkewOption = SkewOption,
    Mean = cmoms[1,1],
    StandDev = sqrt(cmoms[2,1]),
    Skew = cmoms[3,1],
    AtSiteMean = cmoms[1,2],
    AtSiteStandDev = sqrt(cmoms[2,2]),
    AtSiteSkew = cmoms[3,2],
    AtSiteMSEG = EMAout[[24]],
    AtSiteMSEG_GagedOnly = EMAout[[25]],
    #PRL = EMAout[[26]],
    RegSkew = rG,
    RegMSEG = rGmse,
    HistPeaks = length(QT$dtype[QT$dtype == 1]), #Count historic peaks
    MGBTPeaks = EMAout[[19]],
    PILF_Method = PILF_Method,
    PILF_Thresh = PILFthresh,
    PILFs = EMAout[[21]],
    PILF_0s = EMAout[[20]], 
    WeightOpt = weightOpt,
    WeightCo = EMAout[[32]], 
    RunTime = as.numeric(endT - startT))
  
  if("peak_WY" %in% colnames(QT)){
    EXPout$BegYear <- min(QT$peak_WY)
    EXPout$EndYear <- max(QT$peak_WY)
  }
  
  Qout <- data.frame(EXC_Prob = AEPs,
                     Estimate = 10^EMAout[[28]],
                     Variance = EMAout[[31]],
                     Conf_Low = 10^EMAout[[29]],
                     Conf_Up = 10^EMAout[[30]]#, 
                     #nu = EMAout[[37]] #Don't publish nu values to output; nu may not be well-defined due to CI weighting
                     )
  
  dataout <- data.frame(ql_ema = 10^EMAout[[33]][1:n], 
                        qu_ema = 10^EMAout[[34]][1:n], 
                        tl_ema = 10^EMAout[[35]][1:n], 
                        tu_ema = 10^EMAout[[36]][1:n]
  )
  
  
  nPILFs <- EMAout[[21]]
  
  if(nPILFs > 0){
    
    #Need to do this by indices due to complications with interval peaks 
    #which may be less than PILF threshold.
    PILFinds <- which(EMAout[[23]][1:nPILFs] < PILFthresh) #Get indices of PILFs
    
    
    MGBTout <- data.frame(
      peak_va = EMAout[[23]][PILFinds], #Vector of MGBT PILF values
      p_value = EMAout[[22]][PILFinds] #Vector of MGBT PILF values
    )
    
    #For fixed PILF, all pvals are zero and shouldn't be output
    if(PILF_Method == "FIXED"){
      MGBTout$p_value <- as.numeric(NA)
    }
    
  }
  else{
    MGBTout <- data.frame(
      peak_va = as.numeric(NA), #Vector of MGBT PILF values
      p_value = as.numeric(NA) #Vector of MGBT PILF values
    )
  }
  
  
  
  
  
  return(list(EXPout, Qout, dataout, MGBTout))
  
  
  
}

#' Hirsch-Stedinger Plotting Positions
#'
#' Function to compute Hirsch-Stedinger Plotting Positions
#'
#' @param QT dataframe - dataframe of flow intervals and perception thresholds. 
#'           Must have the following columns:
#'           - ql - lower flow interval values
#'           - qu - upper flow interval values
#'           - tl - lower perception thresholds
#'           - tu - upper perception thresholds
#'           - peak_WY - water year
#'           - dtype (optional) - indicator if value is historic peak. 1 = historic peak, 0 = not historic peak. 
#' @param ppalpha numeric - alpha value for plotting positions. Set to zero for Weibull. 
#'
#' @return dataframe with empirical plotting positions
#' 
#' @examples
#' #Initialize all flow intervals to discrete value
#' qu <- ql <- c(0.9, 1.4, 2.5, 3.0, 1.2, 1.4, 2.1, 2.6, 3.2, 2.7)
#'  
#' ql[1] <- 1e-20 #Censor first data value
#' # Set up data frame with flow intervals and perception thresholds
#' # Use same lower perception threshold for all observations
#' QT <- data.frame(ql = ql, 
#'                  qu = qu, 
#'                  tl = rep(0.9, length(qu)), 
#'                  tu = rep(1e20, length(qu)),
#'                  peak_WY = c(1:10)
#'                  )
#' 
#' pp <- plotPosHS(QT, ppalpha = 0)
#' 
#' @useDynLib peakfq, .registration = TRUE
#' 
#' @export
#'
plotPosHS <- function(QT, ppalpha = 0){
  
  if(!all(c("ql", "qu", "tl", "tu", "peak_WY") %in% colnames(QT))){
    stop("QT dataframe is missing on of the required columns")
  }
  
  if(any(QT[, c("ql", "qu", "tl", "tu")] <= 0)){
    stop("All flow intervals and perception thresholds in QT dataframe must be positive")
  }
  
  n <- nrow(QT)
  
  #Address greater than peaks in same manner as PeakFQ 7.4
  #If upper interval is really large/infinity, set to lower interval 
  QT$qu_adj <- QT$qu
  QT$qu_adj[QT$qu_adj >= Qmax] <- QT$ql[QT$qu_adj >= Qmax]
  
  
  
  
  
  
  #Log-transform and convert data for FORTRAN
  QT_log <- data.frame(ql = as.double(log10(QT$ql)), 
                       qu = as.double(log10(QT$qu_adj)), 
                       tl = as.double(log10(QT$tl)),
                       tu = as.double(log10(QT$tu))
                       )
  
  if(any(is.na(QT_log))){
    stop("Invalid data in QT dataframe")
  }
  
  #Lengths of arrays in FORTRAN are a mess. Need to fix there later. 
  outVec <- as.double(rep(-99, n))
  outVecLong <- as.double(rep(-99, 20000))
  
  
  
  PPout <- .Fortran("PLOTPOSHS",
                    n, #Number of observations, censored or uncensored
                    QT_log$ql, #Vector of lower interval values
                    QT_log$qu, #Vector of adjusted upper interval values
                    QT_log$tl, #Vector of lower threshold values
                    QT_log$tu, #Vector of upper threshold values
                    as.double(ppalpha), #Alpha value for plotting position. Zero for Weibull
                    outVecLong, #Sorted discharges
                    outVec, #Estimated exceedance probabilities
                    integer(0), #Number of censoring thresholds
                    outVecLong, #Vector of censoring thresholds, need to trim end of vector
                    outVec, #Estimated exceedance probabilities of thresholds, need to trim end of vector
                    as.integer(rep(0, n))   #Number of observations censored at each threshold
  )
  
  
  EMP <- data.frame(peak_WY = QT$peak_WY,
                    peak_va = 10^PPout[[7]][1:n],
                    plot_pos = PPout[[8]], 
                    ql = QT$ql, 
                    qu = QT$qu, 
                    tl = QT$tl, 
                    tu = QT$tu, 
                    historic = NA)
  
  if("dtype" %in% colnames(QT)){
    EMP$historic <- as.logical(QT$dtype)
  }
  
  #Set peak values with an upper interval value of infinity to NA
  EMP$peak_va[EMP$peak_va >= Qmax/2] <- NA 
  
  
  return(EMP)
  
}

#' Expected Value of Central Moments Derivative
#'
#' Wrapper for EXPMOMCDERIV subroutine
#'
#' @param m numeric - mean
#' @param sd numeric - standard deviation
#' @param g numeric - skew
#' @param tl numeric - lower censoring threshold
#' @param tu numeric - upper censoring threshold
#' 
#' @useDynLib peakfq, .registration = TRUE
#'
#' @keywords internal
#'
EXPMOMCDERIV <- function(m, sd, g, tl, tu){
  
  #Compute lp3 parameters
  alpha <- 4/g^2
  beta <- g/abs(g)*(sd^2/alpha)^0.5
  tau <- m - alpha*beta

  r <- .Fortran("EXPMOMCDERIV", as.double(c(tau, alpha, beta)),
                     as.double(tl),
                     as.double(tu),
                     #Ensure proper dimensions for outputs
                     as.double(rep(-99, 3)), #Noncentral moments
                     matrix(as.double(rep(-99, 9)), nrow=3, ncol=3)
  )
  
  return(r[4:5])
  
}

#' Expected Value of Noncentral Moments Derivative
#'
#' Wrapper for DEXPECT subroutine
#'
#' @param m numeric - mean
#' @param sd numeric - standard deviation
#' @param g numeric - skew
#' @param tl numeric - lower censoring threshold
#' @param tu numeric - upper censoring threshold
#' 
#' @useDynLib peakfq, .registration = TRUE
#'
#' @keywords internal
#'
DEXPECT <- function(m, sd, g, tl, tu){
  
  #Compute lp3 parameters
  alpha <- 4/g^2
  beta <- g/abs(g)*(sd^2/alpha)^0.5
  tau <- m - alpha*beta
  
  r <- .Fortran("DEXPECT", as.double(tau),
                     as.double(alpha), 
                     as.double(beta),
                     as.double(tl),
                     as.double(tu),
                     #Ensure proper dimensions for outputs
                     as.double(rep(-99, 3)), #Noncentral moments
                     matrix(as.double(rep(-99, 9)), nrow=3, ncol=3)
  )
  
  return(r[6:7])
  
}


#' Compute Skew Weighting Determinate Ratio
#'
#' Computes the determinate ratio used to correct for censored data when weighting
#' at-site and regional skew coefficients. Is a wrapper for the DEXPECT subroutine. 
#'
#' @param m numeric - mean
#' @param sd numeric - standard deviation
#' @param g numeric - skew
#' @param n integer - number of observations
#' @param nthresh integer - number of censoring thresholds
#' @param nobs integer vector - number of observations in each censoring threshold
#' @param tl numeric vector - lower censoring thresholds
#' @param tu numeric vector - upper censoring thresholds
#' 
#' @return Determinate ratio used as a correction factor when computing weighted skew coefficient
#' 
#' 
#' 
#' @useDynLib peakfq, .registration = TRUE
#'
#' @keywords internal
#'
detrat <- function(m, sd, g, n, nthresh, nobs, tl, tu){
  
  #Call EXPMOMCDERIV subroutine
  r <- .Fortran("detratsub", as.double(c(m, sd^2, g)),
                     as.integer(n), 
                     as.integer(nthresh), 
                     as.double(nobs),
                     as.double(tl),
                     as.double(tu),
                     #Ensure proper dimensions for outputs
                     as.double(-99) 
  )
  
  return(r[[7]])
  
}


#' Pearson Type III Expected Moments
#'
#' Computes the expected value of a Pearson Type III variate with parameters given by mc_old(), with observed values known to lie in the interval (ql,qu)
#'
#' @param ql numeric - vector of lower bounds on (log) floods
#' @param qu numeric - vector of upper bounds on (log) floods
#' @param rG numeric - regional skew value (needed for B17C EQ 7-10)
#' @param nG numeric - equivalent years of record for regional skew
#' @param mc_old numeric vector - p3 parameters
#' 
#' @return First 3 central moments (mean, variance, skew)
#' 
#' @useDynLib peakfq, .registration = TRUE
#'
#' @keywords internal
#'
moms_p3_R <- function(ql, qu, rG, nG, mc_old){
  
  ql <- as.double(ql)
  qu <- as.double(qu)
  
  if(any(is.na(ql)) | any(is.na(qu))){
    stop("Invalid value in ql or qu")
  }
  
  if(length(ql) != length(qu)){
    stop("ql and qu must be the same length")
  }
  
  n <- length(ql)
  
  #Call FORTRAN subroutine
  r <- .Fortran("moms_p3",
                n,
                ql, 
                qu, 
                as.double(rG),
                as.double(nG),
                as.double(mc_old),
                #Ensure proper dimensions for outputs
                as.double(c(-99, -99, -99)) 
  )
  
  return(r[[7]])
  
}





#' Pearson Type III EMA Moments
#'
#' Fits Pearson Type III distribution to censored data using expected moments algorithm and method of moments estimators
#'
#' @param ql numeric - vector of lower bounds on (log) floods
#' @param qu numeric - vector of upper bounds on (log) floods
#' @param as_M_mse numeric - mean square error of at-site mean (M)
#' @param as_S2_mse numeric - mean square error of at-site variance (S2)
#' @param as_G_mse numeric - mean square error of at-site skew (G)
#' @param r_M numeric - regional mean (M)
#' @param r_S2 numeric - regional variance (std.dev^2) (S2)
#' @param r_G numeric - regional skew (G)
#' @param r_M_mse numeric - mean square error of regional mean; set to -99 to turn off weighting the mean
#' @param r_S2_mse numeric - mean square error of regional variance; set to -99 to turn off weighting the variance
#' @param r_G_mse numeric - mean square error of regional skew; set to -99 to turn off weighting the skew
#' @param Wd numeric weighting factor for regional skew effective record length with Halloween method
#' 
#' @return first 3 central moments (mean, variance, skew)
#' 
#' @useDynLib peakfq, .registration = TRUE
#'
#' @keywords internal
#'
p3est_ema_R <- function(ql, qu, as_M_mse, as_S2_mse, as_G_mse, r_M = -99, r_S2 = -99, r_G = -99, r_M_mse = -99, r_S2_mse = -99, r_G_mse = -99, Wd = 1){
  
  ql <- as.double(ql)
  qu <- as.double(qu)
  
  if(any(is.na(ql)) | any(is.na(qu))){
    stop("Invalid value in ql or qu")
  }
  
  if(length(ql) != length(qu)){
    stop("ql and qu must be the same length")
  }
  
  n <- length(ql)
  
  #Call FORTRAN subroutine
  r <- .Fortran("p3est_ema",
                n,
                ql, 
                qu, 
                as.double(as_M_mse),
                as.double(as_S2_mse),
                as.double(as_G_mse),
                as.double(r_M),
                as.double(r_S2),
                as.double(r_G),
                as.double(r_M_mse),
                as.double(r_S2_mse),
                as.double(r_G_mse), 
                as.double(Wd),
                #Ensure proper dimensions for outputs
                as.double(c(-99, -99, -99)) 
  )
  
  return(r[[14]])
  
}

#' Pearson Type III EMA Variance
#'
#' Computes variance of EMA estimator
#'
#' @param nthresh integer - number of distinct censoring thresholds
#' @param nobs integer vector - number of observations in each censoring threshold
#' @param tl numeric vector - lower censoring thresholds
#' @param tu numeric vector - upper censoring thresholds
#' @param mc numeric vector - estimates of first 3 central moments (mean, variance, skew)
#' 
#' @return variance-covariance matrix for the central moment estimators
#' 
#' @useDynLib peakfq, .registration = TRUE
#'
#' @keywords internal
#'
var_mom_R <- function(nthresh, nobs, tl, tu, mc){
  
  tl <- as.double(tl)
  tu <- as.double(tu)
  
  if(any(is.na(tl)) | any(is.na(tu))){
    stop("Invalid value in tl or tu")
  }
  
  
  #Call FORTRAN subroutine
  r <- .Fortran("var_mom",
                as.integer(nthresh),
                as.double(nobs), #Yes, this is always an integer but needs to be converted to a double for FORTRAN
                tl, 
                tu,
                as.double(mc),
                #Ensure proper dimensions for outputs
                matrix(rep(as.double(-99), 9), nrow = 3, ncol=3)
  )
  
  return(r[[6]])
  
}

#' Pearson Type III Quantiles
#'
#' Computes the inverse CDF of a Pearson Type III variate with moments m=(mu,sigma^2,g)
#'
#' @param q numeric - probability
#' @param m numeric - vector moments (mean, variance, skew) defining the Pearson Type III distribution
#' 
#' @return quantile value
#' 
#' @useDynLib peakfq, .registration = TRUE
#'
#' @keywords internal
#'
qP3_R <- function(q, m){
  
  q <- as.double(q)
  m <- as.double(m)
  
  if(is.na(q) | q < 0 | q > 1){
    stop("Invalid q")
  }
  
  if(any(is.na(m))){
    stop("Invalid m")
  }
  
  #Call FORTRAN subroutine
  r <- .Fortran("qP3sub",
                q,
                m, 
                #Ensure proper dimensions for outputs
                as.double(-99) 
  )
  
  return(r[[3]])
  
}


#' Skew Mean-Square Error
#'
#' Computes the mean-square error of the at-site skew coefficient
#'
#' @param nthresh integer - number of distinct censoring thresholds
#' @param nobs integer vector - number of observations in each censoring threshold
#' @param tl numeric vector - lower censoring thresholds
#' @param tu numeric vector - upper censoring thresholds
#' @param mc numeric vector - estimates of first 3 central moments (mean, variance, skew)
#' 
#' @return mean-square error of the at-site skew coefficient
#' 
#' @useDynLib peakfq, .registration = TRUE
#'
#' @keywords internal
#'
mseg_all_R <- function(nthresh, nobs, tl, tu, mc){
  
  tl <- as.double(tl)
  tu <- as.double(tu)
  
  if(any(is.na(tl)) | any(is.na(tu))){
    stop("Invalid value in tl or tu")
  }
  
  
  #Call FORTRAN subroutine
  r <- .Fortran("mseg_all_sub",
                as.integer(nthresh),
                as.double(nobs), #Yes, this is always an integer but needs to be converted to a double for FORTRAN
                tl, 
                tu,
                as.double(mc),
                #Ensure proper dimensions for outputs
                as.double(-99)
  )
  
  return(r[[6]])
  
}
