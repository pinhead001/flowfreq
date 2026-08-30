testthat::context("Included Tests for Skew Weighting")
# ############## Tests for Greg Schwarz's Halloween Skew Weighting ##############
#
# This test compares the Halloween determinate ratio weighting factor
# from peakfq to independent output from Greg Schwarz's SAS program. 
# Output from SAS is contained in the results_WHIST.csv file. 
# Greg described this output in detail in an email to Seth Siefken and 
# Peter McCarthy on November 15, 2022. 
#
# Description of the main fields in the results_WHIST.csv file are as follows: 
#
# SIGMA2 = variance of log peak annual flow
# SKEW = skew of log peak annual flow
# FREQ = frequency of censoring (for the censored share of records)
# TAU =  implied normalized censoring threshold 
# NU = share of record not subject to censoring
# PROBUT = the probability of a being an uncensored value for the share of the record (1 – NU) subjected to censoring
# HERL_factor = the Halloween ERL scaling factor; this column is computed in Excel by hand using the equations for determinant applied to the elements.
# JERL_factor = the implied ERL scaling factor for Jery’s method; this column is computed in Excel from the ACOV_SKEWT and acov_SKEW_SKEW columns.

test_that("Skew Weighting", {
  
  WHIST <- read.csv(system.file("testdata",
                                "results_WHIST.csv",
                                package="peakfq"))
  
  #Compute inputs to FORTRAN code
  
  d <- WHIST[, c("NU",	"SIGMA2",	"FREQ",	"TAU",	"SKEW", "PROBUT")]
  d$MU <- 2 #Choose mean of 2
  d$SIGMA <- sqrt(d$SIGMA2)
  
  n <- 100 #Total number of observations, doesn't matter for testing
  
  nobs_list <- mapply(c,  n*(1-d$NU), n*d$NU, SIMPLIFY = F)
  
  #Construct lists of pereception thresholds for censored and uncensored portions of the record
  tl_list <- mapply(c, d$SIGMA*d$TAU+d$MU, rep(-20, nrow(d)), SIMPLIFY = F) 
  tu_list <- mapply(c, rep(20, nrow(d)), rep(20, nrow(d)), SIMPLIFY = F) 
  
  
  r <- mapply(detrat, 
              m = d$MU, 
              sd = d$SIGMA, 
              g = d$SKEW, 
              n = rep(n, nrow(d)),
              nthresh = rep(2, nrow(d)), 
              nobs = nobs_list, 
              tl = tl_list,
              tu = tu_list)
  
  #Get Wd results from FORTRAN output
  Wd <- unlist(r)
  
  
  expect_equal(Wd, WHIST$HERL_factor, tolerance=0.001)
})