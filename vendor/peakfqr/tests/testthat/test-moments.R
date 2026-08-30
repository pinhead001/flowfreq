testthat::context("Included Moment Tests")
#Make temporary directory for test output files

inPSF <- system.file("testdata",
                     "wymt_ffa_2022A.psf",
                     package="peakfq")

inData <- system.file("testdata",
                     "wymt_ffa_2022A_WATSTORE.TXT",
                     package="peakfq")

r <- peakfq::peakfq(inPSF)

d <- r[[1]]
q <- r[[2]] #quantile output
e <- r[[3]] #empirical plotting positions
mgb <- r[[4]] #MGBT results

d <- d[order(d$site_no),]
q <- q[order(q$site_no),]

test_that("Site Info", {
  
  
  
  v7_4_output <- read.csv(system.file("testdata",
                         "wymt_ffa_2022A_EXPinfo_7_4.csv",
                         package="peakfq"), header = TRUE, 
                         colClasses=c("site_no"="character"))
  
  v7_4_output$site_no <- as.character(v7_4_output$site_no)
  v7_4_output <- v7_4_output[order(v7_4_output$site_no),]
  
  

  
  testCols <- c("site_no", "station_nm", "BegYear", "EndYear", "HistPeaks",  "SkewOption")
  
  
  expect_equal(d[,testCols], v7_4_output[, testCols])
})

test_that("Moments", {
  

  v7_4_output <- read.csv(system.file("testdata",
                                      "wymt_ffa_2022A_EXPinfo_7_4.csv",
                                      package="peakfq"), header = TRUE,
                                      colClasses=c("site_no"="character"))
  
  v7_4_output$site_no <- as.character(v7_4_output$site_no)
  v7_4_output <- v7_4_output[order(v7_4_output$site_no),]
  
  
  testCols <- c("site_no", "Mean", "StandDev", "AtSiteSkew")

  expect_equal(d[,testCols], v7_4_output[, testCols], tolerance=0.03)
})

test_that("MGBT", {
  

  
  v7_4_output <- read.csv(system.file("testdata",
                                      "wymt_ffa_2022A_EXPinfo_7_4.csv",
                                      package="peakfq"), header = TRUE,
                          colClasses=c("site_no"="character"))
  
  v7_4_output$site_no <- as.character(v7_4_output$site_no)
  v7_4_output <- v7_4_output[order(v7_4_output$site_no),]
  
  
  testCols <- c("site_no", "PILF_Method", "PILF_Thresh", "PILFs", "PILF_0s")
  
  expect_equal(d[,testCols], v7_4_output[, testCols])
})

test_that("MGBT_pvals", {
  
  
  
  v7_4_output <- read.csv(system.file("testdata",
                                      "wymt_ffa_2022A_MGBT_7_5_1.csv",
                                      package="peakfq"), header = TRUE,
                          colClasses=c("site_no"="character"))
  
  v7_4_output$site_no <- as.character(v7_4_output$site_no)
  v7_4_output <- v7_4_output[order(v7_4_output$site_no),]
  
  #Need to exclude data from 06326960 because PeakFQ 7.X treats 2005
  #peak of 0.5 cfs with code 4 like an exactly known peak even though 
  #it's shown as an interval in the .PRT file
  mgb <- mgb[mgb$site_no != "06326960.00",]
  v7_4_output <- v7_4_output[v7_4_output$site_no != "06326960.00",]
  
  mgb <- mgb[!is.na(mgb$peak_va),] #Remove NA's from sites without PILFs
  mgb <- mgb[mgb$peak_va > 0,] #Remove zero peaks not included in PeakFQ 7.X
  
  #Renumber rows for test
  rownames(mgb) <- 1:nrow(mgb) 
  rownames(v7_4_output) <- 1:nrow(v7_4_output) 
  
  
  testCols <- c("site_no", "peak_va", "p_value")
  
  expect_equal(mgb[,testCols], v7_4_output[, testCols])
})

test_that("Quantiles", {
  
  
  v7_4_output <- read.csv(system.file("testdata",
                                      "wymt_ffa_2022A_EXPdata_7_4.csv",
                                      package="peakfq"), header = TRUE,
                          colClasses=c("site_no"="character"))
  
  v7_4_output$site_no <- as.character(v7_4_output$site_no)
  v7_4_output <- v7_4_output[order(v7_4_output$site_no),]
  
  #Only test station skew sites due to change in skew algorithm
  v7_4_output <- v7_4_output[v7_4_output$site_no %in% d$site_no[d$SkewOption == "Station"],]
  q <- q[q$site_no %in% d$site_no[d$SkewOption == "Station"],]
  
  
  testCols <- c("site_no", "EXC_Prob",	"Estimate",	"Variance",	"Conf_Low",	"Conf_Up")
  
  expect_equal(q[,testCols], v7_4_output[, testCols], tolerance=0.001)
})

test_that("Plotting Positions", {


  v7_4_output <- read.csv(system.file("testdata",
                                      "wymt_ffa_2022A_EMPdata_7_4.csv",
                                      package="peakfq"), header = TRUE,
                          colClasses=c("site_no"="character"))

  v7_4_output$site_no <- as.character(v7_4_output$site_no)
  
  #Sort by peak value, then year
  v7_4_output <- v7_4_output[order(v7_4_output$peak_va, decreasing = T),]
  
  

  e <- e[order(e$peak_va, decreasing = T),]
  
  #e <- e[complete.cases(e$peak_va),]
  
  #For testing, select only rows with peaks and site numbers included in v 7.4 output
  e_test <- e[paste(e$site_no, e$peak_WY, sep="_") %in% paste(v7_4_output$site_no, abs(v7_4_output$peak_WY), sep="_"),]
  
  
  
  v7_4_output$peak_WY <- abs(v7_4_output$peak_WY) #Get rid of minus signs
  
  #Now sort by plotting position, then site number 
  v7_4_output <- v7_4_output[order(v7_4_output$plot_pos, decreasing = F),]
  e_test <- e_test[order(e_test$plot_pos, decreasing = F),]
  
  v7_4_output <- v7_4_output[order(v7_4_output$site_no),]
  e_test <- e_test[order(e_test$site_no),]
  
  #Remove sites affected by rounding of perception thresholds
  #This is a bug in version 7.4
  badPPSites <- c("06326960.00", "06328100.00")
  v7_4_output <- v7_4_output[v7_4_output$site_no %!in% badPPSites,]
  e_test <- e_test[e_test$site_no %!in% badPPSites,]
  

  
  #Only test water year and plot position
  testCols <- c("site_no", "peak_WY", "plot_pos")
  v7_4_output <- v7_4_output[!duplicated(v7_4_output[, testCols]), c("site_no", "peak_WY",	"peak_va", "plot_pos")]
  rownames(v7_4_output) <- 1:nrow(v7_4_output)
  rownames(e_test) <- 1:nrow(e_test)
  

  expect_equal(e_test[,testCols], v7_4_output[,testCols], tolerance=0.001)
})

