
#' Plot Flood Frequency Analysis Curve
#'
#' Function to plot results of peak-flow frequency analysis
#'
#' @param plotPositions dataframe containing observed peak values in a column named "peak_va",corresponding plotting positions in a column named "plot_pos", 
#'              and an optional column named "historic" containing a logical TRUE for peaks which are historic
#' @param LP3fit dataframe containing LP3 curve fit for a site with exceedance probabilities in a column named "EXC_Prob", and corresponding peak estimates in a column named "Estimate",
#'        lower confidence limits in a column named "Conf_Low", and upper confidence limits in a column named "Conf_Up"
#' @param eps numeric - confidence interval coverage, 0.90 indicates 90% confidence intervals with lower 5% and upper 95 confidence limits
#' @param siteNum character - site number to print on plot
#' @param siteName character - site name to print on plot
#' @param PILF numeric - PILF threshold for analysis
#' @param infoText character - text to display in the plot information text box
#' 
#' @returns Base R graphics object with a plot of the frequency curve
#' 
#' @details
#' The plot generated is optimized for viewing and saving in a window
#' 9 inches wide and 6 inches high. The plot can be resized
#' but excessive resizing may result in overlapping text or other display
#' issues.
#' 
#' @examples
#' \dontrun{
#' #Location of analysis specifications file
#' psf_path <- system.file("extdata", "ExampleSpecifications.psf", package="peakfq")
#' 
#' #Compute plotting positions and LP3fit
#' results <- peakfq(psf_path)
#' 
#' #Get plotting positions for site to plot
#' site <- "03554000"
#' pp_all <- results[[3]]
#' pp <- pp_all[pp_all$site_no == site,]
#' 
#' #Get LP3 curve for site to plot
#' lp3_all <- results[[2]]
#' lp3 <- lp3_all[lp3_all$site_no == site,]
#'  
#' #Create plot
#' p <- ffaPlot(pp, lp3, eps = unique(lp3$Conf_Interval), siteNum = site)
#'  
#' }
#'
#'
#' @export
#' 
ffaPlot <- function(plotPositions, LP3fit, eps = NA, siteNum="", siteName="", PILF = 0, infoText = ""){
  #Set logical for whether each peak is historic
  if(!is.null(plotPositions$historic)){
    plotPositions$historic[is.na(plotPositions$historic)] <- FALSE
    plotPositions$historic[plotPositions$historic != TRUE] <- FALSE
  }else{
    plotPositions$historic <- FALSE
  }
  
  
  #Apply normal transformation to plotting positions
  #Negative is needed as these are exceedance (not nonexceedance) probabilities
  plotPositions$plot_pos_qnorm <- -stats::qnorm(plotPositions$plot_pos)
  
  #Handle censored and interval peaks
  plotPositions$interval <- FALSE
  plotPositions$censored <- FALSE
  plotPositions$interval[plotPositions$ql != plotPositions$qu] <- TRUE
  #Censored peaks are intervals either a zero lower bound or infinite upper bound
  plotPositions$censored[plotPositions$interval & (plotPositions$ql <= Qmin | plotPositions$qu >= Qmax)] <- TRUE
  
  #Pull interval and censored peaks into their own dataframe
  intervalPeaks <- plotPositions[plotPositions$interval & !plotPositions$censored, c("plot_pos_qnorm", "ql", "qu")] 
  cPeaks <- plotPositions[plotPositions$censored, c("plot_pos_qnorm", "ql", "qu")] 
  
  #Need some mathematical trickery to plot the censored peaks in one series
  if(nrow(cPeaks) > 0){
    cPeaks$q_bound <- NA
    cPeaks$q_inf <- NA
    cPeaks$q_bound[cPeaks$ql <= Qmin] <- cPeaks$qu[cPeaks$ql <= Qmin]
    cPeaks$q_bound[cPeaks$qu >= Qmax] <- cPeaks$ql[cPeaks$qu >= Qmax]
    
    cPeaks$q_inf[cPeaks$ql <= Qmin] <- cPeaks$ql[cPeaks$ql <= Qmin]
    cPeaks$q_inf[cPeaks$qu >= Qmax] <- cPeaks$qu[cPeaks$qu >= Qmax]
  }
  
  
  plotPositions <- plotPositions[!plotPositions$interval, ]
  
  #Separate systematic peaks from PILFs
  #plotPositions$gagedPeaks <- NA 
  gagedPeaks <- plotPositions[plotPositions$peak_va >= PILF & plotPositions$historic != TRUE & !plotPositions$interval, c("plot_pos_qnorm", "peak_va"), drop=FALSE]
  PILFs <- plotPositions[plotPositions$peak_va < PILF, c("plot_pos_qnorm", "peak_va"), drop=FALSE]
  PILFs <- PILFs[PILFs$peak_va > Qmin,] #Remove zeros from plotting
  historicPeaks <- plotPositions[plotPositions$peak_va >= PILF & plotPositions$historic == TRUE & !plotPositions$interval, c("plot_pos_qnorm", "peak_va"), drop=FALSE]
  
  
  #Remove portion of frequency curve less than PILF threshold
  LP3fit <- LP3fit[LP3fit$Estimate >= PILF,]
  
  
  
  #Set up AEPs and AEP labels for the plot
  AEPs <- c(0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2,  0.3, 0.4, 0.5, 0.6, 0.7,  0.8, 0.9, 0.95, 0.98, 0.99, 0.995)
  AEPlabels <- c("0.2", "", "1", "", "5", "10", "20",  "", "40", "", "60", "",  "80", "90", "95", "98", "", "99.5")
  AEPs_qnorm <- -stats::qnorm(AEPs)
  
  
  #Calculate the number of axis breaks needed
  plotVals <- log10(c(plotPositions$peak_va, LP3fit$Estimate, LP3fit$Conf_Low, LP3fit$Conf_Up)) #Combine values and log transform
  plotVals[is.infinite(plotVals) | plotVals == log10(Qmin) | plotVals == log10(Qmax)] <- NA #Remove any infinite values, can occur from zero peaks

  
  maxBreak <- ceiling(max(plotVals, na.rm = TRUE)) #Get max break on the axis
  minBreak <- floor(min(plotVals, na.rm = TRUE)) #Get min break on the axis
  
  #minBreak <- max(minBreak, -1, na.rm=TRUE) #If a zero value is in the plotted data, set minBreak to -1
  
  breakPoints <- 10^seq(minBreak, maxBreak) #Set breaks at integer powers of 10
  
  #use outer product to get minor tick locations
  minorTicks <- as.vector(outer(1:9, breakPoints))
  minorTicks <- minorTicks[minorTicks < max(breakPoints)] #Remove ticks greater than greatest break point
  
  #Set number of decimal places to show on scale
  if(minBreak < 0){
    scaleAccuracy <- 0.1
  }else{
    scaleAccuracy <- 1
  }
  
  #Prepare text for confidence interval legend item
  
  
  if(!is.na(eps) & is.numeric(eps)){
    #Compute upper and lower confidence interval values in percent
    CI_low <- ((1-eps)/2)*100
    CI_up <- (1 - (1-eps)/2)*100
  }else{
    CI_low <- ""
    CI_up <- ""
  }
  
  CI_legend_text <- paste("Confidence limits:", CI_low, "percent lower,", CI_up, "percent upper")
  
  xTitle <- paste("Annual exceedance probability, in percent \n Station - ", siteNum, siteName)
  
  
  ### Prepare the legend items ###
  leg_labels <- c("Gaged peak discharge",
                  "PILF",
                  "Historic peak discharge", 
                  "Interval peak discharge", 
                  "Censored peak discharge", 
                  "Fitted frequency curve", 
                  CI_legend_text)
  leg_points <- c(16, 1, 24, NA, NA, NA, NA)
  leg_lines <- c(0, 0, 0, 1, 1, 1, 1)
  leg_bg <- c("cyan", NA, "purple", NA, NA, NA, NA)
  leg_col <- c("cyan", "black", "black", "cyan", "grey", "red", "blue")
  leg_lwd <- c(1,1.5,1,1.5,1,2,2)
  
  #Find out which potential legend items are actually used
  leg_used <- as.logical(c(nrow(gagedPeaks) > 0,
                           nrow(PILFs) > 0,
                           nrow(historicPeaks) > 0, 
                           nrow(intervalPeaks) > 0,
                           nrow(cPeaks) > 0,
                           TRUE, 
                           TRUE))
  leg_rows <- sum(leg_used, na.rm = TRUE) #Count number of items in legend
  
  #Set up to record plot without displaying
  cur_dev <- grDevices::dev.cur()   # store current device
  pdf(NULL, width = 7.5, height = 5)  # open null device
  null_dev <- grDevices::dev.cur()  # store null device
  
  # make sure we always clean up properly, even if something causes an error
  on.exit({
    grDevices::dev.off(null_dev)
    if (cur_dev > 1) grDevices::dev.set(cur_dev) # only set cur device if not null device
  })
  
  ### Do a dummy plot to get the size of the legend ###
  
  # plot(gagedPeaks$plot_pos_qnorm, gagedPeaks$peak_va, col="black", pch=1, 
  #      ylim = c(min(breakPoints), max(breakPoints)),
  #      yaxs = "i", #Needed to remove 4% inner margin
  #      yaxt='n', 
  #      xaxt='n')
  # 
  # #Add legend to dummy plot and get its height
  # legend_size <- legend("topleft",
  #                       legend=leg_labels[leg_used],
  #                       pch=leg_points[leg_used],
  #                       lty=leg_lines[leg_used],
  #                       pt.bg=leg_bg[leg_used],
  #                       col = leg_col[leg_used],
  #                       lwd = leg_lwd[leg_used],
  #                       title="Explanation",
  #                       ncol = 2)$rect$h
  # 
  # grDevices::dev.off(null_dev) #close null device
  
  # #re-adjust y axis to accommodate legend
  # maxBreak <- maxBreak + legend_size*1.16
  # tu$thresh_plot <- pmin(tu$thresh, maxBreak)
  
  ### Now do the actual plot ###
  
  
  pdf(NULL, width = 9, height = 6)  # open a new null device
  grDevices::dev.control("enable")  # turn on recording for the null device
  
  
  
  #Adjust bottom margin (first vector element) to give room for explanation
  #Order is bottom, left, top, right
  par(mar=c(12, 7, 1, 1) + 0.1, #margin, units are lines of text
      cex = 0.8, #Scale fonts by 75% 
      las=1, #Make all axis labels horizontal
      mgp = c(2, 0.1, 0), #Move axis labels closer to axes
      xpd=TRUE, 
      bg="white") 
  
  
  
  plot(gagedPeaks$plot_pos_qnorm, gagedPeaks$peak_va, col="cyan", pch=16, 
       ylim = c(min(breakPoints), max(breakPoints)),
       xlim = c(min(AEPs_qnorm), max(AEPs_qnorm)),
       #main = ptitle,
       xlab = xTitle, 
       ylab = "", 
       #yaxs = "i", #Needed to remove 4% inner margin
       yaxt='n', 
       xaxt='n',
       log="y")
  
  title(ylab = "Annual peak discharge, in cubic feet per second", mgp = c(4.5, 0.1, 0))
  
  #Take care of x-axis formatting
  #Bottom (side 1)
  axis(side=1, at=AEPs_qnorm,tck = 0.02, labels=AEPlabels)
  #Top (side 3)
  axis(side=3, at=AEPs_qnorm,tck = 0.02, labels=FALSE)
  
  
  #Take care of y-axis formatting
  #Left side (side 2)
  axis(side=2, at=breakPoints,tck = 0.02, labels=formatC(breakPoints,
                                                         format="fg",
                                                         big.mark=','))
  rug(minorTicks, ticksize = 0.01, side = 2)
  
  #Right side (side 4)
  axis(side=4, at=breakPoints, tck = 0.02, labels= FALSE)
  rug(minorTicks, ticksize = 0.01, side = 4)
  
  #Plot PILFs and historic peaks
  points(PILFs$plot_pos_qnorm, PILFs$peak_va, col="black", pch=1, lwd=1.5)
  points(historicPeaks$plot_pos_qnorm, historicPeaks$peak_va, bg="purple", col="black", pch=24)
  
  #Plot interval and censored peaks
  
  #Base R plots don't support error bars, but do support arrows
  #with arrow heads that look like error bars
  arrows(intervalPeaks$plot_pos_qnorm, intervalPeaks$ql, 
         intervalPeaks$plot_pos_qnorm, pmin(intervalPeaks$qu, 10^maxBreak),
         length=0.05, angle=90, code=3, col="cyan", lwd=2, xpd=FALSE)
  
  if(nrow(cPeaks) > 0){
    arrows(cPeaks$plot_pos_qnorm, cPeaks$q_bound,
           cPeaks$plot_pos_qnorm, cPeaks$q_inf,
           length=0.05, angle=90, code=1, col="grey", xpd=FALSE)
  }
  
  
  #Plot fitted distribution
  lines(-qnorm(LP3fit$EXC_Prob), LP3fit$Estimate, col="red", lwd=2)
  lines(-qnorm(LP3fit$EXC_Prob), LP3fit$Conf_Low, col="blue", lwd=2)
  lines(-qnorm(LP3fit$EXC_Prob), LP3fit$Conf_Up, col="blue", lwd=2)
  
  legend("bottomleft", 
         inset = c(0.02, -0.22 - 0.04*leg_rows), #Hold top of legend approximately constant. Need to fine-tune math to keep truly fixed.
         legend=leg_labels[leg_used],
         pch=leg_points[leg_used],
         lty=leg_lines[leg_used],
         pt.bg=leg_bg[leg_used],
         col = leg_col[leg_used],
         lwd = leg_lwd[leg_used],
         title="EXPLANATION", #Needed for legend under plot
         bg="white",
         bty="n", #Don't draw box around legend
         ncol = 1,
         box.lty = 1, # Line type of the box
         box.lwd = 1, # Width of the line of the box
         box.col = "black")
  
  mtext(infoText, 
        side = 1, #bottom margin
        adj = 0, #left align
        line = 10, #line from inner, counting outwards
        at = 1, #position on x axis (which is is normal transformed)
        cex = 0.8)
  
  #Save the plot to variable
  p <- recordPlot()
  dev.off()
  return(p)
}
  
