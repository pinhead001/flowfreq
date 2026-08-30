
#' Plot Peak Flow Data
#'
#' Function to plot peak-flow input data
#'
#' @param QT dataframe - dataframe containing water years, flow intervals, perception thresholds, and identifier for historic peaks in the columns:
#'   1. peak_WY integer - water year
#'   2. ql numeric - lower flow interval values
#'   3. qu numeric - upper flow interval values
#'   4. tl numeric - lower perception threshold values
#'   5. tu numeric - upper perception threshold values
#'   6. dtype integer - indicator if value is historic peak. 1 = historic peak, 0 = not historic peak
#'   OR
#'      historic logical - indicator if value is historic peak. TRUE = historic peak, FALSE = not historic peak
#'      If dtype and historic are both supplied, the historic column will be used. 
#'      
#'      
#' @param siteNum character - site number to print on plot
#' @param siteName character - site name to print on plot
#' @param ytrans character - y-axis transformation. Set to "log10" for log axis
#' @param ptitle character - text to use for plot title
#' @param nodatyears integer - vector of user-specified no data years to color differently from other missing data years
#' 
#' @returns Base R graphics object with plot of peak-flow data
#' 
#' @details
#' The plot generated is optimized for viewing and saving in a window
#' 9 inches wide and 6 inches high. The plot can be resized
#' but excessive resizing may result in overlapping text or other display
#' issues.
#' 
#' 
#' 
#' @examples
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
#'                  dtype = rep(0, length(qu)), #No historic peaks
#'                  peak_WY = c(1:10)
#'                  )
#'  
#' #Create plot
#' p <- pkPlot(QT, siteNum = "12345678", siteName="Example Site",
#'             ytrans="identity", ptitle="Some Input Data")
#'  
#' 
#' 
#'
#' @export
#' 
pkPlot <- function(QT, siteNum="", siteName="", ytrans="identity", ptitle="", nodatyears=c()){
  
  if(nrow(QT) == 0){
    
    warning("Supplied dataframe has zero rows. Returning empty plot.")
    
    plot(0, 0, col="black", pch=1,
         xlim = c(1950, 2000),
         ylim = c(0, 1000),
         xlab = "Water Year", 
         ylab = "Annual peak discharge, in cubic feet per second", 
         family = "sans",
         yaxs = "i", #Needed to remove 4% inner margin
         yaxt="n", 
         xaxt="n")
    
    p <- recordPlot()
    dev.off()
    return(p)
    
  }
  
  #Handle censored peaks
  QT$censored <- FALSE
  QT$censored[QT$ql != QT$qu] <- TRUE
  
  cPeaks <- QT[QT$ql != QT$qu, c("peak_WY", "ql", "qu")]
  #Separate greater than peaks
  gtPeaks <- cPeaks[cPeaks$qu >= Qmax,]
  cPeaks <- cPeaks[cPeaks$qu < Qmax,]
  
  #Set logical for whether each peak is historic
  QT$historicPeaks <- NA
  QT$gagedPeaks <- NA
  
  if("historic" %!in% colnames(QT) & "dtype" %in% colnames(QT)){
    QT$historic <- as.logical(QT$dtype)
  }
  else if("historic" %!in% colnames(QT) & "dtype" %!in% colnames(QT)){
    QT$historic <- FALSE #No info to designate historic peaks
  }
  
  QT[!QT$historic & !QT$censored, "gagedPeaks"] <- QT[ !QT$historic & !QT$censored, "ql"]
  #QT[QT$historic & !QT$censored, "historicPeaks"] <- QT[QT$historic & !QT$censored, "ql"]
  hPeaks <- QT[QT$historic & !QT$censored, c("peak_WY", "ql", "qu")] #Get historic peaks
  
  #Get gap years
  gaps <- min(QT$peak_WY):max(QT$peak_WY)
  gaps <- gaps[gaps %!in% QT$peak_WY]
  
  PTs <- QT[, c("peak_WY", "tl", "tu")]
  
  if(length(gaps) > 0){
    
    #If gaps are present, add these years to QT dataframe
    PTs <- rbind(PTs,data.frame(peak_WY = gaps,
                               tl = Qmax,
                               tu = Qmax)
                 )
    
    PTs <- PTs[order(PTs$peak_WY),]
	
	#Special handling of coded missing Data
	nodat <- gaps[gaps %in% nodatyears]
	gaps_plot <- gaps[gaps %!in% nodatyears]
	

  }
  else{
    nodat <- integer(0)
	  gaps_plot <- integer(0)
  }
  
  
  #Process PTs
  
  #Get years where PT changes
  tlStarts <- c(1, diff(PTs$tl))
  tuStarts <- c(1, diff(PTs$tu))
  tlEnds <- c(diff(PTs$tl), 1)
  tuEnds <- c(diff(PTs$tu), 1)
  
  
  tl <- data.frame(startYear = PTs$peak_WY[tlStarts != 0], 
                   endYear = PTs$peak_WY[tlEnds != 0], 
                   thresh = PTs$tl[tlStarts != 0])
  
  tu <- data.frame(startYear = PTs$peak_WY[tuStarts != 0], 
                   endYear = PTs$peak_WY[tuEnds != 0], 
                   thresh = PTs$tu[tuStarts != 0])
  
  
  
  #Get max break on the axis from pretty()
  #Need 1.1 to make sure max value doesn't end up on top edge
  maxBreak <- max(QT$qu[QT$qu < Qmax]) 
  maxBreak <- max(pretty(c(0,maxBreak*1.05))) 
  
  
  
   
  
  xTitle <- paste("Water year \n Station - ", siteNum, siteName)
  
  ### Prepare the legend items ###
  leg_labels <- c("Gaged peak discharge",
                  "Historic peak discharge", 
                  "Censored peak discharge", 
                  "Lower perception threshold", 
                  "Upper perception threshold", 
                  "No data (program default)", 
				          "No data (user specified)")
  leg_points <- c(1, 24, NA, NA, NA, NA, NA)
  leg_lines <- c(0, 0, 1, 1, 1, 1, 1)
  leg_bg <- c(NA, "black", NA, NA, NA, NA, NA)
  leg_col <- c("black", "black", "black", "darkgoldenrod1", "darkorchid1", "indianred", "lightgoldenrod1")
  leg_lwd <- c(1,1,1,4,4,4,4)
  
  #Find out which potential legend items are actually used
  leg_used <- as.logical(c(nrow(QT) > 0,
                           nrow(hPeaks) > 0,
                           nrow(cPeaks) > 0, 
                           TRUE, 
                           TRUE, 
                           length(gaps_plot) > 0, 
						               length(nodat) > 0))
  
  
  
  
  
  ### Do a dummy plot to get the size of the legend ###
  
  #Set up to record plot without displaying
  cur_dev <- grDevices::dev.cur()   # store current device
  pdf(NULL, width = 7.5, height = 5)  # open null device
  null_dev <- grDevices::dev.cur()  # store null device
  
  # make sure we always clean up properly, even if something causes an error
  on.exit({
    grDevices::dev.off(null_dev)
    if (cur_dev > 1) grDevices::dev.set(cur_dev) # only set cur device if not null device
  })
  
  
  plot(QT$peak_WY, QT$gagedPeaks, col="black", pch=1, 
       ylim = c(0, maxBreak),
       yaxs = "i", #Needed to remove 4% inner margin
       yaxt='n', 
       xaxt='n')
  
  #Add legend to dummy plot and get its height
  legend_size <- legend("topleft",
                        legend=leg_labels[leg_used],
                        pch=leg_points[leg_used],
                        lty=leg_lines[leg_used],
                        pt.bg=leg_bg[leg_used],
                        col = leg_col[leg_used],
                        lwd = leg_lwd[leg_used],
                        title="Explanation",
                        ncol = 2)$rect$h
  
  grDevices::dev.off(null_dev) #close null device
  
  #re-adjust y axis to accommodate legend
  maxBreak <- maxBreak + legend_size*1.16
  tu$thresh_plot <- pmin(tu$thresh, maxBreak)
  
  ### Now do the actual plot ###
  
  #Set up to record plot without displaying

  pdf(NULL, width = 9, height = 6)  # open a new null device
  grDevices::dev.control("enable")  # turn on recording for the null device

  
  #Setting fonts cannot be simple
  #windowsFonts(Arial = windowsFont("Arial"))
  
  #Adjust bottom margin (first vector element) to give room for explanation
  par(mar=c(4, 7, 2, 2) + 0.1, #units are lines of text
      las=1, #Make all axis labels horizontal
      mgp = c(3, 0.5, 0), 
      bg="white"
      ) 
  

  
  plot(QT$peak_WY, QT$gagedPeaks, col="black", pch=1, 
       ylim = c(0, maxBreak),
       main = ptitle,
       xlab = xTitle, 
       ylab = "", 
       #family = "Arial",
       yaxs = "i", #Needed to remove 4% inner margin
       yaxt='n', 
       xaxt='n')
  
  if(ytrans == "log10"){
    
    #Need to find lower bound for y axis
    
    
    
    plot(QT$peak_WY, QT$gagedPeaks, col="black", pch=1, 
         ylim = c(10^min(floor(log10(QT$ql[QT$ql > Qmin]))), maxBreak),
         main = ptitle,
         xlab = xTitle, 
         ylab = "", 
         #family = "Arial",
         yaxs = "i", #Needed to remove 4% inner margin
         yaxt='n', 
         xaxt='n',
         log="y")
  }
  
  title(ylab = "Annual peak discharge, in cubic feet per second", mgp = c(4.5, 0.1, 0))
  
  #Take care of x-axis formatting
  #Bottom (side 1)
  axis(side=1, at=axTicks(1),tck = 0.01, labels=formatC(axTicks(1),
                                                        format="fg"))
  #Top (side 3)
  axis(side=3, at=axTicks(1),tck = 0.01, labels=FALSE)
  
  
  #Take care of y-axis formatting
  #Left side (side 2)
  axis(side=2, at=axTicks(2),tck = 0.01, labels=formatC(axTicks(2),
                                             format="d",
                                             big.mark=','))
  #Right side (side 4)
  axis(side=4, at=axTicks(2), tck = 0.01, labels= FALSE)

  #Plot perception thresholds
  segments(tl$startYear - 0.5, tl$thresh, tl$endYear + 0.5, tl$thresh, col="darkgoldenrod1", lwd=4)
  segments(tu$startYear - 0.5, tu$thresh_plot, tu$endYear + 0.5, tu$thresh_plot, col="darkorchid1", lwd=4)
  
  segments(gaps_plot, rep(maxBreak - legend_size*1.16, length(gaps_plot)), gaps_plot, rep(0, length(gaps_plot)), col="indianred", lwd=4)
  segments(nodat, rep(maxBreak - legend_size*1.16, length(nodat)), nodat, rep(0, length(nodat)), col="lightgoldenrod1", lwd=4)
  
  points(hPeaks$peak_WY, as.numeric(hPeaks$ql), 
         col="black", bg = "black", pch=24)
  
  #Base R plots don't support error bars, but do support arrows
  #with arrow heads that look like error bars
  arrows(cPeaks$peak_WY, cPeaks$ql, cPeaks$peak_WY, pmin(cPeaks$qu, maxBreak), length=0.05, angle=90, code=3, xpd=FALSE)
  arrows(gtPeaks$peak_WY, gtPeaks$ql, gtPeaks$peak_WY, pmin(gtPeaks$qu, maxBreak), length=0.05, angle=90, code=1, xpd=FALSE)
  
  

  legend("topleft", 
         inset = c(0.02, 0.02),
         legend=leg_labels[leg_used],
         pch=leg_points[leg_used],
         lty=leg_lines[leg_used],
         pt.bg=leg_bg[leg_used],
         col = leg_col[leg_used],
         lwd = leg_lwd[leg_used],
         title="EXPLANATION", #Needed for legend under plot
         bg="white",
         box.col = "white",
         #text.width = 60, #Used to set width of legend columns in x-axis coordinates
         #bty="n",
         ncol = 2)
  
  #Save the plot to variable
  p <- recordPlot()
  dev.off()
  return(p)
  
}



