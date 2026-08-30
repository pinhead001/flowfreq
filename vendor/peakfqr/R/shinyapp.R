################################################################################
#           
#                          PeakFQ through R Shiny
#                   
#                              December 2023
#
################################################################################

################################# SET UP #######################################

#options(scipen = 3)

# declare global variables for check tool
utils::globalVariables(c(
  "psf_int_init", 
  "psf_peaks_init", 
  "psf_int_init",
  "psf_opt_init",
  "psf_peaks_init",
  "psf_thresh_init",
  "data_flow_int_df_init",
  "station_specs_df_init",
  "Qmin",
  "Qmax"
))

################################################################################

#' peakfq Shiny App
#'
#' Launch Shiny app user interface for peakfq
#' 
#' @examples
#' \dontrun{
#' #Shiny app will launch in web browser
#' peakfq_app() 
#' }
#' 
#'
#' @export
#' 
#' @import shiny
#' @import tidyverse
#' @import tibble
#' @import shinythemes
#' @import dataRetrieval
#' @import dplyr
#' @import rhandsontable
#' @import shinyFiles
#' @import shinyBS
#' @importFrom stringr str_replace_all
#' @importFrom shinyalert shinyalert
#' @importFrom shinyjs useShinyjs
#' @importFrom shinyjs toggle
#' @importFrom fs path_home
#' @importFrom rlang .data
#' @importFrom shinyBS bsButton bsPopover
#' @importFrom utils globalVariables
peakfq_app <- function(){
  if (interactive()) {
    library(peakfq) #Needed to get peakfq functions for app
    shiny::runApp(appDir = system.file("app", package = "peakfq"),
                  launch.browser = TRUE)
  } else {
    shiny::shinyAppDir(appDir = system.file("app", package = "peakfq"))
  }
}




