##########        PSF dataframe creation functions

########  Function to create a siteQT dataframe for all sites 
#' Create SiteQT Dataframe 
#' @description
#' Function to run the siteQT function for multiple sites at once, integrating
#' user specifications. See documentation for "siteQT" for additional information.   
#' @param pf_data dataframe of peak flow data for one or multiple sites 
#' @param psf_int optional, dataframe of interval specifications for sites included in pf_data
#' @param psf_thresh  optional, dataframe of threshold specifications for sites included in pf_data
#' @param psf_peaks optional, dataframe of peak specifications for sites included in pf_data
#' @param psf_opt optional, dataframe of processing specifications for sites included in pf_data
#' @param keepNoInfo optional, logical value of whether rows with no information (flow interval of zero to infinity) should be kept
#'
#' @keywords internal
#' 
#' @import dplyr
#' @importFrom rlang .data
#'
#' 
all_siteQT <- function(pf_data, psf_int = data.frame(), psf_thresh = data.frame(), psf_peaks = data.frame(), psf_opt = data.frame(), keepNoInfo = FALSE) {
  
  allsites_qt <- data.frame(matrix(ncol = 9, nrow = 0, dimnames = list(NULL, c("site_no", "peak_WY", "peak_va", "peak_cd", "ql", "qu", "tl", "tu", "dtype"))))
  
  #Add message if dataframe does not have rows. Not error since need blank siteQT.
  if(nrow(pf_data) == 0){
    message("No data rows provided.")
  }
  
  if(nrow(pf_data) > 0){
    #Generate psfs if not provided, as needed 
    
    #If no intervals: Leave blank (nrows = 0)
    if(nrow(psf_int) == 0) {
      psf_int <- psf_int_init[0,]
    }
    
    #If no thresholds: populate with defaults
    if(nrow(psf_thresh) == 0) {
      psf_thresh <- pf_data %>% 
        group_by(.data$site_no) %>% 
        summarise(label = "Interval",
                  start = min(.data$peak_WY), 
                  end = max(.data$peak_WY)) %>% 
        mutate(min = 0, 
               max = Qmax, 
               comment = "Default")
    }  
    
    # Limnotech created the bit below, apparently to always show the full data
    # view on input plots. Don't think it's a desireable behavior, so turned off
    # #Handle the possibility that the psf_thresh from a spec file upload would have an end year less than the data end year, 
    # # which would prevent siteQT from generating an output
    # if(nrow(psf_thresh) > 0){
    #   
    #   sites <- unique(pf_data$site_no)
    #   
    #   for(site in sites) {
    #     sel_site <- site
    #     data_max_yr <- pf_data %>%
    #       filter(site_no == sel_site) %>%
    #       select(peak_WY) %>%
    #       max()
    #     
    #     psf_end_yrs <- psf_thresh %>%
    #       filter(site_no == sel_site) %>%
    #       select(end) 
    #     
    #     psf_max_yr <- ifelse(nrow(psf_end_yrs) == 0, NA, max(psf_end_yrs))
    #     
    #     
    #     #If data end year is later than psf end year, update 'end' on default row to match data end year
    #     if(nrow(psf_end_yrs) != 0 & (data_max_yr > psf_max_yr)){
    #       row_i <- which(psf_thresh$site_no == sel_site & tolower(psf_thresh$comment) == "default")
    #       psf_thresh[row_i, 'end'] <- as.integer(data_max_yr)
    #     }
    #   }
    # }
    
    #If no peaks: Leave blank (nrows = 0)
    if(nrow(psf_peaks) == 0) {
      psf_peaks <- psf_peaks_init[0,]
    }
    
    #If no options/specs: Leave blank
    if(nrow(psf_opt) == 0) {
      
      site_qt_list <- lapply(unique(pf_data$site_no), siteQT, pf_data, psf_int, psf_thresh, psf_peaks)
      
    }else{ #If psf_opt is provided, ensure that Urb/Reg option is used in siteQT()
      
      #Need to handle Limnotech not using 'Urb/Reg' column same as siteQT()
      
      urb_reg_specs <- data.frame(site_no = psf_opt$site_no)
      urb_reg_specs$`Urb/Reg` <- "YES"
      urb_reg_specs$`Urb/Reg`[is.na(psf_opt$`Urb/Reg`) | psf_opt$`Urb/Reg` == FALSE] <- "NO"
      
      site_qt_list <- lapply(unique(pf_data$site_no), siteQT, pf_data, psf_int, psf_thresh, psf_peaks, TRUE, urb_reg_specs, keepNoInfo)
    }
    
    
    #Add site numbers to each datafames
    site_qt_list <-lapply(site_qt_list,
                          function(x){add_column(x, site_no = attributes(x)$site_no, .before = 1)})
    
    allsites_qt <-  dplyr::bind_rows(site_qt_list)
    
    
  }  
  
  return(allsites_qt)
}  


###   Function to create psf_thresh from user inputs
#' Create PSF Thresholds 
#' @description
#' Create dataframe of threshold specifications from perceptible ranges table 
#' @param perc_ranges_df Populated dataframe of perceptible ranges
#'
#' @return dataframe of threshold specifications
#' @keywords internal 
#' 
#' @import dplyr
#' @importFrom rlang .data
#'
#' 
create_psf_thresh <- function(perc_ranges_df) {
  
  psf_thresh <- perc_ranges_df %>% 
    filter(!(is.na(.data$site_no)| is.na(.data$start) | is.na(.data$end) | is.na(.data$min) | is.na(.data$max) | is.na(.data$comment))) %>% 
    filter(!(.data$start == "" | .data$end == "" | .data$min == "" | .data$max == "" | .data$comment == "")) %>%
    mutate(label = "PCPT_Thresh", 
           start = as.integer(.data$start),
           end = as.integer(.data$end)) %>% 
    select(names(psf_thresh_init)) %>%
    
    
    return(psf_thresh)
  
}

###   function to create psf_peaks from user inputs
#' Create PSF Peaks 
#' @description
#' Create dataframe of peak specifications from data/flow intervals table
#' @param data_flow_int_df Populated dataframe of data/flow intervals
#'
#' @return dataframe of peak specifications
#' @keywords internal
#' 
#' @import dplyr
#' @importFrom rlang .data
#'
#' 
create_psf_peaks <- function(data_flow_int_df) {
  
  psf_peaks <- data_flow_int_df %>%
               filter(!(is.na(.data$site_no)| is.na(.data$peak_WY) | is.na(.data$peak_va) | is.na(.data$interval_low) | is.na(.data$interval_up) | is.na(.data$comment))) %>% 
               filter(!(.data$site_no == "" | .data$peak_WY == "" | .data$peak_va == "" | .data$interval_low == "" | .data$interval_up == "" | .data$comment == "")) %>%
               mutate(label = "Peak", 
                      peak_WY = as.integer(.data$peak_WY)) %>%
               select(names(psf_peaks_init)) %>%
               filter(!is.na(.data$site_no)) %>% #Remove blank rows from rhot
               filter(comment != "") #Remove all standard entries.
  
  return(psf_peaks)
  
}

###   function to create psf_int from user inputs
#' Create PSF Intervals
#' @description
#' Create dataframe of interval specifications from data/flow intervals table
#' @param data_flow_int_df Populated dataframe of data/flow intervals
#'
#' @return dataframe of interval specifications
#' @keywords internal
#' 
#' @import dplyr
#' @importFrom rlang .data
#'
#' 
create_psf_int <- function(data_flow_int_df) {
  
  psf_int <- data_flow_int_df %>%
             filter(!(is.na(.data$site_no)| is.na(.data$peak_WY) | is.na(.data$interval_low) | is.na(.data$interval_up) | is.na(.data$comment))) %>% 
             filter(!(.data$site_no == "" | .data$peak_WY == "" | .data$interval_low == "" | .data$interval_up == "" | .data$comment == "")) %>%
             mutate(label = "Interval", 
                    peak_WY = as.integer(.data$peak_WY)) %>%
             select(names(psf_int_init)) %>%
             filter(!is.na(.data$site_no)) %>% 
             filter(comment != "") %>% #Remove all standard entries.
             filter(!is.na(.data$comment))
  
  return(psf_int)
  
}

###   Function to create psf_opt df from user inputs
#' Create PSF Opt
#' @description
#' Create dataframe of processing option specifications from populated station specifications, perceptible ranges, and data/flow intervals dataframes.
#' @param station_specs_df Populated dataframe of processing option specifications
#' @param perc_ranges_df Populated dataframe of perceptible ranges
#' @param data_flow_int_df Populated dataframe of data/flow intervals
#' @param psf_peaks Peak specifications dataframe 
#' 
#'
#' @return dataframe  of processing option specifications
#' @keywords internal
#' 
#' @import dplyr
#' @importFrom rlang .data
#'
#' 
create_psf_opt <- function(station_specs_df, perc_ranges_df, data_flow_int_df, psf_peaks = data.frame()) {
  
  #Reformat all data to 
  psf_thresh <- create_psf_thresh(perc_ranges_df)
  # psf_peaks <- create_psf_peaks(data_flow_int_df) #Cannot be changed in app; capturing causes problems.
  psf_peaks <- psf_peaks
  psf_int <- create_psf_int(data_flow_int_df)
  
  from_station_specs <- station_specs_df %>%
    select(names(psf_opt_init)[names(psf_opt_init) %!in% c("Interval", "Peak", "PCPT_Thresh")]) #Retain the columns that belong in the psf_opt df
  
  thresh_concat <- psf_thresh %>% 
    mutate(PCPT_Thresh = ifelse(is.na(comment), NA, paste(.data$label, .data$start, .data$end, .data$min, .data$max, .data$comment))) %>%
    select(.data$site_no, .data$PCPT_Thresh) %>%
    group_by(.data$site_no) %>%
    mutate(PCPT_Thresh = paste(.data$PCPT_Thresh, collapse = ";")) %>% 
    unique() %>%
    ungroup()
  
  int_concat <- psf_int %>% 
    mutate(Interval = ifelse(is.na(.data$comment), NA, paste(.data$label, .data$peak_WY, .data$interval_low, .data$interval_up, .data$comment))) %>%
    select(.data$site_no, .data$Interval) %>%
    group_by(.data$site_no) %>%
    mutate(Interval = paste(.data$Interval, collapse = ";")) %>% 
    unique()  %>%
    ungroup()
  
  if(nrow(int_concat) == 0){
    
    int_concat[1:nrow(from_station_specs),] <- NA
    int_concat$site_no <- from_station_specs$site_no
    
  }
  
  peak_concat <- psf_peaks %>%
    mutate(Peak = ifelse(is.na(.data$comment), NA, paste(.data$label, .data$peak_WY, .data$peak_va, .data$comment))) %>%
    select(.data$site_no, .data$Peak) %>%
    group_by(.data$site_no) %>%
    mutate(Peak = paste(.data$Peak, collapse = ";")) %>% 
    unique() %>%
    ungroup()
  
  if(nrow(peak_concat) == 0){
    
    peak_concat[1:nrow(from_station_specs),] <- NA
    peak_concat$site_no <- from_station_specs$site_no
    
  }
  
  psf_opt <- from_station_specs %>% 
    merge(thresh_concat, by = "site_no", all = TRUE) %>%
    merge(int_concat, by = "site_no", all = TRUE) %>%
    merge(peak_concat, by = "site_no", all = TRUE) %>%
    select(names(psf_opt_init)) %>%
    unique() %>%
    filter(!is.na(.data$site_no))
  
  return(psf_opt)
  
}
