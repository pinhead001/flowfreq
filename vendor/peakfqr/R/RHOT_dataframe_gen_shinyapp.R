####   Utility function to add comments to data frame containing remark codes
#' Add Remark Comments 
#' @description
#' Add comments to data frame based on included remark codes. See the following
#' link for additional information on remark codes: 
#' <a href="https://nwis.waterdata.usgs.gov/nwis/peak?help#flow_qual_cd">"https://nwis.waterdata.usgs.gov/nwis/peak?help#flow_qual_cd" </a>
#' @param df A dataframe containing a column with remark codes 
#'
#' @return Input dataframe with the addition of descriptive comments based on remark codes. 
#' @keywords internal
#'
#' @import dplyr
#' @importFrom rlang .data
#' 
add_cd_comments <- function(df){
  
  commented_df <- df %>% mutate(comment4 = ifelse(grepl("4", as.character(.data$peak_cd)) == TRUE, "PEAK < STATED VALUE", ""),
                                comment8 = ifelse(grepl("8", as.character(.data$peak_cd)) == TRUE, "PEAK > STATED VALUE", ""),
                                commentO = ifelse(grepl("O", as.character(.data$peak_cd)) == TRUE, "OPPORTUNISTIC PEAK", ""),
                                comment3 = ifelse(grepl("3", as.character(.data$peak_cd)) == TRUE, "DAM BREAK", ""),
                                comment7 = ifelse(grepl("7", as.character(.data$peak_cd)) == TRUE, "HISTORIC PEAK", "")) %>%
    mutate(comment = paste(.data$comment4, .data$comment8, .data$commentO, .data$comment3, .data$comment7, sep = ", ") %>%
             gsub(pattern = " ,", replacement = "") %>% #drop strings of just commas (i.e., blank comments)
             trimws() %>% #Drop whitespace
             gsub(pattern = "^, ", replacement = "") %>% #Replace leading commas
             gsub(pattern = ",$", replacement = "")) %>% #Replace trailing commas
    select(names(df), comment)
  
  return(commented_df)  
  
}    

####   Utility function to identify and group contiguous values within perc ranges
#' Group Contiguous Values
#' @description
#' Identify and consolidate contiguous values within perceptible ranges
#' @param perc_ranges_df dataframe of perceptible ranges 
#'
#' @return dataframe of perceptible ranges for a set of given years 
#' @keywords internal
#'
#' @import dplyr
#' @importFrom rlang .data
#' 
group_contiguous_perc_ranges <- function(perc_ranges_df) {
  
  contig_perc_ranges <- data.frame(matrix(nrow = 0, ncol = 6, dimnames = list(NULL, c("site_no", "start", "end", "min", "max", "comment"))))  #format for tacking onto perc ranges output
  
  for(site in unique(perc_ranges_df$site_no)) {
    
    station_ranges_update <- data.frame(matrix(nrow = 0, ncol = 6, dimnames = list(NULL, c("site_no", "start", "end", "min", "max", "comment"))))  #format for tacking onto perc ranges output
    
    site_pr_df <- perc_ranges_df %>% filter(.data$site_no == site)
    
    #Find contiguous values to group - check if contiguous with preceding row
    contig_check <- site_pr_df %>% 
      group_by(min, max) %>% 
      mutate(lag = start - lag(start, default = 1),
             contig = ifelse(lag == 1 | lag == (start-1), "contiguous", "not contiguous")) #identify if contiguous with preceding group
    
    contig_check$group <- NA
    group = 1
    
    for(row in 1:nrow(contig_check)) {
      
      sel_row = contig_check[row,]
      group = ifelse(sel_row$contig == "not contiguous", group + 1, group) # increment if not contiguous
      
      contig_check$group[row] <- group
    }
    
    #consolidate contiguous rows
    station_ranges_update <- contig_check %>% 
                            group_by(group) %>%
                            summarise(site_no = first(.data$site_no),
                                      start = min(.data$start),
                                      end = max(.data$end),
                                      min = first(.data$min),
                                      max = first(.data$max),
                                      comment = first(.data$comment))
    
  } 
  
  contig_perc_ranges <- contig_perc_ranges %>% rbind(station_ranges_update) %>% ungroup()
  
  return(contig_perc_ranges)
}

####  Process rows with 'largest since' values and return the appropriate perceptible ranges
#' Process 'Largest Since' Values
#' @description
#' Process rows with 'largest since' values and return the 
#' appropriate perceptible ranges. 
#' 
#' Identify all discharges with the "year_last_pk" 
#' column populated and sort these from largest to smallest based on peak 
#' discharge ("peak_va" column). Starting with the largest peak discharge (Q1) , 
#' apply a lower perception threshold of Q1 and upper perception threshold of infinity 
#' to all years without discharge values going back to the year AFTER the year associated 
#' with Q1 in the "year_last_pk" column.  
#' Repeat step 2 with the next largest peak discharge with a "year_last_pk" value, 
#' overwriting the perception thresholds from the previous discharge if necessary. 
#' Continue until the all years with a "year_last_peak" value have been processed.
#' @param pf_data dataframe of formatted peak flow data
#'
#' @return Dataframe of perceptible ranges based on "largest since" values
#' @keywords internal
#'
#' @import dplyr
#' @importFrom rlang .data
#' 
process_largest_since <- function(pf_data) {
  
  # Identify all discharges with the "year_last_pk" column populated and sort these from largest to smallest based on peak discharge ("peak_va" column). (pull this back into pf_data!)   
  # Starting with the largest peak discharge (Q1) , apply a lower perception threshold of Q1 and upper perception threshold of infinity to all years without discharge values going back to the year AFTER the year associated with Q1 in the "year_last_pk" column.    
  # Repeat step 2 with the next largest peak discharge with a "year_last_pk" value, overwriting the perception thresholds from the previous discharge if necessary. Continue until the all years with a "year_last_peak" value have been processed. 
  
  #Dataframe to gather all station rows
  perc_ranges_largest_since <- data.frame(matrix(nrow = 0, ncol = 6, dimnames = list(NULL, c("site_no", "start", "end", "min", "max", "comment"))))  #format for tacking onto perc ranges output
  
  year_last_rows_all <- pf_data %>% 
                        filter(!is.na(.data$year_last_pk)) 
  
  if(nrow(year_last_rows_all) == 0) {
    
    return(perc_ranges_largest_since)
    
  } else { 
    
    for(site in unique(year_last_rows_all$site_no)){ 
      
      sel_data <- pf_data %>% filter(.data$site_no == site)
      
      year_last_rows <- sel_data %>% 
        filter(!is.na(.data$year_last_pk)) %>%
        arrange(desc(as.numeric(.data$peak_va)))
      
      
      #Dataframe to gather all instance rows
      station_ranges_update <- data.frame(matrix(nrow = 0, ncol = 6, dimnames = list(NULL, c("site_no", "start", "end", "min", "max", "comment"))))  #format for tacking onto perc ranges output
      
      for(i in 1:nrow(year_last_rows)) {
        
        row_ranges_update <- data.frame(matrix(nrow = 0, ncol = 6, dimnames = list(NULL, c("site_no", "start", "end", "min", "max", "comment"))))  #format for tacking onto perc ranges output
        
        sel_row <- year_last_rows[i,]
        
        Q1 = sel_row$peak_va #Largest peak discharge
        year_Q1 = as.integer(sel_row$peak_WY)
        year_referenced = as.integer(sel_row$year_last_pk)
        year_after = year_referenced + 1 
        
        year_span <- year_after:year_Q1 #all possible years going back to year after year associated with Q1
        present_years <- sel_data %>% filter(!is.na(.data$peak_va))
        present_years <- present_years$peak_WY #years within range above that have a peak value in the provided data
        missing_years <- year_span[year_span %!in% present_years]#years within range above that DO NOT have a peak value in the provided data. These require manual handling. 
        
        if(length(missing_years) == 0) {next} #Patch based on sites where "year last" is same as peak_WY
        
        row_ranges_update[1:length(missing_years),] <- NA #Ensure sufficient rows
        
        row_ranges_update$site_no = site 
        row_ranges_update$start = missing_years 
        row_ranges_update$end = missing_years
        row_ranges_update$min = Q1
        row_ranges_update$max = Qmax
        row_ranges_update$comment = paste0(year_Q1, " LARGEST SINCE ", year_referenced) 
        
        station_ranges_update <- station_ranges_update %>% rbind(row_ranges_update)
        
      } #end row loop
      
      if(nrow(station_ranges_update) == 0) {next}
      
      #Keep only the lower min value for a given year
      station_ranges_update <- station_ranges_update %>% 
        group_by(start, end) %>% 
        arrange(min) %>%
        summarise(site_no = first(.data$site_no),
                  start = first(.data$start),
                  end = first(.data$end),
                  min = min(.data$min),
                  max = max(.data$max),
                  comment = first(.data$comment)) %>%
        group_contiguous_perc_ranges() %>% 
        select(!"group")
      
      perc_ranges_largest_since <- perc_ranges_largest_since %>% rbind(station_ranges_update) 
      
    }   #end station loop
    
    return(perc_ranges_largest_since)
  }
} 

####   Process code 7 periods and return the appropriate perceptible ranges
#' Process Code Seven
#' @description
#' Process code 7 periods and return the appropriate perceptible ranges
#' @param pf_data A dataframe of peak flow data containing a remarks code column 
#'
#' @return Dataframe of perceptible ranges corresponding to data marked with code 7. 
#' @keywords internal
#' 
#' @import dplyr
#' @importFrom rlang .data
#'
#' 
process_cd_7 <- function(pf_data) {
  
  # Identify each period of missing systematic data with a code 7 peak in that period
  # If one or more code 7 peaks occurs before the first systematic peak (peak with no code 7), the period is from the first code 7 peak to the year before the first systematic peak
  # If one or more code 7 peaks occurs after the last systematic peak (peak with no code 7), the period is from the year after the last systematic peak to the last code 7 peak
  # If one or more code 7 peaks occurs during a gap in systematic data, the period is from the year after the last systematic peak before the code 7 peak to the year before the first systematic peak after the code 7 peak
  # Gaps in systematic data that DO NOT have a code 7 peak in the gap are ignored
  # For each period identified, set the lower perception threshold for the period to the lowest code 7 peak discharge during that period and the upper perception threshold to infinity. 
  
  #Note on terminology - peaks without a code 7 or code O are called "Systematic" data. 
  #Since code O peaks have already had their discharge values removed, they shouldn't affect the determination of what data are Systematic for handling the code 7s. 
  
  pf_data_sys <- pf_data %>% 
    filter(grepl("O", as.character(.data$peak_cd)) == FALSE) %>% #Remove discharge values for code O peaks
    mutate(systematic = ifelse(grepl("7", as.character(.data$peak_cd)) == TRUE, FALSE, TRUE))
  
  cd7_periods <- data.frame(matrix(nrow = 0, ncol = 6, dimnames = list(NULL, c("site_no", "start", "end", "min", "max", "comment"))))  #format for tacking onto perc ranges output
  
  cd7_sites <- unique(pf_data_sys$site_no[grepl("7", as.character(pf_data_sys$peak_cd))])
  
  for(site in cd7_sites){ 
    
    site_periods <- data.frame(matrix(nrow = 0, ncol = 7, dimnames = list(NULL, c("site_no", "start", "end", "min", "max", "comment", "index"))))  #format for tacking onto perc ranges output
    n = 0
    
    sel_data <- pf_data_sys %>% filter(.data$site_no == site)
    
    cd7_rows <- sel_data %>%
      filter(grepl("7", as.character(.data$peak_cd)) == TRUE)
    
    cd7_indices <- which(grepl("7", as.character(sel_data$peak_cd)))
    
    systematic_indices <- which(!grepl("7", as.character(sel_data$peak_cd)))
    
    
    for(i in cd7_indices) { 
      
      #increment number of historic intervals
      n = n + 1
      
      sel_row <- sel_data[i,]
      period <- data.frame(matrix(nrow = 0, ncol = 7, dimnames = list(NULL, c("site_no", "start", "end", "min", "max", "comment", "index"))))  #format for tacking onto perc ranges output
      
      if (i < min(systematic_indices)){
        # If one or more code 7 peaks occurs before the first systematic peak (peak with no code 7), the period is from the first code 7 peak to the year before the first systematic peak
        
        interval <- sel_row$peak_WY:(sel_data[min(systematic_indices),]$peak_WY-1)
        
      } else if(i > max(systematic_indices)){
        # If one or more code 7 peaks occurs after the last systematic peak (peak with no code 7), the period is from the year after the last systematic peak to the last code 7 peak
        
        interval <- (sel_data[max(systematic_indices),]$peak_WY+1):(sel_row$peak_WY)
        
      } else { 
        # If one or more code 7 peaks occurs during a gap in systematic data, the period is from the year after the last systematic peak before the code 7 peak to the year before the first systematic peak after the code 7 peak
        
        i_sys_before <- systematic_indices[systematic_indices < i] %>% max()
        i_sys_after <- systematic_indices[systematic_indices > i] %>% min()
        
        interval <- (sel_data[i_sys_before,]$peak_WY+1):(sel_data[i_sys_after,]$peak_WY-1)
        
      }
      
      period <- period[1,] %>% mutate(site_no = site,
                                      start = min(interval),
                                      end = max(interval),
                                      min = sel_row$peak_va,
                                      max = Qmax,
                                      comment = paste0("HISTORIC PLACEHOLDER"),
                                      index = i) 
      
      site_periods <- site_periods %>% rbind(period)
      
    }  #end index loop
    
    #For each period identified, set the lower perception threshold for the period to the lowest code 7 peak discharge during that period and the upper perception threshold to infinity. "
    #Find contiguous values to group - check if contiguous with preceeding row
    contig_check <- site_periods %>% 
      mutate(lag = .data$index - lag(.data$index),
             contig = ifelse(lag == 1 | is.na(lag), "contiguous", "not contiguous")) #identify if contiguous with preceding group
    
    contig_check$group <- NA
    group = 1
    
    for(row in 1:nrow(contig_check)) {
      
      sel_row = contig_check[row,]
      group = ifelse(sel_row$contig == "not contiguous", group + 1, group) # increment if not contiguous
      
      contig_check$group[row] <- group
    }
    
    #consolidate contiguous rows - For each period identified, set the lower perception threshold for the period to the lowest code 7 peak discharge during that period and the upper perception threshold to infinity. 
    site_periods <- contig_check %>% 
      group_by(group) %>%
      summarise(site_no = first(.data$site_no),
                start = min(.data$start),
                end = max(.data$end),
                min = min(.data$min),
                max = Qmax,
                comment = paste0("HISTORIC ", first(group))) %>% 
      select(!"group") %>%
      ungroup
    
    cd7_periods <- cd7_periods %>% rbind(site_periods)
    row.names(cd7_periods) <- NULL
    
  } #end site loop
  
  return(cd7_periods)
  
} 

####   Function to generate station specifications df from input data and, where provided, specs
#' Generate Station Specifics
#' @description
#' Generate station specifications dataframe from input data and, 
#' where provided, specifications. 
#' @param pf_data A dataframe of peak flow data
#' @param psf_opt Optional, a dataframe containing processing specifications 
#' @param site_info Optional, a dataframe containing site name, latitude, and longitude
#'
#' @return A formatted dataframe of station specifications for display on the user interface
#' @keywords internal
#'
#' @import dplyr
#' @importFrom rlang .data
#' 
generate_station_specifications <- function(pf_data, psf_opt = data.frame(), site_info = data.frame()){ 
  
  if (nrow(pf_data) == 0) {
    stop(paste("No data rows have been provided"))
  }
  
  #Process data if provided
  if (nrow(pf_data) != 0) {
    
    #Summarise dataset by site (this provides the years & historical period)
    site_data_cols <- pf_data %>% 
                      group_by(.data$site_no) %>% 
                      summarise(BegYear = min(.data$peak_WY), 
                                EndYear = max(.data$peak_WY), 
                                HistPeriod = sum(!is.na(.data$peak_va))) %>% 
                      select(.data$site_no, .data$HistPeriod, .data$BegYear, .data$EndYear)
    
    collected_cols <- site_data_cols
    
    #pull in information from specifications, if provided 
    if (nrow(psf_opt) != 0){ 
      psf_cols <- psf_opt[names(psf_opt) %in% names(station_specs_df_init)]
      collected_cols <- merge(collected_cols, psf_cols, all.x = TRUE)
    } 
    
    #pull in information from site information (file or NWIS), if provided 
    if (nrow(site_info) != 0){ 
      site_info_cols <- site_info
      collected_cols <- merge(collected_cols, site_info_cols, all.x = TRUE)
    }
    
    station_specs_df <- merge(collected_cols, station_specs_df_init, all.x = TRUE) %>%
                        select(names(station_specs_df_init)) %>%
                        mutate(LOType = ifelse(is.na(.data$LOType), "MGBT", .data$LOType))#PILF test option on the first tab should default to the MGBT option.
    
    #Ensure that the order of sites is consistent with the input data
    site_order <- pf_data$site_no %>% unique()
    station_specs_df <- station_specs_df[match(site_order, site_data_cols$site_no),] %>%
                        unique()
    
    return(station_specs_df)
    
  }
}

####   Function to generate perceptible ranges df from input data and, where provided, specs
#' Generate Perceptible Ranges
#' @description
#' Generate perceptible ranges dataframe from input data and, where 
#' provided, specifications
#' @param pf_data A dataframe of peak flow data
#' @param psf_thresh Optional, a dataframe containing threshold specifications 
#'
#' @return A formatted dataframe of perceptible ranges for display on the user interface
#' @keywords internal
#'
#' @import dplyr
#' @importFrom rlang .data
#' 
generate_perceptible_ranges <- function(pf_data, psf_thresh = data.frame()){
  
  #If psf is provided, use that as the (sole) source of truth and skip subsequent processing
  if (nrow(psf_thresh) != 0){ 
    
    perc_ranges_df <- psf_thresh %>% select(!.data$label)
    
  } else { #If PSF is not provided, handle as below
    
    #Handle 'Largest Since' rows
    largest_since_perc_ranges <- process_largest_since(pf_data)
    if (nrow(largest_since_perc_ranges) > 0){
      warning('Perception thresholds have been automatically populated based on historic peak-flow data from NWIS. The user should carefully review these perception thresholds to verify they correctly represent historic peak streamflow.')
    }
    
    #Generate default perceptible ranges based on input data
    default_perc_ranges <- pf_data %>% 
      group_by(.data$site_no) %>% 
      summarise(start = min(.data$peak_WY), 
                end = max(.data$peak_WY),
                HistPeriod = as.integer(sum(!is.na(.data$peak_va)))) %>% 
      select(.data$site_no, .data$start, .data$end) %>%
      mutate(min = 0, 
             max = Qmax, 
             comment = "Default")
    
    
    #Update to reflect that 'largest since' may indicate an earlier start year 
    if(nrow(largest_since_perc_ranges) > 0) {
      
      for(site in unique(largest_since_perc_ranges$site_no)) {
        
        sel_largest_since <- largest_since_perc_ranges %>% filter(.data$site_no == site)
        sel_default <- default_perc_ranges %>% filter(.data$site_no == site)
        
        min_year_largest_since <- min(sel_largest_since$start)
        max_year_largest_since <- max(sel_largest_since$end)
        
        min_year_default <- min(sel_default$start)
        max_year_default <- max(sel_default$end)
        
        min_year <- min(min_year_largest_since, min_year_default)
        max_year <- max(max_year_largest_since, max_year_default)
        
        default_perc_ranges$start[(default_perc_ranges$site_no == site)] = min_year
        default_perc_ranges$end[default_perc_ranges$site_no == site] = max_year
        
      }
      
    }
    
    # Modify perc ranges defaults based on peak_cd values, where relevant. 
    remark_cds <- c("4", "8", "O", "3") #Remark codes requiring special handling. Code 7 handled separately
    
    remark_cd_perc_ranges <- pf_data %>% 
      filter(!is.na(.data$peak_cd)) %>%
      select(.data$site_no, .data$peak_WY, .data$peak_va, .data$peak_cd)
    
    if(nrow(remark_cd_perc_ranges) != 0) {
      remark_cd_perc_ranges <- remark_cd_perc_ranges %>% 
        mutate(min = case_when(grepl("4", as.character(.data$peak_cd)) == TRUE ~ .data$peak_va,
                               grepl("8", as.character(.data$peak_cd)) == TRUE ~ 0,
                               grepl("O", as.character(.data$peak_cd)) == TRUE ~ Qmax,
                               grepl("3", as.character(.data$peak_cd)) == TRUE ~ Qmax),
               max = case_when(grepl("4", as.character(.data$peak_cd)) == TRUE ~ Qmax,
                               grepl("8", as.character(.data$peak_cd)) == TRUE ~ .data$peak_va,
                               grepl("O", as.character(.data$peak_cd)) == TRUE ~ Qmax,
                               grepl("3", as.character(.data$peak_cd)) == TRUE ~ Qmax)) %>%
        add_cd_comments() %>%
        filter(!(is.na(.data$min) & is.na(.data$max))) %>% #remove cd7 rows since handled separately for perc ranges
        select(.data$site_no, start = .data$peak_WY, end = .data$peak_WY, .data$min, .data$max, .data$comment) 
      
    }
    
    #Handle code 7 perceptible ranges
    cd7_perc_ranges <- pf_data %>% process_cd_7()
    if (nrow(cd7_perc_ranges) > 0){
      warning('Perception thresholds have been automatically populated based on historic peak-flow data from NWIS. The user should carefully review these perception thresholds to verify they correctly represent historic peak streamflow.')
    }
    
    #Historic thresholds (cd7) should supercede 'largest since' thresholds: "Thresholds assigned by a historic peak should overwrite those from a  largest-since year." 
    #Identify all code 7 years 
    cd7_yrs <- c()
    
    if(nrow(cd7_perc_ranges) != 0){
      
      #List with vectors of code 7 years at all sites
      cd7_yrs <- mapply(seq, cd7_perc_ranges$start, cd7_perc_ranges$end, SIMPLIFY = FALSE) 
      
      #Create dataframe matching years to site numbers
      cd7s <- data.frame(site_no = rep(cd7_perc_ranges$site_no, sapply(cd7_yrs, length)), 
                 cd7_yrs = unlist(cd7_yrs))
        
    }
    
    #Identify all largest_since years that should drop off because they are superseded by a historic perception threshold
    if(nrow(largest_since_perc_ranges) != 0){
      for(j in 1:nrow(largest_since_perc_ranges)){
        
        sel_row <- largest_since_perc_ranges[j,]
        
        sel_site <- sel_row$site_no
        
        st_yr <- sel_row$start
        end_yr <- sel_row$end
        
        yr_seq <- seq(st_yr, end_yr)
        
        if(nrow(cd7_perc_ranges) == 0){next} #No processing needed no code 7 peaks
        
        #Get code 7 years for this site
        cd7s_site <- cd7s$cd7_yrs[cd7s$site_no == sel_site]
        
        overlap_yrs <- yr_seq[yr_seq %in% cd7s_site]
        
        if(length(overlap_yrs) == 0){next} #No processing needed if it doesn't overlap with a historic perception threshold
        
        #Figure out the start and end of the overlap
        st_olp <- min(overlap_yrs)
        end_olp <- max(overlap_yrs)
        
        if(st_olp == st_yr & end_olp == end_yr){ #If there's a complete overlap: 
          largest_since_perc_ranges[j,] <- NA
        } else { #If there's incomplete overlap: 
          st_out <- min(yr_seq[yr_seq %!in% overlap_yrs])
          end_out <- max(yr_seq[yr_seq %!in% overlap_yrs])
          
          #Revise the row accordingly
          largest_since_perc_ranges$start[j] <- st_out
          largest_since_perc_ranges$end[j] <- end_out
        }
      }
      
      #Remove the rows where largest since overlapped with (was superseded by) historic peak data
      largest_since_perc_ranges <- na.omit(largest_since_perc_ranges)
    }
    
    #Build: Default rows + remark code rows
    perc_ranges_df <- rbind(default_perc_ranges, largest_since_perc_ranges, remark_cd_perc_ranges, cd7_perc_ranges) %>%
      group_by(.data$site_no) %>% 
      arrange(start, desc(end), .by_group = TRUE) #SORT BY STATION, Beg YEAR AFTER DEFAULT
    
  }
  
  perc_ranges_df <- perc_ranges_df %>% 
                    ungroup() %>% 
                    unique()
  
  return(perc_ranges_df)
  
}

####   Function to generate data/flow intervals df from input data and, where provided, specs
#' Generate Data/Flow Intervals
#' @description
#' Generate data/flow intervals dataframe from input data and, 
#' where provided, specifications 
#' @param pf_data A dataframe of peak flow data
#' @param psf_peaks Optional, a dataframe containing peak specifications
#' @param psf_int Optional, a dataframe containing interval specifications
#' @param psf_thresh Optional, a dataframe containing threshold specifications
#' @param psf_opt Optional, a dataframe containing processing specifications
#'
#' @return
#' A list containing:
#' 1. A dataframe with data/flow intervals based on the peak flow data, and, where provided, specifications.
#' 2. A list of site IDs for which all data are urban/regulated peaks that have not been selected for inclusion. 
#' @keywords internal
#' 
#' @import dplyr
#' @importFrom rlang .data
#'
#' 
generate_data_flow_intervals <- function(pf_data, psf_peaks = data.frame(), psf_int = data.frame(), psf_thresh = data.frame(), psf_opt = data.frame()){
  
  #Need to get list of sites that only have urban/regulated peaks and don't 
  #have the use urban/regulated option set to TRUE
  
  urb_reg_sites <- unique(pf_data$site_no[grepl("6|C", as.character(pf_data$peak_cd))])
  unreg_sites <- unique(pf_data$site_no[!grepl("6|C", as.character(pf_data$peak_cd))])
  urb_reg_only_sites <- urb_reg_sites[urb_reg_sites %!in% unreg_sites]
  
  #Select sites that don't have Urb/Reg set to TRUE in psf
  urb_reg_only_sites <- urb_reg_only_sites[urb_reg_only_sites %!in% psf_opt$site_no[psf_opt$`Urb/Reg`]]
  
  
  
  #Filter only to sites that are not urb/reg only:
  pf_data <- pf_data %>% 
             filter(.data$site_no %!in% urb_reg_only_sites)
  
  #Columns year and comments come from input data. Remark codes, peak values, and intervals come from siteQT
  data_rows <- pf_data %>% 
               select(.data$site_no, .data$peak_WY, .data$peak_cd) %>%
               add_cd_comments() %>%
               select(.data$site_no, .data$peak_WY, .data$comment)
  
  #If psf is provided, ensure it is integrated with the data cols.
  #psf_int and psf_peaks contributions get handled through siteQT, but need to pull in comments.
  psf_rows <- data.frame(matrix(ncol = 3, nrow = 0, dimnames = list(NULL, c("site_no", "peak_WY", "comment"))))
  
  #pull in psf_int comments, where provided
  if (nrow(psf_int) != 0){
    psf_int_rows <- psf_int %>% 
      select(.data$site_no, .data$peak_WY, .data$comment)
  } else {
    psf_int_rows <- data.frame(matrix(ncol = 3, nrow = 0, dimnames = list(NULL, c("site_no", "peak_WY", "comment"))))
  }
  
  #pull in psf_peaks comments, where provided
  if (nrow(psf_peaks) != 0){
    psf_peaks_rows <- psf_peaks %>% 
                      select(.data$site_no, .data$peak_WY, .data$comment)
  } else {
    psf_peaks_rows <- data.frame(matrix(ncol = 3, nrow = 0, dimnames = list(NULL, c("site_no", "peak_WY", "comment"))))
  }
  
  #collect comments coming from the specifications
  psf_rows <- psf_rows %>% 
              rbind(psf_int_rows) %>% 
              rbind(psf_peaks_rows) %>% 
              unique()
  
  # #drop rows from data_rows that are covered by psf inputs to avoid duplication
  if(nrow(psf_rows) > 0){
    dup_rows <- inner_join(data_rows[,1:2], psf_rows[,1:2]) %>%
                mutate(unique_id = paste0(.data$site_no, "_", .data$peak_WY))
  
    data_rows <- data_rows %>%
                 mutate(unique_id = paste0(.data$site_no, "_", .data$peak_WY)) %>%
                 filter(.data$unique_id %!in% dup_rows$unique_id) %>%
                 select(!.data$unique_id)
  }
  
  #Combine comments from data, and, if provided, psf. 
  comment_cols <- merge(data_rows, psf_rows, all = TRUE) %>% 
                  arrange(.data$peak_WY) 
  
  #Columns upper bound and lower bound come from site QT output
  qt_df <- all_siteQT(pf_data = pf_data, psf_int = psf_int, psf_thresh = psf_thresh, psf_peaks = psf_peaks, psf_opt = psf_opt, keepNoInfo = TRUE)
  
  qt_cols <- qt_df %>% 
             select(.data$site_no, .data$peak_WY, .data$peak_va, remark_cd = .data$peak_cd, interval_low = .data$ql, interval_up = .data$qu) 
  
  
  #Combine columns from data and siteQT
  data_flow_int_df <- merge(qt_cols, comment_cols, by = c("site_no", "peak_WY"), all = TRUE) %>%
                      select(names(data_flow_int_df_init)) %>%
                      ungroup()

  # Handle Gage Height Only rows
    #If a year in the input data only has a gage height (gage_ht column) but no peak discharge value (peak_va), then in the data/intervals table on the second tab add a comment "Gage height only" to the comments column
    
    #Identify and format rows
    gho_data_rows <- pf_data %>% 
                     filter(is.na(.data$peak_va) & !is.na(.data$gage_ht))
    
    for(i in 1:nrow(gho_data_rows)){
      
      sel_row <- gho_data_rows[i,]
      sel_site <- sel_row$site_no
      sel_WY <- sel_row$peak_WY
      
      #Identify index of corresponding row to replace
      i_to_rep <- which(data_flow_int_df$site_no == sel_site & data_flow_int_df$peak_WY == sel_WY & (is.na(data_flow_int_df$comment) | data_flow_int_df$comment == "GAGE HEIGHT ONLY" | data_flow_int_df$comment == ""))
      
      data_flow_int_df[i_to_rep, "comment"] <- "GAGE HEIGHT ONLY"
      
    }
  
  #Remove rows that have neither a recorded peak nor comment nor code. 
  data_flow_int_df <- data_flow_int_df %>% 
                      filter(!is.na(.data$peak_va) | !is.na(.data$comment)) %>%
                      unique()
    
  list(data_flow_int_df = data_flow_int_df, urb_reg_only_sites = urb_reg_only_sites)
  
}

