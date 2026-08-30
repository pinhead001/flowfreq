#############       RHOT renderer functions

#Renderer for read-only colums - show as light grey
read_only_style <- "function(instance, td, row, col, prop, value, cellProperties) {
                      Handsontable.renderers.TextRenderer.apply(this, arguments);
                        td.style.background = '#ECECEC';
                      }"


#Renderer for incomplete comment cells - show as light red
incomplete_style <- "function (instance, td, row, col, prop, value, cellProperties) {
                    Handsontable.renderers.TextRenderer.apply(this, arguments);

                    if (instance.params) {
                       hrows = instance.params.row_highlight
                       hrows = hrows instanceof Array ? hrows : [hrows]
                       hcols = instance.params.col_highlight
                       hcols = hcols instanceof Array ? hcols : [hcols]

                       if (hrows.includes(row) && hcols.includes(col) && value != 0 && (!value || value == '')) {
                         td.style.background = '#ff8888'
                       }

                     }

                    }"

safeHtmlRenderer <- "function (instance, td, row, col, prop, value, cellProperties) {
  td.innerHTML = value;
}"

#renderRHandsontable requires an input whose output is an rhandsontable. This allows for a dataframe input with not function nesting.
#' Render R HandsOnTable
#' @description
#' Function that allows the user to pass a dataframe that displays as an RHandsOnTable 
#' without function nesting
#' @param df Any dataframe 
#'
#' @return UI component - RHandsOnTable 
#' @keywords internal
#'
#' 
render_RHOT = function(df){renderRHandsontable({rhandsontable(df)})}  


####   config_station_spec_RHOT - generate editable table for specification inputs in first tab
#' Configure Station Specification Inputs
#' @description
#' Generate a formatted and editable table for input of processing 
#' specifications in the UI
#' @param station_specs_df A formatted dataframe of processing specifications, 
#' created by function 'generate_station_specifications'
#'
#' @return RHandsonTable object
#' @keywords internal
#'
#' @import dplyr
#' @import rhandsontable
#' @importFrom rlang .data
#' 
config_station_spec_RHOT <- function(station_specs_df){
  
  col_widths <- c(100, 100, 100, 125, 90, 75, 100, 85, 85, 80, 80, 80, 90, 85, 85)

  # Build the data frame
  df_for_rhot <- station_specs_df %>% mutate('Station ID' = as.character(.data$site_no),
                                             "Input Data Start Year" = as.integer(.data$BegYear), #Based on 8/21 input from  Seth - "populated from first... year of input data for the site"
                                             "Input Data End Year" = as.integer(.data$EndYear), #Based on 8/21 input from  Seth - "populated from... last year of input data for the site"
                                             "Input Data Record Length" = as.integer(.data$HistPeriod), #Based on 8/21 input from  Seth - "populated from the record length computed as the number of rows of data for the site (with valid peak discharge values, some may only have peak gage height, which we don't want to count)." 
                                             "Skew Option" = factor(.data$SkewOpt, levels = c("Station", "Weighted", "Regional")),
                                             "Use Map Skew" = as.logical(.data$MapSkew),
                                             "Map Skew Source" = as.character(.data$MapSkewSource),
                                             "Regional Skew" = as.numeric(.data$GenSkew),
                                             "Reg Skew Std Error" = as.numeric(.data$SkewSE),
                                             "Mean Sqr Err" = round(.data$`Reg Skew Std Error`^2, 3), 
                                             "PILF (LO) Test" = factor(.data$LOType, levels = c("MGBT", "FIXED", "NONE")),
                                             "PILF (LO) Threshold" = as.numeric(ifelse(.data$LOType == "MGBT", 0, .data$LoThresh)), 
                                             "Urban/ Reg Peaks" = as.logical(.data$`Urb/Reg`),
                                             "Latitude" = as.numeric(.data$Latitude),
                                             "Longitude" = as.numeric(.data$Longitude))
  df_for_rhot <- df_for_rhot[,c("Station ID","Input Data Start Year", "Input Data End Year", "Input Data Record Length","Skew Option","Use Map Skew", "Map Skew Source", "Regional Skew", "Reg Skew Std Error",
                                "Mean Sqr Err", "PILF (LO) Test", "PILF (LO) Threshold", "Urban/ Reg Peaks",
                                "Latitude", "Longitude")]
  
  
  #If 'Weighted' or 'Regional' is selected, highlight Skew and Std Error cells missing values
  regional_weighted_rows <- which((df_for_rhot$`Skew Option` == "Weighted" | df_for_rhot$`Skew Option` == "Regional") & (is.na(df_for_rhot$`Regional Skew`) | is.na(df_for_rhot$`Reg Skew Std Error`)))
  row_highlight = regional_weighted_rows - 1 #Zero-based indexing
  col_highlight = which(colnames(df_for_rhot) %in% c("Regional Skew", "Reg Skew Std Error")) - 1 #Zero-based indexing
  
  
  #Render RHOT of data frame, with input controls 
  rhot <- df_for_rhot %>% rhandsontable(height = 300, rowHeaders = FALSE, stretchH = "none", row_highlight = row_highlight, col_highlight = col_highlight) %>%
    hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE) %>%
    hot_col("Station ID", readOnly = TRUE, renderer = read_only_style, halign = "htCenter") %>%
    hot_col("Input Data Start Year", readOnly = TRUE, renderer = read_only_style, halign = "htCenter") %>%
    hot_col("Input Data End Year", readOnly = TRUE, renderer = read_only_style, halign = "htCenter") %>%
    hot_col("Input Data Record Length", readOnly = TRUE, renderer = read_only_style, halign = "htCenter") %>%
    hot_col("Use Map Skew", halign = "htCenter") %>%
    hot_col("Map Skew Source", readOnly = TRUE, renderer = safeHtmlRenderer, halign = "htCenter") %>%
    hot_col("Regional Skew", format = "0.000", halign = "htCenter", renderer = incomplete_style) %>% 
    hot_col("Reg Skew Std Error", format = "0.000", halign = "htCenter", renderer = incomplete_style) %>% 
    hot_col("Mean Sqr Err", readOnly = TRUE, renderer = read_only_style, format = "0.000", halign = "htCenter") %>%
    hot_col("PILF (LO) Threshold", format = "0.00", halign = "htCenter") %>% 
    hot_col("Latitude", readOnly = TRUE, renderer = read_only_style, halign = "htCenter") %>%
    hot_col("Longitude", readOnly = TRUE, renderer = read_only_style, halign = "htCenter") %>%
    hot_col("Urban/ Reg Peaks", halign = "htCenter") %>%
    hot_cols(colWidths = col_widths)
  
  #Set PILF (LO) Threshold to readOnly = T iff the PLIF (LO) Test is set to MGBT for that row
  mgbt_rows <- which(df_for_rhot$"PILF (LO) Test"=="MGBT")
  
  for (r in mgbt_rows){
    rhot <-  rhot %>% hot_cell(r, "PILF (LO) Threshold", readOnly = TRUE, 
                               comment = "Multiple Grubbs-Beck Test sets a default threshold of 0. Select 'FIXED' to enter a different value.") %>%
      hot_cols(colWidths = col_widths)
  }
  
  #Set editablility of skew  inputs based on selected skew option
  
  #If 'Station' is selected, map skew, regional skew, and reg skew std. err. should be FALSE/NA and read only
  station_rows <- which(df_for_rhot$`Skew Option` == "Station")
  
  for (r in station_rows){
    rhot <-  rhot %>% hot_cell(r, "Use Map Skew", readOnly = TRUE,   #! set to false
                               comment = "Map Skew cannot be used when Station Skew is selected.") %>% 
      hot_cell(r, "Map Skew Source", readOnly = TRUE,   #! set to false
               comment = "Map Skew source is not needed when Station skew is selected.") %>% 
      hot_cell(r, "Regional Skew", readOnly = TRUE, 
               comment = "Regional skew is not needed when Station skew is selected.") %>% 
      hot_cell(r, "Reg Skew Std Error", readOnly = TRUE, 
               comment = "Regional skew standard error is not needed when Station skew is selected.") %>% 
      hot_cell(r, "Mean Sqr Err", readOnly = TRUE, 
               comment = "Mean square error is not needed when Station skew is selected.") %>%
      hot_cols(colWidths = col_widths)
  }
  
  
  
  
  #If 'Use map skew' is checked, regional skew and reg skew std. error will be pulled from map, and set to read only
  map_skew_rows <- which(df_for_rhot$`Use Map Skew` == TRUE)
  
  for (r in map_skew_rows){
    rhot <-  rhot %>% hot_cell(r, "Regional Skew", readOnly = TRUE, 
                               comment = "Regional skew is pulled from map when 'Map Skew' is selected.") %>% 
      hot_cell(r, "Reg Skew Std Error", readOnly = TRUE, 
               comment = "Regional skew standard error is pulled from map when 'Map Skew' is selected.") %>% 
      hot_cell(r, "Mean Sqr Err", readOnly = TRUE) %>%
      hot_cell(r, "Map Skew Source", readOnly = TRUE, 
               comment = "Map Skew Source is pulled when 'Map Skew' is selected.") %>%
      hot_cols(colWidths = col_widths)
  }
  
  # Return the R Handsontable
  rhot 
  
} 

####   Perceptible ranges RHOT - generate editable table for perceptible range inputs in second tab
#' Configure Perceptible Range Inputs
#' @description
#' Generate a formatted and editable table for perceptible range inputs 
#' in the UI
#' @param perc_ranges_df A formatted dataframe of perceptible ranges, 
#' created by function'generate_perceptible_ranges'
#'
#' @return RHandsonTable object
#' @keywords internal
#' 
#' @import dplyr
#' @import rhandsontable
#' 
config_perc_ranges_RHOT <- function(perc_ranges_df){ 
  
  
  #Filtered in the argument to retain ability to reactively reference input$station. Rename for clarity.
  station_perc_ranges_df <- perc_ranges_df
  #Ensure a blank row at end
  
  #if every field in the last row is not NA, add another blank row beneath. This appears to be the logic used in the original
  if(sum(!is.na(last(station_perc_ranges_df))) == 6){ #
    station_perc_ranges_df[nrow(station_perc_ranges_df)+1,] <- NA
  }
  
  # Build the data frame
  df_for_rhot <- station_perc_ranges_df %>%
    mutate('Start Year' = as.integer(start),
           'End Year' = as.integer(end),
           'Lower Bound' = ifelse(min == Qmax, "inf", #Show/allow 'Inf' as equal to Qmax (global env var).
                                  ifelse(min == Qmin, 0, #Show/allow Qmin as 0
                                         as.numeric(min))), 
           'Upper Bound' = ifelse(max == Qmax, "inf", #Show/allow 'Inf' as equal to Qmax (global env var).
                                  ifelse(max == Qmin, 0, #Show/allow Qmin as 0
                                         as.numeric(max))), 
           'Comment (Required)' = as.character(comment)) %>% #! Comments need to be required. 
    select('Start Year', 'End Year', 'Lower Bound', 'Upper Bound', 'Comment (Required)')
  
  
  #Identify rows that have at least one entry - https://stackoverflow.com/questions/58142090/color-a-whole-row-in-rhandsontable-based-on-a-string-value-in-one-column
  data_rows <- which(rowSums(is.na(df_for_rhot)) < 5) 
  #uncommented_rows <- which((is.na(df_for_rhot$'Comment (Required)')|nchar(df_for_rhot$'Comment (Required)') == 0)) 
  #uncommented_data_rows <- uncommented_rows[uncommented_rows %in% data_rows]
  
  
  
  row_highlight = data_rows - 1 #Only highlight rows with some data entered
  col_highlight = 5 - 1
  
  #Render RHOT of data frame, with input controls 
  rhot <- df_for_rhot %>%
    rhandsontable(height = 200, rowHeaders = FALSE, row_highlight = row_highlight, col_highlight = col_highlight) %>% 
    hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE) %>%
    hot_cols(colWidths = c(60, 60, 80, 80, 120), 
             renderer = incomplete_style) 
  
  rhot
}

###    Data/Flow Intervals RHOT - generate editable table for Data/Flow interval inputs in second tab
#' Configure Data/Flow Interval Inputs
#' @description
#' Generate a formatted and editable table for Data/Flow interval 
#' inputs in the UI
#' @param data_flow_int_df A formatted dataframe of Data/Flow intervals, 
#' created by function 'generate_data_flow_intervals'
#'
#' @return RHandsonTable object
#' @keywords internal
#' 
#' @import dplyr
#' @import rhandsontable
#' @importFrom rlang .data
#'
#' 
config_data_flow_int_RHOT <- function(data_flow_int_df) {
  
  station_data_flow_int_df <- data_flow_int_df #Filtered in the argument to retain ability to reactively reference input$station. Rename for clarity. 
  
  if(nrow(station_data_flow_int_df) > 0) { 
    #Ensure a blank row at end
    #if every field in the last row is not NA, or every field except Peak is not NA, add another blank row beneath. This appears to be the logic used in the original
    if(sum(!is.na(last(station_data_flow_int_df)$peak_WY), !is.na(last(station_data_flow_int_df)$interval_low), !is.na(last(station_data_flow_int_df)$interval_up), !is.na(last(station_data_flow_int_df)$comment)) == 4) {
      station_data_flow_int_df[nrow(station_data_flow_int_df)+1,] <- NA
    }
  }
  
  
  # Build the data frame
  df_for_rhot <- station_data_flow_int_df %>%
    mutate('Year' = as.integer(.data$peak_WY),
           'Recorded Peak' = as.character(.data$peak_va),
           'Peak Codes' = as.character(.data$remark_cd),
           'Lower Bound' = as.character(ifelse(as.numeric(.data$interval_low) == Qmax, "inf", #Show/allow 'Inf' as equal to Qmax (global env var).
                                               ifelse(as.numeric(.data$interval_low) == Qmin, 0,  #Show/allow Qmin as 0
                                                      .data$interval_low))),
           'Upper Bound' = as.character(ifelse(as.numeric(.data$interval_up) == Qmax, "inf", #Show/allow 'Inf' as equal to Qmax (global env var).
                                               ifelse(as.numeric(.data$interval_up) == Qmin, 0, #Show/allow Qmin as 0
                                                      .data$interval_up))),
           'Comment (Required)' = as.character(.data$comment)) %>%
    select('Year', 'Recorded Peak', 'Peak Codes', 'Lower Bound', 'Upper Bound', 'Comment (Required)')
  
  #Identify rows for highlighting
  #Identify read-only columns, to display in grey - recorded peak & Peak Codes
  read_only_cols <- c(2, 3) - 1 #0-based indexing
  
  #Identify rows that have at least one user entry but are incomplete/lacking a comment. Shade comment cell. - https://stackoverflow.com/questions/58142090/color-a-whole-row-in-rhandsontable-based-on-a-string-value-in-one-column
  started_rows <- which(rowSums(is.na(df_for_rhot)) < 6) %>% unique() #Rows that are not fully NA
  added_rows <- which(is.na(df_for_rhot$`Recorded Peak`) & is.na(df_for_rhot$`Peak Codes`)) #proxy. If missing both peak and remark code, assume user entered. 
  user_added_rows <- added_rows[added_rows %in% started_rows] #Rows where the user has begun to enter data 
  
  uncommented_rows <- which((is.na(df_for_rhot$'Comment (Required)')|nchar(df_for_rhot$'Comment (Required)') == 0)) 
  
  uncommented_data_rows <- uncommented_rows[uncommented_rows %in% user_added_rows]
  
  row_highlight = uncommented_data_rows - 1 #0-based indexing
  col_highlight = 6 - 1 #0-based indexing
  
  #Identify rows with a code O or 3. Highlight entire row - these are not used.
  cd_O_3_rows <- which(grepl("3|O", df_for_rhot$'Peak Codes'))
  
  unused_rows <- cd_O_3_rows - 1 #0-based indexing
  
  # Identify row where comment is needed due to change to existing. 
  # Proxy (since otherwise every click on table may trigger): Require comment if recorded peak != lower bound != upper bound. Pre-programmed 'TRUES' due to Peak Codes should already have comments
  comment_needed_rows <- which(((df_for_rhot$`Recorded Peak` != df_for_rhot$`Lower Bound`) | (df_for_rhot$`Recorded Peak` != df_for_rhot$`Upper Bound`) | (df_for_rhot$`Lower Bound` != df_for_rhot$`Upper Bound`)) & #if ! (peak == lower bound == upper bound)
                                 (df_for_rhot$`Comment (Required)` == "")) #not yet commented
  
  add_comment_rows <- comment_needed_rows - 1 #0-based indexing
  
  row_highlight = row_highlight %>% append(add_comment_rows)
  
  #Define the renderer for data/flow intervals table
  dfi_style <- "function (instance, td, row, col, prop, value, cellProperties) {
                    Handsontable.renderers.TextRenderer.apply(this, arguments);

                    if (instance.params) {
                       rocols = instance.params.read_only_cols
                       rocols = rocols instanceof Array ? rocols : [rocols]
                       hrows = instance.params.row_highlight
                       hrows = hrows instanceof Array ? hrows : [hrows]
                       hcols = instance.params.col_highlight
                       hcols = hcols instanceof Array ? hcols : [hcols]
                       urows = instance.params.unused_rows
                       urows = urows instanceof Array ? urows : [urows]

                       if(rocols.includes(col)){
                         td.style.background = '#ECECEC';
                       }

                       if (hrows.includes(row) && hcols.includes(col)) {
                         td.style.background = '#ff8888'
                       }
                       
                       if (urows.includes(row)){
                         td.style.background = '#ff8888'
                       }

                     }

                    }"  
  
  #Render RHOT of data frame, with input controls 
  rhot <- df_for_rhot %>% 
    rhandsontable(height = 200, rowHeaders = FALSE,
                  row_highlight = row_highlight, col_highlight = col_highlight, #highlight rows with incomplete comments
                  unused_rows = unused_rows, #highlight rows not used in analysis (code O or 3)
                  read_only_cols = read_only_cols) %>% #shade read=only columns gray
    hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE) %>%
    hot_cols(colWidths = c(60, 80, 80, 80, 80, 120),
             renderer = dfi_style, 
             columnSorting = TRUE) %>%
    hot_col('Recorded Peak', readOnly = TRUE) %>%
    hot_col('Peak Codes', readOnly = TRUE)     
  
  rhot
}

