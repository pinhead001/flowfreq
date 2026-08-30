
#' Write PSF
#' @description
#' This function generates the specification file output. 
#' 
#' @param psf_opt dataframe of processing specifications for sites returned from 
#' function 'readPSF'
#' @param psf_thresh dataframe of threshold specifications for included sites 
#' @param psf_int dataframe of interval specifications for included sites 
#' @param psf_peaks dataframe of peak specifications for included sites 
#' @param basename basename for output files as specified in user interface
#' @param full_o_path full output path where specification file will be saved, 
#' based on inputs in user interface
#' @param o_plot_format output plot format based on input from user interface
#' @param o_plot_pos output plot position based on input from user interface
#' @param o_conf symmetric confidence interval based on input from user 
#' interface
#' @param o_extend_analysis extended analysis option based on checkbox in user
#' interface
#' @param tab_file optional, full filepath of tab-delimited file linked in
#' specification file
#' @param site_file optional, full filepath of site information file linked in
#' specification file
#'
#' @return This function returns the text that will ultimately be written into 
#' the .psf file 
#' @keywords internal
#' 
#' @importFrom  R.utils getRelativePath
#'
#' 
writePSF <- function(psf_opt, 
                     psf_thresh, 
                     psf_int, 
                     psf_peaks, 
                     basename, 
                     full_o_path, 
                     o_plot_format, 
                     o_plot_pos, 
                     o_conf, 
                     o_extend_analysis,
                     tab_file = NULL, 
                     site_file = NULL){

  # Get Attribute values
  peakFormat <- attributes(psf_opt)$peakFormat
  peakFile <- attributes(psf_opt)$peakFile
  siteFile <- attributes(psf_opt)$siteinfoFile
  outFile <- attributes(psf_opt)$outFile
  PlotFormat <- attributes(psf_opt)$PlotFormat
  

  df <- data.frame(matrix(nrow = 0, ncol = 1))
  
  # Spec File Header Fields 
  if(!is.na(peakFile) && !is.null(peakFile)) {
    df <- rbind(df, paste0("I ", peakFormat, " ", peakFile))
    }
  else if(!is.na(tab_file) && !is.null(tab_file)){
    
    rel_tab <- R.utils::getRelativePath(tab_file, relativeTo=full_o_path) # Find relative paths

    # Confirm file name is not "NA" (chack for last three characters in filepath)
    if(substr(rel_tab, nchar(rel_tab)-2, nchar(rel_tab)) != "/NA") {
      
      rel_tab <- sub("^\\.{2}", ".", rel_tab) # If there are two dots at beginning of relative path, change it to one dot
      rel_tab <- sub("^\\.\\/", "", rel_tab)  # If "./" is at beginning of relative path, remove it
      
      df <- rbind(df, paste0("I RDB ", rel_tab))
      }
    }
  
  if(!is.na(siteFile) && !is.null(siteFile)){
    df <- rbind(df, paste0("I INFO ", siteFile))
  }
  else if(!is.na(site_file) && !is.null(site_file)) {
    
    rel_site <- R.utils::getRelativePath(site_file, relativeTo=full_o_path) # Find relative paths
    
    # Confirm file name is not "NA" (chack for last three characters in filepath)
    if(substr(rel_site, nchar(rel_site)-2, nchar(rel_site)) != "/NA") {
      
      rel_site <- sub("^\\.{2}", ".", rel_site)   # If there are two dots at beginning of relative path, change it to one dot
      rel_site <- sub("^\\.\\/", "", rel_site)  # If "./" is at beginning of relative path, remove it
      
      df <- rbind(df, paste0("I INFO ", rel_site))
      }
    }

  
  df <- rbind(df, paste0("O File ", basename, ".txt"))
  df <- rbind(df, paste0("O Plot Format ", o_plot_format))
  df <- rbind(df, paste0("O Plot Position ", o_plot_pos))
  df <- rbind(df, paste0("O ConfInterval ", o_conf))
  
  if(!is.na(o_extend_analysis)) {
    if(o_extend_analysis == TRUE){
      df <- rbind(df, paste0("O EXTENDED YES"))}}
  
  
  # Produce station specific information
  for (i in 1:nrow(psf_opt)){
    
    df <- rbind(df, paste0("Station ", psf_opt$site_no[i]))
    
    #USGS was not responsible for the nested for loop
    for (j in 1:nrow(psf_thresh)) {
      
      if(psf_opt$site_no[i] == psf_thresh$site_no[j]){
        df <- rbind(df, paste0("     ", psf_thresh$label[j], " ", psf_thresh$start[j], " ", psf_thresh$end[j], " ", psf_thresh$min[j], " ", psf_thresh$max[j], " ", psf_thresh$comment[j]))
      }
    }
    
    #Get intervals for site
    psf_int_site <- psf_int[psf_int$site_no == psf_opt$site_no[i],]
    
    if(nrow(psf_int_site) > 0){
      introws <- t(rbind(paste0("     ", psf_int_site$label, " ", psf_int_site$peak_WY, " ", " ", psf_int_site$interval_low, " ", psf_int_site$interval_up, " ", psf_int_site$comment)))
      colnames(introws) <- colnames(df)
      df <- rbind(df, introws)
    }
    

    
    if(!is.na(psf_opt$SkewOpt[i])) {
      df <- rbind(df, paste0("     SkewOpt ", psf_opt$SkewOpt[i]))}
    if(!is.na(psf_opt$GenSkew[i])) {
      df <- rbind(df, paste0("     GenSkew ", psf_opt$GenSkew[i]))}
    if(!is.na(psf_opt$SkewSE[i])) {
      df <- rbind(df, paste0("     SkewSE ", psf_opt$SkewSE[i]))}
    if(!is.na(psf_opt$LOType[i])) {
      df <- rbind(df, paste0("     LOType ", psf_opt$LOType[i]))}
    if(!is.na(psf_opt$LoThresh[i])) {
      df <- rbind(df, paste0("     LoThresh ", psf_opt$LoThresh[i]))}
    if(!is.na(psf_opt$`Urb/Reg`[i])) {
      if(psf_opt$`Urb/Reg`[i] == TRUE){
      df <- rbind(df, paste0("     Urb/Reg YES"))}}

  }
  
  # only write the table if the app is being run locally
  write.table(x = df, file = full_o_path, row.names = F, col.names = F, quote = F,sep = '\t')
  
  return(df)
  
} 





#' Write NWIS Peak
#' 
#' Writes out the data returned by an NWIS query. Using the provided site 
#' numbers, reads peak flow from NWISweb. Data is retrieved from
#' https://waterdata.usgs.gov/nwis. Returns a tab delimited 
#' table, written to the file path, that includes both the NWIS peak flow
#' data and comments
#' 
#' @param site_no_list  USGS site number (or multiple sites). This is usually an
#'  8 digit number.
#' @param peak_file_path  a character string naming a file path to which the output 
#' peak data file should be written   
#' @param info_file_path  a character string naming a file path to which the output 
#' site infofile should be written   
#'
#' @keywords internal
#'
#' 
write_NWIS_peak <- function (site_no_list, peak_file_path, info_file_path){
  
  df <- data.frame(matrix(nrow = 0, ncol = 1))
  
  nwis_data <- readNWISpeak(site_no_list, asDateTime = FALSE, convertType = FALSE)
  cmnt <- attributes(nwis_data)$comment
  site_info <- attributes(nwis_data)$siteInfo
  
  col_def <- c("5s","15s","10d","6s","8s","33s","8s","27s","4s","10d","6s","8s","27s")
  col_names <- names(nwis_data)
  
  for (i in 1:length(cmnt)){
    df <- rbind(df, cmnt[i])
  }
  
  col_name_row <- paste(col_names, collapse = "\t")
  df <- rbind(df, col_name_row)
  
  col_def_row <- paste(col_def, collapse = "\t")
  df <- rbind(df, col_def_row)
  
  # Match number of columns in nwis_data
  n_blank_cols <- 12
  for (i in 1:n_blank_cols) {
    new_col_name <- paste("BlankColumn", i, sep = "")
    df[[new_col_name]] <- NA
  }
  
  # Make column names consistent
  new_col_names <- letters[1:13]
  colnames(df) <- new_col_names
  colnames(nwis_data) <- new_col_names
  
  df_bind <- rbind(df, nwis_data)
  
  write.table(df_bind, file = peak_file_path,
              na = "", row.names = F, col.names = F, quote = F,sep = '\t')
  
  #######   Now write out site metadata   ######
  
  site_info_df <- data.frame(matrix(nrow = 0, ncol = 1))
  site_info_col_def <- c("5s", "15s", "50s", "7s", "16s", "16s", "16s", "16s", "1s", "1s", "10s", "10s", "3s", "2s", "3s", "2s", "23s", "20s", "7s", "8s", "1s", "3s", "10s", "16s", "2s", "1s", "30s", "8s", "8s", "8s", "8s", "6s", "1s", "1s", "30s", "10s", "8s", "1s", "8s", "8s", "1s", "12s")
  site_info_col_names <- names(site_info)
  
  site_info_col_name_row <- paste(site_info_col_names, collapse = "\t")
  site_info_df <- rbind(site_info_df, site_info_col_name_row)
  
  site_info_col_def_row <- paste(site_info_col_def, collapse = "\t")
  site_info_df <- rbind(site_info_df, site_info_col_def_row)
  
  # Match number of columns in nwis_data
  site_info_n_blank_cols <- 41
  for (i in 1:site_info_n_blank_cols) {
    new_col_name <- paste("BlankColumn", i, sep = "")
    site_info_df[[new_col_name]] <- NA
  }
  
  # Make column names consistent
  site_info_new_col_names <- 1:42
  colnames(site_info_df) <- site_info_new_col_names
  colnames(site_info) <- site_info_new_col_names
  
  site_info_df_bind <- rbind(site_info_df, site_info)
  
  write.table(site_info_df_bind, file = info_file_path,
              na = "", row.names = F, col.names = F, quote = F,sep = '\t')
  
  return(list(df_bind, site_info_df_bind))
  
}




#' Write Log 
#'
#' Writes a log file of any errors and warnings generated during PeakFQ runs 
#' 
#' @param x any code passed to the function for which warnings and errors should
#' be recorded
#'
#' @keywords internal
#'
#' 
write_log <- function (x){
  track <- list()
  
  tryCatch(withCallingHandlers(
    x,
    
    # Handle the warnings.
    warning = function(w) {
      # warns <<- c(warns, list(w))
      track <<- c(track, list(w))
      invokeRestart("muffleWarning")
    }),
    
    # Handle the errors.
    error = function(e) {
      # errs <<- c(errs, list(e))
      track <<- c(track, list(e))
    }
  )
  
  # Remove ANSI escape codes to make messages more readable
  # Have to add extra double slash (\\) before each special character
  codes_slash <- c("\\\033\\[1m", "\\\033\\[33m", "\\\033\\[39m", "\\\033\\[22m", "\\\033\\[38\\;5\\;232m", "\\\n")
  
  if (length(track) > 0){
    for (i in 1:length(track)){
      
      message <- as.character(track[[i]])
      for (j in 1:length(codes_slash)){
        message <- gsub(codes_slash[j], "", message, ignore.case = TRUE)
        
      }
      
      #track_log <<- rbind(track_log, message)
    }
    
  }
  
}

