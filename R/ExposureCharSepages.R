# quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("."))
}

# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c("mean_expo", "t2", "t3"))


#' Calculates exposure characteristics
#'
#' @param input_data A data.frame containing information on exposure concentrations and <LOD
#' @param lod_values A double defining LOD value for a compound
#' @param path A string defining the path to save the output
#' @param file_name A string defining specific part of the file name
#'
#' @return A data.frame with exposure descriptive statistics
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @importFrom naniar replace_with_na_all
#' @importFrom stats quantile
#' @importFrom utils write.csv
#'
ExposureCharSepages <- function(input_data,
                                lod_values,
                                path,
                                file_name = NULL) {
  
  input_data_long <- input_data %>% 
    
    # Transform to long format
    tidyr::pivot_longer(cols = -id,
                        names_to = c("Exposure", "Trim"),
                        values_to = "Conc",
                        names_pattern = "(.*)_(.*)")
  
  exposure_char <- input_data_long %>%
    
    # Add LOD values for each compound
    merge(lod_values, by = "Exposure") %>% 
    dplyr::group_by(Exposure, Trim) %>% 
    
    # Calculate percentage of exposure samples >LOD and basic statistics
    dplyr::summarise(n = sum(!is.na(Conc)),
                     
                     "pct_det" = round((sum(Conc >= LOD, na.rm = TRUE) * 100) / n, 1),
                     "p5" = round(stats::quantile(Conc, na.rm = TRUE, probs = 0.05), 6),
                     "median" = round(stats::quantile(Conc, na.rm = TRUE, probs = 0.5), 6),
                     "p95" = round(stats::quantile(Conc, na.rm = TRUE, probs = 0.95), 3)
    ) %>% 
    merge(lod_values, by = "Exposure") %>% 
    # dplyr::mutate(p5 = ifelse(p5 < LOD, "<LOD", p5),
    #               median = ifelse(median < LOD, "<LOD", median)) %>%
    tidyr::pivot_wider(id_cols = -n, 
                       names_from = Trim,
                       values_from = c("pct_det", "p5", "median", "p95")) %>%
    dplyr::select(Exposure, 
                  LOD, 
                  dplyr::contains("pct"),
                  dplyr::contains("_t2"),
                  dplyr::contains("_t3"))
  
  exposure_char_aver <- input_data_long %>%
    dplyr::group_by(id, Exposure) %>%
    
    # Calculate concentration averaged over T2 and T3
    dplyr::summarize(aver_expo_conc_Std = mean(Conc, na.rm = TRUE)) %>%
    dplyr::group_by(Exposure) %>%
    
    # Calculate basic statistics
    dplyr::summarise("p5_av" = round(stats::quantile(aver_expo_conc_Std, na.rm = TRUE, probs = 0.05), 2),
                     "median_av" = round(stats::quantile(aver_expo_conc_Std, na.rm = TRUE, probs = 0.5), 2),
                     "p95_av" = round(stats::quantile(aver_expo_conc_Std, na.rm = TRUE, probs = 0.95), 1))
  
  exposure_characteristics <- dplyr::left_join(exposure_char, exposure_char_aver, by = "Exposure") %>% 
    dplyr::mutate(Expo_family = ifelse(stringr::str_detect(Exposure, "BPA|BPS|BP3|BUPA|ETPA|MEPA|PRPA"), "Phenol", "Phthalate")) %>% 
                  # p5_av = ifelse(p5_av < LOD, "<LOD", p5_av),
                  # median_av = ifelse(median_av < LOD, "<LOD", median_av)) %>% 
    dplyr::select(Exposure, Expo_family, everything()) %>% 
    dplyr::arrange(Expo_family, Exposure)
  
  if (!is.null(file_name)) {
    write.csv(exposure_characteristics, here::here(path, file_name), row.names = FALSE)
  }
  
  return(exposure_characteristics)
}
