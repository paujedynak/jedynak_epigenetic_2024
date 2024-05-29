# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c("ident", "mother_bmi", "mother_edu",
                         "delivery_mode", "parity", "ch_sex", "corr_lmp",
                         "date_delivery", "date_lmp", "mean_ddg_echo", 
                         "mo_age", "mo_he", "mo_par",
                         "mo_we_bepr", "mother_height", "mother_weight",
                         "po_datedel", "po_datelmp", "mother_edu",
                         "season_conc", "child_sex", "mother_age", "gest_age"))

#' Create covariates data
#'
#' @param cov_expo_data A .sas7bdat file containing exposure and covariate data
#' @param technical_covariates_data A .csv file containing technical covariates data
#' @param exposures A character vector defining original names of exposures (from t2 and t3)
#' @param cell_mix A data.frame containing cell proportions for each individual
#'
#' @return A data.frame containing covariates data
#' @export
#' @import dplyr
#' @importFrom lubridate month
#' @importFrom tidyselect all_of
#' @importFrom tidyr replace_na

PrepareDatasetsSepagesAllComp <- function(cov_expo_data,
                                          technical_covariates_data,
                                          cell_mix,
                                          exposures) {
  
  # Process phenol exposures
  expo_log <- cov_expo_data %>%
    
    # Select exposures
    dplyr::select(ident, tidyselect::all_of(exposures)) %>%
    dplyr::mutate(id = as.character(ident), .keep = "unused") %>% 
    dplyr::select(id, dplyr::everything()) %>% 
    
    # Log2 imputed standardized exposures
    dplyr::mutate(dplyr::across(dplyr::contains("_cor_"), 
                                .fns = list(log2 = ~log2(.)),
                                .names = "{fn}_{col}")) %>% 
    
    # Change categorized standardized exposure concentrations to factors
    # 1 == < LOD
    # 2 > LOD & < LOQ
    # 3 > LOQ
    dplyr::mutate_at(dplyr::vars(contains("_cat_")), factor) %>% 
    dplyr::mutate(dplyr::across(dplyr::contains(c("_cat_")), ~ factor(.x, 
                                                                      levels = c("1", "2", "3"),
                                                                      labels = c("<LOD", "LOD-LOQ", ">LOQ")))) %>% 
    
    # Change "Non-detected" (ND) in the non-standardized BPS, BUPA into 0 and change to numeric
    dplyr::mutate_at(dplyr::vars(dplyr::contains("_string_")),
                     function(x) stringr::str_replace(x, 
                                                      pattern = "ND", 
                                                      replacement = "0")) %>% 
    dplyr::mutate_at(dplyr::vars(dplyr::contains("_string_")), as.numeric) %>% 
    
    # Change non-detected records in the non-standardized BPS, BUPA into values: LOD / sqrt(2)
    dplyr::mutate_at(c("mo_BPS_total_string_t2", "mo_BPS_total_string_t3"), function(x) ifelse(x == 0, 0.1 / sqrt(2), x)) %>%
    dplyr::mutate_at(c("mo_BUPA_total_string_t2", "mo_BUPA_total_string_t3"), function(x) ifelse(x == 0, 0.07 / sqrt(2), x))
  
  # Clean colnames
  colnames(expo_log) <- gsub('i_cor_|mo_|total_|ms_', '', colnames(expo_log))
  colnames(expo_log) <- gsub("OXBE", "BP3", colnames(expo_log))
  
  # Add mean value of two standardized log2 exposure concentrations for T2 and T3 for each compound
  # with continuous outcome
  expo_log_av <- expo_log %>%
    
    # Select log2 exposures
    dplyr::select(id, dplyr::contains("log2_")) %>%
    
    # Transform to long format
    tidyr::pivot_longer(cols = -id, 
                        names_to = c("Expo", "Trim"),
                        values_to = "Conc_log2",
                        names_pattern = "(.*)_(.*)") %>%
    dplyr::group_by(id, Expo) %>%
    
    # Calculate mean standardized log2 exposure concentration value over trimester T2 and T3 
    # for each compound
    dplyr::summarize(mean_expo_log2 = mean(Conc_log2, na.rm = TRUE)) %>% 
    tidyr::pivot_wider(names_from = "Expo",
                       values_from = "mean_expo_log2",
                       names_prefix = "aver_")
  
  # Merge raw and averaged concentration values
  expo_raw_aver <- expo_log %>% 
    dplyr::left_join(expo_log_av, by = "id") %>% 
    
    # Add summarizing value of categorized exposures (BUPA, BPS) for T2 and T3
    
    # if t2 == <LOD or NA and t3 == <LOD or NA assign <LOD_T2&3 (<LOD_T2&3 for trimesters 2 and 3) 
    # if t2 != <LOD and t3 == <LOD or NA or t3 != <LOD and t2 == <LOD or NA assign >LOD_T2orT3 (> LOD for trimester 2 or 3)
    # else (if t2 and t3 != 1) assign >LOD_T2&T3 (>LOD for trimesters 2 and 3)
    dplyr::rowwise() %>% 
    dplyr::mutate(aver_BPS = dplyr::case_when((BPS_cat_t2 == "<LOD" | is.na(BPS_cat_t2)) & (BPS_cat_t3 == "<LOD" | is.na(BPS_cat_t3)) ~ "<LOD_T2&3",
                                              (BPS_cat_t2 != "<LOD" & (BPS_cat_t3 == "<LOD" | is.na(BPS_cat_t3))) | (BPS_cat_t3 != "<LOD" & (BPS_cat_t2 == "<LOD" | is.na(BPS_cat_t2))) ~ ">LOD_T2orT3",
                                              TRUE ~ ">LOD_T2&T3"),
                  aver_2cat_BPS = dplyr::case_when((BPS_cat_t2 == "<LOD" | is.na(BPS_cat_t2)) & (BPS_cat_t3 == "<LOD" | is.na(BPS_cat_t3)) ~ "<LOD_T2&3",
                                                    TRUE ~ ">LOD_T2orT3"),
                  aver_BUPA = dplyr::case_when((BUPA_cat_t2 == "<LOD" | is.na(BUPA_cat_t2)) & (BUPA_cat_t3 == "<LOD" | is.na(BUPA_cat_t3)) ~ "<LOD_T2&3",
                                               (BUPA_cat_t2 != "<LOD" & (BUPA_cat_t3 == "<LOD" | is.na(BUPA_cat_t3))) | (BUPA_cat_t3 != "<LOD" & (BUPA_cat_t2 == "<LOD" | is.na(BUPA_cat_t2))) ~ ">LOD_T2orT3",
                                               TRUE ~ ">LOD_T2&T3"), 
                  aver_2cat_BUPA = dplyr::case_when((BUPA_cat_t2 == "<LOD" | is.na(BUPA_cat_t2)) & (BUPA_cat_t3 == "<LOD" | is.na(BUPA_cat_t3)) ~ "<LOD_T2&3",
                                                     TRUE ~ ">LOD_T2orT3")) %>% 
    dplyr::mutate_at(c("aver_BPS", "aver_2cat_BPS", "aver_BUPA", "aver_2cat_BUPA"), as.factor)
  
  
  # Process covariates
  cov <- cov_expo_data %>%
    dplyr::transmute(
      id = as.character(ident),
      
      # Maternal age
      mother_age = as.numeric(mo_age),
      
      # Maternal active smoking
      maternal_smoke_bef_pregn = dplyr::case_when(mo_tob_avgr == 0 ~ 0,
                                                  mo_tob_avgr > 0 ~ 1),
      maternal_smoke_gest_yn = mo_tob_grstt1_yn,
      maternal_smoke_anytime_yn = mo_tob_gr_anyt_yn_n2,
      
      # Maternal education
      mother_edu =
        dplyr::case_when(
          mo_dipl %in% c(1, 2, 3) ~ "bac+2",
          mo_dipl == 4 ~ "bac+3,4",
          mo_dipl == 5 ~ "bac>5"
        ),
      mother_edu = factor(mother_edu, 
                          levels = c("bac+2", "bac+3,4", "bac>5")),
      
      # Parity
      parity = factor(mo_par,
                      
                      # Combine two last levels One and more
                      labels = c("Nulliparous", "One_and_more", "One_and_more")),
      
      # Maternal weight
      mother_weight = mo_we_bepr,
      
      # Maternal height
      mother_height = mo_he,
      
      # Child sex
      child_sex = factor(ch_sex,
                         labels = c("Male", "Female")),
      
      # Delivery date
      date_delivery = as.Date(po_datedel),
      
      # Echo based conception date (mean out of 3 Usounds)
      mean_ddg_echo = po_datestartpreg_us,
      
      # LMP date
      date_lmp = as.Date(po_datelmp)) %>%
    
    # Gestational duration based on LMP and corrected for US
    dplyr::rowwise() %>%
    dplyr::mutate(
      # Calculate difference in days between conception date based on the LMP date reported by the woman and 
      # the US-based conception date
      diff = as.numeric(mean_ddg_echo - (date_lmp + 14)),
      
      # Based on the date of the LMP or gestational duration assessed by the obstetrician if it differed
      # from the LMP-based estimate by more than 2 weeks
      corr_lmp = dplyr::case_when(abs(diff) > 14 ~ mean_ddg_echo,
                                  TRUE ~ date_lmp),
      
      # Calculate gestational duration based on the LMP reported by a woman
      # corrected with US data if differed by more than 2 weeks
      gest_age = date_delivery - corr_lmp,
      gest_age = as.numeric(gest_age / 7),
      
      # Calculate season of conception
      season_conc = case_when(
        lubridate::month(corr_lmp + 14) %in% c(1:3) ~ "Jan_March",
        lubridate::month(corr_lmp + 14) %in% c(4:6) ~ "April_June",
        lubridate::month(corr_lmp + 14) %in% c(7:9) ~ "July_Sept",
        lubridate::month(corr_lmp + 14) %in% c(10:12) ~ "Oct_Dec"),
      season_conc = factor(season_conc, levels = c("Jan_March", "April_June", "July_Sept", "Oct_Dec")),
      
      # Calculate pre-pregnancy BMI [weight (kg) / height (cm) / height (cm)] x 10,000
      mother_bmi = (mother_weight / mother_height / mother_height) * 10000,
      mother_bmi = dplyr::case_when(mother_bmi < 18.5 ~ "Underweight",
                                    mother_bmi >= 18.5 & mother_bmi < 25 ~ "Normal_weight",
                                    mother_bmi >= 25 & mother_bmi < 30 ~ "Overweight_obesity",
                                    mother_bmi >= 30 ~ "Overweight_obesity"),
      mother_bmi = as.factor(mother_bmi),
      mother_bmi = stats::relevel(mother_bmi, ref = "Normal_weight"),
      
      # Generate maternal active smoking variable
      maternal_smoke_sum = sum(maternal_smoke_bef_pregn, maternal_smoke_gest_yn, maternal_smoke_anytime_yn, na.rm = TRUE),
      maternal_smoke_sum = dplyr::case_when(is.na(maternal_smoke_bef_pregn) & is.na(maternal_smoke_gest_yn) &
                                              is.na(maternal_smoke_anytime_yn) ~ NA_real_,
                                            TRUE ~ maternal_smoke_sum),
      mother_act_smoke = case_when(maternal_smoke_sum == 0 ~ "Didn't_smoke",
                                   maternal_smoke_sum != 0 ~ "Smoked_before_OrAnd_in_pregn"), 
      mother_act_smoke = factor(mother_act_smoke)
    )
  
  cov <- select(cov,
                !matches("echo|date|maternal"),
                -mother_weight,
                -mother_height,
                -diff,
                -corr_lmp
  )
  
  # create covariates data
  
  # Merge exposures with covariates
  expo_cov <- merge(expo_raw_aver, cov, by = "id") # dim 395
  
  # Process technical confounders
  technical_conf <- technical_covariates_data %>% # dim 395   6
    
    # Select variables of interest. "Plaque" is the correct variable coding the plate number
    # Change chip, plate, and batch to factors
    dplyr::transmute(id = as.character(ident),
                     chip = factor(chip),
                     plate = factor(plate),
                     batch = as.factor(batch)) %>% 
    
    # Add cellmix
    dplyr::left_join(cell_mix, by = "id")
  
  # Merge exposures-covariates with technical factors
  expo_cov_tech_conf <- dplyr::left_join(expo_cov, technical_conf, by = "id") # dim 395  37
  
  return(expo_cov_tech_conf)
}

#' Impute missing data in categorical covariates (with mode)
#'
#' @param data_w_missing A data.frame with covariates containing missing values
#'
#' @return A data.frame with covariates with missing values imputed
#' @export

ImputeMissingCov <- function(data_w_missing) {
  
  # data_w_missing %>% dplyr::mutate_if(is.numeric, funs(replace(.,is.na(.), mean(., na.rm = TRUE)))) %>%
  #   dplyr::mutate_if(is.factor, funs(replace(.,is.na(.), Mode(na.omit(.)))))
  
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  data_complete <- data_w_missing %>%
    dplyr::mutate_if(is.factor, ~tidyr::replace_na(., Mode(na.omit(.))))
  
  return(data_complete)
}
