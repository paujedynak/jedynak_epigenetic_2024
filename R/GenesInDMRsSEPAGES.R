# quiets concerns of R CMD check when variables appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(".", "Exposure_name", "Estimate", "Estimate_CI", "p_value_FDR", "Location_in_gene"))
}

#' List genes within DMRs or any other dataset that contains genes names grouped within one cell
#'
#' @param data A dataframe containing a column Gene containing strings of genes separated by ';'
#' @param var_to_distinct A string or a vector specifying the names of grouping variables 
#'
#' @return A dataframe containing a list of genes (one gene per row)
#' @export

GenesInDMRsSEPAGES <- function(data, var_to_distinct) {
  
  split_genes <- data %>% 
    dplyr::mutate(
      Gene = stringr::str_split(Gene, ";")) %>% 
    tidyr::unnest_longer(Gene) %>% 
    dplyr::distinct_at(vars(!!dplyr::sym(var_to_distinct[1]), 
                      !!dplyr::sym(var_to_distinct[2]),
                                     !!dplyr::sym(var_to_distinct[3])), 
                    .keep_all = TRUE) %>% 
    stats::na.omit() %>% 
    dplyr::filter(Gene != "Unknown") %>% 
    dplyr::ungroup()
  
  return(split_genes)
  
}