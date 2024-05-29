# quiets concerns of R CMD check when variables appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(".", "Exposure_name", "Estimate", "Estimate_CI", "p_value_FDR", "Location_in_gene"))
}

#' List CpGs within DMRs
#'
#' @param meth_data A matrix containing methylation data
#' @param probes A scalar defining to how many CpGs within a DMR the result should be restricted to
#' @param annotated_EWAS_result A data.frame with the annotated EWAS result
#' @param path_to_combp A string defining the path where the *.regions-t.bed files (results of the comb-p procedure) are stored
#' @param path A string defining the path to save the .bed files
#' @param file_name A string defining file name without extension
#'
#' @return A data.frame with a CpGs list for each identified DMR
#' @export
#' @import dplyr
#' @importFrom here here
#' @importFrom stringr word
#' @importFrom data.table fread
#' @importFrom tibble rownames_to_column
#' @importFrom utils write.csv
#'

CpGsInDMRsSEPAGES <- function(meth_data,
                       probes,
                       annotated_EWAS_result,
                       path_to_combp,
                       path,
                       file_name) {

  # Create an empty data.frame where the results will be collected
  CpGs_DMR <- data.frame()

  # Process the DMR analysis results files (ending with *.regions-t.bed)
  for (DMR_file in list.files(path = here::here(path_to_combp),
                              pattern = "\\.regions-t.bed$")) {

    # Subset the exposure name
    expo <- stringr::word(DMR_file, 1:3, sep = "_") %>%
      paste0(collapse = "_")

    print(expo)

    # List the DMRs with corresponding chromosomal location
    regions_CpGs <- data.table::fread(here::here(path_to_combp, DMR_file)) %>%

      # Filter regions with at leat 5 probes
      dplyr::filter(n_probes >= probes)

    if (nrow(regions_CpGs) == 0) {
      next
    }

    regions_CpGs <- regions_CpGs %>%
      dplyr::mutate(Exposure = expo) %>%
      dplyr::select(Exposure, '#chrom', start, end)

    CpG_list <- data.frame()

    for (i in seq_len(nrow(regions_CpGs))) {

      # List all CpGs within each DMR and the corresponding estimate
      list_CpGs_within_DMR <- dplyr::filter(annotated_EWAS_result,
                                            Exposure == unique(regions_CpGs$Exposure),

                                            Chr == as.character(regions_CpGs$'#chrom')[i] &
                                              Position >= regions_CpGs$start[i] &
                                              Position <= regions_CpGs$end[i]) %>%
        dplyr::select(Exposure,
                      Exposure_name,
                      CpG:Estimate,
                      Estimate_CI,
                      raw_p_value,
                      p_value_FDR) %>%

        # Categorize the estimate
        dplyr::mutate(Direction = dplyr::case_when(Estimate > 0 ~ "+",
                                                   TRUE ~ "-"))

      CpG_list <- dplyr::bind_rows(CpG_list, list_CpGs_within_DMR) %>%
        dplyr::arrange(Gene)
    }

    CpGs_DMR <- dplyr::bind_rows(CpGs_DMR, CpG_list) %>%
      dplyr::arrange(Gene, Exposure)

  }

  # Calculate mean methylation values for each CpG
  mean_meth <- data.frame("mean_meth" = colMeans(meth_data, na.rm = TRUE)) %>%
    tibble::rownames_to_column(var = "CpG")

  # Add mean methylation level to the result
  CpGs_in_DMRs <- CpGs_DMR %>%
    merge(mean_meth, by = "CpG") %>%
    dplyr::select(CpG,
                  Exposure,
                  Exposure_name,
                  Chr:Location_in_gene,
                  mean_meth,
                  Estimate,
                  Estimate_CI,
                  raw_p_value,
                  p_value_FDR) %>%
    dplyr::arrange(Exposure,
                   CpG)

  # Save the result
  utils::write.csv(CpGs_in_DMRs, here::here(path, paste0(file_name, ".csv")), row.names = FALSE)

  return(CpGs_in_DMRs)
}
