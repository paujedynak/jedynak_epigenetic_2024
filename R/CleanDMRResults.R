# quiets concerns of R CMD check when variables appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(".", "Exposure", "Gene", "Hg19", "start", "chrom", "end", "n_probes",
                           "raw_p_value", "Chr", "CpG", "chrom", "start", "end"))
}

#' Organizes DMR result in a neat form
#'
#' @param annotated_EWAS_result A data.frame with the annotated EWAS result
#' @param path_to_combp A string defining the path where the *.regions-t.bed files (results of the comb-p procedure) are stored
#' @param expo_levels A string defining names of the exposures in the dataset
#' @param expo_labels A string defining short neat labels for exposures
#' @param annotation_object Annotation object (Illumina)
#' @param sepages A logical to walk around different naming patterns between EDEN and SEPAGES
#'
#' @return A data frame with a summarised result of the DMR analysis
#' @export
#' @import dplyr
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @importFrom data.table fread
#' @importFrom here here
#' @importFrom stringr word
#' @importFrom minfi getAnnotation
#'
CleanDMRResults <- function(annotated_EWAS_result,
                            path_to_combp,
                            expo_levels,
                            expo_labels,
                            annotation_object = IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19,
                            sepages = FALSE) {

  # Create an annotation table
  annot_table <- minfi::getAnnotation(annotation_object)


  # Create a data.frame to contain the final result
  DMR_results <- data.frame()

  # Decompress the DMR analysis results files (ending with *.regions-t.bed)
  for (DMR_file in list.files(path = here::here(path_to_combp), pattern = "\\.regions-t.bed$")) {

    print(DMR_file)
    # Subset the exposure name

    if (isFALSE(sepages)) {
      expo <- stringr::word(DMR_file, 1, sep = "_log2_") %>%
        paste0("_log2")

    } else {
      expo <- stringr::word(DMR_file, 1, 3, sep = "_")

    }

    # Create a subset of annotated data specific for the exposure being analysed
    annot_subset <- dplyr::filter(annotated_EWAS_result, stringr::str_detect(Exposure, expo))

    # List the DMRs with corresponding p values and direction of the estimate
    regions_p <- data.table::fread(here::here(path_to_combp, DMR_file)) %>%

      # Add gene name
      dplyr::mutate(Gene = annot_table$UCSC_RefGene_Name[match(start, annot_table$pos)],

                    # Add exposure name
                    Exposure = unique(annot_subset$Exposure),

                    chrom = as.character(.$'#chrom'),

                    # Add chromosomal coordinates of the DMR
                    Hg19 = stringr::str_c("chr", chrom, ":", start, "-", end),

                    'No. of probes' = n_probes,

                    'Slk p-value' = formatC(.$z_p, format = "e", digits = 2),

                    'Sidak p-value' = formatC(.$z_sidak_p, format = "e", digits = 2),

                    # Add estimate value corresponding to the Exposure and position in the gene
                    dir_effect = annot_subset$Estimate[match(start, annotated_EWAS_result$Position)],

                    # Add categorical variable defining the direction of the estimate
                    'Direction of effect' = factor(dplyr::case_when(
                      dir_effect > 0 ~ "+",
                      TRUE ~ "-"))) %>%

      dplyr::select(Gene, Hg19, Exposure, 'No. of probes':'Sidak p-value', 'Direction of effect') %>%
      droplevels()

    # Append results for a given exposure
    DMR_results <- dplyr::bind_rows(DMR_results, regions_p) %>%
      dplyr::arrange(Gene, Hg19, Exposure)
  }

   DMR_res <- DMR_results %>%
    dplyr::mutate(Exposure = factor(Exposure, levels = expo_levels, labels = expo_labels))

  return(DMR_res)

}
