# Epigenetic footprints: Investigating placental DNA methylation in the context of prenatal exposure to phenols and phthalates

## Article info

Repository contains reproducible analyses for the paper: [Jedynak et al., 2024](https://www.sciencedirect.com/science/article/pii/S0160412024003490).

Corresponding author: [Claire Philippat](claire.philippat@inserm.fr).

## Analysis overview

Analyses depend on R (>= 4.3.0).

Code used to produce the files is found in the corresponding file: `jedynak_epigenetic_2024.R` (executable file).

This analysis was performed under Windows 11 x64 (build 22621) using:    
* [R 4.4.0](https://cran.r-project.org/bin/windows/base) (2024-04-24)    

### Packages

All packages used in the analyses are listed in the `sessionInfo.txt` and below. All "in-house" packages and functions are available in the R/ folder.


### To run

Re-running the analysis requires an additional `data/raw_data` folder that is not shared here. These data can only be provided upon request and after approval by the SEPAGES consortium (contact: [Sarah Lyon-Caen](sarah.lyon-caen@univ-grenoble-alpes.fr)). Running the script: `jedynak_epigenetic_2024.R` will allow to fully reproduce the analyses, figures and tables.


## Repo organization

### data/raw_data/ folder

Analysis input data-files are not made available as they contain sensitive information. The analysis input data files would be:

* `meth_clean.RDS` = pre-processed DNA methylation data (same as in the previous work Jedynak et al. 2023 on SEPAGES), n = 395
* `data_pe_sepages_pau_221010.sas7bdat` = Data provided by Anne Boudier with covariates, standardized and non-standardized exposures, n = 479
* `bSEPAGES_SampleSheet_nodup_20211208.csv` = dataset containing DNA methylation measurement technical factors (batch, plate, chip) info, provided by Lucile Broséus, n = 395
* `SEPAGES_CC.planet_rpc.csv` = data on the cell mix provided by Lucile Broséus, n = 395

Additional file `EWAS_EDEN.RDS` containing data on EWAS results for previous EDEN study on phenols and phthalates (Jedynak et al. 2021, 2022) is not available in the repo due to its large size.


### R/ folder

This folder contains all the in-house functions and packages used for the analyses.


## Session info

```
R version 4.4.0 (2024-04-24)
Platform: x86_64-pc-linux-gnu
Running under: Linux Mint 21.3

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8          LC_NUMERIC=C                  LC_TIME=en_US.UTF-8           LC_COLLATE=en_US.UTF-8       
 [5] LC_MONETARY=es_ES.UTF-8       LC_MESSAGES=en_US.UTF-8       LC_PAPER=es_ES.UTF-8          LC_NAME=es_ES.UTF-8          
 [9] LC_ADDRESS=es_ES.UTF-8        LC_TELEPHONE=es_ES.UTF-8      LC_MEASUREMENT=es_ES.UTF-8    LC_IDENTIFICATION=es_ES.UTF-8

time zone: America/Lima
tzcode source: system (glibc)

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] VennDiagram_1.7.3         futile.logger_1.4.3       lubridate_1.9.3           dplyr_1.0.2              
 [5] purrr_1.0.2               readr_2.1.5               tidyr_1.3.1               tibble_3.2.1             
 [9] ggplot2_3.5.1             tidyverse_2.0.0           stringr_1.5.1             rio_1.1.0                
[13] RobustRegrMod_1.0         readxl_1.4.3              qgcomp_2.15.2             mycor_0.1.1              
[17] Hmisc_5.1-3               Helpers_0.1.0             forcats_1.0.0             DMR_0.2                  
[21] DescriptiveStatistics_2.1 corrplot_0.92             pacman_0.5.1              here_1.0.1               

loaded via a namespace (and not attached):
  [1] matrixStats_1.3.0                                  bitops_1.0-7                                      
  [3] rms_6.8-1                                          doParallel_1.0.17                                 
  [5] naniar_1.1.0                                       httr_1.4.7                                        
  [7] RColorBrewer_1.1-3                                 tools_4.4.0                                       
  [9] doRNG_1.8.6                                        backports_1.5.0                                   
 [11] utf8_1.2.4                                         R6_2.5.1                                          
 [13] HDF5Array_1.32.0                                   jomo_2.7-6                                        
 [15] rhdf5filters_1.16.0                                withr_3.0.0                                       
 [17] gridExtra_2.3                                      base64_2.0.1                                      
 [19] preprocessCore_1.66.0                              quantreg_5.98                                     
 [21] cli_3.6.2                                          Biobase_2.64.0                                    
 [23] textshaping_0.4.0                                  formatR_1.14                                      
 [25] fastDummies_1.7.3                                  sandwich_3.1-0                                    
 [27] officer_0.6.6                                      mvtnorm_1.2-5                                     
 [29] polspline_1.1.25                                   genefilter_1.86.0                                 
 [31] askpass_1.2.0                                      Rsamtools_2.20.0                                  
 [33] systemfonts_1.1.0                                  foreign_0.8-86                                    
 [35] R.utils_2.12.3                                     gfonts_0.2.0                                      
 [37] siggenes_1.78.0                                    illuminaio_0.46.0                                 
 [39] svglite_2.1.3                                      pscl_1.5.9                                        
 [41] scrime_1.3.5                                       limma_3.60.2                                      
 [43] rstudioapi_0.16.0                                  RSQLite_2.3.7                                     
 [45] httpcode_0.3.0                                     generics_0.1.3                                    
 [47] shape_1.4.6.1                                      BiocIO_1.14.0                                     
 [49] xlsx_0.6.5                                         zip_2.3.1                                         
 [51] Matrix_1.6-5                                       fansi_1.0.6                                       
 [53] S4Vectors_0.42.0                                   abind_1.4-5                                       
 [55] R.methodsS3_1.8.2                                  lifecycle_1.0.4                                   
 [57] multcomp_1.4-25                                    yaml_2.3.8                                        
 [59] SummarizedExperiment_1.34.0                        rhdf5_2.48.0                                      
 [61] SparseArray_1.4.8                                  blob_1.2.4                                        
 [63] promises_1.3.0                                     crayon_1.5.2                                      
 [65] mitml_0.4-5                                        lattice_0.22-5                                    
 [67] cowplot_1.1.3                                      GenomicFeatures_1.56.0                            
 [69] annotate_1.82.0                                    xlsxjars_0.6.1                                    
 [71] KEGGREST_1.44.0                                    pillar_1.9.0                                      
 [73] knitr_1.47                                         beanplot_1.3.1                                    
 [75] GenomicRanges_1.56.0                               rjson_0.2.21                                      
 [77] boot_1.3-30                                        bigstatsr_1.5.12                                  
 [79] codetools_0.2-19                                   compareGroups_4.8.0                               
 [81] bigassertr_0.1.6                                   pan_1.9                                           
 [83] glue_1.7.0                                         fontLiberation_0.1.0                              
 [85] data.table_1.15.4                                  vctrs_0.6.5                                       
 [87] png_0.1-8                                          cellranger_1.1.0                                  
 [89] gtable_0.3.5                                       IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1
 [91] cachem_1.1.0                                       xfun_0.44                                         
 [93] S4Arrays_1.4.1                                     mime_0.12                                         
 [95] survey_4.4-2                                       coda_0.19-4.1                                     
 [97] survival_3.5-8                                     rJava_1.0-11                                      
 [99] iterators_1.0.14                                   statmod_1.5.0                                     
[101] TH.data_1.1-2                                      nlme_3.1-163                                      
[103] flock_0.7                                          bit64_4.0.5                                       
[105] fontquiver_0.2.1                                   GenomeInfoDb_1.40.1                               
[107] rprojroot_2.0.4                                    nor1mix_1.3-3                                     
[109] rpart_4.1.23                                       colorspace_2.1-0                                  
[111] BiocGenerics_0.50.0                                DBI_1.2.2                                         
[113] nnet_7.3-19                                        tidyselect_1.2.1                                  
[115] bigparallelr_0.3.2                                 bit_4.0.5                                         
[117] compiler_4.4.0                                     curl_5.2.1                                        
[119] chron_2.3-61                                       glmnet_4.1-8                                      
[121] htmlTable_2.4.2                                    HardyWeinberg_1.7.8                               
[123] SparseM_1.83                                       flextable_0.9.6                                   
[125] mice_3.16.0                                        xml2_1.3.6                                        
[127] fontBitstreamVera_0.1.1                            DelayedArray_0.30.1                               
[129] rtracklayer_1.64.0                                 sfsmisc_1.1-18                                    
[131] checkmate_2.3.1                                    scales_1.3.0                                      
[133] quadprog_1.5-8                                     digest_0.6.35                                     
[135] minqa_1.2.7                                        rmarkdown_2.27                                    
[137] GEOquery_2.72.0                                    XVector_0.44.0                                    
[139] htmltools_0.5.8.1                                  pkgconfig_2.0.3                                   
[141] base64enc_0.1-3                                    lme4_1.1-35.3                                     
[143] sparseMatrixStats_1.16.0                           MatrixGenerics_1.16.0                             
[145] fastmap_1.2.0                                      rlang_1.1.3                                       
[147] htmlwidgets_1.6.4                                  UCSC.utils_1.0.0                                  
[149] shiny_1.8.1.1                                      DelayedMatrixStats_1.26.0                         
[151] zoo_1.8-12                                         jsonlite_1.8.8                                    
[153] BiocParallel_1.38.0                                mclust_6.1.1                                      
[155] R.oo_1.26.0                                        RCurl_1.98-1.14                                   
[157] magrittr_2.0.3                                     kableExtra_1.4.0                                  
[159] Formula_1.2-5                                      GenomeInfoDbData_1.2.12                           
[161] Rhdf5lib_1.26.0                                    munsell_0.5.1                                     
[163] Rcpp_1.0.12                                        gdtools_0.3.7                                     
[165] visdat_0.6.0                                       stringi_1.8.4                                     
[167] zlibbioc_1.50.0                                    MASS_7.3-60                                       
[169] plyr_1.8.9                                         bumphunter_1.46.0                                 
[171] minfi_1.50.0                                       parallel_4.4.0                                    
[173] Biostrings_2.72.0                                  splines_4.4.0                                     
[175] multtest_2.60.0                                    hms_1.1.3                                         
[177] locfit_1.5-9.9                                     uuid_1.2-0                                        
[179] rngtools_1.5.2                                     stats4_4.4.0                                      
[181] futile.options_1.0.1                               crul_1.4.2                                        
[183] XML_3.99-0.16.1                                    evaluate_0.23                                     
[185] mitools_2.4                                        lambda.r_1.2.4                                    
[187] nloptr_2.0.3                                       tzdb_0.4.0                                        
[189] foreach_1.5.2                                      httpuv_1.6.15                                     
[191] MatrixModels_0.5-3                                 openssl_2.2.0                                     
[193] reshape_0.8.9                                      broom_1.0.6                                       
[195] xtable_1.8-4                                       restfulr_0.0.15                                   
[197] Rsolnp_1.16                                        later_1.3.2                                       
[199] viridisLite_0.4.2                                  ragg_1.3.2                                        
[201] truncnorm_1.0-9                                    arm_1.14-4                                        
[203] memoise_2.0.1                                      AnnotationDbi_1.66.0                              
[205] GenomicAlignments_1.40.0                           IRanges_2.38.0                                    
[207] writexl_1.5.0                                      cluster_2.1.6                                     
[209] timechange_0.3.0                                  

```