## Code to recreate the analyses for the paper: Epigenetic footprints: 
# Investigating placental DNA methylation in the context of prenatal exposure to 
# phenols and phthalates"
#
# Authors: Paulina Jedynak, Valérie Siroux, Lucile Broséus, Jörg Tost, Florence Busato,
# Stephan Gabet, Cathrine Thomsen, Amrit K. Sakhi, Azemira Sabaredzovic, Sarah Lyon-Caen, 
# Sam Bayat, Rémy Slama, Claire Philippat, Johanna Lepeule"
#
# DOI: https://doi.org/10.1016/j.envint.2024.108763

## -------------------------------------------------------------------------------------------------

## General setup

# Create a path to save output files
path_main <- "results/main_analyses" # main analyses of the paper
path_suppl <- "results/supplementary_analyses" # supplementary analyses of the paper (not directly
# used in the paper but useful)


## -------------------------------------------------------------------------------------------------

## Install and load packages and functions
install.packages("pacman")
library(pacman)

p_load(here)
install.packages(here("R", "packages_tarballs", "DescriptiveStatistics_2.1.tar.gz"), repos = NULL, type = "source")
install.packages(here("R", "packages_tarballs", "DMR_0.2.tar.gz"), repos = NULL, type = "source")
install.packages(here("R", "packages_tarballs", "Helpers_0.1.0.tar.gz"), repos = NULL, type = "source")
install.packages(here("R", "packages_tarballs", "RobustRegrMod_1.0.tar.gz"), repos = NULL, type = "source")

p_load(corrplot,
       DescriptiveStatistics,
       DMR,
       forcats,
       grDevices,
       Helpers,
       Hmisc,
       mycor,
       qgcomp,
       readxl,
       RobustRegrMod,
       rio,
       stats,
       stringr,
       tidyverse, # NOTE: Keep dplyr package version 1.0. as newer versions will not work (changes in pivot_wider function) with ExposureCharacteristicsSepages.R
       VennDiagram)

install.packages(here("R", "packages_tarballs", "dplyr_1.0.2.tar.gz"), repos = NULL, type = "source")
library(dplyr)

source(here("R", "ExposureCharSepages.R"))
source(here("R", "PrepareDatasetsSepagesAllComp.R"))
source(here("R", "CpGsInDMRsSEPAGES.R"))
source(here("R", "GenesInDMRsSEPAGES.R"))
source(here("R", "CleanDMRResults.R"))
source(here("R", "ReadExcel.R"))


## -------------------------------------------------------------------------------------------------

## Load datasets

# Note: re-running the analyses requires an additional data/raw_data folder not shared here, containing the different 
# data sets presented hereafter.
# These data can only be provided upon request.

# Load the methylation dataset with removed outliers and excessive missing values.
# This dataset was created for the previous SEPAGES study on triclosan (Jedynak et al. 2023)
meth_clean <- import(here("data", "raw_data", "meth_clean.RDS")) # 395 x 752,577

# Data on exposure and covariates
# Data provided by Anne Boudier with covariates, standardized and non-standardized exposures
cov_expo_data <- read_sas(here("data", "raw_data", "data_pe_sepages_pau_221010.sas7bdat"), here("data", "raw_data", "formats.sas7bcat"))

# Load the dataset containing DNA methylation measurement technical factors (batch, plate, chip) info, n = 395, provided by Lucile Broséus 02-11-2021
tech_conf <- import(here("data", "raw_data", "SEPAGES_SampleSheet_nodup_20211208.csv"))

# Load the cell mix provided by Lucile Broséus 07-02-2022
SEPAGES_CC_planet_rpc <- import(here("data", "raw_data", "SEPAGES_CC.planet_rpc.csv")) %>%
  mutate(id = as.character(id))

# Load EWAS results for previous EDEN study on phenols and phthalates (Jedynak et al. 2021, 2022)
EWAS_EDEN <- import(here("data", "EWAS_EDEN.RDS"))

# Load LOD values for all exposures
lod_values <- import(here("data", "LODs.csv"))

# Load the list of genes identified as differentially methylated in association with triclosan in the DMR analyses (Appendix Table 4 from Jedynak et al. 2021)
genes_DMR_2probes <- import(here("data", "genes_DMR_2probes.csv"))

# Load the list of imprinted genes (same as in Jedynak et al. 2021)
imprinted_genes <- import(here("data", "201028_imprinted_genes_list_cp200514_corrected.csv")) %>%
  arrange(gene_symbol) # in my study, only 297 of these genes are represented (in the total of 20,283 Illumina genes available for this study)

# Load list of genes identified in other studies as associated with phthalates
phthalates_genes <- import(here("data", "genes_phthalates.csv"))

# Enriched terms based on at least 20 genes provided by Lucile Broséus
dbgap_filtered_FDR001_20GenesPerTerm <- import(here("data", "nl_EWAS_agnostic_quant_res_0.001_enrichR_dbgap_filtered_FDR001_20GenesPerTerm.rds")) 


## -------------------------------------------------------------------------------------------------

## Create variables names

# Manually created lists of compounds, models variables, etc.

# Original names for non-standardized exposure concentrations

# Name patterns:
# mo_CCC_total_i_t2 or  mo_CCC_total_i_t3 for continuous phenols
# mo_CCC_total_string_t2 or mo_CCC_total_string_t3 for BUPA and BPS [contains measured concentrations, if 0 (non detected) replaced with LOD / sqrt(2)]
# mo_CCC_i_t2 or mo_CCC_i_t3 for phthalates
# mo_CCC_ms_i_t2 or mo_CCC_ms_i_t3 for phthalate molar sum
# where CCC stands for the name of the compound

exposures_ns <- c(
  # Phenols
  "mo_BPA_total_i_t2", "mo_BPA_total_i_t3",
  "mo_BPS_total_string_t2", "mo_BPS_total_string_t3",
  "mo_BUPA_total_string_t2", "mo_BUPA_total_string_t3",
  "mo_ETPA_total_i_t2", "mo_ETPA_total_i_t3",
  "mo_MEPA_total_i_t2", "mo_MEPA_total_i_t3",
  "mo_OXBE_total_i_t2", "mo_OXBE_total_i_t3",
  "mo_PRPA_total_i_t2", "mo_PRPA_total_i_t3",
  
  # Phthalates
  "mo_MBzP_i_t2", "mo_MBzP_i_t3",
  "mo_MEP_i_t2", "mo_MEP_i_t3",
  "mo_MiBP_i_t2", "mo_MiBP_i_t3",
  "mo_MnBP_i_t2", "mo_MnBP_i_t3",
  
  "mo_ohMPHP_i_t2", "mo_ohMPHP_i_t3",
  
  "mo_MECPP_i_t2", "mo_MECPP_i_t3",
  "mo_MEHHP_i_t2", "mo_MEHHP_i_t3",
  "mo_MEHP_i_t2", "mo_MEHP_i_t3",
  "mo_MEOHP_i_t2", "mo_MEOHP_i_t3",
  "mo_MMCHP_i_t2", "mo_MMCHP_i_t3",
  "mo_DEHP_ms_i_t2", "mo_DEHP_ms_i_t3",
  
  "mo_cxMiNP_i_t2", "mo_cxMiNP_i_t3",
  "mo_ohMiNP_i_t2", "mo_ohMiNP_i_t3",
  "mo_oxoMiNP_i_t2", "mo_oxoMiNP_i_t3",
  "mo_DiNP_ms_i_t2", "mo_DiNP_ms_i_t3",
  
  "mo_ohMINCH_i_t2", "mo_ohMINCH_i_t3",
  "mo_oxoMINCH_i_t2", "mo_oxoMINCH_i_t3",
  "mo_DINCH_ms_i_t2", "mo_DINCH_ms_i_t3"
)

# Original names for exposure concentrations standardized acc. to Mortamais et al. 2012

# Name patterns:
# mo_CCC_total_i_cor_t2 or mo_CCC_total_i_cor_t3 for continuous phenols
# mo_CCC_total_cat_t2 or mo_CCC_total_cat_t3 for BUPA and BPS
# mo_CCC_i_cor_t2 or mo_CCC_i_cor_t3 for phthalates
# mo_CCC_ms_i_cor_t2 or mo_CCC_ms_i_cor_t2 for phthalate molar sums
# where CCC stands for the name of the compound

exposures_sd <- c(
  
  # Phenols
  "mo_BPA_total_i_cor_t2", "mo_BPA_total_i_cor_t3",
  "mo_BPS_total_cat_t2", "mo_BPS_total_cat_t3",
  "mo_BUPA_total_cat_t2", "mo_BUPA_total_cat_t3",
  "mo_ETPA_total_i_cor_t2", "mo_ETPA_total_i_cor_t3",
  "mo_MEPA_total_i_cor_t2", "mo_MEPA_total_i_cor_t3",
  "mo_OXBE_total_i_cor_t2", "mo_OXBE_total_i_cor_t3",
  "mo_PRPA_total_i_cor_t2", "mo_PRPA_total_i_cor_t3",
  
  # Phthalates
  "mo_MBzP_i_cor_t2", "mo_MBzP_i_cor_t3",
  "mo_MEP_i_cor_t2", "mo_MEP_i_cor_t3",  
  "mo_MiBP_i_cor_t2", "mo_MiBP_i_cor_t3",  
  "mo_MnBP_i_cor_t2", "mo_MnBP_i_cor_t3",
  
  "mo_ohMPHP_i_cor_t2", "mo_ohMPHP_i_cor_t3",
  
  "mo_MECPP_i_cor_t2", "mo_MECPP_i_cor_t3",
  "mo_MEHHP_i_cor_t2", "mo_MEHHP_i_cor_t3",
  "mo_MEHP_i_cor_t2", "mo_MEHP_i_cor_t3",
  "mo_MEOHP_i_cor_t2", "mo_MEOHP_i_cor_t3",
  "mo_MMCHP_i_cor_t2", "mo_MMCHP_i_cor_t3",
  "mo_DEHP_ms_i_cor_t2", "mo_DEHP_ms_i_cor_t3",
  
  "mo_cxMiNP_i_cor_t2", "mo_cxMiNP_i_cor_t3",
  "mo_ohMiNP_i_cor_t2", "mo_ohMiNP_i_cor_t3",
  "mo_oxoMiNP_i_cor_t2", "mo_oxoMiNP_i_cor_t3",
  "mo_DiNP_ms_i_cor_t2", "mo_DiNP_ms_i_cor_t3",
  
  "mo_ohMINCH_i_cor_t2", "mo_ohMINCH_i_cor_t3",
  "mo_oxoMINCH_i_cor_t2", "mo_oxoMINCH_i_cor_t3",
  "mo_DINCH_ms_i_cor_t2", "mo_DINCH_ms_i_cor_t3")


# Exposure names to be used in regressions in candidate approach, restricted to those overlapping with EDEN

exposures_overlap <- c(
  
  # Phenols
  "aver_log2_BPA", 
  "aver_log2_BP3", 
  "aver_2cat_BUPA",
  "aver_log2_ETPA", 
  "aver_log2_MEPA", 
  "aver_log2_PRPA",
  
  # Phthalates
  "aver_log2_MEP", 
  "aver_log2_MiBP", 
  "aver_log2_MnBP", 
  "aver_log2_MBzP", 
  "aver_log2_MECPP", 
  "aver_log2_MEHHP", 
  "aver_log2_MEHP", 
  "aver_log2_MEOHP",
  "aver_log2_DEHP_red")

# Exposure full names to be used in regressions in candidate approach, restricted to those overlapping with EDEN
expo_labels_overlap <- c(
  
  # Phenols
  "Bisphenol A (BPA)",
  "Benzophenone-3 (BP-3)", 
  "Butylparaben (BUPA)",
  "Ethylparaben (ETPA)", 
  "Methylparaben (MEPA)", 
  "Propylparaben (PRPA)",
  
  # Phthalates
  "Monoethyl phthalate (MEP)",
  "Mono-iso-butyl phthalate (MiBP)",                     
  "Mono-n-butyl phthalate (MnBP)",
  "Monobenzyl phthalate (MBzP)",
  "Mono-2-ethyl-5-carboxypentyl phthalate (MECPP)",
  "Mono-2-ethyl-5-hydroxyhexyl phthalate (MEHHP)",
  "Mono-2-ethylhexyl phthalate (MEHP)",
  "Mono-2-ethyl-5-oxohexyl phthalate (MEOHP)",
  "Molar sum of di(2-ethylhexyl) phthalate (DEHP) metabolites (ΣDEHP)"
)


# Exposure names to be used in regressions in exploratory study
exposures <- c(
  
  # Phenols
  "aver_log2_BPA",
  "aver_log2_BP3", 
  "aver_log2_ETPA", 
  "aver_log2_MEPA", 
  "aver_log2_PRPA",
  
  # Phthalates
  "aver_log2_MEP", 
  "aver_log2_MiBP", 
  "aver_log2_MnBP", 
  "aver_log2_MBzP", 
  "aver_log2_ohMPHP", 
  "aver_log2_MECPP", 
  "aver_log2_MEHHP", 
  "aver_log2_MEHP", 
  "aver_log2_MEOHP", 
  "aver_log2_MMCHP", 
  "aver_log2_DEHP", 
  "aver_log2_cxMiNP", 
  "aver_log2_ohMiNP", 
  "aver_log2_oxoMiNP", 
  "aver_log2_DiNP", 
  "aver_log2_ohMINCH" , 
  "aver_log2_oxoMINCH", 
  "aver_log2_DINCH")


# Categorical exposure names to be used in regressions in exploratory study
exposures_cat <- c("aver_2cat_BPS",
                   "aver_2cat_BUPA")


# Exposure full names to be used in regressions in exploratory study
expo_labels <- c(
  
  # Phenols
  "Bisphenol A (BPA)", 
  "Benzophenone-3 (BP-3)", 
  "Ethylparaben (ETPA)", 
  "Methylparaben (MEPA)", 
  "Propylparaben (PRPA)",
  
  # Phthalates
  "Monoethyl phthalate (MEP)",
  "Mono-iso-butyl phthalate (MiBP)",                     
  "Mono-n-butyl phthalate (MnBP)",
  "Monobenzyl phthalate (MBzP)",
  "Mono-6-hydroxy-propylheptyl phthalate (OH-MPHP)",
  "Mono-2-ethyl-5-carboxypentyl phthalate (MECPP)",
  "Mono-2-ethyl-5-hydroxyhexyl phthalate (MEHHP)",
  "Mono-2-ethylhexyl phthalate (MEHP)",
  "Mono-2-ethyl-5oxo-hexyl phthalate (MEOHP)",
  "Mono-2-carboxymethylhexyl phthalate (MMCHP)", 
  "Molar sum of di(2-ethylhexyl) phthalate (DEHP) metabolites (ΣDEHP)",
  "Mono-4-methyl-7-carboxyoctyl phthalate (cx-MiNP)",
  "Mono-4-methyl-7-hydroxyoctyl phthalate (OH-MiNP)",
  "Mono-4-methyl-7-oxooctyl phthalate (oxo-MiNP)",
  "Molar sum of di-isononyl phthalate (DiNP) metabolites (ΣDiNP)",
  "2-(((hydroxy-4-methyloctyl)oxy)carbonyl) cyclohexanecarboxylic acid (OH-MINCH)",
  "2-(((4-methyl-7-oxooctyl)oxy)carbonyl) cyclohexanecarboxylic acid (oxo-MINCH)",
  "Molar sum of 1,2-cyclohexane dicarboxylic acid di-isononyl ester (DINCH) metabolites (ΣDiNCH)"
)

# Categorical exposure full names to be used in regressions in exploratory study
expo_labels_cat <- c("Bisphenol S (BPS)",
                     "Butylparaben (BUPA)")


# List of confounders used as adjustment factors in the main statistical analyses
confounders <- c("season_conc", 
                 "mother_act_smoke",
                 "parity",
                 "mother_edu",
                 "mother_bmi",
                 "child_sex",
                 "mother_age",
                 "gest_age")

# List of technical confounders for the EWAS and DMR study
tech_conf_EWAS_DMR <- c("batch", "plate", "chip")

# List cell types obtained by planet package, to be added to the models in the sensitivity analyses. Since cell types sum up to 1, one needs to be removed to avoid singularity. Remove nRBC type as least represented 
cell_mix_names <- c("Endothelial", "Hofbauer", "Stromal", "Syncytiotrophoblast", "Trophoblasts")

# Create a list of CpGs showing linear associations with exposure in EDEN 
linear_CpGs_EDEN <- c("cg02366345", "cg04214090", "cg02630105", "cg03341641")


## -------------------------------------------------------------------------------------------------

## Data preparation

## Methylation data

# Calculate mean methylation and its standard deviation for each CpG in the sample
mean_meth <- data.frame(Mean_meth = colMeans(meth_clean, na.rm = TRUE),
                        SD_meth = colSds(meth_clean, na.rm = TRUE)) %>%
  rownames_to_column("CpG")


## Exposures, covariates, technical factors

# Add DEHP molar sum containing only the compounds used in EDEN: MEHP, MEHHP, MEOHP, and MECPP
cov_expo_data <- cov_expo_data %>%
  mutate(
    mo_DEHP_ms_i_cor_red_t2 = mo_MEOHP_i_cor_t2/292 + mo_MECPP_i_cor_t2/308 + mo_MEHP_i_cor_t2/278 + mo_MEHHP_i_cor_t2/294,
    
    mo_DEHP_ms_i_cor_red_t3 = mo_MEOHP_i_cor_t3/292 + mo_MECPP_i_cor_t3/308 + mo_MEHP_i_cor_t3/278 + mo_MEHHP_i_cor_t3/294)

# Prepare dataset containing all individuals (merge confounders, technical factors, exposures, and cell mix)
expo_cov_tech_conf_all <- PrepareDatasetsSepagesAllComp(
  cov_expo_data = cov_expo_data,
  technical_covariates_data = tech_conf,
  exposures = c(c(exposures_sd, "mo_DEHP_ms_i_cor_red_t2", "mo_DEHP_ms_i_cor_red_t3"), exposures_ns),
  cell_mix = SEPAGES_CC_planet_rpc) %>%
  
  # Remove women (n = 5) that dropped out before delivery
  drop_na("child_sex")

# Save dataset
export(expo_cov_tech_conf_all, here(path_main, "expo_cov_tech_conf_all.RDS"))

# Restrict covariates data to methylation IDs
expo_cov_tech_conf_miss <- expo_cov_tech_conf_all[expo_cov_tech_conf_all$id %in% rownames(meth_clean), ]

# Save dataset with missing values
export(expo_cov_tech_conf_miss, here(path_main, "expo_cov_tech_conf_miss.RDS"))


## -------------------------------------------------------------------------------------------------

## Impute missing covariates

expo_cov_tech_conf <- ImputeMissingCov(data_w_missing = expo_cov_tech_conf_miss) %>%
  
  # Select variables of interest
  select(id, aver_log2_BP3:Trophoblasts)

# Drop rownames
rownames(expo_cov_tech_conf) <- NULL

# Save datset
# Confounders, technical factors, exposures, and cell mix - imputed values for confounders
export(expo_cov_tech_conf, here(path_main, "expo_cov_tech_conf.RDS"))


## -------------------------------------------------------------------------------------------------

## Population characteristics and descriptive statistics

# Population characteristics for the mother-child pairs included in the study (n = 395)

DescriptiveStatistics::PopulationCharacteristics(
  input_data = expo_cov_tech_conf_miss,
  variable_list = confounders,
  path = path_main,
  file_name = "population_charactersitics_incl")


## -------------------------------------------------------------------------------------------------

# Comparison of the population characteristics between included (395) and excluded (84) individuals

# Restrict IDs to those present in the methylation dataset
cov_expo_data_excl <- expo_cov_tech_conf_all[expo_cov_tech_conf_all$id %nin% rownames(meth_clean), ] %>%
  
  # Remove women (n = 5) that dropped out before delivery
  drop_na("child_sex")

# Calculate descriptive statistics for the excluded individuals
DescriptiveStatistics::PopulationCharacteristics(
  input_data = cov_expo_data_excl,
  variable_list = confounders,
  path = path_main,
  file_name = "population_charactersitics_excl")


## -------------------------------------------------------------------------------------------------
# Test if there is a sign. difference between included and excluded individuals
cov_comp <- expo_cov_tech_conf_all %>%
  select(id, all_of(confounders)) %>%
  mutate(Inclusion = factor(case_when(id %in% cov_expo_data_excl$id ~ "Excl",
                                      TRUE ~ "Incl")))

# Comparison of categorical variables
chi_res <- data.frame()
for (cname in colnames(select(cov_comp, season_conc:child_sex))) {
  Chi <- chisq.test(table(cov_comp$Inclusion, cov_comp[, cname]))
  Chi <- data.frame(Comp = cname, p_val = Chi$p.value)
  chi_res <- rbind(chi_res, Chi)
}
chi_res # There are no differences between included and excluded women (chisq test p-value > 0.05)


## -------------------------------------------------------------------------------------------------

## Exposure descriptive statistics

# Maternal exposure descriptive statistics for phenols and phthalates

expo_cov_tech_conf_miss_sel <- expo_cov_tech_conf_miss %>% 
  select(!matches("log2_|_cat_|aver_|_i_|_red_"), 
         -c(mother_age:Trophoblasts)) %>% 
  rename_at(.vars = vars(contains("_string")),
            .funs = funs(sub("_string", "", .)))

exposure_char <- ExposureCharSepages(input_data = expo_cov_tech_conf_miss_sel,
                                     lod_values = lod_values,
                                     path = path_main)


## -------------------------------------------------------------------------------------------------

# Comparison of the averaged exposure levels between included (395) and excluded (84) individuals

# Test if there is a sign. difference between included and excluded individuals for exposure concentrations
expo_comp <- expo_cov_tech_conf_all %>%
  select(id, contains("aver_log2"), contains("aver_2cat")) %>%
  mutate(Inclusion = factor(case_when(id %in% cov_expo_data_excl$id ~ "Excl",
                                      TRUE ~ "Incl")))

# Comparison of categorical exposures
chi_res_expo <- data.frame()
for (cname in colnames(select(expo_comp, contains("aver_2cat")))) {
  Chi <- chisq.test(table(expo_comp$Inclusion, expo_comp[, cname]))
  Chi <- data.frame(Comp = cname, Chi = Chi$p.value)
  chi_res_expo <- rbind(chi_res_expo, Chi)
}
chi_res_expo # There are no differences between included and excluded women (chisq test p-value > 0.05)


## -------------------------------------------------------------------------------------------------

# Comparison of numerical variables
expo_comp %>%
  select(id, Inclusion, !contains("aver_2cat")) %>%
  pivot_longer(cols = -c(id, Inclusion),
               names_to = "Conf",
               values_to = "Value") %>%
  group_by(Conf) %>%
  do(w = wilcox.test(Value ~ Inclusion,
                     data = .,
                     paired = FALSE)) %>%
  summarise(Conf, Wilcox = w$p.value) # There are no differences between included and excluded women (Wilcoxon test p-value > 0.05)


## -------------------------------------------------------------------------------------------------
# Comparison of numerical variables
wilcox <- cov_comp %>%
  select(id, Inclusion, contains("age")) %>%
  pivot_longer(cols = -c(id, Inclusion),
               names_to = "Comp",
               values_to = "Value") %>%
  group_by(Comp) %>%
  do(w = wilcox.test(Value ~ Inclusion,
                     data = .,
                     paired = FALSE)) %>%
  summarise(Comp, p_val = w$p.value) # There is a difference in gestational age between included and excluded women (Wilcoxon test p-value = 0.015)


## -------------------------------------------------------------------------------------------------
# Calculate Spearman's correlations between averaged exposure concentrations
corr_data <- expo_cov_tech_conf |> 
  select(starts_with("aver")) |> 
  select(-matches("BUPA|BPS|DEHP_red")) |> 
  rename_all(~str_remove(., "aver_log2_")) |> 
  select(BP3, BPA, ETPA, MEPA, PRPA, 
         MEP, MnBP, MiBP, MBzP, 
         ohMPHP, 
         DiNP, cxMiNP, oxoMiNP, ohMiNP, 
         DEHP, MMCHP, MECPP, MEHHP, MEHP, MEOHP,
         DINCH, ohMINCH, oxoMINCH) # + BUPA and BPS

corr <- corr_data |> 
  mycor(method = "spearman", na.action = na.omit)

col2 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3",
                           "#92C5DE", "#D1E5F0", "#FFFFFF",
                           "#FDDBC7", "#F4A582", "#D6604D",
                           "#B2182B", "#67001F"))

jpeg(width = 3600,
     height = 3600, 
     filename = here("results", "correlation_compounds.jpeg"), 
     res = 300)

corrplot(corr = corr$r,
         p.mat = corr$P,
         tl.srt = 45,
         tl.col = "black",
         type = "lower", 
         insig = "pch", 
         sig.level = 0.1, 
         pch.cex = 0.9,
         col = col2(200), 
         addCoef.col = 'black')

dev.off()

# Calculate the effective number of independent exposures (taking into account the correlation structure between the exposures)

# Calculate the correlation matrix
cormat <- cor(corr_data, use = "pairwise.complete.obs", method = "spearman")

# Calculate the effective number of independent exposures ("meff") based on Li 2012 method (see Jedynak et al. 2020 on HELIX cohort) - this will be useful for the answer to the reviewers
lambdas <- eigen(cormat)$values
ncol(cormat) - sum((lambdas > 1) * (lambdas - 1)) # 11.04922 + BUPA and BPS, 13 in total


## -------------------------------------------------------------------------------------------------

# Remove upper and lower 1% of extreme concentrations from each continuous exposure (influential values)

expo_cov_tech_conf_miss_long <- expo_cov_tech_conf_miss_sel %>%     
  pivot_longer(cols = -id,
               names_to = c("Exposure", "Trim"),
               values_to = "Conc",
               names_pattern = "(.*)_(.*)") 

# Create the dataset for plotting later
exposure_char_red_plot <- expo_cov_tech_conf_miss_long %>%
  group_by(id, Exposure) %>%
  
  # Calculate concentration averaged over T2 and T3
  summarize(aver_expo_conc_Std = mean(Conc, na.rm = TRUE)) %>%
  ungroup() %>% 
  group_by(Exposure) %>% 
  mutate_at(vars(aver_expo_conc_Std),
            function(x) x <- ifelse(x < quantile(x, 0.01) |
                                      x > quantile(x, 0.99),
                                    NA,
                                    x)) %>%
  
  mutate(expo_group = case_when(Exposure %in% c("BPA",
                                                "BPS",
                                                "BP3",
                                                "MEPA",
                                                "ETPA",
                                                "PRPA",
                                                "BUPA") ~ "Phenol",
                                Exposure %in% c("OhMINCH",
                                                "oxoMINCH",
                                                "DINCH") ~ "DINCH",
                                .default = "Phthalate"),
         
         expo_group = factor(expo_group, levels = c("Phenol", "Phthalate", "DINCH")))

exposure_char_red <- exposure_char_red_plot %>% 
  
  # Calculate basic statistics for the dataset with influential values removed
  summarise("p5_av_red" = round(quantile(aver_expo_conc_Std, na.rm = TRUE, probs = 0.05), 2),
            "median_av_red" = round(quantile(aver_expo_conc_Std, na.rm = TRUE, probs = 0.5), 2),
            "p95_av_red" = round(quantile(aver_expo_conc_Std, na.rm = TRUE, probs = 0.95), 1)) %>% 
  right_join(lod_values, by = "Exposure")


## -------------------------------------------------------------------------------------------------
# Create new dataset (Confounders, technical factors, exposures, and cell mix - imputed values for confounders) with removed individuals that have influential values for exposure
quant_expo <- expo_cov_tech_conf %>%
  mutate_at(vars(contains("aver_log2_")),
            function(x) x <- ifelse(x < quantile(x, 0.01) |
                                      x > quantile(x, 0.99),
                                    NA,
                                    x))

# Save dataset -  THIS WILL BE THE DATASET USED IN ALL MAIN ANALYSES
export(quant_expo, here(path_main, "expo_cov_tech_conf_quant.RDS"))


## -------------------------------------------------------------------------------------------------

## I. Agnostic approach (exploratory analysis)

# All CpGs and compounds available in the SEPAGES would be considered in the EWAS (linear + 
# non-linear associations). Two complementary analyses will be run: including all available 
# individuals and population reduced for the extreme upper and lower 1% of exposure values (8 women for each compound).


# Continuous exposures
EWAS_agnostic_cont_sex <- data.frame()

for (sex in c("Male", "Female")) {
  
  expo_cov_tech_conf_sex <- filter(expo_cov_tech_conf, child_sex == sex)
  meth_clean_restr_sex <- meth_clean[rownames(meth_clean) %in% expo_cov_tech_conf_sex$id, ]
  
  EWAS_sex <- RobustRegressions::SpecificLociMethylationRegressions(
    meth_data = meth_clean_restr_sex,
    expo_cov = expo_cov_tech_conf_sex,
    exposures = exposures,
    confounders = confounders[confounders != "child_sex"],
    technical_confounders = c(tech_conf_EWAS_DMR, cell_mix_names),
    maxit = 400,
    annotation_object = IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19,
    expo_labels = expo_labels,
    path = path_main,
    file_name = paste0("1.EWAS_agnostic_sex/EWAS_agnostic_cont_", sex),
  ) %>%
    mutate(child_sex = sex)
  
  EWAS_agnostic_cont_sex <- bind_rows(EWAS_agnostic_cont_sex, EWAS_sex)
}

# Save EWAS results
export(EWAS_agnostic_cont_sex, here(path_main, "1.EWAS_agnostic_sex", "EWAS_agnostic_cont_sex.RDS"))


## -------------------------------------------------------------------------------------------------
# Categorical exposures
EWAS_agnostic_cat_sex <- data.frame()

for (sex in c("Male", "Female")) {
  
  expo_cov_tech_conf_sex <- filter(expo_cov_tech_conf, child_sex == sex)
  meth_clean_restr_sex <- meth_clean[rownames(meth_clean) %in% expo_cov_tech_conf_sex$id, ]
  
  EWAS_sex <- RobustRegressions::SpecificLociMethylationRegressions(
    meth_data = meth_clean_restr_sex,
    expo_cov = expo_cov_tech_conf_sex,
    exposures = exposures_cat,
    confounders = confounders[confounders != "child_sex"],
    technical_confounders = c(tech_conf_EWAS_DMR, cell_mix_names),
    maxit = 400,
    annotation_object = IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19,
    expo_labels = expo_labels_cat,
    path = path_main,
    file_name = paste0("1.EWAS_agnostic_sex/EWAS_agnostic_cat_", sex)
  ) %>%
    mutate(child_sex = sex)
  
  EWAS_agnostic_cat_sex <- bind_rows(EWAS_agnostic_cat_sex, EWAS_sex)
}

export(EWAS_agnostic_cat_sex, here(path_main, "1.EWAS_agnostic_sex", "EWAS_agnostic_cat_sex.RDS"))


## -------------------------------------------------------------------------------------------------
# Merge results for continuous and categorical exposures and save the result
EWAS_agnostic_sex <- bind_rows(EWAS_agnostic_cont_sex, EWAS_agnostic_cat_sex)
rm(EWAS_agnostic_cont_sex, EWAS_agnostic_cat_sex)

export(EWAS_agnostic_sex, here(path_main, "1.EWAS_agnostic_sex", "EWAS_agnostic_sex.RDS"))


## -------------------------------------------------------------------------------------------------
# Print how many linear associations were detected
EWAS_agnostic_sex_res <- EWAS_agnostic_sex %>%
  
  # Select only significant associations
  filter(p_value_FDR < 0.05) %>%
  group_by(child_sex, Exposure) %>%
  summarise(n_CpGs = n()) %>%
  arrange(child_sex, desc(n_CpGs))


## -------------------------------------------------------------------------------------------------

## 2. Non-linear agnostic sex-stratified EWAS on whole population (n = 395)

# Check if any of the continuous exposure-DNA methylation associations are non-linear.

nl_assoc_EWAS_agnostic_sex <- data.frame()

for (sex in c("Male", "Female")) {
  
  expo_cov_tech_conf_sex <- filter(expo_cov_tech_conf, child_sex == sex)
  meth_clean_restr_sex <- meth_clean[rownames(meth_clean) %in% expo_cov_tech_conf_sex$id, ]
  
  nl_assoc_EWAS_sex <-
    RobustRegressions::NonLinearityCheck(meth_data = meth_clean_restr_sex,
                                         expo_cov = expo_cov_tech_conf_sex,
                                         exposures = exposures,
                                         confounders = confounders[confounders != "child_sex"],
                                         technical_confounders = c(tech_conf_EWAS_DMR, cell_mix_names),
                                         maxit = 400,
                                         path = path_main,
                                         file_name = paste0("2.nl_EWAS_agnostic_sex/nl_assoc_EWAS_agnostic_", sex)) %>%
    mutate(child_sex = sex)
  
  nl_assoc_EWAS_agnostic_sex <- bind_rows(nl_assoc_EWAS_agnostic_sex, nl_assoc_EWAS_sex)
}

# Save the result
export(nl_assoc_EWAS_agnostic_sex, here(path_main, "2.nl_EWAS_agnostic_sex", "nl_assoc_EWAS_agnostic_sex.RDS"))


## -------------------------------------------------------------------------------------------------

# After identification of the non-linearly associated CpGs, the nominal p-values in the agnostic sex-stratified
# EWAS for the CpGs that showed non-linear associations should be replaced by the nominal p-values of the non-linear 
# associations. The FDR-corrected p-values will be replaced by the FDR-corrected p-values for the non-linear associations.

# Replace FDR-corrected p-values in the main EWAS for the CpGs that showed non-linear associations (spline total FDR-corrected p-value < 0.05) and for these CpGs replace nominal p-values with spline total p-value
nl_EWAS_agnostic_sex <- full_join(EWAS_agnostic_sex,
                                  select(nl_assoc_EWAS_agnostic_sex, CpG, Exposure, child_sex, spline_tot_FDR, spline_tot_pval),
                                  by = c("CpG", "Exposure", "child_sex")) %>%
  mutate(p_value_FDR = case_when(spline_tot_FDR < 0.05 ~ spline_tot_FDR,
                                 TRUE ~ p_value_FDR),
         raw_p_value = case_when(spline_tot_FDR < 0.05 ~ spline_tot_pval,
                                 TRUE ~ raw_p_value),
         Estimate_CI = case_when(spline_tot_FDR < 0.05 ~ "Non-linear",
                                 TRUE ~ Estimate_CI)) %>%
  select(-c(spline_tot_FDR, spline_tot_pval))

# Save the result
export(nl_EWAS_agnostic_sex, here(path_main, "2.nl_EWAS_agnostic_sex", "nl_EWAS_agnostic_sex.RDS"))


# Results: we observed a major inflation of significant hits for some compounds for the initial population 
# (especially for non-linear associations with MECPP in females).
# Therefore, an additional analysis will be run after removal of 1% upper and 1% lower extreme exposure 
# values (calculated separately for each compound).

# Only the results that overlap between the whole (n = 395) and reduced (n = 387) population will be considered \
# and discussed.


## -------------------------------------------------------------------------------------------------

## 3. Linear agnostic sex-stratified EWAS on reduced population (n = 387)

# Continuous exposures
EWAS_agnostic_quant_cont_sex <- data.frame()

for (sex in c("Male", "Female")) {
  
  expo_cov_tech_conf_sex <- filter(quant_expo, child_sex == sex)
  meth_clean_restr_sex <- meth_clean[rownames(meth_clean) %in% expo_cov_tech_conf_sex$id, ]
  
  EWAS_quant <- RobustRegressions::SpecificLociMethylationRegressions(
    meth_data = meth_clean_restr_sex,
    expo_cov = expo_cov_tech_conf_sex,
    exposures = exposures,
    confounders = confounders[confounders != "child_sex"],
    technical_confounders = c(tech_conf_EWAS_DMR, cell_mix_names),
    maxit = 400,
    annotation_object = IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19,
    expo_labels = expo_labels,
    path = path_main,
    file_name = paste0("EWAS_agnostic_quant_cont_", sex),
  ) %>%
    mutate(child_sex = sex)
  
  EWAS_agnostic_quant_cont_sex <- bind_rows(EWAS_agnostic_quant_cont_sex, EWAS_quant)
}

export(EWAS_agnostic_quant_cont_sex, here(path_main, "3.EWAS_agnostic_quant_sex", "EWAS_agnostic_quant_cont_sex.RDS"))

# For categorical exposures there is no need to run this analysis (as no influential values for dichotomised variables)


## -------------------------------------------------------------------------------------------------
# Merge results for continuous with extreme values removed and original categorical exposures (as for categorized exposures there was no need to remove influential values) and save the result
EWAS_agnostic_quant_sex <- bind_rows(EWAS_agnostic_quant_cont_sex, EWAS_agnostic_cat_sex)
rm(EWAS_agnostic_quant_cont_sex, EWAS_agnostic_cat_sex)

export(EWAS_agnostic_quant_sex, here(path_main, "3.EWAS_agnostic_quant_sex", "EWAS_agnostic_quant_sex.RDS"))


## -------------------------------------------------------------------------------------------------

## 4. Non-linear agnostic sex-stratified EWAS on reduced population (n = 387)

nl_assoc_EWAS_agnostic_quant_sex <- data.frame()

for (sex in c("Male", "Female")) {
  
  expo_cov_tech_conf_sex <- filter(quant_expo, child_sex == sex)
  meth_clean_restr_sex <- meth_clean[rownames(meth_clean) %in% expo_cov_tech_conf_sex$id, ]
  
  nl_EWAS_sex_quant <-
    RobustRegressions::NonLinearityCheck(meth_data = meth_clean_restr_sex,
                                         expo_cov = expo_cov_tech_conf_sex,
                                         exposures = exposures,
                                         confounders = confounders[confounders != "child_sex"],
                                         technical_confounders = c(tech_conf_EWAS_DMR, cell_mix_names),
                                         maxit = 400,
                                         path = path_main,
                                         file_name = paste0("nl_assoc_EWAS_agnostic_quant_", sex)) %>%
    mutate(child_sex = sex)
  
  nl_assoc_EWAS_agnostic_quant_sex <- bind_rows(nl_assoc_EWAS_agnostic_quant_sex, nl_EWAS_sex_quant)
}

export(nl_assoc_EWAS_agnostic_quant_sex, here(path_main, "4.nl_EWAS_agnostic_quant_sex", "nl_assoc_EWAS_agnostic_quant_sex.RDS"))


## -------------------------------------------------------------------------------------------------
# Count the non-linear associations
nl_assoc_EWAS_agnostic_quant_sex %>%
  filter(spline_tot_FDR < 0.05) %>%
  nrow() # There are 70 non-linear associations

## -------------------------------------------------------------------------------------------------

# After identification of the non-linearly associated CpGs, the nominal p-values in the agnostic sex-stratified 
# EWAS (with 2% of extreme exposure concentrations removed) for the CpGs that showed non-linear associations 
# should be replaced by the nominal p-values of the non-linear associations.
# The FDR-corrected p-values will be replaced by the FDR-corrected p-values for the non-linear associations.
# This EWAS result will be used in the further analyses.

# Replace FDR-corrected p-values in the main EWAS for the CpGs that showed non-linear associations (spline total FDR-corrected p-value < 0.05) and for these CpGs replace nominal p-values with spline total p-value
nl_EWAS_agnostic_quant_sex <- full_join(
  EWAS_agnostic_quant_sex,
  select(nl_assoc_EWAS_agnostic_quant_sex, child_sex, CpG, Exposure, spline_tot_FDR, spline_tot_pval),
  by = c("CpG", "Exposure", "child_sex")) %>%
  mutate(p_value_FDR = case_when(spline_tot_FDR < 0.05 ~ spline_tot_FDR,
                                 TRUE ~ p_value_FDR),
         raw_p_value = case_when(spline_tot_FDR < 0.05 ~ spline_tot_pval,
                                 TRUE ~ raw_p_value),
         Estimate_CI = case_when(spline_tot_FDR < 0.05 ~ "Non-linear",
                                 TRUE ~ Estimate_CI)) %>%
  select(-c(spline_tot_FDR, spline_tot_pval))

# Save the result
export(nl_EWAS_agnostic_quant_sex, here(path_main, "4.nl_EWAS_agnostic_quant_sex", "nl_EWAS_agnostic_quant_sex.RDS"))


## -------------------------------------------------------------------------------------------------
# Select only significant associations
nl_EWAS_agnostic_quant_sex_res <- nl_EWAS_agnostic_quant_sex %>%
  filter(p_value_FDR < 0.05) # There were 114 associations in a sex-stratified analysis: males (68 CpGs) and females (46 CpGs)


## -------------------------------------------------------------------------------------------------

## 5. Linear agnostic EWAS on both sexes and whole population (n = 395)

# Continuous exposures
EWAS_agnostic_cont <-
  RobustRegressions::SpecificLociMethylationRegressions(
    meth_data = meth_clean,
    expo_cov = expo_cov_tech_conf,
    exposures = exposures,
    confounders = confounders,
    technical_confounders = c(tech_conf_EWAS_DMR, cell_mix_names),
    maxit = 400,
    annotation_object = IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19,
    expo_labels = expo_labels,
    path = path_main,
    file_name = "5.EWAS_agnostic/EWAS_agnostic_cont") %>%
  mutate(child_sex = "both_sexes")


## -------------------------------------------------------------------------------------------------
# # Categorical exposures
EWAS_agnostic_cat <-
  RobustRegressions::SpecificLociMethylationRegressions(
    meth_data = meth_clean,
    expo_cov = expo_cov_tech_conf,
    exposures = exposures_cat,
    confounders = confounders,
    technical_confounders = c(tech_conf_EWAS_DMR, cell_mix_names),
    maxit = 400,
    annotation_object = IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19,
    expo_labels = expo_labels_cat,
    path = path_main,
    file_name = "5.EWAS_agnostic/EWAS_agnostic_cat") %>%
  mutate(child_sex = "both_sexes")


## -------------------------------------------------------------------------------------------------
# Merge results for continuous and categorical exposures and save the result
EWAS_agnostic <- bind_rows(EWAS_agnostic_cont, EWAS_agnostic_cat)

export(EWAS_agnostic, here(path_main, "5.EWAS_agnostic", "EWAS_agnostic.RDS"))


## -------------------------------------------------------------------------------------------------

## 6. Non-linear agnostic EWAS on both sexes and whole population (n = 395)

nl_assoc_EWAS_agnostic <-
  RobustRegressions::NonLinearityCheck(meth_data = meth_clean,
                                       expo_cov = expo_cov_tech_conf,
                                       exposures = exposures,
                                       confounders = confounders,
                                       technical_confounders = c(tech_conf_EWAS_DMR, cell_mix_names),
                                       maxit = 400,
                                       path = path_main,
                                       file_name = "6.nl_EWAS_agnostic/nl_assoc_EWAS_agnostic") %>%
  mutate(child_sex = "both_sexes")

# For categorical exposures there is no need to run this analysis (as no influential values for dichotomised variables)


## -------------------------------------------------------------------------------------------------

# After identification of the non-linearly associated CpGs, the nominal p-values in the main agnostic sex-stratified 
# EWAS for the CpGs that showed non-linear associations should be replaced by the nominal p-values of the non-linear associations.
# The FDR-corrected p-values will be replaced by the FDR-corrected p-values for the non-linear associations.

# Replace FDR-corrected p-values in the main EWAS for the CpGs that showed non-linear associations (spline total FDR-corrected p-value < 0.05) and for these CpGs replace nominal p-values with spline total p-value
nl_EWAS_agnostic <- full_join(EWAS_agnostic,
                              select(nl_assoc_EWAS_agnostic, CpG, Exposure, spline_tot_FDR, spline_tot_pval),
                              by = c("CpG", "Exposure")) %>%
  mutate(p_value_FDR = case_when(spline_tot_FDR < 0.05 ~ spline_tot_FDR,
                                 TRUE ~ p_value_FDR),
         raw_p_value = case_when(spline_tot_FDR < 0.05 ~ spline_tot_pval,
                                 TRUE ~ raw_p_value),
         Estimate_CI = case_when(spline_tot_FDR < 0.05 ~ "Non-linear",
                                 TRUE ~ Estimate_CI)) %>%
  select(-c(spline_tot_FDR, spline_tot_pval))

export(nl_EWAS_agnostic, here(path_main, "6.nl_EWAS_agnostic", "nl_EWAS_agnostic.RDS"))


## -------------------------------------------------------------------------------------------------

## 7. Linear agnostic EWAS on both sexes and reduced population (n = 387)

# # Continuous exposures 
EWAS_agnostic_quant_cont <- RobustRegressions::SpecificLociMethylationRegressions(
  meth_data = meth_clean,
  expo_cov = quant_expo,
  exposures = exposures,
  confounders = confounders,
  technical_confounders = c(tech_conf_EWAS_DMR, cell_mix_names),
  maxit = 400,
  annotation_object = IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19,
  expo_labels = expo_labels,
  path = path_main,
  file_name = "7.EWAS_agnostic_quant/EWAS_agnostic_quant_cont",
) %>%
  mutate(child_sex = "both_sexes")

# For categorical exposures there is no need to run this analysis (as no influential values for dichotomised variables)


## -------------------------------------------------------------------------------------------------
# Merge results for continuous with extreme values removed and original categorical exposures (as for categorized exposures there was no need to remove influential values) and save the result
EWAS_agnostic_quant <- bind_rows(EWAS_agnostic_quant_cont, EWAS_agnostic_cat)
rm(EWAS_agnostic_quant_cont, EWAS_agnostic_cat)

export(EWAS_agnostic_quant, here(path_main, "7.EWAS_agnostic_quant", "EWAS_agnostic_quant.RDS"))


## -------------------------------------------------------------------------------------------------

## 8. Non-linear agnostic EWAS on both sexes and reduced population (n = 387)

nl_assoc_EWAS_agnostic_quant <-
  RobustRegressions::NonLinearityCheck(meth_data = meth_clean,
                                       expo_cov = quant_expo,
                                       exposures = exposures,
                                       confounders = confounders,
                                       technical_confounders = c(tech_conf_EWAS_DMR, cell_mix_names),
                                       maxit = 400,
                                       path = path_main,
                                       file_name = "8.nl_EWAS_agnostic_quant/nl_assoc_EWAS_agnostic_quant")


## -------------------------------------------------------------------------------------------------

# After identification of the non-linearly associated CpGs, the nominal p-values in the main agnostic 
# EWAS for the CpGs that showed non-linear associations should be replaced by the nominal p-values of the 
# non-linear associations.

# The FDR-corrected p-values will be replaced by the FDR-corrected p-values for the non-linear associations.

# Replace FDR-corrected p-values in the main EWAS for the CpGs that showed non-linear associations (spline total FDR-corrected p-value < 0.05) and for these CpGs replace nominal p-values with spline total p-value
nl_EWAS_agnostic_quant <- full_join(
  EWAS_agnostic_quant,
  select(nl_assoc_EWAS_agnostic_quant, CpG, Exposure, spline_tot_FDR, spline_tot_pval),
  by = c("CpG", "Exposure")) %>%
  mutate(p_value_FDR = case_when(spline_tot_FDR < 0.05 ~ spline_tot_FDR,
                                 TRUE ~ p_value_FDR),
         raw_p_value = case_when(spline_tot_FDR < 0.05 ~ spline_tot_pval,
                                 TRUE ~ raw_p_value),
         Estimate_CI = case_when(spline_tot_FDR < 0.05 ~ "Non-linear",
                                 TRUE ~ Estimate_CI)) %>%
  select(-c(spline_tot_FDR, spline_tot_pval))

# Save the result
export(nl_EWAS_agnostic_quant, here(path_main, "8.nl_EWAS_agnostic_quant", "nl_EWAS_agnostic_quant.RDS"))


## -------------------------------------------------------------------------------------------------
# Select only significant associations
nl_EWAS_agnostic_quant_res <- nl_EWAS_agnostic_quant %>%
  filter(p_value_FDR < 0.05) # There are no associations when both sexes were studied together


# As there were no associations for the agnostic approach, we selected CpGs with nominal p-value <0.001 to run 
# further enrichment analyses on these CpGs
# Save CpGs list for further enrichment analyses

## -------------------------------------------------------------------------------------------------

nl_EWAS_agnostic_quant_res_0.001 <- nl_EWAS_agnostic_quant_res %>% 
  filter(raw_p_value < 0.001)

export(nl_EWAS_agnostic_quant_res_0.001, here(path_main, "nl_EWAS_agnostic_quant_res_0.001.RDS"))

# Count the number of CpGs with nominal p-value <0.001
nl_EWAS_agnostic_quant_res_0.001 %>% 
  group_by(child_sex, Exposure) %>%
  summarise(n_CpGs = n())

# Save genes encompassed by the CpGs with the nominal p-value <0.001 for further enrichment analyses
genes_CpGs_agnostic_quant_res_0.001 <- nl_EWAS_agnostic_quant_res %>% 
  filter(raw_p_value < 0.001) |> 
  GenesInDMRsSEPAGES(var_to_distinct = c("Gene", "child_sex", "Exposure"))

export(genes_CpGs_agnostic_quant_res_0.001, here(path_main, "genes_CpGs_agnostic_quant_res_0.001.RDS"))

# Count the number of genes encompassed by the CpGs with nominal p-value <0.001
genes_CpGs_agnostic_quant_res_0.001 %>% 
  group_by(child_sex, Exposure) %>% 
  summarise(n_genes = n()) |> 
  arrange(n_genes)


## -------------------------------------------------------------------------------------------------

## Summary of the results for the whole and reduced populations, for both sexes and sex-stratified analyses

# Display significant associations for the whole (n = 395) and reduced (n = 387) populations

# Agnostic sex-stratified, whole population
nl_EWAS_agnostic_sex_res <- nl_EWAS_agnostic_sex %>%
  filter(p_value_FDR < 0.05)

# Agnostic, both sexes, whole population
nl_EWAS_agnostic_res <- nl_EWAS_agnostic %>%
  filter(p_value_FDR < 0.05)

# Agnostic sex-stratified, reduced population
nl_EWAS_agnostic_quant_sex_res <- nl_EWAS_agnostic_quant_sex %>%
  filter(p_value_FDR < 0.05)

# Agnostic, both sexes, reduced population
nl_EWAS_agnostic_quant_res <- nl_EWAS_agnostic_quant %>%
  filter(p_value_FDR < 0.05)


## -------------------------------------------------------------------------------------------------

## Check which CpGs from the exploratory study (out of 114 identified CpGs) were present in the 450k study

114 - length(setdiff(pull(nl_EWAS_quant_res, CpG), unique(pull(ungroup(EWAS_EDEN), CpG))))
# 55 out of 114 CpGs were present in 450k study - this may be the partial explanation of the low overlap between the studies


## -------------------------------------------------------------------------------------------------

## DMR study

## 9. Agnostic sex-stratified DMR on linear and non-linear associations on reduced population (n = 387)

# This DMR analysis will be based on the associations identified for the reduced dataset (1% upper and 1% lower 
# extreme exposure values removed) for a sex-stratified analysis.
# Standard p-value corrections will be applied (the same setting as in EDEN).

# Create bed files for each exposure so they can be processed by the comb-p
for (sex in c("Male", "Female")) {
  
  DMR::bedCreator(annotated_EWAS_result = filter(nl_EWAS_agnostic_quant_sex, child_sex == sex),
                  prefix = paste0("agnostic_quant_", sex),
                  path = paste0(path_main, "/9.DMR_agnostic_quant_sex/bed_files"))
}

# Setup for DMR analysis: since it is an exploratory study, p-value for starting a DMR would be kept as in 
# EDEN at 0.001 level.

# After running the combp analysis externally, copy output `*.regions-t.bed` files to /comb_p_results folder.

## -------------------------------------------------------------------------------------------------
# Clean and annotate the DMR results
DMR_agnostic_quant_sex <- data.frame()

# Clean and unify the names of the exposures (variable levels)
expo_levels <- c(exposures, exposures_cat)

# Remove exposures that were not identified as associate4d with DMRs
expo_levels_females <- expo_levels[!expo_levels %in% c("aver_log2_BPA", "aver_log2_BP3")]
expo_levels_males <- expo_levels[!expo_levels %in% c("aver_log2_oxoMINCH", "aver_log2_BP3", "aver_log2_MMCHP")]

# Clean and unify the names of the exposures (variable labels)
expo_lab <- c(expo_labels, expo_labels_cat)
expo_labels_females <- expo_lab[!expo_lab %in% c("Bisphenol A (BPA)", "Benzophenone-3 (BP-3)")]
expo_labels_males <- expo_lab[!expo_lab %in% c("2-(((4-methyl-7-oxooctyl)oxy)carbonyl) cyclohexanecarboxylic acid (oxo-MINCH)", "Benzophenone-3 (BP-3)", "Mono-2-carboxymethylhexyl phthalate (MMCHP)" )]

for (sex in c("Male", "Female")) {
  
  if (sex == "Male") {
    
    expo_lv <- expo_levels_males
    expo_lb <- expo_labels_males
    
  } else {
    
    expo_lv <- expo_levels_females
    expo_lb <- expo_labels_females
  }
  
  annot_data <- filter(nl_EWAS_agnostic_quant_sex, child_sex == sex)
  DMR_sex <- DMR::CleanDMRResults(
    annotated_EWAS_result = annot_data,
    expo_levels = expo_lv,
    expo_labels = expo_lb,
    path_to_combp = paste0(path_main, paste0("/9.DMR_agnostic_quant_sex/comb_p_results/", sex)),
    annotation_object = IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19,
    sepages = TRUE) %>%
    mutate(child_sex = sex)
  
  DMR_agnostic_quant_sex <- bind_rows(DMR_agnostic_quant_sex, DMR_sex)
}

# Save the results
export(DMR_agnostic_quant_sex, here(path_main, "9.DMR_agnostic_quant_sex", "DMR_agnostic_quant_sex.RDS"))


## -------------------------------------------------------------------------------------------------

## 10. Agnostic DMR on linear and non-linear associations on both sexes on reduced population (n = 387)

# Categorized variables (BUPA and BPS) are ignored as no individuals were removed

# Create bed files for each exposure so they can be processed by the comb-p
DMR::bedCreator(annotated_EWAS_result = nl_EWAS_agnostic_quant,
                prefix = "agnostic_quant_both_sexes",
                path = paste0(path_main, "/10.DMR_agnostic_quant/bed_files"))

# Setup for DMR analysis: since it is an exploratory study, p-value for starting a DMR would be kept as in 
# EDEN at 0.001 level.

# After running the combp analysis, copy output `*.regions-t.bed` files to /comb_p_results folder.


## -------------------------------------------------------------------------------------------------
# DMR results cleaning
DMR_agnostic_quant <- DMR::CleanDMRResults(
  annotated_EWAS_result = nl_EWAS_agnostic_quant,
  expo_levels = c(exposures, exposures_cat),
  expo_labels = c(expo_labels, expo_labels_cat),
  path_to_combp = paste0(path_main, "/10.DMR_agnostic_quant/comb_p_results"),
  annotation_object = IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19,
  sepages = TRUE) %>%
  mutate(child_sex = "both_sexes")

# Save the results
export(DMR_agnostic_quant, here(path_main, "10.DMR_agnostic_quant", "DMR_agnostic_quant.RDS"))


## -------------------------------------------------------------------------------------------------
# Merge both sexes and sex-stratified analyses
DMR_agnostic_quant_res <- rbind(DMR_agnostic_quant, DMR_agnostic_quant_sex) %>%
  select(child_sex, Exposure, everything()) %>%
  arrange(child_sex, Exposure, Gene)

export(DMR_agnostic_quant_res, here(path_main, "DMR_agnostic_quant_res.RDS"))


## -------------------------------------------------------------------------------------------------
# Clean the names
DMR_agnostic_quant_res <- DMR_agnostic_quant_res %>%
  mutate(Expo_family = ifelse(str_detect(Exposure, "BPA|BPS|BP-3|BUPA|ETPA|MEPA|PRPA"), "Phenol", "Phthalate"),
         child_sex = ifelse(child_sex == "both_sexes", "Females and males", child_sex),
         child_sex = factor(child_sex,
                            levels = c("Females and males", "Female", "Male"),
                            labels = c("Females and males", "Females", "Males"))) %>%
  select(child_sex, Expo_family, Exposure, everything())

# Select DMRs with at least 5 CpGs
DMR_agnostic_quant_res_5 <- DMR_agnostic_quant_res %>%
  filter(`No. of probes` >= 5)

# Count the number of DMRs with at least 5 CpGs
DMR_agnostic_quant_res_5_tab <- DMR_agnostic_quant_res_5 %>%
  group_by(child_sex, Expo_family, Exposure) %>%
  summarise(n = n()) %>%
  arrange(child_sex, desc(n))

# List genes encompassed by the DMRs with at least 5 CpGs
DMR_agnostic_quant_res_5 <- DMR_agnostic_quant_res_5 %>%
  mutate(Gene = ifelse(Gene == "", "Unknown", Gene))

DMR_agnostic_quant_res_5_genes <- DMR_agnostic_quant_res_5 %>%
  GenesInDMRsSEPAGES(var_to_distinct = c("Gene", "child_sex", "Exposure")) %>%
  distinct(Gene, .keep_all = TRUE) # males and females 23 genes (188 CpGs); males 32 genes; females 28 genes (603 CpGs in total)


## -------------------------------------------------------------------------------------------------

## Check the correlation between estimates within the overlapping DMRs between the whole population, males and females

# List CpGs within DMRs
# Sex-stratified analysis
CpGs_in_DMRs_5probes_agnostic_sex <- data.frame()

for (sex in c("Female", "Male")) {
  
  CpGs_in_DMRs <- CpGsInDMRsSEPAGES(meth_data = meth_clean,
                                    probes = 5,
                                    annotated_EWAS_result = filter(nl_EWAS_agnostic_quant_sex, child_sex == sex),
                                    path_to_combp = paste0(path_main, "/9.DMR_agnostic_quant_sex/comb_p_results/", sex),
                                    path = here(path_main, "9.DMR_agnostic_quant_sex"),
                                    file_name = paste0("CpGs_in_DMRs_5probes_agnostic_sex_", sex)) %>%
    mutate(child_sex = sex)
  
  CpGs_in_DMRs_5probes_agnostic_sex <- rbind(CpGs_in_DMRs_5probes_agnostic_sex, CpGs_in_DMRs)
}

# Save the results
export(CpGs_in_DMRs_5probes_agnostic_sex, here(path_main, "9.DMR_agnostic_quant_sex", "CpGs_in_DMRs_5probes_agnostic_sex.csv"))


# Both sexes analysis
CpGs_in_DMRs_5probes_agnostic <- CpGsInDMRsSEPAGES(meth_data = meth_clean,
                                                   probes = 5,
                                                   annotated_EWAS_result = nl_EWAS_agnostic_quant,
                                                   path_to_combp = paste0(path_main, "/10.DMR_agnostic_quant/comb_p_results"),
                                                   path = here(path_main, "10.DMR_agnostic_quant"),
                                                   file_name = "CpGs_in_DMRs_5probes_agnostic") %>%
  mutate(child_sex = "Female and male")

# Save the results
export(CpGs_in_DMRs_5probes_agnostic, here(path_main, "10.DMR_agnostic_quant/CpGs_in_DMRs_5probes_agnostic.csv"))


## -------------------------------------------------------------------------------------------------
# Merge results listing probes within DMRs
CpGs_in_DMRs_5probes <- rbind(CpGs_in_DMRs_5probes_agnostic_sex, CpGs_in_DMRs_5probes_agnostic) %>%
  
  # Clean exposure names
  mutate(Expo_family = ifelse(str_detect(Exposure, "BPA|BPS|BP3|BUPA|ETPA|MEPA|PRPA"), "Phenol", "Phthalate"),
         Expo_subfamily = case_when(str_detect(Exposure_name, "MINCH|DINCH") ~ "DINCH metabolites",
                                    str_detect(Exposure, "BP3") ~ "Benzophenone-3",
                                    str_detect(Exposure_name, "MiNP|DiNP") ~ "DiNP metabolites",
                                    str_detect(Exposure, "BPA|BPS") ~ "Bisphenols",
                                    str_detect(Exposure, "BUPA|ETPA|MEPA|PRPA") ~ "Parabens",
                                    str_detect(Exposure, "MECPP|MEHHP|MEHP|MEOHP|MMCHP|DEHP") ~ "DEHP metabolites",
                                    TRUE ~ "Other phthalates"))

# Save the results
export(CpGs_in_DMRs_5probes, here(path_main, "CpGs_in_DMRs_5probes_agnostic_annotated.csv"))


## -------------------------------------------------------------------------------------------------

## Check if the CpGs identified in exploratory and candidate appraoch are present on lists available 
# in Delahaye et al. and Deyssenroth et al. 

# Load the list of eQTMs identified by Delahaye et al. and Deyssenroth et al. 
list_cpgs <- read_excel_allsheets(here("data", "CpGs_for_eQTM.xlsx"))

# Delahaye
intersect(list_cpgs$exploratory_EWAS$CpG, list_cpgs$Delahaye$CpG) # 0 with Delahaye
intersecting_cpgs_delhaye <- intersect(list_cpgs$exploratory_DMR$CpG, list_cpgs$Delahaye$CpG) # 12 with Delahaye
intersecting_cpgs_delhaye
intersect(list_cpgs$candidate_CpGs$CpG, list_cpgs$Delahaye$CpG) # 0 with Delahaye

# Deyssenroth
# Lucile Broséus did not find any intersection

# For the DMR study, the following CpGs were overlapping with Delahaye et al.:
  
# [1] "cg05881566" "cg06137032" "cg07119472" "cg07486017" "cg18610205" "cg05952543"
# [7] "cg09437135" "cg16131766" "cg20769842" "cg20792895" "cg00446235" "cg13466383"

# No CpGs were intersecting with Deyssenroth et al.


## -------------------------------------------------------------------------------------------------

## Check the correlation between estimates within the overlapping DMRs between phenols and phthalates

CpG_in_DMRs_5probes_males <- CpGs_in_DMRs_5probes_agnostic_sex %>%
  filter(child_sex == "Male") %>%
  select(CpG, Exposure, Gene, Estimate, raw_p_value, p_value_FDR) %>%
  left_join(filter(select(nl_EWAS_agnostic_quant_sex, Exposure, CpG, Estimate, child_sex, raw_p_value, p_value_FDR), child_sex == "Female"), by = c("Exposure", "CpG")) %>%
  select(-child_sex) %>%
  rename(Est_males = Estimate.x,
         Est_females = Estimate.y,
         raw_p_value_males = raw_p_value.x,
         raw_p_value_females = raw_p_value.y,
         p_value_FDR_males = p_value_FDR.x,
         p_value_FDR_females = p_value_FDR.y)


# Plot the CpGs identified within DMRs for males and females
tiff(here(path_suppl, "figures", "5DMR_est_males_females.tiff"),
     width = 800, height = 600, res = 150)
plot(CpG_in_DMRs_5probes_males$Est_males, 
     CpG_in_DMRs_5probes_males$Est_females,
     xlab = "Beta estimate males",
     ylab = "Beta estimate females",
     main = paste0("Significant DMRs (>=5 probes) in males \ncorr. = ", round(cor(CpG_in_DMRs_5probes_males$Est_males, 
                                                                                  CpG_in_DMRs_5probes_males$Est_females), 1)))
dev.off()

CpG_in_DMRs_5probes_females <- CpGs_in_DMRs_5probes_agnostic_sex %>%
  filter(child_sex == "Female") %>%
  select(CpG, Exposure, Gene, Estimate, raw_p_value, p_value_FDR) %>%
  left_join(filter(select(nl_EWAS_agnostic_quant_sex, Exposure, CpG, Estimate, child_sex, raw_p_value, p_value_FDR), child_sex == "Male"), by = c("Exposure", "CpG")) %>%
  select(-child_sex) %>%
  rename(Est_females = Estimate.x,
         Est_males = Estimate.y,
         raw_p_value_females = raw_p_value.x,
         raw_p_value_males = raw_p_value.y,
         p_value_FDR_females = p_value_FDR.x,
         p_value_FDR_males = p_value_FDR.y)

tiff(here(path_suppl, "figures", "5DMR_est_females_males.tiff"),
     width = 800, height = 600, res = 150)
plot(CpG_in_DMRs_5probes_females$Est_females, 
     CpG_in_DMRs_5probes_females$Est_males,
     xlab = "Beta estimate females",
     ylab = "Beta estimate males",
     main = paste0("Significant DMRs (>=5 probes) in females \n corr. = ", round(cor(CpG_in_DMRs_5probes_females$Est_females, 
                                                                                     CpG_in_DMRs_5probes_females$Est_males), 1)))
dev.off()

# Export regression estimates comparison between males and females for CpGs identified within DMRs
estimates_DMRs_5 = list("reference_males" = CpG_in_DMRs_5probes_males, "reference_females" = CpG_in_DMRs_5probes_females)
export(estimates_DMRs_5, here(path_main, "estimates_DMRs_5.RDS"))


## -------------------------------------------------------------------------------------------------

## Imprinted genes

# Check which genes idenitified within DMRs are imprinted
DMR_agnostic_quant_res_5_genes[DMR_agnostic_quant_res_5_genes %in% imprinted_genes$Gene_id] %>% 
  sort()
# APC, GNAS, GNASAS, MIMT1, MKRN3, PEG3, PEG10, PPAP2C, SGCE, ZIM2 

# Save imprinted genes list
DMR_agnostic_quant_res_5 %>% 
  filter(str_detect(Gene, paste(imprinted_genes$Gene_id, collapse = "|"))) %>% 
  export(here(path_main, "imp_genes_5probes.csv"))


## -------------------------------------------------------------------------------------------------
# Save lists of genes identified within DMRs for phenols and phthtalates separately
# Phthalates
DMR5_agnostic_quant_res_phth %>% 
  select(child_sex:Gene) %>% 
  group_by(child_sex) %>% 
  Helpers::CpGsInGenes(var_to_distinct = "Gene") %>% 
  arrange(Gene, child_sex, Exposure) %>% 
  export(here(path_main, "DMR5_per_expo_genes_phth.csv"))

# Phenols
DMR5_agnostic_quant_res_phenols %>% 
  select(child_sex:Gene) %>% 
  group_by(child_sex) %>% 
  Helpers::CpGsInGenes(var_to_distinct = "Gene") %>% 
  arrange(Gene, child_sex, Exposure) %>% 
  export(here(path_main, "DMR5_per_expo_genes_phen.csv"))

#  Split genes within the identified DMRs so each gene occupies one row
# Phthalates
DMR5_agnostic_quant_res_phth_genes <- DMR5_agnostic_quant_res_phth %>% 
  group_by(child_sex) %>% 
  Helpers::CpGsInGenes(var_to_distinct = "Gene") %>% 
  pull(Gene)

# Phenols
DMR5_agnostic_quant_res_phen_genes <- DMR5_agnostic_quant_res_phenols %>% 
  group_by(child_sex) %>% 
  Helpers::CpGsInGenes(var_to_distinct = "Gene") %>% 
  pull(Gene)

# Combine phenols and phthalates genes lists
DMR5_agnostic_quant_res_genes <- c(DMR5_agnostic_quant_res_phth_genes, DMR5_agnostic_quant_res_phen_genes) %>% 
  unique()


## -------------------------------------------------------------------------------------------------

## II. Concept-driven approach - overlapping compounds and CpGs

# Only CpGs and compounds overlapping between SEPAGES and EDEN cohorts will be considered. Only DMR analysis 
# will be run with less restrictive threshold to start a region. A region will be considered replicated if it 
# encompasses genes identified as associated with phenols/phthalates in our previous studies (Jedynak et al. 2021, 2022)

## 11. Sex-stratified DMR on linear and non-linear associations on reduced population (n = 387) and compounds (15) and 
## CpGs (337,722) overlapping with EDEN (note: DEHP sum, in contrary to the agnostic approach, is composed of 4 DEHP 
## metabolites)

# First run EWAS for overlapping compounds and CpGs (n = 387, 337,722 CpGs). These results will be used in further 
# DMR analysis but not discussed.

EWAS_overlap_quant_sex <- data.frame()

for (sex in c("Male", "Female")) {
  
  expo_cov_tech_conf_sex <- filter(quant_expo, child_sex == sex) %>%
    select(id, aver_log2_DEHP_red, mother_age:Trophoblasts)
  
  meth_clean_restr_sex <- meth_clean[rownames(meth_clean) %in% expo_cov_tech_conf_sex$id, colnames(meth_clean) %in% unique(EWAS_EDEN$CpG)]
  
  EWAS_quant_DEHP_red <- RobustRegressions::SpecificLociMethylationRegressions(
    meth_data = meth_clean_restr_sex,
    expo_cov = expo_cov_tech_conf_sex,
    exposures = exposures_overlap,
    confounders = confounders[confounders != "child_sex"],
    technical_confounders = c(tech_conf_EWAS_DMR, cell_mix_names),
    maxit = 400,
    annotation_object = IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19,
    expo_labels = expo_labels_overlap,
    path = path_main,
    file_name = paste0("EWAS_overlap_quant_", sex),
  ) %>%
    mutate(child_sex = sex)
  
  EWAS_overlap_quant_sex <- bind_rows(EWAS_overlap_quant_sex, EWAS_quant_DEHP_red)
}

export(EWAS_overlap_quant_sex, here(path_main, "11.DMR_overlap_quant_sex", "EWAS_overlap_quant_sex.RDS"))


## -------------------------------------------------------------------------------------------------

## Check if there are any non-linear associations

nl_assoc_EWAS_overlap_quant_sex <- data.frame()

for (sex in c("Male", "Female")) {
  
  expo_cov_tech_conf_sex <- filter(quant_expo, child_sex == sex) %>%
    select(id, aver_log2_DEHP_red, mother_age:Trophoblasts)
  meth_clean_restr_sex <- meth_clean[rownames(meth_clean) %in% expo_cov_tech_conf_sex$id, colnames(meth_clean) %in% unique(EWAS_EDEN$CpG)]
  
  nl_EWAS_quant <-
    RobustRegressions::NonLinearityCheck(meth_data = meth_clean_restr_sex,
                                         expo_cov = expo_cov_tech_conf_sex,
                                         exposures = exposures_overlap,
                                         confounders = confounders[confounders != "child_sex"],
                                         technical_confounders = c(tech_conf_EWAS_DMR, cell_mix_names),
                                         maxit = 400,
                                         path = path_main,
                                         file_name = paste0("nl_assoc_EWAS_quant_", sex)) %>%
    mutate(child_sex = sex)
  
  nl_assoc_EWAS_overlap_quant_sex <- bind_rows(nl_assoc_EWAS_overlap_quant_sex, nl_EWAS_quant)
} # There were no non-linear associations detected

export(nl_assoc_EWAS_overlap_quant_sex, here(path_main, "11.DMR_overlap_quant_sex", "nl_assoc_EWAS_overlap_quant_sex.RDS"))


## -------------------------------------------------------------------------------------------------

## First select the CpGs and compounds overlapping between SEPAGES and EDEN from the agnostic EWAS results.

nl_EWAS_overlap_quant_sex <- nl_EWAS_agnostic_quant_sex %>%
  filter(CpG %in% unique(EWAS_EDEN$CpG)) %>%
  bind_rows(EWAS_overlap_quant_sex) %>%
  filter(Exposure %in% exposures_overlap)
# Adjusting of the FDR-corrected p-values for the new number of CpGs is not necessary, as only nominal p-values will be used

export(nl_EWAS_overlap_quant_sex, here(path_main, "11.DMR_overlap_quant_sex", "nl_EWAS_overlap_quant_sex.RDS"))


## -------------------------------------------------------------------------------------------------

# Setup for DMR analysis: since it is a replication study, p-value for starting a DMR would be increased from 0.001 to 0.05.

# After running the combp analysis, copy output `*.regions-t.bed` files to /comb_p_results folder.

# Create bed files for each exposure so they can be processed by the comb-p
for (sex in c("Male", "Female")) {
  
  DMR::bedCreator(annotated_EWAS_result = filter(nl_EWAS_overlap_quant_sex, child_sex == sex),
                  prefix = paste0("overlap_quant_sex_", sex),
                  path = paste0(path_main, "/11.DMR_overlap_quant_sex/bed_files"))
}


## -------------------------------------------------------------------------------------------------
# Sex-stratified DMR results cleaning
# Change name of the DEHP exposure so it fits the naming pattern from DMR result cleaning function
nl_EWAS_overlap_quant_sex <- nl_EWAS_overlap_quant_sex %>%
  mutate(Exposure = ifelse(Exposure == "aver_log2_DEHP_red", "aver_log2_DEHP", Exposure))

exposures_overlap_mod <- replace(exposures_overlap, exposures_overlap == "aver_log2_DEHP_red", "aver_log2_DEHP")

DMR_overlap_quant_sex <- data.frame()

for (sex in c("Male", "Female")) {
  
  DMR_quant_sex <- DMR::CleanDMRResults(
    annotated_EWAS_result = filter(nl_EWAS_overlap_quant_sex, child_sex == sex),
    expo_levels = exposures_overlap_mod,
    expo_labels = expo_labels_overlap,
    path_to_combp = paste0(path_main, paste0("/11.DMR_overlap_quant_sex/comb_p_results/", sex)),
    annotation_object = IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19,
    sepages = TRUE) %>%
    mutate(child_sex = sex)
  
  DMR_overlap_quant_sex <- bind_rows(DMR_overlap_quant_sex, DMR_quant_sex)
}

export(DMR_overlap_quant_sex, here(path_main, "11.DMR_overlap_quant_sex", "DMR_overlap_quant_sex.RDS"))


## -------------------------------------------------------------------------------------------------

## 12. DMR on both sexes on linear and non-linear associations on reduced population (n = 387) and compounds (15) 
## and CpGs (337,722) overlapping with EDEN

# First run EWAS for overlapping compounds and CpGs (n = 387, 337,722 CpGs). These results will be used in 
# further DMR analysis but not discussed.

meth_clean_restr <- meth_clean[, colnames(meth_clean) %in% unique(EWAS_EDEN$CpG)]
export(meth_clean_restr, here(path_main, "meth_clean_restr.RDS"))

EWAS_overlap_quant <- RobustRegressions::SpecificLociMethylationRegressions(
  meth_data = meth_clean_restr,
  expo_cov = quant_expo,
  exposures = exposures_overlap,
  confounders = confounders,
  technical_confounders = c(tech_conf_EWAS_DMR, cell_mix_names),
  maxit = 400,
  annotation_object = IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19,
  expo_labels = expo_labels_overlap,
  path = path_main,
  file_name = "12.DMR_overlap_quant/EWAS_overlap_quant",
) %>%
  mutate(child_sex = "both_sexes")


## -------------------------------------------------------------------------------------------------

## Check if there are any non-linear associations

nl_assoc_EWAS_overlap_quant <-
  RobustRegressions::NonLinearityCheck(meth_data = meth_clean_restr,
                                       expo_cov = quant_expo,
                                       exposures = exposures_overlap,
                                       confounders = confounders,
                                       technical_confounders = c(tech_conf_EWAS_DMR, cell_mix_names),
                                       maxit = 400,
                                       path = path_main,
                                       file_name = "12.DMR_overlap_quant/nl_assoc_EWAS_overlap_quant") # There were no non-linear associations


## -------------------------------------------------------------------------------------------------

## First select the CpGs and compounds overlapping between SEPAGES and EDEN from the agnostic EWAS results.

nl_EWAS_overlap_quant <- nl_EWAS_agnostic_quant %>%
  filter(CpG %in% unique(EWAS_EDEN$CpG)) %>%
  bind_rows(EWAS_overlap_quant) %>%
  filter(Exposure %in% exposures_overlap)

export(nl_EWAS_overlap_quant, here(path_main, "12.DMR_overlap_quant", "nl_EWAS_overlap_quant.RDS"))


## -------------------------------------------------------------------------------------------------

# Setup for DMR analysis: since it is a replication study, p-value for starting a DMR would be increased from 0.001 to 0.05.

# After running the combp analysis, copy output `*.regions-t.bed` files to /comb_p_results folder.

# Create bed files for each exposure so they can be processed by the comb-p
DMR::bedCreator(annotated_EWAS_result = nl_EWAS_overlap_quant,
                prefix = "overlap_quant",
                path = paste0(path_main, "/12.DMR_overlap_quant/bed_files"))


## -------------------------------------------------------------------------------------------------
# Change name of the DEHP exposure so it fits the naming pattern from DMR result cleaning function
nl_EWAS_overlap_quant <- nl_EWAS_overlap_quant %>%
  mutate(Exposure = ifelse(Exposure == "aver_log2_DEHP_red", "aver_log2_DEHP", Exposure))

# DMR results cleaning
DMR_overlap_quant <- DMR::CleanDMRResults(
  annotated_EWAS_result = nl_EWAS_overlap_quant,
  expo_levels = exposures_overlap_mod,
  expo_labels = expo_labels_overlap,
  path_to_combp = paste0(path_main, "/12.DMR_overlap_quant/comb_p_results"),
  annotation_object = IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19,
  sepages = TRUE) %>%
  mutate(child_sex = "both_sexes")

export(DMR_overlap_quant, here(path_main, "12.DMR_overlap_quant", "DMR_overlap_quant.RDS"))


## -------------------------------------------------------------------------------------------------
# Merge both sexes and sex-stratified analyses
DMR_overlap_quant_res <- rbind(DMR_overlap_quant, DMR_overlap_quant_sex) %>%
  select(child_sex, Exposure, everything()) %>%
  arrange(child_sex, Exposure, Gene)

export(DMR_overlap_quant_res, here(path_main, "DMR_overlap_quant_res.RDS"))


## -------------------------------------------------------------------------------------------------

## ## Genes overlapping between EDEN and SEPAGES studies

# Display how many single genes (splitting e.g. SGCE/PEG10 in two) are in the result and how many replicate with EDEN study
DMR_overlap_res <- DMR_overlap_quant_res %>%
  separate(Exposure, c("Exp", "Exposure"), sep = "\\(") %>% 
  mutate(Exposure = str_replace_all(DMR_overlap_res$Exposure, "\\)", ""))

DMR_overlap_res %>%
  filter(child_sex == "both_sexes") %>% 
  Helpers::CpGsInGenes(var_to_distinct = "Gene") %>% 
  select(Exposure, Gene) %>% # 154 genes, 185 DMRs
  inner_join(genes_DMR_2probes, by = c("Exposure", "Gene")) # 0 genes overlapping

DMR_overlap_res %>%
  filter(child_sex == "Female") %>% 
  Helpers::CpGsInGenes(var_to_distinct = "Gene") %>% 
  select(Exposure, Gene) %>% # 263 genes, 394 DMRs
  inner_join(genes_DMR_2probes, by = c("Exposure", "Gene")) # 0 genes overlapping

DMR_overlap_res %>%
  filter(child_sex == "Male") %>% 
  Helpers::CpGsInGenes(var_to_distinct = "Gene") %>% 
  select(Exposure, Gene) %>% # 339 genes, 428 DMRs
  inner_join(genes_DMR_2probes, by = c("Exposure", "Gene")) # 0 genes overlapping

# There were no overlapping DMRs.

## -------------------------------------------------------------------------------------------------

## **III. Concept-driven approach - CpGs associated with phenol and phthalate exposure in EDEN**

# There were 24 CpGs associated with phenol and phthalate exposure in males in EDEN.
# 20 of them were available in SEPAGES.
# We will check if any of these CpGs shows association in SEPAGES (on the level of 0.05 for the nominal p-value)

CpGs_overlap_res <- nl_EWAS_agnostic_quant_res %>%
  
  # Select variables of interest
  select(child_sex, Exposure_name, CpG, Estimate, Estimate_CI, raw_p_value) %>%
  
  # Add suffix _SEPAGES to the columns of interest to distinguish from the EDEN study
  rename_at(vars(-c(CpG, child_sex, Exposure_name)), function(x) x = paste0(x, "_SEPAGES"))

# Select CpGs associated with exposure in EDEN (23 for phenols and 1 for phthalates)
EWAS_EDEN_sel <- EWAS_EDEN %>%
  
  # Restrict EWAS results in EDEN to compounds common to EDEN and SEPAGES
  filter(Exposure %nin% c("DCP_24_log2", "DCP_25_log2", "Sum_DCPs_log2", "TCS_log2", "MCPP_log2", "MCOP_log2", "MCNP_log2")) %>%
  droplevels() %>%
  
  # Select significant CpGs
  filter(p_value_FDR < 0.05) %>%
  ungroup() %>%
  
  mutate(Expo_family = ifelse(Exposure_name %in% c("Bisphenol A", "Benzophenone-3", "Butylparaben", "Ethylparaben", "Methylparaben", "Propylparaben"), "Phenol", "Phthalate")) %>%
  
  # Select variables of interest
  select(Expo_family, Exposure_name, CpG:Location_in_gene, Estimate, Estimate_CI, raw_p_value) %>%
  
  # Replace linear associations with non-linear
  mutate(Estimate_CI = case_when(CpG %in% linear_CpGs_EDEN ~ Estimate_CI,
                                 TRUE ~ "Non-linear")) %>%
  
  # Add suffix _EDEN to the columns of interest to distinguish from the EDEN study
  rename_at(vars(c(Estimate_CI, raw_p_value)), function(x) x = paste0(x, "_EDEN"))

# Merge SEPAGES and EDEN datasets for comparison
EWAS_EDEN_SEPAGES <- EWAS_EDEN_sel %>%
  
  # Merge with SEPAGES EWAS results
  merge(CpGs_overlap_res, by = c("CpG", "Exposure_name")) %>%
  
  # Order the columns and arrange by gene
  select(child_sex, Expo_family, Exposure_name, CpG, Chr:Location_in_gene, contains("SEPAGES"), contains("EDEN")) %>%
  arrange(child_sex, Expo_family, Exposure_name, Gene)

# There are 20 CpGs merged with SEPAGES data (males, females, both, so the same 20 CpGs are repeated 3 times) 

# Select overlapping CpGs from SEPAGES with nominal p-value < 0.05
EWAS_EDEN_SEPAGES %>% 
  filter(raw_p_value_SEPAGES < 0.05) %>% 
  arrange(Gene) # There are 3 overlapping CpGs and 1 gene: LGALS8 (protein-coding, imprinted)
# Male Ethylparaben	cg13444964 LGALS8
# Female Propylparaben	cg27455890 Unknown
# both_sexes Propylparaben	cg27455890 Unknown


## -------------------------------------------------------------------------------------------------

## 13. Mixture analysis - quantile g-computation

# qgcomp is a package to implement g-computation for analyzing the effects of exposure mixtures. 
# Quantile g-computation yields estimates of the effect of increasing all exposures by one quantile, 
# simultaneously. This, it estimates a “mixture effect” useful in the study of exposure mixtures such as air pollution, 
# diet, and water contamination.

## Both sexes

# Transform beta-values to m-values so they approach normality
meth_clean_M <- beta2m(meth_clean)
export(meth_clean_M, here("data", "meth_clean_M.RDS")) # 395 x 752,577


## -------------------------------------------------------------------------------------------------
# Run the g-comp analysis
results_gq <- data.frame()

system.time(for (i in seq_along(meth_clean_M)) {
  
  y <- data.frame(y = meth_clean_M[, i])
  data_comp <- bind_cols(y,
                         select(quant_expo, all_of(c(exposures[!exposures %in% c("aver_log2_DEHP", "aver_log2_DiNP", "aver_log2_DINCH")],
                                                     confounders, 
                                                     tech_conf_EWAS_DMR, 
                                                     cell_mix_names)))) |> 
    na.omit()
  
  
  qc.fit <- qgcomp.glm.noboot(f = y ~ ., 
                              data = data_comp, 
                              expnms = exposures[!exposures %in% c("aver_log2_DEHP", "aver_log2_DiNP", "aver_log2_DINCH")], 
                              family = gaussian())
  
  res_qg_fit <- data.frame(CpG = colnames(meth_clean_M)[i], 
                           est_psi1 = unname(qc.fit$psi),
                           CI_psi1_lower = unname(qc.fit$ci)[1],
                           CI_psi1_upper = unname(qc.fit$ci)[2],
                           p_val_psi1 = unname(qc.fit$pval[2]))
  
  results_gq <- bind_rows(results_gq, res_qg_fit)
})


# Save the results
export(results_gq, here(path_main, "13.Mixture_analysis", "results_mixture_both_sexes.RDS"))


## -------------------------------------------------------------------------------------------------

## Add FDR correction of p-value\

results_gq <- results_gq |> 
  mutate(FDR_p_value = p.adjust(p_val_psi1, method = "BH"))

min(results_gq$FDR_p_value) # 0.098 (no significant effect of the mixture)


## -------------------------------------------------------------------------------------------------

## Add CpGs features needed to create BED files for further DMR analysis

results_gq <- nl_EWAS_agnostic_quant |> 
  select(CpG, Chr, Position) |> 
  distinct(CpG, .keep_all = TRUE) |> 
  merge(results_gq, by = "CpG") |> 
  mutate(Exposure = "mixture",
         raw_p_value = p_val_psi1) |> 
  select(Exposure, CpG:Position, raw_p_value)


## -------------------------------------------------------------------------------------------------
# Create bed files
bedCreator(annotated_EWAS_result = results_gq,
           prefix = "quant_both_sexes",
           path = paste0(path_main, "/13.Mixture_analysis/bed_files"))


## -------------------------------------------------------------------------------------------------

## Sex-stratified analysis

results_gq <- 
  results_gq_sex <- data.frame()
for (sex in c("Male", "Female")) {
  
  quant_expo_sex <- filter(quant_expo, child_sex == sex)
  meth_clean_M_sex <- meth_clean_M[quant_expo_sex$id, ]
  
  results_gq <- data.frame()
  for (i in seq_along(meth_clean_M)) { 
    
    y <- data.frame(y = meth_clean_M_sex[, i])
    
    data_comp <- bind_cols(y,
                           select(quant_expo_sex, all_of(c(exposures[!exposures %in% c("aver_log2_DEHP", "aver_log2_DiNP", "aver_log2_DINCH")],
                                                           confounders[confounders != "child_sex"], 
                                                           tech_conf_EWAS_DMR, 
                                                           cell_mix_names)))) |> 
      na.omit()
    
    
    qc.fit <- qgcomp.glm.noboot(f = y ~ ., 
                                data = data_comp, 
                                expnms = exposures[!exposures %in% c("aver_log2_DEHP", "aver_log2_DiNP", "aver_log2_DINCH")], 
                                family = gaussian())
    
    res_qg_fit <- data.frame(child_sex = sex,
                             CpG = colnames(meth_clean_M)[i], 
                             est_psi1 = unname(qc.fit$psi),
                             CI_psi1_lower = unname(qc.fit$ci)[1],
                             CI_psi1_upper = unname(qc.fit$ci)[2],
                             p_val_psi1 = unname(qc.fit$pval[2]))
    
    results_gq <- bind_rows(results_gq, res_qg_fit)
  }
  
  results_gq <- results_gq |> 
    mutate(child_sex == sex)
  
  export(results_gq, here(path_main, str_c("results_mixture_", sex, ".RDS")))
  
  results_gq_sex <- bind_rows(results_gq_sex, results_gq)
}

export(results_gq_sex, here(path_main, "13.Mixture_analysis", "results_mixture_per_sex.RDS"))


## -------------------------------------------------------------------------------------------------

## Add FDR correction of p-value

results_gq_sex <- import(here(path_main, "13.Mixture_analysis", "results_mixture_per_sex.RDS"))
results_gq_sex <- results_gq_sex |> 
  
  group_by(child_sex) |> 
  mutate(FDR_p_value = p.adjust(p_val_psi1, method = "BH"))

min(results_gq_sex$FDR_p_value) # 0.06 (no significant effect of the mixture)


## -------------------------------------------------------------------------------------------------

## Add CpGs features needed to create BED files for further DMR analysis

results_gq_sex <- nl_EWAS_agnostic_quant |> 
  select(CpG, Chr, Position) |> 
  distinct(CpG, .keep_all = TRUE) |> 
  merge(results_gq_sex, by = "CpG") |> 
  mutate(Exposure = "mixture",
         raw_p_value = p_val_psi1) |> 
  select(child_sex, Exposure, CpG:Position, raw_p_value)


## -------------------------------------------------------------------------------------------------
# Create bed files for each exposure so they can be processed by the comb-p
for (sex in c("Male", "Female")) {
  
  bedCreator(annotated_EWAS_result = filter(results_gq_sex, child_sex == sex),
             prefix = paste0("agnostic_quant_", sex),
             path = paste0(path_main, "/13.Mixture_analysis/bed_files"))
}

# In the comb-p analysis no DMRs for mixtures were identified, neither for both sexes and in sex-stratified analysis.

## -------------------------------------------------------------------------------------------------

## 14. Enrichment anaysis - results provided by Lucile Broséus

# The following analyses will be performed:
#   
# -the nb of phenotypes by exposure / family of exposure 
# -look at whether the phenotypes overlap between some exposures/some flimsy of exposure
# -look at whether the phenotypes overlap between boys and girls
# 
# 
# More restrictive conditions: 
#   
# - filtering significant terms with a more stringent threshold (FDR < 0.01);
# - considering only those terms associated with at least 20 genes (because tests for smaller gene sets are likely inconsistent).
# 
# As is, the table is ordered by exposure, child_sex and decreasing combined score.
# With these selection criteria, we find 65 distinct significant terms.


# Create figures that will help to interpret the results of the enrichment analysis
dbgap_filtered_FDR001_20GenesPerTerm |> 
  ggplot(aes(x = fct_infreq(Term))) +
  facet_wrap("child_sex",
             nrow = 3) +
  geom_bar() +
  scale_x_discrete(guide = guide_axis(angle = 45))  +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 1.2, "cm"),
        axis.text = element_text(size = 10))

ggsave(filename = "enrichment_phenotypes.jpg", path = here(path_suppl, "figures"), width = 2400, height = 1500, units = "px", dpi = 200)


dbgap_filtered_FDR001_20GenesPerTerm |> 
  mutate(expo_family = ifelse(str_detect(exposure, c("BPA|BPS|BUPA|BP3|ETPA|MEPA|PRPA")), "phenol", "phthalate")) |> 
  ggplot(aes(x = fct_infreq(Term))) +
  facet_grid(rows = vars(child_sex),
             cols = vars(expo_family)) +
  geom_bar() +
  scale_x_discrete(guide = guide_axis(angle = 45))  +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 1.2, "cm"),
        axis.text = element_text(size = 10))

ggsave(filename = "enrichment_phenotypes_expo_family.jpg", path = here(path_suppl, "figures"), width = 7400, height = 1500, units = "px", dpi = 200)

dbgap_filtered_FDR001_20GenesPerTerm |> 
  mutate(expo_subfamily = case_when(str_detect(exposure, c("BPA|BPS")) ~ "bisphenol", 
                                    str_detect(exposure, c("BUPA|ETPA|MEPA|PRPA")) ~ "paraben",
                                    str_detect(exposure, c("BP3")) ~ "BP3",
                                    str_detect(exposure, c("MECPP|MEHHP|MEHP|MEOHP|MMCHP|DEHP")) ~ "DEHP",
                                    str_detect(exposure, c("ohMINCH|oxoMINCH|DINCH")) ~ "DINCH",
                                    str_detect(exposure, c("cxMiNP|ohMiNP|oxoMiNP|DiNP")) ~ "DiNP",
                                    .default = "other_phthalate")) |> 
  ggplot(aes(x = fct_infreq(Term))) +
  facet_grid(rows = vars(child_sex),
             cols = vars(expo_subfamily)) +
  geom_bar() +
  scale_x_discrete(guide = guide_axis(angle = 45))  +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 1.2, "cm"),
        axis.text = element_text(size = 10))

ggsave(filename = "enrichment_phenotypes_expo_subfamily.jpg", path = here(path_suppl, "figures"), width = 2400, height = 1500, units = "px", dpi = 200)


dbgap_filtered_FDR001_20GenesPerTerm |> 
  ggplot(aes(x = fct_infreq(Term))) +
  facet_wrap("exposure") +
  geom_bar() +
  scale_x_discrete(guide = guide_axis(angle = 45))  +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 1.2, "cm"),
        axis.text = element_text(size = 10))

ggsave(filename = "enrichment_phenotypes_exposure.jpg", path = here(path_main, "figures"), width = 2400, height = 1500, units = "px", dpi = 200)


## -------------------------------------------------------------------------------------------------
enrichR_dbgap_filter <- dbgap_filtered_FDR001_20GenesPerTerm |> 
  mutate(expo_family = ifelse(str_detect(exposure, c("BPA|BPS|BUPA|BP3|ETPA|MEPA|PRPA")), "phenol", "phthalate"),
         expo_subfamily = case_when(str_detect(exposure, c("BPA|BPS")) ~ "bisphenol", 
                                    str_detect(exposure, c("BUPA|ETPA|MEPA|PRPA")) ~ "paraben",
                                    str_detect(exposure, c("BP3")) ~ "",
                                    str_detect(exposure, c("MECPP|MEHHP|MEHP|MEOHP|MMCHP|DEHP")) ~ "DEHP",
                                    str_detect(exposure, c("ohMINCH|oxoMINCH|DINCH")) ~ "DINCH",
                                    str_detect(exposure, c("cxMiNP|ohMiNP|oxoMiNP|DiNP")) ~ "DiNP",
                                    .default = "other_phthalate"))

tab_phen_per_expo <- enrichR_dbgap_filter |> 
  group_by(child_sex, exposure) |> 
  count(Term)

# By exposure, one phenotype per exposure (no overlap)
# =======
tab_phen_per_family_both_sexes <- enrichR_dbgap_filter |> 
  group_by(child_sex, expo_family, expo_subfamily) |> 
  count(Term) |> 
  filter(child_sex == "both_sexes")

tab_phen_per_family_females <- enrichR_dbgap_filter |> 
  group_by(child_sex, expo_family, expo_subfamily) |> 
  count(Term) |> 
  filter(child_sex == "Female")

tab_phen_per_family_males <- enrichR_dbgap_filter |> 
  group_by(child_sex, expo_family, expo_subfamily) |> 
  count(Term) |> 
  filter(child_sex == "Male")

tab_phen_per_family <- enrichR_dbgap_filter |> 
  group_by(child_sex, expo_family, expo_subfamily) |> 
  count(Term) |> 
  ungroup()

export(tab_phen_per_family, here(path_main, "14.Enrichment", "tab_phen_per_family_reduced.csv"))

# Both sexes
both <- tab_phen_per_family |> 
  rowwise() |> 
  mutate(term_expo_subfamily = paste(expo_subfamily, Term)) |> 
  filter(child_sex == "both_sexes") |> 
  select(term_expo_subfamily)

# Females
f <- tab_phen_per_family |> 
  rowwise() |> 
  mutate(term_expo_subfamily = paste(expo_subfamily, Term)) |> 
  filter(child_sex == "Female") |> 
  select(term_expo_subfamily)

# Males
m <- tab_phen_per_family |> 
  rowwise() |> 
  mutate(term_expo_subfamily = paste(expo_subfamily, Term)) |> 
  filter(child_sex == "Male") |> 
  select(term_expo_subfamily)

venn.diagram(
  x = list(pull(both, term_expo_subfamily), 
           pull(f, term_expo_subfamily), 
           pull(m, term_expo_subfamily)),
  category.names = c("both" , "female" , "Male"),
  filename = 'venn_diagramm_sex_reduced.png',
  output=TRUE
)

overlap_males_females <- intersect(f, m)
export(overlap_males_females, here(path_main, "14.Enrichment", "overlap_males_females_reduced.csv"))

bisphenols <- tab_phen_per_family |> 
  filter(expo_subfamily == "bisphenol") |> 
  select(Term)

parabens <- tab_phen_per_family |> 
  filter(expo_subfamily == "paraben") |> 
  select(Term)

bp3 <- tab_phen_per_family |> 
  filter(expo_subfamily == "BP3") |> 
  select(Term)

dehp <- tab_phen_per_family |> 
  filter(expo_subfamily == "DEHP") |> 
  select(Term)

dinp <- tab_phen_per_family |> 
  filter(expo_subfamily == "DiNP") |> 
  select(Term)

dinch <- tab_phen_per_family |> 
  filter(expo_subfamily == "DINCH") |> 
  select(Term)

other_phth <- tab_phen_per_family |> 
  filter(expo_subfamily == "other_phthalate") |> 
  select(Term)

venn.diagram(
  x = list(pull(bisphenols, Term), 
           pull(parabens, Term), 
           pull(bp3, Term)),
  
  category.names = c("bisphenols" , "parabens" , "bp3"),
  filename = 'venn_diagramm_expo_phenols_reduced.png',
  output=TRUE
)

venn.diagram(
  x = list(pull(dehp, Term),
           pull(dinp, Term),
           pull(dinch, Term),
           pull(other_phth, Term)),
  
  category.names = c("dehp", "dinp", "dinch", "other_phth"),
  filename = 'venn_diagramm_expo_phthalates_reduced.png',
  output=FALSE
)

overlap_phenols <- intersect(intersect(bisphenols, parabens), bp3)
export(overlap_phenols, here(path_main, "14.Enrichment", "overlap_phenols_reduced.csv"))


overlap_phthalates <- intersect(intersect(intersect(dehp, dinp), dinch), other_phth)
export(overlap_phthalates, here(path_main, "14.Enrichment", "overlap_phthalates_reduced.csv"))


## CREATE TABLES AND FIGURES FOR THE PAPER

## TABLES

## TABLE 1

# Population characteristics for the mother–child pairs included (n = 395) and
# excluded (n = 84) from the study.

# Load results in excel format for included individuals
Table_1_1 <- import(here(path_main, "population_charactersitics_incl.csv")) %>% 
  as.data.frame()

colnames(Table_1_1) <- c("Covariate", "Distribution_included", "N_included")

# remove repeated entries for missing values
Table_1_1$Covariate[10] <- "    'Missing_smoke'"
Table_1_1$Covariate[18] <- "    'Missing_edu'"
Table_1_1$Covariate[23] <- "    'Missing_bmi'"

# Load results in excel format for excluded individuals
Table_1_2 <- import(here(path_main, "population_charactersitics_excl.csv")) %>% 
  as.data.frame()

colnames(Table_1_2) <- c("Covariate", "Distribution_excluded", "N_excluded")

# remove repeated entries for missing values
Table_1_2$Covariate[10] <- "    'Missing_smoke'"

# Merge included and excluded individuals
Table_1 <- full_join(Table_1_1, Table_1_2, by = "Covariate")

export(Table_1, here(path_main, "tables", "Table_1.csv"))
Table_1


## -------------------------------------------------------------------------------------------------

## TABLE 2

# Candidate CpGs selected from the study on males from the French EDEN cohort (2003–2006) (Jedynak et al., 2021,2022) 
# and differentially methylated in the current study (nominal p-value <0.05) in association with pregnancy exposure 
# biomarker concentrations (n = 387).

Table_2 <- EWAS_EDEN_SEPAGES %>% 
  filter(raw_p_value_SEPAGES < 0.05) %>% 
  merge(mean_meth, by = "CpG") %>% 
  select(child_sex:Exposure_name, CpG, Chr:Location_in_gene, Mean_meth, SD_meth, everything(), -Estimate_SEPAGES) %>% 
  mutate_at(vars(contains("meth"), contains("p_value")), round, 3)

export(Table_2, here(path_main, "tables", "Table_2.csv"))


## -------------------------------------------------------------------------------------------------

## FIGURES

## FIGURE 1

# Workflow of the study.
# 
# Figure made externally.


## FIGURE 2

# Maternal urinary concentrations of synthetic phenols, phthalates, and DINCH metabolites 
# (concentrations standardized for analytical batch and sampling conditions and log2-transformed), 
# restricted to the reduced population (n = 387 pregnant women), where 2% of the extreme continuous exposure 
# biomarker values (n = 8 for each biomarker) were removed.


exposure_char_red_plot %>% 
  filter(!Exposure_long %in% c("Bisphenol S", "Butylparaben")) %>% 
  mutate(aver_expo_conc_Std_log2 = log2(aver_expo_conc_Std)) %>% 
  ggplot(aes(Exposure_long, aver_expo_conc_Std_log2, fill = expo_group)) +
  geom_boxplot(outlier.size = 1) +
  scale_fill_manual(values = c("blue", "gray", "green"), name = "Exposure family") +
  theme_bw() +
  xlab("") + 
  ylab("Standardized averaged log2-transformed concentration") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title.y = element_text(size = 9))

ggsave(filename = "Figure_2.jpg", path = here(path_main, "figures"), width = 1400, height = 900, units = "px", dpi = 200)

## -------------------------------------------------------------------------------------------------

# FIGURE 3

# Heat map summarizing the enrichment analysis on differentially methylated CpGs associated with pregnancy exposure 
# in the exploratory study (p-value <0.001), when both sexes were studied together and separately (n = 387, 752,577 CpGs).

# Figure provided by Lucile Broséus.


## -------------------------------------------------------------------------------------------------

## FIGURE 4

# DMRs (≥5 probes) associated with pregnancy exposure (Šidák-corrected p-value <0.05, n = 387, 752,577 CpGs) in females (yellow), males (orange), 
# and when both sexes were studied together (red). 

# Figure made externally.

## -------------------------------------------------------------------------------------------------

## FIGURE 5

# Volcano plot for CpGs within DMRs (≥5 probes) associated with pregnancy exposure in the exploratory study 
# (Šidák-corrected p-value <0.05, n = 387, 752,577 CpGs), when both sexes were studied together and separately.

# Clean the names
CpG_in_DMRs5_plot <- CpGs_in_DMRs_5probes %>%
  select(child_sex, CpG, Expo_subfamily, Exposure_name, Gene, Estimate, Estimate_CI, raw_p_value) %>% 
  mutate(child_sex = factor(child_sex),
         child_sex = relevel(child_sex, ref = "Female and male"),
         Expo_subfamily = factor(Expo_subfamily, levels = c("Bisphenols", "Parabens", "Benzophenone-3", "DEHP metabolites", "DINCH metabolites", "DiNP metabolites", "Other phthalates")))

jpeg(here(path_main, "figures", "Figure_5.jpeg"),
     width = 2200, height = 1200, res = 250)

ggplot(CpG_in_DMRs5_plot, aes(x = Estimate,
                              y = -log10(raw_p_value),
                              col = Expo_subfamily)) +
  facet_wrap(vars(child_sex), ncol = 3) +
  geom_point() +
  labs(colour = "") +
  scale_color_manual(values = c("red1",
                                "darkorange1",
                                "deeppink",
                                "#02e302",
                                "forestgreen",
                                "blue3",
                                "blueviolet")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  ylab(bquote(~-log[10]~"(nominal p-value)")) +
  xlab(expression(beta))

dev.off()



## -------------------------------------------------------------------------------------------------

## SUPPLEMENTARY TABLES

## SUPPLEMENTARY TABLE 1

# List of the metabolites investigated with the name of their parent compounds, if relevant.

# Table made externally.

## -------------------------------------------------------------------------------------------------

## SUPPLEMENTARY TABLE 2

# Number of differentially methylated placental CpGs associated with maternal exposure in the exploratory study 
# on the whole (n = 395) and reduced (n = 387) populations, when both sexes were analyzed together and in the 
# sex-stratified analysis (752,577 CpGs).

# Merge both sexes and sex-stratified analyses for the whole population
nl_EWAS_res <- rbind(nl_EWAS_agnostic_res, nl_EWAS_agnostic_sex_res)

nl_EWAS_res <- nl_EWAS_res %>% 
  mutate(Exposure_name = factor(Exposure_name, levels = levels(nl_EWAS_res$Exposure_name), labels = c(expo_labels, expo_labels_cat)))

# Select CpGs
nl_EWAS_res_CpGs <- nl_EWAS_res %>% 
  select(Exposure, child_sex, CpG)


# Merge both sexes and sex-stratified analyses for the reduced population
nl_EWAS_quant_res <- rbind(nl_EWAS_agnostic_quant_res, nl_EWAS_agnostic_quant_sex_res) 

nl_EWAS_quant_res <- nl_EWAS_quant_res %>% 
  mutate(Exposure_name = factor(Exposure_name, levels = levels(nl_EWAS_quant_res$Exposure_name), labels = c(expo_labels, expo_labels_cat)))

# From reduced population, select only those CpGs that were identified in the whole population
nl_EWAS_quant_res <- nl_EWAS_quant_res %>% 
  inner_join(nl_EWAS_res_CpGs)

# Display significant hits for the whole population
nl_EWAS_res_tab <- nl_EWAS_res %>% 
  group_by(child_sex, Exposure_name) %>% 
  summarise(n_CpGs = n()) %>% 
  arrange(child_sex, desc(n_CpGs))

# Display significant hits for the reduced population
nl_EWAS_quant_res_tab <- nl_EWAS_quant_res %>% 
  group_by(child_sex, Exposure_name) %>% 
  summarise(n_CpGs = n()) %>% 
  arrange(child_sex, desc(n_CpGs))

Suppl_Table_2 <- nl_EWAS_res_tab %>% 
  full_join(nl_EWAS_quant_res_tab, by = c("child_sex", "Exposure_name")) %>%
  mutate(Expo_family = ifelse(str_detect(Exposure_name, "BPA|BPS|BP-3|BUPA|ETPA|MEPA|PRPA"), "Phenol", "Phthalate")) %>% 
  arrange(child_sex, Expo_family, desc(n_CpGs.x)) %>% 
  select(child_sex, Expo_family, Exposure_name, everything())

export(Suppl_Table_2, here(path_main, "tables", "Suppl_Table_2.csv"))


## -------------------------------------------------------------------------------------------------

## SUPPLEMENTARY TABLE 3

# Average maternal urinary exposure concentrations assessed in weekly pools collected twice during pregnancy 
# in the whole (n = 395 pregnant women) and reduced (n = 387) population where 2% of the extreme continuous 
# exposure values were removed (n = 8 for each compound except for categorized bisphenol S and butylparaben).

Suppl_Table_3 <- exposure_char %>% 
  merge(select(exposure_char_red, -LOD), by = "Exposure") %>% 
  arrange(Expo_family, Exposure) %>% 
  select(Exposure_long, Expo_family, everything())

expo_char_all <- expo_char_all[c(2, 3, 1, 6, 5, 7, 4, 16, 19, 18, 11, 22, 10, 23, 25, 8, 12, 13, 14, 15, 17, 9, 21, 24), ]

export(Suppl_Table_3, here(path_main, "tables", "Suppl_Table_3.csv"))


## -------------------------------------------------------------------------------------------------

## SUPPLEMENTARY TABLE 4

# Differentially methylated CpGs associated with pregnancy exposure in the exploratory study, when both sexes were 
# studied together and separately (n = 387, 752,577 CpGs).

Suppl_Table_4 <- nl_EWAS_quant_res %>% 
  merge(mean_meth, by = "CpG") %>% 
  
  select(child_sex, Exposure_name, CpG:Location_in_gene, Mean_meth, SD_meth, Estimate_CI, raw_p_value, p_value_FDR, -Exposure) %>% 
  arrange(child_sex, Exposure_name, CpG) %>% 
  mutate_at(vars(Mean_meth, SD_meth, p_value_FDR), round, 3)

export(Suppl_Table_4, here(path_main, "tables", "Suppl_Table_4.csv"))


## -------------------------------------------------------------------------------------------------

## SUPPLEMENTARY TABLE 5

# Enrichment analysis on differentially methylated CpGs associated (p-value < 0.001) with pregnancy exposure 
# in the exploratory study, when both sexes were studied together and separately (n = 387, 752,577 CpGs).

enrichR_dbgap_filter %>% 
  select(expo_family, expo_subfamily, exposure, child_sex, Term:Adjusted.P.value, Odds.Ratio:Genes) %>% 
  arrange(expo_family, expo_subfamily, exposure, child_sex) %>% 
  export(here(path_main, "tables", "Suppl_Table_5.csv"))


## -------------------------------------------------------------------------------------------------

## SUPPLEMENTARY TABLE 6

# Number of DMRs (with at least five CpGs) associated with pregnancy exposure in the exploratory study, 
# when both sexes were studied together and separately (Šidák-corrected p-value <0.05; n = 387, 752,577 CpGs).

Suppl_Table_6 <- DMR_agnostic_quant_res_5_tab %>% 
  pivot_wider(names_from = child_sex,
              values_from = n) %>% 
  replace(is.na(.), 0)

export(Suppl_Table_6, here(path_main, "tables", "Suppl_Table_6.csv"))


## -------------------------------------------------------------------------------------------------

## SUPPLEMENTARY TABLE 7

# DMRs associated with pregnancy concentrations of phenols and phthalates in the exploratory study, when both 
# sexes studied together and separately (Šidák-corrected p-value <0.05, n = 387, 752,577 CpGs).

export(DMR_agnostic_quant_res, here(path_main, "tables", "Suppl_Table_7.csv"))


## -------------------------------------------------------------------------------------------------
# Placenta
intersect(DMR5_agnostic_quant_res_genes, phthalates_genes$Gene_Jedynak) # FGF12 MECPP (+) Jedynak 2022 males; OH-MINCH (-) females Jedynak 2023
intersect(DMR5_agnostic_quant_res_genes, phthalates_genes$Gene_Grindler) # RNF39 23 phthalates mix (-); cx-MiNP (+) females Jedynak 2023
# RIMS1 phthalates mix Grindler (-); ETPA (+) males and females Jedynak 2023
# GNAS phthalates mix Grindler (+); BPS (+) males Jedynak 2023
intersect(DMR5_agnostic_quant_res_genes, phthalates_genes$Gene_petroff) # FOXA1 (+) Miura; cx-MiNP (+) Jedynak 2023 
intersect(DMR5_agnostic_quant_res_genes, phthalates_genes$Gene_LaRocca_Zhao) # none

# cord blood
intersect(DMR5_agnostic_quant_res_genes, phthalates_genes$Gene_solomon) # RNF39 MEHP (-)
intersect(DMR5_agnostic_quant_res_genes, phthalates_genes$Gene_Montrose) # none
intersect(DMR5_agnostic_quant_res_genes, phthalates_genes$Gene_Tindula) # studied SGCE, PEG10, PEG3 but no signal for those
intersect(DMR5_agnostic_quant_res_genes, phthalates_genes$Gene_Chen) # none
intersect(DMR5_agnostic_quant_res_genes, phthalates_genes$Gene_Miura) # none


## -------------------------------------------------------------------------------------------------

## SUPPLEMENTARY TABLE 8

# CpGs identified within the exploratory DMR study and overlapping with the eQTMs associating methylation levels 
# with the expression of the nearby genes (n = 387, 752,577 CpGs).

Suppl_Table_8 <- filter(CpGs_in_DMRs_5probes, CpG %in% intersecting_cpgs_delhaye) |> 
  arrange(Gene, Exposure) %>% 
  export(here(path_main, "tables", "Suppl_Table_8.csv"))


## -------------------------------------------------------------------------------------------------

## SUPPLEMENTARY FIGURES

## SUPPLEMENTARY FIGURE 1

# Figure made externally.
# Study flow chart.

## -------------------------------------------------------------------------------------------------

## Save session info
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

