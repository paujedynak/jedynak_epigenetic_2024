rm(list = ls())
################################################################################
# SEPAGES - Preprocessing methylation data 
################################################################################
# Author: Lucile
# Date: Finalized in April 2021 (EDEN)
# Update: 2nd November 2021 (SEPAGES)
# Update: 8th November 2021 (check for chip)
# Last update: 8th December 2021 (separate datasets for chrX and chr Y + plate recoding)
################################################################################
# Design:
# - merge several raw files to generate an complete design file
# Pre-processing of DNAme data:
# 1. BMIQ normalisation (already done by the platform)
# 1'. Keep only one technical replicate per individual (randomly picked)
# 2. Removal of CpGs close to known SNP
# 3. Removal of cross-reactive CpGs 
# 4. Three datasets: autosomes/chrX/chrY
#------------------------------------------------------------------------------#
# Notes:
#
################################################################################
# Paths and parameters
#------------------------------------------------------------------------------#

# BMIQ data provided by Johanna, contains duplicates
# (from CNRGH platform)
normFile <- "~/Data/SEPAGES/RAW_DATA/myNorm.RData"

array_type <- "EPIC"

# Design
# (from CNRGH platform)
designFile <- "~/Data/SEPAGES/IDAT/SampleSheet_CANCAIR_P31-35.csv"

# Hand-made design (chip+plate)
# Using file test_visualisation_SEPAGES_final_versionExcel_20211014.xlsx
# (from Johanna and Emie)
hmdesignFile <- "~/Data/SEPAGES/SEPAGES_design_2021-10-15.xlsx"
  
# correspondance sample <-> id 
# (provided by Sarah - 22/10/2021)
correspFile <- "~/Data/SEPAGES/RAW_DATA/placenta_SEC_crb_jointure_20211021.dta"

# Seed for reproducing random picking (replicates)
seed <- 38

#------------------------------------------------------------------------------#
# Outputs
#------------------------------------------------------------------------------#

# BMIQ without duplicates 
bvalFile <- "~/Work/SEPAGES/Data/bmiq_no_duplicates.rds"

# New design, all samples
newDesignFile.all <- paste0("~/Work/SEPAGES/Data/SEPAGES_SampleSheet_", format(Sys.time(), "%Y%m%d"),".csv")
# New design, without replicates
newDesignFile.nodup <- paste0("~/Work/SEPAGES/Data/SEPAGES_SampleSheet_nodup_", format(Sys.time(), "%Y%m%d"),".csv")

#Output file: meth. data after removing SNPs, cross-reactive and CHR XY probes
bvalFile.final <- "~/Work/SEPAGES/Data/bmiq_processed.rds"
bvalFile.chrX <- "~/Work/SEPAGES/Data/bmiq_processed_chrX.rds"
bvalFile.chrY <- "~/Work/SEPAGES/Data/bmiq_processed_chrY.rds"

################################################################################
# R packages
#------------------------------------------------------------------------------#

library(magrittr)

################################################################################
# 1. BMIQ normalisation 
#------------------------------------------------------------------------------#
# (cf: platform)

################################################################################
# 2. Load norm. data + remove duplicates
#------------------------------------------------------------------------------#

load(normFile) # myNorm
Bvalues <- myNorm; rm(myNorm)
dim(Bvalues)
#814481    432

# Load design
design <- read.csv(file = designFile)
dim(design)

stopifnot(identical(colnames(Bvalues), design$Sample_Name))

idsInOrder <- design$Sample_Name

# Add batch (cf: mail forwarded by Johanna)
table(design$Sample_Plate)
design$batch <- tapply(design$Sample_Plate, 
                       IND = seq_along(design$Sample_Plate),
                       FUN = function(x) stringr::str_split(x, pattern = "_", simplify = T)[3])
design <- design %>% dplyr::mutate(batch = ifelse(batch %in% c("P31", "P32"), 1,
                                                  ifelse(batch == "P35", 3, 2)))
table(design$batch)
# 1   2   3 
#192 192  48 

# Define plate
design$plate <- tapply(design$Sample_Plate, 
                       IND = seq_along(design$Sample_Plate),
                       FUN = function(x) 
                         stringr::str_split(x, pattern = "_", simplify = T)[3] %>%
                         stringr::str_sub(start = 3))
design$plate <- as.numeric(design$plate)
table(design$plate)
# 1  2  3  4  5 
# 96 96 96 96 48 


# Stats of technical replicates
replicates <- table(colnames(Bvalues))
replicates <- replicates[replicates>1]
replicates
# 14 indiv.

# Processing of replicates: pick only on sample 
# (needed for adjusting for technical effects)
indices <- tapply(X = names(replicates), 
                  INDEX = seq_along(replicates),
                  FUN = function(x) grep(pattern = x, x = colnames(Bvalues)),
                  simplify = T)

# Pick only one replicate, randomly
indices.kept <- tapply(X = names(replicates), 
                  INDEX = seq_along(replicates),
                  FUN = function(x){
                    set.seed(seed = seed)
                    ind <- grep(pattern = x, x = colnames(Bvalues))
                    return(ind[sample(x = seq_along(ind), size = 1)])
                    },
                  simplify = T)

indices.kept
 
# Methylation data matrix without replicates
indices <- unlist(indices)
Bvalues <- cbind(Bvalues[,-indices], Bvalues[,as.vector(indices.kept)])
Bvalues <- Bvalues[,order(colnames(Bvalues))]
dim(Bvalues)
# 814481    395

# Extended design file for all samples
# Load hand-made design
# design and hmdesign should match with the platform design file on (Barcode_Stock x Batch) - (Barcode_Stock x plate)
hmdesign <- readxl::read_xlsx(hmdesignFile) %>% data.frame()

stopifnot(identical(hmdesign$Barcode_Stock, design$Barcode_Stock))

hmdesign <- cbind.data.frame(hmdesign, design[,c("Sample_Name", "Sex", "batch", "plate")])
table(hmdesign$batch, hmdesign$plate)

hmdesign <- hmdesign %>% dplyr::select(Sample_Name, Sex, batch, plate, chip = Chip)

# odd plates (left side) are recoded as 1
# even plates (right side) are recoded as 2
hmdesign$plate <- ifelse(hmdesign$plate%%2==0, 2, 1)
table(hmdesign$batch, hmdesign$plate)

# Data about technical effect for kept sample (only one technical replicate per indiv.)
subdesign <- rbind.data.frame(hmdesign[-indices,], hmdesign[indices.kept,])
table(subdesign$batch, subdesign$plate)

################################################################################
# Assign indiv names instead of sample names
#------------------------------------------------------------------------------#

corresp <- haven::read_dta(correspFile) %>% data.frame()
colnames(corresp) <- c("Sample_Name", "ident")

stopifnot(nrow(corresp) == ncol(Bvalues))

corresp$Sample_Name <- tapply(corresp$Sample_Name, 
                             seq_along(corresp$Sample_Name),
                             FUN = function(x) stringr::str_split(x, pattern = "\\.", simplify = T)[1])
corresp <- corresp %>% dplyr::arrange(Sample_Name)
length(unique(corresp$Sample_Name))
# must be 395
colnames(Bvalues) <- tapply(colnames(Bvalues), seq_along(colnames(Bvalues)),
                            FUN = function(x) stringr::str_split(x, pattern = "-", simplify = T)[1])
Bvalues <- Bvalues[,order(colnames(Bvalues))]
length(unique(colnames(Bvalues)))
# must be 395

stopifnot(identical(as.character(corresp$Sample_Name),colnames(Bvalues)))
colnames(Bvalues) <- corresp$ident
Bvalues[1:5,1:5]

# Design with all samples
hmdesign$Sample_Name <- tapply(hmdesign$Sample_Name, seq_along(hmdesign$Sample_Name),
                                FUN = function(x) stringr::str_split(x, pattern = "-", simplify = T)[1])
hmdesign <- hmdesign %>% dplyr::arrange(Sample_Name)

hmdesign <- merge(hmdesign, corresp, by = "Sample_Name")
dim(hmdesign)
#432  13

# Design without technical replicates
subdesign$Sample_Name <- tapply(subdesign$Sample_Name, seq_along(subdesign$Sample_Name),
                                FUN = function(x) stringr::str_split(x, pattern = "-", simplify = T)[1])
subdesign <- subdesign %>% dplyr::arrange(Sample_Name)

stopifnot(identical(as.character(corresp$Sample_Name), as.character(subdesign$Sample_Name)))
subdesign$ident <- corresp$ident
dim(subdesign)
#395  13

################################################################################
# Save in bvalFile / new design files (methylation data without replicates)
#------------------------------------------------------------------------------#

saveRDS(Bvalues, file = bvalFile)

hmdesign <- hmdesign[, c("ident", "Sample_Name", "Sex", "batch", "plate", "chip")]
write.csv(x = hmdesign, file = newDesignFile.all, row.names = F)

subdesign <- subdesign[, c("ident", "Sample_Name", "Sex", "batch", "plate", "chip")]
write.csv(x = subdesign, file = newDesignFile.nodup, row.names = F)
 
################################################################################
# 2. Removal of CpGs close to SNPs (dist<2 bp) with MAF < 0.05
#------------------------------------------------------------------------------#

Bvalues <- readRDS(bvalFile)

Bvalues <- DMRcate::rmSNPandCH(Bvalues, 
                               dist = 2, 
                               mafcut = 0.05,
                               rmcrosshyb = FALSE, 
                               rmXY = FALSE)
dim(Bvalues)

# 802 182 CpGs 

################################################################################
# 3. Removal of cross-reactive probes
#------------------------------------------------------------------------------#

xloci <- maxprobes::xreactive_probes(array_type = array_type)
xloci <- unlist(unique(xloci))
length(xloci)
# 87 464 known cross-reactive probes

x <- which(rownames(Bvalues) %in% xloci)
length(x)
# 33 151 known cross-reactive probes still in the data set
round(length(x)/nrow(Bvalues),2)

Bvalues <- Bvalues[-x,]
dim(Bvalues)
# 769 031  CpGs 

################################################################################
# 4. Split meth. dataset into three sets: CpGs on autosomes, ChrX and ChrY
#------------------------------------------------------------------------------#

suppressPackageStartupMessages( library( IlluminaHumanMethylationEPICanno.ilm10b4.hg19 ) ) 
annotation_object <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19

annot <- minfi::getAnnotation(annotation_object)
annot <- data.frame(CpG = annot$Name, chr = annot$chr)

N <- nrow(Bvalues)

Bvalues.chrX <- Bvalues[which(rownames(Bvalues) %in% annot$CpG[annot$chr == "chrX"]),]
Bvalues.chrY <- Bvalues[which(rownames(Bvalues) %in% annot$CpG[annot$chr == "chrY"]),]
Bvalues <- Bvalues[-which(rownames(Bvalues) %in% annot$CpG[annot$chr %in% c("chrY", "chrX")]),]

stopifnot(N == (nrow(Bvalues) + nrow(Bvalues.chrX) + nrow(Bvalues.chrY)))

nrow(Bvalues.chrX)
# 16 412

# for ChrY, keep only male children ids
Bvalues.chrY <- Bvalues.chrY[, colnames(Bvalues.chrY) %in% subdesign$ident[subdesign$Sex=="M"]]
dim(Bvalues.chrY)
# 42 209

################################################################################
# Save
#------------------------------------------------------------------------------#

saveRDS(Bvalues, file = bvalFile.final)
saveRDS(Bvalues.chrX, file = bvalFile.chrX)
saveRDS(Bvalues.chrY, file = bvalFile.chrY)

################################################################################