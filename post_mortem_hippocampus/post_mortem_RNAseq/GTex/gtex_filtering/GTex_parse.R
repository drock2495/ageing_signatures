library(cmapR)
library(tidyverse)

setwd("~/OneDrive - King's College London/Data/Bioinformatic/External/Ageing Datasets/Hippocampus")

# Load GTex dataset, sample attributes and subject phenotypes
ds_path <- "data/GTex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"
GTex <- parse_gctx(ds_path)

sample_annot <- read_tsv("data/GTex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

subj_pd <- read_tsv("data/GTex/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")

# Filter to include only samples from hippocampus

brain_hip_samples <- sample_annot %>% 
  filter(SMTSD == "Brain - Hippocampus") %>% 
  separate(SAMPID, into = c("SUBJID", "code"), sep = "-0011")

# Filter subject data to include participants below 40 and above 60 years old, and cause of death 1 or 2 on
# hardy scale

selected_subj <- subj_pd %>% 
  filter(DTHHRDY %in% c(1, 2)) %>% 
  filter(AGE %in% c("20-29", "30-39", "60-69", "70-79"))

# 69 RNA-seq samples from hippocampus in correct age range and cause of death
# Note that we needed the stem ID (GTEX-XXXXX) to merge with tbls together
# but to later filter the columns of the expression matrix we need to reassemble the full
# samples identifiers


phenotype_data <- left_join(selected_subj, brain_hip_samples, by = "SUBJID") %>% 
  filter(!is.na(code)) %>% 
  filter(SMGEBTCHT == "TruSeq.v1") %>%
  separate(code, into = c("code1", "code2"), sep = "R1[:alpha:]-") %>% 
  unite(., c("SUBJID", "code2"), col = "SUBJID", sep = "-0011-") %>% 
  view()


## Accessing GCT object components

# The components of a `GCT` object can be accessed or modified using a set of accessor functions.
# access the data matrix and subset to include only our samples

# raw counts
raw_counts <- mat(GTex)[, phenotype_data$SUBJID]

# metadata
phenotype_data <- as.data.frame(phenotype_data)
row.names(phenotype_data) <- phenotype_data$SUBJID 
metadata <- phenotype_data[, -1]  

# save raw data

raw_counts <- cbind(row.names(raw_counts), raw_counts)
colnames(raw_counts)[1] <- "ensembl"

write_tsv(as_tibble(raw_counts), "data/GTex/raw_counts.txt")

metadata <- metadata %>% 
  rownames_to_column(var = "SUBJID")
write_tsv(metadata, "data/GTex/metadata.txt")
