library(Biobase)
library(GEOquery)
library(tidyverse)
library(limma)

----------## READ IN TARGETS FILE ###-----------------------------------------------


setwd("~/OneDrive - King's College London/Data/Bioinformatic/External/Ageing Datasets/Hippocampus")


# Using the GEO2Query, a targets.txt file was created with first column as filenames and sample info then saved in the 
# directory of each GSE containing the raw CEL files


targets_25219 <- readTargets(file = "targets.txt", path = "data/GSEseries_raw/GSEseries_raw/GSE25219_RAW")
targets_25219 <-  targets_25219[1:21, -14] # cleaning



---------------### READ IN RAW DATA, RMA CORRECTION AND NORMALISATION ###----------------------

# package is required for annotation of 25219 array
library(pd.huex.1.0.st.v2) 
library(oligo)

setwd("~/OneDrive - King's College London/Data/Bioinformatic/External/Ageing Datasets/Hippocampus/data/GSEseries_raw/GSEseries_raw/GSE25219_RAW")

# output of code below is an ExonFeatureSet
rawData_25219 <- read.celfiles(filenames = targets_25219$Filename, pkgname = "pd.huex.1.0.st.v2")
# Now we run RMA background correction, normalization and probe summarisation in one step
# package for this array automatically downloads
eset25219 <- rma(rawData_25219) 
# rename samples to GSM idenitifiers
colnames(eset25219) <- targets_25219$Accession

-------- ### CREATE PHENOTYPE, FEATURE DATA AND EXPRESSION MATRIX ###------------------

## FEATURE DATA ##

# Gene anotation is not available for this platform so the GPL file was downloaded from GEO as a .txt file
# on the GPL platform page. Note: first few rows explaining the column headers were deleted in excel before
# parsing into R
setwd("~/OneDrive - King's College London/Data/Bioinformatic/External/Ageing Datasets/Hippocampus")

# read in platform file for HuEx-1_0-st (affymetrix) was downloaded from GEO
GPL <- read_tsv("data/Kang2011_GSE25219/GPL5175-3188.txt")

# prepare expressed probes data frame for merge with GPL annotation
mat <- exprs(eset25219) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "ID")
mat$ID <- as.numeric(mat$ID)

# this gives us a data frame containing annotation for probes that were expressed in the experiment
annot_probes <- as_tibble(left_join(mat, GPL, by = "ID"))[, c(1, 23:33)]

# we still have the control probes lurking in there even after RMA. Remove these
controlprobes <- which(str_detect(annot_probes$SPOT_ID, "control"))
annot_probes <- annot_probes[-controlprobes, ]

# split gene assignment column on regex
annot_split <- str_split(annot_probes$gene_assignment, pattern = " // | /// ")
# give list element probe names to keep track
names(annot_split) <- annot_probes$ID

# grab gene info for list sub-elements (eg. entrez is always 5th subelement of the list)
entrez <- map_dbl(annot_split, ~ as.numeric(.x[5]))
symbol <- map_chr(annot_split, ~ as.character(.x[2]))
chr_pos <- map_chr(annot_split, ~ as.character(.x[4]))
full_name <- map_chr(annot_split, ~ as.character(.x[3]))

# Ensembl transcript ID is not same subelement index for each probe, so we subset our list using regex
ENST <- map(annot_split, ~ str_subset(.x, "ENST"))
# Take 1st sub-element of each element (probe). This ENSTXXX transcript corresponds to gene info above
ENST <- map_chr(ENST, ~ as.character(.x[1]))


# Create data frame with everything to be feature data

feature_data <- data.frame(entrez, symbol, chr_pos, full_name, ENST) %>% 
  rownames_to_column(var = "probe_id")
  
# Now we must deal with one-to-many, and many-to-one gene mappings, taking the median expression value for 
# multiple probes mapping to one gene and discarding probes mapping to multiple genes

exprs <- exprs(eset25219) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "probe_id")
feature_data$probe_id <- as.numeric(feature_data$probe_id)
feature_data$symbol <- as.character(feature_data$symbo)
exprs$probe_id <- as.numeric(exprs$probe_id)
feature_data$full_name <- as.character(feature_data$full_name)
feature_data$chr_pos <- as.character(feature_data$chr_pos)
feature_data$ENST <- as.character(feature_data$ENST)



# merge on probe id and summarise many-to-one gene mappings by taking the median expression value
# Note that some probes in GPL file are not mapped to ANY idenitifier! (bizzare). These are removed in the last line

summed_probes <- left_join(feature_data, exprs, by = "probe_id") %>% 
  as_tibble() %>% 
  dplyr::filter(!is.na(entrez)) %>% 
  group_by(entrez) %>% 
  dplyr::summarise_at(vars(starts_with("GSM")), median)

# Group_by entrez but keep rest of annotation data to merge
annotation <- left_join(feature_data, exprs, by = "probe_id") %>% 
  as_tibble() %>% 
  filter(!is.na(entrez)) %>% 
  dplyr::group_by(entrez) %>% 
  summarise(symbol = min(symbol), chr_pos = min(chr_pos), full_name = min(full_name),
            probe_id = min(probe_id), ENST = min(ENST))

# join annotation with probe summary, then split into expression matrix and final feature data
final_annot <- left_join(summed_probes, annotation, by = "entrez")

# Final feature data
feature_data <- as.data.frame(final_annot[, c(1, 23:27)])
row.names(feature_data) <- feature_data$entrez

## EXPRESSION MATRIX ##

assay_data <- as.data.frame(final_annot[, c(1:22)])
row.names(assay_data) <- assay_data$entrez
assay_data <- as.matrix(assay_data[, -1])


## PHENOTYPE DATA ##


# save sample names for use as rownames
samples <- targets_25219$Accession

pheno_data <- targets_25219[, -c(1, 2)]
row.names(pheno_data) <- samples

--------### CREATE NEW EXPRESSION SET OBJECT ### ----------------------------

setwd("~/OneDrive - King's College London/Data/Bioinformatic/External/Ageing Datasets/Hippocampus/data/Kang2011_GSE25219")

first_col_to_rownames <- function(df) {
  names <- df[[1]]
  df <- as.data.frame(df)
  rownames(df) <- names
  df2 <- df[, -1]
  df2
}

x <- as.matrix(first_col_to_rownames(read_tsv("expression matrix.txt")))

f <- first_col_to_rownames(read_tsv("feature_data.txt"))

p <- first_col_to_rownames(read_tsv("phenotype_data.txt"))

eset <- ExpressionSet(assayData = x, phenoData = AnnotatedDataFrame(p), 
                      featureData = AnnotatedDataFrame(f))

# Save expression matrix, fdata and pdata 

assay_data %>% 
  as.data.frame() %>%
  rownames_to_column(var = "entrez") %>%
  write_tsv(., "data/Kang2011_GSE25219/expression matrix.txt")

write_tsv(feature_data, "data/Kang2011_GSE25219/feature_data.txt")

pheno_data %>% 
  rownames_to_column(var = "sample_id") %>% 
  write_tsv(., "data/Kang2011_GSE25219/phenotype_data.txt")


# Check intensity curves for each sample. They should be homogenous, indicating they have been normalised correctly
plotDensities(exprs(eset), legend = F); abline(v=5)


# Do we detecting any clustering by PCA?

plotPCA.any <- function (eset, pc.x, pc.y, labels = FALSE, intgroup, ntop = 500, returnData = FALSE) 
{
  require(genefilter)
  rv <- rowVars(exprs(eset))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(exprs(eset)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(pData(eset)))) {
    stop("the argument 'intgroup' should specify columns of phenotype data in expression set object")
  }
  intgroup.df <- as.data.frame(pData(eset)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    pData(eset)[[intgroup]]
  }
  
  ## Select the PCAs and percentVar that you like instead of 1 and 2
  d <- data.frame(PCX = pca$x[, pc.x], PCY = pca$x[, pc.y], group = group, 
                  intgroup.df, sample = rownames(pData(eset)))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pc.x:pc.y]
    return(d)
  }
  
  ## CREATE PLOTS ###
  require(ggplot2)
  require(ggrepel)
  
  if (labels) {
    g <- ggplot(data = d, aes_string(x = "PCX", y = "PCY", color = "group")) + 
      geom_point(size = 3) + 
      xlab(paste0("PC", pc.x, ": ", round(percentVar[pc.x] * 100), "% Variance")) + 
      ylab(paste0("PC", pc.y,  ": ", round(percentVar[pc.y] * 100), "% Variance")) + 
      coord_fixed() + 
      geom_text_repel(size=3, aes_string(label = "sample"), color = "black") +
      scale_color_discrete(name = intgroup) +
      theme_bw() +
      theme(legend.text = element_text(size = 16), 
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14),
            legend.key.size = unit(1.5, "cm")) 
    return (g)
  } else {
    p <- ggplot(data = d, aes_string(x = "PCX", y = "PCY", color = "group")) + 
      geom_point(size = 3) + 
      xlab(paste0("PC", pc.x, ": ", round(percentVar[pc.x] * 100), "% Variance")) + 
      ylab(paste0("PC", pc.y,  ": ", round(percentVar[pc.y] * 100), "% Variance")) + 
      coord_fixed() + 
      scale_color_discrete(name = intgroup) +
      theme_bw() +
      theme(legend.text = element_text(size = 16), 
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14),
            legend.key.size = unit(1.5, "cm"))
    return(p)
  }
}
# Answer is no, not really
PCA1v2_with_outlier <- plotPCA.any(eset = eset, pc.x = 1, pc.y = 2, intgroup = "Group", labels = TRUE)

# Sample GSM704691 appears to be an outlier as it clusters entirely separately from the rest


ggsave("PCAplot_kang2011.tiff", path = "figures", plot = PCA1v2_with_outlier, device = "tiff",
       width = 8, height = 8)

# Remove outlier before running limma pipeline

assay_data <- assay_data[, -10]
pheno_data <- pheno_data[-10, ]

eset <- ExpressionSet(assayData = assay_data, phenoData = AnnotatedDataFrame(pheno_data), 
                      featureData = AnnotatedDataFrame(feature_data))

# We see that the samples cluster strongly by subject (brain code) and less strongly by sex
PCA1v2_col_by_sample <- plotPCA.any(eset = eset, pc.x = 1, pc.y = 2, intgroup = "Brain.code", labels = T)
PCA1v2_col_by_sex <- plotPCA.any(eset = eset, pc.x = 1, pc.y = 2, intgroup = "Sex", labels = T)



ggsave("PCA1v2_col_by_sex.tiff", path = "figures", plot = PCA1v2_col_by_sex, device = "tiff",
       width = 8, height = 8)

--------### DESIGN MATRIX ### ----------------------------

## Create variables to input into the design matrix ##

# Note that sex has high degree of colinearity with subject as so some coefficients can't be estimated
# Also I tried to account for subject in the design formula but again one subject coefficient was not esimatable
# Hence, we shall treat subject as a random effect and calculate this intra-subject correlation first
# This intra_subject correlation will then get account for in the linear model


Group <- factor(pData(eset)[, "Group"], levels = c("old", "young"))

design <- model.matrix(~0 + Group)
colnames(design) <- c("Old", "Young")


intra_subject <- pData(eset)[, "Brain.code"]
intra_subject <- factor(intra_subject)

corfit <- duplicateCorrelation(eset, design, block = intra_subject)
corfit$consensus


--------### FIT THE MODEL ### ----------------------------

## Then this inter-subject correlation is input into the linear model fit ##

fit <- lmFit(eset,design,block=intra_subject,correlation=corfit$consensus)

----# FILTER LOW EXPRESSED GENES
  
# As we don't have background control probe files...
# We can filter out low expressed probes based on their average expression computed by lmFit
# Two diagnostic plots to help us decide intensity cutoff

# Histogram of average expression values
hist(fit$Amean, main = "Averaged Gene Expression Distribution", 
     xlab = "Expression")

# variance (std dev residual) vs average expression
# generally higher expr values have more variance
# there are genes with some high variation even at intensity =4
# therefore we will filter genes below 4
plotSA(fit)

# Define filter and subset fitted object to remove genes
keep <- fit$Amean > 4
fit <- fit[keep,]


-----# CONTRASTS AND EBAYES 
  
# Set contrasts we wish to make 
  
cm <- makeContrasts(OldvYoung = Old - Young,
                      levels = design)

# fit the contrasts
fit2 <- contrasts.fit(fit, cm)

# Run empirical bayes function for differential expression testing

fit2 <- eBayes(fit2)

## find top DE genes between contrast coefficients ##

results <- decideTests(fit2)
summary(results)

## obtain DEGs for full fitted model ranked by p value ##

all_genes <- topTable(fit2, number = nrow(fit2))

write_tsv(all_genes, "results/diffExpr_all_genes_kang2011.txt")


