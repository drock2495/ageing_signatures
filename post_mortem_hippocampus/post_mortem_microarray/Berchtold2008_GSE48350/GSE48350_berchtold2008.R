library(Biobase)
library(GEOquery)
library(tidyverse)
library(limma)
library(genefilter)
library(SummarizedExperiment)

----------## READ IN TARGETS FILE ###-----------------------------------------------


setwd("~/OneDrive - King's College London/Data/Bioinformatic/External/Ageing Datasets/Hippocampus")


# Using the GEO2Query, a targets.txt file was created with first column as filenames and sample info then saved in the 
# directory of each GSE containing the raw CEL files

targets_48350 <- readTargets(file = "targets.txt", path = "data/GSEseries_raw/GSEseries_raw/GSE48350_RAW")
targets_48350 <- targets_48350[1:35, -c(7:11)] # cleaning


---------------### READ IN RAW DATA, RMA CORRECTION AND NORMALISATION ###----------------------

# package is required for annotation of 25219 array

library(oligo)


setwd("~/OneDrive - King's College London/Data/Bioinformatic/External/Ageing Datasets/Hippocampus/data/GSEseries_raw/GSEseries_raw/GSE48350_RAW")

# output of code below is an ExpressionFeatureSet
rawData_48350 <- read.celfiles(filenames = targets_48350$Filename)
# Now we run RMA background correction, normalization and probe summarisation in one step
eset48350 <- rma(rawData_48350)
# rename samples to GSM idenitifiers
colnames(eset48350) <- targets_48350$Accession


-------- ### CREATE PHENOTYPE, FEATURE DATA AND EXPRESSION MATRIX ###------------------

## FEATURE DATA ##

# map Affy IDs to entrez and symbols

library(biomaRt)
ensembl <-useMart("ensembl",dataset="hsapiens_gene_ensembl")

feature_data <- getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol", "entrezgene_id"), 
                      filters = "affy_hg_u133_plus_2", values = row.names(exprs(eset48350)), mart = ensembl)
biomart_ids <- featuredata

probesdf <- data.frame(affy_hg_u133_plus_2 = row.names(exprs(eset48350)))

# 36,367 unique probes that have either symbol, entrez or both
fData_mapped_probes <- left_join(probesdf, biomart_ids) %>%
  as_tibble() %>%  # 60,093 probes
  rename(affy_hg_u133_plus_2 = "probe_id") %>% 
  filter(isUnique(probe_id)) %>% # 49,211 unique probes, filtering probes that map to >1 gene (one-to-many probes)
  mutate_at(.vars = c("hgnc_symbol", "entrezgene_id"), .funs = list(~ifelse(.=="", NA, as.character(.)))) %>% 
  as_tibble() %>% 
  filter(!is.na(hgnc_symbol) | !is.na(entrezgene_id)) # filter out probes which did not map to neither symbol or entrez (totally unmapped)

# For probes mapping to same gene (many-to-one), take the median expression value
probeSum_exprs <- exprs(eset48350) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "probe_id") %>% 
  as_tibble() %>% # 54,675 unique probes measured on array
  right_join(., fData_mapped_probes) %>%
  group_by(hgnc_symbol,entrezgene_id) %>% 
  summarise_at(.vars = c(2:36), median) # 18,431 unique genes


# We have more probes mapping to symbols than entrez ids, so we will use symbols to annotate our expression matrix
# There are 19 genes that have an entrez but no symbol. We'll assign pseudo_IDs to these using row numbers so 
# that we can keep them for analysis and work out what they are later on (through feature data)

idx <- which(is.na(probeSum_exprs$hgnc_symbol))
probeSum_exprs[idx, ]$hgnc_symbol <- paste0("uid_", 1:nrow(probeSum_exprs[idx, ]))

feature_data <- as.data.frame(probeSum_exprs[, 1:2])
row.names(feature_data) <- feature_data$hgnc_symbol


## PHENOTYPE DATA ##


# clean up targets df to become phenotype data

targets_48350$Age <- str_extract(targets_48350$Age, "[0-9][0-9]+")
# save sample names for rownames later
samples <- targets_48350$Accession

targets_48350 <- targets_48350[, -c(1, 2)]
targets_48350$Age <- as.numeric(targets_48350$Age)

# create group column idenifying young or old
phenotype_data <- targets_48350 %>% 
  mutate(Group = ifelse(Age <= 40, "young", ifelse(Age >= 60, "old", "")))

row.names(phenotype_data) <- samples


## EXPRESSION MATRIX ##

gene_symbols <- probeSum_exprs$hgnc_symbol
probeSum_exprs <- as.data.frame(probeSum_exprs[, str_detect(colnames(probeSum_exprs), "GSM")])
row.names(probeSum_exprs) <- gene_symbols
expression_matrix <- as.matrix(probeSum_exprs)

# Save expression matrix, fdata and pdata 

# Save expression matrix (must move rownames to new column before save)
expression_data <- as.data.frame(cbind(rownames(expression_matrix), expression_matrix))
colnames(expression_data)[1] <- "gene_id"
write_tsv(expression_data, "data/expression_data.txt")

#save feature data
write_tsv(feature_data, "data/feature_data.txt")

#save phenotype data (again move sample rownames to new column)
pdata <- phenotype_data %>% 
  rownames_to_column(var = "sample_id")
write_tsv(pdata, path = "data/phenotype_data.txt")


--------### CREATE NEW EXPRESSION SET OBJECT ### ----------------------------

setwd("~/OneDrive - King's College London/Data/Bioinformatic/External/Ageing Datasets/Hippocampus/data/Berchtold_GSE48350")

first_col_to_rownames <- function(df) {
  names <- df[[1]]
  df <- as.data.frame(df)
  rownames(df) <- names
  df2 <- df[, -1]
  df2
}

x <- as.matrix(first_col_to_rownames(read_tsv("expression_data.txt")))

f <- as.data.frame(read_tsv("feature_data.txt"))
row.names(f) <- f$hgnc_symbol


p <- first_col_to_rownames(read_tsv("phenotype_data.txt"))


eset <- ExpressionSet(assayData = x, phenoData = AnnotatedDataFrame(p), 
                      featureData = AnnotatedDataFrame(f))


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
PCA1v2 <- plotPCA.any(eset = eset, pc.x = 1, pc.y = 2, intgroup = "Group", labels = TRUE)

ggsave("PCAplot_berchtold2008.tiff", path = "figures", plot = PCA1v2, device = "tiff",
       width = 8, height = 8)

--------### DESIGN MATRIX ### ----------------------------

## Create variables to input into the design matrix ##

Group <- as.factor(pData(eset)[, "Group"])
Sex <- as.factor(pData(eset)[, "Sex"])

# Build design matrix, looking at effect of Young vs Old whilst controlling for Sex

design <- model.matrix(~0 + Group + Sex)

colnames(design) <- c("Old", "Young", "Male")

--------### FIT THE MODEL ### ----------------------------

----# FILTER LOW EXPRESSED GENES

# As we don't have background control probe files...
# We can filter out low expressed probes based on their average expression computed by lmFit

fit <- lmFit(eset, design)

# Two diagnostic plots to help us decide intensity cutoff

# Histogram of average expression values
hist(fit$Amean)

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
# Trend = TRUE takes account of the mean -variance relationship when calculating log2 fold changes
fit2 <- eBayes(fit2, trend=TRUE)
plotSA(fit2) # blue line is the mean-variance fit



## find top DE genes between contrast coefficients ##

results <- decideTests(fit2)
summary(results)

## obtain DEGs for full fitted model ranked by p value ##

all_genes <- topTable(fit2, number = nrow(fit2))

write_tsv(all_genes, "results/diffExpr_all_genes_berchtold.txt")

