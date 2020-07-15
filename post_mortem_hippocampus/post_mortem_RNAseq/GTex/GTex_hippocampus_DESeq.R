

###----------------CREATE DESEQ DATASET OBJECT--------------------------------------

setwd("~/OneDrive - King's College London/Data/Bioinformatic/External/Ageing Datasets/Hippocampus")

library(DESeq2)
library(ggplot2)
library(ggfortify)
library(tidyverse)


first_col_to_rownames <- function(df) {
  names <- df[[1]]
  df <- as.data.frame(df)
  rownames(df) <- names
  df2 <- df[, -1]
  df2
}

metadata <- first_col_to_rownames(read_tsv("data/GTex/metadata.txt"))

metadata$Sex <- as.factor(metadata$Sex)

raw_counts <- first_col_to_rownames(read_tsv("data/GTex/raw_counts.txt"))

# create the deseq object and specify experimental design

dds <- DESeqDataSetFromMatrix(countData =  raw_counts,
                              colData =  metadata,
                              design = ~ Group)

###-----UNSUPERVISED HIERARCHICAL CLUSTERING--------------------------------------------------------------------------


## Heatmap ##

dds <- estimateSizeFactors(dds) # adjust for library sizes

hip_norm_counts <- counts(dds, normalized = TRUE)

hip_vsd_cor <- vst(dds, blind = TRUE) %>% #log transform
  assay() %>% # extract matrix of transformed counts
  cor() # compute correlation values between samples
library(pheatmap)
# plot correlation matrix as heatmap
pheatmap(hip_vsd_cor, annotation = select(metadata, Group))

---------------------
  
  ## PCA analysis ##
  
# log transform the normalised counts two ways
hip_counts_log <- vst(dds, blind = TRUE) # first way with vst()
dds_log <- rlog(dds, blind = T ) # second way with rlog()

# define function to plot any two PCs against each other. Object is log transform dds. Specify PCs you want by
# defining pc.x and pc.y. Remember to set intgroup

library(genefilter)
library(ggrepel)

plotPCA.any <- function (object, pc.x, pc.y, labels = FALSE, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  require(genefilter)
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  
  ## Select the PCAs and percentVar that you like instead of 1 and 2
  d <- data.frame(PCX = pca$x[, pc.x], PCY = pca$x[, pc.y], group = group, 
                  intgroup.df, sample = rownames(colData(object)))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pc.x:pc.y]
    return(d)
  }
  if (labels) {
    g <- ggplot(data = d, aes_string(x = "PCX", y = "PCY", color = "group")) + 
      geom_point(size = 3, alpha = 0.8) + 
      xlab(paste0("PC", pc.x, ": ", round(percentVar[pc.x] * 100), "% Variance")) + 
      ylab(paste0("PC", pc.y,  ": ", round(percentVar[pc.y] * 100), "% Variance")) + 
      coord_fixed() + 
      geom_text_repel(size=3, aes_string(label = "sample"), color = "black") +
      scale_color_discrete(name = intgroup)
    return (g)
  } else {
    p <- ggplot(data = d, aes_string(x = "PCX", y = "PCY", color = "group")) + 
      geom_point(size = 3, alpha = 0.8) + 
      xlab(paste0("PC", pc.x, ": ", round(percentVar[pc.x] * 100), "% Variance")) + 
      ylab(paste0("PC", pc.y,  ": ", round(percentVar[pc.y] * 100), "% Variance")) + 
      coord_fixed() + 
      scale_color_discrete(name = intgroup)
    return(p)
  }
}


# call plotPCA.any function to get ggplot objects
pc_12_sex <- plotPCA.any(hip_counts_log, pc.x = 1, pc.y = 2, intgroup="Sex", labels = F)
pc_12_age <- plotPCA.any(hip_counts_log, pc.x = 1, pc.y = 2, intgroup="AGE", labels = F)

# set returndata = TRUE to obtain data frames
p12_data <- plotPCA.any(hip_counts_log, pc.x = 1, pc.y = 2, intgroup="genotype", labels = T, returnData = T)
p12age_data <- plotPCA.any(hip_counts_log, pc.x = 1, pc.y = 2, intgroup="AGE", labels = T, returnData = T)


library(ggforce) # package for adding shapes around groups of interest
library(ggrepel) # for text labelling on graphs

# PCA plot pc1 vs pc2

gtex_PCA_sex <- pc_12_sex +
  theme_bw() +
  theme(legend.text = element_text(size = 16), 
        legend.title = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.key.size = unit(0.8, "cm")) 

ggsave("PCA_GTex_Sex.tiff", plot = gtex_PCA_sex, dpi = 300, device = "tiff")


###-----ADJUST DESIGN FORMULA & RUN DESEQ ANALYSIS-----------------------------------------------------------------

# Adjust design formula to account for main sources of variation identified through heatmap/PCA

dds <- DESeqDataSetFromMatrix(countData =  raw_counts,
                              colData =  metadata,
                              design = ~ Sex + Group)
# run DEseq function

dds_analysed <- DESeq(dds)

## Check dispersions ie. how well do raw counts fit to model##

# save HQ tiff in work dir
tiff("dispersions.tiff", units="in", width=5, height=5, res=300)
plotDispEsts(dds_analysed) # plot should show decreasing disp with increasing mean
dev.off()

--------------------### EXTRACT RESULTS ###------------------------------------

res <- results(dds_analysed, # extract results from dds object run through deseq() function
               name = "Group_young_vs_old", # use name arg to specify coefficient. Can find this through resultsNames(dds_analysed)
               alpha = 0.05, # sig level P < 0.05
               lfcThreshold = 0.32) # logFC threshold at > 1.25 raw fold change

## Shrink LogFCs for better DE estimation
res <- lfcShrink(dds = dds_analysed, # pass deseq-analysed dds object
                        coef = "Group_young_vs_old", # coef should be same as name arg above
                        type = 'apeglm', res = res, svalue = FALSE) # pass extracted results object 
#to res. Note svalue = F keeps p and padj values

###----------SUMMARISE RESULTS DATA ###----------------------------------------------

summary <- summary(res) # get number of up/down genes

# filter out DE genes at adjusted p val < 0.05
data.frame(res) # convert results to data frame

# get rid of version number in ensembl ids as it seems to be reducing mapping ability of biomart
# 88 ensembl ids are now duplicates. Summarise these. Leaves us with 56,156 unique ensembl ids

res_df <- as.data.frame(res) %>% 
  rownames_to_column(var = "ensembl") %>% 
  separate(ensembl, into = c("ensembl", "version"), sep = "\\.[:digit:]+") %>% 
  dplyr::select(-version) %>% 
  group_by(ensembl) %>% 
  summarise(baseMean = round(mean(baseMean)), log2FoldChange = median(log2FoldChange),
            lfcSE = max(lfcSE), pvalue = min(pvalue), padj = min(padj))

mapped_ids <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"), 
                filters = "ensembl_gene_id", values = res_df$ensembl, mart = ensembl)

colnames(mapped_ids) <- c("ensembl", "symbol", "entrez")

# Summarise duplicate ensembl ids in mapped ids

res_annotated <- mapped_ids %>% 
  group_by(ensembl) %>% 
  summarise(symbol = min(symbol), entrez = min(entrez)) %>% 
  left_join(res_df, ., by = "ensembl")
  
res_annotated <- read_tsv("results_all_genes.txt")
write_tsv(res_annotated, "results_all_genes.txt")

res_sig <- res_annotated %>% 
  filter(padj < 0.05) %>% # filter gene list
  arrange(padj)

write_tsv(res_sig, "results_P_0.05.txt")


#-----###EXPLORE RESULTS DATA###-------------

# Create MA plot
plotMA(res)

# save HQ tiff in work dir
tiff("MA_plot_GTex.tiff", units="in", width=5, height=5, res=300)
plotMA(res_shrunk) # plot should show significant genes across a wide range of mean values
dev.off()


# create new logical column for padj < 0.05 and column for labelling genes on volcano plot



volcano_data <- res_annotated %>% 
  arrange(padj) %>%
  mutate(threshold = padj < 0.05) %>% # tibble deletes rownames so must create column to store
  mutate(genelabs = ifelse(.$symbol %in% .$symbol[1:60],
                           as.character(.$symbol), "")) %>% # must convert to character class!
  filter(!is.na(padj)) %>% 
  mutate(up_down = ifelse(log2FoldChange < 0 & padj < 0.05, "Down", 
                          ifelse(log2FoldChange > 0 & padj < 0.05, "Up", "Not Sig"))) #column for labelling sig up or down genes




# CREATE VOLCANO PLOT

library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(ggforce)

volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), 
                                          col = up_down)) + 
  geom_point(show.legend = T, alpha = 0.7) +
  theme_bw() +
  scale_color_brewer(palette = 2, type = "qual", name="FDR < 0.05",
                     breaks=c("Up", "Down", "Not Sig"),
                     labels=c("Upregulated", "Downregulated", "Not Significant")) +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = 2, size = 0.3) +
  geom_vline(xintercept = -0.32, color = "black", linetype = 2, size = 0.3) +
  geom_vline(xintercept = 0.32, color = "black", linetype = 2, size = 0.3) +
  geom_text_repel(aes(label = genelabs), size = 3.5, color = "black", point.padding = 0.22) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15)) +
  labs(x = "Log2 Fold Change", y = "-log10 (Adjusted P Value)") 

ggsave("volcano_plot.tiff", plot = volcano_plot, dpi = 300, device = "tiff")


### CREATE HEATMAP OF SIGNIFICANT GENES ###

library(pheatmap)

# Subset normalized counts to significant genes

hip_norm_counts <- as.data.frame(hip_norm_counts) # coerce to df

hip_norm_counts_0.05 <- hip_norm_counts[res_sig$ensembl, ] # subset rownames of normalised counts with significant genes

# Where there is no symbol or "NA" fill in with ENSG identifier
gaps <- which(res_sig$symbol[1:41] == "")
na <- which(is.na(res_sig$symbol))

res_sig$symbol[gaps] <- res_sig$ensembl[gaps]
res_sig$symbol[na] <- res_sig$ensembl[na]


row.names(hip_norm_counts_0.05) <- res_sig$symbol[1:41]

nrow(hip_norm_counts_0.05) == nrow(res_sig) # check this evaluates to true

hip_norm_counts_0.05 <- as.matrix(hip_norm_counts_0.05) # get rid of genes column

# Choose heatmap color palette

heat_colors <- brewer.pal(n = 6, name = "PuBu") # or go with default palette in pheatmap

ann_colors = list(Genotype = c(Progerin = "#F1948A", Wild_Type ="#16A085"), # pass list of colors for groups to annotation colors
                  Age = c(Old = "#F4D03F", Young = "#5D6D7E"))

# Plot heatmap
heatmap <- pheatmap(hip_norm_counts_0.05, 
                    cluster_cols = T, 
                    show_rownames = TRUE,
                    show_colnames = F,
                    border_color = FALSE,
                    annotation = dplyr::select(metadata, Group, Sex), 
                    scale = "row",
                    angle_col = 45,
                    legend = T,
                    annotation_legend = T,
                    annotation_names_col = F,
                    fontsize_col = 10,
                    main = "Heatmap of Normalised Counts of Significant DE genes", 
                    treeheight_row = 100,
                    treeheight_col = 30)

tiff("DE_genes_heatmap.tiff", units="in", width=8, height=8, res=300)
heatmap # plot should show decreasing disp with increasing mean
dev.off()




