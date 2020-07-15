# ----------### PROCESSING OF RAW COUNTS ###-------------


  # load packages

library(readr)
library(tidyverse)

### Creation of Metadata ###

genotype <- c("Wt", "Wt", "Progerin", "Progerin", 
              "Wt", "Wt", "Progerin", "Progerin")
age <- c("Old", "Old", "Old", "Old", "Young", "Young", "Young", "Young")
sample_id <- c("GSM1266739", "GSM1266740", "GSM1266741", "GSM1266742", 
               "GSM1266743", "GSM1266744", "GSM1266745", "GSM1266746")
names <- c("wt_old1", "wt_old2", "progerin_old1", "progerin_old2",
           "wt_young1", "wt_young2", "progerin_young1", "progerin_young2")

metadata <- data.frame(genotype,age,sample_id, stringsAsFactors = F)

rownames(metadata) <- names

### Read in Files ###

setwd("~/OneDrive - King's College London/Data/Bioinformatic/External/RNA-seq/Progerin_signatures/Miller_GSE524a31/GSE524a31_raw_counts")

read_tsv_plus <- function(flnm) { #function to label list ouput from map() with filenames
  read_tsv(flnm, col_names = F) %>% 
    mutate(filename = flnm) #adds a column to each df in each list element stating filename
}

files <- list.files(pattern = "*.txt") %>% # creates chr vector of txt files from work dir
  map(~read_tsv_plus(.)) %>% # iterates function over each file
  set_names(map_chr(., ~ .x[[3]][[1]])) # creates list names from first sub element of each element ie. the filename
  
### Rename and filter tibbles in list ###

files <- map(files, ~{names(.)[names(.) == 'X1'] <- 'gene_id'; .}) # rename first column of each tibble to gene_id

files <- purrr::map2(files, as.list(names), ~{names(.x)[names(.x) == 'X2'] <- .y; .x}) # apply sample names to second column
                                                                                # of each df in correct order
files <- map(files, ~ {.x[["filename"]] <- NULL; .x}) # get rid of filename column in each df

### Merge tibbles in list together ###

merged_tibs <- files %>% 
  purrr::reduce(left_join, by = "gene_id") %>% # reduce will iteratively merge tibbles together element-wise
  as.data.frame() # need as data frame  as tibbles not allowed rownames

rownames(merged_tibs) <- merged_tibs$gene_id # set genes as rownames in raw counts
raw_counts <- merged_tibs[, -1] # remove gene_id column from raw counts

if(!all(rownames(metadata) == colnames(raw_counts))) {
  stop("rownames of metadata and colnames of raw counts do not match")
} #should evauluate to TRUE

###----------------CREATE DESEQ DATASET OBJECT--------------------------------------


library(DESeq2)
library(ggplot2)
library(ggfortify)

  # create the deseq object and specify experimental design
dds <- DESeqDataSetFromMatrix(countData =  raw_counts,
                              colData =  metadata,
                              design = ~ genotype)

###-----UNSUPERVISED HIERARCHICAL CLUSTERING--------------------------------------------------------------------------


  ## Heatmap ##

dds <- estimateSizeFactors(dds) # adjust for library sizes

prog_norm_counts <- counts(dds, normalized = TRUE)

vsd_cor_prog <- vst(dds, blind = TRUE) %>% #log transform
  assay() %>% # extract matrix of transformed counts
  cor() # compute correlation values between samples
library(pheatmap)
# plot correlation matrix as heatmap
pheatmap(vsd_cor_prog, annotation = select(metadata, genotype))

---------------------

  ## PCA analysis ##

# log transform the normalised counts two ways
norm_counts_log <- vst(dds, blind = TRUE) # first way with vst()
dds_log <- rlog(dds, blind = T ) # second way with rlog()

# define function to plot any two PCs against each other. Object is log transform dds. Specify PCs you want by
# defining pc.x and pc.y. Remember to set intgroup

library(genefilter)

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
      geom_point(size = 3) + 
      xlab(paste0("PC", pc.x, ": ", round(percentVar[pc.x] * 100), "% Variance")) + 
      ylab(paste0("PC", pc.y,  ": ", round(percentVar[pc.y] * 100), "% Variance")) + 
      coord_fixed() + 
      geom_text_repel(size=3, aes_string(label = "sample"), color = "black") +
      scale_color_discrete(name = intgroup)
    return (g)
    } else {
      p <- ggplot(data = d, aes_string(x = "PCX", y = "PCY", color = "group")) + 
        geom_point(size = 3) + 
        xlab(paste0("PC", pc.x, ": ", round(percentVar[pc.x] * 100), "% Variance")) + 
        ylab(paste0("PC", pc.y,  ": ", round(percentVar[pc.y] * 100), "% Variance")) + 
        coord_fixed() + 
        scale_color_discrete(name = intgroup)
      return(p)
    }
}


# call plotPCA.any function to get ggplot objects
pc_12 <- plotPCA.any(dds_log, pc.x = 1, pc.y = 2, intgroup="genotype", labels = F)
pc_23 <- plotPCA.any(dds_log, pc.x = 2, pc.y = 3, intgroup="genotype", labels = F)

# set returndata = TRUE to obtain data frames
p12_data <- plotPCA.any(dds_log, pc.x = 1, pc.y = 2, intgroup="genotype", labels = T, returnData = T)
p23_data <- plotPCA.any(dds_log, pc.x = 2, pc.y = 3, intgroup="genotype", labels = T, returnData = T)


# split sample column to extract age column
p23_data <- p23_data %>%
  separate(sample, into = c(NA, "age"), sep = "_") %>% 
  separate(age, into = c("age", NA), sep = "(?<=[A-Za-z])(?=[0-9])")
p12_data <- p12_data %>% 
  separate(sample, into = c(NA, "age"), sep = "_") %>% 
  separate(age, into = c("age", NA), sep = "(?<=[A-Za-z])(?=[0-9])")

library(ggforce) # package for adding shapes around groups of interest
library(ggrepel) # for text labelling on graphs

# PCA plot pc1 vs pc2

PCA_12_final <- pc_12 +
  theme_bw() +
  geom_mark_ellipse(aes(fill = genotype, col = genotype), alpha = 0.3, show.legend = T, color = "white") +
  geom_text_repel(size=5, aes(label = age), color = "black", point.padding = 0.2) +
  guides(color = FALSE) +
  labs(fill = "Genotype") +
  theme(legend.text = element_text(size = 16), 
        legend.title = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.key.size = unit(1.5, "cm")) +
  scale_fill_discrete(labels = c("Progerin", "Wild Type"))

ggsave("PC1_vs_PC2.tiff", plot = PCA_12_final, dpi = 300, device = "tiff")


# PCA plot pc2 vs pc3

PCA_23_final <- pc_23 +
  theme_bw() +
  geom_mark_ellipse(data = p23_data, aes(group = age, fill = age), 
                    alpha = 0.3, show.legend = T, color = "white") + # pass data containing age column to this geom
  geom_text_repel(size=4, aes(label = sample), color = "black", point.padding = 0.25) +
  guides(color = FALSE) +
  labs(fill = "Age") +
  theme(legend.text = element_text(size = 16), 
        legend.title = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.key.size = unit(1.5, "cm")) +
  scale_fill_discrete(labels = c("Old", "Young"))

ggsave("PC2_vs_PC3.tiff", plot = PCA_23_final, dpi = 300, device = "tiff")



###-----ADJUST DESIGN FORMULA & RUN DESEQ ANALYSIS-----------------------------------------------------------------

# Adjust design formula to account for main sources of variation identified through heatmap/PCA

dds <- DESeqDataSetFromMatrix(countData =  raw_counts,
                              colData =  metadata,
                              design = ~ age + genotype)
# run DEseq function

dds_analysed <- DESeq(dds)
  
## Check dispersions ie. how well do raw counts fit to model##

    # save HQ tiff in work dir
tiff("dispersions.tiff", units="in", width=5, height=5, res=300)
plotDispEsts(dds_analysed) # plot should show decreasing disp with increasing mean
dev.off()

--------------------### EXTRACT RESULTS ###------------------------------------
 
res <- results(dds_analysed, # extract results from dds object run through deseq() function
                     name = "genotype_Wt_vs_Progerin", # use name arg to specify coefficient
                     alpha = 0.05, # sig level P < 0.05
                     lfcThreshold = 0.32) # logFC threshold at > 1.25 raw fold change

  ## Shrink LogFCs for better DE estimation
res_shrunk <- lfcShrink(dds = dds_analysed, # pass deseq-analysed dds object
                     coef = "genotype_Wt_vs_Progerin", # coef should be same as name arg above
                     type = 'apeglm', res = res, svalue = FALSE) # pass extracted results object 
                                                #to res. Note svalue = F keeps p and padj values

###----------SUMMARISE RESULTS DATA ###----------------------------------------------

summary(res_shrunk) # get number of up/down genes

  # filter out DE genes at adjusted p val < 0.05
res_shrunk_df <- data.frame(res_shrunk) # convert results to data frame

res_shrunk_df2 <- res_shrunk_df %>% 
  rownames_to_column(var = "Gene")

genes <- res_shrunk_df2$Gene

entrez <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), 
      filters = "external_gene_name", values = genes, mart = ensembl)

colnames(entrez) <- c("Gene", "Ensembl")

merged <- left_join(res_shrunk_df2, entrez, by = "Gene")

idx <- which(isUnique(merged$Gene))

merged <- view(merged[idx, ])

table(isUnique(merged$Gene))

write_tsv(merged, "results_all_genes.txt")

res_sig <- subset(res_shrunk_df2, padj < 0.05) # filter gene list

genes <- rownames(res_sig) # grab all significant genes

res_sig <- res_sig %>% 
  mutate(Gene = genes) %>% # create column for sig genes as rownames are lost with tibble
  arrange(padj) # order adj p values smallest to largest

res_sig <- res_sig[, c(6, 1:5)] # reorder columns for clarity

head(res_sig, 20) # look at top 20 DE genes



#-----###EXPLORE RESULTS DATA###-------------

  # Create MA plot
plotMA(res_shrunk)

# save HQ tiff in work dir
tiff("MA_plot.tiff", units="in", width=5, height=5, res=300)
plotMA(res_shrunk) # plot should show significant genes across a wide range of mean values
dev.off()


# create new logical column for padj < 0.05 and column for labelling genes on volcano plot

genes_all <- rownames(res)

prog_res_all <- data.frame(res_shrunk) %>% 
  mutate(threshold = padj < 0.05, 
         gene = genes_all) %>% # tibble deletes rownames so must create column to store
  arrange(padj)

table(prog_res_all$gene == "alignment_not_unique") # find how many genes display this message (only 1 in this case)
prog_res_all <- prog_res_all[-3, ] # get rid of this unidentified gene in third row

prog_res_all <- prog_res_all %>% # create column for labelling top 20 DE genes
  mutate(genelabs = ifelse(prog_res_all$gene %in% prog_res_all$gene[1:10]|prog_res_all$gene == "LMNA"|
                             prog_res_all$gene == "LIN28A"| prog_res_all$gene == "FEZF2",
                           as.character(prog_res_all$gene), "")) # must convert to character class!

prog_res_all <- prog_res_all %>% # filter out rows with no multiple corrected p value
  filter(!is.na(padj))


prog_res_all <- prog_res_all %>%
  mutate(up_down = ifelse(log2FoldChange < 0 & padj < 0.05, "Down", 
                          ifelse(log2FoldChange > 0 & padj < 0.05, "Up", "Not Sig")))



# CREATE VOLCANO PLOT

library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(ggforce)

volcano_plot2 <- ggplot(prog_res_all, aes(x = log2FoldChange, y = -log10(padj), 
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
  geom_mark_ellipse(aes(filter = log2FoldChange < -3, 
                        description = "Only 3 genes are found to be significantly downregulated"), 
                    color = "#2DC59E", size = 0.8, label.fontsize = 10, linetype = 2, expand = unit(7.8, "mm")) +
  geom_text_repel(aes(label = genelabs), size = 3, color = "black", point.padding = 0.22) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15)) +
  labs(x = "Log2 Fold Change", y = "-log10 (Adjusted P Value)") 
  
ggsave("volcano_plot2.tiff", plot = volcano_plot2, dpi = 300, device = "tiff")


### CREATE HEATMAP OF SIGNIFICANT GENES ###

library(pheatmap)

# Subset normalized counts to significant genes

prog_norm_counts <- as.data.frame(prog_norm_counts) # coerce to df

prog_norm_sig_counts <- prog_norm_counts[prog_res_sig$Gene, ] # subset rownames of normalised counts with significant genes

nrow(prog_norm_sig_counts) == nrow(prog_res_sig) # check this evaluates to true

prog_norm_sig_counts <- as.matrix(prog_norm_sig_counts[, 1:8]) # get rid of genes column

# Choose heatmap color palette

heat_colors <- brewer.pal(n = 6, name = "PuBu") # or go with default palette in pheatmap

ann_colors = list(Genotype = c(Progerin = "#F1948A", Wild_Type ="#16A085"), # pass list of colors for groups to annotation colors
                  Age = c(Old = "#F4D03F", Young = "#5D6D7E"))

heat_metadata <- metadata # create new metadata df to edit names (can't edit labels of pheatmap easily)

colnames(heat_metadata) <- c("Genotype", "Age", "Sample ID")

heat_metadata$Genotype[heat_metadata$Genotype == "Wild Type"] <- "Wild_Type"

samples <- rownames(heat_metadata)

samples <- c("Wild Type | Old 1", "Wild Type | Old 2", "Progerin | Old 1", 
             "Progerin | Old 2", "Wild Type | Young 1", "Wild Type | Young 2",
             "Progerin | Young 1", "Progerin | Young 2")

rownames(heat_metadata) <- samples # give readable sample names to metadata

colnames(prog_norm_sig_counts) <- samples # give matching sample names to filtered norm counts

# Plot heatmap
heatmap <- pheatmap(prog_norm_sig_counts, 
         cluster_rows = TRUE, 
         show_rownames = FALSE,
         border_color = FALSE,
         annotation = select(heat_metadata, Genotype, Age), 
         scale = "row",
         angle_col = 45,
         legend = T,
         annotation_legend = T,
         annotation_names_col = F,
         fontsize_col = 10,
         main = "Heatmap of Normalised Counts of Significant DE genes", 
         treeheight_row = 100,
         treeheight_col = 30,
         annotation_colors = ann_colors)

tiff("DE_genes_heatmap.tiff", units="in", width=5, height=5, res=300)
heatmap # plot should show decreasing disp with increasing mean
dev.off()


