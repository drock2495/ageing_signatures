

----------# COMBINED AGEING SIGNATURE FOR BERCHTOLD & KANG ----------


library(tidyverse)

setwd("~/OneDrive - King's College London/Data/Bioinformatic/External/Ageing Datasets/Hippocampus/results")

# Load DE genes from all ageing datasets

berchtold <- read_tsv("diffExpr_all_genes_berchtold.txt")

kang <- read_tsv("diffExpr_all_genes_kang2011.txt")


# Find common entrez ids between all studies

common <- intersect(berchtold$entrezgene_id, kang$entrez)

# There are 13893 common entrez genes between datasets


berchtold <- berchtold[match(common, berchtold$entrezgene_id),]
kang <- kang[match(common, kang$entrez),]


# Now we have same entrez in same order for each dataset we can cbind them together

colnames(berchtold)[c(3, 6, 7)] <- paste0("berch", "_" ,colnames(berchtold)[c(3, 6, 7)])
colnames(kang)[c(7,10,11)] <- paste0("kang", "_" ,colnames(kang)[c(7,10,11)])

ber_kang <- cbind(berchtold[, c("hgnc_symbol", "entrezgene_id")], berchtold[, c(3, 6, 7)],
                       kang[c(7,10,11)])

# Save combined signature for IPA analysis

setwd("~/OneDrive - King's College London/Data/Bioinformatic/IPA/Data")

ber_kang %>% 
  filter(berch_logFC > 0 & kang_logFC > 0 |
           berch_logFC < 0 & kang_logFC < 0) %>%
  dplyr::select(hgnc_symbol, entrezgene_id, berch_logFC, kang_logFC) %>%
  pivot_longer(cols = -c(hgnc_symbol, entrezgene_id), names_to = "study", values_to = "logFC") %>% 
  group_by(hgnc_symbol) %>% 
  mutate(meanlogFC = mean(logFC)) %>% 
  pivot_wider(names_from = "study", values_from = "logFC") %>% 
  ungroup() %>% 
  arrange(desc(meanlogFC)) %>% 
  select(hgnc_symbol,entrezgene_id,meanlogFC) %>% 
  write_tsv(.,"Berch_Kang_Combined_Signature.txt")


