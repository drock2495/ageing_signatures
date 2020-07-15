--------### MODEL FITTING & GET RESULTS###-------------------

## Load packages ##

library(Biobase)
library(limma)
library(statmod)

## Read in data ##

source("~/OneDrive - King's College London/Dan Array Data/Data for Limma/TYTUS/Data/probeset_summarisation_data/design_matrix.R")

## Estimate correlation between measurements made on the same subject ##

intra_subject <- pData(eset)[, "Sample"]
intra_subject <- factor(intra_subject)
corfit <- duplicateCorrelation(eset, design, block = intra_subject)
corfit$consensus

## Then this inter-subject correlation is input into the linear model fit ##

fit <- lmFit(eset,design,block=intra_subject,correlation=corfit$consensus)

## set contrasts we wish to make ##

cm <- makeContrasts(OldVYoungforTP1 = Old.1 - Young.1, OldVYoungforTP6 = Old.6 - Young.6, 
                    OldVYoungforTP24 = Old.24 - Young.24, OldVYoungforTP72 = Old.72 - Young.72, 
                    OldVYoungforTP144 = Old.144 - Young.144,
                    levels = design)

## compute these contrasts and moderated t tests ##

fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)

## find top DE genes between contrast coefficients ##

results <- decideTests(fit2)
summary(results) # only 3 DEGs detected at padj < 0.05!

## obtain DEGs for full fitted model ranked by p value ##

all_genes <- topTable(fit2, number = nrow(fit2))

## results for top DE genes at each timepoint ##

OldvYoungTP144 <- topTable(fit2, coef="OldVYoungforTP144", number = nrow(fit2))
OldvYoungTP72 <- topTable(fit2, coef="OldVYoungforTP72", number = nrow(fit2))
OldvYoungTP24 <- topTable(fit2, coef="OldVYoungforTP24", number = nrow(fit2))
OldvYoungTP6 <- topTable(fit2, coef="OldVYoungforTP6", number = nrow(fit2))
OldvYoungTP1 <- topTable(fit2, coef="OldVYoungforTP1", number = nrow(fit2))

setwd("~/OneDrive - King's College London/Dan Array Data/Data for Limma/TYTUS/TopDE genes/TopDEGs/bycoeff_contrast/probeSet_summarised_before_LmFit")

write_tsv(OldvYoungTP1, "OldvYoungTP1.txt", na = "NA")
write_tsv(OldvYoungTP6, "OldvYoungTP6.txt", na = "NA")
write_tsv(OldvYoungTP24, "OldvYoungTP24.txt", na = "NA")
write_tsv(OldvYoungTP72, "OldvYoungTP72.txt", na = "NA")
write_tsv(OldvYoungTP144, "OldvYoungTP144.txt", na = "NA")

## filter genes at p < 0.05 

p0.05 <- all_genes %>% 
  filter(P.Value < 0.05) 
write_tsv(p0.05, "DEGs_p0.05.txt", na = "NA")


