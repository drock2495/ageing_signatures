## Load packages ##

library(Biobase)
library(limma)
library(Rcpp)
library(tidyverse)

## Read in data ##

setwd("~/OneDrive - King's College London/Dan Array Data/Data for Limma/TYTUS/Data/probeset_summarisation_data")

first_col_to_rownames <- function(df) {
  names <- df[[1]]
  df <- as.data.frame(df)
  rownames(df) <- names
  df2 <- df[, -1]
  df2
}

p <- read_tsv("phenotype_data.txt")

p <- first_col_to_rownames(p)

f <- read_tsv("feature_data.txt")

f <- first_col_to_rownames(f)

x <- read_tsv("expression_matrix.txt")

x <- first_col_to_rownames(x)


## convert expression matrix to matrix class ##

x <- as.matrix(x)


## Create espressionset object ##

eset <- ExpressionSet(assayData = x, phenoData = AnnotatedDataFrame(p), featureData = AnnotatedDataFrame(f))

## combine age and timepoint into combined factor ##

YngOld_tmpt_comb <- with(pData(eset), paste(Type, Timepoint, sep = "."))

YngOld_tmpt_comb <- factor(YngOld_tmpt_comb)

## Create variables to input into the design matrix ##

Batch <- as.factor(pData(eset)[, "Batch"])
Gender <- as.factor(pData(eset)[, "Gender"])


# Create design matrix without intercept but include combined age.timepoint variable and batch variable in the linear model ##

design <- model.matrix(~0 + YngOld_tmpt_comb + Gender + Batch)

## rename colnames of design matrix to more interpretable names ##

names <- c("Old.1", "Old.144", "Old.24", "Old.6", "Old.72", "Young.1", 
           "Young.144", "Young.24", "Young.6", "Young.72", "Gender.1", "Batch.2", "Batch.3")

colnames(design) <- names


# A little info regarding batch effects - source https://www.biostars.org/p/266507/

#It should be pointed out that there is a difference between modelling a batch effect and 
#directly modifying data in order to counteract a batch effect. 
#Programs like ComBat aim to directly modify your data in an attempt to 
#eliminate batch effects (it literally 'subtracts' out the modeled effect, 
#which can result in the infamous negative values after employing ComBat). 
#After employing ComBat, statistical tests are conducted on the modified data,
#with batch not appearing in the design formula.

#However, by including 'batch' as a covariate in a design formula for the purpose of differential expression analysis, 
#one is simply modeling the effect size of the batch, without actually modifying your data. 
#The modeled effect size is then taken into account when statistical inferences are made. 
#In DESeq2, however, the vst() and rld() transformations can adjust the actual transformed counts based on batch 
##(and anything else in your design formula) by setting blind=FALSE, and this is the recommended procedure by 
#Michael Love when the aim is to use the normalised+transformed expression levels for downstream analyses.

