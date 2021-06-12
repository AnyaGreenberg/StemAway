setwd('D:/stemaway/bioinformatics/level2')

library(affy)
library(affyPLM)
library(arrayQualityMetrics)
library(ggplot2)
library(pheatmap)

gse <- readRDS('data/affybatch.rds')


### AFFYPLM
pset <- fitPLM(gse, background=TRUE, normalize=TRUE)

# Relative log expression
par(mai=c(1.5,1,0.5,0.5))
RLE(pset, main='RLE', ylab='Probe Intensities', las=2)

# Normalized unscaled standard error
par(mai=c(1.5,1,0.5,0.5))
NUSE(pset, main='NUSE', ylab='Probe Intensities', las=2)


### ARRAYQUALITYMETRICS
arrayQualityMetrics(gse, 'output/arrayQualityMetrics/', force=TRUE, do.logtransform = TRUE)


### NORMALIZATION
rma <- rma(gse)
rma_df <- exprs(rma)
