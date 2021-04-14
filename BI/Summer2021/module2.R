setwd('/Volumes/SamsungUSB/STEM-Away/BI-mentoring/Level1/2021Summer')

library(GEOquery)
library(affy)

gse <- getGEO(filename='metadata/GSE19804_series_matrix.txt')

meta <- gse@phenoData@data

meta <- data.frame(Sample=rownames(meta),
                   Tissue=meta$`tissue:ch1`)

for (i in 1:length(meta$Tissue)) {
  meta$Tissue[i] <- strsplit(meta$Tissue, ',')[[i]][1]
  if (meta$Tissue[i] == 'lung cancer') {
    meta$Tissue[i] <- 'CANCER'
  } else {
    meta$Tissue[i] <- 'NORMAL'
  }
}

# write.csv(meta, 'metadata/metadata.txt')
meta <- read.csv('metadata/metadata.txt', row.names=1)

gse <- ReadAffy(compress=T, celfile.path='Data/GSE19804/')
# saveRDS(gse, 'data/gse.rds')
gse <- readRDS('data/gse.rds')
