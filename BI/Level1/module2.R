setwd("/Volumes/SamsungUSB/STEM-Away/BI-mentoring/Level1/SA-BI-Project")

library(GEOquery)
library(affy)

gse32323 <- getGEO(filename="Metadata/GSE32323_series_matrix.txt")
gse8671 <- getGEO(filename="Metadata/GSE8671_series_matrix.txt")

meta32323 <- gse32323@phenoData@data
meta8671<- gse8671@phenoData@data

meta <- data.frame(Sample=c(rownames(meta32323)[1:34], rownames(meta8671)),
                   Tissue=c(meta32323$`tissue:ch1`[1:34], meta8671$`Tissue:ch1`),
                   Batch=c(rep(1,34), rep(2,64)))

for (i in 1:length(meta$Tissue)) {
  meta$Tissue[i] <- strsplit(meta$Tissue, ",")[[i]][1]
  if (meta$Tissue[i] != "normal") {
    meta$Tissue[i] <- "cancer"
  }
}

meta$Group <- as.factor(c(rep("B1_Normal", 17), rep("B1_Cancer", 17), rep("B2_Normal", 32), rep("B2_Cancer", 32)))

# write.csv(meta, "Metadata/metadata.txt")
meta <- read.csv("Metadata/metadata.txt", row.names=1)


gse32323 <- ReadAffy(compress=T, celfile.path="Data/GSE32323/")
gse8671 <- ReadAffy(compress=T, celfile.path="Data/GSE8671/")
gse <- merge(gse32323, gse8671)
# saveRDS(gse, "Data/gse.rds")
gse <- readRDS("Data/gse.rds")
