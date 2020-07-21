setwd("C:/Users/durbe/OneDrive/Documents/ReposExtras/STEM-Away/July2020/")

library(affy)
library(affyPLM)
library(simpleaffy)
library(arrayQualityMetrics)
library(affyQCReport)
library(gcrma)
library(sva)
library(ggplot2)
library(pheatmap)


##### WEEK TWO - WEEK THREE #####
gse32323 <- ReadAffy(compress=T, celfile.path="./Data/Raw/GSE32323/")
gse8671 <- ReadAffy(compress=T, celfile.path="./Data/Raw/GSE8671/") 
gse <- merge(gse32323, gse8671)
row.names(gse@phenoData@data) <- sub(".CEL.gz", "", row.names(gse@phenoData@data))
row.names(gse@phenoData@data) <- sub("_chip_array", "", row.names(gse@phenoData@data))
row.names(gse@protocolData@data) <- sub(".CEL.gz", "", row.names(gse@protocolData@data))
row.names(gse@protocolData@data) <- sub("_chip_array", "", row.names(gse@protocolData@data))
saveRDS(gse, "./Data/gse.rds")


### QUALITY CONTROL - before normalization
# simpleaffy
sa <- qc(gse)
plot(sa)
# arrayQualityMetrics
arrayQualityMetrics(gse, "./QC/arrayQualityMetrics-Raw/", force=T, do.logtransform=T)
# affyQCReport
QCReport(gse, "./QC/affyQCReport-Raw.pdf")
# affyPLM
pset <- fitPLM(gse, background=F, normalize=T)
rle <- RLE(pset, type="stats")
rle <- data.frame(Median=rle[1,])
ggplot(rle, aes(Median))+ geom_histogram()+ xlab("Median")+ ylab("Frequency")+ ggtitle("RLE Histogram (Raw)")+ theme_bw()
par(mai=c(3.5,1,1,1))
RLE(pset, main="RLE (Raw)", las=2, ce.lab=0.5)
nuse <- NUSE(pset, type="stats")
nuse <- data.frame(Median=nuse[1,])
ggplot(nuse, aes(Median))+ geom_histogram()+ xlab("Median")+ ylab("Frequency")+ ggtitle("NUSE Histogram (Raw)")+ theme_bw()
par(mai=c(3.5,1,1,1))
NUSE(pset, main="NUSE (Raw)", las=2)


### NORMALIZATION
# raw
raw <- as.data.frame(exprs(gse))
par(mai=c(3.5,1,1,1))
boxplot(gse, main="Boxplot Raw Data", ylab="Probe Intensities", las=2)
# mas5
mas5 <- mas5(gse)
par(mai=c(3.5,1,1,1))
boxplot(mas5, main="Boxplot mas5 Data", ylab="Probe Intensities", las=2)
# rma
rma <- rma(gse)
par(mai=c(3.5,1,1,1))
boxplot(exprs(rma), main="Boxplot rma Data", ylab="Probe Intensities", las=2)
# gcrma
gcrma <- gcrma(gse)
par(mai=c(3.5,1,1,1))
boxplot(exprs(gcrma), main="Boxplot gcrma Data", ylab="Probe Intensities", las=2)


### QUALITY CONTROL - after normalization
# arrayQualityMetrics
arrayQualityMetrics(mas5, "./QC/arrayQualityMetrics-Mas5/", force=T, do.logtransform=T)
arrayQualityMetrics(rma, "./QC/arrayQualityMetrics-Rma/", force=T, do.logtransform=T)
arrayQualityMetrics(gcrma, "./QC/arrayQualityMetrics-Gcrma/", force=T, do.logtransform=T)
# affyQCReport
QCReport(mas5, "./QC/affyQCReport-Mas5.pdf")
QCReport(rma, "./QC/affyQCReport-Rma.pdf")
QCReport(gcrma, "./QC/affyQCReport-Gcrma.pdf")
# affyPLM
pset <- fitPLM(gse, background=TRUE, normalize=TRUE) #mas5
pset <- fitPLM(gse, background=TRUE, normalize=TRUE, background.method="RMA")
pset <- fitPLM(gse, background=TRUE, normalize=TRUE, background.method="GCRMA")
rle <- RLE(pset, type="stats")
rle <- data.frame(Median=rle[1,])
ggplot(rle, aes(Median))+ geom_histogram()+ xlab("Median")+ ylab("Frequency")+ ggtitle("RLE Histogram (Mas5)")+ theme_bw()
par(mai=c(3.5,1,1,1))
RLE(pset, main="RLE (Mas5)", las=2)
nuse <- NUSE(pset, type="stats")
nuse <- data.frame(Median=nuse[1,])
ggplot(nuse, aes(Median))+ geom_histogram()+ xlab("Median")+ ylab("Frequency")+ ggtitle("NUSE Histogram (Mas5)")+ theme_bw()
par(mai=c(3.5,1,1,1))
NUSE(pset, main="NUSE (Mas5)", las=2)


### BATCH CORRECTION
pheno <- read.csv("./Data/metadata.csv")
type <- pheno$CN
row.names(pheno) <- pheno$X.Sample_geo_accession
design <- model.matrix(~0+type, pheno)

combatMas5 <- ComBat(exprs(mas5), pheno$BATCH, design[,1])
combatRma <- ComBat(exprs(rma), pheno$BATCH, design[,1])
combatGcrma <- ComBat(exprs(gcrma), pheno$BATCH, design[,1])
write.csv(combatMas5, "./Data/BatchCorrected/mas5.csv")
write.csv(combatRma, "./Data/BatchCorrected/rma.csv")
write.csv(combatGcrma, "./Data/BatchCorrected/gcrma.csv")


### VISUALIZATION
## PCA
group <- as.factor(pheno$mygroup)
# raw
pca <- prcomp(raw, scale=F, center=F)
pca <- as.data.frame(pca$rotation)
ggplot(pca, aes(x=PC1, y=PC2, color=group))+
  geom_point()+
  scale_colour_manual(values=c("hotpink1", "cyan4", "coral1", "aquamarine4"))+
  ggtitle("PCA before Normalization")+ theme_bw()
# mas5
pca <- prcomp(exprs(mas5), scale=F, center=F)
pca <- as.data.frame(pca$rotation)
ggplot(pca, aes(x=PC1, y=PC2, color=group))+
  geom_point()+
  scale_colour_manual(values=c("hotpink1", "cyan4", "coral1", "aquamarine4"))+
  ggtitle("PCA mas5 Normalization")+ theme_bw()
# rma
pca <- prcomp(exprs(rma), scale=F, center=F)
pca <- as.data.frame(pca$rotation)
ggplot(pca, aes(x=PC1, y=PC2, color=group))+
  geom_point()+
  scale_colour_manual(values=c("hotpink1", "cyan4", "coral1", "aquamarine4"))+
  ggtitle("PCA rma Normalization")+ theme_bw()
# gcrma
pca <- prcomp(exprs(gcrma), scale=F, center=F)
pca <- as.data.frame(pca$rotation)
ggplot(pca, aes(x=PC1, y=PC2, color=group))+
  geom_point()+
  scale_colour_manual(values=c("hotpink1", "cyan4", "coral1", "aquamarine4"))+
  ggtitle("PCA gcrma Normalization")+ theme_bw()
# mas5 - batch correction
pca <- prcomp(combatMas5, scale=F, center=F)
pca <- as.data.frame(pca$rotation)
ggplot(pca, aes(x=PC1, y=PC2, color=group))+
  geom_point()+
  scale_colour_manual(values=c("hotpink1", "cyan4", "coral1", "aquamarine4"))+
  ggtitle("PCA mas5 Batch Correction")+ theme_bw()
# rma - batch correction
pca <- prcomp(combatRma, scale=F, center=F)
pca <- as.data.frame(pca$rotation)
ggplot(pca, aes(x=PC1, y=PC2, color=group))+
  geom_point()+
  scale_colour_manual(values=c("hotpink1", "cyan4", "coral1", "aquamarine4"))+
  ggtitle("PCA rma Batch Correction")+ theme_bw()
# gcrma - batch correction
pca <- prcomp(combatGcrma, scale=F, center=F)
pca <- as.data.frame(pca$rotation)
ggplot(pca, aes(x=PC1, y=PC2, color=group))+
  geom_point()+
  scale_colour_manual(values=c("hotpink1", "cyan4", "coral1", "aquamarine4"))+
  ggtitle("PCA gcrma Batch Correction")+ theme_bw()

## HEATMAP
# mas5
dismat <- 1-cor(exprs(mas5))
colnames(dismat) <- pheno$CN
row.names(dismat) <- sub(".CEL.gz", "", row.names(dismat))
row.names(dismat)[1:34] <- sub("_chip_array", "", row.names(dismat)[1:34])
pheatmap(dismat, main="Hierarchical Clustering Heatmap mas5 Normalization")
# rma
dismat <- 1-cor(exprs(rma))
colnames(dismat) <- pheno$CN
row.names(dismat) <- sub(".CEL.gz", "", row.names(dismat))
row.names(dismat)[1:34] <- sub("_chip_array", "", row.names(dismat)[1:34])
pheatmap(dismat, main="Hierarchical Clustering Heatmap rma Normalization")
# gcrma
dismat <- 1-cor(exprs(gcrma))
colnames(dismat) <- pheno$CN
row.names(dismat) <- sub(".CEL.gz", "", row.names(dismat))
row.names(dismat)[1:34] <- sub("_chip_array", "", row.names(dismat)[1:34])
pheatmap(dismat, main="Hierarchical Clustering Heatmap gcrma Normalization")
# mas5 - batch correction
dismat <- 1-cor(combatMas5)
colnames(dismat) <- pheno$CN
row.names(dismat) <- sub(".CEL.gz", "", row.names(dismat))
row.names(dismat)[1:34] <- sub("_chip_array", "", row.names(dismat)[1:34])
pheatmap(dismat, main="Hierarchical Clustering Heatmap mas5 Batch Correction")
# rma - batch correction
dismat <- 1-cor(combatRma)
colnames(dismat) <- pheno$CN
row.names(dismat) <- sub(".CEL.gz", "", row.names(dismat))
row.names(dismat)[1:34] <- sub("_chip_array", "", row.names(dismat)[1:34])
pheatmap(dismat, main="Hierarchical Clustering Heatmap rma Batch Correction")
# gcrma - batch correction
dismat <- 1-cor(combatGcrma)
colnames(dismat) <- pheno$CN
row.names(dismat) <- sub(".CEL.gz", "", row.names(dismat))
row.names(dismat)[1:34] <- sub("_chip_array", "", row.names(dismat)[1:34])
pheatmap(dismat, main="Hierarchical Clustering Heatmap gcrma Batch Correction")


















