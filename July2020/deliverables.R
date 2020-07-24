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


### QUALITY CONTROL - before normalization
arrayQualityMetrics(gse, "./QC/arrayQualityMetrics/", force=T, do.logtransform=T)

# simpleaffy
sa <- qc(gse)
plot(sa)

# affyQCReport
QCReport(gse, "./QC/affyQCReport.pdf")

# affyPLM
pset <- fitPLM(gse, background=T, normalize=T)
rle <- RLE(pset, type="stats")
rle <- data.frame(Median=rle[1,])
ggplot(rle, aes(Median))+ geom_histogram()+ xlab("Median")+ ylab("Frequency")+ ggtitle("RLE Histogram")+ theme_bw()
par(mai=c(3.5,1,1,1))
RLE(pset, main="RLE", las=2, ce.lab=0.5)

nuse <- NUSE(pset, type="stats")
nuse <- data.frame(Median=nuse[1,])
ggplot(nuse, aes(Median))+ geom_histogram()+ xlab("Median")+ ylab("Frequency")+ ggtitle("NUSE Histogram")+ theme_bw()
par(mai=c(3.5,1,1,1))
NUSE(pset, main="NUSE", las=2)


### NORMALIZATION
box_plot <- function(d, t) {
  par(mai=c(3.5,1,1,1))
  boxplot(d, main=t, ylab="Probe Intensities", las=2)
}

box_plot(gse, "Boxplot Raw Data")

mas5 <- mas5(gse)
par(mai=c(3.5,1,1,1))
box_plot(exprs(mas5), "Boxplot mas5 Data")

rma <- rma(gse)
box_plot(exprs(rma), "Boxplot rma Data")

gcrma <- gcrma(gse)
box_plot(exprs(gcrma), "Boxplot gcrma Data")

### BATCH CORRECTION
pheno <- read.csv("./Data/metadata.csv")
type <- pheno$CN
row.names(pheno) <- pheno$X.Sample_geo_accession
design <- model.matrix(~0+type, pheno)

combatMas5 <- ComBat(exprs(mas5), pheno$BATCH, design[,1])
combatRma <- ComBat(exprs(rma), pheno$BATCH, design[,1])
combatGcrma <- ComBat(exprs(gcrma), pheno$BATCH, design[,1])


### VISUALIZATION
## PCA
group <- as.factor(pheno$mygroup)

pca_plot <- function(d, t) {
  pca <- prcomp(d, scale=F, center=F)
  pca <- as.data.frame(pca$rotation)
  ggplot(pca, aes(x=PC1, y=PC2, color=group))+
    geom_point()+
    scale_colour_manual(values=c("hotpink1", "navy", "goldenrod2", "aquamarine4"))+
    ggtitle(t)+ theme_bw()
}

pca_plot(exprs(gse), "PCA before Normalization")
pca_plot(t(mas5), "PCA mas5 Normalization")
pca_plot(t(rma), "PCA rma Normalization")
pca_plot(t(gcrma), "PCA gcrma Normalization")
pca_plot(combatMas5, "PCA mas5 Batch Correction")
pca_plot(combatRma, "PCA rma Batch Correction")
pca_plot(combatGcrma, "PCA gcrma Batch Correction")


## HEATMAP
heatmap <- function(d,t) {
  dismat <- 1-cor(d)
  colnames(dismat) <- pheno$CN
  row.names(dismat) <- sub(".CEL.gz", "", row.names(dismat))
  row.names(dismat)[1:34] <- sub("_chip_array", "", row.names(dismat)[1:34])
  pheatmap(dismat, main=t)
}

heatmap(t(mas5), "Hierarchical Clustering Heatmap mas5 Normalization")
heatmap(t(rma), "Hierarchical Clustering Heatmap rma Normalization")
heatmap(t(gcrma), "Hierarchical Clustering Heatmap gcrma Normalization")
heatmap(combatMas5, "Hierarchical Clustering Heatmap mas5 Batch Correction")
heatmap(combatRma, "Hierarchical Clustering Heatmap rma Batch Correction")
heatmap(combatGcrma, "Hierarchical Clustering Heatmap gcrma Batch Correction")


