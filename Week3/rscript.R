setwd("C:/Users/durbe/Documents/repositories/StemAway/Week3")
library(affy)
library(simpleaffy)
library(gcrma)
library(affyPLM)
library(ggplot2)

data <- ReadAffy(compress=TRUE, celfile.path="./data/raw/")

# Quality Control
## simpleAffy
report_simple <- qc(data)
plot(report_simple)

data_norm <- call.exprs(data, "mas5")
report_simple_norm <- qc(data, data_norm)
plot(report_simple_norm)

wd# Data Preprocessing
norm <- gcrma(data, normalize=TRUE)
norm <- as.data.frame(exprs(norm))
colnames(norm) <- sub(".CEL.gz", "", colnames(norm))
write.csv(norm, "./data/norm.csv")

# Visualization
pset_norm <- fitPLM(data, background=FALSE, background.method="GCRMA")
pset_raw <- fitPLM(data, background=FALSE, normalize=FALSE)

## RLE
RLE(pset_norm, main="RLE Boxplot (Normalized)")
rle_norm <- RLE(pset_norm, type="stats")
rle_norm <- as.data.frame(t(rle_norm))
ggplot(rle_norm, aes(median))+
  geom_histogram()+
  xlab("Median")+ ylab("Frequency")+ggtitle("RLE Histogram (Normalized)")+
  theme_bw()
  
RLE(pset_raw, main="RLE Boxplot (Raw)")
rle_raw <- RLE(pset_raw, type="stats")
rle_raw <- as.data.frame(t(rle_raw))
ggplot(rle_raw, aes(median))+
  geom_histogram()+
  xlab("Median")+ ylab("Frequency")+ggtitle("RLE Histogram (Raw)")+
  theme_bw()

## NUSE
NUSE(pset_norm, main="NUSE Boxplot (Normalized)")
nuse_norm <- NUSE(pset_norm, type="stats")
nuse_norm <- as.data.frame(t(nuse_norm))
ggplot(nuse_norm, aes(median))+
  geom_histogram()+
  xlab("Median")+ ylab("Frequency")+ggtitle("NUSE Histogram (Normalized)")+
  theme_bw()

NUSE(pset_raw, main="NUSE Boxplot (Raw)")
nuse_raw <- NUSE(pset_raw, type="stats")
nuse_raw <- as.data.frame(t(nuse_raw))
ggplot(nuse_raw, aes(median))+
  geom_histogram()+
  xlab("Median")+ ylab("Frequency")+ ggtitle("NUSE Histogram (Raw)")+
  theme_bw()

## PCA
pca <- prcomp(norm, center=FALSE, scale=FALSE)
plot(pca, type="l", main="Variation by PC")

pc <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2])
ggplot(pc)+
  geom_point(aes(x=PC1, y=PC2))+
  xlab("PC1")+ ylab("PC2")+ ggtitle("PCA")+
  theme_bw()
