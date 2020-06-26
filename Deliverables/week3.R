setwd("C:/Users/durbe/Documents/repositories/StemAway/Deliverables")
library(affy)
library(simpleaffy)
library(gcrma)
library(affyPLM)
library(ggplot2)
library(pheatmap)

data <- ReadAffy(compress=TRUE, celfile.path="./data/raw/")

# Quality Control
## simpleAffy
report_simple <- qc(data)
plot(report_simple)

data_norm <- call.exprs(data, "mas5")
report_simple_norm <- qc(data, data_norm)
plot(report_simple_norm)


# Data Preprocessing
raw <- as.data.frame(exprs(data))
colnames(raw) <- sub(".CEL.gz", "", colnames(raw))
write.csv(raw, "./data/raw.csv")

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
norm <- read.csv("./data/norm.csv")
row.names(norm) <- norm$X
norm <- norm[-1]

raw <- read.csv("./data/raw.csv")
row.names(raw) <- raw$X
raw <- raw[-1]

cancer <- factor(c(rep("control",32), rep("cancer",32)))

ta <- t(norm) 
pca <- prcomp(ta, center=FALSE, scale=FALSE)
pc <- as.data.frame(pca$x)
ggplot(pc, aes(x=PC1, y=PC2, color=cancer))+
  geom_point()+
  ggtitle("PCA after GCRMA Normalization")+ theme_bw()

ta <- t(raw) 
pca <- prcomp(ta, center=FALSE, scale=FALSE)
pc <- as.data.frame(pca$x)
ggplot(pc, aes(x=PC1, y=PC2, color=cancer))+
  geom_point()+
  ggtitle("PCA before GCRMA Normalization")+ theme_bw()

## heatmap
dismat <- 1-cor(norm)
colnames(dismat) <- factor(c(rep("control",32), rep("cancer",32)))
pheatmap(dismat, main="Hierarchical Clustering Heatmap after GCRMA Normalization")

dismat <- 1-cor(raw)
colnames(dismat) <- factor(c(rep("control",32), rep("cancer",32)))
pheatmap(dismat, main="Hierarchical Clustering Heatmap before GCRMA Normalization")

