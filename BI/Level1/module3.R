setwd("/Volumes/SamsungUSB/STEM-Away/BI-mentoring/Level1/SA-BI-Project")

library(affy)
library(affyPLM)
library(simpleaffy)
library(arrayQualityMetrics) 
library(affyQCReport)
library(sva)
library(ggplot2)
library(pheatmap)

# remove file extensions from labels
gse <- readRDS("Data/gse.rds")
# for (i in 1:length(colnames(exprs(gse)))) {
#   if (i < 35) {
#     row.names(gse@phenoData@data)[i] <- strsplit(row.names(gse@phenoData@data)[i],"_")[[1]][1]
#     row.names(gse@protocolData@data)[i] <- strsplit(row.names(gse@protocolData@data)[i],"_")[[1]][1]
#     colnames(exprs(gse))[i] <- strsplit(colnames(exprs(gse))[i], "_")[[1]][1] #RLE/NUSE
#   } else {
#     row.names(gse@phenoData@data)[i] <- strsplit(row.names(gse@phenoData@data)[i],".CEL.gz")[[1]][1]
#     row.names(gse@protocolData@data)[i] <- strsplit(row.names(gse@protocolData@data)[i],".CEL.gz")[[1]][1]
#     colnames(exprs(gse))[i] <- strsplit(colnames(exprs(gse))[i], ".CEL.gz")[[1]][1] #RLE/NUSE
#   }
# }
# 
# saveRDS(gse, "Data.gse.rds")

### QUALITY CONTROL
##### SIMPLEAFFY
sa <- qc(gse)
plot(sa)


##### ARRAYQUALITYMETRICS
# arrayQualityMetrics(gse, "QC/arrayQualityMetrics/", force=T, do.logtransform=T)


##### AFFYQCREPORT
# QCReport(gse, "QC/QCReport.pdf")


##### AFFYPLM
pset <- fitPLM(gse, background=T, normalize=T)
par(mai=c(1.5, 1, 0.5, 0.5))
RLE(pset, main="RLE", ylab="Probe Intensities", las=2)
par(mai=c(1.5, 1, 0.5, 0.5))
NUSE(pset, main="NUSE", ylab="Probe Intensities", las=2)



### BACKGROUND CORRECTION AND NORMALIZATION
##### MAS5
# mas5 <- mas5(gse)
# write.csv(log2(exprs(mas5)), "Data/mas5.csv")
mas5 <- read.csv("Data/mas5.csv", row.names=1)


##### RMA
# rma <- rma(gse)
# write.csv(exprs(rma), "Data/rma.csv")
rma <- read.csv("Data/rma.csv", row.names=1)


##### GCRMA
# gcrma <- gcrma(gse)
# write.csv(exprs(gcrma), "Data/gcrma.csv")
gcrma <- read.csv("Data/gcrma.csv", row.names=1)



### BATCH CORRECTION
# meta <- read.csv("Metadata/metadata.txt", row.names=1)
# design <- model.matrix(~Tissue, meta)
# bc <- ComBat(rma, meta$Batch, design)
# write.csv(bc,"Data/bc.csv")
bc <- read.csv("Data/bc.csv", row.names=1)



### VISUALIZATIONS
##### BOXPLOTS
boxplot(gse, main="Raw", ylab="Probe Intensities", las=2, col=c("red", "orange", "yellow", "green", "blue", "purple"))
boxplot(rma, main="RMA", ylab="Probe Intensities", las=2, col=c("red", "orange", "yellow", "green", "blue", "purple"))
boxplot(bc, main="RMA Batch Corrected", ylab="Probe Intensities", las=2, col=c("red", "orange", "yellow", "green", "blue", "purple"))


##### PCA
group <- as.factor(meta$Group)

####### RAW
pca_raw <- prcomp(exprs(gse), scale=F, center=F)
pca_raw <- as.data.frame(pca_raw$rotation)
ggplot(pca_raw, aes(x=PC1, y=PC2, color=group))+
  geom_point()+ stat_ellipse()+
  scale_colour_manual(values=c("hotpink1", "navy", "goldenrod2", "aquamarine4"))+
  geom_text(aes(label=row.names(meta)),hjust=0, vjust=0)+
  ggtitle("Raw")+ theme_bw()

####### NORM
pca_norm <- prcomp(rma, scale=F, center=F)
pca_norm <- as.data.frame(pca_norm$rotation)
ggplot(pca_norm, aes(x=PC1, y=PC2, color=group))+
  geom_point()+ stat_ellipse()+
  scale_colour_manual(values=c("hotpink1", "navy", "goldenrod2", "aquamarine4"))+
  geom_text(aes(label=row.names(meta)),hjust=0, vjust=0)+
  ggtitle("RMA")+ theme_bw()

####### BC
pca_bc <- prcomp(bc, scale=F, center=F)
pca_bc <- as.data.frame(pca_bc$rotation)
ggplot(pca_bc, aes(x=PC1, y=PC2, color=group))+
  geom_point()+ stat_ellipse()+
  scale_colour_manual(values=c("hotpink1", "navy", "goldenrod2", "aquamarine4"))+
  geom_text(aes(label=row.names(meta)),hjust=0, vjust=0)+
  ggtitle("RMA Batch Corrected")+ theme_bw()


##### HEATMAP
group <- data.frame(Tissue=meta$Group)
row.names(group) <- meta$Sample

####### RAW
dismat_raw <- 1-cor(exprs(gse))
pheatmap(dismat_raw, annotation_col=group, annotation_row=group, main="Raw", labels_col=meta$Tissue, 
         annotation_colors=list(Tissue=c(B1_Normal="hotpink1", B1_Cancer="navy", B2_Normal="goldenrod2", B2_Cancer="aquamarine4")))

####### NORM
dismat_norm <- 1-cor(rma)
pheatmap(dismat_norm, annotation_col=group, annotation_row=group, main="RMA", labels_col=meta$Tissue, 
         annotation_colors=list(Tissue=c(B1_Normal="hotpink1", B1_Cancer="navy", B2_Normal="goldenrod2", B2_Cancer="aquamarine4")))

####### BC
dismat_bc <- 1-cor(bc)
pheatmap(dismat_bc, annotation_col=group, annotation_row=group, main="RMA Batch Corrected", labels_col=meta$Tissue, 
         annotation_colors=list(Tissue=c(B1_Normal="hotpink1", B1_Cancer="navy", B2_Normal="goldenrod2", B2_Cancer="aquamarine4")))
