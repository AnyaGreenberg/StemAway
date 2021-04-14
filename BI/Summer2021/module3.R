setwd('/Volumes/SamsungUSB/stemaway/bi-mentoring/level1/2021summer')

library(affy)
library(affyPLM)
library(arrayQualityMetrics) 
library(affyQCReport)
# library(sva)
library(ggplot2)
library(pheatmap)

# remove file extensions from labels
gse <- readRDS('Data/gse.rds')
# for (i in 1:length(colnames(exprs(gse)))) {
#     row.names(gse@phenoData@data)[i] <- strsplit(row.names(gse@phenoData@data)[i],'.CEL')[[1]][1]
#     row.names(gse@protocolData@data)[i] <- strsplit(row.names(gse@protocolData@data)[i], '.CEL')[[1]][1]
#     colnames(exprs(gse))[i] <- strsplit(colnames(exprs(gse))[i], '.CEL')[[1]][1] #RLE/NUSE
# }
# #
# saveRDS(gse, 'Data/gse.rds')

### QUALITY CONTROL
##### ARRAYQUALITYMETRICS
# arrayQualityMetrics(gse, 'QC/array-quality-metrics/', force=T, do.logtransform=T)


##### AFFYQCREPORT
# QCReport(gse, 'QC/qc-report.pdf')


##### AFFYPLM
pset <- fitPLM(gse, background=T, normalize=T)
par(mai=c(1.5, 1, 0.5, 0.5))
RLE(pset, main='RLE', ylab='Probe Intensities', las=2)
par(mai=c(1.5, 1, 0.5, 0.5))
NUSE(pset, main='NUSE', ylab='Probe Intensities', las=2)



### BACKGROUND CORRECTION AND NORMALIZATION
##### MAS5
# mas5 <- mas5(gse)
# write.csv(log2(exprs(mas5)), 'data/mas5.csv')
mas5 <- read.csv('data/mas5.csv', row.names=1)


##### RMA
# rma <- rma(gse)
# write.csv(exprs(rma), 'data/rma.csv')
rma <- read.csv('data/rma.csv', row.names=1)


##### GCRMA
# gcrma <- gcrma(gse)
# write.csv(exprs(gcrma), 'data/gcrma.csv')
gcrma <- read.csv('data/gcrma.csv', row.names=1)

### VISUALIZATIONS
# paper <- rma[,c(1,2,7,16,19:24,26:32,34:36,38,40,43,45:47,55:59,
#                 61,62,67,76,79:84,86:92,94:96,98,100,103,105:107,115:119)]
# write.csv(paper, 'data/selected-rma.csv')
paper <- read.csv('data/selected-rma.csv', row.names=1)

##### BOXPLOTS
boxplot(gse, main='Raw', ylab='Probe Intensities', las=2, col=c('red', 'orange', 'yellow', 'green', 'blue', 'purple'))
boxplot(rma, main='RMA', ylab='Probe Intensities', las=2, col=c('red', 'orange', 'yellow', 'green', 'blue', 'purple'))


##### PCA
meta <- read.csv('metadata/metadata.txt', row.names=1)
paper.meta <- read.csv('metadata/selected-metadata.csv', row.names=1)
# paper.meta <- meta[c(1,2,7,16,19:24,26:32,34:36,38,40,43,45:47,55:59,
#                   61,62,67,76,79:84,86:92,94:96,98,100,103,105:107,115:119),]
# write.csv(paper.meta, 'metadata/selected-metadata.csv')
group <- as.factor(meta$Tissue)

####### RAW
pca_raw <- prcomp(exprs(gse), scale=F, center=F)
pca <- as.data.frame(pca_raw$rotation)
ggplot(pca, aes(x=PC1, y=PC2, color=group))+
  geom_point()+ stat_ellipse()+
  scale_colour_manual(values=c('hotpink1', 'navy'))+
  geom_text(aes(label=row.names(meta)),hjust=0, vjust=0)+
  ggtitle('PCA - Raw')+ theme_bw()

####### NORM
pca_norm <- prcomp(rma, scale=F, center=F)
pca <- as.data.frame(pca_norm$rotation)
ggplot(pca, aes(x=PC1, y=PC2, color=group))+
  geom_point()+ stat_ellipse()+
  scale_colour_manual(values=c('hotpink1', 'navy'))+
  geom_text(aes(label=row.names(meta)),hjust=0, vjust=0)+
  ggtitle('PCA - RMA')+ theme_bw()

####### NORM - SELECTED
group <- as.factor(paper.meta$Tissue)
pca_norm <- prcomp(paper, scale=F, center=F)
pca <- as.data.frame(pca_norm$rotation)
ggplot(pca, aes(x=PC1, y=PC2, color=group))+
  geom_point()+ stat_ellipse()+
  scale_colour_manual(values=c('hotpink1', 'navy'))+
  geom_text(aes(label=row.names(paper.meta)),hjust=0, vjust=0)+
  ggtitle('PCA - RMA (selected)')+ theme_bw()


##### HEATMAP
group <- data.frame(Tissue=meta$Tissue)
row.names(group) <- meta$Sample

####### RAW
dismat_raw <- 1-cor(exprs(gse))
pheatmap(dismat_raw, annotation_col=group, annotation_row=group, main='Heatmap - Raw', labels_col=meta$Tissue, 
         annotation_colors=list(Tissue=c(CANCER='hotpink1', NORMAL='navy')))

####### NORM
dismat_norm <- 1-cor(rma)
pheatmap(dismat_norm, annotation_col=group, annotation_row=group, main='Heatmap - RMA', labels_col=meta$Tissue, 
         annotation_colors=list(Tissue=c(CANCER='hotpink1', NORMAL='navy')))

####### NORM
group <- data.frame(Tissue=paper.meta$Tissue)
row.names(group) <- paper.meta$Sample
dismat_norm <- 1-cor(paper)
pheatmap(dismat_norm, annotation_col=group, annotation_row=group, main='Heatmap - RMA (selected)', labels_col=paper.meta$Tissue, 
         annotation_colors=list(Tissue=c(CANCER='hotpink1', NORMAL='navy')))
