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

library(hgu133plus2.db)
library(WGCNA)
library(limma)
library(EnhancedVolcano)
library(ggplot2)
library(pheatmap)

row.names(filtMas5) <- filtMas5$X
filtMas5 <- filtMas5[-1]

row.names(filtGcrma) <- filtGcrma$X
filtGcrma <- filtGcrma[-1]

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
box_plot(log2(exprs(mas5)), "Boxplot mas5 Data")

rma <- rma(gse)
box_plot(exprs(rma), "Boxplot rma Data")

gcrma <- gcrma(gse)
box_plot(exprs(gcrma), "Boxplot gcrma Data")

### BATCH CORRECTION
pheno <- read.csv("./Data/metadata.csv")
type <- pheno$CN
row.names(pheno) <- pheno$X.Sample_geo_accession
design <- model.matrix(~0+type, pheno)

combatMas5 <- ComBat(log2(exprs(mas5)), pheno$BATCH, design[,1])
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



##### WEEK FOUR #####
mas5 <- read.csv("./Data/BatchCorrected/mas5.csv")
rma <- read.csv("./Data/BatchCorrected/rma.csv")
gcrma <- read.csv("./Data/BatchCorrected/gcrma.csv")

row.names(gcrma) <- gcrma$X
gcrma <- gcrma[-1]

### ANNOTATION
annot <- function(d) {
  symbols <- select(hgu133plus2.db, keys=rownames(d), columns=c("SYMBOL"), keytype="PROBEID")
  symbols <- symbols[!duplicated(symbols$PROBEID),]
  
  row.names(symbols) <- symbols$PROBEID
  d <- merge(symbols, d, by="row.names")
  row.names(d) <- d$Row.names
  d <- d[-c(1,2)]
  
  d <- na.omit(d)
  temp <- collapseRows(d[-1], rowGroup=d$SYMBOL, rowID=rownames(d))
  return(temp)
}

annotMas5 <- annot(mas5)
annotRma <- annot(rma)
annotGcrma <- annot(gcrma)

temp <- annotGcrma$datETcollapsed

### GENE FILTERING
filt <- function(d) {
  means <- rowMeans(d)
  perc2 <- as.numeric(quantile(means, probs=0.02, na.rm=T))
  temp <- d[which(means > perc2),]
  return(temp)
}

filtMas5 <- filt(annotMas5$datETcollapsed)
filtRma <- filt(annotRma$datETcollapsed)
filtGcrma <- filt(annotGcrma$datETcollapsed)


### LIMMA
pheno <- read.csv("./Data/metadata.csv")
row.names(pheno) <- pheno$X.Sample_geo_accession
type <- pheno$CN
names(type) <- pheno$X.Sample_geo_accession
design <- model.matrix(~0+type, pheno)

max <- 50000

limmaFit <- function(d) {
  lm <- lmFit(d, design)
  contrast <- makeContrasts(typeCancer-typeNormal, levels=design)
  colnames(contrast) <- c("Cancer-Control")
  fit <- contrasts.fit(lm, contrast)
  fit <- eBayes(fit)
  return(fit)
}

geneMas5 <- rownames(filtMas5)
geneRma <- rownames(filtRma)
geneGcrma <- rownames(filtGcrma)

fitMas5 <- limmaFit(filtMas5)
fitRma <- limmaFit(filtRma)
fitGcrma <- limmaFit(filtGcrma)

tTMas5 <- topTable(fitMas5, genelist=geneMas5, adjust.method="fdr", sort.by="P", number=max)
tTRma <- topTable(fitRma, genelist=geneRma, adjust.method="fdr", sort.by="P", number=max)
tTGcrma <- topTable(fitGcrma, genelist=geneGcrma, adjust.method="fdr", sort.by="P", number=max)

tTGcrma <- topTable(fitGcrma, coef=2, genelist=geneGcrma, adjust.method="fdr", sort.by="P", number=max)

length(rownames(tTMas5[which(tTMas5$adj.P.Val < 0.05),]))
length(rownames(tTRma[which(tTRma$adj.P.Val < 0.05),]))
length(rownames(tTGcrma[which(tTGcrma$adj.P.Val < 0.05),]))

## volcano plot
EnhancedVolcano(tTMas5, lab=row.names(tTMas5), x="logFC", y="adj.P.Val", 
                pointSize=1, legendLabSize=10, labSize=3.0,
                title="Volcano Plot", subtitle="mas5 Normalization")
EnhancedVolcano(tTRma, lab=row.names(tTRma), x="logFC", y="adj.P.Val", 
                pointSize=1, legendLabSize=10, labSize=3.0,
                title="Volcano Plot", subtitle="rma Normalization")
EnhancedVolcano(tTGcrma, lab=row.names(tTGcrma), x="logFC", y="adj.P.Val", 
                pointSize=1, legendLabSize=10, labSize=3.0,
                title="Volcano Plot", subtitle="gcrma Normalization")

## heatmaps
hm <- function(dFit, dFilt, gl, t) {
  temp <- topTable(dFit, genelist=gl, adjust.method="fdr", sort.by="P", number=50)
  input <- merge(temp, dFilt, by="row.names")
  row.names(input) <- input$Row.names
  input <- input[-c(1,2,3,4,5,6,7,8)]
  test <- strsplit(colnames(input), "_")
  for (i in 1:98) {
    colnames(input)[i] <- test[[i]][1]
  }
  colnames(input) <- sub(".CEL.gz", "", colnames(input))
  par(2,2,2,2)
  pheatmap(input, main=t, cluster_rows=F, legend=T, annotation_col=as.data.frame(type))
}

hm(fitMas5, filtMas5, geneMas5, "Heatmap mas5")
hm(fitRma, filtRma, geneRma, "Heatmap rma")
hm(fitGcrma, filtGcrma, geneGcrma, "Heatmap gcrma")





































