setwd("C:/Users/durbe/OneDrive/Documents/ReposExtras/STEM-Away/July2020/final")

##### PHENO DATA #####
library(GEOquery)
gse <- getGEO(filename="./GEODatasets/AML/GSE97485_series_matrix.txt.gz")

# 2. copy phenotype data and transpose into a new sheet
pheno <- gse@phenoData@data

# 3. remove non-relevant columns/rows
meta <- factor(rep(c("cancer", "normal"), c(20,10)))
names(meta) <- row.names(pheno)



##### WEEK TWO - WEEK THREE #####
library(affy)
library(affyPLM)
library(simpleaffy)
library(arrayQualityMetrics)
library(affyQCReport)
library(gcrma)
library(sva)
library(ggplot2)
library(pheatmap)

gse <- ReadAffy(compress=T, celfile.path="./GEODatasets/AML/GSE97485/")

### QUALITY CONTROL - before normalization
arrayQualityMetrics(gse, "./QC/arrayQualityMetrics/", force=T, do.logtransform=T)

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
par(mai=c(3.5,1,1,1))
boxplot(gse, main="Boxplot Raw Data", ylab="Probe Intensities", las=2, col=c("red", "orange", "yellow", "green", "blue", "purple"))

rma <- rma(gse)
saveRDS(rma, "./QC/data/norm/rma.rds")
par(mai=c(3.5,1,1,1))
boxplot(exprs(rma), main="Boxplot rma Data", ylab="Probe Intensities", las=2, col=c("red", "orange", "yellow", "green", "blue", "purple"))

mas5 <- mas5(gse)
gcrma <- gcrma(gse)


#format colnames
temp <- strsplit(colnames(rma), "_")
for (i in 1:98) {
  colnames(rma)[i] <- temp[[i]][1]
}
colnames(rma) <- sub(".CEL.gz", "", colnames(rma))


### VISUALIZATION
## PCA
pca_plot <- function(d, t) {
  group <- as.factor(meta$mygroup)
  pca <- prcomp(d, scale=F, center=F)
  pca <- as.data.frame(pca$rotation)
  ggplot(pca, aes(x=PC1, y=PC2, color=group))+
    geom_point()+ stat_ellipse()+
    scale_colour_manual(values=c("hotpink1", "navy", "goldenrod2", "aquamarine4"))+
    ggtitle(t)+ theme_bw()
}

pca_plot(exprs(gse), "PCA before Normalization")
pca_plot(exprs(rma), "PCA rma Normalization")


## HEATMAP
heatmap <- function(d,t) {
  dismat <- 1-cor(d)
  group <- data.frame(sample=factor(meta$CN))
  row.names(group) <- colnames(d)
  pheatmap(dismat, cluster_rows=F, annotation_col=group, main=t)
}

heatmap(exprs(rma), "Hierarchical Clustering Heatmap rma Normalization")
heatmap(combatRma, "Hierarchical Clustering Heatmap rma Batch Correction")



##### WEEK FOUR #####
library(hgu133plus2.db)
library(WGCNA)
library(limma)
library(EnhancedVolcano)
library(ggplot2)
library(pheatmap)


rma <- read.csv("./QC/data/bc/rma.csv")

row.names(rma) <- rma$X
rma <- rma[-1]


### ANNOTATION
annot <- function(d) {
  symbols <- AnnotationDbi::select(hgu133plus2.db, keys=rownames(d), columns=c("SYMBOL"), keytype="PROBEID")
  symbols <- symbols[!duplicated(symbols$PROBEID),]
  row.names(symbols) <- symbols$PROBEID
  d <- merge(symbols, d, by="row.names")
  row.names(d) <- d$Row.names
  d <- d[-c(1,2)]
  d <- na.omit(d)
  d <- collapseRows(d[-1], rowGroup=d$SYMBOL, rowID=rownames(d))
  return(d)
}

annotRma <- annot(rma)


### GENE FILTERING
filt <- function(d) {
  means <- rowMeans(d)
  perc2 <- as.numeric(quantile(means, probs=0.02, na.rm=T))
  temp <- d[which(means > perc2),]
  return(temp)
}

filtRma <- filt(annotRma$datETcollapsed)


### LIMMA
type <- factor(meta$CN, levels = c('Normal','Cancer'), ordered = F)
row.names(meta) <- meta$X.Sample_geo_accession
design <- model.matrix(~0+type, meta)

max <- 50000

limmaFit <- function(d) {
  lm <- lmFit(d, design)
  contrast <- makeContrasts(typeCancer-typeNormal, levels=design)
  colnames(contrast) <- c("Cancer-Control")
  fit <- contrasts.fit(lm, contrast)
  fit <- eBayes(fit)
  return(fit)
}

geneRma <- rownames(filtRma)
fitRma <- limmaFit(filtRma)

tTRma <- topTable(fitRma, coef=1, p.value=0.05, sort="P", genelist=geneRma, number=max)
write.table(tTRma, "./DGE/data/rma.txt")

length(rownames(tTRma[which(tTRma$adj.P.Val < 0.05),]))


## VOLCANO PLOT
EnhancedVolcano(tTRma, lab=row.names(tTRma), x="logFC", y="adj.P.Val", 
                pointSize=1, legendLabSize=10, labSize=3.0,
                title="Volcano Plot", subtitle="rma Normalization")


## HEATMAP
hm <- function(dFit, dFilt, gl, t) {
  t50 <- topTable(dFit, genelist=gl, adjust.method="fdr", sort.by="P", number=50)
  input <- dFilt[(rownames(dFilt) %in% rownames(t50)),]
  group <- data.frame(sample=factor(meta$CN))
  row.names(group) <- colnames(dFilt)
  par(2,2,2,2)
  pheatmap(input, main=t, cluster_rows=F, annotation_col=group)
}

hm(fitRma, filtRma, geneRma, "Heatmap rma")


##### WEEK FIVE #####library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(magrittr)
library(tidyr)
library(msigdbr)
library(enrichplot)


rma <- read.table("./DGE/data/rma.txt")

### FUNCTIONAL ANALYSIS
## DEG VECTOR
genelist <- function(d) {
  gene <- d$logFC
  names(gene) <- d$ID
  DEG.gene <- gene[gene > 1.5]
  DEG.gene <- sort(DEG.gene, decreasing=T)
  return(DEG.gene)
}

rGene <- genelist(rma)
DEG.r <- AnnotationDbi::select(org.Hs.eg.db, keys=names(rGene), columns=c("ENTREZID"), keytype="SYMBOL")


# GENE ONTOLOGY
eGoD <- function(d, o, f) {
  ego <- enrichGO(d$ENTREZID, org.Hs.eg.db, ont=o, readable=T)
  write.csv(ego, f)
}

eGoD(DEG.r, "CC", "./FA/data/GO/rmaCC.csv")
eGoD(DEG.r, "BP", "./FA/data/GO/rmaBP.csv")
eGoD(DEG.r, "MF", "./FA/data/GO/rmaMF.csv")

egoPlot  <- function(d, o, t) {
  ego <- enrichGO(d$ENTREZID, org.Hs.eg.db, ont=o, readable=F)
  egoR <- setReadable(ego, org.Hs.eg.db, keyType="ENTREZID")
  barplot(egoR, title=t)
}

egoPlot(DEG.r, "CC", "rma CC")
egoPlot(DEG.r, "BP", "rma BP")
egoPlot(DEG.r, "MF", "rma MF")

# OPTIONAL***
ggoPlot <- function(d, o, l, t) {
  ggo <- groupGO(d$ENTREZID, org.Hs.eg.db, keyType="ENTREZID", ont=o, level=l, readable=F)
  barplot(ggo, title=t)
}


### KEGG ANALYSIS
keggPlot <- function(d, gene, c, t) {
  ek <- enrichKEGG(d$ENTREZID)
  dotplot(ek, title=t)
}

keggPlot(DEG.r, rGene, T, "rma KEGG")


### GENE CONCEPT NETWORK
cnetPlot <- function(d, gene, c) {
  # ed <- enrichDGN(d$ENTREZID)
  ek <- enrichKEGG(d$ENTREZID)
  edR <- setReadable(ek, org.Hs.eg.db, keyType="ENTREZID")
  cnetplot(edR, foldChange=gene, categorySize="pvalue", colorEdge=T, circular=c)
}

cnetPlot(DEG.r, rGene, F)
cnetPlot(DEG.r, rGene, T)


## TRANSCRIPTIONAL FACTOR ANALYSIS
tfa <- function(d, gene) {
  c3 <- read.gmt("./FA/data/c3.gmt")
  e <- enricher(DEG.m$ENTREZID, TERM2GENE=c3)
  temp <- setReadable(e, org.Hs.eg.db, keyType="ENTREZID")
  cnetplot(temp, foldChange=gene, categorySize="pvalue", colorEdge=T)
}

tfa(DEG.r, rGene)


# EXTERNAL TOOLS

write(names(rGene), "./FA/data/genelist/rma_up.txt")
write(names(rGeneA), "./FA/data/genelist/rma_all.txt")
