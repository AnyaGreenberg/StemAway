setwd("C:/Users/durbe/OneDrive/Documents/ReposExtras/STEM-Away/July2020/")

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


gse32323 <- ReadAffy(compress=T, celfile.path="./data/raw/GSE32323/")
gse8671 <- ReadAffy(compress=T, celfile.path="./data/raw/GSE8671/") 
gse <- merge(gse8671, gse32323)
saveRDS(gse, "./data/gse.rds")


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
par(mai=c(3.5,1,1,1))
boxplot(gse, main="Boxplot Raw Data", ylab="Probe Intensities", las=2, col=c("red", "orange", "yellow", "green", "blue", "purple"))

mas5 <- mas5(gse)
saveRDS(mas5, "./QC/data/norm/mas5.rds")
par(mai=c(3.5,1,1,1))
boxplot(log2(exprs(mas5)), main="Boxplot mas5 Data", ylab="Probe Intensities", las=2, col=c("red", "orange", "yellow", "green", "blue", "purple"))

rma <- rma(gse)
saveRDS(rma, "./QC/data/norm/rma.rds")
par(mai=c(3.5,1,1,1))
boxplot(exprs(rma), main="Boxplot rma Data", ylab="Probe Intensities", las=2, col=c("red", "orange", "yellow", "green", "blue", "purple"))

gcrma <- gcrma(gse)
saveRDS(gcrma, "./QC/data/norm/gcrma.rds")
par(mai=c(3.5,1,1,1))
boxplot(exprs(gcrma), main="Boxplot gcrma Data", ylab="Probe Intensities", las=2, col=c("red", "orange", "yellow", "green", "blue", "purple"))


#format colnames
temp <- strsplit(colnames(gcrma), "_")
for (i in 1:98) {
  colnames(gcrma)[i] <- temp[[i]][1]
}
colnames(gcrma) <- sub(".CEL.gz", "", colnames(gcrma))


### BATCH CORRECTION
meta <- read.csv("./Data/metadata.csv")
design <- model.matrix(~CN, meta)

combatMas5 <- ComBat(log2(exprs(mas5)), meta$BATCH, design)
write.csv(combatMas5,"./QC/data/bc/mas5.csv")

combatRma <- ComBat(exprs(rma), meta$BATCH, design)
write.csv(combatRma,"./QC/data/bc/rma.csv")

combatGcrma <- ComBat(exprs(gcrma), meta$BATCH, design)
write.csv(combatGcrma,"./QC/data/bc/gcrma.csv")

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

pca_plot(exprs(mas5), "PCA mas5 Normalization")
pca_plot(exprs(rma), "PCA rma Normalization")
pca_plot(exprs(gcrma), "PCA gcrma Normalization")

pca_plot(combatMas5, "PCA mas5 Batch Correction")
pca_plot(combatRma, "PCA rma Batch Correction")
pca_plot(combatGcrma, "PCA gcrma Batch Correction")


## HEATMAP
heatmap <- function(d,t) {
  dismat <- 1-cor(d)
  group <- data.frame(sample=factor(meta$CN))
  row.names(group) <- colnames(d)
  pheatmap(dismat, cluster_rows=F, annotation_col=group, main=t)
}

heatmap(exprs(mas5), "Hierarchical Clustering Heatmap mas5 Normalization")
heatmap(exprs(rma), "Hierarchical Clustering Heatmap rma Normalization")
heatmap(exprs(gcrma), "Hierarchical Clustering Heatmap gcrma Normalization")

heatmap(combatMas5, "Hierarchical Clustering Heatmap mas5 Batch Correction")
heatmap(combatRma, "Hierarchical Clustering Heatmap rma Batch Correction")
heatmap(combatGcrma, "Hierarchical Clustering Heatmap gcrma Batch Correction")



##### WEEK FOUR #####
library(hgu133plus2.db)
library(WGCNA)
library(limma)
library(EnhancedVolcano)
library(ggplot2)
library(pheatmap)


mas5 <- read.csv("./QC/data/bc/mas5.csv")
rma <- read.csv("./QC/data/bc/rma.csv")
gcrma <- read.csv("./QC/data/bc/gcrma.csv")

row.names(mas5) <- mas5$X
mas5 <- mas5[-1]

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

annotMas5 <- annot(mas5)
annotRma <- annot(rma)
annotGcrma <- annot(gcrma)

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
meta <- read.csv("./Data/metadata.csv")
type <- meta$CN
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

geneMas5 <- rownames(filtMas5)
geneRma <- rownames(filtRma)
geneGcrma <- rownames(filtGcrma)

fitMas5 <- limmaFit(filtMas5)
fitRma <- limmaFit(filtRma)
fitGcrma <- limmaFit(filtGcrma)

tTMas5 <- topTable(fitMas5, coef=1, p.value=0.05, sort="P", genelist=geneMas5, number=max)
write.table(tTMas5, "./DGE/data/mas5.txt")
tTRma <- topTable(fitRma, coef=1, p.value=0.05, sort="P", genelist=geneRma, number=max)
write.table(tTRma, "./DGE/data/rma.txt")
tTGcrma <- topTable(fitGcrma, coef=1, p.value=0.05, sort="P", genelist=geneGcrma, number=max)
write.table(tTGcrma, "./DGE/data/gcrma.txt")

length(rownames(tTMas5[which(tTMas5$adj.P.Val < 0.05),]))
length(rownames(tTRma[which(tTRma$adj.P.Val < 0.05),]))
length(rownames(tTGcrma[which(tTGcrma$adj.P.Val < 0.05),]))

## VOLCANO PLOT
EnhancedVolcano(tTMas5, lab=row.names(tTMas5), x="logFC", y="adj.P.Val", 
                pointSize=1, legendLabSize=10, labSize=3.0,
                title="Volcano Plot", subtitle="mas5 Normalization")
EnhancedVolcano(tTRma, lab=row.names(tTRma), x="logFC", y="adj.P.Val", 
                pointSize=1, legendLabSize=10, labSize=3.0,
                title="Volcano Plot", subtitle="rma Normalization")
EnhancedVolcano(tTGcrma, lab=row.names(tTGcrma), x="logFC", y="adj.P.Val", 
                pointSize=1, legendLabSize=10, labSize=3.0,
                title="Volcano Plot", subtitle="gcrma Normalization")

## HEATMAP
hm <- function(dFit, dFilt, gl, t) {
  t50 <- topTable(dFit, genelist=gl, adjust.method="fdr", sort.by="P", number=50)
  input <- dFilt[(rownames(dFilt) %in% rownames(t50)),]
  group <- data.frame(sample=factor(meta$CN))
  row.names(group) <- colnames(dFilt)
  par(2,2,2,2)
  pheatmap(input, main=t, cluster_rows=F, annotation_col=group)
}

hm(fitMas5, filtMas5, geneMas5, "Heatmap mas5")
hm(fitRma, filtRma, geneRma, "Heatmap rma")
hm(fitGcrma, filtGcrma, geneGcrma, "Heatmap gcrma")



##### WEEK FIVE #####library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(magrittr)
library(tidyr)
library(msigdbr)
library(enrichplot)


mas5 <- read.table("./DGE/data/mas5.txt")
rma <- read.table("./DGE/data/rma.txt")
gcrma <- read.table("./DGE/data/gcrma.txt")

### FUNCTIONAL ANALYSIS
## DEG VECTOR
genelist <- function(d) {
  gene <- d$logFC
  names(gene) <- d$ID
  DEG.gene <- gene[gene > 1.5]
  DEG.gene <- sort(DEG.gene, decreasing=T)
  return(DEG.gene)
}

mGene <- genelist(mas5)
rGene <- genelist(rma)
gGene <- genelist(gcrma)

DEG.m <- AnnotationDbi::select(org.Hs.eg.db, keys=names(mGene), columns=c("ENTREZID"), keytype="SYMBOL")
DEG.r <- AnnotationDbi::select(org.Hs.eg.db, keys=names(rGene), columns=c("ENTREZID"), keytype="SYMBOL")
DEG.g <- AnnotationDbi::select(org.Hs.eg.db, keys=names(gGene), columns=c("ENTREZID"), keytype="SYMBOL")


## GENE ONTOLOGY
# eGoD <- function(d, o, f) {
#   ego <- enrichGO(d$ENTREZID, org.Hs.eg.db, ont=o, readable=T)
#   write.csv(ego, f)
# }
# 
# eGoD(DEG.m, "CC", "./FA/data/GO/mas5CC.csv")
# eGoD(DEG.m, "BP", "./FA/data/GO/mas5BP.csv")
# eGoD(DEG.m, "MF", "./FA/data/GO/mas5MF.csv")
# 
# eGoD(DEG.r, "CC", "./FA/data/GO/rmaCC.csv")
# eGoD(DEG.r, "BP", "./FA/data/GO/rmaBP.csv")
# eGoD(DEG.r, "MF", "./FA/data/GO/rmaMF.csv")
# 
# eGoD(DEG.g, "CC", "./FA/data/GO/gcrmaCC.csv")
# eGoD(DEG.g, "BP", "./FA/data/GO/gcrmaBP.csv")
# eGoD(DEG.g, "MF", "./FA/data/GO/gcrmaMF.csv")

egoPlot  <- function(d, o, t) {
  ego <- enrichGO(d$ENTREZID, org.Hs.eg.db, ont=o, readable=F)
  egoR <- setReadable(ego, org.Hs.eg.db, keyType="ENTREZID")
  barplot(egoR, title=t)
}

egoPlot(DEG.m, "CC", "mas5 CC")
egoPlot(DEG.m, "BP", "mas5 BP")
egoPlot(DEG.m, "MF", "mas5 MF")

egoPlot(DEG.r, "CC", "rma CC")
egoPlot(DEG.r, "BP", "rma BP")
egoPlot(DEG.r, "MF", "rma MF")

egoPlot(DEG.g, "CC", "gcrma CC")
egoPlot(DEG.g, "BP", "gcrma BP")
egoPlot(DEG.g, "MF", "gcrma MF")

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

keggPlot(DEG.m, mGene, T, "mas5 KEGG")
keggPlot(DEG.r, rGene, T, "rma KEGG")
keggPlot(DEG.g, gGene, T, "gcrma KEGG")


### GENE CONCEPT NETWORK
cnetPlot <- function(d, gene, c) {
  # ed <- enrichDGN(d$ENTREZID)
  ek <- enrichKEGG(d$ENTREZID)
  edR <- setReadable(ek, org.Hs.eg.db, keyType="ENTREZID")
  cnetplot(edR, foldChange=gene, categorySize="pvalue", colorEdge=T, circular=c)
}

cnetPlot(DEG.m, mGene, F)
cnetPlot(DEG.m, mGene, T)

cnetPlot(DEG.r, rGene, F)
cnetPlot(DEG.r, rGene, T)

cnetPlot(DEG.g, gGene, F)
cnetPlot(DEG.g, gGene, T)


### GSEA
gsealist <- function(d) {
  gene <- d$logFC
  names(gene) <- d$ID
  gene <- sort(gene, decreasing=T)
  return(gene)
}

mGeneA <- gsealist(mas5)
rGeneA <- gsealist(rma)
gGeneA <- gsealist(gcrma)


gsea <- function(gene, d) {
  # h <- read.gmt("./FA/data/hallmark.gmt")
  DEG.entrez <- gene
  names(DEG.entrez) <- d$ENTREZID
  # gsea <- GSEA(DEG.entrez, TERM2GENE=h)
  
  msig <- msigdbr(species="Homo sapiens", category="H")
  msigS <- msig %>% select(gs_name, entrez_gene)
  gsea <- GSEA(DEG.entrez, TERM2GENE=msigS)
  
  return(gsea)
}

gseaM <- gsea(mGeneA, DEG.m)
gseaR <- gsea(rGeneA, DEG.r)
gseaG <- gsea(gGeneA, DEG.g)

enrich <- function(d) {
  h <- read.gmt("./FA/data/hallmark.gmt")
  e <- enricher(d$ENTREZID, TERM2GENE=h)
  eR <- setReadable(e, org.Hs.eg.db, keyType="ENTREZID")
}

eM <- enrich(DEG.m)
eR <- enrich(DEG.r)
eG <- enrich(DEG.g)

#weird plots
gseaplot2(gseaM, geneSetID=2, pvalue_table=T)
gseaplot2(gseaR, geneSetID=1, pvalue_table=T)
gseaplot2(gseaG, geneSetID=5, pvalue_table=T)

heatplot(eM, foldChange=mGeneA)
heatplot(eR, foldChange=rGeneA)
heatplot(eG, foldChange=gGeneA)

## TRANSCRIPTIONAL FACTOR ANALYSIS
tfa <- function(d, gene) {
  c3 <- read.gmt("./FA/data/c3.gmt")
  e <- enricher(DEG.m$ENTREZID, TERM2GENE=c3)
  temp <- setReadable(e, org.Hs.eg.db, keyType="ENTREZID")
  cnetplot(temp, foldChange=gene, categorySize="pvalue", colorEdge=T)
}

tfa(DEG.m, mGene)
tfa(DEG.r, rGene) #error - dataframe should contain at least two columns 
tfa(DEG.g, gGene)

## EXTERNAL TOOLS
# write(names(mGene), "./FA/data/genelist/mas5_up.txt")
# write(names(mGeneA), "./FA/data/genelist/mas5_all.txt")
# 
# write(names(rGene), "./FA/data/genelist/rma_up.txt")
# write(names(rGeneA), "./FA/data/genelist/rma_all.txt")
# 
# write(names(gGene), "./FA/data/genelist/gcrma_up.txt")
# write(names(gGeneA), "./FA/data/genelist/gcrma_all.txt")
