setwd("C:/Users/durbe/OneDrive/Documents/ReposExtras/STEM-Away/June2020")
library(affy)
library(affyPLM)
library(sva) # 'genefilter' package cannot be loaded
library(AnnotationDbi)
library(hgu133plus2.db)
library(simpleaffy)
library(arrayQualityMetrics)
library(affyQCReport)
library(gcrma)
library(ggplot2)
library(pheatmap)

library(limma)
library(EnhancedVolcano)

library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(magrittr)
library(tidyr)

library(GEOquery)

data <- ReadAffy(compress=TRUE, celfile.path="./Data/Raw/")


##### QUALITY CONTROL #####
### assessing quality of RNA samples (RNA degradation, hybridization, spike-in)
## SIMPLEAFFY - acquire analysis report 
report_sa <- qc(data)
plot(report_sa)
## ARRAYQUALITYMETRICS - find correlation between two samples
arrayQualityMetrics(data, "./QC/arrayQyalityMetrics-Raw/", force=TRUE, do.logtransform=TRUE)
# !do this one again with normalized data
## AFFYQCREPORT
QCReport(data, "./QC/affyQCReport.pdf")
## AFFYPLM - compute RLSE and NUSE scores of microarray samples
pset <- fitPLM(data, background=TRUE, normalize=TRUE)
rle <- RLE(pset, type="stats")
rle <- data.frame(Median=rle[1,])
ggplot(rle, aes(Median))+ geom_histogram()+ xlab("Median")+ ylab("Frequency")+ ggtitle("RLE Histogram (Raw)")+ theme_bw()
nuse <- NUSE(pset, type="stats")
nuse <- data.frame(Median=nuse[1,])
ggplot(nuse, aes(Median))+ geom_histogram()+ xlab("Median")+ ylab("Frequency")+ ggtitle("NUSE Histogram (Raw)")+ theme_bw()


##### DATA PREPROCESSING: BACKGROUND CORRECTION + NORMALIZATION #####
## MAS5
mas5N <- mas5(data)
arrayQualityMetrics(mas5N, "./QC/arrayQyalityMetrics-mas5/", force=TRUE, do.logtransform=TRUE)
mas5Ndf <- as.data.frame(exprs(mas5N))
colnames(mas5Ndf) <- sub(".CEL.gz", "", colnames(mas5Ndf))
write.csv(mas5Ndf, "./Data/Normalized/mas5.csv")
## RMA
rmaN <- rma(data)
arrayQualityMetrics(rmaN, "./QC/arrayQyalityMetrics-rma/", force=TRUE, do.logtransform=TRUE)
rmaNdf <- as.data.frame(exprs(rmaN))
colnames(rmaNdf) <- sub(".CEL.gz", "", colnames(rmaNdf))
write.csv(rmaNdf, "./Data/Normalized/rma.csv")
## GCRMA
gcrmaN <- gcrma(data)
arrayQualityMetrics(gcrmaN, "./QC/arrayQyalityMetrics-gcrma/", force=TRUE, do.logtransform=TRUE)
gcrmaNdf <- as.data.frame(exprs(gcrmaN))
colnames(gcrmaNdf) <- sub(".CEL.gz", "", colnames(gcrmaNdf))
write.csv(gcrmaNdf, "./Data/Normalized/gcrma.csv")


##### DATA VISUALIZATION #####
type <- factor(c(rep("control",32), rep("cancer",32)))
## RAW
raw <- as.data.frame(exprs(data))
t <- t(raw)
pca <- prcomp(t, scale=FALSE, center=FALSE)
pca <- as.data.frame(pca$x)
ggplot(pca, aes(x=PC1, y=PC2, color=type))+
  geom_point()+
  ggtitle("PCA before Normalization")+ theme_bw()
dismat <- 1-cor(raw)
colnames(dismat) <- factor(c(rep("control",32), rep("cancer",32)))
pheatmap(dismat, main="Hierarchical Clustering Heatmap before Normalization")
## MAS5
t <- t(mas5Ndf)
pca <- prcomp(t, scale=FALSE, center=FALSE)
pca <- as.data.frame(pca$x)
ggplot(pca, aes(x=PC1, y=PC2, color=type))+
  geom_point()+
  ggtitle("PCA after mas5 Normalization")+ theme_bw()
dismat <- 1-cor(mas5Ndf)
colnames(dismat) <- factor(c(rep("control",32), rep("cancer",32)))
pheatmap(dismat, main="Hierarchical Clustering Heatmap after mas5 Normalization")
## RMA
t <- t(rmaNdf)
pca <- prcomp(t, scale=FALSE, center=FALSE)
pca <- as.data.frame(pca$x)
ggplot(pca, aes(x=PC1, y=PC2, color=type))+
  geom_point()+
  ggtitle("PCA after rma Normalization")+ theme_bw()
dismat <- 1-cor(rmaNdf)
colnames(dismat) <- factor(c(rep("control",32), rep("cancer",32)))
pheatmap(dismat, main="Hierarchical Clustering Heatmap after rma Normalization")
## GCRMA
t <- t(gcrmaNdf)
pca <- prcomp(t, scale=FALSE, center=FALSE)
pca <- as.data.frame(pca$x)
ggplot(pca, aes(x=PC1, y=PC2, color=type))+
  geom_point()+
  ggtitle("PCA after gcrma Normalization")+ theme_bw()
dismat <- 1-cor(gcrmaNdf)
colnames(dismat) <- factor(c(rep("control",32), rep("cancer",32)))
pheatmap(dismat, main="Hierarchical Clustering Heatmap after gcrma Normalization")


##### ANNOTATION #####
mas5Ndf <- mas5Ndf[-c(13,43)]
rmaNdf <- rmaNdf[-c(13,43)]
gcrmaNdf <- gcrmaNdf[-c(13,43)]
## MAS5
mas5P <- select(hgu133plus2.db, keys=rownames(mas5Ndf), columns=c("SYMBOL"), keytype="PROBEID")
mas5P <- mas5P[!duplicated(mas5P$PROBEID),]
mas5P <- mas5P[!duplicated(mas5P$SYMBOL),]
rownames(mas5P) <- mas5P$PROBEID
mas5S <- merge(mas5Ndf, mas5P["SYMBOL"], by="row.names")
mas5S <- na.omit(mas5S)
## RMA
rmaP <- select(hgu133plus2.db, keys=rownames(rmaNdf), columns=c("SYMBOL"), keytype="PROBEID")
rmaP <- rmaP[!duplicated(rmaP$PROBEID),]
rmaP <- rmaP[!duplicated(rmaP$SYMBOL),]
rownames(rmaP) <- rmaP$PROBEID
rmaS <- merge(rmaNdf, rmaP["SYMBOL"], by="row.names")
rmaS <- na.omit(rmaS)
## GCRMA
gcrmaP <- select(hgu133plus2.db, keys=rownames(gcrmaNdf), columns=c("SYMBOL"), keytype="PROBEID")
gcrmaP <- gcrmaP[!duplicated(gcrmaP$PROBEID),]
gcrmaP <- gcrmaP[!duplicated(gcrmaP$SYMBOL),]
rownames(gcrmaP) <- gcrmaP$PROBEID
gcrmaS <- merge(gcrmaNdf, gcrmaP["SYMBOL"], by="row.names")
gcrmaS <- na.omit(gcrmaS)


##### GENE FILTERING #####
## MAS5
means <- rowMeans(mas5S[,2:63])
perc4 <- as.numeric(quantile(means, probs=0.04, na.rm=TRUE))
mas5F <- mas5S[which(means > perc4),]
colnames(mas5F)[1] <- "PROBEID"
write.csv(mas5F, "./Data/Filtered/mas5F.csv")
## RMA
means <- rowMeans(rmaS[,2:63])
perc4 <- as.numeric(quantile(means, probs=0.04, na.rm=TRUE))
rmaF <- rmaS[which(means > perc4),]
colnames(rmaF)[1] <- "PROBEID"
write.csv(rmaF, "./Data/Filtered/rmaF.csv")
## GCRMA
means <- rowMeans(gcrmaS[,2:63])
perc4 <- as.numeric(quantile(means, probs=0.04, na.rm=TRUE))
gcrmaF <- gcrmaS[which(means > perc4),]
colnames(gcrmaF)[1] <- "PROBEID"
write.csv(gcrmaF, "./Data/Filtered/gcrmaF.csv")


##### LIMMA #####
type <- factor(paste(rep(c("Control", "Cancer"), each=31)))
design <- model.matrix(~0+type)
max <- length(mas5Ndf$GSM215051)
## MAS5
genelist <- mas5F$SYMBOL
mas5Fnum <- mas5F[-c(1,64)]
### compute pooled variance
row.names(design) <- colnames(mas5Fnum)
write.csv(design, "./Data/DifferentialExpression/model_matrix.csv")
### compute log2 values
mas5Flog <- log2(mas5Fnum)
### analysis
lm <- lmFit(mas5Flog, design)
contrast <- makeContrasts(typeCancer-typeControl, levels=design)
colnames(contrast) <- c("Cancer-Control")
fit <- contrasts.fit(lm, contrast)
fit <- eBayes(fit)
### compute top 100 gene list
tT100 <- topTable(fit, genelist=genelist, adjust.method="fdr", sort.by="B", number=100)
write.table(tT100, file="./Data/DifferentialExpression/mas5100.txt", sep="\t")
### volcano plot with all genes
tTv <- topTable(fit, genelist=genelist, adjust.method="fdr", sort.by="B", number=max)
write.table(tTv, file="./Data/DifferentialExpression/mas5.txt", sep="\t")
EnhancedVolcano(tTv, lab=tTv$ID, x="logFC", y="adj.P.Val", 
                pCutoff=0.01, FCcutoff=1.0, pointSize=1, legendLabSize=10, labSize=3.0,
                title="Volcano Plot", subtitle="mas5 Normalization, log transform")
length(rownames(tTv[which(tTv$adj.P.Val < 0.01),]))
length(rownames(tTv[which(tTv$adj.P.Val < 0.05),]))
## RMA
genelist <- rmaF$SYMBOL
rmaFnum <- rmaF[-c(1,64)]
row.names(design) <- colnames(rmaFnum)
rmaFlog <- log2(rmaFnum)
lm <- lmFit(rmaFnum, design)
contrast <- makeContrasts(typeCancer-typeControl, levels=design)
colnames(contrast) <- c("Cancer-Control")
fit <- contrasts.fit(lm, contrast)
fit <- eBayes(fit)
tT100 <- topTable(fit, genelist=genelist, adjust.method="fdr", sort.by="B", number=100)
write.table(tT100, file="./Data/DifferentialExpression/rma100.txt", sep="\t")
tTv <- topTable(fit, genelist=genelist, adjust.method="fdr", sort.by="B", number=max)
write.table(tTv, file="./Data/DifferentialExpression/rma.txt", sep="\t")
EnhancedVolcano(tTv, lab=tTv$ID, x="logFC", y="adj.P.Val", 
                pCutoff=0.01, FCcutoff=1.0, pointSize=1, legendLabSize=10, labSize=3.0,
                title="Volcano Plot", subtitle="rma Normalization, no log transform")
length(rownames(tTv[which(tTv$adj.P.Val < 0.01),]))
length(rownames(tTv[which(tTv$adj.P.Val < 0.05),]))
## GCRMA
genelist <- gcrmaF$SYMBOL
gcrmaFnum <- gcrmaF[-c(1,64)]
row.names(design) <- colnames(gcrmaFnum)
gcrmaFlog <- log2(gcrmaFnum)
lm <- lmFit(gcrmaFnum, design)
contrast <- makeContrasts(typeCancer-typeControl, levels=design)
colnames(contrast) <- c("Cancer-Control")
fit <- contrasts.fit(lm, contrast)
fit <- eBayes(fit)
tT100 <- topTable(fit, genelist=genelist, adjust.method="fdr", sort.by="B", number=100)
write.table(tT100, file="./Data/DifferentialExpression/gcrma100.txt", sep="\t")
tTv <- topTable(fit, genelist=genelist, adjust.method="fdr", sort.by="B", number=max)
write.table(tTv, file="./Data/DifferentialExpression/gcrma.txt", sep="\t")
EnhancedVolcano(tTv, lab=tTv$ID, x="logFC", y="adj.P.Val", 
                pCutoff=0.01, FCcutoff=1.0, pointSize=1, legendLabSize=10, labSize=3.0,
                title="Volcano Plot", subtitle="gcrma Normalization, log transform")
length(rownames(tTv[which(tTv$adj.P.Val < 0.01),]))
length(rownames(tTv[which(tTv$adj.P.Val < 0.05),]))


##### GENE VECTOR #####
## MAS5
mas5 <- read.table("./Data/DifferentialExpression/mas5.txt", header=TRUE, sep="\t")
mas5G <- mas5$logFC
names(mas5G) <- mas5$ID
mas5U <- sort(mas5G[mas5G > 2], decreasing=TRUE)
mas5D <- sort(mas5G[mas5G < -2], decreasing=TRUE)
mas5idU <- bitr(names(mas5U), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db, drop=TRUE)
mas5idD <- bitr(names(mas5D), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db, drop=TRUE)
row.names(mas5idU) <- mas5idU$SYMBOL
row.names(mas5idD) <- mas5idD$SYMBOL
## RMA
rma <- read.table("./Data/DifferentialExpression/rma.txt", header=TRUE, sep="\t")
rmaG <- rma$logFC
names(rmaG) <- rma$ID
rmaU <- sort(rmaG[rmaG > 2], decreasing=TRUE)
rmaD <- sort(rmaG[rmaG < -2], decreasing=TRUE)
rmaidU <- bitr(names(rmaU), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db, drop=TRUE)
rmaidD <- bitr(names(rmaD), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db, drop=TRUE)
row.names(rmaidU) <- rmaidU$SYMBOL
row.names(rmaidD) <- rmaidD$SYMBOL
## GCRMA
gcrma <- read.table("./Data/DifferentialExpression/gcrma.txt", header=TRUE, sep="\t")
gcrmaG <- gcrma$logFC
names(gcrmaG) <- gcrma$ID
gcrmaU <- sort(gcrmaG[gcrmaG > 2], decreasing=TRUE)
gcrmaD <- sort(gcrmaG[gcrmaG < -2], decreasing=TRUE)
gcrmaidU <- bitr(names(gcrmaU), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db, drop=TRUE)
gcrmaidD <- bitr(names(gcrmaD), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db, drop=TRUE)
row.names(gcrmaidU) <- gcrmaidU$SYMBOL
row.names(gcrmaidD) <- gcrmaidD$SYMBOL


##### GENE ONTOLOGY ANALYSIS #####
## groupGO 
i <- 5
ont <- "CC"
### MAS5
ggoU <- groupGO(mas5idU$ENTREZID, org.Hs.eg.db, readable=FALSE, level=i, ont=ont)
ggoD <- groupGO(mas5idD$ENTREZID, org.Hs.eg.db, readable=FALSE, level=i, ont=ont)
barplot(ggoU, title=paste("mas5 Normalization (Up-Regulated), Level=", i, sep=""))
barplot(ggoD, title=paste("mas5 Normalization (Down-Regulated), Level=", i, sep=""))
### RMA
ggoU <- groupGO(rmaidU$ENTREZID, org.Hs.eg.db, readable=FALSE, level=i, ont=ont)
ggoD <- groupGO(rmaidD$ENTREZID, org.Hs.eg.db, readable=FALSE, level=i, ont=ont)
barplot(ggoU, title=paste("rma Normalization (Up-Regulated), Level=", i, sep=""))
barplot(ggoD, title=paste("rma Normalization (Down-Regulated), Level=", i, sep=""))
### GCRMA
ggoU <- groupGO(gcrmaidU$ENTREZID, org.Hs.eg.db, readable=FALSE, level=i, ont=ont)
ggoD <- groupGO(gcrmaidD$ENTREZID, org.Hs.eg.db, readable=FALSE, level=i, ont=ont)
barplot(ggoU, title=paste("gcrma Normalization (Up-Regulated), Level=", i, sep=""))
barplot(ggoD, title=paste("gcrma Normalization (Down-Regulated), Level=", i, sep=""))
## enrichGO
ont <- "MF"
### MAS5
egoU <- enrichGO(mas5idU$ENTREZID, org.Hs.eg.db, readable=TRUE, ont=ont)
write.csv(as.data.frame(egoU), "./Data/GeneOntology/mas5U.csv")
egoU <- enrichGO(mas5idU$ENTREZID, org.Hs.eg.db, readable=FALSE, ont=ont)
egoU <- setReadable(egoU, org.Hs.eg.db, keyType="ENTREZID")
dotplot(egoU, title="mas5 Normailzation (Up-Regulated)")
barplot(egoU, title="mas5 Normailzation (Up-Regulated)")
plotGOgraph(egoU)
egoD <- enrichGO(mas5idD$ENTREZID, org.Hs.eg.db, readable=TRUE, ont=ont)
write.csv(as.data.frame(egoD), "./Data/GeneOntology/mas5D.csv")
egoD <- enrichGO(mas5idD$ENTREZID, org.Hs.eg.db, readable=FALSE, ont=ont)
egoD <- setReadable(egoD, org.Hs.eg.db, keyType="ENTREZID")
dotplot(egoD, title="mas5 Normailzation (Down-Regulated)")
barplot(egoD, title="mas5 Normailzation (Down-Regulated)")
plotGOgraph(egoD)
### RMA
egoU <- enrichGO(rmaidU$ENTREZID, org.Hs.eg.db, readable=TRUE, ont=ont)
write.csv(as.data.frame(egoU), "./Data/GeneOntology/rmaU.csv")
egoU <- enrichGO(rmaidU$ENTREZID, org.Hs.eg.db, readable=FALSE, ont=ont)
egoU <- setReadable(egoU, org.Hs.eg.db, keyType="ENTREZID")
dotplot(egoU, title="rma Normailzation (Up-Regulated)")
barplot(egoU, title="rma Normailzation (Up-Regulated)")
plotGOgraph(egoU)
egoD <- enrichGO(rmaidD$ENTREZID, org.Hs.eg.db, readable=TRUE, ont=ont)
write.csv(as.data.frame(egoD), "./Data/GeneOntology/rmaD.csv")
egoD <- enrichGO(rmaidD$ENTREZID, org.Hs.eg.db, readable=FALSE, ont=ont)
egoD <- setReadable(egoD, org.Hs.eg.db, keyType="ENTREZID")
dotplot(egoD, title="rma Normailzation (Down-Regulated)")
barplot(egoD, title="rma Normailzation (Down-Regulated)")
plotGOgraph(egoD)
### GCRMA
egoU <- enrichGO(gcrmaidU$ENTREZID, org.Hs.eg.db, readable=TRUE, ont=ont)
write.csv(as.data.frame(egoU), "./Data/GeneOntology/gcrmaU.csv")
egoU <- enrichGO(gcrmaidU$ENTREZID, org.Hs.eg.db, readable=FALSE, ont=ont)
egoU <- setReadable(egoU, org.Hs.eg.db, keyType="ENTREZID")
dotplot(egoU, title="gcrma Normailzation (Up-Regulated)")
barplot(egoU, title="gcrma Normailzation (Up-Regulated)")
plotGOgraph(egoU)
egoD <- enrichGO(gcrmaidD$ENTREZID, org.Hs.eg.db, readable=TRUE, ont=ont)
write.csv(as.data.frame(egoD), "./Data/GeneOntology/gcrmaD.csv")
egoD <- enrichGO(gcrmaidD$ENTREZID, org.Hs.eg.db, readable=FALSE, ont=ont)
egoD <- setReadable(egoD, org.Hs.eg.db, keyType="ENTREZID")
dotplot(egoD, title="gcrma Normailzation (Down-Regulated)")
barplot(egoD, title="gcrma Normailzation (Down-Regulated)")
plotGOgraph(egoD)


##### KEGG ANALYSIS #####
## MAS5
keggU <- enrichKEGG(mas5idU$ENTREZID)
write.csv(as.data.frame(keggU), "./Data/KEGG/mas5U.csv")
dotplot(keggU, title="mas5 Normailzation (Up-Regulated)")
keggD <- enrichKEGG(mas5idD$ENTREZID)
write.csv(as.data.frame(keggD), "./Data/KEGG/mas5D.csv")
dotplot(keggD, title="mas5 Normailzation (Down-Regulated)")
## RMA
keggU <- enrichKEGG(rmaidU$ENTREZID)
write.csv(as.data.frame(keggU), "./Data/KEGG/rmaU.csv")
dotplot(keggU, title="rma Normailzation (Up-Regulated)")
keggD <- enrichKEGG(rmaidD$ENTREZID)
write.csv(as.data.frame(keggD), "./Data/KEGG/rmaD.csv")
dotplot(keggD, title="rma Normailzation (Down-Regulated)")
## GCRMA
keggU <- enrichKEGG(gcrmaidU$ENTREZID)
write.csv(as.data.frame(keggU), "./Data/KEGG/gcrmaU.csv")
dotplot(keggU, title="gcrma Normailzation (Up-Regulated)")
keggD <- enrichKEGG(gcrmaidD$ENTREZID)
write.csv(as.data.frame(keggD), "./Data/KEGG/gcrmaD.csv")
dotplot(keggD, title="gcrma Normailzation (Down-Regulated)")


##### WIKIPATHWAYS ANALYSIS #####
extdata <- system.file("extdata/wikipathways-20180810-gmt-Homo_sapiens.gmt", package="clusterProfiler")
hsp <- read.gmt(extdata)
hsp <- hsp %>% tidyr::separate(term, c("name", "version", "wpid", "org"), "%")
wpgenes <- data.frame(wpid=hsp$wpid, gene=hsp$gene)
wpnames <- data.frame(wpid=hsp$wpid, name=hsp$name)
## enricher
### MAS5
enrich <- enricher(mas5idU$ENTREZID, TERM2GENE=wpgenes, TERM2NAME=wpnames)
enrich <- setReadable(enrich, org.Hs.eg.db, keyType="ENTREZID")
write.table(enrich@result, "./Data/wikiPathways/enricher/Umas5.txt")
enrich <- enricher(mas5idD$ENTREZID, TERM2GENE=wpgenes, TERM2NAME=wpnames)
enrich <- setReadable(enrich, org.Hs.eg.db, keyType="ENTREZID")
write.table(enrich@result, "./Data/wikiPathways/enricher/Dmas5.txt")
### RMA
enrich <- enricher(rmaidU$ENTREZID, TERM2GENE=wpgenes, TERM2NAME=wpnames)
enrich <- setReadable(enrich, org.Hs.eg.db, keyType="ENTREZID")
write.table(enrich@result, "./Data/wikiPathways/enricher/Urma.txt")
enrich <- enricher(rmaidD$ENTREZID, TERM2GENE=wpgenes, TERM2NAME=wpnames)
enrich <- setReadable(enrich, org.Hs.eg.db, keyType="ENTREZID")
write.table(enrich@result, "./Data/wikiPathways/enricher/Drma.txt")
### GCRMA
enrich <- enricher(gcrmaidU$ENTREZID, TERM2GENE=wpgenes, TERM2NAME=wpnames)
enrich <- setReadable(enrich, org.Hs.eg.db, keyType="ENTREZID")
write.table(enrich@result, "./Data/wikiPathways/enricher/Ugcrma.txt")
enrich <- enricher(gcrmaidD$ENTREZID, TERM2GENE=wpgenes, TERM2NAME=wpnames)
enrich <- setReadable(enrich, org.Hs.eg.db, keyType="ENTREZID")
write.table(enrich@result, "./Data/wikiPathways/enricher/Dgcrma.txt")
## GSEA
### MAS5
temp <- merge(mas5idU, mas5U, by="row.names")
genelist <- temp$y
names(genelist) <- as.character(temp$ENTREZID)
genelist <- sort(genelist, decreasing=TRUE)
gseaU <- GSEA(genelist, TERM2GENE=wpgenes, TERM2NAME=wpnames, minGSSize=1)
gseaU <- setReadable(gseaU, org.Hs.eg.db, keyType="ENTREZID")
saveRDS(gseaU, "./Data/wikiPathways/GSEA/Umas5GSEA.rds")
temp <- merge(mas5idD, mas5D, by="row.names")
genelist <- temp$y
names(genelist) <- as.character(temp$ENTREZID)
genelist <- sort(genelist, decreasing=TRUE)
gseaD <- GSEA(genelist, TERM2GENE=wpgenes, TERM2NAME=wpnames, minGSSize=1)
gseaD <- setReadable(gseaD, org.Hs.eg.db, keyType="ENTREZID")
saveRDS(gseaD, "./Data/wikiPathways/GSEA/Dmas5GSEA.rds")
### RMA
temp <- merge(rmaidU, rmaU, by="row.names")
genelist <- temp$y
names(genelist) <- as.character(temp$ENTREZID)
genelist <- sort(genelist, decreasing=TRUE)
gseaU <- GSEA(genelist, TERM2GENE=wpgenes, TERM2NAME=wpnames, minGSSize=1)
gseaU <- setReadable(gseaU, org.Hs.eg.db, keyType="ENTREZID")
saveRDS(gseaU, "./Data/wikiPathways/GSEA/UrmaGSEA.rds")
temp <- merge(rmaidD, rmaD, by="row.names")
genelist <- temp$y
names(genelist) <- as.character(temp$ENTREZID)
genelist <- sort(genelist, decreasing=TRUE)
gseaD <- GSEA(genelist, TERM2GENE=wpgenes, TERM2NAME=wpnames, minGSSize=1)
gseaD <- setReadable(gseaD, org.Hs.eg.db, keyType="ENTREZID")
saveRDS(gseaD, "./Data/wikiPathways/GSEA/DrmaGSEA.rds")
### GCRMA
temp <- merge(gcrmaidU, gcrmaU, by="row.names")
genelist <- temp$y
names(genelist) <- as.character(temp$ENTREZID)
genelist <- sort(genelist, decreasing=TRUE)
gseaU <- GSEA(genelist, TERM2GENE=wpgenes, TERM2NAME=wpnames, minGSSize=1)
gseaU <- setReadable(gseaU, org.Hs.eg.db, keyType="ENTREZID")
saveRDS(gseaU, "./Data/wikiPathways/GSEA/UgcrmaGSEA.rds")
temp <- merge(gcrmaidD, gcrmaD, by="row.names")
genelist <- temp$y
names(genelist) <- as.character(temp$ENTREZID)
genelist <- sort(genelist, decreasing=TRUE)
gseaD <- GSEA(genelist, TERM2GENE=wpgenes, TERM2NAME=wpnames, minGSSize=1)
gseaD <- setReadable(gseaD, org.Hs.eg.db, keyType="ENTREZID")
saveRDS(gseaD, "./Data/wikiPathways/GSEA/DgcrmaGSEA.rds")


##### GET IDS #####
## MAS5
write(as.vector(mas5idU$SYMBOL), "./Data/IDs/mas5Usymbols.txt")
write(as.vector(mas5idU$ENTREZID), "./Data/IDs/mas5Uentrezid.txt")
write(as.vector(mas5idD$SYMBOL), "./Data/IDs/mas5Dsymbols.txt")
write(as.vector(mas5idD$ENTREZID), "./Data/IDs/mas5Dentrezid.txt")
## RMA
write(as.vector(rmaidU$SYMBOL), "./Data/IDs/rmaUsymbols.txt")
write(as.vector(rmaidU$ENTREZID), "./Data/IDs/rmaUentrezid.txt")
write(as.vector(rmaidD$SYMBOL), "./Data/IDs/rmaDsymbols.txt")
write(as.vector(rmaidD$ENTREZID), "./Data/IDs/rmaDentrezid.txt")
## GCRMA
write(as.vector(gcrmaidU$SYMBOL), "./Data/IDs/gcrmaUsymbols.txt")
write(as.vector(gcrmaidU$ENTREZID), "./Data/IDs/gcrmaUentrezid.txt")
write(as.vector(gcrmaidD$SYMBOL), "./Data/IDs/gcrmaDsymbols.txt")
write(as.vector(gcrmaidD$ENTREZID), "./Data/IDs/gcrmaDentrezid.txt")
