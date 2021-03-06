setwd("C:/Users/durbe/OneDrive/Documents/ReposExtras/STEM-Away/June2020")
library(affy)
library(affyPLM)
library(sva)
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
## AFFYPLM - compute RLSE and NUSE scores of microarray samples
pset <- fitPLM(data, background=FALSE, normalize=TRUE)
rle <- RLE(pset, type="stats")
rle <- data.frame(Median=rle[1,])
ggplot(rle, aes(Median))+ geom_histogram()+ xlab("Median")+ ylab("Frequency")+ ggtitle("RLE Histogram (Raw)")+ theme_bw()
nuse <- NUSE(pset, type="stats")
nuse <- data.frame(Median=nuse[1,])
ggplot(nuse, aes(Median))+ geom_histogram()+ xlab("Median")+ ylab("Frequency")+ ggtitle("NUSE Histogram (Raw)")+ theme_bw()


##### DATA PREPROCESSING: BACKGROUND CORRECTION + NORMALIZATION #####
gcrmaN <- gcrma(data)
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
gcrmaNdf <- gcrmaNdf[-c(13,43)]
gcrmaP <- select(hgu133plus2.db, keys=rownames(gcrmaNdf), columns=c("SYMBOL"), keytype="PROBEID")
gcrmaP <- gcrmaP[!duplicated(gcrmaP$PROBEID),]
gcrmaP <- gcrmaP[!duplicated(gcrmaP$SYMBOL),]
rownames(gcrmaP) <- gcrmaP$PROBEID
gcrmaS <- merge(gcrmaNdf, gcrmaP["SYMBOL"], by="row.names")
gcrmaS <- na.omit(gcrmaS)


##### GENE FILTERING #####
means <- rowMeans(gcrmaS[,2:63])
perc4 <- as.numeric(quantile(means, probs=0.04, na.rm=TRUE))
gcrmaF <- gcrmaS[which(means > perc4),]
colnames(gcrmaF)[1] <- "PROBEID"
write.csv(gcrmaF, "./Data/Filtered/gcrmaF.csv")


##### LIMMA #####
type <- factor(paste(rep(c("Control", "Cancer"), each=31)))
design <- model.matrix(~0+type)
max <- length(mas5Ndf$GSM215051)

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
ggoU <- groupGO(gcrmaidU$ENTREZID, org.Hs.eg.db, readable=FALSE, level=i, ont=ont)
ggoD <- groupGO(gcrmaidD$ENTREZID, org.Hs.eg.db, readable=FALSE, level=i, ont=ont)
barplot(ggoU, title=paste("gcrma Normalization (Up-Regulated), Level=", i, sep=""))
barplot(ggoD, title=paste("gcrma Normalization (Down-Regulated), Level=", i, sep=""))
## enrichGO
ont <- "MF"
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
enrich <- enricher(gcrmaidU$ENTREZID, TERM2GENE=wpgenes, TERM2NAME=wpnames)
enrich <- setReadable(enrich, org.Hs.eg.db, keyType="ENTREZID")
write.table(enrich@result, "./Data/wikiPathways/enricher/Ugcrma.txt")
enrich <- enricher(gcrmaidD$ENTREZID, TERM2GENE=wpgenes, TERM2NAME=wpnames)
enrich <- setReadable(enrich, org.Hs.eg.db, keyType="ENTREZID")
write.table(enrich@result, "./Data/wikiPathways/enricher/Dgcrma.txt")
## GSEA
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
write(as.vector(gcrmaidU$SYMBOL), "./Data/IDs/gcrmaUsymbols.txt")
write(as.vector(gcrmaidU$ENTREZID), "./Data/IDs/gcrmaUentrezid.txt")
write(as.vector(gcrmaidD$SYMBOL), "./Data/IDs/gcrmaDsymbols.txt")
write(as.vector(gcrmaidD$ENTREZID), "./Data/IDs/gcrmaDentrezid.txt")
