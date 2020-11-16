setwd("C:/Users/durbe/OneDrive/Documents/ReposExtras/STEM-Away/June2020/Final/")

library(affy)
library(arrayQualityMetrics)
library(gcrma)
library(GEOquery)
library(ggplot2)
library(pheatmap)

library(limma)
library(EnhancedVolcano)

library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)


gse <- ReadAffy(compress=T, celfile.path="./Data/")
arrayQualityMetrics(gse, "./arrayQualityMetrics/", force=T, do.logtransform=T)

# Normalization
gcrma <- gcrma(gse)
gcrma <- as.data.frame(exprs(gcrma))
for (i in 1:44) {
  colnames(gcrma)[i] <- strsplit(colnames(gcrma), "_")[[i]][1]
}
gcrma <- gcrma[, -c(1,2,8,25)]
write.csv(gcrma, "./Data/gcrma.csv")

# Metadata
meta <- getGEO(filename="./Data/metadata.txt")
meta <- meta@phenoData@data
meta <- meta[40]
meta <- meta[105:148,]
meta <- meta[-c(1,2,8,25)]
for (i in 1:40) {
  meta[i] <- strsplit(meta[i], ",")[[1]][1]
}
names(meta) <- colnames(gcrma)
write.table(meta, "./Data/metadata-cleaned.txt", sep="\t")

# Annotation
annot <- AnnotationDbi::select(hgu133plus2.db, keys=rownames(gcrma), columns=c("SYMBOL"), keytype="PROBEID")
annot <- annot[!duplicated(annot$PROBEID),]
annot <- annot[!duplicated(annot$SYMBOL),]
rownames(annot) <- annot$PROBEID
annot <- merge(gcrma, annot["SYMBOL"], by="row.names")
row.names(annot) <- annot$Row.names
annot <- annot[-1]
annot <- na.omit(annot)
write.csv(annot, "./Data/annotated.csv")

# Filtering
means <- rowMeans(annot[,1:40])
perc4 <- as.numeric(quantile(means, probs=0.04))
filt <- annot[which(means > perc4),]
write.csv(filt, "./Data/filtered.csv")

# DEG Analysis
design <- model.matrix(~0+meta)
max <- length(gcrma$GSM549101)

genelist <- filt$SYMBOL
deg <- filt[-41]
row.names(design) <- names(meta)

lm <- lmFit(deg, design)

contrast <- makeContrasts(metacancer-metanormal, levels=design)
colnames(contrast) <- c("Cancer-Normal")

fit <- contrasts.fit(lm, contrast)
fit <- eBayes(fit)

tT100 <- topTable(fit, genelist=genelist, adjust.method="fdr", sort.by="B", number=100)
tT <- topTable(fit, genelist=genelist, adjust.method="fdr", sort.by="B", number=max)
write.table(tT, "./Data/DEG.txt")

# Volcano Plot
tT$DEG <- ifelse(tT$logFC > 0, "Up-Regulated", "Down-Regulated")
tT$NegP <- ifelse(tT$adj.P.Val>0, 0-log(tT$adj.P.Val), abs(log(tT$adj.P.Val)))

vp <- tT[which(abs(tT$logFC) > 1),]
vp <- vp[which(abs(tT$adj.P.Val) < 0.0000000001),]
ggplot(vp, aes(x=logFC, NegP, color=DEG))+
  geom_point()+
  geom_vline(xintercept=1, linetype="dashed")+
  geom_vline(xintercept=-1, linetype="dashed")+
  scale_color_manual(values=c("springgreen4", "mediumpurple"))+
  ggtitle("Volcano Plot")+
  labs(x="LogFC", y="-LogAdj.P")+ ylim(20,35)+
  theme_bw()

# Heatmap
hm <- as.data.frame(tT100$ID)
colnames(hm) <- c("SYMBOL")
hm <- merge(filt, hm, by="SYMBOL")
row.names(hm) <- hm$SYMBOL
hm <- hm[-1]
colnames(hm) <- as.factor(meta)

thm <- as.data.frame(t(hm))
col.pal <- RColorBrewer::brewer.pal(9, "Purples")
pheatmap(thm, color=col.pal, cluster_rows=F,  
         main="Hierarchical Clustering Map of Expression Patterns")

# Gene Vector
deg <- read.table("./Data/DEG.txt")
genes <- deg$logFC
names(genes) <- deg$ID

up <- sort(genes[genes > 1], decreasing=T)
upID <- bitr(names(up), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

down <- sort(genes[genes < -1], decreasing=T)
downID <- bitr(names(down), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# Gene Ontology
eGO <- function(d, o, t) {
  ego <- enrichGO(d$ENTREZID, org.Hs.eg.db, readable=F, ont=o, pool=T)
  barplot(ego, title=t)
}

eGO(upID, "CC", "Up-Regulated: CC")
eGO(upID, "MF", "Up-Regulated: MF")
eGO(upID, "BP", "Up-Regulated: BP")

eGO(downID, "CC", "Down-Regulated: CC")
eGO(downID, "MF", "Down-Regulated: MF")
eGO(downID, "BP", "Down-Regulated: BP")

# ID Stuff
write(as.vector(upID$ENTREZID), "./Data/upID.txt")
write(as.vector(downID$ENTREZID), "./Data/downID.txt")


