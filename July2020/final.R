setwd("C:/Users/durbe/OneDrive/Documents/ReposExtras/STEM-Away/Final/")

##### METADATA
# library(GEOquery)
# geo <- getGEO(filename="./data/GSE68417_series_matrix.txt.gz")
# pheno <- geo@phenoData@data
# meta <- pheno$`analysis group:ch1`
# names(meta) <- rownames(pheno)
# meta <- sub("normal", "Normal", meta)
# write.table(meta, "./data/meta.txt")
meta <- read.table("./data/meta.txt", header=T, row.names=1)


##### QUALITY CONTROL
library(oligo)
library(ggplot2)
library(pheatmap)

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68417
gse68417 <- read.celfiles(list.files("./data/raw/GSE68417/", full.names=T)) 

norm <- rma(gse68417)
colnames(norm) <- rownames(meta)


## VISUALIZATION
boxplot(gse68417, which="all")
boxplot(norm, which="all")

plm <- fitProbeLevelModel(gse68417)
RLE(plm)
NUSE(plm)

Tissue <- as.factor(meta$x)
pca <- prcomp(norm, scale=F, center=F)
pca <- as.data.frame(pca$rotation)
ggplot(pca, aes(x=PC1, y=PC2, color=Tissue))+
  geom_point()+ stat_ellipse()+
  scale_colour_manual(values=c("hotpink1", "navy", "goldenrod2", "aquamarine4"))+
  ggtitle("PCA")+ theme_bw()

dismat <- 1-cor(exprs(norm))
Tissue <- data.frame(Tissue=factor(meta$x))
rownames(Tissue) <- rownames(meta)
pheatmap(dismat, cluster_rows=F, annotation_col=Tissue, main="Heatmap")


##### DGE
library(hugene10sttranscriptcluster.db)
library(WGCNA)
library(limma)
library(EnhancedVolcano)

## ANNOTATION
symbols <- AnnotationDbi::select(hugene10sttranscriptcluster.db, keys=rownames(norm), columns=c("SYMBOL"), keytype="PROBEID")
symbols <- symbols[!duplicated(symbols$PROBEID),]
rownames(symbols) <- symbols$PROBEID
annot <- merge(symbols, exprs(norm), by="row.names")
rownames(annot) <- annot$Row.names
annot <- annot[-c(1,2)]
annot <- na.omit(annot)
annot <- collapseRows(annot[-1], rowGroup=annot$SYMBOL, rowID=rownames(annot))
annot <- annot$datETcollapsed


## FILTERING
means <- rowMeans(annot)
perc2 <- as.numeric(quantile(means, probs=0.02, na.rm=T))
filt <- annot[which(means > perc2),]
# write.csv(filt, "./data/filt.csv")


## LIMMA
meta$Type <- c(rep("Normal", 6), rep("Cancer", 29), rep("Normal",14))
Tissue <- as.factor(meta$Type)
design <- model.matrix(~0+Tissue, meta)

lm <- lmFit(filt, design)
contrast <- makeContrasts(TissueCancer-TissueNormal, levels=design)
fit <- contrasts.fit(lm, contrast)
fit <- eBayes(fit)
tT <- topTable(fit, adjust="fdr", sort.by="B", number=max(rownames(filt)))
write.table(tT, "./data/DGE.txt")


## VISUALIZATION
EnhancedVolcano(tT, lab=rownames(tT), x="logFC", y="adj.P.Val",
                pointSize=1, legendLabSize=10, labSize=3.0,
                title="Volcano Plot")

t50 <- topTable(fit, adjust="fdr", sort.by="B", number=50)
input <- filt[(rownames(filt) %in% rownames(t50)),]
Tissue <- data.frame(Tissue=factor(meta$x))
rownames(Tissue) <- colnames(filt)
pheatmap(input, culster_rows=F, annotation_col=Tissue, main="Top 50 Genes")


##### FUNCTIONAL ANALYSIS
library(org.Hs.eg.db)

geneLF <- tT$logFC
names(geneLF) <- rownames(tT)
DEG.up <- geneLF[geneLF > 1.0]
DEG.up <- sort(DEG.up, decreasing=T)

DEG <- AnnotationDbi::select(org.Hs.eg.db, keys=names(DEG.up), columns=c("ENTREZID"), keytype="SYMBOL")


## GENE ONTOLOGY
library(clusterProfiler)

ego <- enrichGO(DEG$ENTREZID, org.Hs.eg.db, ont="CC", readable=F)
egoR <- setReadable(ego, org.Hs.eg.db, keyType="ENTREZID")
barplot(egoR, title="CC")

ego <- enrichGO(DEG$ENTREZID, org.Hs.eg.db, ont="BP", readable=F)
egoR <- setReadable(ego, org.Hs.eg.db, keyType="ENTREZID")
barplot(egoR, title="BP")

ego <- enrichGO(DEG$ENTREZID, org.Hs.eg.db, ont="MF", readable=F)
egoR <- setReadable(ego, org.Hs.eg.db, keyType="ENTREZID")
barplot(egoR, title="MF")


## KEGG
ek <- enrichKEGG(DEG$ENTREZID)
dotplot(ek, title="KEGG")


## GENE-CONCEPT NETWORK
library(enrichplot)
ekR <- setReadable(ek, org.Hs.eg.db, keyType="ENTREZID")
cnetplot(ekR, foldchange=DEG.up, categorySize="pvalue", colorEdge=T)


## TF ANALYSIS
library(msigdbr)
msig <- msigdbr(species="Homo sapiens", category="C3")
c3 <- msig %>% select(gs_name, entrez_gene)

e <- enricher(DEG$ENTREZID, TERM2GENE=c3)
eR <- setReadable(e, org.Hs.eg.db, keyType="ENTREZID")
cnetplot(eR, foldchange=DEG.up, categorySize="pvalue", colorEdge=T)


## EXTERNAL TOOLS
write(names(DEG.up), "./data/genelist.txt")


## GSEA
library(enrichplot)
genes <- sort(geneLF, decreasing=T)
DEG <- AnnotationDbi::select(org.Hs.eg.db, keys=names(genes), columns=c("ENTREZID"), keytype="SYMBOL")
DEG <- DEG[!duplicated(DEG$SYMBOL),]
names(genes) <- DEG$ENTREZID
  
msig <- msigdbr(species="Homo sapiens", category="H")
h <- msig %>% select(gs_name, entrez_gene)

gsea <- GSEA(genes, TERM2GENE=h)
gseaplot2(gsea, geneSetID=1:5, pvalue_table=T)

e <- enricher(names(genes), TERM2GENE=h)
eR <- setReadable(e, org.Hs.eg.db, keyType="ENTREZID")
heatplot(eR, foldChange=genes) #[0x9]
