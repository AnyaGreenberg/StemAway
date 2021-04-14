setwd("/Volumes/SamsungUSB/STEM-Away/BI-mentoring/Level1/SA-BI-Project")

library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library(msigdbr)
# library(topGO)
# library(pathview)
# library(magrittr)
# library(tidyr)




### GENE VECTOR
tT <- read.table("DGE/rma.txt", row.names=1)
DEG.genes <- tT$logFC
names(DEG.genes) <- tT$ID
DEG.genes <- DEG.genes[DEG.genes > 1.5]
DEG.genes <- sort(DEG.genes, decreasing=T)

DEG.entrez <- AnnotationDbi::select(org.Hs.eg.db, keys=names(DEG.genes), columns=c("ENTREZID"), keytype="SYMBOL")



### FUNCTIONAL ENRICHMENT ANALYSIS
##### GENE ONTOLOGY
egoCC <- enrichGO(DEG.entrez$ENTREZID, org.Hs.eg.db, ont="CC", readable=T)
egoBP <- enrichGO(DEG.entrez$ENTREZID, org.Hs.eg.db, ont="BP", readable=T)
egoMF <- enrichGO(DEG.entrez$ENTREZID, org.Hs.eg.db, ont="MF", readable=T)

barplot(egoCC, title="CC")
barplot(egoBP, title="BP")
barplot(egoMF, title="MF")


##### KEGG PATHWAYS
ek <- enrichKEGG(DEG.entrez$ENTREZID)
dotplot(ek, title="Enriched KEGG Pathways")



### GENE-CONCEPT NETWORKS
edR <- setReadable(ek, org.Hs.eg.db, keyType="ENTREZID")
cnetplot(edR, foldChange=DEG.genes, categorySize="pvalue", colorEdge=T)



### GSEA
##### GENE VECTOR
GSEA.entrez <- AnnotationDbi::select(org.Hs.eg.db, keys=tT$ID, columns=c("ENTREZID"), keytype="SYMBOL")
GSEA.entrez <- GSEA.entrez[!duplicated(GSEA.entrez$SYMBOL),]

row.names(GSEA.entrez) <- GSEA.entrez$SYMBOL
GSEA.genes <- merge(GSEA.entrez, tT, by="row.names")

genelist <- GSEA.genes$logFC
names(genelist) <- GSEA.genes$ENTREZID
genelist <- sort(genelist, decreasing=T)



### GSE ANALYSIS
msig <- msigdbr(species="Homo sapiens", category="H")
h <- msig %>% select(gs_name, entrez_gene)

gsea <- GSEA(genelist, TERM2GENE=h, eps=0)
gseaplot2(gsea, geneSetID=1:5, pvalue_table=T)



## TRANSCRIPTIONAL FACTOR ANALYSIS
msig <- msigdbr(species="Homo sapiens", category="C3")
c3 <- msig %>% select(gs_name, entrez_gene)

e <- enricher(names(genelist[genelist > 1]), TERM2GENE=c3)
eR <- setReadable(e, org.Hs.eg.db, keyType="ENTREZID")
cnetplot(eR, foldChange=DEG.genes, categorySize="pvalue", colorEdge=T)



# EXTERNAL TOOLS
write(names(DEG.genes), "FA-1/upDEG.txt")
