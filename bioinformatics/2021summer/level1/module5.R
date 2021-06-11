setwd('/Volumes/SamsungUSB/stemaway/bi-mentoring/level1/2021summer')

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
tT <- read.table("dge/rma-tt.txt", row.names=1)
DEG.genes <- tT$logFC
names(DEG.genes) <- tT$ID
DEG.genes <- DEG.genes[DEG.genes > 1.5]
DEG.genes <- sort(DEG.genes, decreasing=T)

DEG.entrez <- AnnotationDbi::select(org.Hs.eg.db, keys=names(DEG.genes), columns=c("ENTREZID"), keytype="SYMBOL")

###

tT.selected <- read.table("dge/rma-tt-selected.txt", row.names=1)
DEG.genes.selected <- tT.selected$logFC
names(DEG.genes.selected) <- tT.selected$ID
DEG.genes.selected <- DEG.genes.selected[DEG.genes.selected > 1.5]
DEG.genes.selected <- sort(DEG.genes.selected, decreasing=T)

DEG.entrez.selected <- AnnotationDbi::select(org.Hs.eg.db, keys=names(DEG.genes.selected), columns=c("ENTREZID"), keytype="SYMBOL")




### FUNCTIONAL ENRICHMENT ANALYSIS
##### GENE ONTOLOGY
egoCC <- enrichGO(DEG.entrez$ENTREZID, org.Hs.eg.db, ont="CC", readable=T)
egoBP <- enrichGO(DEG.entrez$ENTREZID, org.Hs.eg.db, ont="BP", readable=T)
egoMF <- enrichGO(DEG.entrez$ENTREZID, org.Hs.eg.db, ont="MF", readable=T)

barplot(egoCC, title="CC")
barplot(egoBP, title="BP")
barplot(egoMF, title="MF")

###

egoCC.selected <- enrichGO(DEG.entrez.selected$ENTREZID, org.Hs.eg.db, ont="CC", readable=T)
egoBP.selected <- enrichGO(DEG.entrez.selected$ENTREZID, org.Hs.eg.db, ont="BP", readable=T)
egoMF.selected <- enrichGO(DEG.entrez.selected$ENTREZID, org.Hs.eg.db, ont="MF", readable=T)

barplot(egoCC.selected, title="CC (selected samples)")
barplot(egoBP.selected, title="BP (selected samples)")
barplot(egoMF.selected, title="MF (selected samples)")



##### KEGG PATHWAYS
ek <- enrichKEGG(DEG.entrez$ENTREZID)
dotplot(ek, title="Enriched KEGG Pathways")

###

ek.selected <- enrichKEGG(DEG.entrez.selected$ENTREZID)
dotplot(ek.selected, title="Enriched KEGG Pathways (selected samples)")



### GENE-CONCEPT NETWORKS
edR <- setReadable(ek, org.Hs.eg.db, keyType="ENTREZID")
cnetplot(edR, foldChange=DEG.genes, categorySize="pvalue", colorEdge=T)

###

edR.selected <- setReadable(ek.selected, org.Hs.eg.db, keyType="ENTREZID")
cnetplot(edR.selected, foldChange=DEG.genes.selected, categorySize="pvalue", colorEdge=T)



### GSEA
##### GENE VECTOR
GSEA.entrez <- AnnotationDbi::select(org.Hs.eg.db, keys=tT$ID, columns=c("ENTREZID"), keytype="SYMBOL")
GSEA.entrez <- GSEA.entrez[!duplicated(GSEA.entrez$SYMBOL),]

row.names(GSEA.entrez) <- GSEA.entrez$SYMBOL
GSEA.genes <- merge(GSEA.entrez, tT, by="row.names")

genelist <- GSEA.genes$logFC
names(genelist) <- GSEA.genes$ENTREZID
genelist <- sort(genelist, decreasing=T)

###

GSEA.entrez.selected <- AnnotationDbi::select(org.Hs.eg.db, keys=tT.selected$ID, columns=c("ENTREZID"), keytype="SYMBOL")
GSEA.entrez.selected <- GSEA.entrez.selected[!duplicated(GSEA.entrez.selected$SYMBOL),]

row.names(GSEA.entrez.selected) <- GSEA.entrez.selected$SYMBOL
GSEA.genes.selected <- merge(GSEA.entrez.selected, tT.selected, by="row.names")

genelist.selected <- GSEA.genes.selected$logFC
names(genelist.selected) <- GSEA.genes.selected$ENTREZID
genelist.selected <- sort(genelist.selected, decreasing=T)



### GSE ANALYSIS
msig <- msigdbr(species="Homo sapiens", category="H")
h <- msig %>% select(gs_name, entrez_gene)

gsea <- GSEA(genelist, TERM2GENE=h, eps=0)
gseaplot2(gsea, geneSetID=1:5, pvalue_table=T)

###

gsea.selected <- GSEA(genelist.selected, TERM2GENE=h, eps=0)
gseaplot2(gsea.selected, geneSetID=1:5, pvalue_table=T)



## TRANSCRIPTIONAL FACTOR ANALYSIS
msig <- msigdbr(species="Homo sapiens", category="C3")
c3 <- msig %>% select(gs_name, entrez_gene)

e <- enricher(names(genelist[genelist > 1]), TERM2GENE=c3)
eR <- setReadable(e, org.Hs.eg.db, keyType="ENTREZID")
cnetplot(eR, foldChange=DEG.genes, categorySize="pvalue", colorEdge=T)

###

e.selected <- enricher(names(genelist.selected[genelist.selected > 1]), TERM2GENE=c3)
eR.selected <- setReadable(e.selected, org.Hs.eg.db, keyType="ENTREZID")
cnetplot(eR.selected, foldChange=DEG.genes.selected, categorySize="pvalue", colorEdge=T)



# EXTERNAL TOOLS
write(names(DEG.genes), "fa-1/up-deg.txt")

###

write(names(DEG.genes.selected), "fa-1/up-deg-selected.txt")
