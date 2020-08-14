gsealist <- function(d) {
  gene <- d$logFC
  names(gene) <- d$ID
  
  DEG <- AnnotationDbi::select(org.Hs.eg.db, keys=names(gene), columns=c("ENTREZID"), keytype="SYMBOL")
  DEG <- DEG[!duplicated(DEG$SYMBOL),]
  row.names(DEG) <- DEG$SYMBOL
  DEG <- merge(DEG, gene, by="row.names")
  
  genelist <- DEG$y
  names(genelist) <- DEG$ENTREZID
  genelist <- sort(genelist, decreasing=T)
  
  return(genelist)
}

mGene <- gsealist(mas5)
rGene <- gsealist(rma)
gGene <- gsealist(gcrma)

msig <- msigdbr(species="Homo sapiens", category="H")
msigS <- msig %>% select(gs_name, entrez_gene)

gseaM <- GSEA(mGene, TERM2GENE=msigS)
gseaR <- GSEA(rGene, TERM2GENE=msigS)
gseaG <- GSEA(gGene, TERM2GENE=msigS)

gseaplot2(gseaM, geneSetID=1:5, pvalue_table=T)
gseaplot2(gseaR, geneSetID=1:5, pvalue_table=T)
gseaplot2(gseaG, geneSetID=1:5, pvalue_table=T)


e <- enricher(names(mGene), TERM2GENE=msigS)
eM <- setReadable(e, org.Hs.eg.db, keyType="ENTREZID")
heatplot(eM, foldChange=mGene)

e <- enricher(names(rGene), TERM2GENE=msigS)
eR <- setReadable(e, org.Hs.eg.db, keyType="ENTREZID")
heatplot(eR, foldChange=rGene)

e <- enricher(names(gGene), TERM2GENE=msigS)
eG <- setReadable(e, org.Hs.eg.db, keyType="ENTREZID")
heatplot(eG, foldChange=gGene)


### EXTERNAL TOOLS
write(names(mGene), "./FA/data/genelist/mas5_all.txt")
write(names(rGene), "./FA/data/genelist/rma_all.txt")
write(names(gGene), "./FA/data/genelist/gcrma_all.txt")
