setwd('D:/stemaway/bioinformatics/level2')

library(GEOquery)
library(affy)

# Get expression data from .CEL.gz files and save as RDS
gse <- ReadAffy(celfile.path='data/raw/', compress=TRUE)

# Reformat sample names
for (i in 1:length(colnames(exprs(gse)))) {
    row.names(gse@phenoData@data)[i] <- strsplit(row.names(gse@phenoData@data)[i],'.CEL')[[1]][1]
    row.names(gse@protocolData@data)[i] <- strsplit(row.names(gse@protocolData@data)[i], '.CEL')[[1]][1]
    colnames(exprs(gse))[i] <- strsplit(colnames(exprs(gse))[i], '.CEL')[[1]][1] #RLE/NUSE
}
  # #
saveRDS(gse, 'data/affybatch.rds')

# Get metadata from series matrix file
series <- getGEO('GSE19804', GSEMatrix=TRUE)
meta_df <- series$GSE19804_series_matrix.txt.gz@phenoData@data
colnames(meta_df)

# Select features of interest
meta <- meta_df[c('tissue:ch1', 'stage:ch1')]
colnames(meta) <- c('TISSUE', 'STAGE')

# Reformat text
meta$TISSUE <- replace(meta$TISSUE, meta$TISSUE != 'lung cancer', 'NORMAL')
meta$TISSUE <- replace(meta$TISSUE, meta$TISSUE == 'lung cancer', 'CANCER')
meta$STAGE <- replace(meta$STAGE, meta$STAGE == 'n/a', '0')
View(meta)

write.csv(meta, 'data/metadata.csv', quote=FALSE)
