setwd("C:/Users/durbe/Documents/repositories/StemAway/Deliverables")

# 1. download series matrix file from GEO
gse <- getGEO(filename="./data/GSE8671_series_matrix.txt.gz")

filt <- read.csv("./data/filt.csv")
row.names(filt) <- filt$X
filt <- filt[-1]

gset <- ExpressionSet(assayData = as.matrix(filt))

gset@phenoData@data <- gse@phenoData@data

# 2. copy meta data
meta <- as.data.frame(gset@phenoData@data)

# 3. remove non-relevant columns
meta <- meta[-c(13,43),]

# 4. rename rest of columns with meaningful names

# 5. remove unwanted strings from cells

# 6. save sheet as tab-delimited
write.csv(meta, "./meta/meta.csv")
write.table(meta, "./meta/meta.tsv", sep="\t")

# 7. read file into R using read.table() function

# 8. summarize(phenodata)

# 9. 