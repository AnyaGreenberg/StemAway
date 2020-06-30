setwd("C:/Users/durbe/Documents/repositories/StemAway/Deliverables")
library(GEOquery)
library(affy)

# 1. download series matrix file from GEO 
gse <- getGEO(filename="./meta/GSE8671_series_matrix.txt.gz")

# 2. copy meta data and transpose into a new sheet
meta <- gse@phenoData@data
meta <- meta[-c(13,43),]

# 3. remove non-relevant columns/rows
meta <- meta[,c(1,2,42,44,45)]

# 4. rename rest of columns with meaningful names (avoid spaces)
colnames(meta) <- c("ptID", "GEOassession", "location","size", "tissue")

# 5. remove unwanted strings from cells (replace empty strings)
for (i in 1:length(meta$ptID)) {
  temp <- strsplit(meta$ptID[i], split="#")
  meta$ptID[i] <- temp[[1]][2]
}

# 6. save sheet as tab-delimited (.txt)
## write.csv(meta, "./meta/meta.csv")
write.table(meta, "./meta/meta.txt", sep="\t")

# 7. read file into R using read.table() function
meta <- read.table("./meta/meta.txt", sep="\t", header=TRUE)

# 8. have a look into the object - summarize(phenodata)
s <- summary(meta)

# 9. read in the .CEL files, introduce phenotypic data into the created expression set using the argument pData(expressionSet) <- new("AnnotatedDataFrame", data=object)
data <- ReadAffy(compress=TRUE, celfile.path="./data/raw/")
raw <- exprs(data)

pData(data) <- AnnotatedDataFrame(data=meta)

# 10. save the expressionSet object (including both expression matrix and phenodata) using saveRDS() function
saveRDS(phenoData, "./meta/pheno.csv")
