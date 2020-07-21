setwd("C:/Users/durbe/OneDrive/Documents/ReposExtras/STEM-Away/June2020/Final/")

library(affy)
library(arrayQualityMetrics)
library(gcrma)
library(ggplot2)
library(pheatmap)

gse <- ReadAffy(compress=T, celfile.path="./Data/")
arrayQualityMetrics(gse, "./arrayQualityMetrics-Raw/", force=T, do.logtransform=T)
