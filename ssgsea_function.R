#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("4 arguments must be supplied ()",
   call.=FALSE)
}


####### Parameter of main function
path_Rlibrary <- args[1]
countscsv <- args[2]
genesetrds <- args[3]
output <- args[4]

.libPaths(path_Rlibrary, FALSE)
print(.libPaths())

library(readr)
library(data.table)
library(GSVA)
library(dplyr)

print("GSVA version is")
packageVersion("GSVA")

## read csv file as counts
counts = readRDS(countscsv)
geneset = readRDS(genesetrds)

if (!is.matrix(counts) || !is.numeric(counts)) {
  stop("counts_matrix should be a numeric matrix")
}
if (!is.list(geneset)) {
  stop("gene_sets should be a list")
}

ssgsea_results <- gsva(counts, geneset, method = "ssgsea")


saveRDS(ssgsea_results, file = output)
print("Results saved!")
