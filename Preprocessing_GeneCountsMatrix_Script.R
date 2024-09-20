###### Preprocessing of Gene Counts Matrix

#### Input: Raw Gene Counts matrix obtained from reads processing and metadata

### Load libraries
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(openxlsx)
library(pheatmap)
#BiocManager::install("EnsDb.Hsapiens.v79")
library(EnsDb.Hsapiens.v79)


### Load the input files
counts_data_sra94 <- read.table("feature_counts_GSE107994.txt", header = TRUE, row.names = 1)
dim(counts_data_sra94)
head(counts_data_sra94)
srr_list94 <- readLines("SRR_Acc_List_GSE107994.txt") ######## Optional: Need to modify because the samples are named by SRR runs 

### Modify the raw matrix according to your needs

################ Optional: remove decimals from the Ensembl gene ids 
### changing the rownames
library(stringr)
rownames(data94_filtered) <- str_replace(rownames(data94_filtered),
                                         pattern = "\\..+$",
                                         replacement = "")
head(rownames(data94_filtered))   ## confirmation
#################

### Check for duplicate genes in rows
duplicatedGenes <- which(duplicated(rownames(counts_data_sra_clean94)) == T)
if (length(duplicatedGenes) > 0) {
  print("Duplicate ENSGs. Deleted entries:")
  print(counts_data_sra_clean94[duplicatedGenes,1:2]) ## examine duplicates
  counts_data_sra_clean94 <- counts_data_sra_clean94[-duplicatedGenes,] ## remove duplicates
}

dim(counts_data_sra_clean94) ## check if you lost any duplicates

### Filtering 1: Remove all the genes with null expression across all samples
row_counts <- rowSums(counts_data_sra_clean94)
data94_filtered <- counts_data_sra_clean94[row_counts > 0, ]
dim(counts_data_sra_clean94)
dim(data94_filtered)

### Filtering 2: Only protein coding genes to be used 

############ Optional step: using gene annotation for Homo_sapiens.GRCh38.94 obatined to get linking information between ensembl ids and gene biotype
library(dplyr)
gene_annotation <- read.table("Homo_sapiens.GRCh38.94_gene_annotation_table.txt", header = TRUE, sep = "\t") ### Input your gene annotation table 
dim(gene_annotation)
gene_annotation <- gene_annotation %>%
  mutate(
    gene_id = str_replace(gene_id, "gene_id ", ""),
    GeneSymbol = str_replace(GeneSymbol, "gene_name ", ""),
    Class = str_replace(Class, "gene_biotype ", ""),
    Strand = str_replace(Strand, "strand ", "")
  )
protein_coding_genes_table <- gene_annotation %>% ### selecting only protein coding genes
  filter(Class == "protein_coding") %>%
  dplyr::select(GeneSymbol, gene_id)
head(protein_coding_genes_table)
#############

### filtering on the basis of protein coding genes
protein_coding_genes <- protein_coding_genes_table$gene_id
filtered_matrix94 <- data94_filtered[rownames(data94_filtered) %in% protein_coding_genes, ]
dim(filtered_matrix94)



##### Normalization and Log trasnformation using edgeR
dge <- DGEList(counts = filtered_matrix94)
dim(dge)
keep <- filterByExpr(dge) ### Filtering 3: By low gene expression
dge <- dge[keep, ]
dim(dge)
dge <- calcNormFactors(dge, method = "TMM")  ## Normalization using TMM
z <- cpm(dge, log = TRUE)  #### applies the TMM factors for normalization and log transforms the counts for downstream analysis 
normalized_tmm_edger <- z
normalized_tmm_edger <- as.matrix(normalized_tmm_edger)
dim(normalized_tmm_edger)


###Plot histogram to visualize the chnange in distribution before and after normalization
par(mfrow = c(1, 2))
hist(as.matrix(data94_filtered), main = "Raw Data",
     xlab = "Raw Values", col = "lightgreen")
hist(normalized_tmm_edger, main = "TMM EdgeR", 
     xlab = "Normalized values", col = "lightpink")


