###### Continued Exploratory Data Analysis: Differential Gene Expression Analysis 

### Inputs required: raw matrix(filtered only-not normalized) and the group information containing dataframe all in rds format


### Load libraries
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(pheatmap)
library(VennDiagram)
library(UpSetR)

###### Functions used in the script

### 1. Function to create volcano plot
create_volcano_plot <- function(res, title) {
  res$significant <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "Significant", "Not Significant")
  num_sig_genes <- sum(res$significant == "Significant", na.rm = TRUE)
  
  
  ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("grey30", "red2")) +
    labs(title = title, 
         x = "Log2 Fold Change", 
         y = "-Log10 Adjusted P-value",
         color = paste("Number of Significant\n Genes:", num_sig_genes )) + # Update legend title
    theme_minimal()
  
}

### Load datasets: filtered gene counts matrix (not normalized) and group information
filtered_matrix94 <- readRDS("lowcounts_filtered94.rds")
dim(filtered_matrix94)
group_info94 <- readRDS("Group_information.rds")
dim(group_info94)


### Differential Expression Analysis using DESeq2
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = filtered_matrix94,
                                colData = group_info94,
                                design = ~ Group, 
                                tidy = FALSE,
                                ignoreRank = FALSE,
                                rownames(filtered_matrix94))
print(dds)
### Perform DGE
dds <- DESeq(dds)
### Result generation for different comparisons
res_active_tb <- results(dds, contrast = c("Group", "Active_TB", "Control"))
res_ltbi <- results(dds, contrast = c("Group", "LTBI_Non_Progressor", "Control"))
res_ltbi_progressor <- results(dds, contrast = c("Group", "LTBI_Progressor", "Control"))
res_active_ltbinonprog <- results(dds, contrast = c("Group", "Active_TB", "LTBI_Non_Progressor"))
### Total number of significant genes for each comparison: based on FDR adjusted pvalue < 0.05 and logFC > 1 or logFC < -1
sum(res_active_tb$padj < 0.05 & abs(res_active_tb$log2FoldChange) > 1, na.rm = TRUE)
sum(res_ltbi$padj < 0.05 & abs(res_ltbi$log2FoldChange) > 1, na.rm = TRUE)
sum(res_ltbi_progressor$padj < 0.05 & abs(res_ltbi_progressor$log2FoldChange) > 1, na.rm = TRUE)
sum(res_active_ltbinonprog$padj < 0.05 & abs(res_active_ltbinonprog$log2FoldChange) > 1, na.rm = TRUE)
### create volcano plot for all
create_volcano_plot(res_active_tb, "Active TB vs Control")
create_volcano_plot(res_ltbi, "LTBI Non-Progressor vs Control")
create_volcano_plot(res_ltbi_progressor, "LTBI Progressor vs Control")
create_volcano_plot(res_active_ltbinonprog, "Active TB vs LTBI Non Progressor")



### More DGE on subset of the dataset
subset_colData <- group_info94[ (group_info94$Group == "LTBI_Progressor" & group_info94$Timepoint == 0) | group_info94$Group == "LTBI_Non_Progressor" | group_info94$Group == "Control", ]
head(subset_colData)
subset_raw_ltbprog_baseline <- filtered_lc_data94[, rownames(subset_colData)]
dim(subset_raw_ltbprog_baseline)
subset_dds <- DESeqDataSetFromMatrix(countData = subset_raw_ltbprog_baseline,
                                     colData = subset_colData,
                                     design = ~ Group,
                                     tidy = FALSE,
                                     ignoreRank = FALSE,
                                     rownames(subset_raw_ltbprog_baseline))

subset_dds <- DESeq(subset_dds)
res_ltbprog_baseline_LTBI <- results(subset_dds, contrast = c("Group", "LTBI_Progressor", "LTBI_Non_Progressor"))
sig_genes_ltbprog <- rownames(res_ltbprog_baseline_LTBI)[
  res_ltbprog_baseline_LTBI$padj < 0.05 & abs(res_ltbprog_baseline_LTBI$log2FoldChange) > 1
]
sig_genes_active <- rownames(res_active_ltbinonprog)[
  res_active_ltbinonprog$padj < 0.05 & abs(res_active_ltbinonprog$log2FoldChange) > 1
]


### Create Upset plot to compare between LTBI Prog Vs LTBI Non prog and Active TB vs LTBI Non Prog
sig_list <- list(
  "ActiveTB_LTBINP" = sig_genes_active,
  "LTBIProg_LTBINP" = sig_genes_ltbprog
)
upset_data <- fromList(sig_list)
upset(upset_data, order.by = "freq")
