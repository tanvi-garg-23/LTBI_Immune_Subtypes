#### Pathway Enrichment Analysis using clusterProfiler package: 
## Created using this: https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/


### Input required: DGE results for all comparisons 

### Load libraries
# BiocManager::install("clusterProfiler")
# BiocManager::install("pathview")
# BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
organism = "org.Hs.eg.db"
# BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)


### Load the input data 
df_ltb1 <- readRDS("ltb1_Control.rds")
head(df_ltb1)
df_ltb2 <- readRDS("ltb2_Control.rds")
head(df_ltb2)
df_ltb3 <- readRDS("ltb3_Control.rds")
head(df_ltb3)
df_ltb4 <- readRDS("ltb4_Control.rds")
head(df_ltb4)


#### Create gene list for all the comparisons

df_ltb1$x <- rownames(df_ltb1)
logfc_df_ltb1 <- df_ltb1$log2FoldChange
names(logfc_df_ltb1) <- df_ltb1$x
gene_list_ltb1 <- na.omit(logfc_df_ltb1) # remove na
length(gene_list_ltb1)
gene_list_ltb1 <- sort(gene_list_ltb1, decreasing = TRUE) # order in descending order 
head(gene_list_ltb1)

df_ltb2$x <- rownames(df_ltb2)
logfc_df_ltb2 <- df_ltb2$log2FoldChange
names(logfc_df_ltb2) <- df_ltb2$x
gene_list_ltb2 <- na.omit(logfc_df_ltb2)
length(gene_list_ltb2)
gene_list_ltb2 <- sort(gene_list_ltb2, decreasing = TRUE)
head(gene_list_ltb2)

df_ltb3$x <- rownames(df_ltb3)
logfc_df_ltb3 <- df_ltb3$log2FoldChange
names(logfc_df_ltb3) <- df_ltb3$x
gene_list_ltb3 <- na.omit(logfc_df_ltb3)
length(gene_list_ltb3)
gene_list_ltb3 <- sort(gene_list_ltb3, decreasing = TRUE)
head(gene_list_ltb3)

df_ltb4$x <- rownames(df_ltb4)
logfc_df_ltb4 <- df_ltb4$log2FoldChange
names(logfc_df_ltb4) <- df_ltb4$x
gene_list_ltb4 <- na.omit(logfc_df_ltb4)
length(gene_list_ltb4)
gene_list_ltb4 <- sort(gene_list_ltb4, decreasing = TRUE)
head(gene_list_ltb4)


#### Visualize the gene list distributions using histogram
par(mfrow = c(2, 2))
hist(gene_list_ltb1, breaks = 50, main = "Distribution of Gene List Values LTB1", xlab = "Gene Expression Values")
hist(gene_list_ltb2, breaks = 50, main = "Distribution of Gene List Values LTB2", xlab = "Gene Expression Values")
hist(gene_list_ltb3, breaks = 50, main = "Distribution of Gene List Values LTB3", xlab = "Gene Expression Values")
hist(gene_list_ltb4, breaks = 50, main = "Distribution of Gene List Values LTB4", xlab = "Gene Expression Values")




### KEGG GSEA
#### gseKEGG
library(AnnotationDbi)

ensembl_gene_ids <- names(gene_list_ltb1)
gene_info <- select(org.Hs.eg.db, keys = ensembl_gene_ids, keytype = "ENSEMBL", columns = c("ENSEMBL", "ENTREZID"))
print(gene_info)
gene_info <- gene_info[!is.na(gene_info$ENTREZID), ]
gene_list_entrez <- gene_list_ltb1[gene_info$ENSEMBL]
names(gene_list_entrez) <- gene_info$ENTREZID
head(gene_list_entrez)

gse_ltb1_KEGG <- gseKEGG(gene_list_entrez,
                         organism = "hsa",
                         keyType = "kegg",
                         exponent = 1,
                         minGSSize = 10,
                         maxGSSize = 200,
                         eps = 1e-10,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH",
                         verbose = TRUE,
                         by = "fgsea")


ensembl_gene_ids <- names(gene_list_ltb2)
gene_info <- select(org.Hs.eg.db, keys = ensembl_gene_ids, keytype = "ENSEMBL", columns = c("ENSEMBL", "ENTREZID"))
print(gene_info)
gene_info <- gene_info[!is.na(gene_info$ENTREZID), ]
gene_list_entrez <- gene_list_ltb2[gene_info$ENSEMBL]
names(gene_list_entrez) <- gene_info$ENTREZID
head(gene_list_entrez)

gse_ltb2_KEGG <- gseKEGG(gene_list_entrez,
                         organism = "hsa",
                         keyType = "kegg",
                         exponent = 1,
                         minGSSize = 10,
                         maxGSSize = 200,
                         eps = 1e-10,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH",
                         verbose = TRUE,
                         by = "fgsea")


ensembl_gene_ids <- names(gene_list_ltb3)
gene_info <- select(org.Hs.eg.db, keys = ensembl_gene_ids, keytype = "ENSEMBL", columns = c("ENSEMBL", "ENTREZID"))
print(gene_info)
gene_info <- gene_info[!is.na(gene_info$ENTREZID), ]
gene_list_entrez <- gene_list_ltb3[gene_info$ENSEMBL]
names(gene_list_entrez) <- gene_info$ENTREZID
head(gene_list_entrez)

gse_ltb3_KEGG <- gseKEGG(gene_list_entrez,
                         organism = "hsa",
                         keyType = "kegg",
                         exponent = 1,
                         minGSSize = 10,
                         maxGSSize = 200,
                         eps = 1e-10,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH",
                         verbose = TRUE,
                         by = "fgsea")



ensembl_gene_ids <- names(gene_list_ltb4)
gene_info <- select(org.Hs.eg.db, keys = ensembl_gene_ids, keytype = "ENSEMBL", columns = c("ENSEMBL", "ENTREZID"))
print(gene_info)
gene_info <- gene_info[!is.na(gene_info$ENTREZID), ]
gene_list_entrez <- gene_list_ltb4[gene_info$ENSEMBL]
names(gene_list_entrez) <- gene_info$ENTREZID
head(gene_list_entrez)

gse_ltb4_KEGG <- gseKEGG(gene_list_entrez,
                         organism = "hsa",
                         keyType = "kegg",
                         exponent = 1,
                         minGSSize = 10,
                         maxGSSize = 200,
                         eps = 1e-10,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH",
                         verbose = TRUE,
                         by = "fgsea")



### Plotting the results
require(DOSE)
library(gridExtra)
p1 <- dotplot(gse_ltb1_KEGG, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("LTB1 KEGG Pathways")
p2 <- dotplot(gse_ltb2_KEGG, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("LTB2 KEGG Pathways")
p3 <- dotplot(gse_ltb3_KEGG, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("LTB3 KEGG Pathways")
p4 <- dotplot(gse_ltb4_KEGG, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("LTB4 KEGG Pathways")

combined_plot <- grid.arrange(p1, p2, p3, p4, nrow = 2)
ggsave("combined_dotplot_gseKEGG.png", plot = combined_plot, width = 16, height = 20)



### enrichgo analysis
ora_result_ltb4 <- enrichGO(
  gene = names(gene_list_ltb4)[abs(gene_list_ltb4) > 1],  # Select significant genes
  OrgDb = organism,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

ora_result_ltb1 <- enrichGO(
  gene = names(gene_list_ltb1)[abs(gene_list_ltb1) > 1],  # Select significant genes
  OrgDb = organism,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

ora_result_ltb2 <- enrichGO(
  gene = names(gene_list_ltb2)[abs(gene_list_ltb2) > 1],  # Select significant genes
  OrgDb = organism,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

ora_result_ltb3 <- enrichGO(
  gene = names(gene_list_ltb3)[abs(gene_list_ltb3) > 1],  # Select significant genes
  OrgDb = organism,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)



