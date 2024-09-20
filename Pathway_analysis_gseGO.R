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



#### Performing the pathway enrichment analysis using Gene Ontology (gseGO)
gse_ltb1 <- gseGO(geneList = gene_list_ltb1, 
                  ont = "BP",
                  keyType = "ENSEMBL",
                  minGSSize = 3,
                  maxGSSize = 200,
                  pvalueCutoff = 0.05,
                  verbose = TRUE,
                  OrgDb = organism,
                  pAdjustMethod = "BH")
gse_ltb1 <- clusterProfiler::simplify(gse_ltb1,
                                      by = "p.adjust",
                                      cutoff = 0.7,
                                      select_fun = min,
                                      measure = "Wang",
                                      semData = NULL)

gse_ltb2 <- gseGO(geneList = gene_list_ltb2, 
                  ont = "BP",
                  keyType = "ENSEMBL",
                  minGSSize = 3,
                  maxGSSize = 200,
                  pvalueCutoff = 0.05,
                  verbose = TRUE,
                  OrgDb = organism,
                  pAdjustMethod = "BH")
gse_ltb2 <- clusterProfiler::simplify(gse_ltb2,
                                      by = "p.adjust",
                                      cutoff = 0.7,
                                      select_fun = min,
                                      measure = "Wang",
                                      semData = NULL)

gse_ltb3 <- gseGO(geneList = gene_list_ltb3, 
                  ont = "BP",
                  keyType = "ENSEMBL",
                  minGSSize = 3,
                  maxGSSize = 200,
                  pvalueCutoff = 0.05,
                  verbose = TRUE,
                  OrgDb = organism,
                  pAdjustMethod = "BH")
gse_ltb3 <- clusterProfiler::simplify(gse_ltb3,
                                      by = "p.adjust",
                                      cutoff = 0.7,
                                      select_fun = min,
                                      measure = "Wang",
                                      semData = NULL)


gse_ltb4 <- gseGO(geneList = gene_list_ltb4, 
                  ont = "BP",
                  keyType = "ENSEMBL",
                  minGSSize = 5,
                  maxGSSize = 200,
                  pvalueCutoff = 0.05,
                  verbose = TRUE,
                  OrgDb = organism,
                  pAdjustMethod = "BH") 
gse_ltb4 <- clusterProfiler::simplify(gse_ltb4,
                                      by = "p.adjust",
                                      cutoff = 0.7,
                                      select_fun = min,
                                      measure = "Wang",
                                      semData = NULL)



### Plotting the results 
require(DOSE)
library(gridExtra)

g1 <-dotplot(gse_ltb4, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("LTB1 GO Pathways")
g2 <- dotplot(gse_ltb1, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("LTB2 GO Pathways")
g3 <- dotplot(gse_ltb2, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("LTB3 GO Pathways")
g4 <- dotplot(gse_ltb3, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("LTB4 GO Pathways")

combined_plot_g <- grid.arrange(g1, g2, g3, g4, nrow = 2)
ggsave("combined_dotplot_gseGOBP.png", plot = combined_plot_g, width = 16, height = 20)


