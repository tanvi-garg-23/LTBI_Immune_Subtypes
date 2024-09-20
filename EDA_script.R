####### Exploratory Data Analysis of Gene Counts Matrix

### Inputs required: raw matrix and normalized matrix and the group information containing dataframe all in rds format

### Load libraries
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(pheatmap)
library(VennDiagram)
library(UpSetR)


### Load datasets: raw gene counts matrix and normalized gene counts matrix
raw_matrix94 <- readRDS("rawcounts_proteinfiltered94.rds")
dim(raw_matrix94)

norm_matrix94 <- readRDS("counts_protein_exprs94.rds")
dim(norm_matrix94)

group_info94 <- readRDS("Group_information.rds")


### Plotting datasets distributions using histogram
par(mfrow = c(1, 2))
hist(as.matrix(raw_matrix94), main = "Raw Data",
     xlab = "Raw Counts Values", col = "lightgreen")
hist(as.matrix(norm_matrix94), main = "TMM-Normalized Data", 
     xlab = "Normalized values", col = "lightpink")


######### Unsupervised Clustering 

### 1. Heatmap using normalized gene counts matrix
hc_row <- hclust(as.dist(1-cor(t(norm_matrix94))),method="ward.D2") # genes/rows by correlation
hc_col <- hclust(dist(t(norm_matrix94)),method="ward.D2")   ### samples/cols by euclidean distance
annot <- data.frame(Group = as.factor(group_info94$Group)) ### setting annotation for coloring wrt to the groups 
rownames(annot) <- colnames(norm_matrix94)
unique(group_info94$Group)
# Define the colors for different disease groups
group_colors <- c(
  "Active_TB" = "#F8766D",     # Light red/pink
  "Control" = "#00BA38",       # Green
  "LTBI_Non_Progressor" = "#00BFC4",  # Cyan/Light blue
  "LTBI_Progressor" = "#C77CFF"       # Light purple
)
annotCol <- list(Group = group_colors)
colGradient <- colorRampPalette(c("blue", "white", "red"))
heatmaptitle <- "Heatmap: Normalized Gene Counts"
pheatmap(norm_matrix94,
         main = heatmaptitle,
         color = colGradient(40),
         annotation_col = annot,
         annotation_colors = annotCol,
         cluster_rows = hc_row,
         cluster_cols = hc_col,
         show_rownames = FALSE,
         show_colnames = FALSE,
         scale = "row",
         use_raster = FALSE)


### 2. PCA using normalized gene counts matrix
pca_res94 <- prcomp(t(norm_matrix94), scale. = TRUE)
pca_data94 <- as.data.frame(pca_res94$x)
pca_data94$label <- group_info94$Group
ggplot(pca_data94, aes(x = PC1, y = PC2, color = label)) +
  geom_point() + 
  ggtitle("PCA: Normalized Gene Counts") +
  xlab(paste("PC1 - ", round(summary(pca_res94)$importance[2, 1] * 100, 2), "% variance", sep="")) +
  ylab(paste("PC2 - ", round(summary(pca_res94)$importance[2, 2] * 100, 2), "% variance", sep=""))

