library(clusterProfiler)
library(enrichplot)
library(europepmc)
library(ggupset)
library(wordcloud)

library(ggplot2)
library(Seurat)
library(patchwork)
BiocManager::install("org.Mm.eg.db", character.only = TRUE)
library("org.Hs.eg.db")


#Question 1
df <- readRDS('exam1.rds')

#Making PErcent MT col 
df[['percent.mt']]<-PercentageFeatureSet(df, pattern = "^mt-")
VlnPlot(df, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(df, feature1 = "nCount_RNA", feature2 = "percent.mt")

#subsetting out dead cells 
df <-subset(df, subset = nFeature_RNA > 200 & nFeature_RNA < 15000 & percent.mt < 5)

#Question 2 
#normalize data 
#Normalizing the DATA 
df <- NormalizeData(df, normalization.method = "LogNormalize", scale.factor = 1e4)

#Question 3 
#Highly Variable Features: values with high varience between cells. usually select 200 genes 
df<- FindVariableFeatures(df, selection.method = "vst", nfeatures = 2000)
top10<- head(VariableFeatures(df))

#Step 4 Scaling (standarization) method 1 
all.genes <- rownames(df)
#list of all genes present 
df<-ScaleData(df, features = all.genes)

df <-RunPCA(df, features = VariableFeatures(df))
#print x # of PCA vectors, and names of # of genes/PCA vector 
print(df[["pca"]], dims = 1:5, nfeatures = 6)


df<-JackStraw(df, num.replicate = 100)
df <-ScoreJackStraw(df, dims = 1:20) #choose 20 vectors 
JackStrawPlot(df, dims = 1:20 ) #based on P value the best is PC13

library(reticulate)
library(umap)

#For plots to have clusters these commands are necessary first. Dont forget 
df <- FindNeighbors(df, dims = 1:10)
df <- FindClusters(df)
Idents(df)

df <-RunUMAP(df, dims = 1:13) #using PCA 10 vectors 


DimPlot(df, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(df, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)



df <- RunTSNE(df, seed.use = 1, dims = 1:10)
DimPlot(df, reduction = "tsne", label = T, repel = F)

#Question 4 
# the first step is already completed with the clustering. Now it is time to ID the clusters 
#this is done usign FindMarkers 

markers <- FindAllMarkers(df, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers)
library(tidyverse)
#markers is going to contain distinct markers for each group indicidually 
markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)


#By identifying ident.1 and iden.t2 we can either do 1 specific cluster to another. Or id ident.2 is blank it is all against 1 
cluster5.markers <- FindMarkers(df, ident.1 = 0, ident.2 = 2, min.pct = 0.25)
dim(cluster5.markers) # 1.3 k genes
head(cluster5.markers, n = 5)


#Plotting 
library(dplyr)

top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE)
top5_list <- split(top5$gene, top5$cluster)

FeaturePlot(df, features = top5_list[["0"]])
FeaturePlot(df, features = top5_list[["1"]])
FeaturePlot(df, features = top5_list[["2"]])
FeaturePlot(df, features = top5_list[["3"]])
#Automagic Method 
for (clust in names(top5_list)) {
  print(FeaturePlot(df, features = top5_list[[clust]]) + ggtitle(paste("Cluster", clust)))
}
VlnPlot(df, features = top5_list[["0"]])
VlnPlot(df, features = top5_list[["1"]])
VlnPlot(df, features = top5_list[["2"]])
VlnPlot(df, features = top5_list[["3"]])
for (clust in names(top5_list)) {
  print(VlnPlot(df, features = top5_list[[clust]]) + ggtitle(paste("Cluster", clust)))
}


#Question 6: GSEA 
#Step 1
markers_all <- FindMarkers(
  df,
  ident.1 = largest_cluster,
  ident.2 = NULL,   # compares vs all other cells
  logfc.threshold = 0,  # IMPORTANT: keep all genes
  min.pct = 0
)

gene_list <- markers_all$avg_log2FC
names(gene_list) <- rownames(markers_all)

# Only remove truly bad values
gene_list <- gene_list[is.finite(gene_list)]

# DO NOT threshold by logFC
gene_list <- sort(gene_list, decreasing = TRUE)

library(clusterProfiler)
library('org.Mm.eg.db')
head(names(gene_list), 20)

# Convert gene symbols to Entrez IDs make sure using correct species 
#GPT says the data is Mouse 

gene_df <- bitr(names(gene_list),
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)

gene_list_entrez <- gene_list[gene_df$SYMBOL]
names(gene_list_entrez) <- gene_df$ENTREZID

# Minimal cleaning only
gene_list_entrez <- gene_list_entrez[!duplicated(names(gene_list_entrez))]
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)
set.seed(123)
gene_list_entrez <- gene_list_entrez + rnorm(length(gene_list_entrez), 0, 1e-6)
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

#This is a check. Greater than 70% is a pass. we have 95 

#Run Bio Proces 


gsea_bp <- gseGO(
  geneList = gene_list_entrez,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  pvalueCutoff = 0.1,
  verbose = FALSE
)
length(gene_list_entrez)
summary(gene_list_entrez)

#molecular Functions 
gsea_mf <- gseGO(
  geneList = gene_list_entrez,
  OrgDb = org.Mm.eg.db,
  ont = "MF",
  keyType = "ENTREZID",
  pvalueCutoff = 0.1,
  verbose = FALSE
)

#Cellular Component 
gsea_cc <- gseGO(geneList = gene_list_entrez,
                 OrgDb = org.Mm.eg.db,
                 ont = "CC",
                 keyType = "ENTREZID",
                 verbose = FALSE, pvalueCutoff = 0.05)

#Visualize 
dotplot(gsea_bp, showCategory = 10)
dotplot(gsea_mf, showCategory = 10)
dotplot(gsea_cc, showCategory = 10) 

library(enrichplot)

cnetplot(
  gsea_bp,
  showCategory = 5,
  foldChange = gene_list_entrez
)


gsea_results <- list(
  markers_all = markers_all,
  gene_list = gene_list,
  gene_list_entrez = gene_list_entrez,
  gsea_bp = gsea_bp,
  gsea_mf = gsea_mf, 
  gsea_cc = gsea_cc
)

other_results <- list(df, gene_df, genes_df, markers, markers_all, top5, top5_list)

saveRDS(other_results, file = 'other.rds')
saveRDS(gsea_results, file = "gsea_results.rds")

