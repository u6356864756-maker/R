#### Irina S. Moreira
### System Biology class 1
#####

#install.packages("ggplot2")
#install.packages("patchwork")
#install.packages("Seurat")
#install.packages("Seurat", dependencies = T)
# Install the remotes package
#if (!requireNamespace("remotes", quietly = TRUE)) {
#  install.packages("remotes")
#}
#install.packages('Signac')
#remotes::install_github("satijalab/seurat-data", quiet = TRUE)
#remotes::install_github("satijalab/azimuth", quiet = TRUE)
#remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)
library(Seurat)
library(ggplot2)
library(patchwork)

# sets the working directory or
# "Session" on the top, and choose "Set working directory" 
# to whatever you want
setwd("~/Desktop/Irina_UC/Professora_FCTUC/MBC/2025-2026-Edicao6/BSM")

?Read10X # ?keyword opens a help menu

BaseDir <- getwd() # Create a variable called BaseDir that provides your directory
BaseDir
#is the path to your working directory
data_dir <- file.path(BaseDir, "hg19") #hg19 was the provided zip file
data_dir
pbmc.data <- Read10X(data.dir = data_dir)

# cells should have more than 200 genes and each gene should be expressed
# in at least 3 cells
pbmc <- CreateSeuratObject(
  counts = pbmc.data,
  project = "pbmc3k",
  min.cells = 3,
  min.features = 200
)

pbmc # 13714 genes and 2700 cells; RNA
?CreateSeuratObject # if you need help to know how to use command
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

pbmc.data@Dimnames # @ sign allows you to access a few number of different
#data

1:3 # vector with positions 1, 2 and 3
c(2, 5, 7) # vector
# before comma we have rows and, after comma we have columns
# i am selecting the last ten columns
pbmc.data[, 1:10]  # I am selecting the first ten columns
pbmc.data[1:10, ] # I am selecting the first ten rows
pbmc.data[1:3, c(2, 4)] # select the first 3 rows and columns 2 and 4
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:5] # select the 3 different genes in the 
#first 5 columns/cells

# we are back to the object seurat = pbmc

### STEP 1 - QUALITY CONTROL 
# we are going to create in this object a new variable
# called "percent.mt"
# for that I use keyword PercentageFeatureSet
# calculate the number of time that something occurs
#pattern = "^MT-" we are going to calculate the % of gene that have
# the pattern/the symbol MT in its name
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc@meta.data
#head(dataset, 8) this selects the first 8 rows of a dataset called dataset
head(pbmc@meta.data, 5) # using head(, 5) I am selecting the first 5 rows of my
# dataset

# Violin plot (VlnPlot is the keyword
# features = to decide which columns you want to plot
# ncol = 3
p1 = VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # side by side ; 3 columns
p1
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1) # on top of each other; 1 column

 # instead of plotting directly to the screen, I am creating/saving a variable called p1 and p2
# where I have my plots
# FeatureScatter is the Seurat keyword to save scatter plots
# feature1 = input the variable that you want for x axis
# feature2 = input the variable that you want for y axis
p1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1 # just calls plot p1
p1 + p2 # this calls both plots at the same time

# From graphical analysis we take out outliers
# & is used to impose that all conditions are fulfilled at the same time
# & is the same AND
# subset(seuratobject, subset = condition1 & condition2)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc
#with QC we took out 62 problematic cells/samples

## STEP 2 - Normalize data
?NormalizeData # *10000 and then make logarithm 
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize")

## STEP 3- FIND HVF
?FindVariableFeatures
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000) # 2000 the typical valuable
# it can be played around

# we are creating a variable called "top10" with the first 10 HVF
top10 <- head(VariableFeatures(pbmc), 10)
top10
p1 <-VariableFeaturePlot(pbmc)
p2 <- LabelPoints(plot = p1, points = top10)
#The LabelPoints() function is particularly useful when you want to highlight specific points or annotate them with labels in a plot. By specifying the points argument as top10, you are indicating that you want to label the top 10 points on the plot. When repel = TRUE, the geom_text_repel function is used to create nicely repelled labels, which helps prevent label overlap. However, it's important to note that using geom_text_repel can be slow when plotting a large number of points
p1
p2
p1 + p2
# p2 is a more interesting plot

# STEP 4 - Scale
all.genes <- rownames(pbmc) # a way to get the names of the rows/in
# this case it is a list of genes
all.genes # call the list of genes
length(all.genes) # provides the length of the vector with the list of genes

pbmc <- ScaleData(pbmc, features = all.genes)

# Step 5- PCA for dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))

# printing x principal components/vector 1:5 you decide that you wanted
# the first 5 and with nfeatures = number of genes to print
print(pbmc[["pca"]], dims = 1:5, nfeatures = 6)

# the keyword for the Seurat plot of PCA
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# Scatter plot for the PCA reduction
DimPlot(pbmc, reduction = "pca") # are very useful
# for pattern recognition, to be seen in a moment

#Heatmaps
#heatmap for PC= 1 (dims =1)
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

## we have to decide of vectors to keep
# they should be as low as possible while still representing
# the majority of the variability of the data

# METHOD 1- simpler, quicker but less reliable
ElbowPlot(pbmc)
# from graphical analysis choose the PCs after which you reach 
# a plateau, adding more PCs does not increase the explainability
# In this case probably around 10

# Method 2 - more reliable but computationally demanding
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20) # 20 vectors are chosen
JackStrawPlot(pbmc, dims = 1:20) #plot and choose based on p-value
# that appears on the right
# 13 PCs

# ANOTHER METHOD FOR REDUCTION OF DIMENSIONALITY
# UMAP AND tSNE
#install.packages("reticulate")
#reticulate::py_install(packages = "umap-learn")
pbmc <- RunUMAP(pbmc, dims = 1:10) # I decide to use 10 as provided by ELbowMethod 
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE) 

# Another - third and last method to perform reduction of dimensionality is tSNE
# is a probablistic method, so have to set seed to whatever you want
pbmc <- RunTSNE(pbmc, seed.use = 1,  dims = 1:10)
DimPlot(pbmc, reduction = "tsne", label = TRUE, repel = TRUE) 


p10 <- DimPlot(pbmc, reduction = "pca", label = TRUE, repel = TRUE)
p11 <- DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE)
p12 <- DimPlot(pbmc, reduction = "tsne", label = TRUE, repel = TRUE)
p10 + p11 
p11 + p12

saveRDS(pbmc, file = "pbmc3k_final.rds")

pbmc
pbmc <- FindNeighbors(pbmc, dims = 1:10)
# By specifying dims = 1:10, the code performs the SNN graph construction 
#using the first 10 principle components (PCs) obtained from a dimensionality 
#reduction step, such as PCA. Computing the SNN graph based on a dimensionally reduced form of the data, like the first few PCs, is a common approach to reduce computation time and focus on the most informative features of the dataset.

pbmc <- FindClusters(pbmc, resolution = 0.5)
# By specifying resolution = 0.5, the code sets the granularity or sensitivity of the clustering analysis. A higher resolution value leads to a greater number of clusters, while a lower resolution value results in fewer clusters but with higher cell-to-cell similarity within each cluster
pbmc
Idents(pbmc) # identity default is the cluster
head(Idents(pbmc), 10)

# This is to compare one cluster in relation all the others
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, 5)

# Ro make the calculation for all cluster at the same time:
pbmc.markers <- FindAllMarkers(pbmc, min.pct = 0.25)
head(pbmc.markers, n = 5)
?FindAllMarkers
# object = seurat_object == pbmc
# features = you can choose the genes of interest
# group.by = default is you "ident" cluster
# logfc.threshold = deafult is 0.1
#test.use = wilcox is the default (MAST, t.test, roc)
#min.pct = minimum expression in the different groups
#min.diff.pct = minimum expression difference
# only.pos = positive markers (default is F)
dim(pbmc.markers) # almost 6k genes DE is the cells
library(tidyverse)
head(pbmc.markers)
pbmc.markers %>% group_by(cluster) %>% 
  slice_max(n = 2, order_by = avg_log2FC) %>%
  as.data.frame()

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers_2 <- FindAllMarkers(pbmc, only.pos = TRUE, 
                                 min.pct = 0.25, logfc.threshold = 0.25)
dim(pbmc.markers_2) # 3.5 k genes instead of almost 6k

# find all markers distinguishing cluster 5 from clusters 0 
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = 0, min.pct = 0.25)
dim(cluster5.markers) # 1.3 k genes
head(cluster5.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0  and 4
cluster5.markers_2 <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0,4), min.pct = 0.25)
dim(cluster5.markers_2) # 1.3 k genes
head(cluster5.markers_2, n = 5)


cluster0.markers_1 <- FindMarkers(pbmc, ident.1 = 0, 
                                  logfc.threshold = 0.25, 
                                  only.pos = TRUE)
dim(cluster0.markers_1)
cluster0.markers_2 <- FindMarkers(pbmc, ident.1 = 0, 
                                logfc.threshold = 0.25, test.use = "roc", 
                                only.pos = TRUE)
dim(cluster0.markers_2)
# ident.1 = 0 specifies the first group or cluster to compare. In this case, it refers to Cluster 0.
# logfc.threshold = 0.25 sets the log-fold change threshold, indicating the 
#minimum magnitude of differential expression to consider a gene as a marker.
# test.use = "roc" specifies the statistical test to be used for differential expression analysis. In this case, the ROC test is used.
# only.pos = TRUE indicates that only positive markers (upregulated in Cluster 0) will be returned.

#Different plots answer different questions:
#  - **VlnPlot**: distribution of expression per cluster  
#- **FeaturePlot**: spatial pattern over UMAP  
#- **DotPlot**: compact view of many markers across clusters  
#- **Heatmap**: top markers per cluster


# VlnPlot() is a function in Seurat that generates violin plots to visualize the expression 
#probability distributions across clusters or groups.
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# features = c("MS4A1", "CD79A") specifies the features (genes) of interest to be plotted. In this case, the expression of two genes, MS4A1 and CD79A, will be visualized.
# The resulting violin plots will show the expression distribution of the selected genes (MS4A1 and CD79A) across clusters or groups within the Seurat object. These plots can provide insights into the expression patterns and heterogeneity of the chosen genes across different cell populations.

?VlnPlot
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = 'counts', log = TRUE)
# features = c("NKG7", "PF4") specifies the features (genes) of interest to be plotted. In this case, the expression of two genes, NKG7 and PF4, will be visualized.
# slot = 'counts' indicates that the violin plots will be based on the raw count data stored in the 'counts' slot of the Seurat object.
# log = TRUE specifies that the expression values will be log-transformed before plotting.
# The resulting violin plots will show the expression distribution of the selected genes (NKG7 and PF4) in the raw count data of the Seurat object. By setting slot = 'counts' and log = TRUE, the plots will visualize the log-transformed expression values from the raw counts.

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A",
                               "PPBP"), ncol = 2)
# The resulting feature plots will show the expression patterns of the selected genes (MS4A1, GNLY, CD3E, CD14, FCER1A, FCGR3A, LYZ, PPBP, CD8A) on a dimensional reduction plot, allowing for the visualization of their expression patterns across the cells in the Seurat object.


### DotPlot with canonical PBMC markers 

DotPlot(
  pbmc,
  features = c("IL7R","CCR7","CD3D","CD3E","CD14","LYZ","MS4A1","CD79A","NKG7","GNLY","FCER1A","CST3","PPBP")
) + RotatedAxis()


### Heatmap of top markers per cluster
top5 <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%as.data.frame()
dim(top5)
head(top5)
top5$gene
DoHeatmap(pbmc, features = top5$gene) # without legend
DoHeatmap(pbmc, features = top5$gene) + NoLegend() # with legend

# Assign labels to different clusters according to our canonical markers' cell types
# create a vector with the names of the 9 different types of cels
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T",
                     "FCGR3A+ Mono", "NK", "DC", "Platelet")
levels(pbmc)
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
levels(pbmc)
DimPlot(pbmc, reduction = 'umap', label = TRUE, pt.size = 0.5) + NoLegend()



# Let's create a UMAP plot of the cells in the pbmc object, with customized axis labels, 
#text sizes, and legend appearance.
plot1 <-DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + 
  xlab("UMAP 1") + ylab("UMAP 2") + 
  theme(axis.title = element_text(size = 18),
        legend.text = element_text(size = 18)) + 
  guides(colour = guide_legend(override.aes = list(size = 10)))

ggsave(filename = "pbmc3k_umap.jpg", height = 7, width = 12, plot = plot1,
       quality = 50)


##### EXERCISE
ex1 <- readRDS("./Ex1.RDS")
#Create a metadata column indicating the percentage of mitochondrial genes. 
ex1 # 5070 cells and 17k genes

ex1[["percent.mt"]] <- PercentageFeatureSet(ex1, pattern = "^MT-")

plot1 <- FeatureScatter(ex1, feature1 = "nCount_RNA", 
                        feature2 = "percent.mt") + geom_hline(yintercept = 25)
plot1
# The code is used to create a scatter plot to visualize the relationship between the number of RNA counts
#(nCount_RNA) and the percentage of mitochondrial genes (percent.mt) in the Seurat object pbmc

plot2 <- FeatureScatter(ex1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# The code is used to create a scatter plot to visualize the relationship between the number of RNA counts (nCount_RNA)
#and the number of RNA features (nFeature_RNA) in the Seurat object pbmc

plot1 + plot2

VlnPlot(ex1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# A good filter might be around 25
ex1 <- subset(ex1, subset = percent.mt < 25)

#Mitochondrial genes can be markers for apoptosis or other irregulatiries where mitochondral DNA has escaped the organelle, 
#so it is a good marker for contamination

#b.	How can you normalize the data and remove genes that are likely uninformative for the overall structure? 
# Generate the necessary calculation and graphs to ensure quality control for each object. If required, filter out cells.

FeatureScatter(ex1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_vline(xintercept = 1e4) + 
  geom_hline(yintercept = 200) + 
  geom_hline(yintercept = 2700)

# There are few cells lower that the lower bound, so this removal might be necessary. Upper bounds also seem reasonable
ex1 <- subset(ex1, subset = nFeature_RNA > 200 & nFeature_RNA < 2700)
# we lost 34 cells of low quality
34/5070 *100
# 0.7 % of low quality cells were taken out
ex1
ex1 <- NormalizeData(ex1, normalization.method = "LogNormalize", scale.factor = 1e4)
ex1 <- FindVariableFeatures(ex1, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(ex1) %>% 
  LabelPoints(points = head(VariableFeatures(ex1), 10), repel = TRUE)

#c.	Could you scale up the data and perform dimensionality reduction? 
# Scale the data, perform PCA to estimate the optimal number of PCs that would best represent the data. Afet performing clustering, produce t-SNE and UMAP plots.

all.genes <- rownames(ex1)
ex1 <- ScaleData(ex1, features = all.genes)

ex1 <- RunPCA(ex1, features = VariableFeatures(object = ex1))
VizDimLoadings(ex1, dims = 1:4, reduction = "pca")
DimPlot(ex1, reduction = "pca")
ElbowPlot(ex1)

ex1 <- JackStraw(ex1, num.replicate = 100)
ex1 <- ScoreJackStraw(ex1, dims = 1:20)
JackStrawPlot(ex1, dims = 1:20) 

ex1 <- FindNeighbors(ex1, dims = 1:20)
ex1 <- FindClusters(ex1, resolution = 0.5)

ex1 <- RunTSNE(ex1, dims = 1:20)
p1 <-DimPlot(ex1, reduction = "tsne", label = T)

ex1 <- RunUMAP(ex1, dims = 1:20)
p2<- DimPlot(ex1, reduction = "umap", label = T)
p3<- DimPlot(ex1, reduction = "pca", label = T)
p1+p2+p3

#d.	Is it possible to define cell clusters based on the processed data? 

# Conduct a clustering analysis of cells. Interpret the clustering using initial differential gene analysis

## Individually
cluster0.markers <- FindMarkers(ex1, ident.1 = 0, min.pct = 0.25)
cluster0.markers %>% 
  ggplot(aes(x = avg_log2FC, y = -log(p_val_adj), col = ifelse(p_val_adj < 0.05, 1, 0))) + 
  geom_point()
VlnPlot(ex1, features = rownames(cluster0.markers)[1:2])

cluster1.markers <- FindMarkers(ex1, ident.1 = 1, min.pct = 0.25)
cluster1.markers %>% ggplot(aes(x = avg_log2FC, y = -log(p_val_adj), 
                                col = ifelse(p_val_adj < 0.05, 1, 0))) + geom_point()
VlnPlot(ex1, features = rownames(cluster1.markers)[1:2])

cluster2.markers <- FindMarkers(ex1, ident.1 = 2, min.pct = 0.25)
cluster2.markers %>% 
  ggplot(aes(x = avg_log2FC, y = -log(p_val_adj), col = ifelse(p_val_adj < 0.05, 1, 0))) + 
  geom_point()
VlnPlot(ex1, features = rownames(cluster2.markers)[1:2])

cluster3.markers <- FindMarkers(ex1, ident.1 = 3, min.pct = 0.25)
cluster3.markers %>% ggplot(aes(x = avg_log2FC, y = -log(p_val_adj), 
                                col = ifelse(p_val_adj < 0.05, 1, 0))) + geom_point()
VlnPlot(ex1, features = rownames(cluster3.markers)[1:2])

## Grouped
ex1.markers <- FindAllMarkers(ex1, only.pos = TRUE, min.pct = 0.25, 
                              logfc.threshold = 0.25)
ex1.markers %>% ggplot(aes(x = avg_log2FC, y = -log(p_val_adj), col = ifelse(p_val_adj < 0.05, 1, 0))) + geom_point()

group_markers <- c()
group_markers <- c(group_markers,ex1.markers$gene[ex1.markers$cluster==0] %>% head(1))
group_markers <- c(group_markers,ex1.markers$gene[ex1.markers$cluster==1] %>% head(1))
group_markers <- c(group_markers,ex1.markers$gene[ex1.markers$cluster==2] %>% head(1))
group_markers <- c(group_markers,ex1.markers$gene[ex1.markers$cluster==3] %>% head(1))


VlnPlot(ex1, features = group_markers)

ex1.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5

DoHeatmap(ex1, features = top5$gene) + NoLegend()

# Display the top marker gene of each cluster in FeaturePlot and RidgePlot
ex1.markers %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = avg_log2FC) -> top1

FeaturePlot(ex1, features = top1$gene)
RidgePlot(ex1, features = top1$gene[1:2])
RidgePlot(ex1, features = top1$gene[3:4])

#e.	Can you identify the markers for cluster 0 using at least two statistical tests?
  
# Install MAST, in case you do not have it
#install.packages("BiocManager") # Needed to install all Bioconductor packages
#BiocManager::install("MAST")

#?FindMarkers
cluster0.markers_1 <- FindMarkers(ex1, ident.1 = 0, 
                                  logfc.threshold = 0.25,
                                  test.use = "roc", only.pos = TRUE)
cluster0.markers_2 <- FindMarkers(ex1, ident.1 = 0, 
                                  logfc.threshold = 0.25, 
                                  test.use = "MAST", only.pos = TRUE)

head(cluster0.markers_1); head(cluster0.markers_2)




