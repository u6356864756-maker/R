### Class 2 continuation
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# The following initializes usage of Bioc devel
#BiocManager::install("clusterProfiler")
#BiocManager::install("Mfuzz")
#BiocManager::install("ComplexHeatmap")

library(Seurat)
library(cowplot)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(pheatmap)
library(e1071)
library(Mfuzz)
library(clValid)
library(clusterProfiler)
# install dataset
ctrl.data <- read.table(file = "immune_control_expression_matrix.txt.gz", sep = "\t")
stim.data <- read.table(file = "immune_stimulated_expression_matrix.txt.gz", sep = "\t")
dim(ctrl.data)
dim(stim.data)

# Set up control object
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "IMMUNE_CTRL", 
                           min.cells = 5)
ctrl$stim <- "CTRL"
ctrl <- subset(ctrl, subset = nFeature_RNA > 500)
ctrl <- NormalizeData(ctrl, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)

# Set up stimulated object
stim <- CreateSeuratObject(counts = stim.data, project = "IMMUNE_STIM", min.cells = 5)
stim$stim <- "STIM"
stim <- subset(stim, subset = nFeature_RNA > 500)
stim <- NormalizeData(stim, verbose = FALSE)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000)

#calculate anchors 
immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim),
                                         dims = 1:20)
#compare the two sets 
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
saveRDS(immune.combined, "Immune.RDS")


##### Class March 24, first day of linux 
#we hard started to integrate two different conditions
#need to call the libraries and the data. Immune.RDS
readRDS("Immune.RDS")
DefaultAssay(immune.combined) <- "integrated"


immune.combined<-ScaleData(immune.combined, verbose = FALSE)
immune.combined<-RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined<-RunUMAP(immune.combined, reduction = 'pca', dims = 1:20)
immune.combined<-FindNeighbors(immune.combined, reduction = 'pca', dims = 1:20)
immune.combined<-FindClusters(immune.combined, resolution = 0.5)

#Visualization 
#Label by whether sample is immune stimulated 
p1<- DimPlot(immune.combined, reduction = 'umap', group.by = 'stim') 
p2<- DimPlot(immune.combined, redution = 'umap', label = TRUE)
plot_grid(p1, p2) 
#this is garbage, really bad seperation 
#ID by cluster and seperates the two conditions 
DimPlot(immune.combined, reduction = "umap", split.by = "stim") #split.by command is the important one 


#To identify canonical cell type marker genes that are conserved across conditions, we provide the `FindConservedMarkers()` function. 
#This function performs differential gene expression testing for each dataset/group and combines the p-values 
#using meta-analysis methods from the MetaDE R package. For example, we can calculated the 
#genes that are conserved markers irrespective of stimulation condition in cluster 6 (NK cells) 
#lookinlg for differentially expressed genes in the dataset 
install.packages('BiocManager') 
BiocManager::install('multtest')
install.packages('metap')
DefaultAssay(immune.combined) <-"RNA"

nk.markers <-FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "stim")
#find conserved genes across conditions. What genes are in both in CLUSTER 6, ident.1=6

#using the idents functions, we can change the cell labels to be the actual cell types 
#We have 14 x 2 cells, 14 of each type and 2 conditions 
Idents(immune.combined) <- factor(
  Idents(immune.combined),
  levels = c("HSPC", "Mono/Mk Doublets", "pDC", "Eryth","Mk", "DC", "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated", "CD4 Naive T", "CD4 Memory T"))
markers.to.plot <- c("CD3D","CREM","HSPH1","SELL","GIMAP5","CACYBP","GNLY","NKG7","CCL5","CD8A","MS4A1","CD79A","MIR155HG","NME1","FCGR3A","VMO1","CCL2","S100A9","HLA-DQA1","GPR183","PPBP","GNG11","HBA2","HBB","TSPAN13","IL3RA","IGJ","PRSS57")
DotPlot(immune.combined, features = markers.to.plot, cols = c(‘blue’, ‘red’), dot.scale = 8, split.by = "stim") + RotatedAxis()

