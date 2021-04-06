library(Seurat)
library(tidyverse)
library(readxl)
library(viridis)
library(clustree)

#load functions and colors
source(file.path("R", "functions.R"))

#load data
load(file.path("data", "seurat_all_without_doublets.RData"))

#subset for microglia cells
micr <- subset(all_singl, cells=grep("Micr", colnames(all_singl)))

#re-run seurat
# run sctransform
micr <- SCTransform(micr, vars.to.regress = "percent.mt", verbose = FALSE)

#pca etc
micr <- RunPCA(micr, verbose = FALSE)

#run elbow plot to find most relevant PCs
ElbowPlot(micr)

#run UMAP and clustering on chosen PCs
micr <- RunUMAP(micr, dims = 1:15, verbose = FALSE)
micr <- FindNeighbors(micr, dims = 1:15, verbose = FALSE)
micr <- FindClusters(micr, resolution = seq(from=0.2,to=2,by=0.2))

#run clustree to find a cluster resolution
#url https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html#seurat-objects
clustree(micr)

#re cluster based on clustree output
micr <- FindClusters(micr, resolution = 1.2)

#plot clusters
DimPlot(micr, label = TRUE) + NoLegend()

#save data
save(micr, file = file.path("data", "seurat_micr.RData"))
