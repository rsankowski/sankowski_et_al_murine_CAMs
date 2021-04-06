library(Seurat)
library(tidyverse)
library(readxl)
library(viridis)
library(clustree)

#load functions and colors
source(file.path("R", "functions.R"))

#load data
load(file.path("data", "seurat_all_without_doublets.RData"))

#subset for menoglia cells
men <- subset(all_singl, cells=grep("Men", colnames(all_singl)))

#re-run seurat
# run sctransform
men <- SCTransform(men, vars.to.regress = "percent.mt", verbose = FALSE)

#pca etc
men <- RunPCA(men, verbose = FALSE)

#run elbow plot to find most relevant PCs
ElbowPlot(men)

#run UMAP and clustering on chosen PCs
men <- RunUMAP(men, dims = 1:15, verbose = FALSE)
men <- FindNeighbors(men, dims = 1:15, verbose = FALSE)
men <- FindClusters(men, resolution = seq(from=0.2,to=2,by=0.2))

#run clustree to find a cluster resolution
#url https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html#seurat-objects
clustree(men)

#re cluster based on clustree output
men <- FindClusters(men, resolution = 1.2)

#plot clusters
DimPlot(men, label = TRUE) + NoLegend()

#save data
save(men, file = file.path("data", "seurat_men.RData"))
