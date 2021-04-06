library(Seurat)
library(tidyverse)
library(readxl)
library(viridis)
library(clustree)

#load functions and colors
source(file.path("R", "functions.R"))

#load data
load(file.path("data", "seurat_all_without_doublets.RData"))

#subset for cpoglia cells
cp <- subset(all_singl, cells=grep("CP", colnames(all_singl)))

#re-run seurat
# run sctransform
cp <- SCTransform(cp, vars.to.regress = "percent.mt", verbose = FALSE)

#pca etc
cp <- RunPCA(cp, verbose = FALSE)

#run elbow plot to find most relevant PCs
ElbowPlot(cp)

#run UMAP and clustering on chosen PCs
cp <- RunUMAP(cp, dims = 1:15, verbose = FALSE)
cp <- FindNeighbors(cp, dims = 1:15, verbose = FALSE)
cp <- FindClusters(cp, resolution = seq(from=0.2,to=2,by=0.2))

#run clustree to find a cluster resolution
#url https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html#seurat-objects
clustree(cp)

#re cluster based on clustree output
cp <- FindClusters(cp, resolution = 1.2)

#plot clusters
DimPlot(cp, label = TRUE) + NoLegend()

#save data
save(cp, file = file.path("data", "seurat_cp.RData"))
