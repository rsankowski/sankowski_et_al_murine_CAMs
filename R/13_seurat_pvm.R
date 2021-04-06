library(Seurat)
library(tidyverse)
library(readxl)
library(viridis)
library(clustree)

#load functions and colors
source(file.path("R", "functions.R"))

#load data
load(file.path("data", "seurat_all_without_doublets.RData"))

#subset for pvmoglia cells
pvm <- subset(all_singl, cells=grep("PVM", colnames(all_singl)))

#re-run seurat
# run sctransform
pvm <- SCTransform(pvm, vars.to.regress = "percent.mt", verbose = FALSE)

#pca etc
pvm <- RunPCA(pvm, verbose = FALSE)

#run elbow plot to find most relevant PCs
ElbowPlot(pvm)

#run UMAP and clustering on chosen PCs
pvm <- RunUMAP(pvm, dims = 1:15, verbose = FALSE)
pvm <- FindNeighbors(pvm, dims = 1:15, verbose = FALSE)
pvm <- FindClusters(pvm, resolution = seq(from=0.2,to=2,by=0.2))

#run clustree to find a cluster resolution
#url https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html#seurat-objects
clustree(pvm)

#re cluster based on clustree output
pvm <- FindClusters(pvm, resolution = 1.2)

#plot clusters
DimPlot(pvm, label = TRUE) + NoLegend()

#save data
save(pvm, file = file.path("data", "seurat_pvm.RData"))
