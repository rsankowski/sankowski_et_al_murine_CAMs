library(Seurat)
library(DoubletFinder)

#load data
load(file.path("data", "seurat_all_with_doublets.RData"))

#run doublet finder
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.025*nrow(all@meta.data))  ## Assuming 2.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
#low confidence
all <- doubletFinder_v3(all, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
#high confidence
all <- doubletFinder_v3(all, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_70", sct = TRUE)

#plot duplet prediction
all%>%
  DimPlot(group.by = "DF.classifications_0.25_0.09_70")
  

all_singl <- subset(all, cells = rownames(all@meta.data[all@meta.data$DF.classifications_0.25_0.09_70 == "Singlet",]))
save(all_singl, file = file.path("data", "seurat_all_without_doublets.RData"))
