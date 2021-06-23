library(Seurat)
library(tidyverse)
library(readxl)
library(viridis)

#create plotting subfolders
dir.create(file.path("plots", "umap", "all"))
dir.create(file.path("plots", "heatmaps", "all"))
dir.create(file.path("plots", "others", "all"))

#load functions and colors
source(file.path("R", "functions.R"))

#load data
load(file.path("data", "seurat_all_without_doublets.RData"))

#set cluster order
order_clusters <- data.frame(seurat_clusters= all_singl@meta.data[,"seurat_clusters"], row.names = rownames(all_singl@meta.data)) %>%
  bind_cols(as.data.frame(t(all_singl[["SCT"]]@scale.data))) %>%
  group_by(seurat_clusters) %>%
  summarize_all(.funs=mean) %>%
  as.data.frame()
  
rownames(order_clusters) <- order_clusters$seurat_clusters
order_clusters <- order_clusters$seurat_clusters[hclust(dist(order_clusters[,-1]))$order]

#reorder the clusters
levels(all_singl) <- order_clusters[c(10:8,15, 7, 13:14, 16:17, 12, 11, 3:6, 1:2)]

#extract metadata
metadata <- all_singl@meta.data
metadata$seurat_clusters <- factor(metadata$seurat_clusters, levels = levels(all_singl))

#export metadata
write.csv(metadata, file.path("data", "all_cell_metadata.csv"))
tools::md5sum(file.path("data", "all_cell_metadata.csv"))

# fig 1a - conditions umap
all_singl %>%
DimPlot(group.by = "Condition", pt.size = 3) +
  theme_void() +
  scale_color_manual(values=colors_many[3:2])
ggsave(file.path("plots", "umap", "all", "conditions_all_umap.pdf"), useDingbats=F)

# fig 1b - cell signature umaps
signature_genes <-  read_excel(file.path("data","cell_signatures.xlsx"), "Core signature", skip = 2)

for (i in colnames(signature_genes)) {
  plt <- plot_expmap_seurat(na.omit(signature_genes[[i]]), object = all_singl, point_size = 3, .retain_cl = levels(all_singl)) + labs(subtitle=paste0(as.character(i), ' Signature'))
  print(plt)
  ggsave(file.path("plots", "umap", "all", paste0(i,"_signature_all.pdf")), useDingbats=F)
} 

# fig 1c - compartments umap
all_singl %>%
  DimPlot(group.by = "Compartment", pt.size = 3) +
  theme_void() +
  scale_color_brewer(palette="Dark2", direction = -1)

ggsave(file.path("plots", "umap", "all", "compartments_all_umap.pdf"), useDingbats=F)

# fig 1d - clusters umap
all_singl %>%
  DimPlot(label=T, pt.size = 3) +
  theme_void() +
  scale_color_manual(values=colors_many[-14])
ggsave(file.path("plots", "umap", "all", "clusters_all_umap.pdf"), useDingbats=F)

#fig 1e - gene heatmap
if (!file.exists(file.path("data", "markers_all.RData"))) {
  markers_all <- FindAllMarkers(all_singl)          
  save(markers_all, file = file.path("data", "markers_all.RData"))
  write_csv(markers_all, file.path("data", "markers_all.csv"))
} else {
  load(file.path("data", "markers_all.RData"))
}

top10 <- markers_all %>% 
  #remove non annotated genes
  filter(!gene %in% grep("(^Gm|^Rp|Rik|mt-|RP)", .$gene, value = T)) %>%
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

heat <- DoHeatmap(all_singl,features = top10$gene, group.colors = colors_many[-14]) 
heat + 
  scale_fill_viridis(option = "B")

ggsave(file.path("plots", "heatmaps", "all", "clusters_all_heatmap.pdf"), useDingbats=F)

DoHeatmap(all_singl,features = top10$gene, group.bar = F)  + 
  theme_void() + 
  scale_fill_viridis(option = "B") +
  NoLegend()
ggsave(file.path("plots", "heatmaps", "all", "clusters_all_heatmap_panel.png"))

#fig 1f - barplot
mosaicGG2(metadata, "seurat_clusters", "Condition", colors = colors_many[3:2]) 

ggsave(file.path("plots", "others", "all", "marimekko_all.pdf"), useDingbats=F)

hyper_test_n(metadata, "seurat_clusters", "Condition") %>%
  write_csv(file.path("data", "marimekko_clusters_condition_all_cells.csv"))

#fig 1f - barplot compartments
mosaicGG2(metadata, "seurat_clusters", "Compartment", colors = colors_pat)  +
  scale_fill_brewer(palette="Dark2", direction = -1)

ggsave(file.path("plots", "others", "all", "marimekko_compartment_all.pdf"), useDingbats=F)

hyper_test_n(metadata, "seurat_clusters", "Compartment") %>%
  write_csv(file.path("data", "marimekko_clusters_compartment_all_cells.csv"))

#fig 1 g - gene expression umaps
genes <- c("P2ry12", "Tmem119", "Hexb", "Apoe", "Gpr34", "Mrc1", "Ms4a7", "Stab1", "Ifitm2", "Ccr2", "Nr4a1", "H2-Aa", "Cd209a", "Sox9", "Mki67", "Clec7a", "Lpl", "Gzmb", "Nkg7", "Klrb1c")

for (i in genes) {
  plt <- plot_expmap_seurat(features=i, object=all_singl, point_size = 3,.retain_cl = levels(all_singl))
  print(plt)
  ggsave(file.path("plots", "umap", "all", paste0(i,"_all.pdf")), useDingbats=F)
}  

#export suppl tables
a <- as.data.frame(table(metadata$seurat_clusters)) 
colnames(a) <- c("Cluster", "Freq")
assertthat::assert_that(sum(a$Freq) == nrow(metadata))
write.csv(a, file.path("data", "Tbl_S1_cluster_cell_count.csv"))
