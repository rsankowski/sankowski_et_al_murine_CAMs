library(Seurat)
library(tidyverse)
library(readxl)
library(viridis)
library(ggrepel)
library(vegan)
library(ggpubr)

#create plotting subfolders
dir.create(file.path("plots", "umap", "pvm"))
dir.create(file.path("plots", "heatmaps", "pvm"))
dir.create(file.path("plots", "others", "pvm"))

#load functions and colors
source(file.path("R", "functions.R"))

#load data
load(file.path("data", "seurat_pvm.RData"))

#set cluster order
order_clusters <- data.frame(seurat_clusters= pvm@meta.data[,"seurat_clusters"], row.names = rownames(pvm@meta.data)) %>%
  bind_cols(as.data.frame(t(pvm[["SCT"]]@scale.data))) %>%
  group_by(seurat_clusters) %>%
  summarize_all(.funs=mean) %>%
  as.data.frame()

rownames(order_clusters) <- order_clusters$seurat_clusters
order_clusters <- order_clusters$seurat_clusters[hclust(dist(order_clusters[,-1]))$order]

#reorder the clusters
levels(pvm) <- order_clusters[c(8:7,3,6:4,10:9,2,1)]

#extract metadata
metadata <- pvm@meta.data
metadata$seurat_clusters <- factor(metadata$seurat_clusters, levels = order_clusters[c(8:7,3,6:4,10:9,2,1)])
metadata$cell_type <- case_when(
  metadata$seurat_clusters %in% c("0", "2", "6", "9")   ~ "CAMs",
  metadata$seurat_clusters %in% c("5")  ~ "cDC2",
  metadata$seurat_clusters %in% c("7") ~ "Micr",
  metadata$seurat_clusters %in% c("7", "3") ~ "Monocytes",
  metadata$seurat_clusters %in% c("8") ~ "Lymphocytes",
  TRUE ~ "other"
)

#export metadata
write.csv(metadata, file.path("data", "pvm_cell_metadata.csv"))
tools::md5sum(file.path("data", "pvm_cell_metadata.csv"))

#export suppl tables
a <- as.data.frame(table(metadata$seurat_clusters)) 
colnames(a) <- c("Cluster", "Freq")
assertthat::assert_that(sum(a$Freq) == nrow(metadata))
write.csv(a, file.path("data", "Tbl_S13_cluster_cell_count_pvm.csv"))
  
# fig 4a - cell signature umaps
signature_genes <-  read_excel(file.path("data","cell_signatures.xlsx"), "Core signature", skip = 2)

for (i in colnames(signature_genes)) {
  plt <- plot_expmap_seurat(na.omit(signature_genes[[i]]), point_size = 6, object = pvm, .retain_cl = levels(pvm)) + labs(subtitle= paste0(i,' Signature'))
  print(plt)
  ggsave(file.path("plots", "umap", "pvm", paste0(i,"_signature_pvm.pdf")), useDingbats=F)
}

# fig 4b - conditions umap
pvm %>%
  DimPlot(group.by = "Condition", pt.size = 6) +
  theme_void() +
  scale_color_manual(values=colors_many[3:2])
ggsave(file.path("plots", "umap", "pvm", "conditions_pvm_umap.pdf"), useDingbats=F)

# fig 4c - clusters umap
pvm %>%
  DimPlot(label=T, pt.size = 6) +
  theme_void() +
  scale_color_manual(values=colors_pat)
ggsave(file.path("plots", "umap", "pvm", "clusters_pvm_umap.pdf"), useDingbats=F)

#fig 4d - barplot
mosaicGG2(metadata, "seurat_clusters", "Condition", colors = colors_many[3:2]) 

ggsave(file.path("plots", "others", "pvm", "marimekko_pvm.pdf"), useDingbats=F)

hyper_test_n(metadata, "seurat_clusters", "Condition") %>%
  write_csv(file.path("data", "marimekko_clusters_condition_pvm_cells.csv"))

#fig 4e - gene heatmap
if (!file.exists(file.path("data", "markers_pvm.RData"))) {
  markers_pvm <- FindAllMarkers(pvm)          
  save(markers_pvm, file = file.path("data", "markers_pvm.RData"))
  write_csv(markers_pvm, file.path("data", "markers_pvm.csv"))
} else {
  load(file.path("data", "markers_pvm.RData"))
}

top10 <- markers_pvm %>% 
  #remove non annotated genes
  filter(!gene %in% grep("(^Gm|^Rp|Rik|mt-|RP)", .$gene, value = T)) %>%
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

heat <- DoHeatmap(pvm,features = top10$gene, group.colors = colors_pat) 
heat + 
  scale_fill_viridis(option = "B")

ggsave(file.path("plots", "heatmaps", "pvm", "pvm_heatmap.pdf"), useDingbats=F)

DoHeatmap(pvm,features = top10$gene, group.bar = F) +
  theme_void() + 
  scale_fill_viridis(option = "B") +
  NoLegend()

ggsave(file.path("plots", "heatmaps", "pvm", "heatmap_color_panel.png"))


#volcano plot
pvm2 <- pvm

Idents(pvm2) <- paste(metadata$Condition, metadata$cell_type, sep = "_")

#epiplexus
if (!file.exists(file.path("data", "pvm_diffgenes_spf_vs_gf.RData"))) {
  pvm2_genes_gf <- FindMarkers(pvm2, 
                           ident.1 = "GF_CAMs",
                           ident.2 = "SPF_CAMs",
                          logfc.threshold = 0.01,
                          min.pct = 0.01) %>%
    rownames_to_column(var="gene") %>%
    mutate(Condition = "GF")
  
  pvm2_genes_spf <- FindMarkers(pvm2, 
                              ident.1 = "SPF_CAMs",
                              ident.2 = "GF_CAMs",
                              logfc.threshold = 0.01,
                              min.pct = 0.01)  %>%
    rownames_to_column(var="gene") %>%
    mutate(Condition = "SPF")

  pvm2_genes_all <- bind_rows(pvm2_genes_gf,
                             pvm2_genes_spf)
  
  save(pvm2_genes_all, file = file.path("data", "pvm_diffgenes_spf_vs_gf.RData"))
} else {
  load(file.path("data", "pvm_diffgenes_spf_vs_gf.RData"))
}

write.csv(pvm2_genes_all, file.path("data", "Table_S12_pvm_DEGs.csv"))

#remove unannotated genes
pvm2_genes_all <- pvm2_genes_all %>%
  filter(!grepl("(Gm|Rp|Rik|mt-|RP)", .$gene))

pvm2_genes_all <- pvm2_genes_all %>%
  filter(avg_log2FC > 0) %>%
  mutate(
    genes_sig = ifelse(p_val_adj < .05 & avg_log2FC > .2, "sig.", "not sig."),
    show_genes = ifelse(genes_sig == "sig.", gene, NA),
    avg_log2FC = case_when(
      Condition == "SPF"  ~ -1 * avg_log2FC,
      TRUE ~ avg_log2FC
    )
  ) 


pvm_volcano <- ggplot(pvm2_genes_all, aes(x=avg_log2FC, y= -log10(p_val_adj), label=show_genes, color=genes_sig)) +
  annotate("rect", xmin=-Inf,xmax=0,ymin=-Inf,ymax=Inf,
            fill=colors_many[2], alpha =.4)+
  annotate("rect" ,xmin=0,xmax=Inf,ymin=-Inf,ymax=Inf,
            fill=colors_many[3], alpha =.4) +
  geom_point(size=5) + 
  geom_text_repel(size=7, box.padding=1.15, max.overlaps = 20) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size=25)) +
  scale_color_manual(values = c("light grey", "black")) +
  NoLegend() +
  labs(x="avg. log2FC", y="-log10 transf. adj. p-value") 

pvm_volcano 

ggsave(file.path("plots", "others", "pvm", "pvm_cams_volcano.pdf"), useDingbats=F)


#fig 4 g - gene expression umaps
genes <- c("Mrc1", "Stab1", "P2ry12", "Apoe", "Nr4a1", "Ccr2", "Enpp2", "H2-Aa", "Cd209a", "Hexb", "Fos", "Nkg7", "Tmem119", "Olfml3","Sall1")

for (i in genes) {
  plt <- plot_expmap_seurat(features=i, object=pvm,  point_size = 6,.retain_cl = levels(pvm))
  print(plt)
  ggsave(file.path("plots", "umap", "pvm", paste0(i,"_pvm.pdf")), useDingbats=F)
}  

# violin plots
genes <- c("Ly86",
           "Fcrls",
           "Cd209b")

map(genes, function(x) {
  a <- data.frame(gene = pvm[["SCT"]]@counts[x,], "ID" = colnames(pvm[["SCT"]]@counts))
  a <- metadata[,c('Condition', 'seurat_clusters')] %>%
    rownames_to_column(var="ID") %>%
    right_join(a) %>%
    na.omit %>%
    mutate(Condition = factor(.$Condition, levels = c("SPF", "GF")))
  
  plt <- ggplot(a, aes(seurat_clusters, gene, fill = Condition)) +
    geom_violin(scale="width") +
    #geom_jitter(pch=21, size = 3, width = 0.1) +
    labs(title = x) +
    scale_fill_manual("Status", values = colors_many[2:4]) +
    theme_minimal() 
  
  print(plt)
  ggsave(file.path('plots','others', 'pvm', paste0('Cluster-violin-plot-pvm-', x,'.pdf')), useDingbats=F)
})

#other analyses
#bray-curtis dissimilarity
set.seed(79106)

#compare based on variable genes - CAMs
a_spf <- pvm[["SCT"]]@counts[VariableFeatures(pvm), which(metadata$seurat_clusters %in% c("3", "4"))] %>%
  as.matrix() %>%
  as.data.frame() %>%
  dplyr::select(contains("SPF")) %>%
  t() 

a_spf <- a_spf %>%
  vegdist(method="bray", binary=FALSE, diag=FALSE, upper=FALSE,
          na.rm = FALSE) %>%
  as.numeric()

a_gf <- pvm[["SCT"]]@counts[VariableFeatures(pvm), which(metadata$seurat_clusters %in% c("3", "4"))] %>%
  as.matrix() %>%
  as.data.frame() %>%
  dplyr::select(contains("GF")) %>%
  t() 

a_gf <- a_gf %>%
  vegdist(method="bray", binary=FALSE, diag=FALSE, upper=FALSE,
          na.rm = FALSE) %>%
  as.numeric()

all <- bind_rows(data.frame(condition="GF", dist=a_gf),
                 data.frame(condition="SPF", dist=a_spf))
all$condition <- factor(all$condition, levels = c("SPF", "GF"))

all %>%
  ggplot(aes(condition, dist, fill=condition)) +
  geom_violin(draw_quantiles = 0.5) +
  scale_fill_manual(values = colors_many[2:4]) +
  stat_compare_means() +
  theme_linedraw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=20)) +
  labs(title = "pvm Bray-Curtis Dissimilarity", y="Bray-Curtis Dissimilarity Coefficient")

ggsave(file.path("plots", "others", "pvm", "pvm_epiplexus_bray_curtis_violin_plots.pdf"), useDingbats=F)

 