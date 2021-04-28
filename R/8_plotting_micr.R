library(Seurat)
library(tidyverse)
library(readxl)
library(viridis)
library(ggrepel)
library(vegan)
library(ggpubr)

#create plotting subfolders
dir.create(file.path("plots", "umap", "micr"))
dir.create(file.path("plots", "heatmaps", "micr"))
dir.create(file.path("plots", "others", "micr"))

#load functions and colors
source(file.path("R", "functions.R"))

#load data
load(file.path("data", "seurat_micr.RData"))

#set cluster order
order_clusters <- data.frame(seurat_clusters= micr@meta.data[,"seurat_clusters"], row.names = rownames(micr@meta.data)) %>%
  bind_cols(as.data.frame(t(micr[["SCT"]]@scale.data))) %>%
  group_by(seurat_clusters) %>%
  summarize_all(.funs=mean) %>%
  as.data.frame()

rownames(order_clusters) <- order_clusters$seurat_clusters
order_clusters <- order_clusters$seurat_clusters[hclust(dist(order_clusters[,-1]))$order]

#reorder the clusters
levels(micr) <- order_clusters

#extract metadata
metadata <- micr@meta.data
metadata$seurat_clusters <- factor(metadata$seurat_clusters, levels = levels(micr))

#export suppl tables
a <- as.data.frame(table(metadata$seurat_clusters)) 
colnames(a) <- c("Cluster", "Freq")
assertthat::assert_that(sum(a$Freq) == nrow(metadata))
write.csv(a, file.path("data", "Tbl_S5_cluster_cell_count.csv"))

# fig 2a - cell signature umaps
signature_genes <-  read_excel(file.path("data","cell_signatures.xlsx"), "Core signature", skip = 2)

plot_expmap_seurat(na.omit(signature_genes[["Microglia"]]), point_size = 6, object = micr, .retain_cl = levels(micr)) + labs(subtitle= 'Microglia Signature')

ggsave(file.path("plots", "umap", "micr", "microglia_signature_micr.pdf"), useDingbats=F)

# fig 2b - conditions umap
micr %>%
  DimPlot(group.by = "Condition", pt.size = 6) +
  theme_void() +
  scale_color_manual(values=colors_many[3:2])
ggsave(file.path("plots", "umap", "micr", "conditions_micr_umap.pdf"), useDingbats=F)

# fig 2c - clusters umap
micr %>%
  DimPlot(label=T, pt.size = 6) +
  theme_void() +
  scale_color_manual(values=colors_pat)
ggsave(file.path("plots", "umap", "micr", "clusters_micr_umap.pdf"), useDingbats=F)

#fig 2d - barplot
mosaicGG2(metadata, "seurat_clusters", "Condition", colors = colors_many[3:2]) 

ggsave(file.path("plots", "others", "micr", "marimekko_micr.pdf"), useDingbats=F)

hyper_test_n(metadata, "seurat_clusters", "Condition") %>%
  write_csv(file.path("data", "marimekko_clusters_condition_micr_cells.csv"))

#fig 2e - gene heatmap
if (!file.exists(file.path("data", "markers_micr.RData"))) {
  markers_micr <- FindAllMarkers(micr)          
  save(markers_micr, file = file.path("data", "markers_micr.RData"))
  write_csv(markers_micr, file.path("data", "markers_micr.csv"))
} else {
  load(file.path("data", "markers_micr.RData"))
}

top10 <- markers_micr %>% 
  #remove non annotated genes
  filter(!gene %in% grep("(^Gm|^Rp|Rik|mt-|RP)", .$gene, value = T)) %>%
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

heat <- DoHeatmap(micr,features = top10$gene, group.colors = colors_pat) 
heat + 
  scale_fill_viridis(option = "B")

ggsave(file.path("plots", "heatmaps", "micr", "clusters_micr_heatmap.pdf"), useDingbats=F)

DoHeatmap(micr,features = top10$gene, group.bar = F) +
  theme_void() + 
  scale_fill_viridis(option = "B") +
  NoLegend()

ggsave(file.path("plots", "heatmaps", "micr", "clusters_micr_heatmap_color_panel.png"))

#fig 2 - bray-curtis dissimilarity
set.seed(79106)

#copare based on variable genes
a_spf <- micr[["SCT"]]@counts[VariableFeatures(micr),] %>%
  as.matrix() %>%
  as.data.frame() %>%
  dplyr::select(contains("SPF")) %>%
  t() 

a_spf <- a_spf %>%
  vegdist(method="bray", binary=FALSE, diag=FALSE, upper=FALSE,
          na.rm = FALSE) %>%
  as.numeric()

a_gf <- micr[["SCT"]]@counts[VariableFeatures(micr),] %>%
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
  labs(title = "MG Bray-Curtis Dissimilarity", y="Bray-Curtis Dissimilarity Coefficient")

ggsave(file.path("plots", "others", "micr", "micr_bray_curtis_violin_plots.pdf"), useDingbats=F)

#fig 2 - volcano plot

micr2 <- micr
Idents(micr2) <- micr2@meta.data$Condition

if (!file.exists(file.path("data", "micr_diffgenes_spf_vs_gf.RData"))) {
  micr2_genes <- FindAllMarkers(micr2, 
                       logfc.threshold = 0.01,
                      min.pct = 0.01) 

  save(micr2_genes, file = file.path("data", "micr_diffgenes_spf_vs_gf.RData"))
} else {
  load(file.path("data", "micr_diffgenes_spf_vs_gf.RData"))
}

#export suppl table 8
write.csv(micr2_genes, file.path("data", "Table_S8_micr.csv"))
#remove unannotated genes
micr2_genes <- micr2_genes %>%
  filter(!grepl("(Gm|Rp|Rik|mt-|RP)", .$gene))

micr2_genes <- micr2_genes %>%
  filter(avg_log2FC > 0) %>%
  mutate(
    genes_sig = ifelse(p_val_adj < .05 & avg_log2FC > .2, "sig.", "not sig."),
    show_genes = ifelse(genes_sig == "sig.", gene, NA),
    avg_log2FC = case_when(
      cluster == "SPF"  ~ -1 * avg_log2FC,
      TRUE ~ avg_log2FC
    )
  ) 


micr_volcano <- ggplot(micr2_genes, aes(x=avg_log2FC, y= -log10(p_val_adj), label=show_genes, color=genes_sig)) +
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

micr_volcano 

ggsave(file.path("plots", "others", "micr", "micr_volcano.pdf"), useDingbats=F)

#fig 2 g - gene expression umaps
genes <- c("Hexb", "Cst3", "P2ry12", "Apoe")

for (i in genes) {
  plt <- plot_expmap_seurat(features=i, object=micr,  point_size = 6,.retain_cl = levels(micr))
  print(plt)
  ggsave(file.path("plots", "umap", "micr", paste0(i,"_micr.pdf")), useDingbats=F)
}  

# plot violin plots
genes <- c("Ddit4", 
           "C1qb",
           "Rnase4",
           "Apoe",
           "Cd180",
           "Ly86",
           "P2ry12",
           "Fcgr2b")

map(genes, function(x) {
a <- data.frame(gene = micr[["SCT"]]@counts[x,], "ID" = colnames(micr[["SCT"]]@counts))
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
ggsave(file.path('plots','others', 'micr', paste0('Cluster-violin-plot-micr-', x,'.pdf')), useDingbats=F)
})

