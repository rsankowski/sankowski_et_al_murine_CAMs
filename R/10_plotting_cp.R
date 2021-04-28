library(Seurat)
library(tidyverse)
library(readxl)
library(viridis)
library(ggrepel)
library(vegan)
library(ggpubr)

#create plotting subfolders
dir.create(file.path("plots", "umap", "cp"))
dir.create(file.path("plots", "heatmaps", "cp"))
dir.create(file.path("plots", "others", "cp"))

#load functions and colors
source(file.path("R", "functions.R"))

#load data
load(file.path("data", "seurat_cp.RData"))

#set cluster order
order_clusters <- data.frame(seurat_clusters= cp@meta.data[,"seurat_clusters"], row.names = rownames(cp@meta.data)) %>%
  bind_cols(as.data.frame(t(cp[["SCT"]]@scale.data))) %>%
  group_by(seurat_clusters) %>%
  summarize_all(.funs=mean) %>%
  as.data.frame()

rownames(order_clusters) <- order_clusters$seurat_clusters
order_clusters <- order_clusters$seurat_clusters[hclust(dist(order_clusters[,-1]))$order]

#reorder the clusters
levels(cp) <- order_clusters[c(6,3:2,4,7:8,5,1)]

#extract metadata
metadata <- cp@meta.data
metadata$seurat_clusters <- factor(metadata$seurat_clusters, levels = levels(cp))
metadata$cell_type <- case_when(
  metadata$seurat_clusters %in% c("2", "3", "4") ~ "CP_epi",
  metadata$seurat_clusters %in% c("0", "1") ~ "CP_MF",
  TRUE ~ "other"
)

#export suppl tables
a <- as.data.frame(table(metadata$seurat_clusters)) 
colnames(a) <- c("Cluster", "Freq")
assertthat::assert_that(sum(a$Freq) == nrow(metadata))
write.csv(a, file.path("data", "Tbl_S9_cluster_cell_count.csv"))

# fig 3a - cell signature umaps
signature_genes <-  read_excel(file.path("data","cell_signatures.xlsx"), "Core signature", skip = 2)

for (i in colnames(signature_genes)) {
  plt <- plot_expmap_seurat(na.omit(signature_genes[[i]]), point_size = 6, object = cp, .retain_cl = levels(cp)) + labs(subtitle= paste0(i,' Signature'))
  print(plt)
  ggsave(file.path("plots", "umap", "cp", paste0(i,"_signature_cp.pdf")), useDingbats=F)
}

# fig 3b - conditions umap
cp %>%
  DimPlot(group.by = "Condition", pt.size = 6) +
  theme_void() +
  scale_color_manual(values=colors_many[3:2])
ggsave(file.path("plots", "umap", "cp", "conditions_cp_umap.pdf"), useDingbats=F)

# fig 3c - clusters umap
cp %>%
  DimPlot(label=T, pt.size = 6) +
  theme_void() +
  scale_color_manual(values=colors_pat)
ggsave(file.path("plots", "umap", "cp", "clusters_cp_umap.pdf"), useDingbats=F)

#fig 3d - barplot
mosaicGG2(metadata, "seurat_clusters", "Condition", colors = colors_many[3:2]) 

ggsave(file.path("plots", "others", "cp", "marimekko_cp.pdf"), useDingbats=F)

hyper_test_n(metadata, "seurat_clusters", "Condition") %>%
  write_csv(file.path("data", "marimekko_clusters_condition_cp_cells.csv"))

#fig 3e - gene heatmap
if (!file.exists(file.path("data", "markers_cp.RData"))) {
  markers_cp <- FindAllMarkers(cp)          
  save(markers_cp, file = file.path("data", "markers_cp.RData"))
  write_csv(markers_cp, file.path("data", "markers_cp.csv"))
} else {
  load(file.path("data", "markers_cp.RData"))
}

top10 <- markers_cp %>% 
  #remove non annotated genes
  filter(!gene %in% grep("(^Gm|^Rp|Rik|mt-|RP)", .$gene, value = T)) %>%
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

heat <- DoHeatmap(cp,features = top10$gene, group.colors = colors_pat) 
heat + 
  scale_fill_viridis(option = "B")

ggsave(file.path("plots", "heatmaps", "cp", "cp_heatmap.pdf"), useDingbats=F)

DoHeatmap(cp,features = top10$gene, group.bar = F) +
  theme_void() + 
  scale_fill_viridis(option = "B") +
  NoLegend()

ggsave(file.path("plots", "heatmaps", "cp", "clusters_micr_heatmap_color_panel.png"))


#fig 3 - bray-curtis dissimilarity
set.seed(79106)

#copare based on variable genes - epiplexus
a_spf <- cp[["SCT"]]@counts[VariableFeatures(cp), which(metadata$seurat_clusters %in% c("2", "3", "4"))] %>%
  as.matrix() %>%
  as.data.frame() %>%
  dplyr::select(contains("SPF")) %>%
  t() 

a_spf <- a_spf %>%
  vegdist(method="bray", binary=FALSE, diag=FALSE, upper=FALSE,
          na.rm = FALSE) %>%
  as.numeric()

a_gf <- cp[["SCT"]]@counts[VariableFeatures(cp), which(metadata$seurat_clusters %in% c("2", "3", "4"))] %>%
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
  labs(title = "CP epiplexus Bray-Curtis Dissimilarity", y="Bray-Curtis Dissimilarity Coefficient")

ggsave(file.path("plots", "others", "cp", "cp_epiplexus_bray_curtis_violin_plots.pdf"), useDingbats=F)

#bray curtis for CP MFs
#copare based on variable genes - MF
a_spf <- cp[["SCT"]]@counts[VariableFeatures(cp), which(metadata$seurat_clusters %in% c("0", "1"))] %>%
  as.matrix() %>%
  as.data.frame() %>%
  dplyr::select(contains("SPF")) %>%
  t() 

a_spf <- a_spf %>%
  vegdist(method="bray", binary=FALSE, diag=FALSE, upper=FALSE,
          na.rm = FALSE) %>%
  as.numeric()

a_gf <- cp[["SCT"]]@counts[VariableFeatures(cp), which(metadata$seurat_clusters %in% c("0", "1"))] %>%
  as.matrix() %>%
  as.data.frame() %>%
  dplyr::select(contains("GF")) %>%
  t() 
#a_gf <- a_gf[sample(1:n_cells),]

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
  labs(title = "CP MF Bray-Curtis Dissimilarity", y="Bray-Curtis Dissimilarity Coefficient")

ggsave(file.path("plots", "others", "cp","cp_MF_bray_curtis_violin_plots.pdf"), useDingbats=F)


#fig 3 - volcano plot
cp2 <- cp

Idents(cp2) <- paste(metadata$Condition, metadata$cell_type, sep = "_")

#epiplexus
if (!file.exists(file.path("data", "cp_epi_diffgenes_spf_vs_gf.RData"))) {
  cp2_genes_gf <- FindMarkers(cp2, 
                           ident.1 = "GF_CP_epi",
                           ident.2 = "SPF_CP_epi",
                          logfc.threshold = 0.01,
                          min.pct = 0.01) %>%
    rownames_to_column(var="gene") %>%
    mutate(Condition = "GF")
  
  cp2_genes_spf <- FindMarkers(cp2, 
                              ident.1 = "SPF_CP_epi",
                              ident.2 = "GF_CP_epi",
                              logfc.threshold = 0.01,
                              min.pct = 0.01)  %>%
    rownames_to_column(var="gene") %>%
    mutate(Condition = "SPF")

  cp2_genes_all <- bind_rows(cp2_genes_gf,
                             cp2_genes_spf)
  
  save(cp2_genes_all, file = file.path("data", "cp_epi_diffgenes_spf_vs_gf.RData"))
} else {
  load(file.path("data", "cp_epi_diffgenes_spf_vs_gf.RData"))
}

write.csv(cp2_genes_all, file.path("data", "Table_S12a_cp_epi_DEGs.csv"))

#remove unannotated genes
cp2_genes_all <- cp2_genes_all %>%
  filter(!grepl("(Gm|Rp|Rik|mt-|RP)", .$gene))

cp2_genes_all <- cp2_genes_all %>%
  filter(avg_log2FC > 0) %>%
  mutate(
    genes_sig = ifelse(p_val_adj < .05 & avg_log2FC > .2, "sig.", "not sig."),
    show_genes = ifelse(genes_sig == "sig.", gene, NA),
    avg_log2FC = case_when(
      Condition == "SPF"  ~ -1 * avg_log2FC,
      TRUE ~ avg_log2FC
    )
  ) 


cp_volcano <- ggplot(cp2_genes_all, aes(x=avg_log2FC, y= -log10(p_val_adj), label=show_genes, color=genes_sig)) +
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

cp_volcano 

ggsave(file.path("plots", "others", "cp", "cp_epi_volcano.pdf"), useDingbats=F, height = 9.57, width = 6.25)


#cp MF
if (!file.exists(file.path("data", "cp_MF_diffgenes_spf_vs_gf.RData"))) {
  cp2_genes_gf <- FindMarkers(cp2, 
                              ident.1 = "GF_CP_MF",
                              ident.2 = "SPF_CP_MF",
                              logfc.threshold = 0.01,
                              min.pct = 0.01) %>%
    rownames_to_column(var="gene") %>%
    mutate(Condition = "GF")
  
  cp2_genes_spf <- FindMarkers(cp2, 
                               ident.1 = "SPF_CP_MF",
                               ident.2 = "GF_CP_MF",
                               logfc.threshold = 0.01,
                               min.pct = 0.01)  %>%
    rownames_to_column(var="gene") %>%
    mutate(Condition = "SPF")
  
  cp2_genes_all <- bind_rows(cp2_genes_gf,
                             cp2_genes_spf)
  
  save(cp2_genes_all, file = file.path("data", "cp_MF_diffgenes_spf_vs_gf.RData"))
} else {
  load(file.path("data", "cp_MF_diffgenes_spf_vs_gf.RData"))
}

write.csv(cp2_genes_all, file.path("data", "Table_S12b_cp_MF_DEGs.csv"))

#remove unannotated genes
cp2_genes_all <- cp2_genes_all %>%
  filter(!grepl("(Gm|Rp|Rik|mt-|RP)", .$gene))

cp2_genes_all <- cp2_genes_all %>%
  filter(avg_log2FC > 0) %>%
  mutate(
    genes_sig = ifelse(p_val_adj < .05 & avg_log2FC > .2, "sig.", "not sig."),
    show_genes = ifelse(genes_sig == "sig.", gene, NA),
    avg_log2FC = case_when(
      Condition == "SPF"  ~ -1 * avg_log2FC,
      TRUE ~ avg_log2FC
    )
  ) 


cp_volcano <- ggplot(cp2_genes_all, aes(x=avg_log2FC, y= -log10(p_val_adj), label=show_genes, color=genes_sig)) +
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

cp_volcano 

ggsave(file.path("plots", "others", "cp", "cp_MF_volcano.pdf"), useDingbats=F, height = 9.57, width = 6.25)

#fig 3 g - gene expression umaps
genes <- c("Hexb", "Cst3", "P2ry12", "Apoe", "Mki67", "Pcp4l1", "Stab1", "Ttr", "Mrc1", "Tmem119")

for (i in genes) {
  plt <- plot_expmap_seurat(features=i, object=cp,  point_size = 6,.retain_cl = levels(cp))
  print(plt)
  ggsave(file.path("plots", "umap", "cp", paste0(i,"_cp.pdf")), useDingbats=F)
}  

# violin plots
genes <- c("Ly86",
           "Apoe",
           "Cd52")

map(genes, function(x) {
  a <- data.frame(gene = cp[["SCT"]]@counts[x,], "ID" = colnames(cp[["SCT"]]@counts))
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
  ggsave(file.path('plots','others', 'cp', paste0('Cluster-violin-plot-cp-', x,'.pdf')), useDingbats=F)
})

