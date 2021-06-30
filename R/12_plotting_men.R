library(Seurat)
library(tidyverse)
library(readxl)
library(viridis)
library(ggrepel)
library(vegan)
library(ggpubr)

#create plotting subfolders
dir.create(file.path("plots", "umap", "men"))
dir.create(file.path("plots", "heatmaps", "men"))
dir.create(file.path("plots", "others", "men"))

#load functions and colors
source(file.path("R", "functions.R"))

#load data
load(file.path("data", "seurat_men.RData"))

#set cluster order
order_clusters <- data.frame(seurat_clusters= men@meta.data[,"seurat_clusters"], row.names = rownames(men@meta.data)) %>%
  bind_cols(as.data.frame(t(men[["SCT"]]@scale.data))) %>%
  group_by(seurat_clusters) %>%
  summarize_all(.funs=mean) %>%
  as.data.frame()

rownames(order_clusters) <- order_clusters$seurat_clusters
order_clusters <- order_clusters$seurat_clusters[hclust(dist(order_clusters[,-1]))$order]

#reorder the clusters
levels(men) <- order_clusters[c(10:9,2,4:6,3,7,8,1)]

#extract metadata
metadata <- men@meta.data
metadata$seurat_clusters <- factor(metadata$seurat_clusters, levels = order_clusters[c(10:9,2,7:3,8,1)])
metadata$cell_type <- case_when(
  metadata$seurat_clusters %in% c("0", "2", "7")  ~ "CAMs",
  metadata$seurat_clusters %in% c("4") & colSums(as.matrix(men[["SCT"]]@counts[c("Cd209a", "Tnfsf9", "Tnip3", "Kcne3"),])) > 2  ~ "cDC2",
  metadata$seurat_clusters %in% c("1", "3", "5", "6", "8") ~ "Micr",
  metadata$seurat_clusters %in% c("4") ~ "Ly6chi Monocytes",
  metadata$seurat_clusters %in% c("4") ~ "Ly6clo Monocytes",
  TRUE ~ "other"
)
write.csv(metadata, file.path("data", "men_cell_metadata.csv"))
tools::md5sum(file.path("data", "men_cell_metadata.csv"))

#export suppl tables
a <- as.data.frame(table(metadata$seurat_clusters)) 
colnames(a) <- c("Cluster", "Freq")
assertthat::assert_that(sum(a$Freq) == nrow(metadata))
write.csv(a, file.path("data", "Tbl_S17_cluster_cell_count.csv"))

# fig 5a - cell signature umaps
signature_genes <-  read_excel(file.path("data","cell_signatures.xlsx"), "Core signature", skip = 2)

for (i in colnames(signature_genes)) {
  plt <- plot_expmap_seurat(na.omit(signature_genes[[i]]), point_size = 6, object = men, .retain_cl = levels(men)) + labs(subtitle= paste0(i,' Signature'))
  print(plt)
  ggsave(file.path("plots", "umap", "men", paste0(i,"_signature_men.pdf")), useDingbats=F)
}

# fig 5b - conditions umap
men %>%
  DimPlot(group.by = "Condition", pt.size = 6) +
  theme_void() +
  scale_color_manual(values=colors_many[3:2])
ggsave(file.path("plots", "umap", "men", "conditions_men_umap.pdf"), useDingbats=F)

# fig 5c - clusters umap
men %>%
  DimPlot(label=T, pt.size = 6) +
  theme_void() +
  scale_color_manual(values=colors_pat)
ggsave(file.path("plots", "umap", "men", "clusters_men_umap.pdf"), useDingbats=F)

#fig 5d - barplot
mosaicGG2(metadata, "seurat_clusters", "Condition", colors = colors_many[3:2]) 

ggsave(file.path("plots", "others", "men", "marimekko_men.pdf"), useDingbats=F)

hyper_test_n(metadata, "seurat_clusters", "Condition") %>%
  write_csv(file.path("data", "marimekko_clusters_condition_men_cells.csv"))

#fig 5e - gene heatmap
if (!file.exists(file.path("data", "markers_men.RData"))) {
  markers_men <- FindAllMarkers(men)          
  save(markers_men, file = file.path("data", "markers_men.RData"))
  write_csv(markers_men, file.path("data", "markers_men.csv"))
} else {
  load(file.path("data", "markers_men.RData"))
}

top10 <- markers_men %>% 
  #remove non annotated genes
  filter(!gene %in% grep("(^Gm|^Rp|Rik|mt-|RP)", .$gene, value = T)) %>%
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

heat <- DoHeatmap(men,features = top10$gene, group.colors = colors_pat) 
heat + 
  scale_fill_viridis(option = "B")

ggsave(file.path("plots", "heatmaps", "men", "men_heatmap.pdf"), useDingbats=F)


DoHeatmap(men,features = top10$gene, group.bar = F) +
  theme_void() + 
  scale_fill_viridis(option = "B") +
  NoLegend()

ggsave(file.path("plots", "heatmaps", "men", "heatmap_color_panel.png"))

#volcano plot
men2 <- men

Idents(men2) <- paste(metadata$Condition, metadata$cell_type, sep = "_")

#calculate differentially expressed genes
if (!file.exists(file.path("data", "men_diffgenes_spf_vs_gf.RData"))) {
  men2_genes_gf <- FindMarkers(men2, 
                           ident.1 = "GF_CAMs",
                           ident.2 = "SPF_CAMs",
                          logfc.threshold = 0.01,
                          min.pct = 0.01) %>%
    rownames_to_column(var="gene") %>%
    mutate(Condition = "GF")
  
  men2_genes_spf <- FindMarkers(men2, 
                              ident.1 = "SPF_CAMs",
                              ident.2 = "GF_CAMs",
                              logfc.threshold = 0.01,
                              min.pct = 0.01)  %>%
    rownames_to_column(var="gene") %>%
    mutate(Condition = "SPF")

  men2_genes_all <- bind_rows(men2_genes_gf,
                             men2_genes_spf)
  
  save(men2_genes_all, file = file.path("data", "men_diffgenes_spf_vs_gf.RData"))
} else {
  load(file.path("data", "men_diffgenes_spf_vs_gf.RData"))
}
#remove unannotated genes
men2_genes_all <- men2_genes_all %>%
  filter(!grepl("(Gm|Rp|Rik|mt-|RP)", .$gene))

men2_genes_all <- men2_genes_all %>%
  filter(avg_log2FC > 0) %>%
  mutate(
    genes_sig = ifelse(p_val_adj < .05 & avg_log2FC > .2, "sig.", "not sig."),
    show_genes = ifelse(genes_sig == "sig.", gene, NA),
    avg_log2FC = case_when(
      Condition == "SPF"  ~ -1 * avg_log2FC,
      TRUE ~ avg_log2FC
    )
  ) 


men_volcano <- ggplot(men2_genes_all, aes(x=avg_log2FC, y= -log10(p_val_adj), label=show_genes, color=genes_sig)) +
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

men_volcano 

ggsave(file.path("plots", "others", "men", "men_cams_volcano.pdf"), useDingbats=F)

#fig 5 g - gene expression umaps
genes <- c("Mrc1", "Stab1", "P2ry12", "Apoe", "Nr4a1", "Ccr2", "Enpp2", "H2-Aa", "Cd209a", "Hexb", "Fos")

for (i in genes) {
  plt <- plot_expmap_seurat(features=i, object=men,  point_size = 6,.retain_cl = levels(men))
  print(plt)
  ggsave(file.path("plots", "umap", "men", paste0(i,"_men.pdf")), useDingbats=F)
}  


#other analyses

#bray-curtis dissimilarity
set.seed(79106)

#copare based on variable genes - CAMs
a_spf <- men[["SCT"]]@counts[VariableFeatures(men), which(metadata$seurat_clusters %in% c("1", "3", "4"))] %>%
  as.matrix() %>%
  as.data.frame() %>%
  dplyr::select(contains("SPF")) %>%
  t() 

a_spf <- a_spf %>%
  vegdist(method="bray", binary=FALSE, diag=FALSE, upper=FALSE,
          na.rm = FALSE) %>%
  as.numeric()

a_gf <- men[["SCT"]]@counts[VariableFeatures(men), which(metadata$seurat_clusters %in% c("1", "3", "4"))] %>%
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
  labs(title = "men Bray-Curtis Dissimilarity", y="Bray-Curtis Dissimilarity Coefficient")

ggsave(file.path("plots", "others", "men", "men_bray_curtis_violin_plots.pdf"), useDingbats=F)

