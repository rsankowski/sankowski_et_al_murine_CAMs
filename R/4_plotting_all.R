library(Seurat)
library(tidyverse)
library(readxl)
library(viridis)
library(biomaRt)
library(ggpubr)

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

lst <- list()
#plot the correlations of Apoe with cell activation signatures
fls <- list.files(file.path("data", "genesets"))
act_genes <- map(fls, function(x) {
  lst[[x]] <- read_table(file.path("data", "genesets", x), skip = 2, col_names = F)
})

names(act_genes) <- fls
act_genes <- lapply(act_genes, as.data.frame) 

#find murine orthologues
ensembl <- useMart("ensembl")
ensembl.human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#save into new list
act_genes2 <- list()
act_genes2 <- map(2:length(act_genes), function(x) {
  orthologs <- getLDS(attributes = c("mgi_symbol"),
                      filters = "mgi_symbol", values = act_genes[[x]][[1]], mart = ensembl.mouse,
                      attributesL = c("hgnc_symbol") , martL = ensembl.human)
  act_genes2[x] <- data.frame(X1=orthologs$MGI.symbol)
})

#append the dam signature
dam <- data.frame(X1=act_genes[[1]])
act_genes2 <- c(list(dam), act_genes2)

lst_micr <- list()
lst_micr <- map(1:4, function(x) {
  lst_micr[[x]] <- data.frame(cum_expr=colSums(all_singl[["SCT"]]@counts[act_genes2[[x]][[1]][act_genes2[[x]][[1]] != "Apoe" & act_genes2[[x]][[1]] %in% rownames(all_singl[["SCT"]]@counts)], all_singl@active.ident %in% c("0", "5", "10", "2")]),
                              Apoe_expr= c(all_singl[["SCT"]]@counts["Apoe", all_singl@active.ident %in% c("0", "5", "10", "2")])) 
  
})
names(lst_micr) <- c("DAM Signature", "GENESET GO IMMUNE RESPONSE", "GENESET GO IMMUNE SYSTEM PROC.", "GENESET CTRL vs IFNG Micr")


lst_micr <- lst_micr %>% 
  bind_rows(.id="Geneset")

broom::glance(lm(cum_expr ~ Apoe_expr + Geneset, data=lst_micr))

lst_micr %>% 
  ggplot(aes(Apoe_expr, cum_expr)) +
  geom_point(size=2.5) +
  geom_smooth(method="lm", se=F, size=2) +
  theme_pubclean() +
  facet_wrap(~ Geneset) +
  stat_cor(size=5) +
  theme(text = element_text(size=18))

ggsave(file.path("plots", "others", "all", "Apoe_geneset_correlation_micr.pdf"), useDingbats=F)
ggsave(file.path("plots", "others", "all", "Apoe_geneset_correlation_micr.png"))

#same for cams
lst_cams <- list()
lst_cams <- map(1:4, function(x) {
  lst_cams[[x]] <- data.frame(cum_expr=colSums(all_singl[["SCT"]]@counts[act_genes2[[x]][[1]][act_genes2[[x]][[1]] != "Apoe" & act_genes2[[x]][[1]] %in% rownames(all_singl[["SCT"]]@counts)], all_singl@active.ident %in% c("7", "1", "13", "4")]),
                              Apoe_expr=c(all_singl[["SCT"]]@counts["Apoe", all_singl@active.ident %in% c("7", "1", "13", "4")])) 
  
})

names(lst_cams) <- c("DAM Signature", "GENESET GO IMMUNE RESPONSE", "GENESET GO IMMUNE SYSTEM PROC.", "GENESET CTRL vs IFNG micr")


lst_cams <- lst_cams %>% 
  bind_rows(.id="Geneset")

broom::glance(lm(cum_expr ~ Apoe_expr + Geneset, data=lst_cams))

lst_cams %>% 
  ggplot(aes(Apoe_expr, cum_expr)) +
  geom_point(size=2.5) +
  geom_smooth(method="lm", se=F, size=2) +
  theme_pubclean() +
  facet_wrap(~ Geneset) +
  stat_cor(size=5) +
  theme(text = element_text(size=18))

ggsave(file.path("plots", "others", "all", "Apoe_geneset_correlation_cams.pdf"), useDingbats=F)
ggsave(file.path("plots", "others", "all", "Apoe_geneset_correlation_cams.png"))

