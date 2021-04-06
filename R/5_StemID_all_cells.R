library(RaceID)
library(Seurat)
library(tidyverse)
library(assertthat)

date <- Sys.Date()

#load functions and colors
source(file.path("R", "functions.R"))

#load data
load(file.path("data", "seurat_all_without_doublets.RData"))

# 1.Run RaceID (control cell filtering with mintotal, i.e. total UMI count per cell) -> sc object
if (!file.exists(file.path("data", "sc.RData"))) {
  
  sc <- SCseq(all_singl[["SCT"]]@counts)
  
  #filter data
  sc <- filterdata(sc, 
                   mintotal=100,
                   minnumber = 1,
                   knn=10,
                   minexpr = 1)
  
  # 2.Run Seurat with filtering such that the same cells are retained
  assert_that(length(colnames(all_singl)) == length(colnames(sc@ndata)))
  
  # 3.Re-initialize RaceID output with Seurat data:
  
  part <- as.numeric(as.character(all_singl@meta.data$seurat_clusters))
  d <- as.matrix( dist(all_singl@reductions$pca@cell.embeddings) )
  tsne <- as.data.frame(all_singl@reductions$umap@cell.embeddings)
  names(part) <- colnames(sc@ndata)
  
  n <- colnames(sc@ndata)
  part <- part[n]
  
  # partition
  sc@cpart <- sc@cluster$kpart <- part
  # distances
  sc@distances <- d[n,n]
  # tsne
  sc@tsne <- tsne[n,]
  rm(d, tsne)
  
  sc@medoids <- compmedoids(sc, sc@cpart)
  
  # colors for clusters
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
  idx <- which(order_clusters %in% unique(as.character(sc@cpart)))
  sc@cpart <- factor(sc@cpart, levels = order_clusters)
  
  sc@fcol <- c(colors_pat, colors_many, 'steelblue', colors_many[-14])[idx] #sample(rainbow(max(sc@cpart)))
  
  save(sc, file="data/sc.RData")
  
} else {
  load("data/sc.RData")
}


#StemID
library(FateID)
library(readxl)
library(limma)
library(MASS)

begin <- Sys.time()

if (!file.exists("data/ltr_all.RData")){ #data/ltr-larger-clusters.RData
  ltr <- Ltree(sc)
  
  #convert clusters in integers
  ltr@sc@cpart <- as.numeric(ltr@sc@cpart) +1
  names(ltr@sc@cpart) <- colnames(ltr@sc@ndata)
  
  ltr <- compentropy(ltr)
  ltr <- projcells(ltr,nmode=TRUE,fr=FALSE) #400
  ltr <- projback(ltr,pdishuf=100)
  ltr <- lineagegraph(ltr)
  ltr <- comppvalue(ltr,pthr=0.01)
  
  save(ltr, file = 'data/ltr_all.RData')
} else {
  load('data/ltr_all.RData')
}

plotgraph(ltr,scthr=0.9,showCells=FALSE)

svg(paste0('plots/umap/lineage-graph-tsne-plot.svg'), width = 8.57, height = 5.79)
plotgraph(ltr,scthr=0.9,showCells=FALSE)
dev.off()

pdf(paste0('plots/umap/lineage-graph-tsne-plot.pdf'), width = 8.57, height = 5.79)
plotgraph(ltr,scthr=0.9,showCells=FALSE,)
dev.off()

x <- compscore(ltr,scthr=0.9)
plotdistanceratio(ltr)
plotspantree(ltr)
plotprojections(ltr)

#lineage tree for moDCs
ccr7 <- c(1,4,2)
ifn_cl <- c(1,7,3)

#pseudotemporal ccr7
n <- cellsfromtree(ltr,ccr7)
x <- getfdata(ltr@sc)

fs  <- filterset(x,n=n$f, minexpr = 1, minnumber = 1)

if (!file.exists("data/s1d-ccr7-new.Robj")) {
  s1d <- getsom(fs,nb=1000,alpha=.5)
  save(s1d, file = "data/s1d-ccr7-new.Robj")
} else {
  load("data/s1d-ccr7-new.Robj")
}
ps  <- procsom(s1d,corthr=.85,minsom=3)
y    <- ltr@sc@cpart[n$f]
fcol <- sc@fcol
plotheatmap(ps$all_singl_ea.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)

pdf(paste0('plots/heatmaps/', date, '-ccr7-trajectory-heatmap.pdf'), width = 8.57, height = 5.79)
plotheatmap(ps$all_singl_ea.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
dev.off()

png(paste0('plots/heatmaps/', date, '-ccr7-trajectory-heatmap.png'), width = 8.57, height = 5.79, res = 600, units = "in")
plotheatmap(ps$all_singl_ea.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
dev.off()

#export node genes
modules <- data.frame('Node' = NA, 'Genes' = NA)
for (i in 1: max(ps$nodes)) {
  gene_names <- names(ps$nodes)[ps$nodes == i]
  gene_names <- gsub('_.*', '', gene_names)
  modules2 <- data.frame('Node' = NA, 'Genes' = gene_names)
  modules2$Node <- rep(as.character(i), nrow(modules2))
  modules <- rbind(na.omit(modules), modules2)
}

write_csv(modules, paste0('data/',date,'-nodes-stemid-vector-ccr7-new.csv'))

#find TFs in node genes
tfs <- read_csv("http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv")

tfs_found <- modules[modules$Genes %in% tfs[[3]][tfs$`Binding mode` != "Not a DNA binding protein"],]

#export TFs
write_csv(tfs_found, paste0('data/',date,'-nodes-stemid-TFs-ccr7-new.csv'))

#plot expression of TFs
#adjust gene names
lst = list()
apply(data.frame(genes=as.character(tfs_found[,2])), 1, function(i) 
{
  .g=name2id(paste0("^(", i , ")"), rownames(sc@ndata))
  .n=n$f
  pdf(paste0("plots/others/", i, ".pdf"))
  plotexpression(fs,y,.g,.n,col=fcol,name=i,cluster=FALSE,alpha=.5,types=NULL)
  dev.off()
  return(lst[[i]] <- unlist(x[.g,.n]))
})

lst <- list()
lst <- apply(data.frame(genes=as.character(tfs_found[,2])), 1, function(i) {
  .g=name2id(paste0("^(", i , ")"), rownames(sc@ndata))
  .n=n$f
  return(x[.g,.n])
})
colnames(lst) <- tfs_found[,2]
library(broom)
library(MASS)

trajectory <- 1:nrow(lst)
lst <- as.data.frame(lst)
a <- apply(lst, 2, function(gene) tidy(glm.nb(gene ~ trajectory))) %>%
  bind_rows(.id = "id") %>%
  filter(term=="trajectory") %>% 
  mutate(padj = p.adjust(p.value, method = "BH"))
write_csv(a, "data/stat-testing-TFs-in-ccr7-Trajectory-new.csv")

#pseudotemporal ordering - ifn_cl
n <- cellsfromtree(ltr,ifn_cl)
x <- getfdata(ltr@sc)


fs  <- filterset(x,n=n$f, minexpr = 1, minnumber = 1)
if (!file.exists("data/s1d-ifn_cl-new.Robj")) {
  s1d <- getsom(fs,nb=1000,alpha=.5)
  save(s1d, file = "data/s1d-ifn_cl-new.Robj")
} else {
  load("data/s1d-ifn_cl.Robj")
}

ps  <- procsom(s1d,corthr=.85,minsom=3)
y    <- ltr@sc@cpart[n$f]
fcol <- ltr@sc@fcol
plotheatmap(ps$all_singl_ea.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)

pdf(paste0('plots/heatmaps/', date, '-ifn_cl-trajectory-heatmap.pdf'), width = 8.57, height = 5.79)
plotheatmap(ps$all_singl_ea.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
dev.off()

png(paste0('plots/heatmaps/', date, '-ifn_cl-trajectory-heatmap.png'), width = 8.57, height = 5.79, res = 600, units = "in")
plotheatmap(ps$all_singl_ea.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
dev.off()

#eport node genes
modules <- data.frame('Node' = NA, 'Genes' = NA)
for (i in 1: max(ps$nodes)) {
  gene_names <- names(ps$nodes)[ps$nodes == i]
  gene_names <- gsub('_.*', '', gene_names)
  modules2 <- data.frame('Node' = NA, 'Genes' = gene_names)
  modules2$Node <- rep(as.character(i), nrow(modules2))
  modules <- rbind(na.omit(modules), modules2)
}

write.csv(modules, paste0('data/',date,'-nodes-stemid-vector-ifn_cl.csv'))

#find TFs in node genes
tfs <- read_csv("http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv")

tfs_found <- modules[modules$Genes %in% tfs[[3]][tfs$`Binding mode` != "Not a DNA binding protein"],]

#export TFs
write_csv(tfs_found, paste0('data/',date,'-nodes-stemid-TFs-ifn_cl.csv'))

#plot fate bias

