rm(list=ls())
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SingleR)
library(R.matlab)
library(pheatmap)
library(monocle)
library(CellChat)
library(ggalluvial)
library(svglite)
library(NMF)
require(data.table)
library(knitr)
library(cowplot)
library(RColorBrewer)
library(org.Mm.eg.db)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GSEABase)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 4000 * 1024^2)


if (T) {
  wd <- "E:/data/work/xinqiao/2021.01/zhureport/"
  setwd(wd)
  if (!file.exists("merge")) {dir.create("merge")}
  wd <- paste0(wd,"merge/")
  setwd(wd)
  if (!file.exists("out")) {dir.create("out")}
  out <- paste0(wd,"out/")
  if (!file.exists("DE")) {dir.create("DE")}
  DE <- paste0(wd,"DE/")
  if (!file.exists("go_kegg")){dir.create("go_kegg")}
  GO <- paste0(wd,"go_kegg/")
  if (!file.exists("gsea")){dir.create("gsea")}
  gsea <- paste0(wd,"gsea/")
}

################have_mt_spot#######################
################have_mt_spot#######################
folders=list.files('E:/data/work/kongjianshujv/',pattern='_beirui$')
#positon <- c("edge","middle","core","normal","normal","edge","middle","core","core","dege","middle","normal")
#sample <- c(rep("sample3",4),rep("sample4",4),rep("sample5",4))
project <- c("D8","E13")

sceList = lapply(folders,function(folder){ 
  x <- Load10X_Spatial(file.path(paste0("E:/data/work/kongjianshujv/",folder,"/out")),slice = folder)
  return(x)
})

for (i in 1:length(sceList)) {
  sceList[[i]]@meta.data$orig.ident = project[i]
}

for (i in 1:length(sceList)) {
  sceList[[i]] <- SCTransform(sceList[[i]], assay = "Spatial", verbose = FALSE)
  sceList[[i]]@meta.data$orig.ident = project[i]
}

pancreas.features <- SelectIntegrationFeatures(object.list = sceList, nfeatures = 3000)
pancreas.list <- PrepSCTIntegration(object.list = sceList, anchor.features = pancreas.features, 
                                    verbose = FALSE)
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT", 
                                           anchor.features = pancreas.features, verbose = FALSE)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)

DefaultAssay(pancreas.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 50, verbose = FALSE)
ElbowPlot(pancreas.integrated,ndims = 50)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
pancreas.integrated <- FindNeighbors(pancreas.integrated, dims = 1:30)
pancreas.integrated <- FindClusters(pancreas.integrated, verbose = FALSE,resolution = 0.2)

DimPlot(pancreas.integrated, reduction = "umap", group.by = c("ident"))
DimPlot(pancreas.integrated, reduction = "umap", group.by = c("orig.ident"))
SpatialDimPlot(pancreas.integrated)
###############################################################
###############################################################
#m=1
de_markers <- FindAllMarkers(pancreas.integrated)
x <- de_markers %>% group_by(cluster)  %>% top_n(2,avg_logFC)
for (i in 1:nrow(x)) {
  plot <- SpatialFeaturePlot(object = pancreas.integrated, features = x$gene[i],ncol = 1)
  #pdf(paste0("cluster_",x$cluster[i],"_",i,".pdf"))
  print(plot)
  #dev.off()
}
#SpatialFeaturePlot(object = pancreas.integrated, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)
top10 <- de_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
Doheatmap_chuhan(pancreas.integrated,top10)

source("E://data/work/xinqiao/2021.01/zhureport/three_enrich.R")
gene_all_pre <- rownames(pancreas.integrated)
gene_all_id <-  bitr(gene_all_pre, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Mm.eg.db)
gene_all <- gene_all_id$ENTREZID
for (n in unique(de_markers$cluster)) {
  #n=1
  tmp <- rownames(de_markers[de_markers$cluster == n,])
  tmp2 <-  lapply(strsplit(tmp,split = "[.]"), function(x){
    return(x[1]) 
  })
  tmp2 <-  unlist(tmp2)
  id <-  bitr(tmp2,
              fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = org.Mm.eg.db)
  gene_diff <- id$ENTREZID

  fit <-  try(three_enrich(genelist = gene_diff,name = paste0(GO,"ST_cluster",n),
                           method = c("BP","MF","CC"),type="mmu",universe = "",showCategory=15))
  if("try-error" %in% class(fit)){next}else{print(paste0("cluster_",n,"_enrich_finished"))
  }
}

