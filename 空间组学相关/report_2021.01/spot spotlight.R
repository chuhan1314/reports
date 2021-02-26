#devtools::install_github("https://github.com/MarcElosua/SPOTlight")
library(SPOTlight)
library(Seurat)
library(dplyr)
#ref
load("ref.rda")
cortex_sc <- ref

cortex =  Load10X_Spatial(file.path(paste0("E:/data/work/kongjianshujv/e13_beirui/out")),slice = "e13")
cortex <- NormalizeData(cortex, normalization.method = "LogNormalize", scale.factor = 10000)
cortex <- FindVariableFeatures(cortex, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(cortex)
cortex <- ScaleData(cortex, features = all.genes)
cortex <- RunPCA(cortex, features = VariableFeatures(object = cortex), npcs = 50)
ElbowPlot(cortex,ndims = 50)
cortex <- FindNeighbors(cortex, dims = 1:50) 
cortex <- RunUMAP(cortex, dims = 1:50,seed.use = 20)
cortex <- FindClusters(cortex, resolution = 0.5)
SpatialDimPlot(cortex,label = T)

library(ggplot2)
x <- cortex@images$e13@coordinates
ggplot2::ggplot(x,aes(x=imagecol,y=imagerow))+
  geom_point()+
  geom_text(aes(rownames(x)))

Seurat::Idents(object = cortex_sc) <- cortex_sc@meta.data$type
cluster_markers_all <- Seurat::FindAllMarkers(object = cortex_sc, 
                                              slot = "data",
                                              verbose = TRUE, 
                                              only.pos = TRUE, 
                                              logfc.threshold = 1,
                                              min.pct = 0.9)

set.seed(123)
spotlight_ls <- spotlight_deconvolution(se_sc = cortex_sc,
                                        counts_spatial = cortex@assays$Spatial@counts,
                                        clust_vr = "type",
                                        cluster_markers = cluster_markers_all,
                                        cl_n = 100,
                                        hvg = 3000,
                                        ntop = NULL,
                                        transf = "uv",
                                        method = "nsNMF",
                                        min_cont = 0.09)

decon_mtrx <- spotlight_ls[[2]]
cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]
cortex@meta.data <- cbind(cortex@meta.data, decon_mtrx)
cortex@meta.data$spot.pre <- apply(decon_mtrx,1,function(x){tmp <- names(which.max(x))
return(tmp)})
cortex@meta.data$score.max <- apply(decon_mtrx,1,function(x){tmp <- x[which.max(x)]
return(tmp)})
write.csv(cortex@meta.data,"spotlight_res.csv")
Seurat::SpatialDimPlot(cortex,group.by = "spot.pre")

SPOTlight::spatial_scatterpie(se_obj = cortex,
                              cell_types_all = cell_types_all,
                              img_path = "E:/data/work/kongjianshujv/e13_beirui/out/spatial/tissue_lowres_image.png",
                              pie_scale = 0.45)
SPOTlight::spatial_scatterpie(se_obj = cortex,
                              cell_types_all = cell_types_all,
                              img_path = "E:/data/work/kongjianshujv/e13_beirui/out/spatial/tissue_lowres_image.png",
                              pie_scale = 1,
                              cell_types_interest = "MC38")
SpatialFeaturePlot(cortex,
                   features = "T.Cell", 
                   pt.size.factor = 1.6,
                   alpha = c(0, 1)) +
  ggplot2::scale_fill_gradientn(
    colours = heat.colors(10, rev = TRUE),
    limits = c(0, 1))

graph_ntw <- get_spatial_interaction_graph(decon_mtrx = decon_mtrx[, cell_types_all])
library(igraph)
library(RColorBrewer)
deg <- degree(graph_ntw, mode="all")
# Get color palette for difusion
edge_importance <- E(graph_ntw)$importance
# Select a continuous palette
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'seq',]
# Create a color palette
getPalette <- colorRampPalette(brewer.pal(9, "YlOrRd"))
# Get how many values we need
grad_edge <- seq(0, max(edge_importance), 0.1)
# Generate extended gradient palette dataframe
graph_col_df <- data.frame(value = as.character(grad_edge),
                           color = getPalette(length(grad_edge)),
                           stringsAsFactors = FALSE)
# Assign color to each edge
color_edge <- data.frame(value = as.character(round(edge_importance, 1)), stringsAsFactors = FALSE) %>%
  dplyr::left_join(graph_col_df, by = "value") %>%
  dplyr::pull(color)
# Open a pdf file
plot(graph_ntw,
     # Size of the edge
     edge.width = edge_importance,
     edge.color = color_edge,
     # Size of the buble
     vertex.size = deg,
     vertex.color = "#cde394",
     vertex.frame.color = "white",
     vertex.label.color = "black",
     vertex.label.family = "Ubuntu", # Font family of the label (e.g.“Times”, “Helvetica”)
     layout = layout.circle)

nmf_mod_ls <- spotlight_ls[[1]]
nmf_mod <- nmf_mod_ls[[1]]
h <- NMF::coef(nmf_mod)
rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
topic_profile_plts <- dot_plot_profiles_fun(h = h,
                                            train_cell_clust = nmf_mod_ls[[2]])
topic_profile_plts[[2]] + theme(axis.text.x = element_text(angle = 90), 
                                axis.text = element_text(size = 12))
