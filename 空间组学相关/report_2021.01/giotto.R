library(Seurat)
library(dplyr)
library(Giotto)
meta.seurat <- read.csv("/data/e13/meta.csv")

data_path ="/data/input/"
workdir="./"
expr_data_path=fs::path(data_path, "raw_feature_bc_matrix")
raw_matrix=get10Xmatrix(path_to_data=expr_data_path, gene_column_index=2)
spatial_locations=data.table::fread(fs::path(data_path, "spatial", "tissue_positions_list.csv"))
spatial_locations = spatial_locations[match(colnames(raw_matrix), V1)]
colnames(spatial_locations) = c('barcode', 'in_tissue', 'array_row', 'array_col', 'col_pxl', 'row_pxl')
myinst=createGiottoInstructions(save_plot=T, show_plot=F, save_dir=workdir, python_path="/usr/bin/python3")
e13 <- createGiottoObject(raw_exprs = raw_matrix, spatial_locs = spatial_locations[,.(row_pxl,-col_pxl)], 
                          instructions = myinst, cell_metadata = spatial_locations[,.(in_tissue, array_row, array_col)])
spatPlot(gobject = e13,  point_size = 2, cell_color = 'in_tissue', 
         cell_color_code = c('0' = 'lightgrey', '1' = 'blue'), save_param=c(save_name="1-spatplot"))

e13 = subsetGiotto(e13, cell_ids = e13@cell_metadata[in_tissue==1]$cell_ID)
# ?? QC ??e13 <- filterGiotto(gobject = e13, expression_threshold = 1, gene_det_in_min_cells = 50, min_det_genes_per_cell = 1000,expression_values = c('raw'),verbose = T)
#spatPlot(gobject = e13)

e13@cell_metadata$seurat.cluster <- meta.seurat$seurat_clusters
e13 <- normalizeGiotto(gobject = e13, scalefactor = 6000, verbose = T)
e13 <- addStatistics(gobject = e13)
#spatPlot(gobject = e13,  point_size = 2, save_param=c(save_name="2-spatplot"))
#spatPlot(gobject = e13, cell_color = 'nr_genes', color_as_factor = F,  point_size = 2, save_param=c(save_name="3-spatplot"))
e13 <- calculateHVG(e13)
featgenes = fDataDT(e13)[hvg == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$gene_ID
e13 <- runPCA(gobject = e13, genes_to_use = featgenes, scale_unit = F, center=T, method="factominer")
signPCA(e13, genes_to_use = featgenes, scale_unit = F)
plotPCA(gobject = e13)
e13 <- runUMAP(e13, dimensions_to_use = 1:10)
plotUMAP(gobject = e13)
e13 <- runtSNE(e13, dimensions_to_use = 1:10)
plotTSNE(gobject = e13)
e13 <- createNearestNetwork(gobject = e13, dimensions_to_use = 1:10, k = 30)
e13 <- doLeidenCluster(gobject = e13, resolution = 0.5, n_iterations = 1000)
plotUMAP(gobject = e13, cell_color = 'leiden_clus', show_NN_network = T, point_size = 2)
spatDimPlot(gobject = e13, cell_color = 'leiden_clus',dim_point_size = 1.5, spat_point_size = 1.5)
spatPlot(gobject = e13, cell_color = 'leiden_clus',point_size = 8)
spatPlot(gobject = e13, cell_color = 'seurat.cluster',point_size = 8)
#DG_subset = subsetGiottoLocs(e13, x_max = 9000, x_min = 6000, y_max = -7000, y_min = -10000, return_gobject = T)
#spatDimPlot(gobject = DG_subset, cell_color = 'leiden_clus', spat_point_size = 5)
gini_markers_subclusters = findMarkers_one_vs_all(gobject = e13,method = 'gini',
                                                  expression_values = 'normalized',cluster_column = 'leiden_clus',min_genes = 20,min_expr_gini_score = 0.5,min_det_gini_score = 0.5)
topgenes_gini = gini_markers_subclusters[, head(.SD, 10), by = 'cluster']$genes
violinPlot(e13, genes = unique(topgenes_gini), cluster_column = 'leiden_clus',strip_text = 8, strip_position = 'right')
topgenes_gini = gini_markers_subclusters[, head(.SD, 2), by = 'cluster']$genes
my_cluster_order = c(1:11)
plotMetaDataHeatmap(e13, selected_genes = topgenes_gini, 
                    custom_cluster_order = my_cluster_order, metadata_cols = c('leiden_clus'), x_text_size = 10, y_text_size = 10)
scran_markers_subclusters = findMarkers_one_vs_all(gobject = e13,method = 'scran',expression_values = 'normalized',cluster_column = 'leiden_clus')
# violinplot
topgenes_scran = scran_markers_subclusters[, head(.SD, 1), by = 'cluster']$genes
violinPlot(e13, genes = unique(topgenes_scran), cluster_column = 'leiden_clus',strip_text = 8, strip_position = 'right')

topgenes_scran = scran_markers_subclusters[, head(.SD, 2), by = 'cluster']$genes
plotMetaDataHeatmap(e13, selected_genes = topgenes_scran, custom_cluster_order = my_cluster_order,metadata_cols = c('leiden_clus'))
# umap plots
dimGenePlot2D(e13, expression_values = 'normalized',genes = scran_markers_subclusters[, head(.SD, 1), by = 'cluster']$genes,cow_n_col = 3, point_size = 1)

e13_sc_markers = read.csv( '/data/e13/sig_matrix.csv',row.names = 2)
sig_matrix <- e13_sc_markers[,-1]


#############PAGE####################
e13 = runSpatialEnrich(e13, sign_matrix = sig_matrix, enrich_method = 'PAGE')
e13 = runSpatialEnrich(e13, sign_matrix = sig_matrix, enrich_method = 'rank')
e13 = runSpatialEnrich(e13, sign_matrix = sig_matrix, enrich_method = 'hypergeometric')
tmp <- 'hypergeometric'
cell_types = colnames(sig_matrix)
plotMetaDataCellsHeatmap(gobject = e13,metadata_cols = 'leiden_clus',
                         value_cols = cell_types,spat_enr_names = tmp,x_text_size = 8, y_text_size = 8)
cell_types_subset = colnames(sig_matrix)[1:10]
spatCellPlot(gobject = e13, spat_enr_names = tmp,cell_annotation_values = cell_types_subset,
             cow_n_col = 4,coord_fix_ratio = NULL, point_size = 2,cell_color_gradient=c("blue","white","red"))
spatCellPlot(gobject = visium_brain, spat_enr_names = tmp ,cell_annotation_values = c("MC38", "Myeloid", "T.Cell", "Kupffer","LSEC"),
             cow_n_col = 2, point_size = 2)
for (i in c('PAGE','rank','hypergeometric')) {
  tmp <- e13@spatial_enrichment[[i]]
  write.csv(tmp,file = paste0(i,"gio_enrich.csv"))
}
##################rank#########################

colnames(PAGEgio_enrich) <- paste0("page,",colnames(PAGEgio_enrich))
colnames(rankgio_enrich) <- paste0("rank,",colnames(rankgio_enrich))
colnames(hypergeometricgio_enrich) <- paste0("hyper,",colnames(hypergeometricgio_enrich))
table(rownames(cortex@meta.data)==rankgio_enrich$`rank,cell_ID`)
cortex@meta.data <- cbind(cortex@meta.data,PAGEgio_enrich[,3:12],rankgio_enrich[,3:12],hypergeometricgio_enrich[,3:12])
name <- colnames(cortex@meta.data)
for (i in 18:47) {
  #i=18
  x <- name[i]
  pdf(paste0(i,".pdf"))
  print(SpatialFeaturePlot(cortex,features = x))
  dev.off()
}
