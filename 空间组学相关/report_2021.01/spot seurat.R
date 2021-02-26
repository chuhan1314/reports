#BiocManager::install("rhdf5")
library(rhdf5)
library(Seurat)
library(dplyr)

if (!file.exists("ref.rda")) {
#count<- h5read("../GSE164430_RAW/GSM5009539_10X_mouse_down_processed.h5ad/GSM5009539_10X_mouse_down_processed.h5ad","X")
barcode <-  h5read("./GSE164430_RAW/GSM5009539_10X_mouse_down_processed.h5ad/GSM5009539_10X_mouse_down_processed.h5ad","obs")
gene <-  h5read("./GSE164430_RAW/GSM5009539_10X_mouse_down_processed.h5ad/GSM5009539_10X_mouse_down_processed.h5ad","var")

#count <- as.data.frame(count)
gene.phe <- as.data.frame(cbind(gene=gene[["_index"]],ens=gene[["gene_ids"]]))
meta <- as.data.frame(unlist(cbind(id=barcode[["_index"]],celltype=barcode[["celltype"]],percent.mt=barcode[["percent_mito"]])))
#rownames(count) <- gene.phe$gene
#colnames(count) <- meta$id
rownames(meta) <- meta$id
####################字典？？？？？？#######################################
dict <- eval(parse(text=paste0("c(",(paste(paste0('"',as.character(c(0:9)),'"','= "',barcode[["__categories"]]$celltype,'"'),collapse = ",")),")")))
meta$type <- unlist(lapply(meta$celltype,function(x){return(dict[x])}))
##########################################################################
count <- Read10X_h5("./GSE164430_RAW/GSM5009539_10X_raw_feature_bc_matrix.h5")
ref_pre <- CreateSeuratObject(counts = count)
ref_pre <- subset(ref_pre,cells=rownames(meta))
ref_pre <- as.data.frame(ref_pre@assays$RNA@counts)
ref_pre <- ref_pre[-grep(rownames(ref_pre),pattern = "^GRCh"),]
rownames(ref_pre) <- unlist(lapply(strsplit(rownames(ref_pre),"---"),function(x){return(x[2])}))

ref <- CreateSeuratObject(counts = ref_pre)
ref@meta.data <- cbind(ref@meta.data,meta)
ref <- NormalizeData(ref, normalization.method = "LogNormalize", scale.factor = 10000)
ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = 3000)
ref <- ScaleData(ref, features = VariableFeatures(ref))
ref <-  RunPCA(ref,features = rownames(ref))
ElbowPlot(ref,ndims = 50)
ref <- RunUMAP(ref, dims = 1:30)
DimPlot(ref, group.by = "type", label = TRUE)
save(ref,file = "ref.rda")}else{load("ref.rda")}

cortex =  Load10X_Spatial(file.path(paste0("E:/data/work/kongjianshujv/e13_beirui/out")),slice = "e13")
cortex <- NormalizeData(cortex, normalization.method = "LogNormalize", scale.factor = 10000)
cortex <- FindVariableFeatures(cortex, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(cortex)
cortex <- ScaleData(cortex, features = all.genes)
cortex <- RunPCA(cortex, features = VariableFeatures(object = cortex), npcs = 50)
ElbowPlot(cortex,ndims = 50)
cortex <- FindNeighbors(cortex, dims = 1:50) 
cortex <- RunUMAP(cortex, dims = 1:50,seed.use = 20)
cortex <- FindClusters(cortex, resolution = 0.8)
SpatialDimPlot(cortex)

anchors <- FindTransferAnchors(reference = ref, query = cortex, normalization.method = "LogNormalize")
predictions.assay <- TransferData(anchorset = anchors, refdata = ref$type, prediction.assay = TRUE)
cortex[["predictions"]] <- predictions.assay

DefaultAssay(cortex) <- "predictions"
#SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)
cortex <- FindSpatiallyVariableFeatures(cortex, assay = "predictions", selection.method = "moransi", 
                                        features = rownames(cortex), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(cortex), 6)
SpatialPlot(object = cortex, features = top.clusters, ncol = 2)
SpatialFeaturePlot(cortex, features = c("Kupffer"), crop = T)
