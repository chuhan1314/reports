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
rm(list=ls())

tmp=read.csv("./spotlight_res.csv",row.names = 1)
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


counts_foucs <- cortex@assays$Spatial@data
data.input <- as.matrix(counts_foucs)
lables <- tmp$spot.pre
identity <- data.frame(group = lables,row.names = rownames(cortex@meta.data))

cellchat <- createCellChat(data = data.input)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
#levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 2) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)

cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
levels(cellchat@idents)
head(cellchat@LR$LRsig)
cellchat@netP$pathways
save(counts_foucs,cellchat,file = "e13_cellchat.rda")


