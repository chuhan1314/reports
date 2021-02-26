x=read.csv("./spotlight_res.csv",row.names = 1)
y=x[,-which(colnames(x) %in% c("orig.ident","nCount_Spatial",
                         "nFeature_Spatial","Spatial_snn_res.0.5","seurat_clusters","res_ss"))]
tmp <- c(paste0("prediction.score.",gsub("[.]","_",colnames(y)))[c(-ncol(y),-(ncol(y)-1))],"predicted.id","prediction.score.max")
colnames(y) <- tmp
write.table(y,file = "label_transfer_bc.csv",sep = "\t",quote = F)
