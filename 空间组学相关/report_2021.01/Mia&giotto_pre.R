library(dplyr)
library(reshape2)
load("ref.rda")

DimPlot(ref,group.by = "type")
Idents(ref) <- ref$type
deg <- FindAllMarkers(ref)
deg.re <- deg %>% filter(p_val_adj < 0.0001,abs(avg_logFC) > 1)
table(deg.re$cluster)
table(duplicated(deg.re$gene))
 x <- data.frame(cluster=deg.re$cluster,gene=deg.re$gene,v=1)
y <- dcast(x,gene~cluster,value.var = "v",length) 

write.csv(y,file = "sig_matrix.csv")

