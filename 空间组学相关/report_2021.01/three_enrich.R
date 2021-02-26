cat('fit <-  try(three_enrich(genelist = gene_diff,name = paste0(GO,"ST_cluster",n),method = c("BP","MF","CC"),type="mmu",universe=gene_all))','\n',
             'if("try-error" %in% class(fit)){next}else{print(paste0("cluster_",n,"_enrich_finished"))')

three_enrich<-function(genelist,name,method=c("BP","MF","CC"),type=c("hsa","mmu"),universe,
                       showCategory){
  library(clusterProfiler)
  library(ReactomePA)
  library(enrichplot)
  library(S4Vectors)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)

if (type=="hsa") {
    go <- org.Hs.eg.db
    rec <- "human"}else{
      go <- org.Mm.eg.db
      rec <-  "mouse"
    }
  
  kegg_enrich<- enrichKEGG(
    gene        = genelist,
    organism    = type,#'mmu', #
    universe   = gene_all,
    pvalueCutoff = 0.05,
    qvalueCutoff =0.2,
    minGSSize = 5,
    maxGSSize = 1000,
  )
  fit <- try(isEmpty(kegg_enrich@result))
  if("try-error" %in% class(fit)){kegg_enrich_df <- data.frame()}else{
    kegg_enrich@result$Description <- paste0("Kegg: ",kegg_enrich@result$Description)
    kegg_enrich_df <- kegg_enrich@result}
  
  
  go_enrich_df <- data.frame()
  for (i in method) { 
    go_enrich <- enrichGO(gene= genelist,
                          universe      = gene_all,
                          OrgDb         = go,
                          ont           = i,
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.2,
                          minGSSize = 5,
                          maxGSSize = 1000,
                          readable      = TRUE
    )
    fit <- try(isEmpty(go_enrich@result))
    if("try-error" %in% class(fit)){go_enrich_df_tmp <- data.frame()}else{
      go_enrich@result$Description <- paste0("Go_",i,": ",go_enrich@result$Description)
      go_enrich_df_tmp <- go_enrich@result}
    go_enrich_df <- rbind(go_enrich_df,go_enrich_df_tmp)
  }
  
  rec_enrich <- enrichPathway(gene  =  genelist,
                              universe = gene_all,
                              organism = rec,
                              pAdjustMethod = "BH", 
                              qvalueCutoff = 0.2,
                              pvalueCutoff=0.05, 
                              minGSSize = 5,
                              maxGSSize = 1000, 
                              readable = FALSE)
  if("try-error" %in% class(fit)){rec_enrich_df <- data.frame()}else{
    rec_enrich@result$Description <- paste0("Reactome: ",rec_enrich@result$Description)
    rec_enrich_df <- rec_enrich@result}
  
  tmp <- rbind(go_enrich_df,kegg_enrich_df,rec_enrich_df)
  write.csv(tmp,paste0(name,"_enriched.csv"))
  go_enrich@result <- tmp
  
  pdf(paste0(name,"_dotplot.pdf"),width = 12)
  print(dotplot(go_enrich,showCategory = showCategory))
  dev.off()
  
  pdf(paste0(name,"_barplot.pdf"), width = 12)
  print(barplot(go_enrich,showCategory = showCategory))
  dev.off()
}

