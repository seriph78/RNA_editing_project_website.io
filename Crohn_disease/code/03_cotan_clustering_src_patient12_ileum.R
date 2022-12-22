devtools::load_all("../COTAN/")
library(stringr)

out_dir = "Crohn_disease/data/03_clustred_data/"

sample.id = "patient12_ileum"

obj <- readRDS(paste0("Crohn_disease/data/02_cleaned_data/",sample.id,"_obj_cotan.RDS"))

clusters <- cellsUniformClustering(obj,
                                   cores = 12
                                   )

obj <- addClusterization(obj,clName = "first_round", clusters = clusters)

saveRDS(obj, file = paste0(out_dir,sample.id,".cotan.clustered.RDS") )


list[clusters.coex,pval.table] <- DEA_on_clusters(obj)
obj <- addClusterizationCoex(obj,clName = "first_round",coexDF = clusters.coex)


saveRDS(pval.table, paste(out_dir,"pval.",sample.id,".clustering.RDS"))
saveRDS(obj, file = paste0(out_dir,sample.id,".cotan.clustered.RDS") )

obj <- merge_cell.clusters(obj = obj,
                           cond = sample.id,
                           cores=15, 
                           #srat = paste0("Seurat_obj_",sample.id,"_with_cotan_clusters.RDS"),
                           out_dir = paste0(out_dir) ,
                           GEO=getMetadataDataset(obj)[getMetadataDataset(obj)[,1]=="GEO:",2],
                           sc.method="10X"
                           )



saveRDS(obj, file = paste0(out_dir,sample.id,".clustered.merged.cotan.RDS") )

out_dir <- "Crohn_disease/data/04_cluster_identification/"

expression.cl <- get.clusters.delta.expression(obj)
write.csv(expression.cl,paste0(out_dir,sample.id,"_expression_clusters.csv"))
