devtools::load_all("../COTAN/")
library(stringr)

out_dir = "Crohn_disease/data/03_clustred_data/"

sample.id = "patient5_ileum"
GEO = "GSM3972010_GSM3972009"
gc()
obj <- cell_homogeneous_clustering(cond = sample.id,
                                   out_dir = out_dir, 
                                   in_dir = "Crohn_disease/data/02_cleaned_data/", 
                                   cores = 12, 
                                   dataset_name = "patient5_ileum_obj_cotan.RDS", 
                                   GEO=GEO, 
                                   sc.method="10X"
)

#saveRDS(obj, file = paste0(out_dir,sample.id,".cotan.clustered.RDS") )


gc()
temp <- DEA_on_clusters(obj)
obj <- temp[[1]]
pval <- temp[[2]]

saveRDS(pval, paste(out_dir,"pval.",sample.id,".clustering.RDS"))
saveRDS(obj, file = paste0(out_dir,sample.id,".cotan.clustered.RDS") )
        
obj <- merge_cell.clusters(obj = obj,cond = sample.id,cores=10, 
                           #out_dir_root = out_dir,
                           srat = paste0("Seurat_obj_",sample.id,"_with_cotan_clusters.RDS"),
                           out_dir = out_dir ,
                           GEO=GEO, sc.method="10X"
                           )

saveRDS(obj, file = paste0(out_dir,sample.id,".clustered.merged.cotan.RDS") )

out_dir <- "Crohn_disease/data/04_cluster_identification/"

expression.cl <- get.clusters.delta.expression(obj)
write.csv(expression.cl,paste0(out_dir,sample.id,"_expression_clusters.csv"))
