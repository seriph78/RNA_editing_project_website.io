devtools::load_all("../COTAN/")

setwd("asuna/4.RNA_Editing_AI/")
out_dir = "data/IBD/clustering/"

obj <- cell_homogeneous_clustering(cond = "control",out_dir = out_dir, 
                                   in_dir = "data/IBD/check_clean/COTAN_clean/", 
                                   cores = 12, 
                                   dataset_type = "COTAN", dataset_name = "control.cotan.RDS", 
                                   GEO="GSE121380", sc.method="10X 5'"
)

#saveRDS(obj, file = file.path("data/IBD/clustering/","control.cotan.clustered.RDS") )


gc()
temp <- DEA_on_clusters(obj)
obj <- temp[[1]]
pval <- temp[[2]]

saveRDS(pval, paste(out_dir,"pval.control.clustering.RDS"))
saveRDS(obj, file = file.path("data/IBD/clustering","control.cotan.clustered.RDS") )
        
obj <- merge_cell.clusters(obj = obj,cond = "control",cores=10, 
                           out_dir_root = out_dir,
                           srat = "Seurat_obj_control_with_cotan_clusters.RDS",
                           out_dir = out_dir ,
                           GEO="GSE121380", sc.method="10X 5'"
                           #mt = TRUE, mt_prefix="^MT"
                           )

saveRDS(obj, file = file.path("data/IBD/clustering","control.cotan.clustered.merged.RDS") )
