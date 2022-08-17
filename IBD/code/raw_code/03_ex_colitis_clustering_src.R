devtools::load_all("../../COTAN/")
out_dir = "data/IBD/clustering/"

obj <- cell_homogeneous_clustering(cond = "colitis",out_dir = out_dir, 
                                   in_dir = "data/IBD/check_clean/COTAN_clean/", 
                                   cores = 12, 
                                   dataset_type = "COTAN", dataset_name = "colitis.cotan.RDS", 
                                   GEO="GSE121380", sc.method="10X 5'"
)

saveRDS(obj, file = file.path("data/IBD/clustering/","colitis.cotan.clustered.RDS") )


gc()
temp <- DEA_on_clusters(obj)
obj <- temp[[1]]
pval <- temp[[2]]

saveRDS(pval, paste(out_dir,"pval.colitis.clustering.RDS"))
saveRDS(obj, file = file.path("data/IBD/clustering","colitis.cotan.clustered.RDS") )
        
obj <- merge_cell.clusters(obj = obj,cond = "colitis",cores=10, 
                           out_dir_root = out_dir,
                           srat = "Seurat_obj_colitis_with_cotan_clusters.RDS",
                           out_dir = out_dir ,
                           GEO="GSE121380", sc.method="10X 5'"
                           #mt = TRUE, mt_prefix="^MT"
                           )

saveRDS(obj, file = file.path("data/IBD/clustering","colitis.cotan.clustered.merged.RDS") )
