
#library

devtools::load_all("../COTAN/")

#parameters

root <- "Crohn_disease/"


for (sample.id in c("patient5_PBMC",  "patient5_Uninvolved","patient5_Involved")) { #
  out_dir <- paste0(root,"data/02_cleaned_data/")
  obj <- readRDS(paste0(out_dir,sample.id,"_obj_cotan.RDS"))
  obj <- estimateDispersionBisection(obj, cores = 12)
  gc()
  obj <- calculateCoex(obj)
  
  saveRDS(obj,paste0(out_dir,sample.id,"_obj_cotan.RDS"))
  
  clusters <- cellsUniformClustering(obj,
                                     cores = 12
  )
  
  obj <- addClusterization(obj,clName = "first_round", clusters = clusters)
  
  list[clusters.coex,pval.table] <- DEAOnClusters(obj)
  
  obj <- addClusterizationCoex(obj,clName = "first_round",coexDF = clusters.coex)
  
  
  saveRDS(pval.table, paste(out_dir,"pval.",sample.id,".clustering.RDS"))
  saveRDS(obj, file = paste0(out_dir,sample.id,".cotan.clustered.RDS") )
  
  #obj <- readRDS(file = paste0(out_dir,sample.id,".cotan.clustered.RDS") )
  
  obj <- merge_cell.clusters(obj = obj,
                             cond = sample.id,
                             cores=15, 
                             out_dir = out_dir ,
                             GEO=getMetadataDataset(obj)[getMetadataDataset(obj)[,1]=="GEO:",2],
                             sc.method="10X"
  )
  
  saveRDS(obj, file = paste0(out_dir,sample.id,".clustered.merged.cotan.RDS") )
  
  out_dir <- paste0(root,"data/04_cluster_identification/")
  
  expression.cl <- get.clusters.delta.expression(obj)
  write.csv(expression.cl,paste0(out_dir,sample.id,"_expression_clusters.csv"))
  
}

