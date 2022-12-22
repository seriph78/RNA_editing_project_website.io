
#library

devtools::load_all("../COTAN/")

#parameters

root <- "Crohn_disease/"
out_dir <- paste0(root,"data/02_cleaned_data/")

for (sample.id in c("patient12_Uninvolved", "patient12_Involved", "PBMC")) {
  obj <- readRDS(paste0(out_dir,sample.id,"_obj_cotan.RDS"))
  
  #obj <- cotan_analysis(obj,cores = 12)
  obj <- estimateDispersionBisection(obj, cores = 12)
  gc()
  obj <- calculateCoex(obj)
  
  saveRDS(obj,paste0(out_dir,sample.id,"_obj_cotan.RDS"))
}

