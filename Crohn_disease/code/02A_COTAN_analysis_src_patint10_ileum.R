
#library

devtools::load_all("../COTAN/")

#parameters

root <- "Crohn_disease/"
out_dir <- paste0(root,"data/02_cleaned_data/")

sample.id = "patient10_ileum"#"sample69"

###code
obj <- readRDS(paste0(out_dir,sample.id,"_obj_cotan.RDS"))

obj <- cotan_analysis(obj,cores = 12)

saveRDS(obj,paste0(out_dir,sample.id,"_obj_cotan.RDS"))

#gc()

#obj = get.coex(obj)


#saveRDS(obj,paste0(out_dir,sample.id,"_obj_cotan.RDS"))
