out_dir <- "IBD/data/04_cluster_identification/"
devtools::load_all("../../COTAN/") 
 
cond <- "colitis"

obj <- readRDS("IBD/data/03_clustered_data/colitis.cotan.clustered.merged.RDS")

# cluster assignment with expression increment
expression.cl <- get.clusters.delta.expression(obj)

write.csv(expression.cl,paste0(out_dir,cond,"_expression_clusters.csv"))

markers <- list("B cells" = c("MS4A1", "CD83", "CD79A", "VPREB3", "CD37",
                              "BANK1", "LY9" ),
                "Plasma" = c( "JCHAIN", "IGHA2", "IGHA1", "IGKC"),
                "T & NK cells"= c("CCL5", "CD7", "CD3E", "CD3D", "IL32", 
                                  "NKG7", "CST7", "FYB", "CD8A", "CD2", 
                                  "HCST"),
                "LTi"=c( "IL22", "CSF2", "IL4I1", "XCL1", "KIT", "BHLHE40",
                         "TNFRSF18", "KRT86", "NCOA7", "IL23R", "AQP3", 
                         "ALDOC", "TNFSF13B"), 
                "Myeloid"= c("LYZ", "CST3", "IL1B", "TPSAB1", "CXCL8", 
                             "CPA3", "TPSB2", "AIF1", "FCER1G", "C1QC", 
                             "C1QA", "C1QB", "S100A8"), 
                "Mitotic B & T cells"= c("STMN1","HMGB2", "HMGN2", "HIST1H4C",
                                         "H2AFZ","TUBB", "TUBA1B", "HMGB1", "UBE2C",                                     "RRM2", "PTTG1", "DUT", "TYMS", "H2AFV",
                                         "MKI67", "TOP2A", "CENPF", "NUSAP1"), 
                
                "Fibro & Endo"= c( "CXCL14", "ADAMDEC1", "COL3A1", "APOE",
                                   "CFD", "IGFBP7", "CCL2", "COL1A2", "COL1A1",
                                   "A2M", "COL6A2", "IFITM3", "LUM", "DCN",
                                   "MFAP4"),
                "Epithelial" = c( "PHGR1", "TFF3", "FABP1", "ELF3", "ZG16",
                                  "KRT8", "SPINK4", "LGALS4", "KRT18", "CLDN4",
                                  "PIGR", "S100A6", "FCGBP", "KRT19", "CLDN3")
)

df <- as.data.frame(matrix(nrow = length(names(markers)),ncol = ncol(obj@cluster_data)))
rownames(df) <- names(markers)
colnames(df) <- colnames(obj@cluster_data)
teta <- -1/0.1 * (log(0.25))
#not_ass_clusters <- NA
for (ro in colnames(df)) {
  for (m in names(markers)) {
    #pv <- p_value[unlist(markers[[m]]),ro]
    #co <- obj@cluster_data[unlist(markers[[m]]),ro]
    ex <- expression.cl[rownames(expression.cl) %in% markers[[m]],ro]
    ex[ex < 0 & !is.na(ex)] <- 0
    
    ex <- 1-exp(- teta * ex)
    n.markers <- sum(rownames(obj@raw) %in% markers[[m]])
    df[m,ro] <- sum(ex,na.rm = T)/n.markers
  }
}
write.csv(df,file = paste(out_dir,cond,"_markers_clusters_expression.csv", sep = ""))

############################
ISG.genes <- read.csv("../4.RNA_Editing_AI/ISG_fixed.txt",header = T)
colnames(ISG.genes) <- "genes"

expression.cl <- read.csv(paste0(out_dir,cond,"_expression_clusters.csv"),row.names = 1)

genes <- list("genes" = ISG.genes$genes)

df2 <- as.data.frame(matrix(nrow = 1,ncol = ncol(obj@cluster_data)))
rownames(df2) <- names(genes)
colnames(df2) <- colnames(obj@cluster_data)
teta <- -1/0.1 * (log(0.25))
#not_ass_clusters <- NA
for (ro in colnames(df2)) {
  for (m in names(genes)) {
    #pv <- p_value[unlist(genes[[m]]),ro]
    #co <- obj@cluster_data[unlist(genes[[m]]),ro]
    ex <- expression.cl[rownames(expression.cl) %in% genes[[m]],ro]
    ex[ex < 0 & !is.na(ex)] <- 0
    
    ex <- 1-exp(- teta * ex)
    n.genes <- sum(rownames(obj@raw) %in% genes[[m]])
    df2[m,ro] <- sum(ex,na.rm = T)/n.genes
  }
}
write.csv(df2,file = paste(out_dir,cond,"_ISG_clusters_expression.csv", sep = ""))


