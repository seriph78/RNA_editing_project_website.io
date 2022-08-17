#script to take the first seurat clustering. For each cluster this check the cluster GDI.
# If there are only few genes with GDI higher then 1.5 the cluster in saved. If not, the cells are re-clustered
# recursively.
library(Rgb)
library(viridis)
library(ggplot2)

#library(COTAN)
library(stringr)
devtools::load_all("../../COTAN/")
library(Seurat)
# load package w/o installing
#load_all('../MN2_ScRNAseq/src/scRNA_analysis/COTAN/')

#library(latex2exp)
library(patchwork)

#cond = "old"

cond = "IBD_colitis"

dir.create(file.path("data/IBD/clustering/"))
out_dir <- "data/IBD/clustering/"
dataset_name <- "Seurat_colitis_data.RDS"


cores <- 15
###################
cond = "IBD_control"

dir.create(file.path("data/IBD/clustering/"))
out_dir <- "data/IBD/clustering/"

cores <- 15
in_dir <- "data/IBD/"
dataset_name <- "Seurat_control_data.RDS"

obj <- cell_homogeneous_clustering(cond = cond,out_dir = out_dir, in_dir = in_dir,cores = cores, datates_name = dataset_name)
###################
out_dir_root <- paste0(out_dir,cond,"/")

srat <- readRDS(paste(out_dir_root,"Seurat_obj_",cond,"_with_cotan_clusters.RDS",sep = ""))


gc()


markers <- list("B cells" = c("MS4A1", "CD79A", "VPREB3", "CD37", "BANK1", "LY9" ),
                "Plasma" = c( "JCHAIN", "IGHA2", "IGHA1", "IGKC"),
                "T & NK cells"= c("CCL5", "CD7", "CD3E", "CD3D", "NKG7", "CST7", "FYB", "CD8A", "HCST"),
                "LTi"=c( "IL22", "CSF2", "IL4I1", "XCL1", "KIT", "BHLHE40", "TNFRSF18", "KRT86", "NCOA7", "IL23R", "AQP3", "ALDOC", "TNFSF13B"), 
                "Myeloid"= c("LYZ", "CST3", "IL1B", "TPSAB1", "CXCL8", "CPA3", "TPSB2", "AIF1", "FCER1G", "C1QC", "C1QA", "C1QB", "S100A8"), 
                "Mitotic"= c( "STMN1", "HMGB2", "HMGN2", "HIST1H4C", "H2AFZ", "TUBB", "TUBA1B", "HMGB1", "UBE2C", "RRM2", "PTTG1", "DUT", "TYMS", "H2AFV", "MKI67", "TOP2A", "CENPF", "NUSAP1"), 
                "Fibro & Endo"= c( "CXCL14", "ADAMDEC1", "COL3A1", "APOE", "CFD", "IGFBP7", "CCL2", "COL1A2", "COL1A1", "A2M", "COL6A2", "IFITM3", "LUM", "DCN", "MFAP4"),
                "Epithelial" = c( "PHGR1", "TFF3", "FABP1", "ELF3", "ZG16", "KRT8", "SPINK4", "LGALS4", "KRT18", "CLDN4", "PIGR", "S100A6", "FCGBP", "KRT19", "CLDN3")
                
)



pdf(paste(out_dir_root,cond,"expression_seurat_cotan_cl.pdf",sep = ""),width = 14,height = 10)
plot(DimPlot(srat, reduction = "umap",group.by = "seurat_clusters",  label = TRUE, pt.size = 0.5) + NoLegend())
plot(DimPlot(srat, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "cotan") + NoLegend())
plot(DoHeatmap(srat,group.by = "seurat_clusters",size=2, hjust=0, angle = 0, features = unlist(markers),slot = 'scale.data')+ scale_fill_viridis() + theme(text = element_text(size = 10)) + NoLegend())
plot(DoHeatmap(srat,group.by = "cotan",size=2, hjust=0, angle = 0, features = unlist(markers),slot = 'scale.data')+ scale_fill_viridis() + theme(text = element_text(size = 10)) + NoLegend())

for (name.list in names(markers)) {
  p <- FeaturePlot(srat, pt.size = 0.5,features = markers[[name.list]])
  print(p + plot_annotation(title = paste0("Feature plot for ",name.list," markers"),
                            theme = theme(plot.title = element_text(size = 18,face = "bold",hjust = 0.5))))
  
}


dev.off()
## After the clustering check the analysis need to be done with:

obj <- drop.genes.cells(obj, cells = names(obj@clusters[is.na(obj@clusters)]))

ttm = clean(obj)

obj = ttm$object
ttm$pca.cell.2

gc()
#da fare!
obj = cotan_analysis(obj,cores = 15)
# saving the structure
saveRDS(obj,file = paste(out_dir_root,cond,".cotan.RDS", sep = ""))

# COEX evaluation and storing
obj = get.coex(obj)
# saving the structure
saveRDS(obj,file = paste(out_dir_root,cond,".cotan.RDS", sep = ""))

### DEA
obj = readRDS(file = paste(out_dir_root,cond,".cotan.RDS", sep = ""))


#obj_list = cotan_on_cluster_DE(obj,p_v = 0.1)
obj_list = DEA_on_clusters(obj)

obj = obj_list[[1]]
p_value = obj_list[[2]]


write.csv(p_value,file = paste(out_dir_root,cond,"p_values_clusters.csv", sep = ""))
write.csv(obj@cluster_data,file = paste(out_dir_root,cond,"coex_clusters.csv", sep = ""))

write.csv(p_value[unlist(markers),],file = paste(out_dir_root,cond,"p_values_clusters_markers.csv", sep = ""))
write.csv(obj@cluster_data[unlist(markers),],file = paste(out_dir_root,cond,"coex_clusters_markers.csv", sep = ""))

saveRDS(obj,file = paste(out_dir_root,cond,".cotan.RDS", sep = ""))
#}

obj.old <- readRDS(file = paste(out_dir_root,cond,".cotan.RDS", sep = ""))
p_value <- read.csv(file = paste(out_dir_root,cond,"p_values_clusters.csv", sep = ""), row.names = 1)

cl1_not_mergiable <- c(49,19,36,3,34,5,12,44,1,28,29)
cl2_not_mergiable <- c(50,54,37,7,35,22,2,45,32,45,39,8)


library(dendextend)

cl1_not_mergiable <- c()
cl2_not_mergiable <- c()
cl1_not_mergiable_old <- 1
round <- 0
while (!identical(cl1_not_mergiable,cl1_not_mergiable_old)) {
  round <- round + 1
  print(paste0("Start merging smallest clusters: round ", round))
  cl1_not_mergiable_old <- cl1_not_mergiable

  #merge small cluster based on distances
  cluster_data <- obj@cluster_data
  
  ######## This is the best: cosine dissimilarity
  Matrix <- as.matrix(t(cluster_data))
  sim <- Matrix / sqrt(rowSums(Matrix * Matrix))
  sim <- sim %*% t(sim)
  D_sim <- as.dist(1 - sim)
  tree <- hclust(D_sim,method = "ward.D2")
  #plot(tree)
  dend <- as.dendrogram(tree)
  
  
  ############### This part check if any little two pair of leaf clusters could be merged
  # based on the tree
  
  
  id <- NA
  for (i in c(1:length(get_nodes_attr(dend,"members")))) {
    if(get_nodes_attr(dend,"members")[i] == 2){
      id <- c(id,i+1,i+2)
    }
  }
  id <- id[2:length(id)]
  print(paste0("Created leafs id form marging: ",
               paste(get_nodes_attr(dend,"label",id = id ),collapse=" ")))
  
  cores <- 15
  
  dir <- paste0(out_dir_root,"leafs_merge/")
  dir.create(dir)
  
  p = 1
  while (p < length(id)) {
    p1 <- p
    p2 <- p+1
    p <- p2+1
    
    cl.1 <- str_remove(get_nodes_attr(dend,"label",id = id[p1]),pattern = "cl.")
    cl.2 <- str_remove(get_nodes_attr(dend,"label",id = id[p2]),pattern = "cl.")
    
    
    cond.merge <- paste0("merge_cl_",cl.1,"_",cl.2)
    print(cond.merge)
    if(cl.1 %in% cl1_not_mergiable){
      print(paste0("Clusters ", cl.1," ", cl.2,"already analyzed and not mergiable: skip."))
      next
    }
    
    mat <- srat@assays$RNA@counts[,colnames(srat@assays$RNA@counts) %in% 
                                    rownames(srat@meta.data[srat@meta.data$cotan %in% c(cl.1,cl.2),])]
    
    merged.obj <- automatic.COTAN.object.creation(df = as.data.frame(mat),out_dir = dir,GEO = "IBD GSE121380" 
                                                  ,sc.method="10X 5'end",cond = cond.merge,cores = 15,mt = F,mt_prefix = "^MT"  
    )
    GDI_data_wt1 = get.GDI(merged.obj)
    
    #Test if the number of genes with GDI > 1.5 is more than 1%
    if (dim(GDI_data_wt1[GDI_data_wt1$GDI >= 1.5,])[1]/dim(GDI_data_wt1)[1] > 0.01) {
      print(paste("Clusters ",cl.1, " and ",cl.2," too high GDI!",sep = " "))
      cl1_not_mergiable <- c(cl1_not_mergiable,cl.1)
      cl2_not_mergiable <- c(cl2_not_mergiable,cl.2)
      #cells_to_cluster = colnames(merged.obj@raw)
      #write.csv(cells_to_cluster, file = paste(out_dir,"to_recluster_",cond,"_cl.",cl,".csv",sep = ""))
      #to_rec = c(to_rec,cells_to_cluster)
    }else{
      print(paste("Clusters ",cl.1, " and ",cl.2," can be merged.",sep = " "))
      write.csv(colnames(mat),paste0(dir,"merged_clusters_",cl.1,"_",cl.2,"cell_ids.csv"))
      min.cl <- min(as.numeric(c(cl.1,cl.2)))
      max.cl <- max(as.numeric(c(cl.1,cl.2)))
      srat@meta.data[srat@meta.data$cotan == max.cl,]$cotan <- min.cl
    }
    
    
    GDI_data_wt1$color = "normal"
    genes.to.label = GDI_data_wt1[order(GDI_data_wt1$GDI,decreasing = T),][1:20,]
    
    #genes.to.label = rbind(genes.to.label,GDI_data_wt1[more_genes,])
    genes.to.label$color = "dif"
    #genes.to.label[more_genes,]$color = "mk"
    GDI_data_wt1[rownames(genes.to.label),]$color = "dif"
    #GDI_data_wt1[more_genes,]$color = "mk"
    
    # from here it is just a plot example
    mycolours <- c("dif" = "#3C5488B2","normal"="#F39B7FE5","hk"="#7E6148B2","mk"="#E64B35B2")
    
    f1 = ggplot( subset(GDI_data_wt1,!rownames(GDI_data_wt1) %in% unique(rownames(genes.to.label))),  aes(x=sum.raw.norm, y=GDI)) +  geom_point(alpha = 0.4, color = "#8491B4B2", size=2)
    
    si=12
    GDI_plot_wt1 = f1 + geom_point(data = subset(GDI_data_wt1,rownames(GDI_data_wt1) %in% c(rownames(genes.to.label),"Lhx1os","5330434G04Rik")),aes(x=sum.raw.norm, y=GDI, color=color),alpha = 1, size=2)+
      #geom_hline(yintercept=quantile(GDI$GDI)[4], linetype="dashed", color = "darkblue") +
      #geom_hline(yintercept=quantile(GDI$GDI)[3], linetype="dashed", color = "darkblue") +
      geom_hline(yintercept=1.5, linetype="dotted", color = "#3C5488B2", size= 0.5) +
      scale_color_manual("color", values = mycolours)  +
      scale_fill_manual("color", values = mycolours)  +
      xlab("log normalized counts")+ylab("GDI")+
      geom_label_repel(data =genes.to.label , aes(x=sum.raw.norm, y=GDI, label = rownames(genes.to.label),
                                                  fill=color),
                       label.size = NA,
                       alpha = 0.5,
                       direction = "both",
                       na.rm=TRUE,
                       seed = 1234) +
      theme(axis.text.x = element_text(size = si, angle = 0, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF" ),
            axis.text.y = element_text( size = si, angle = 0, hjust = 0, vjust = .5, face = "plain", colour ="#3C5488FF"),
            axis.title.x = element_text( size = si, angle = 0, hjust = .5, vjust = 0, face = "plain", colour ="#3C5488FF"),
            axis.title.y = element_text( size = si, angle = 90, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF"),
            legend.title = element_blank(),
            legend.text = element_text(color = "#3C5488FF",face ="italic" ),
            legend.position = "none")  +ggtitle(paste(cond.merge,obj@n_cells,sep = " "))
    
    
    pdf(paste(dir,cond.merge, ".GDI_plots.pdf", sep = ""), onefile=TRUE)
    plot(GDI_plot_wt1)
    graphics.off()
    
    rm(merged.obj)
    gc()
  }
  
  saveRDS(srat,paste(out_dir_root,"Seurat_obj_",cond,"_with_cotan_clusters_merged.RDS",sep = ""))
  #srat <- readRDS(paste(out_dir_root,"Seurat_obj_",cond,"_with_cotan_clusters_merged.RDS",sep = ""))
  gc()
  
  
  # Update the homogenus cluster array in cotan object
  obj.clusters.new <- srat@meta.data$cotan
  names(obj.clusters.new) <- rownames(srat@meta.data)
  if(all(colnames(obj@raw) %in% names(obj.clusters.new))){
    obj@clusters <- obj.clusters.new
  }else{
    print("Problem! not all all(colnames(obj@raw) %in% names(obj.clusters.new)")
    break
  }
  
  
  # New DEA on clusters
  clusters.names = unique(obj@clusters)[!is.na(unique(obj@clusters))]
  list.clusters = list(names(obj@clusters[obj@clusters %in% clusters.names[1]]))
  names(list.clusters)=clusters.names[1]
  for (c in c(2:length(clusters.names))) {
    tmp = list(names(obj@clusters[obj@clusters %in% clusters.names[c]]))
    names(tmp)= clusters.names[c]
    list.clusters = c(list.clusters,tmp)
  }
  
  rm(srat)
  gc()
  
  
  
  obj_list = DEA_on_clusters(obj,list.clusters)
  gc()
  
  srat <- readRDS(paste(out_dir_root,"Seurat_obj_",cond,"_with_cotan_clusters_merged.RDS",sep = ""))
  obj = obj_list[[1]]
  
  p_value = obj_list[[2]]
  
  rm(obj_list)
  gc()
  
  write.csv(p_value,file = paste(out_dir_root,"p_values_clusters_merged.csv", sep = ""))
  write.csv(obj@cluster_data,file = paste(out_dir_root,"coex_clusters_merged.csv", sep = ""))
  
  write.csv(p_value[unlist(markers),],file = paste(out_dir_root,"p_values_clusters_merged_markers.csv", sep = ""))
  write.csv(obj@cluster_data[unlist(markers),],file = paste(out_dir_root,"coex_clusters_merged_markers.csv", sep = ""))
  
  saveRDS(obj,file = paste(out_dir_root,cond,"_merged_cotan.RDS", sep = ""))
  

}

p_value <- read.csv(file = paste(out_dir_root,"p_values_clusters_merged.csv", sep = ""),row.names = 1)


# cluster assignment percentage
expression.cl <- get.clusters.change.expression(obj)
df <- as.data.frame(matrix(nrow = length(names(markers)),ncol = ncol(obj@cluster_data)))
rownames(df) <- names(markers)
colnames(df) <- colnames(obj@cluster_data)
#not_ass_clusters <- NA
for (ro in colnames(df)) {
  for (m in names(markers)) {
    pv <- p_value[unlist(markers[[m]]),ro]
    co <- obj@cluster_data[unlist(markers[[m]]),ro]
    ex <- expression.cl[unlist(markers[[m]]),ro]
    df[m,ro] <- (sum(pv < 0.05 & co > 0 & ex >= 0.1, na.rm = T)/length(unlist(markers[[m]])))*100
  }
}
write.csv(df,file = paste(out_dir_root,cond,"_markers_clusters_percentage_expression.csv", sep = ""))


# cluster assignment with expression increment
expression.cl <- get.clusters.change.expression(obj)
df <- as.data.frame(matrix(nrow = length(names(markers)),ncol = ncol(obj@cluster_data)))
rownames(df) <- names(markers)
colnames(df) <- colnames(obj@cluster_data)
teta <- -1/0.1 * (log(0.25))
#not_ass_clusters <- NA
for (ro in colnames(df)) {
  for (m in names(markers)) {
    #pv <- p_value[unlist(markers[[m]]),ro]
    #co <- obj@cluster_data[unlist(markers[[m]]),ro]
    ex <- expression.cl[unlist(markers[[m]]),ro]
    ex[ex < 0 & !is.na(ex)] <- 0
    
    ex <- 1-exp(- teta * ex)
    n.markers <- sum(unlist(markers[[m]]) %in% rownames(obj@raw))
    df[m,ro] <- sum(ex,na.rm = T)/n.markers
  }
}
write.csv(df,file = paste(out_dir_root,cond,"_markers_clusters_expression_function_exp_markers.csv", sep = ""))


###########################################
#merge small cluster based on distances
cluster_data <- obj@cluster_data

######## This is the best: cosine dissimilarity
Matrix <- as.matrix(t(cluster_data))
sim <- Matrix / sqrt(rowSums(Matrix * Matrix))
sim <- sim %*% t(sim)
D_sim <- as.dist(1 - sim)
tree <- hclust(D_sim,method = "ward.D2")
plot(tree)
library(dendextend)
dend <- as.dendrogram(tree)
cut = cutree(tree, k = 5)
dend =branches_color(dend,k=5,col=c("#4DBBD5FF","#91D1C2FF","#E64B35FF","gray80","#3C5488FF"))#,"#F39B7FFF","gray80" ),groupLabels = T)
dend =color_labels(dend,k=5)#,labels = rownames(pca_1),col=pca_1$colors)
dend %>%
  dendextend::set("labels") %>%
  #  set("branches_k_color", value = c("gray80","#4DBBD5FF","#91D1C2FF" ,"gray80","#F39B7FFF","#E64B35FF","#3C5488FF"), k = 7) %>%
  plot(horiz=F, axes=T)

df$cluster <- str_split(rownames(df),"[.]",simplify = T)[,2]

df$hclust <- cut

seurat_obj <- readRDS("/media/asuna/seriph/scRNAseq_Nothobranchius_furzeri/results/cleaning/old_young_seurat_dataset.RDS")

head(tmp_seurat@meta.data)
df <- as.data.frame(df)
df$type1 <- NA
df[df$hclust == 2,]$type1 <- "neuron_gaba"
df[df$hclust == 1,]$type1 <- "Glia/SMC/ENDO/MICRO/OPC"
df[df$hclust == 4,]$type1 <- "Other mixed"
df[df$hclust == 3,]$type1 <- "neuron_glut"
df[df$hclust == 5,]$type1 <- "neurons"


head(df,n = 15)
tmp_seurat@meta.data$row.names <- rownames(tmp_seurat@meta.data)
meta.data <- tmp_seurat@meta.data
dim(meta.data)
meta.data <- merge(meta.data,df[,c("cluster","type1")],by.x="cotan",
                   by.y ="cluster",all.x = T)
rownames(meta.data) <- meta.data$row.names
tmp_seurat@meta.data <- meta.data
tmp_seurat@meta.data = tmp_seurat@meta.data[names(tmp_seurat@active.ident),]

plot(DimPlot(tmp_seurat, reduction = "umap",group.by = "type1",  label = TRUE, pt.size = 0.5) + NoLegend())
plot(DimPlot(tmp_seurat, reduction = "umap",group.by = "type1",  label = TRUE, pt.size = 0.5) + NoLegend())
#########################################
tmp_seurat@meta.data$Type2 = NA
tmp_seurat@meta.data[tmp_seurat@meta.data$cotan %in% c(12,6),]$Type2 = "A"
tmp_seurat@meta.data[tmp_seurat@meta.data$cotan %in% c(20,16,15,32,42,14),]$Type2 = "B"
tmp_seurat@meta.data[tmp_seurat@meta.data$cotan %in% c(45,56,58),]$Type2 = "C"
tmp_seurat@meta.data[tmp_seurat@meta.data$cotan %in% c(53,60,57,62,61,51,59),]$Type2 = "D"
tmp_seurat@meta.data[tmp_seurat@meta.data$cotan %in% c(13,5,35,44,41),]$Type2 = "E"
tmp_seurat@meta.data[tmp_seurat@meta.data$cotan %in% c(30,49,54,55,47,52,50),]$Type2 = "N1"
tmp_seurat@meta.data[tmp_seurat@meta.data$cotan %in% c(33,34,2,19,29),]$Type2 = "N2"
tmp_seurat@meta.data[tmp_seurat@meta.data$type1 == "neuron_gaba",]$Type2 = "N3"
tmp_seurat@meta.data[tmp_seurat@meta.data$type1 == "neuron_glut",]$Type2 = "N4"

plot(DimPlot(tmp_seurat, reduction = "umap",group.by = "Type2",  label = TRUE, pt.size = 0.5) + NoLegend())

tmp_seurat <- RunTSNE(tmp_seurat, dims = 1:25)
plot(DimPlot(tmp_seurat, reduction = "tsne",group.by = "Type2",  label = TRUE, pt.size = 0.5) + NoLegend())

tmp_seurat@meta.data$Type3 = NA
tmp_seurat@meta.data[tmp_seurat@meta.data$cotan %in% c(12,6),]$Type3 = "A"
tmp_seurat@meta.data[tmp_seurat@meta.data$cotan %in% c(20,16,15),]$Type3 = "B1"
tmp_seurat@meta.data[tmp_seurat@meta.data$cotan %in% c(32,42,14),]$Type3 = "B2"
tmp_seurat@meta.data[tmp_seurat@meta.data$cotan %in% c(45,56,58),]$Type3 = "C"
tmp_seurat@meta.data[tmp_seurat@meta.data$cotan %in% c(53),]$Type3 = "D1"
tmp_seurat@meta.data[tmp_seurat@meta.data$cotan %in% c(60,57),]$Type3 = "D2"
tmp_seurat@meta.data[tmp_seurat@meta.data$cotan %in% c(62),]$Type3 = "D3"
tmp_seurat@meta.data[tmp_seurat@meta.data$cotan %in% c(61,51,59),]$Type3 = "D4"
tmp_seurat@meta.data[tmp_seurat@meta.data$cotan %in% c(13,5,35),]$Type3 = "E1"
tmp_seurat@meta.data[tmp_seurat@meta.data$cotan %in% c(44,41),]$Type3 = "E2"
tmp_seurat@meta.data[tmp_seurat@meta.data$cotan %in% c(30,49,54,55,47,52,50),]$Type3 = "N1"
tmp_seurat@meta.data[tmp_seurat@meta.data$cotan %in% c(33,34,2,19,29),]$Type3 = "N2"
tmp_seurat@meta.data[tmp_seurat@meta.data$type1 == "neuron_gaba",]$Type3 = "N3"
tmp_seurat@meta.data[tmp_seurat@meta.data$type1 == "neuron_glut",]$Type3 = "N4"

#####################
# Colors
library(RColorBrewer)
n <- 11
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
####################################

plot(DimPlot(tmp_seurat,cols = sample(col_vector,14), reduction = "umap",group.by ="Type3",  label = TRUE, pt.size = 0.5))

plot(DimPlot(tmp_seurat,cols = sample(col_vector,14), reduction = "tsne",group.by = "Type3",  label = TRUE, pt.size = 0.5) + NoLegend())

Idents(tmp_seurat) <- "Type3"
seurat.markers <- FindAllMarkers(tmp_seurat,only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)
write.csv(seurat.markers[seurat.markers$gene %in% unlist(markers),],"results/cell_assignment/seurat_markers_type3_cl.csv")


StackedVlnPlot(tmp_seurat,features = unlist(markers),
               ncol = 1,
               pt = 0.1,
               split.by = "orig.ident",
               group.by = "Type3",split.plot = T,combine = T)



library(dplyr)
seurat.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

saveRDS(tmp_seurat,"results/clustering/Seurat_dataset.RDS")



### DEA after assignmet
clusters.names = unique(tmp_seurat@meta.data$Type3)
#clusters.names = unique(obj@clusters)[!is.na(unique(obj@clusters))]
#list.clusters = list(names(obj@clusters[obj@clusters %in% clusters.names[1]]))
list.clusters = NA# array(data = NA,dim = length(clusters.names))
for (c in c(1:length(clusters.names))) {
  tmp = list(rownames(tmp_seurat@meta.data[tmp_seurat@meta.data$Type3 == clusters.names[c],]))
  names(tmp) = clusters.names[c]
  if (length(list.clusters[[1]])==1) {
    list.clusters = tmp
  }else{
    list.clusters = c(list.clusters,tmp)
  }
}

#obj_list = cotan_on_cluster_DE(obj,p_v = 0.1)
obj_list = cotan_on_cluster_DE(obj,list.clusters)
obj = obj_list[[1]]
p_value = obj_list[[2]]
p_value <- p_value[,2:ncol(p_value)]

write.csv(p_value,file = paste(out_dir_root,cond,"p_values_Type3_clusters.csv", sep = ""))
write.csv(obj@cluster_data,file = paste(out_dir_root,cond,"coex_Type3_clusters.csv", sep = ""))

write.csv(p_value[unlist(markers),],file = paste(out_dir_root,cond,"p_values_Type3_clusters_markers.csv", sep = ""))
write.csv(obj@cluster_data[unlist(markers),],file = paste(out_dir_root,cond,"coex_Type3_clusters_markers.csv", sep = ""))

saveRDS(obj,file = paste(out_dir_root,cond,".Type3_clusters.cotan.RDS", sep = ""))

cluster_data_type3 <- obj@cluster_data
dist.cl_type3 <- dist(t(cluster_data_type3),)
plot( hclust(dist.cl_type3))

######## This is the best: coesine dissimilarity
Matrix <- as.matrix(t(cluster_data_type3))
sim <- Matrix / sqrt(rowSums(Matrix * Matrix))
sim <- sim %*% t(sim)
D_sim <- as.dist(1 - sim)
tree_type3 <- hclust(D_sim,method = "ward.D2")
plot(tree_type3)
library(dendextend)
dend_type3 <- as.dendrogram(tree_type3)
cut = cutree(tree_type3, k = 5)
dend_type3 =branches_color(dend_type3,k=5,col=c("#4DBBD5FF","#91D1C2FF","#E64B35FF","gray80","#3C5488FF"))#,"#F39B7FFF","gray80" ),groupLabels = T)
dend_type3 =color_labels(dend_type3,k=5)#,labels = rownames(pca_1),col=pca_1$colors)
dend_type3 %>%
  dendextend::set("labels") %>%
  #  set("branches_k_color", value = c("gray80","#4DBBD5FF","#91D1C2FF" ,"gray80","#F39B7FFF","#E64B35FF","#3C5488FF"), k = 7) %>%
  plot(horiz=F, axes=T)




####### Cell types
tmp_seurat@meta.data$cell_type = NA
tmp_seurat@meta.data[tmp_seurat@meta.data$Type3 %in% c("E1"),]$cell_type = "Olig."
tmp_seurat@meta.data[tmp_seurat@meta.data$Type3 %in% c("N3" ),]$cell_type = "Gaba1"
tmp_seurat@meta.data[tmp_seurat@meta.data$Type3 %in% c("N4" ),]$cell_type = "Glut1"
tmp_seurat@meta.data[tmp_seurat@meta.data$Type3 %in% c("D2" ),]$cell_type = "Endot."
tmp_seurat@meta.data[tmp_seurat@meta.data$Type3 %in% c("B1" ),]$cell_type = "Red Blood c."
tmp_seurat@meta.data[tmp_seurat@meta.data$Type3 %in% c( "A"  ),]$cell_type = "A"
tmp_seurat@meta.data[tmp_seurat@meta.data$Type3 %in% c("N1" ),]$cell_type = "Gaba2"
tmp_seurat@meta.data[tmp_seurat@meta.data$Type3 %in% c( "B2" ),]$cell_type = "B2"
tmp_seurat@meta.data[tmp_seurat@meta.data$Type3 %in% c( "C"  ),]$cell_type = "Glia"
tmp_seurat@meta.data[tmp_seurat@meta.data$Type3 %in% c( "E2" ),]$cell_type = "OPC"
tmp_seurat@meta.data[tmp_seurat@meta.data$Type3 %in% c( "D4" ),]$cell_type = "Microglia2"
tmp_seurat@meta.data[tmp_seurat@meta.data$Type3 %in% c( "D3" ),]$cell_type = "Microglia1"
tmp_seurat@meta.data[tmp_seurat@meta.data$Type3 %in% c( "N2" ),]$cell_type = "Glut2"
tmp_seurat@meta.data[tmp_seurat@meta.data$Type3 %in% c( "D1"),]$cell_type = "SMC"

Idents(tmp_seurat) <- "cell_type"
plot(DimPlot(tmp_seurat,cols = sample(col_vector,14), reduction = "umap",  label = TRUE, pt.size = 0.5) + NoLegend())
plot(DimPlot(tmp_seurat,cols = sample(col_vector,14),group.by = "Type3", reduction = "umap",  label = TRUE, pt.size = 0.5) + NoLegend())

table(tmp_seurat@meta.data[tmp_seurat@meta.data$orig.ident == "Young.5w","cell_type"])/dim(tmp_seurat@meta.data[tmp_seurat@meta.data$orig.ident == "Young.5w",])[1]
table(tmp_seurat@meta.data[tmp_seurat@meta.data$orig.ident == "Old.38w","cell_type"])/dim(tmp_seurat@meta.data[tmp_seurat@meta.data$orig.ident == "Old.38w",])[1]

saveRDS(tmp_seurat,"results/clustering/Seurat_dataset.RDS")

pdf("results/cell_assignment/violin_markers_markers_cell_type.pdf",width = 10, height = 15)
print(StackedVlnPlot(tmp_seurat,features = unlist(markers)[1:16],
                     ncol = 1,
                     pt = 0.1,
                     split.by = "orig.ident",
                     group.by = "cell_type",split.plot = T,combine = T))
print(StackedVlnPlot(tmp_seurat,features = unlist(markers)[17:32],
                     ncol = 1,
                     pt = 0.1,
                     split.by = "orig.ident",
                     group.by = "cell_type",split.plot = T,combine = T))
dev.off()

DoHeatmap(tmp_seurat, features = unique(unlist(markers)))+
  scale_fill_gradient2(low = "#64A7C9", mid = "white", high = "#690321")+
  theme(legend.position="bottom")
