library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)
library(harmony)
#############A
n <- readRDS("H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_Immune.rds")
new.cluster.ids <- c("Plasma",
                     "Plasma",
                     "Plasma",   #"Plasmablast"
                     "CD8+ T",
                     "Plasma",#"Plasmablast",
                     "B",
                     "B",
                     "CD8+ T",
                     "Plasma",#"Plasmablast",
                     "B",
                     "B",
                     "B",
                     "CD4+ T",
                     "B",
                     "Plasma",
                     "Plasma",
                     "unknown",
                     "Macrophage",
                     "CD8+ T",
                     "Macrophage",
                     "Macrophage",
                     "Plasma",
                     "DC",
                     "Plasma",
                     "Macrophage",
                     "Plasma",#"Plasmablast",
                     "DC",
                     "B",
                     "DC",
                     "Mast")
names(new.cluster.ids) <- levels(n)
n <- RenameIdents(n, new.cluster.ids)

DimPlot(n, reduction = "umap",label = TRUE)
n$cell_type <- Idents(n)
n <- subset(n, subset = cell_type != "unknown")
levels <- c("CD4+ T", 
            "CD8+ T",
            "B",
            "Plasma",
            "Macrophage",
            "DC",
            "Mast")
Idents(n) <- factor(Idents(n), levels = levels)
cols<-c("#F8766D","#D59100","#99A800","#00BC56","#06A4FF","#C77CFF","#FE6D8C")
DimPlot(n, reduction = "umap",label = TRUE,cols = cols)
pdf(file = "H:/WRY/data/newRNA/SI/fig4/umap_annotion_Immune.pdf",width = 7,height = 5)
DimPlot(n, reduction = "umap",label = TRUE, label.size = 4, pt.size = 0.15,cols = cols)
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/fig4/umap_annotion_ori_Immune.pdf",width = 21,height = 5)
DimPlot(n, reduction = "umap",label = TRUE, label.size = 4, pt.size = 0.15,split.by = "orig.ident",cols = cols)
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/fig4/umap_annotion_ori_Nolegend_Immune.pdf",width = 21,height = 5)
DimPlot(n, reduction = "umap",label = F, pt.size = 0.15,split.by = "orig.ident",cols = cols)
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/fig4/umap_annotion_Nolegend_Immune.pdf",width = 5,height = 5)
DimPlot(n, reduction = "umap",label = FALSE, label.size = 4, pt.size = 0.15,cols = cols)+NoLegend()
dev.off()
saveRDS(n, file = "H:/WRY/data/newRNA/SI/fig4/combind_clean_harmony_cluster_Immune.rds")
################avg
n.combined <- readRDS("H:/WRY/data/newRNA/SI/fig4/combind_clean_harmony_cluster_Immune.rds")
n <- readRDS("H:/WRY/data/newRNA/SI/new/merge_new.rds")
n<-n[,rownames(n.combined@meta.data)]
n$cell_type<-n.combined@meta.data[rownames(n.combined@meta.data),"cell_type"]
Idents(n) <- "orig.ident"
Idents(n.combined) <- "orig.ident"
DefaultAssay(n.combined)<-"originalexp"
avgn<-function(x){
  datan<-subset(n,subset = cell_type == x)
  gene_data <- data.frame(t(as.matrix(datan@assays$RNA@data)),cluster=datan$orig.ident,check.names = F)
  average_data <- aggregate(.~cluster, gene_data, mean)
  cluster_name <- average_data[,1]
  average_data <- apply(average_data[,2:ncol(average_data)],2,as.numeric)
  rownames(average_data) <- cluster_name
  average_data <- t(average_data)
  average_data <-data.frame(average_data)
  names(average_data)<-paste(x,names(average_data),sep = "")
  avg.new<-cbind(avg.new,average_data)
  return(avg.new)
}
avgt<-function(x){
  datat<-subset(n.combined,subset = cell_type == x)
  gene_data <- data.frame(t(as.matrix(datat@assays$originalexp@data)),cluster=datat$orig.ident,check.names = F)
  average_data <- aggregate(.~cluster, gene_data, mean)
  cluster_name <- average_data[,1]
  average_data <- apply(average_data[,2:ncol(average_data)],2,as.numeric)
  rownames(average_data) <- cluster_name
  average_data <- t(average_data)
  average_data <-data.frame(average_data)
  names(average_data)<-paste(x,names(average_data),sep = "")
  avg.total<-cbind(avg.total,average_data)
  return(avg.total)
}
avg.new<-data.frame(rownames(n))
avg.total<-data.frame(rownames(n.combined))
names(avg.new)<-"gene"
names(avg.total)<-"gene"
rownames(avg.new)<-avg.new$gene
rownames(avg.total)<-avg.total$gene
type <- c("CD4+ T", 
          "CD8+ T",
          "B",
          "Plasma",
          "Macrophage",
          "DC",
          "Mast")
for (i in 1:length(type)) {
  avg.new<-avgn(type[i])
}
for (i in 1:length(type)) {
  avg.total<-avgt(type[i])
}
write.csv(avg.new,file="H:/WRY/data/newRNA/SI/fig4/avgnew_Immune.csv")
write.csv(avg.total,file="H:/WRY/data/newRNA/SI/fig4/avgtotal_Immune.csv")
########
n <- readRDS("H:/WRY/data/newRNA/SI/fig4/combind_clean_harmony_cluster_Immune.rds")
ratio<-read.csv(file = "H:/WRY/data/newRNA/SI/ntr/new_ratio_qc.csv",row.names = 1)
n$new_count<-ratio[rownames(n@meta.data),"new_count"]
n$total_count<-ratio[rownames(n@meta.data),"total_count"]
n$new_ratio<-ratio[rownames(n@meta.data),"new_ratio"]
VlnPlot(n,features = c("new_count","total_count","new_ratio"),group.by="cell_type")

color<-c("#313695","#4575B4","#74ADD1","#ABD9E9","#FDAE61","#F46D43","#D73027","#A50026")
dffs <- cbind(n@reductions$umap@cell.embeddings,n@meta.data[,c("orig.ident","new_ratio")])
dffs_0h<-dffs[dffs$orig.ident == "s00h",]
dffs_0h1<-dffs_0h[order(dffs_0h[,"new_ratio"], decreasing = FALSE),]
dffs_0h1[dffs_0h1$new_ratio > 0.4,"new_ratio"] <- 0.4
pdf("H:/WRY/data/newRNA/SI/fig4/new_ratio_Immune_legend.pdf", width = 5, height = 5)
ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.3)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.4))
dev.off()
pdf("H:/WRY/data/newRNA/SI/fig4/new_ratio_Immune_0h.pdf", width = 5, height = 5)
ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.3)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.4))+NoLegend()
dev.off()
dffs_0h<-dffs[dffs$orig.ident == "s02h",]
dffs_0h1<-dffs_0h[order(dffs_0h[,"new_ratio"], decreasing = FALSE),]
dffs_0h1[dffs_0h1$new_ratio > 0.4,"new_ratio"] <- 0.4
pdf("H:/WRY/data/newRNA/SI/fig4/new_ratio_Immune_2h.pdf", width = 5, height = 5)
ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.3)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.4))+NoLegend()
dev.off()
dffs_0h<-dffs[dffs$orig.ident == "s06h",]
dffs_0h1<-dffs_0h[order(dffs_0h[,"new_ratio"], decreasing = FALSE),]
dffs_0h1[dffs_0h1$new_ratio > 0.4,"new_ratio"] <- 0.4
pdf("H:/WRY/data/newRNA/SI/fig4/new_ratio_Immune_6h.pdf", width = 5, height = 5)
ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.3)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.4))+NoLegend()
dev.off()
dffs_0h<-dffs[dffs$orig.ident == "s72h",]
dffs_0h1<-dffs_0h[order(dffs_0h[,"new_ratio"], decreasing = FALSE),]
dffs_0h1[dffs_0h1$new_ratio > 0.4,"new_ratio"] <- 0.4
pdf("H:/WRY/data/newRNA/SI/fig4/new_ratio_Immune_72h.pdf", width = 5, height = 5)
ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.3)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.4))+NoLegend()
dev.off()

################B
n <- readRDS("H:/WRY/data/newRNA/SI/fig4/combind_clean_harmony_cluster_Immune.rds")
new <- readRDS("H:/WRY/data/newRNA/SI/new/merge_new.rds")
new <- new[,rownames(n@meta.data)]
new$cell_type<-n@meta.data[rownames(new@meta.data),"cell_type"]
Idents(new)<-new$cell_type
new$sample<-paste(new$cell_type,new$orig.ident,sep="_")
cat<-c()
stage<-c("s00h","s02h","s06h","s72h")
type <- c("Macrophage",
          "DC",
          "CD4+ T", 
          "CD8+ T",
          "B",
          "Plasma")
for(i in 1:length(type)){
  for (j in 1:length(stage)) {
    cat<-c(cat,paste(type[i],stage[j],sep="_"))
  }
}
new<-subset(new, subset = sample %in% cat)
Idents(new)<-new$sample
Idents(new) <- factor(Idents(new), levels = cat)
gene<-read.csv(file = "H:/WRY/data/newRNA/SI/fig4/SI_Immune_gene.csv",header = F)
pdf(file = "H:/WRY/data/newRNA/SI/fig4/SI_Immune_dotplot.pdf",width =8,height = 6.5)
DotPlot(new,features = rev(gene[,1]),cols = c("lightgrey","#e20000")) + RotatedAxis() + coord_flip()
dev.off()

##############################D GSEA
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(tidyverse)
n <- readRDS("H:/WRY/data/newRNA/SI/fig4/combind_clean_harmony_cluster_Immune.rds")
new <- readRDS("H:/WRY/data/newRNA/SI/new/merge_new.rds")
new <- new[,rownames(n@meta.data)]
new$cell_type<-n@meta.data[rownames(new@meta.data),"cell_type"]
Idents(new)<-new$cell_type

m_df<-msigdbr(species = "Mus musculus")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
m_df1<-msigdbr(species = "Mus musculus",category = "C5")
fgsea_sets1<- m_df1 %>% split(x = .$gene_symbol, f = .$gs_name)
library(presto)
new_use<-subset(new, subset = orig.ident == "s02h" | orig.ident == "s00h")
cd8<-subset(new_use, subset = cell_type == "CD8+ T")
marker_cd8<-wilcoxauc(cd8,group_by = "orig.ident")
marker_cd81<-marker_cd8 %>% dplyr::filter(group == "s02h") %>% arrange(desc(auc)) %>% dplyr::select(feature,auc)
ranks_cd8<-deframe(marker_cd81)
fgsea_cd8<- fgsea(fgsea_sets1,stats = ranks_cd8,eps=0)
write.csv(fgsea_cd8[,c(1:7)],file="H:/WRY/data/newRNA/SI/fig4/GSEA_M5_CD8T_02h.csv")

fgsea_cd8 <- fgsea_cd8 %>% as_tibble() %>% arrange(desc(NES))

p1<-plotEnrichment(fgsea_sets1[["GOBP_IMMUNE_RESPONSE"]],
                   ranks_cd8) + labs(title="GOBP_IMMUNE_RESPONSE")
p2<-plotEnrichment(fgsea_sets1[["GOBP_T_CELL_ACTIVATION"]],
                   ranks_cd8) + labs(title="GOBP_T_CELL_ACTIVATION")
p3<-plotEnrichment(fgsea_sets1[["GOBP_REGULATION_OF_T_CELL_MEDIATED_CYTOTOXICITY"]],
                   ranks_cd8) + labs(title="GOBP_REGULATION_OF_T_CELL_MEDIATED_CYTOTOXICITY")
p4<-plotEnrichment(fgsea_sets1[["GOBP_CELL_CYCLE"]],
                   ranks_cd8) + labs(title="GOBP_CELL_CYCLE")
pdf(file="H:/WRY/data/newRNA/SI/fig4/GSEA_M5_CD8T_02h_plot.pdf",width = 10,height=6)
plot_grid(p1,p2,p3,p4,ncol = 2)
dev.off()

########################E
plot_density<-function(cell,stage,cell1){
  n <- readRDS("H:/WRY/data/newRNA/SI/fig4/combind_clean_harmony_cluster_Immune.rds")
  info<-n@meta.data
  rm(n)
  n <- readRDS("H:/WRY/data/newRNA/SI/new/merge_new.rds")
  n$cell_type<-info[rownames(n@meta.data),"cell_type"]
  Idents(n)<-n$cell_type
  marker<-FindMarkers(n,subset.ident = cell,ident.1 = stage,ident.2 = "s00h",group.by = "orig.ident",
                      min.pct = 0,logfc.threshold = 0)
  new<-read.csv(file="H:/WRY/data/newRNA/SI/fig4/avgnew_Immune.csv",row.names = 1)
  data.dir<-paste0("H:/WRY/data/newRNA/SI/ntr/Immune/ntr_gene_",stage)
  ntr<-read.delim(file=paste0(data.dir,".tsv"),row.names = 1)
  marker<-marker[rownames(ntr),]
  marker$ntr<-ntr[rownames(marker),cell1]
  marker$exp<-new[rownames(marker),paste0(cell1,stage)]
  marker$exp1<-new[rownames(marker),paste0(cell1,"s00h")]
  marker$exp2<-marker$exp1+marker$exp
  marker<-marker[marker$exp2>0.05,]
  gene<-rownames(marker)
  marker$color <- ifelse(abs(marker$avg_log2FC) >= 0.25,ifelse(marker$avg_log2FC > 0.25 ,'red','blue'),'grey')
  color<- c(red='red', grey = 'grey', blue ='blue')
  p2 <- ggplot(marker, aes(x=ntr, y=avg_log2FC, col=color)) + scale_color_manual(values=color)+
    geom_point(size=0.5) +geom_hline(yintercept = 0, lty=4,col='grey',lwd=0.6)+
    geom_hline(yintercept = -0.25, lty=4,col='grey',lwd=0.6)+
    geom_hline(yintercept = 0.25, lty=4,col='grey',lwd=0.6)+
    theme(panel.grid=element_blank(), panel.background=element_rect(color='black', fill='transparent')) +
    labs(x="New Ratio", y="log2FC",title = "new")
  write.csv(marker,file =paste0(paste0("H:/WRY/data/newRNA/SI/fig4/new_",paste(cell,stage,sep="_")),".csv"))
  
  n <- readRDS("H:/WRY/data/newRNA/SI/fig4/combind_clean_harmony_cluster_Immune.rds")
  marker<-FindMarkers(n,subset.ident = cell,ident.1 = stage,ident.2 = "s00h",group.by = "orig.ident",
                      min.pct = 0,logfc.threshold = 0)
  new<-read.csv(file="H:/WRY/data/newRNA/SI/fig4/avgtotal_Immune.csv",row.names = 1)
  data.dir<-paste0("H:/WRY/data/newRNA/SI/ntr/Immune/ntr_gene_",stage)
  ntr<-read.delim(file=paste0(data.dir,".tsv"),row.names = 1)
  marker<-marker[rownames(ntr),]
  marker$ntr<-ntr[rownames(marker),cell1]
  marker<-marker[gene,]
  marker$color <- ifelse(abs(marker$avg_log2FC) >= 0.25,ifelse(marker$avg_log2FC > 0.25 ,'red','blue'),'grey')
  color<- c(red='red', grey = 'grey', blue ='blue')
  p1 <- ggplot(marker, aes(x=ntr, y=avg_log2FC, col=color)) + scale_color_manual(values=color)+
    geom_point(size=0.5) +geom_hline(yintercept = 0, lty=4,col='grey',lwd=0.6)+
    geom_hline(yintercept = -0.25, lty=4,col='grey',lwd=0.6)+
    geom_hline(yintercept = 0.25, lty=4,col='grey',lwd=0.6)+
    theme(panel.grid=element_blank(), panel.background=element_rect(color='black', fill='transparent')) +
    labs(x="Exp", y="log2FC",title = "total")
  write.csv(marker,file =paste0(paste0("H:/WRY/data/newRNA/SI/fig4/total_",paste(cell,stage,sep="_")),".csv"))
  
  p<-p2+p1
  return(p)
}
pdf(file="H:/WRY/data/newRNA/SI/fig4/volcano_CD8T_02h.pdf",height = 5,width = 10)
p<-plot_density("CD8+ T","s02h","CD8..T")
p
dev.off()

##########G

