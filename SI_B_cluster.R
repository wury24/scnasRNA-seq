library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)
library(harmony)
n <- readRDS("H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_Immune.rds")
########
new.cluster.ids <- c("Plasma",
                     "Plasma",
                     "Plasmablast",
                     "CD8+ T",
                     "Plasmablast",
                     "B",
                     "B",
                     "CD8+ T",
                     "Plasmablast",
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
                     "Plasmablast",
                     "DC",
                     "B",
                     "DC",
                     "Mast")
names(new.cluster.ids) <- levels(n)
n <- RenameIdents(n, new.cluster.ids)

DimPlot(n, reduction = "umap",label = TRUE)
DimPlot(n, reduction = "umap",label = TRUE,group.by = "seurat_clusters")
n$cell_type <- Idents(n)
n<-subset(n, subset = seurat_clusters != "5")
b<-subset(n, subset = cell_type == "B" |
            cell_type == "Plasma" | cell_type == "Plasmablast")
saveRDS(b,file = "H:/WRY/data/newRNA/SI/cluster/B/B_total_data.rds")
new <- readRDS("H:/WRY/data/newRNA/SI/new/merge_new.rds")
new <- new[,rownames(b@meta.data)]
new$cell_type<-b@meta.data[rownames(new@meta.data),"cell_type"]
new <- FindVariableFeatures(new)
new <- ScaleData(new)
new <- RunPCA(new)
new <- RunHarmony(new, group.by.vars = "orig.ident",assay.use = "RNA",plot_convergence = TRUE)
new <- RunUMAP(new, reduction = "harmony", dims = 1:30)
new <- FindNeighbors(new, reduction = "harmony", dims = 1:30) 
new <- FindClusters(new,resolution=0.3)
pdf(file = "H:/WRY/data/newRNA/SI/cluster/B/new_cluster.pdf",width = 6,height = 5)
DimPlot(new, reduction = "umap",label = TRUE)
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/cluster/B/new_cluster_celltype.pdf",width = 6,height = 5)
DimPlot(new, reduction = "umap",label = TRUE,group.by = "cell_type")
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/cluster/B/new_cluster_ori.pdf",width = 21,height = 5)
DimPlot(new, reduction = "umap", split.by = "orig.ident")
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/cluster/B/marker_Igha.pdf",width = 21,height = 5)
FeaturePlot(new, features = c("Igha"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, label.size = 2,pt.size = 0.3)&theme(legend.position = "right")
dev.off()
new$Igha_total<-b@assays$originalexp@data["Igha",rownames(new@meta.data)]
new$Mki67_total<-b@assays$originalexp@data["Mki67",rownames(new@meta.data)]
ratio<-read.csv(file = "H:/WRY/data/newRNA/SI/ntr/new_ratio_qc.csv",row.names = 1)
new$new_count<-ratio[rownames(new@meta.data),"new_count"]
new$total_count<-ratio[rownames(new@meta.data),"total_count"]
new$new_ratio<-ratio[rownames(new@meta.data),"new_ratio"]
pdf(file = "H:/WRY/data/newRNA/SI/cluster/B/marker.pdf",width = 21,height = 20)
FeaturePlot(new, features = c("Igha","Igha_total","Mki67","Mki67_total"), reduction = "umap",split.by = "orig.ident",cols = c("grey", "red"), label = TRUE, label.size = 2,pt.size = 0.3)&theme(legend.position = "right")
dev.off()
#info <- cbind(new@reductions$umap@cell.embeddings,new@meta.data[,c("orig.ident","cell_type")])
#write.csv(info, file = "H:/WRY/data/newRNA/SI/cluster/B/newinfo.csv")

###############
color<-c("#313695","#4575B4","#74ADD1","#ABD9E9","#FDAE61","#F46D43","#D73027","#A50026")
ratio<-read.csv(file = "H:/WRY/data/newRNA/SI/ntr/new_ratio_qc.csv",row.names = 1)
new$new_ratio<-ratio[rownames(new@meta.data),"new_ratio"]
dffs <- cbind(new@reductions$umap@cell.embeddings,new@meta.data[,c("orig.ident","new_ratio")])
dffs_0h<-dffs[dffs$orig.ident == "s00h",]
dffs_0h1<-dffs_0h[order(dffs_0h[,"new_ratio"], decreasing = FALSE),]
dffs_0h1[dffs_0h1$new_ratio > 0.4,"new_ratio"] <- 0.4
p1<-ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.5)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.4))

dffs_0h<-dffs[dffs$orig.ident == "s02h",]
dffs_0h1<-dffs_0h[order(dffs_0h[,"new_ratio"], decreasing = FALSE),]
dffs_0h1[dffs_0h1$new_ratio > 0.4,"new_ratio"] <- 0.4
p2<-ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.5)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.4))

dffs_0h<-dffs[dffs$orig.ident == "s06h",]
dffs_0h1<-dffs_0h[order(dffs_0h[,"new_ratio"], decreasing = FALSE),]
dffs_0h1[dffs_0h1$new_ratio > 0.4,"new_ratio"] <- 0.4
p3<-ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.5)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.4))

dffs_0h<-dffs[dffs$orig.ident == "s72h",]
dffs_0h1<-dffs_0h[order(dffs_0h[,"new_ratio"], decreasing = FALSE),]
dffs_0h1[dffs_0h1$new_ratio > 0.4,"new_ratio"] <- 0.4
p4<-ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.5)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.4))

pdf("H:/WRY/data/newRNA/SI/cluster/B/new_ratio_ori_B.pdf", width = 21, height = 4.5)
plot_grid(p1,p2,p3,p4,ncol = 4)
dev.off()

new.cluster.ids <- c("Plasma",
                     "B",
                     "B",
                     "Plasma",
                     "Plasma",
                     "Plasma",
                     "B",
                     "Plasma",
                     "B")
names(new.cluster.ids) <- levels(new)
new <- RenameIdents(new, new.cluster.ids)

DimPlot(new, reduction = "umap",label = TRUE)
new$cell_type <- Idents(new)
pdf(file = "H:/WRY/data/newRNA/SI/cluster/B/new_cluster_annotation.pdf",width = 6,height = 5)
DimPlot(new, reduction = "umap",label = TRUE)
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/cluster/B/new_cluster_ori.pdf",width = 21,height = 5)
DimPlot(new, reduction = "umap", split.by = "orig.ident")
dev.off()
saveRDS(new,file = "H:/WRY/data/newRNA/SI/cluster/B/combind_clean_harmony_cluster_B_new.rds")

new<-readRDS(file = "H:/WRY/data/newRNA/SI/cluster/B/combind_clean_harmony_cluster_B_new.rds")
n0<-subset(new, subset = orig.ident == "s00h")
marker <- FindAllMarkers(n0, only.pos = TRUE,min.pct = 0.1,logfc.threshold = 0.25)
marker <- marker[marker$p_val_adj < 0.05,]

write.csv(marker,file = "H:/WRY/data/newRNA/SI/cluster/B/B_new_marker0h.csv")
top10 <- marker %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(top10,file = "H:/WRY/data/newRNA/BM/fig1/SI_E_marker_new_top10_n00h.csv")
new <- FindVariableFeatures(new)
new <- ScaleData(new,features = top10$gene)
data <- subset(new, downsample = 300)
pdf("H:/WRY/data/newRNA/BM/fig1/SI_E_Heatmap_new_top10.pdf", width = 12, height = 20)
DoHeatmap(data, features = top10$gene, size = 3, angle = 90, label = TRUE)+
  scale_fill_gradientn(colors = brewer.pal(n = 9, name = "Reds"))
dev.off()
pdf("H:/WRY/data/newRNA/BM/fig1/SI_E_Heatmap_new_top10_nolegend.pdf", width = 12, height = 8)
DoHeatmap(data, features = top10$gene, size = 3, angle = 90, label = FALSE)+NoLegend()+
  scale_fill_gradientn(colors = brewer.pal(n = 9, name = "Reds"))
dev.off()
pdf("H:/WRY/data/newRNA/BM/fig1/SI_E_Heatmap_new_top10_nolegend1.pdf", width = 12, height = 8)
DoHeatmap(data, features = top10$gene, size = 3, angle = 90, label = FALSE)+NoLegend()+
  scale_fill_gradientn(colors = rev(brewer.pal(n = 10, name = "RdBu")))
dev.off()

##############
n<-readRDS(file = "H:/WRY/data/newRNA/SI/cluster/B/combind_clean_harmony_cluster_B_new.rds")
n.combined <- readRDS("H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_Immune.rds")
n.combined<-n.combined[,rownames(n@meta.data)]
n.combined$cell_type<-n@meta.data[rownames(n.combined@meta.data),"cell_type"]
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
type <- c("B",
          "Plasma")
for (i in 1:length(type)) {
  avg.new<-avgn(type[i])
}
for (i in 1:length(type)) {
  avg.total<-avgt(type[i])
}
write.csv(avg.new,file="H:/WRY/data/newRNA/SI/cluster/B/avgnew_Neu.csv")
write.csv(avg.total,file="H:/WRY/data/newRNA/SI/cluster/B/avgtotal_Neu.csv")

################
n<-readRDS(file = "H:/WRY/data/newRNA/SI/cluster/B/combind_clean_harmony_cluster_B_new.rds")
n.combined <- readRDS("H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_Immune.rds")
n.combined<-n.combined[,rownames(n@meta.data)]
n.combined$cell_type<-n@meta.data[rownames(n.combined@meta.data),"cell_type"]
Idents(n.combined)<-n.combined$cell_type
levels <- c("B",
            "Plasma")
stage <- c("s00h","s02h","s06h","s72h")
main_dir <- "H:/WRY/data/newRNA/SI/cluster/B/stage_total"
for (i in 1:length(levels)){
  dir.create(file.path(main_dir,levels[i]))
  for(j in 2:length(stage)){
    marker<-FindMarkers(n.combined,subset.ident = levels[i],ident.1 = stage[j],ident.2 = stage[1],group.by = "orig.ident",min.pct = 0.25,logfc.threshold=0.25,only.pos = TRUE)
    write.csv(marker,file =  paste("H:/WRY/data/newRNA/SI/cluster/B/stage_total/",levels[i],"/","marker_",stage[j],"_",stage[1],"_total.csv",sep=""))
    rm(marker)
  }
}

main_dir <- "H:/WRY/data/newRNA/SI/cluster/B/stage_new"
for (i in 1:length(levels)){
  dir.create(file.path(main_dir,levels[i]))
  for(j in 2:length(stage)){
    marker<-FindMarkers(n,subset.ident = levels[i],ident.1 = stage[j],ident.2 = stage[1],group.by = "orig.ident",min.pct = 0.1,logfc.threshold=0.25,only.pos = TRUE)
    write.csv(marker,file =  paste("H:/WRY/data/newRNA/SI/cluster/B/stage_new/",levels[i],"/","marker_",stage[j],"_",stage[1],"_new.csv",sep=""))
    rm(marker)
  }
}
####################################
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(tidyverse)
n <- readRDS("H:/WRY/data/newRNA/SI/cluster/B/combind_clean_harmony_cluster_B_new.rds")
m_df<-msigdbr(species = "Mus musculus")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
m_df1<-msigdbr(species = "Mus musculus",category = "C5")
fgsea_sets1<- m_df1 %>% split(x = .$gene_symbol, f = .$gs_name)
library(presto)
new_use<-subset(n, subset = orig.ident == "s72h" | orig.ident == "s00h")
cd8<-subset(new_use, subset = cell_type == "Plasma")
marker_cd8<-wilcoxauc(cd8,group_by = "orig.ident")
marker_cd81<-marker_cd8 %>% dplyr::filter(group == "s72h") %>% arrange(desc(auc)) %>% dplyr::select(feature,auc)
ranks_cd8<-deframe(marker_cd81)
fgsea_cd8<- fgsea(fgsea_sets1,stats = ranks_cd8,eps=0)
write.csv(fgsea_cd8[,c(1:7)],file="H:/WRY/data/newRNA/SI/cluster/B/GSEA_M5_Plasma_72h.csv")

new_use<-subset(n, subset = orig.ident == "s72h" | orig.ident == "s00h")
cd8<-subset(new_use, subset = cell_type == "B")
marker_cd8<-wilcoxauc(cd8,group_by = "orig.ident")
marker_cd81<-marker_cd8 %>% dplyr::filter(group == "s72h") %>% arrange(desc(auc)) %>% dplyr::select(feature,auc)
ranks_cd8<-deframe(marker_cd81)
fgsea_cd8<- fgsea(fgsea_sets1,stats = ranks_cd8,eps=0)
write.csv(fgsea_cd8[,c(1:7)],file="H:/WRY/data/newRNA/SI/cluster/B/GSEA_M5_B_72h.csv")
####################
new_use<-subset(n, subset = orig.ident == "s06h" | orig.ident == "s00h")
cd8<-subset(new_use, subset = cell_type == "Plasma")
marker_cd8<-wilcoxauc(cd8,group_by = "orig.ident")
marker_cd81<-marker_cd8 %>% dplyr::filter(group == "s06h") %>% arrange(desc(auc)) %>% dplyr::select(feature,auc)
ranks_cd8<-deframe(marker_cd81)
fgsea_cd8<- fgsea(fgsea_sets1,stats = ranks_cd8,eps=0)
write.csv(fgsea_cd8[,c(1:7)],file="H:/WRY/data/newRNA/SI/cluster/B/GSEA_M5_Plasma_06h.csv")

new_use<-subset(n, subset = orig.ident == "s06h" | orig.ident == "s00h")
cd8<-subset(new_use, subset = cell_type == "B")
marker_cd8<-wilcoxauc(cd8,group_by = "orig.ident")
marker_cd81<-marker_cd8 %>% dplyr::filter(group == "s06h") %>% arrange(desc(auc)) %>% dplyr::select(feature,auc)
ranks_cd8<-deframe(marker_cd81)
fgsea_cd8<- fgsea(fgsea_sets1,stats = ranks_cd8,eps=0)
write.csv(fgsea_cd8[,c(1:7)],file="H:/WRY/data/newRNA/SI/cluster/B/GSEA_M5_B_06h.csv")
####################
new_use<-subset(n, subset = orig.ident == "s02h" | orig.ident == "s00h")
cd8<-subset(new_use, subset = cell_type == "Plasma")
marker_cd8<-wilcoxauc(cd8,group_by = "orig.ident")
marker_cd81<-marker_cd8 %>% dplyr::filter(group == "s02h") %>% arrange(desc(auc)) %>% dplyr::select(feature,auc)
ranks_cd8<-deframe(marker_cd81)
fgsea_cd8<- fgsea(fgsea_sets1,stats = ranks_cd8,eps=0)
write.csv(fgsea_cd8[,c(1:7)],file="H:/WRY/data/newRNA/SI/cluster/B/GSEA_M5_Plasma_02h.csv")

new_use<-subset(n, subset = orig.ident == "s02h" | orig.ident == "s00h")
cd8<-subset(new_use, subset = cell_type == "B")
marker_cd8<-wilcoxauc(cd8,group_by = "orig.ident")
marker_cd81<-marker_cd8 %>% dplyr::filter(group == "s02h") %>% arrange(desc(auc)) %>% dplyr::select(feature,auc)
ranks_cd8<-deframe(marker_cd81)
fgsea_cd8<- fgsea(fgsea_sets1,stats = ranks_cd8,eps=0)
write.csv(fgsea_cd8[,c(1:7)],file="H:/WRY/data/newRNA/SI/cluster/B/GSEA_M5_B_02h.csv")
###########################################################
############
new_use<-subset(n, subset = orig.ident == "s02h" | orig.ident == "s00h")
cd8<-subset(new_use, subset = cell_type == "B")
marker_cd8<-wilcoxauc(cd8,group_by = "orig.ident")
marker_cd81<-marker_cd8 %>% dplyr::filter(group == "s02h") %>% arrange(desc(auc)) %>% dplyr::select(feature,auc)
ranks_cd8<-deframe(marker_cd81)
p1<-plotEnrichment(fgsea_sets1[["GOCC_IMMUNOGLOBULIN_COMPLEX"]],
                   ranks_cd8) + labs(title="GOCC_IMMUNOGLOBULIN_COMPLEX(B)")
p2<-plotEnrichment(fgsea_sets1[["GOBP_B_CELL_MEDIATED_IMMUNITY"]],
                   ranks_cd8) + labs(title="GOBP_B_CELL_MEDIATED_IMMUNITY(B)")
p3<-plotEnrichment(fgsea_sets1[["GOBP_ADAPTIVE_IMMUNE_RESPONSE"]],
                   ranks_cd8) + labs(title="GOBP_ADAPTIVE_IMMUNE_RESPONSE(B)")
pdf(file="H:/WRY/data/newRNA/SI/cluster/B/GSEA_M5_B_02h_plot_1.pdf",width = 10,height=6)
plot_grid(p1,p2,p3,ncol = 2)
dev.off()
plotEnrichment(fgsea_sets1[["GOBP_IMMUNE_RESPONSE"]],
                ranks_cd8) + labs(title="GOBP_IMMUNE_RESPONSE(B)")

new_use<-subset(n, subset = orig.ident == "s06h" | orig.ident == "s00h")
cd8<-subset(new_use, subset = cell_type == "B")
marker_cd8<-wilcoxauc(cd8,group_by = "orig.ident")
marker_cd81<-marker_cd8 %>% dplyr::filter(group == "s06h") %>% arrange(desc(auc)) %>% dplyr::select(feature,auc)
ranks_cd8<-deframe(marker_cd81)
p1<-plotEnrichment(fgsea_sets1[["GOCC_IMMUNOGLOBULIN_COMPLEX"]],
                   ranks_cd8) + labs(title="GOCC_IMMUNOGLOBULIN_COMPLEX(B)")
p2<-plotEnrichment(fgsea_sets1[["GOBP_B_CELL_MEDIATED_IMMUNITY"]],
                   ranks_cd8) + labs(title="GOBP_B_CELL_MEDIATED_IMMUNITY(B)")
p3<-plotEnrichment(fgsea_sets1[["GOBP_ADAPTIVE_IMMUNE_RESPONSE"]],
                   ranks_cd8) + labs(title="GOBP_ADAPTIVE_IMMUNE_RESPONSE(B)")
pdf(file="H:/WRY/data/newRNA/SI/cluster/B/GSEA_M5_B_06h_plot_1.pdf",width = 10,height=6)
plot_grid(p1,p2,p3,ncol = 2)
dev.off()

new_use<-subset(n, subset = orig.ident == "s72h" | orig.ident == "s00h")
cd8<-subset(new_use, subset = cell_type == "B")
marker_cd8<-wilcoxauc(cd8,group_by = "orig.ident")
marker_cd81<-marker_cd8 %>% dplyr::filter(group == "s72h") %>% arrange(desc(auc)) %>% dplyr::select(feature,auc)
ranks_cd8<-deframe(marker_cd81)
p1<-plotEnrichment(fgsea_sets1[["GOCC_IMMUNOGLOBULIN_COMPLEX"]],
                   ranks_cd8) + labs(title="GOCC_IMMUNOGLOBULIN_COMPLEX(B)")
p2<-plotEnrichment(fgsea_sets1[["GOBP_B_CELL_MEDIATED_IMMUNITY"]],
                   ranks_cd8) + labs(title="GOBP_B_CELL_MEDIATED_IMMUNITY(B)")
p3<-plotEnrichment(fgsea_sets1[["GOBP_ADAPTIVE_IMMUNE_RESPONSE"]],
                   ranks_cd8) + labs(title="GOBP_ADAPTIVE_IMMUNE_RESPONSE(B)")
pdf(file="H:/WRY/data/newRNA/SI/cluster/B/GSEA_M5_B_72h_plot_1.pdf",width = 10,height=6)
plot_grid(p1,p2,p3,ncol = 2)
dev.off()


#########################
new_use<-subset(n, subset = orig.ident == "s02h" | orig.ident == "s00h")
cd8<-subset(new_use, subset = cell_type == "Plasma")
marker_cd8<-wilcoxauc(cd8,group_by = "orig.ident")
marker_cd81<-marker_cd8 %>% dplyr::filter(group == "s02h") %>% arrange(desc(auc)) %>% dplyr::select(feature,auc)
ranks_cd8<-deframe(marker_cd81)
p1<-plotEnrichment(fgsea_sets1[["GOBP_HUMORAL_IMMUNE_RESPONSE"]],
                   ranks_cd8) + labs(title="GOBP_HUMORAL_IMMUNE_RESPONSE(Plasma)")
p2<-plotEnrichment(fgsea_sets1[["GOBP_IMMUNOGLOBULIN_PRODUCTION"]],
                   ranks_cd8) + labs(title="GOBP_IMMUNOGLOBULIN_PRODUCTION(Plasma)")
p3<-plotEnrichment(fgsea_sets1[["GOCC_IMMUNOGLOBULIN_COMPLEX"]],
                   ranks_cd8) + labs(title="GOCC_IMMUNOGLOBULIN_COMPLEX(Plasma)")
pdf(file="H:/WRY/data/newRNA/SI/cluster/B/GSEA_M5_Plasma_02h_plot_1.pdf",width = 10,height=6)
plot_grid(p1,p2,p3,ncol = 2)
dev.off()

new_use<-subset(n, subset = orig.ident == "s06h" | orig.ident == "s00h")
cd8<-subset(new_use, subset = cell_type == "Plasma")
marker_cd8<-wilcoxauc(cd8,group_by = "orig.ident")
marker_cd81<-marker_cd8 %>% dplyr::filter(group == "s06h") %>% arrange(desc(auc)) %>% dplyr::select(feature,auc)
ranks_cd8<-deframe(marker_cd81)
p1<-plotEnrichment(fgsea_sets1[["GOBP_HUMORAL_IMMUNE_RESPONSE"]],
                   ranks_cd8) + labs(title="GOBP_HUMORAL_IMMUNE_RESPONSE(Plasma)")
p2<-plotEnrichment(fgsea_sets1[["GOBP_IMMUNOGLOBULIN_PRODUCTION"]],
                   ranks_cd8) + labs(title="GOBP_IMMUNOGLOBULIN_PRODUCTION(Plasma)")
p3<-plotEnrichment(fgsea_sets1[["GOCC_IMMUNOGLOBULIN_COMPLEX"]],
                   ranks_cd8) + labs(title="GOCC_IMMUNOGLOBULIN_COMPLEX(Plasma)")
pdf(file="H:/WRY/data/newRNA/SI/cluster/B/GSEA_M5_Plasma_06h_plot_1.pdf",width = 10,height=6)
plot_grid(p1,p2,p3,ncol = 2)
dev.off()

new_use<-subset(n, subset = orig.ident == "s72h" | orig.ident == "s00h")
cd8<-subset(new_use, subset = cell_type == "Plasma")
marker_cd8<-wilcoxauc(cd8,group_by = "orig.ident")
marker_cd81<-marker_cd8 %>% dplyr::filter(group == "s72h") %>% arrange(desc(auc)) %>% dplyr::select(feature,auc)
ranks_cd8<-deframe(marker_cd81)
p1<-plotEnrichment(fgsea_sets1[["GOBP_HUMORAL_IMMUNE_RESPONSE"]],
                   ranks_cd8) + labs(title="GOBP_HUMORAL_IMMUNE_RESPONSE(Plasma)")
p2<-plotEnrichment(fgsea_sets1[["GOBP_IMMUNOGLOBULIN_PRODUCTION"]],
                   ranks_cd8) + labs(title="GOBP_IMMUNOGLOBULIN_PRODUCTION(Plasma)")
p3<-plotEnrichment(fgsea_sets1[["GOCC_IMMUNOGLOBULIN_COMPLEX"]],
                   ranks_cd8) + labs(title="GOCC_IMMUNOGLOBULIN_COMPLEX(Plasma)")
pdf(file="H:/WRY/data/newRNA/SI/cluster/B/GSEA_M5_Plasma_72h_plot_1.pdf",width = 10,height=6)
plot_grid(p1,p2,p3,ncol = 2)
dev.off()


