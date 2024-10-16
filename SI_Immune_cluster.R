library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)
library(harmony)
cc<-read.csv("H:/WRY/code/ccgene/mouse_cc.csv",row.names = 1)
s.genes <- cc[cc$stage == "S","symbol"]
g2m.genes <- cc[cc$stage == "G2/M","symbol"]

n <- readRDS("H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_cluster.rds")
meta<-n@meta.data
metaI<-meta[meta$cell_type == "Immune",]
rm(n)

n <- readRDS("H:/WRY/data/newRNA/SI/cluster/merge_clean.rds")
n<-n[,rownames(metaI)]
n@meta.data<-n@meta.data[rownames(n@meta.data),c("orig.ident","nCount_originalexp","nFeature_originalexp","percent.mt","percent.ribo","percent.dissociation")]
VlnPlot(n, features = c("nFeature_originalexp", "nCount_originalexp", "percent.mt"), ncol = 3)
n <- NormalizeData(n,normalization.method = "LogNormalize",scale.factor = 10000,assay="originalexp")
n <- FindVariableFeatures(n)
n <- CellCycleScoring(n, s.features = s.genes, g2m.features = g2m.genes)
n$CC.Difference <- n$S.Score - n$G2M.Score
n <- ScaleData(n,vars.to.regress = c('nCount_originalexp', 'percent.mt', 'percent.ribo', 'percent.dissociation', 'CC.Difference'))
n <- RunPCA(n,npcs = 40)
n <- RunHarmony(n, group.by.vars = "orig.ident",assay.use = "originalexp",plot_convergence = TRUE)
n <- RunUMAP(n, reduction = "harmony", dims = 1:30)
n <- FindNeighbors(n, reduction = "harmony", dims = 1:30) 
n <- FindClusters(n,resolution=2)
pdf(file = "H:/WRY/data/newRNA/SI/cluster/umap_harmony_Immune.pdf",width = 6,height = 5)
DimPlot(n, reduction = "umap",label = TRUE)
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/cluster/umap_ori_harmony_Immune.pdf",width = 21,height = 5)
DimPlot(n, reduction = "umap", split.by = "orig.ident")
dev.off()
saveRDS(n, file = "H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_Immune.rds")
###############
n <- readRDS("H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_Immune.rds")
########
new.cluster.ids <- c("Plasma",
                     "Plasma",
                     "Plasma", 
                     "CD8+ T",
                     "Plasma",
                     "B",
                     "B",
                     "CD8+ T",
                     "Plasma",
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
                     "Plasma",
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
pdf(file = "H:/WRY/data/newRNA/SI/cluster/umap_annotion_Immune1.pdf",width = 7,height = 5)
DimPlot(n, reduction = "umap",label = TRUE, label.size = 4, pt.size = 0.15,cols = cols)
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/cluster/umap_annotion_ori_Immune1.pdf",width = 21,height = 5)
DimPlot(n, reduction = "umap",label = TRUE, label.size = 4, pt.size = 0.15,split.by = "orig.ident",cols = cols)
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/cluster/umap_annotion_ori_Nolegend_Immune1.pdf",width = 21,height = 5)
DimPlot(n, reduction = "umap",label = F, pt.size = 0.15,split.by = "orig.ident",cols = cols)
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/cluster/umap_annotion_Nolegend_Immune1.pdf",width = 5,height = 5)
DimPlot(n, reduction = "umap",label = FALSE, label.size = 4, pt.size = 0.15,cols = cols)+NoLegend()
dev.off()
saveRDS(n, file = "H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_cluster_Immune1.rds")

n.combined <- readRDS("H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_cluster_Immune1.rds")
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
write.csv(avg.new,file="H:/WRY/data/newRNA/SI/cluster/avg/avgnew_Immune1.csv")
write.csv(avg.total,file="H:/WRY/data/newRNA/SI/cluster/avg/avgtotal_Immune1.csv")

library(pheatmap)
input<-read.csv("H:/WRY/data/newRNA/SI/cluster/Induced/Ig1.csv",header = FALSE)
total<-read.csv("H:/WRY/data/newRNA/SI/cluster/avg/avgtotal_Immune1.csv",row.names = 1)
stage <- c("s00h","s02h","s06h","s72h")
type <- c("B",
          "Plasma",
          "Induced.plasma")
cat<-c()
for(i in 1:length(type)){
  for (j in 1:length(stage)) {
    cat<-c(cat,paste(type[i],stage[j],sep=""))
  }
}
data<-total[input[,1],cat]
sc<-data.frame(t(scale(t(data))))
annotation_col = data.frame(stage = factor(rep(stage, 3)), cell_type = factor(rep(type,each = 4)))
rownames(annotation_col)<-names(sc)
pdf(file="H:/WRY/data/newRNA/SI/cluster/Induced/Ig1_total.pdf",width=4,height = 3)
p<-pheatmap(sc, scale = "none",cluster_cols = FALSE, cluster_rows = TRUE, fontsize = 8,
            show_rownames = TRUE, show_colnames = TRUE, legend = TRUE,
            color = colorRampPalette(c("#0000FF", "white", "#DD0000"))(100)) 
dev.off()


library(Seurat)
library(dplyr)
library(patchwork)
n <- readRDS("H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_cluster_Immune1.rds")
levels <- c("CD4+ T", 
            "CD8+ T",
            "B",
            "Plasma",
            "Macrophage",
            "DC",
            "Mast")
stage <- c("s00h","s02h","s06h","s72h")
main_dir <- "H:/WRY/data/newRNA/SI/cluster/marker1/stage_total"
for (i in 1:length(levels)){
  dir.create(file.path(main_dir,levels[i]))
  for(j in 2:length(stage)){
    marker<-FindMarkers(n,subset.ident = levels[i],ident.1 = stage[j],ident.2 = stage[1],group.by = "orig.ident",only.pos = FALSE)
    write.csv(marker,file =  paste("H:/WRY/data/newRNA/SI/cluster/marker1/stage_total/",levels[i],"/","marker_",stage[j],"_",stage[1],"_total.csv",sep=""))
  }
}

n <- readRDS("H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_cluster_Immune1.rds")
info<-n@meta.data
rm(n)
n <- readRDS("H:/WRY/data/newRNA/SI/new/merge_new.rds")
n$cell_type<-info[rownames(n@meta.data),"cell_type"]
Idents(n)<-n$cell_type
levels <- c("CD4+ T", 
            "CD8+ T",
            "B",
            "Plasma",
            "Macrophage",
            "DC",
            "Mast")
stage <- c("s00h","s02h","s06h","s72h")
main_dir <- "H:/WRY/data/newRNA/SI/cluster/marker1/stage_new"
for (i in 1:length(levels)){
  dir.create(file.path(main_dir,levels[i]))
  for(j in 2:length(stage)){
    marker<-FindMarkers(n,subset.ident = levels[i],ident.1 = stage[j],ident.2 = stage[1],group.by = "orig.ident",only.pos = FALSE)
    write.csv(marker,file =  paste("H:/WRY/data/newRNA/SI/cluster/marker1/stage_new/",levels[i],"/","marker_",stage[j],"_",stage[1],"_new.csv",sep=""))
  }
}

n <- readRDS("H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_cluster_Immune1.rds")
mycolor<-cols
cell.prop<-data.frame(prop.table(table(Idents(n), n$orig.ident)))
colnames(cell.prop)<-c("cell_type","stage","proportion")
cell_count<-data.frame(table(Idents(n), n$orig.ident))
write.csv(cell_count, file = "H:/WRY/data/newRNA/SI/cluster/cell_count/cell_count_Immune1.csv")
pdf("H:/WRY/data/newRNA/SI/cluster/cell_count/prop_orig_Immune1.pdf", width = 5, height = 5)
ggplot(cell.prop,aes(stage,proportion,fill=cell_type))+
  geom_bar(stat="identity",position="fill")+scale_fill_manual(values=mycolor)+
  ggtitle("")+theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20))+
  theme_bw()+
  guides(fill=guide_legend(title=NULL))
dev.off()

mycolor<-c(brewer.pal(4,"Paired"))
cell.prop<-data.frame(prop.table(table(n$orig.ident,Idents(n))))
colnames(cell.prop)<-c("stage","cell_type","proportion")
pdf("H:/WRY/data/newRNA/SI/cluster/cell_count/prop_celltype_Immune1.pdf", width = 14, height = 5)
ggplot(cell.prop,aes(cell_type,proportion,fill=stage))+
  geom_bar(stat="identity",position="fill")+scale_fill_manual(values=mycolor)+
  ggtitle("")+theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20))+
  theme_bw()+
  guides(fill=guide_legend(title=NULL))
dev.off()

########
n <- readRDS("H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_cluster_Immune1.rds")
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
gene<-read.csv(file = "H:/WRY/data/newRNA/SI/cluster/SI_Immune_gene.csv",header = F)
pdf(file = "H:/WRY/data/newRNA/SI/cluster/SI_Immune_dotplot.pdf",width =8,height = 6.5)
DotPlot(new,features = rev(gene[,1]),cols = c("lightgrey","#e20000")) + RotatedAxis() + coord_flip()
dev.off()
