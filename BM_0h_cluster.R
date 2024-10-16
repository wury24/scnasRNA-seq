library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)
library(harmony)
n <- readRDS("H:/WRY/data/newRNA/BM/cluster/combind_clean_harmony.rds")
new.cluster.ids <- c("Monocyte",
                     "Neutrophil",
                     "Neutrophil",
                     "Monocyte",
                     "Neutrophil",
                     "Neutrophil",
                     "Neutrophil",
                     "Neutrophil",
                     "B",
                     "Monocyte",
                     "Neutrophil",
                     "Neutrophil",
                     "DC",
                     "Erythroid",
                     "Pre B",
                     "Neutrophil",
                     "Neutrophil",
                     "HSPC",
                     "T",
                     "Erythroid",
                     "Erythroid",
                     "Basophil",
                     "Pre B",
                     "NK",
                     "Macrophage",
                     "Megakaryocyte",
                     "B",
                     "Macrophage")
names(new.cluster.ids) <- levels(n)
n <- RenameIdents(n, new.cluster.ids)
n$cell_type <- Idents(n)
levels <- c("HSPC", 
            "Monocyte", 
            "Macrophage", 
            "Neutrophil",
            "Basophil",
            "NK",
            "DC",
            "Pre B",
            "B",
            "T",
            "Megakaryocyte",
            "Erythroid")
Idents(n) <- factor(Idents(n), levels = levels)
n0<-subset(n, subset = orig.ident == "n00h")
info<-n0@meta.data
info<-info[info$cell_type != "Erythroid",]
n0<-subset(n0, subset = cell_type != "Erythroid")
DimPlot(n0, reduction = "umap",label = TRUE, label.size = 4, pt.size = 0.6)
total <- readRDS("H:/WRY/data/newRNA/BM/cluster/merge_clean.rds")
total <- total[,rownames(info)]
total$ori_type <- info[rownames(total@meta.data), "cell_type"]
total <- NormalizeData(total,normalization.method = "LogNormalize",scale.factor = 10000,assay="originalexp")
total <- FindVariableFeatures(total)
total <- ScaleData(total)
total <- RunPCA(total,features = VariableFeatures(object = total))
#total <- JackStraw(total , dims = 30, num.replicate = 50)
#total <- ScoreJackStraw(total , dims = 1:30)
#JackStrawPlot(total , dims = 1:30)
#ElbowPlot(total)
total <- FindNeighbors(total, dims = 1:20) 
total <- FindClusters(total,resolution =3.5)
total <- RunUMAP(total, dims = 1:20)

pdf(file = "H:/WRY/data/newRNA/BM/n00h/umap_n00h_noEry.pdf",width = 6,height = 5)
DimPlot(total, reduction = "umap",label = TRUE, label.size = 4, pt.size = 0.6)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/n00h/umap_n00h_ori_type_noEry.pdf",width = 6,height = 5)
DimPlot(total, reduction = "umap",label = TRUE, label.size = 4, pt.size = 0.6,group.by = "ori_type")
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/n00h/marker_Myloid.pdf",width = 10,height = 15)
FeaturePlot(total, features = c("Adgre4","Adgre1","Csf1r","Itgam","Itgax"), reduction = "umap", cols = c("grey", "blue"), pt.size = 0.3,order=TRUE)
dev.off()

new.cluster.ids <- c("Neutrophil",
                     "Neutrophil",
                     "Neutrophil",
                     "Neutrophil",
                     "Neutrophil",
                     "Neutrophil",
                     "B",
                     "Neutrophil",
                     "Neutrophil",
                     "Monocyte",
                     "Macrophage",
                     "T",
                     "DC",
                     "Neutrophil",
                     "Monocyte",
                     "Pre B",
                     "B",
                     "Neutrophil",
                     "Monocyte",
                     "Monocyte",
                     "HSPC", 
                     "Pre B", 
                     "HSPC",
                     "NK",
                     "B",
                     "Neutrophil",
                     "DC",
                     "T",
                     "Basophil",
                     "Macrophage",
                     "Macrophage",
                     "Pre B",
                     "Pre B")
names(new.cluster.ids) <- levels(total)
total <- RenameIdents(total, new.cluster.ids)
total$cell_type <- Idents(total)
levels <- c("HSPC", 
            "Monocyte", 
            "Macrophage", 
            "Neutrophil",
            "Basophil",
            "NK",
            "DC",
            "Pre B",
            "B",
            "T")
Idents(total) <- factor(Idents(total), levels = levels)
pdf(file = "H:/WRY/data/newRNA/BM/n00h/umap_n00h_annotion_noEry.pdf",width = 6,height = 5)
DimPlot(total, reduction = "umap",label = TRUE, label.size = 4, pt.size = 0.6)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/n00h/umap_n00h_annotion_noEry_nolabel.pdf",width = 5,height = 5)
DimPlot(total, reduction = "umap",label = F, label.size = 4, pt.size = 0.6)+NoLegend()
dev.off()
saveRDS(total, file = "H:/WRY/data/newRNA/BM/n00h/n00h_Seurat_cluster_noEry.rds")

pdf(file = "H:/WRY/data/newRNA/BM/n00h/marker_all.pdf",width = 20,height = 20)
FeaturePlot(total, features = c("Cd34","Mki67","Itgam","Ly6g","Csf1r","Adgre4","S100a8","Prss34","Gata2","Cd3e","Klrb1c","Ighd","Ighm","Jchain","Hba-a1"), reduction = "umap", cols = c("grey", "blue"),  pt.size = 0.3)
dev.off()
##########################################
Idents(total) <- total$cell_type
levels <- c("HSPC", 
            "Monocyte", 
            "Macrophage", 
            "Neutrophil",
            "Basophil",
            "NK",
            "DC",
            "Pre B",
            "B",
            "T")
Idents(total) <- factor(Idents(total), levels = levels)
marker <- FindAllMarkers(total, only.pos = TRUE,assay="originalexp",min.pct = 0.25,logfc.threshold = 0.25)
marker <- marker[marker$p_val_adj < 0.05,]
top10 <- marker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
total <- ScaleData(total,features = top10$gene)
data <- subset(total, downsample = 300)
pdf("H:/WRY/data/newRNA/BM/n00h/Heatmap_total_top10_noEry_nolegend1.pdf", width = 12, height = 8)
DoHeatmap(data, features = top10$gene, size = 3, angle = 90, label = FALSE)+NoLegend()+
  scale_fill_gradientn(colors = rev(brewer.pal(n = 10, name = "RdBu")))
dev.off()

new <- readRDS("H:/WRY/data/newRNA/BM/new/merge_new.rds")
new <- new[,rownames(total@meta.data)]
new$cell_type <- total@meta.data[rownames(new@meta.data), "cell_type"]
Idents(new) <- new$cell_type
levels <- c("HSPC", 
            "Monocyte", 
            "Macrophage", 
            "Neutrophil",
            "Basophil",
            "NK",
            "DC",
            "Pre B",
            "B",
            "T")
Idents(new) <- factor(Idents(new), levels = levels)
marker <- FindAllMarkers(new, only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
marker <- marker[marker$p_val_adj < 0.05,]
top10 <- marker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
new <- FindVariableFeatures(new)
new <- ScaleData(new,features = top10$gene)
data <- subset(new, downsample = 300)
pdf("H:/WRY/data/newRNA/BM/n00h/Heatmap_new_top10.pdf", width = 12, height = 8)
DoHeatmap(data, features = top10$gene, size = 3, angle = 90, label = TRUE)+
  scale_fill_gradientn(colors = brewer.pal(n = 9, name = "Reds"))
dev.off()

new <- readRDS("H:/WRY/data/newRNA/BM/new/merge_new.rds")
new <- new[,rownames(n0@meta.data)]
new <- FindVariableFeatures(new)
new <- ScaleData(new)
new <- RunPCA(new,features = VariableFeatures(object = new))
new <- JackStraw(new , dims = 30, num.replicate = 50)
new <- ScoreJackStraw(new , dims = 1:30)
JackStrawPlot(new , dims = 1:30)
ElbowPlot(new)
new <- FindNeighbors(new, dims = 1:10) 
new <- FindClusters(new,resolution = 2)
new <- RunUMAP(new, dims = 1:10)
pdf(file = "H:/WRY/data/newRNA/BM/n00h/umap_new_n00h.pdf",width = 6,height = 5)
DimPlot(new, reduction = "umap",label = TRUE, label.size = 4, pt.size = 0.6)
dev.off()
new$cell_type <- info[rownames(new@meta.data), "cell_type"]
Idents(new) <- new$cell_type
levels <- c("HSPC", 
            "Monocyte", 
            "Macrophage", 
            "Neutrophil",
            "Basophil",
            "NK",
            "DC",
            "Pre B",
            "B",
            "T",
            "Megakaryocyte",
            "Erythroid")
Idents(new) <- factor(Idents(new), levels = levels)
pdf(file = "H:/WRY/data/newRNA/BM/n00h/umap_new_n00h_celltype.pdf",width = 6,height = 5)
DimPlot(new, reduction = "umap",label = TRUE, label.size = 4, pt.size = 0.6)
dev.off()
ratio<-read.csv(file = "H:/WRY/data/newRNA/BM/ntr/new_ratio.csv",row.names = 1)
new$new_ratio<-ratio[rownames(new@meta.data),"new_ratio"]
pdf("H:/WRY/data/newRNA/BM/n00h/new_ratio_new_n00h.pdf", width = 6, height = 5)
FeaturePlot(new, features = "new_ratio", reduction = "umap", label = TRUE,order = TRUE,pt.size = 0.6)+
  scale_color_gradientn(colors = brewer.pal(n = 9, name = "PuRd"))+
  theme(plot.title = element_blank())
dev.off()


#########################
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)
n0<-readRDS(file = "H:/WRY/data/newRNA/BM/n00h/n00h_Seurat_cluster.rds")
n <- readRDS("H:/WRY/data/newRNA/BM/new/merge_new.rds")
n<-n[,rownames(n0@meta.data)]
n$cell_type<-n0@meta.data[rownames(n@meta.data),"cell_type"]

avg.new<-data.frame(rownames(n))
avg.total<-data.frame(rownames(n0))
names(avg.new)<-"gene"
names(avg.total)<-"gene"
rownames(avg.new)<-avg.new$gene
rownames(avg.total)<-avg.total$gene

gene_data <- data.frame(t(as.matrix(n@assays$RNA@data)),cluster=n$cell_type,check.names = F)
average_data <- aggregate(.~cluster, gene_data, mean)
cluster_name <- average_data[,1]
average_data <- apply(average_data[,2:ncol(average_data)],2,as.numeric)
rownames(average_data) <- cluster_name
average_data <- t(average_data)
average_data <-data.frame(average_data)
avg.new<-cbind(avg.new,average_data)
write.csv(avg.new,file="H:/WRY/data/newRNA/BM/n00h/avgnew_0h.csv")

gene_data <- data.frame(t(as.matrix(n0@assays$originalexp@data)),cluster=n0$cell_type,check.names = F)
average_data <- aggregate(.~cluster, gene_data, mean)
cluster_name <- average_data[,1]
average_data <- apply(average_data[,2:ncol(average_data)],2,as.numeric)
rownames(average_data) <- cluster_name
average_data <- t(average_data)
average_data <-data.frame(average_data)
avg.total<-cbind(avg.total,average_data)
write.csv(avg.total,file="H:/WRY/data/newRNA/BM/n00h/avgtotal_0h.csv")

###########
n0<-readRDS(file = "H:/WRY/data/newRNA/BM/n00h/n00h_Seurat_cluster.rds")
info<-cbind(n0@reductions$umap@cell.embeddings,n0@meta.data)
rownames(info)<-gsub("n00h_","",rownames(info))
write.csv(info,file="H:/WRY/data/newRNA/BM/n00h/dynamo/cellinfo.csv")

########################################
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)
n <- readRDS("H:/WRY/data/newRNA/BM/n00h/n00h_Seurat_cluster_noEry.rds")
info<-n@meta.data
rm(n)
n <- readRDS("H:/WRY/data/newRNA/BM/new/merge_new.rds")
gene_all<-rownames(n@assays$RNA@counts)
rm(n)
info0<-info[info$orig.ident == "n00h",]
bm00t <- read.delim("H:/WRY/data/newRNA/BM/process/BM00h_total.tsv", row.names = 1)
gene<-intersect(gene_all,rownames(bm00t))
names(bm00t)<-paste0("n00h_",names(bm00t))
bm00t<-bm00t[gene,rownames(info0)]
bm00t<-data.frame(t(bm00t))
names(bm00t)<-gene
bm00n <- read.delim("H:/WRY/data/newRNA/BM/process/BM_corrected/BM_0h_batch2.new_corrected_matrix2.tsv", row.names = 1)
names(bm00n)<-paste0("n00h_",names(bm00n))
bm00n<-bm00n[gene,rownames(info0)]
bm00n<-data.frame(t(bm00n))
names(bm00n)<-gene
bm00t$cell_type<-info0[rownames(bm00t),"cell_type"]
datat<-aggregate(bm00t[1:length(gene)],by=list(bm00t$cell_type),FUN = sum)
bm00n$cell_type<-info0[rownames(bm00n),"cell_type"]
datan<-aggregate(bm00n[1:length(gene)],by=list(bm00n$cell_type),FUN = sum)
rownames(datat)<-datat[,1]
datat<-datat[,2:(length(gene)+1)]
datat<-data.frame(t(datat))
rownames(datan)<-datan[,1]
datan<-datan[,2:(length(gene)+1)]
datan<-data.frame(t(datan))
data<-datan/datat
data[is.na(data)]<-0
which(data>1)
write.table(data,file="H:/WRY/data/newRNA/BM/n00h/ntr_gene_0hcluster_n00h.tsv",sep = "\t")
rm(bm00t)
rm(bm00n)
rm(datat)
rm(datan)
rm(data)
