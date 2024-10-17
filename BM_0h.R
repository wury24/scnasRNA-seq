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

#############################plot
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)
library(harmony)
n0<-readRDS(file = "H:/WRY/data/newRNA/BM/n00h/n00h_Seurat_cluster_noEry.rds")
###umap
cols_plot<-c("#F8766D","#DB8E00","#AEA200","#64B200","#00BD5C","#00C1A7","#00BADE",
             "#00A6FF","#B385FF","#EF67EB")
pdf(file = "H:/WRY/data/newRNA/BM/fig1/BM_umap_n00h_label.pdf",width = 6,height = 5)
DimPlot(n0, reduction = "umap",label = TRUE, cols = cols_plot,label.size = 4, pt.size = 0.6)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/fig1/BM_umap_n00h_nolegend.pdf",width = 5,height = 5)
DimPlot(n0, reduction = "umap",label = FALSE, cols = cols_plot, label.size = 4, pt.size = 0.6)+NoLegend()+
  theme(axis.title = element_text(size = 18))+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
dev.off()

###new ratio
ratio<-read.csv(file = "H:/WRY/data/newRNA/BM/ntr/new_ratio_qc.csv",row.names = 1)
n0$new_ratio<-ratio[rownames(n0@meta.data),"new_ratio"]
pdf(file = "H:/WRY/data/newRNA/BM/fig1/BM_Vln_newratio_n00h1.pdf",width = 12,height = 6)
VlnPlot(n0, features = c("new_ratio"),cols = cols_plot,pt.size=0.6)+NoLegend()+xlab("")+
  theme(plot.title = element_blank())+
  theme(axis.title = element_text(size = 18))+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=1.2, colour = 'black'))+
  theme(axis.ticks  = element_line(size=1.2, colour = 'black'))
dev.off()
pdf("H:/WRY/data/newRNA/BM/fig1/BM_new_ratio_n00h.pdf", width = 6, height = 5)
FeaturePlot(n0, features = "new_ratio", reduction = "umap", label = FALSE,order = TRUE,pt.size = 0.6)+
  scale_color_gradientn(colors = brewer.pal(n = 9, name = "PuRd"))
dev.off()
pdf("H:/WRY/data/newRNA/BM/fig1/BM_new_ratio_n00h_nolegend.pdf", width = 5, height = 5)
FeaturePlot(n0, features = "new_ratio", reduction = "umap", label = FALSE,order = TRUE,pt.size = 0.6)+NoLegend()+
  scale_color_gradientn(colors = brewer.pal(n = 9, name = "PuRd"))+
  theme(plot.title = element_blank())+
  theme(axis.title = element_text(size = 18))+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
dev.off()

color<-c("#313695","#4575B4","#74ADD1","#ABD9E9","#FDAE61","#F46D43","#D73027","#A50026")
n0@meta.data[n0@meta.data$new_ratio > 0.3,"new_ratio"] <- 0.3
pdf("H:/WRY/data/newRNA/BM/fig1/BM_new_ratio_n00h_test1.pdf", width = 5.5, height = 5)
FeaturePlot(n0, features = "new_ratio", reduction = "umap", label = FALSE,order = TRUE,pt.size = 0.6)+
  scale_color_gradientn(colors = color)+
  theme(plot.title = element_blank())+
  theme(axis.title = element_text(size = 18))+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
dev.off()
pdf("H:/WRY/data/newRNA/BM/fig1/BM_new_ratio_n00h_test1_nolegend.pdf", width = 5, height = 5)
FeaturePlot(n0, features = "new_ratio", reduction = "umap", label = FALSE,order = TRUE,pt.size = 0.6)+NoLegend()+
  scale_color_gradientn(colors = color)+
  theme(plot.title = element_blank())+
  theme(axis.title = element_text(size = 18))+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
dev.off()
####heatmap
new <- readRDS("H:/WRY/data/newRNA/BM/new/merge_new.rds")
new <- new[,rownames(n0@meta.data)]
new$cell_type <- n0@meta.data[rownames(new@meta.data), "cell_type"]
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
marker <- FindAllMarkers(new, only.pos = TRUE,min.pct = 0.1,logfc.threshold = 0.25)
marker <- marker[marker$p_val_adj < 0.05,]
write.csv(marker,file = "H:/WRY/data/newRNA/BM/fig1/BM_marker_new_n00h.csv")
top10 <- marker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,file = "H:/WRY/data/newRNA/BM/fig1/BM_marker_new_top10_n00h.csv")
new <- FindVariableFeatures(new)
new <- ScaleData(new,features = top10$gene)
data <- subset(new, downsample = 300)
pdf("H:/WRY/data/newRNA/BM/fig1/BM_Heatmap_new_top10.pdf", width = 12, height = 20)
DoHeatmap(data, features = top10$gene, size = 3, group.colors = cols_plot, angle = 90, label = TRUE)+
  scale_fill_gradientn(colors = brewer.pal(n = 9, name = "Reds"))
dev.off()

library(ComplexHeatmap)
library(circlize)
datah<-data@assays$RNA@scale.data
gene<-unique(top10$gene)
cluster_info <- sort(Idents(data))
datah<-datah[gene,names(cluster_info)]
gene_select<-c("Adgrg1","Myb","Ccr2","Ms4a6c","Ctss","Ccl6","Ltf","Mmp8","Mcpt8","Gata2","Gzma","Klra8","Siglech",
                "Runx2","Vpreb3","Bach2","Ighm","Ighd","Trbc2","Tcf7","Car2","Hbb-bs")
gene_select<-intersect(gene_select,top10$gene)
gene_pos<-which(rownames(datah) %in% gene_select)
row_anno<-rowAnnotation(gene=anno_mark(at=gene_pos,labels = gene_select,side="left"))
col_fun <- colorRamp2(seq(-0.5,2.5, length.out =9),brewer.pal(n = 9, name = "Reds"))
p1<-DimPlot(n0, reduction = "umap",label = TRUE, label.size = 4, pt.size = 0.6,combine=FALSE)
color<-p1[[1]]$layers[[2]]$data$color
names(color)<-levels(cluster_info)
top_anno<- HeatmapAnnotation(cluster=anno_block(gp = gpar(fill = color),
                                                labels = levels(cluster_info),
                                                labels_gp = gpar(col = "white", fontsize = 8,fontface="bold",lineheight=1)))

pdf("H:/WRY/data/newRNA/BM/fig1/BM_Heatmap_new_test.pdf", width = 12, height = 8)
Heatmap(datah,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_title = NULL,
        column_split = cluster_info,
        left_annotation = row_anno,
        col = col_fun,
        top_annotation = top_anno,
        use_raster = FALSE)
dev.off()



pdf("H:/WRY/data/newRNA/BM/fig1/BM_Heatmap_new_top10_nolegend.pdf", width = 12, height = 8)
DoHeatmap(data, features = top10$gene, size = 3, angle = 90, group.colors = cols_plot, label = FALSE)+NoLegend()+
  scale_fill_gradientn(colors = brewer.pal(n = 9, name = "Reds"))
dev.off()

total <- readRDS("H:/WRY/data/newRNA/BM/cluster/merge_clean.rds")
total <- total[,rownames(n0@meta.data)]
total$cell_type <- n0@meta.data[rownames(total@meta.data), "cell_type"]
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
total <- NormalizeData(total,normalization.method = "LogNormalize",scale.factor = 10000,assay="originalexp")
marker <- FindAllMarkers(total, only.pos = TRUE,assay="originalexp",min.pct = 0.25,logfc.threshold = 0.25)
marker <- marker[marker$p_val_adj < 0.05,]
write.csv(marker,file = "H:/WRY/data/newRNA/BM/fig1/BM_marker_total_n00h.csv")
top10 <- marker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,file = "H:/WRY/data/newRNA/BM/fig1/BM_marker_total_top10_n00h.csv")
total <- FindVariableFeatures(total)
total <- ScaleData(total,features = top10$gene)
data <- subset(total, downsample = 300)
pdf("H:/WRY/data/newRNA/BM/fig1/BM_Heatmap_total_top10.pdf", width = 12, height = 20)
DoHeatmap(data, features = top10$gene, group.colors = cols_plot, size = 3, angle = 90, label = TRUE)+
  scale_fill_gradientn(colors = rev(brewer.pal(n = 10, name = "RdBu")))
dev.off()
pdf("H:/WRY/data/newRNA/BM/fig1/BM_Heatmap_total_top10_nolegend.pdf", width = 11.5, height = 12)
DoHeatmap(data, features = top10$gene, group.colors = cols_plot, size = 3, angle = 90, label = FALSE)+NoLegend()+
  scale_fill_gradientn(colors = rev(brewer.pal(n = 10, name = "RdBu")))
dev.off()

##########new cluster
new <- readRDS("H:/WRY/data/newRNA/BM/new/merge_new.rds")
new <- new[,rownames(n0@meta.data)]
new <- FindVariableFeatures(new)
#new <- FindVariableFeatures(new, selection.method="mvp",mean.cutoff = c(0.1,Inf),dispersion.cutoff=c(0.5,Inf))
#VariableFeaturePlot(new)
new <- ScaleData(new)
new <- RunPCA(new,features = VariableFeatures(object = new))
#new <- JackStraw(new , dims = 30, num.replicate = 50)
#new <- ScoreJackStraw(new , dims = 1:30)
#JackStrawPlot(new , dims = 1:30)
#ElbowPlot(new)
new <- FindNeighbors(new, dims = 1:20) 
new <- FindClusters(new,resolution = 0.8)
new <- RunUMAP(new, dims = 1:20)
pdf(file = "H:/WRY/data/newRNA/BM/fig1/BM_umap_new_n00h.pdf",width = 6,height = 5)
DimPlot(new, reduction = "umap",label = TRUE, label.size = 4, pt.size = 0.6)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/fig1/BM_umap_new_n00h_nolegend.pdf",width = 5,height = 5)
DimPlot(new, reduction = "umap",label = FALSE,  pt.size = 0.6)+NoLegend()+
  theme(axis.title = element_text(size = 18))+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
dev.off()
new$cell_type <- n0@meta.data[rownames(new@meta.data), "cell_type"]
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
pdf(file = "H:/WRY/data/newRNA/BM/fig1/BM_umap_new_n00h_celltype.pdf",width = 6,height = 5)
DimPlot(new, reduction = "umap",label = TRUE,cols=cols_plot,label.size = 4, pt.size = 0.6)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/fig1/BM_umap_new_n00h_celltype_nolegend.pdf",width = 5,height = 5)
DimPlot(new, reduction = "umap",label = FALSE, cols=cols_plot,label.size = 4, pt.size = 0.6)+NoLegend()+
  theme(axis.title = element_text(size = 18))+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
dev.off()

predictions <- table(new$seurat_clusters,Idents(new))
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
order<-c("7","9","2","11","0","1","3","5","12","8","10","6","4")
predictions<-predictions[order,]
predictions <- as.data.frame(predictions)
pdf(file = "H:/WRY/data/newRNA/BM/fig1/BM_fraction_new_total.pdf",width = 7.5,height = 5)
ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile()+
  scale_fill_gradient(name = "Fraction of cells",low = "#ffffc8", high = "#7d0025") +
  xlab("New cluster") + ylab("Total cluster") +
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

predictions <- table(Idents(new),new$seurat_clusters)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
order<-c("7","9","2","11","0","1","3","5","12","8","10","6","4")
type<-c("HSPC", 
        "Monocyte", 
        "Macrophage", 
        "Neutrophil",
        "Basophil",
        "NK",
        "DC",
        "Pre B",
        "B",
        "T")
predictions<-predictions[,order]
predictions <- as.data.frame(predictions)
pdf(file = "H:/WRY/data/newRNA/BM/fig1/BM_fraction_total_new.pdf",width = 5,height = 5)
ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile()+
  scale_fill_gradient(name = "Fraction of cells",low = "#ffffc8", high = "#7d0025") +
  xlab("Total RNA cluster") + ylab("New RNA cluster") +
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

ratio<-read.csv(file = "H:/WRY/data/newRNA/BM/ntr/new_ratio_qc.csv",row.names = 1)
new$new_ratio<-ratio[rownames(new@meta.data),"new_ratio"]
color<-c("#313695","#4575B4","#74ADD1","#ABD9E9","#FDAE61","#F46D43","#D73027","#A50026")
new@meta.data[new@meta.data$new_ratio > 0.3,"new_ratio"] <- 0.3
pdf("H:/WRY/data/newRNA/BM/fig1/BM_new_ratio_n00h_newcluster.pdf", width = 5.5, height = 5)
FeaturePlot(new, features = "new_ratio", reduction = "umap", label = FALSE,order = TRUE,pt.size = 0.6)+
  scale_color_gradientn(colors = color)+
  theme(plot.title = element_blank())+
  theme(axis.title = element_text(size = 18))+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  scale_y_continuous(breaks = c(-3,0,3,6))
dev.off()
pdf("H:/WRY/data/newRNA/BM/fig1/BM_new_ratio_n00h_newcluster_nolegend.pdf", width = 5.1, height = 5)
FeaturePlot(new, features = "new_ratio", reduction = "umap", label = FALSE,order = TRUE,pt.size = 0.6)+NoLegend()+
  scale_color_gradientn(colors = color)+
  theme(plot.title = element_blank())+
  theme(axis.title = element_text(size = 18))+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  scale_y_continuous(breaks = c(-3,0,3,6))
dev.off()
####################
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(RColorBrewer)
total<-readRDS(file = "H:/WRY/data/newRNA/BM/n00h/n00h_Seurat_cluster_noEry.rds")
new <- readRDS("H:/WRY/data/newRNA/BM/new/merge_new.rds")
new<-new[,rownames(total@meta.data)]
new$cell_type <- total@meta.data[rownames(new@meta.data), "cell_type"]
Idents(new)<-new$cell_type
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
type <- c("HSPC", 
          "Monocyte", 
          "Macrophage", 
          "Neutrophil",
          "Basophil",
          "NK",
          "DC",
          "Pre B",
          "B",
          "T")
stage<-c("Total","New")
cat<-c()
for(i in 1:length(type)){
  for (j in 1:length(stage)) {
    cat<-c(cat,paste(type[i],stage[j],sep=""))
  }
}
Vlntype <- function(x){
  gene_n<-data.frame(new@assays$RNA@data[x,rownames(total@meta.data)],
                     new@meta.data[rownames(total@meta.data),"cell_type"])
  names(gene_n)<-c("Exprssion","cell_type")
  gene_t<-data.frame(total@assays$originalexp@data[x,rownames(total@meta.data)],
                     total@meta.data[rownames(total@meta.data),"cell_type"])
  names(gene_t)<-c("Exprssion","cell_type")
  gene_t$RNA <- "Total"
  gene_t$type <- paste(gene_t$cell_type,gene_t$RNA,sep = "")
  gene_n$RNA <- "New"
  gene_n$type <- paste(gene_n$cell_type,gene_n$RNA,sep = "")
  vln <- rbind(gene_t,gene_n) 
  vln$type <- factor(vln$type, levels = cat)
  print(table(vln$type))
  p <- ggplot(vln,aes(x = type, y = Exprssion))+geom_violin(aes(color=RNA,fill = RNA),trim = TRUE,scale = "width")+
    #    geom_jitter(size=0.1,width = 0.1)+
    scale_fill_manual(values = c("#e60012","#009a3e"))+scale_color_manual(values = c("#e60012","#009a3e"))+
    theme_bw()+theme_classic()+
    theme(axis.text = element_text(size = 27,color="black"))+
    theme(axis.title = element_blank())+
    theme(axis.line=element_line(linetype=1,color="black",size=1.1))+
    theme(axis.ticks=element_line(color="black",size=1.1,lineend = 10))+
    labs(title = x)+theme(plot.title=element_text(hjust=0.5))+
    scale_x_discrete(breaks=cat,
                     labels=c("Total", "New", "Total", "New","Total", "New","Total", "New","Total", "New",
                              "Total", "New", "Total", "New","Total", "New","Total", "New","Total", "New"))
  return(p)
}

p1<-Vlntype("Ptprc")
p2<-Vlntype("Ccr2")
p3<-Vlntype("S100a9")
p4<-Vlntype("Rpl13")

pdf(file="H:/WRY/data/newRNA/BM/fig1/Vln_corgene.pdf",width = 16,height=10)
plot_grid(p1,p2,p3,p4,ncol = 1)
dev.off()

#############GO ntr
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(Seurat)
new<-read.csv(file="H:/WRY/data/newRNA/BM/ntr/gene_celltype_percent_new_0h.csv",row.names = 1)
ntr<-read.csv(file="H:/WRY/data/newRNA/BM/ntr/gene_celltype_NTR_0h.csv",row.names = 1)
type<-c("HSPC", 
        "Monocyte", 
        "Macrophage", 
        "Neutrophil",
        "Basophil",
        "NK",
        "DC",
        "Pre.B",
        "B",
        "T",
        "Erythroid")
go<-read.csv(file="H:/WRY/data/newRNA/BM/fig1/go_select.csv")
GO0h<-data.frame()
gene_use<-c()
for(i in 1:length(type)){
  percent<-data.frame(cbind(rownames(new),new[,type[i]]))
  names(percent)<-c("gene","pct")
  percent<-percent[percent$pct>0.05,]
  gene_use<-c(gene_use,length(rownames(percent)))
  
  data<-ntr[percent$gene,]
  data<-data.frame(cbind(rownames(data),data[,type[i]]))
  names(data)<-c("gene","ntr")
  data1<-data[order(data$ntr, decreasing = TRUE),]
  data1<-data1[c(1:150),]
  input <- bitr(data1$gene, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db",drop = TRUE)
  if(length(rownames(input))>0){
    ego <- enrichGO(gene = input$ENTREZID, 
                    OrgDb = org.Mm.eg.db, 
                    ont = "BP", 
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.05,
                    readable = TRUE)
    if(length(ego) > 0){
      GO <- ego@result
      id<-intersect(GO$Description,go$Description)
      rownames(GO)<-GO$Description
      GO <- GO[id,c("ID","Description","GeneRatio","qvalue","Count")]
      rownames(GO)<-NULL
    }else{
      GO<-NULL
    }
    if(length(rownames(GO))>0){
      GO$cell_type<-type[i]
      GO0h<-rbind(GO0h,GO)
    }
  }
}


GO0h$cell_type <- factor(GO0h$cell_type,levels = type)
GO0h<-na.omit(GO0h)
order_GO<-go[c(15:1),3]
GO0h$Description <- factor(GO0h$Description,levels = order_GO)


mytheme <- theme(axis.title=element_text(face="bold", size=14,colour = 'black'), 
                 axis.text=element_text(face="bold", size=14,colour = 'black'), 
                 axis.line = element_line(size=0.5, colour = 'black'), 
                 panel.background = element_rect(color='black'), 
                 legend.key = element_blank() 
)


p <- ggplot(GO0h,aes(x=cell_type,y=Description,color=-log10(qvalue),size=Count))+
  geom_point()+
  #  scale_size(range=c(2, 20))+
  scale_colour_gradient2(low = "#f8f0f7",mid="#3ca0c8",high = "#003c28",midpoint = 4)+
  #  scale_colour_gradient2(low = "#91a5cd",mid="#ffffff",high = "#b42832",midpoint = 5)+
  #  scale_colour_gradient(low = "white",high = "red")+
  theme_classic()+
  ylab("GO Pathway Terms")+
  xlab("")+
  labs(color=expression(-log[10](qValue)))+
  theme(legend.title=element_text(size=14), legend.text = element_text(size=14))+
  theme(axis.title.y = element_text(margin = margin(r = 50)),axis.title.x = element_text(margin = margin(t = 20)))+
  theme(axis.text.x = element_text(face ="bold",color="black",angle=90))+
  theme(axis.ticks.x = element_blank() )
plot <- p+mytheme
pdf(file="H:/WRY/data/newRNA/BM/fig1/BM_GO_n00h_new.pdf",width=10,height=8)
plot
dev.off()

#################################subtype
n0<-readRDS(file = "H:/WRY/data/newRNA/BM/n00h/n00h_Monocyte_Seurat_cluster.rds")
###umap
pdf(file = "H:/WRY/data/newRNA/BM/fig1/BM_Monocyte_umap_n00h_label.pdf",width = 6,height = 5)
DimPlot(n0, reduction = "umap",label = TRUE, label.size = 4, pt.size = 0.6)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/fig1/BM_Monocyte_umap_n00h_nolegend.pdf",width = 5,height = 5)
DimPlot(n0, reduction = "umap",label = FALSE, label.size = 4, pt.size = 0.6)+NoLegend()+
  theme(axis.title = element_text(size = 18))+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
dev.off()

pdf("H:/WRY/data/newRNA/BM/fig1/BM_Monocyte_Mki67_n00h.pdf", width = 6, height = 5)
FeaturePlot(n0, features = c("Mki67"), reduction = "umap", cols = c("grey", "blue"), pt.size = 0.6)
dev.off()
pdf("H:/WRY/data/newRNA/BM/fig1/BM_Monocyte_Mki67_n00h_nolegend.pdf", width = 5, height = 5)
FeaturePlot(n0, features = c("Mki67"), reduction = "umap", cols = c("grey", "blue"), pt.size = 0.6)+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title = element_text(size = 18))+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
dev.off()
###new ratio
ratio<-read.csv(file = "H:/WRY/data/newRNA/BM/ntr/new_ratio.csv",row.names = 1)
n0$new_ratio<-ratio[rownames(n0@meta.data),"new_ratio"]
pdf("H:/WRY/data/newRNA/BM/fig1/BM_Monocyte_new_ratio_n00h.pdf", width = 6, height = 5)
FeaturePlot(n0, features = "new_ratio", reduction = "umap", label = FALSE,order = TRUE,pt.size = 0.6)+
  scale_color_gradientn(colors = brewer.pal(n = 9, name = "PuRd"))
dev.off()
pdf("H:/WRY/data/newRNA/BM/fig1/BM_Monocyte_new_ratio_n00h_nolegend.pdf", width = 5, height = 5)
FeaturePlot(n0, features = "new_ratio", reduction = "umap", label = FALSE,order = TRUE,pt.size = 0.6)+NoLegend()+
  scale_color_gradientn(colors = brewer.pal(n = 9, name = "PuRd"))+
  theme(plot.title = element_blank())+
  theme(axis.title = element_text(size = 18))+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
dev.off()

#################
n0<-readRDS(file = "H:/WRY/data/newRNA/BM/n00h/n00h_Neutrophil_Seurat_cluster.rds")
###umap
pdf(file = "H:/WRY/data/newRNA/BM/fig1/BM_Neutrophil_umap_n00h_label.pdf",width = 6,height = 5)
DimPlot(n0, reduction = "umap",label = TRUE, label.size = 4, pt.size = 0.6)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/fig1/BM_Neutrophil_umap_n00h_nolegend.pdf",width = 5,height = 5)
DimPlot(n0, reduction = "umap",label = FALSE, label.size = 4, pt.size = 0.6)+NoLegend()+
  theme(axis.title = element_text(size = 18))+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
dev.off()

pdf("H:/WRY/data/newRNA/BM/fig1/BM_Neutrophil_Mki67_n00h.pdf", width = 6, height = 5)
FeaturePlot(n0, features = c("Mki67"), reduction = "umap", cols = c("grey", "blue"), pt.size = 0.6)
dev.off()
pdf("H:/WRY/data/newRNA/BM/fig1/BM_Neutrophil_Mki67_n00h_nolegend.pdf", width = 5, height = 5)
FeaturePlot(n0, features = c("Mki67"), reduction = "umap", cols = c("grey", "blue"), pt.size = 0.6)+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title = element_text(size = 18))+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
dev.off()
###new ratio
ratio<-read.csv(file = "H:/WRY/data/newRNA/BM/ntr/new_ratio.csv",row.names = 1)
n0$new_ratio<-ratio[rownames(n0@meta.data),"new_ratio"]
pdf("H:/WRY/data/newRNA/BM/fig1/BM_Neutrophil_new_ratio_n00h.pdf", width = 6, height = 5)
FeaturePlot(n0, features = "new_ratio", reduction = "umap", label = FALSE,order = TRUE,pt.size = 0.6)+
  scale_color_gradientn(colors = brewer.pal(n = 9, name = "PuRd"))
dev.off()
pdf("H:/WRY/data/newRNA/BM/fig1/BM_Neutrophil_new_ratio_n00h_nolegend.pdf", width = 5, height = 5)
FeaturePlot(n0, features = "new_ratio", reduction = "umap", label = FALSE,order = TRUE,pt.size = 0.6)+NoLegend()+
  scale_color_gradientn(colors = brewer.pal(n = 9, name = "PuRd"))+
  theme(plot.title = element_blank())+
  theme(axis.title = element_text(size = 18))+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
dev.off()

###########################
library(ggrepel)
t<-read.csv(file="H:/WRY/data/newRNA/BM/fig1/var_n00h_1.csv",row.names = 1)
t<-t[,c("NAME","nCells","nCounts","alpha","gamma","half_life")]
t<-na.omit(t)
t<-t[t$nCells >200,]
data<-t
n<-floor(length(rownames(data))*0.1)
data1<-data[order(data$half_life,decreasing = TRUE),]
h_cutoff<-data1[n,"half_life"]
data2<-data[order(data$alpha,decreasing = TRUE),]
a_cutoff<-data2[n,"alpha"]
data$color<-"grey"
data[data$half_life > h_cutoff & data$alpha > a_cutoff,"color"]<-"yellow"
data[data$half_life > h_cutoff & data$alpha < a_cutoff,"color"]<-"green"
data[data$half_life < h_cutoff & data$alpha > a_cutoff,"color"]<-"blue"
color<- c(blue='#4d85bd', grey = 'grey', yellow ='#f7903d',green="#59a95a")

pdf(file="H:/WRY/data/newRNA/BM/fig1/half_alpha_point1.pdf",width=6,height=5)
ggplot(data,aes(log(half_life,base = 2),log(alpha,base = 2),color=color))+geom_point()+
  theme_classic()+
  scale_color_manual(values=color)+
  geom_hline(yintercept = log(a_cutoff,base = 2), lty=4,col='grey',lwd=0.6)+
  geom_vline(xintercept = log(h_cutoff,base = 2), lty=4,col='grey',lwd=0.6)+
  geom_text_repel(data=data[data$NAME == "Ptprc" |
                              data$NAME == "S100a9" |
                              data$NAME == "Ccr2" |
                              data$NAME == "Mki67" |
                              data$NAME == "Camp" |
                              data$NAME == "Rpl13" |
                              data$NAME == "Rpl29", ],
                  aes(label=NAME),size=6,color="black")+
  NoLegend()
dev.off()
