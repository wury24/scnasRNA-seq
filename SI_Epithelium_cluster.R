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
metaE<-meta[meta$cell_type == "E",]
rm(n)
n <- readRDS("H:/WRY/data/newRNA/SI/cluster/merge_clean.rds")
n<-n[,rownames(metaE)]
n@meta.data<-n@meta.data[rownames(n@meta.data),c("orig.ident","nCount_originalexp","nFeature_originalexp","percent.mt","percent.ribo","percent.dissociation")]
#VlnPlot(n, features = c("nFeature_originalexp", "nCount_originalexp", "percent.mt"), ncol = 3)
n <- NormalizeData(n,normalization.method = "LogNormalize",scale.factor = 10000,assay="originalexp")
n <- FindVariableFeatures(n)
n <- CellCycleScoring(n, s.features = s.genes, g2m.features = g2m.genes)
n$CC.Difference <- n$S.Score - n$G2M.Score
n <- ScaleData(n,vars.to.regress = c('nCount_originalexp', 'percent.mt', 'percent.ribo', 'percent.dissociation', 'CC.Difference'))
n <- RunPCA(n)
n <- RunHarmony(n, group.by.vars = "orig.ident",assay.use = "originalexp",plot_convergence = TRUE)
n <- RunUMAP(n, reduction = "harmony", dims = 1:40)
n <- FindNeighbors(n, reduction = "harmony", dims = 1:40) 
n <- FindClusters(n,resolution=2)
pdf(file = "H:/WRY/data/newRNA/SI/cluster/umap_harmony_E.pdf",width = 6,height = 5)
DimPlot(n, reduction = "umap",label = TRUE)
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/cluster/umap_ori_harmony_E.pdf",width = 16,height = 5)
DimPlot(n, reduction = "umap", split.by = "orig.ident")
dev.off()
saveRDS(n, file = "H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_E.rds")

marker27<-FindMarkers(n,ident.1="27",only.pos=TRUE,assay ="originalexp")
write.csv(marker27,file = "H:/WRY/data/newRNA/SI/cluster/marker/E/marker27_all.csv")

#n$cell_type<-meta[rownames(n@meta.data),"seurat_clusters"]
#DimPlot(n, reduction = "umap",label = TRUE,group.by = "cell_type")

n0<-subset(n,subset = orig.ident == "s00h")
pdf(file = "H:/WRY/data/newRNA/SI/cluster/marker/E/marker_Tuft.pdf",width = 10,height = 5)
FeaturePlot(n0, features = c("Sh2d6","Dclk1"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, label.size = 2,pt.size = 0.3)
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/cluster/marker/E/marker_EEC.pdf",width = 10,height = 15)
FeaturePlot(n0, features = c("Chgb","Chga","Reg4","Gip","Cck","Nts"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, label.size = 2,pt.size = 0.3)
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/cluster/marker/E/marker_Goblet.pdf",width = 10,height = 10)
FeaturePlot(n0, features = c("Tff3","Zg16","Fcgbp","Muc2"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, label.size = 2,pt.size = 0.3)
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/cluster/marker/E/marker_Paneth.pdf",width = 10,height = 10)
FeaturePlot(n0, features = c("Lyz1","Defa21","Defa17","Mptx2"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, label.size = 2,pt.size = 0.3)
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/cluster/marker/E/marker_ISC.pdf",width = 10,height = 10)
FeaturePlot(n0, features = c("Olfm4","Lgr5","Ascl2"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, label.size = 2,pt.size = 0.3)
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/cluster/marker/E/marker_TA.pdf",width = 10,height = 10)
FeaturePlot(n0, features = c("Mki67","Pcna","Myc","Cdk6"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, label.size = 2,pt.size = 0.3)
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/cluster/marker/E/marker_Secretory.pdf",width = 10,height = 10)
FeaturePlot(n0, features = c("Hnf4a","Hnf4g","Hes1","Atoh1"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, label.size = 2,pt.size = 0.3)
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/cluster/marker/E/marker_Enterocyte.pdf",width = 20,height = 15)
FeaturePlot(n0, features = c("Krt19","Lct","Cbr1","Ephx2","Fgf15","Krt20","Apoa1","Alpi","Apoa4","Fabp1"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, label.size = 2,pt.size = 0.3)
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/cluster/marker/E/marker_Enterocyte_space.pdf",width = 15,height = 15)
FeaturePlot(n0, features = c("Ada","Reg1","Acad9","Fgd4","Gpx2","Lgals4","Rps12","Rps14","Rps18"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, label.size = 2,pt.size = 0.3)
dev.off()
VlnPlot(n0,features = c("Mki67"))
pdf(file = "H:/WRY/data/newRNA/SI/cluster/marker/E/marker_all.pdf",width = 10,height = 10)
FeaturePlot(n0, features = c("Fabp6"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, label.size = 2,pt.size = 0.3)
dev.off()

n <- readRDS("H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_E.rds")
new.cluster.ids <- c("Early Enterocyte",
                     "Late Enterocyte",
                     "Early Enterocyte",
                     "Late Enterocyte",
                     "TA",
                     "Late Enterocyte",
                     "Late Enterocyte",
                     "Early Enterocyte",
                     "Early Enterocyte",
                     "Late Enterocyte",
                     "ISC",
                     "TA",
                     "Goblet",
                     "TA",
                     "Goblet",
                     "Late Enterocyte",
                     "Goblet",
                     "Late Enterocyte",
                     "Early Enterocyte",
                     "Late Enterocyte",
                     "Early Enterocyte",
                     "Paneth",
                     "Paneth",
                     "Early Enterocyte",
                     "Tuft",
                     "EEC",
                     "EEC",
                     "Early Enterocyte",
                     "Late Enterocyte",
                     "Early Enterocyte")
names(new.cluster.ids) <- levels(n)
n <- RenameIdents(n, new.cluster.ids)
DimPlot(n, reduction = "umap",label = TRUE)
n$cell_type <- Idents(n)

levels <- c("ISC", 
            "TA",
            "Paneth", 
            "Goblet", 
            "EEC",
            "Tuft",
            "Early Enterocyte",
            "Late Enterocyte")
Idents(n) <- factor(Idents(n), levels = levels)
DimPlot(n, reduction = "umap",label = TRUE)
pdf(file = "H:/WRY/data/newRNA/SI/cluster/umap_annotion_E.pdf",width = 7,height = 5)
DimPlot(n, reduction = "umap",label = TRUE, label.size = 4, pt.size = 0.15)
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/cluster/umap_annotion_ori_E.pdf",width = 19,height = 5)
DimPlot(n, reduction = "umap",label = TRUE, label.size = 4, pt.size = 0.15,split.by = "orig.ident")
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/cluster/umap_annotion_ori_E_nolabel.pdf",width = 19,height = 5)
DimPlot(n, reduction = "umap",label = FALSE, pt.size = 0.15,split.by = "orig.ident")
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/cluster/umap_annotion_E_Nolegend.pdf",width = 5,height = 5)
DimPlot(n, reduction = "umap",label = FALSE, label.size = 4, pt.size = 0.15)+NoLegend()
dev.off()
saveRDS(n, file = "H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_cluster_E.rds")


n <- readRDS("H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_cluster_E.rds")
n<- subset(n, subset = cell_type == "Goblet" |
             cell_type == "Paneth" |
             cell_type == "EEC" |
             cell_type == "Tuft" |
             cell_type == "Early Enterocyte" |
             cell_type == "Late Enterocyte")
new <- readRDS("H:/WRY/data/newRNA/SI/new/merge_new.rds")
new <- new[,rownames(n@meta.data)]
new$cell_type<-n@meta.data[rownames(new@meta.data),"cell_type"]
Idents(new)<-new$cell_type
new$sample<-paste(new$cell_type,new$orig.ident,sep="_")
cat<-c()
stage<-c("s00h","s02h","s06h","s72h")
type <- c("Paneth",
          "Goblet",
          "EEC",
          "Tuft",
          "Early Enterocyte",
          "Late Enterocyte")
for(i in 1:length(type)){
  for (j in 1:length(stage)) {
    cat<-c(cat,paste(type[i],stage[j],sep="_"))
  }
}
Idents(new)<-new$sample
Idents(new) <- factor(Idents(new), levels = cat)
pdf(file = "H:/WRY/data/newRNA/SI/cluster/E_dotplot_new.pdf",width =8,height = 7)
DotPlot(new,features = rev(c("Lyz1", "Defa22", "Defa24","Gm15293",
                              "Muc2", "Zg16","Fcgbp","Galnt12","Retnlb","Ccl6", "Ccl9", 
                             "Pyy","Gcg","Cck","Scgn","Zbp1","Znfx1","Oasl1","Ddx60",
                             "Reg3b","Reg3g","Lypd8","H2-D1","H2-K1","Tap1")),cols = c("#ededed","#e20000")) + RotatedAxis() + coord_flip()
dev.off()

n$sample<-paste(n$cell_type,n$orig.ident,sep="_")
Idents(n)<-n$sample
Idents(n) <- factor(Idents(n), levels = cat)
pdf(file = "H:/WRY/data/newRNA/SI/cluster/E_dotplot_total.pdf",width = ,height = 7)
DotPlot(n,features = rev(c("Lyz1", "Defa22", "Defa24","Gm15293",
                           "Muc2", "Zg16","Fcgbp","Galnt12","Retnlb","Ccl6", "Ccl9", 
                           "Pyy","Gcg","Cck","Scgn","Zbp1","Znfx1","Oasl1","Ddx60",
                           "Reg3b","Reg3g","Lypd8","H2-D1","H2-K1","Tap1")),cols = c("#ededed","#e20000")) + RotatedAxis() + coord_flip()
dev.off()

n <- readRDS("H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_cluster_E.rds")
n0<-subset(n,subset = orig.ident == "s00h")
marker<-FindAllMarkers(n0,only.pos=TRUE,assay ="originalexp")
write.csv(marker,file="H:/WRY/data/newRNA/SI/cluster/marker/E/marker_E_cell_type.csv")

n <- readRDS("H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_cluster_E.rds")
meta<-n@meta.data
metaE<-meta[meta$orig.ident == "s00h",]
rm(n)
n <- readRDS("H:/WRY/data/newRNA/SI/cluster/merge_clean.rds")
n<-n[,rownames(metaE)]

n@meta.data<-n@meta.data[rownames(n@meta.data),c("orig.ident","nCount_originalexp","nFeature_originalexp","percent.mt","percent.ribo","percent.dissociation")]
n <- NormalizeData(n,normalization.method = "LogNormalize",scale.factor = 10000,assay="originalexp")
n <- FindVariableFeatures(n)
top20 <- marker %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
n <- ScaleData(n,features=top20$gene)
n$cell_type <- metaE[rownames(n@meta.data), "cell_type"]
Idents(n) <- n$cell_type
levels <- c("ISC", 
            "TA",
            "Paneth", 
            "Goblet", 
            "EEC",
            "Tuft",
            "Early Enterocyte",
            "Late Enterocyte")
Idents(n) <- factor(Idents(n), levels = levels)
data <- subset(n, downsample = 300)
pdf("H:/WRY/data/newRNA/SI/cluster/marker/E/SI0h_Heatmap_total_top20.pdf", width = 12, height = 10)
DoHeatmap(data, features = top20$gene, size = 3, angle = 90, label = TRUE)+
  scale_fill_gradientn(colors = rev(brewer.pal(n = 10, name = "RdBu")))
dev.off()

##########################
n <- readRDS("H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_cluster_E.rds")
ratio<-read.csv(file = "H:/WRY/data/newRNA/SI/ntr/new_ratio_qc.csv",row.names = 1)
n$new_ratio<-ratio[rownames(n@meta.data),"new_ratio"]
library(ggpubr)
cat<-c("s00h","s02h","s06h","s72h")
data<-n@meta.data[,c("orig.ident","cell_type","new_ratio")]

data_plot<-data[data$cell_type == "ISC",]
data_plot$type <- factor(data_plot$orig.ident, levels = cat)
p1<-ggplot(data_plot,aes(type,new_ratio))+
  geom_violin(aes(color=type,fill = type),trim = TRUE,scale = "width")+
  #    geom_jitter(size=0.1,width = 0.1)+
  theme_bw()+NoLegend()+theme_classic()+xlab("")+ylab("")+labs(title = "new_ratio")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  geom_signif(data=data_plot,comparisons = list(c("s00h","s02h"),c("s00h","s06h"),c("s00h","s72h")),
              test = wilcox.test,
              map_signif_level = TRUE,textsize = 3,vjust = 0.2,tip_length = 0,
              y_position = c(max(data_plot$new_ratio)+0.02,max(data_plot$new_ratio)+0.04,max(data_plot$new_ratio)+0.06),
              size=0,color="black")
pdf(file = "H:/WRY/data/newRNA/SI/cluster/Vln_new_ratio_ISC.pdf",width = 5,height=4)
p1
dev.off()

data_plot<-data[data$cell_type == "TA",]
data_plot$type <- factor(data_plot$orig.ident, levels = cat)
p1<-ggplot(data_plot,aes(type,new_ratio))+
  geom_violin(aes(color=type,fill = type),trim = TRUE,scale = "width")+
  #    geom_jitter(size=0.1,width = 0.1)+
  theme_bw()+NoLegend()+theme_classic()+xlab("")+ylab("")+labs(title = "new_ratio")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  geom_signif(data=data_plot,comparisons = list(c("s00h","s02h"),c("s00h","s06h"),c("s00h","s72h")),
              test = wilcox.test,
              map_signif_level = TRUE,textsize = 3,vjust = 0.2,tip_length = 0,
              y_position = c(max(data_plot$new_ratio)+0.02,max(data_plot$new_ratio)+0.04,max(data_plot$new_ratio)+0.06),
              size=0,color="black")
pdf(file = "H:/WRY/data/newRNA/SI/cluster/Vln_new_ratio_TA.pdf",width = 5,height=4)
p1
dev.off()

data_plot<-data[data$cell_type == "Paneth",]
data_plot$type <- factor(data_plot$orig.ident, levels = cat)
p1<-ggplot(data_plot,aes(type,new_ratio))+
  geom_violin(aes(color=type,fill = type),trim = TRUE,scale = "width")+
  #    geom_jitter(size=0.1,width = 0.1)+
  theme_bw()+NoLegend()+theme_classic()+xlab("")+ylab("")+labs(title = "new_ratio")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  geom_signif(data=data_plot,comparisons = list(c("s00h","s02h"),c("s00h","s06h"),c("s00h","s72h")),
              test = wilcox.test,
              map_signif_level = TRUE,textsize = 3,vjust = 0.2,tip_length = 0,
              y_position = c(max(data_plot$new_ratio)+0.02,max(data_plot$new_ratio)+0.04,max(data_plot$new_ratio)+0.06),
              size=0,color="black")
pdf(file = "H:/WRY/data/newRNA/SI/cluster/Vln_new_ratio_Paneth.pdf",width = 5,height=4)
p1
dev.off()

data_plot<-data[data$cell_type == "Goblet",]
data_plot$type <- factor(data_plot$orig.ident, levels = cat)
p1<-ggplot(data_plot,aes(type,new_ratio))+
  geom_violin(aes(color=type,fill = type),trim = TRUE,scale = "width")+
  #    geom_jitter(size=0.1,width = 0.1)+
  theme_bw()+NoLegend()+theme_classic()+xlab("")+ylab("")+labs(title = "new_ratio")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  geom_signif(data=data_plot,comparisons = list(c("s00h","s02h"),c("s00h","s06h"),c("s00h","s72h")),
              test = wilcox.test,
              map_signif_level = TRUE,textsize = 3,vjust = 0.2,tip_length = 0,
              y_position = c(max(data_plot$new_ratio)+0.02,max(data_plot$new_ratio)+0.04,max(data_plot$new_ratio)+0.06),
              size=0,color="black")
pdf(file = "H:/WRY/data/newRNA/SI/cluster/Vln_new_ratio_Goblet.pdf",width = 5,height=4)
p1
dev.off()

data_plot<-data[data$cell_type == "Tuft",]
data_plot$type <- factor(data_plot$orig.ident, levels = cat)
p1<-ggplot(data_plot,aes(type,new_ratio))+
  geom_violin(aes(color=type,fill = type),trim = TRUE,scale = "width")+
  #    geom_jitter(size=0.1,width = 0.1)+
  theme_bw()+NoLegend()+theme_classic()+xlab("")+ylab("")+labs(title = "new_ratio")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  geom_signif(data=data_plot,comparisons = list(c("s00h","s02h"),c("s00h","s06h"),c("s00h","s72h")),
              test = wilcox.test,
              map_signif_level = TRUE,textsize = 3,vjust = 0.2,tip_length = 0,
              y_position = c(max(data_plot$new_ratio)+0.02,max(data_plot$new_ratio)+0.04,max(data_plot$new_ratio)+0.06),
              size=0,color="black")
pdf(file = "H:/WRY/data/newRNA/SI/cluster/Vln_new_ratio_Tuft.pdf",width = 5,height=4)
p1
dev.off()

data_plot<-data[data$cell_type == "EEC",]
data_plot$type <- factor(data_plot$orig.ident, levels = cat)
p1<-ggplot(data_plot,aes(type,new_ratio))+
  geom_violin(aes(color=type,fill = type),trim = TRUE,scale = "width")+
  #    geom_jitter(size=0.1,width = 0.1)+
  theme_bw()+NoLegend()+theme_classic()+xlab("")+ylab("")+labs(title = "new_ratio")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  geom_signif(data=data_plot,comparisons = list(c("s00h","s02h"),c("s00h","s06h"),c("s00h","s72h")),
              test = wilcox.test,
              map_signif_level = TRUE,textsize = 3,vjust = 0.2,tip_length = 0,
              y_position = c(max(data_plot$new_ratio)+0.02,max(data_plot$new_ratio)+0.04,max(data_plot$new_ratio)+0.06),
              size=0,color="black")
pdf(file = "H:/WRY/data/newRNA/SI/cluster/Vln_new_ratio_EEC.pdf",width = 5,height=4)
p1
dev.off()

data_plot<-data[data$cell_type == "Early Enterocyte",]
data_plot$type <- factor(data_plot$orig.ident, levels = cat)
p1<-ggplot(data_plot,aes(type,new_ratio))+
  geom_violin(aes(color=type,fill = type),trim = TRUE,scale = "width")+
  #    geom_jitter(size=0.1,width = 0.1)+
  theme_bw()+NoLegend()+theme_classic()+xlab("")+ylab("")+labs(title = "new_ratio")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  geom_signif(data=data_plot,comparisons = list(c("s00h","s02h"),c("s00h","s06h"),c("s00h","s72h")),
              test = wilcox.test,
              map_signif_level = TRUE,textsize = 3,vjust = 0.2,tip_length = 0,
              y_position = c(max(data_plot$new_ratio)+0.02,max(data_plot$new_ratio)+0.04,max(data_plot$new_ratio)+0.06),
              size=0,color="black")
pdf(file = "H:/WRY/data/newRNA/SI/cluster/Vln_new_ratio_Early_Enterocyte.pdf",width = 5,height=4)
p1
dev.off()

data_plot<-data[data$cell_type == "Late Enterocyte",]
data_plot$type <- factor(data_plot$orig.ident, levels = cat)
p1<-ggplot(data_plot,aes(type,new_ratio))+
  geom_violin(aes(color=type,fill = type),trim = TRUE,scale = "width")+
  #    geom_jitter(size=0.1,width = 0.1)+
  theme_bw()+NoLegend()+theme_classic()+xlab("")+ylab("")+labs(title = "new_ratio")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  geom_signif(data=data_plot,comparisons = list(c("s00h","s02h"),c("s00h","s06h"),c("s00h","s72h")),
              test = wilcox.test,
              map_signif_level = TRUE,textsize = 3,vjust = 0.2,tip_length = 0,
              y_position = c(max(data_plot$new_ratio)+0.02,max(data_plot$new_ratio)+0.04,max(data_plot$new_ratio)+0.06),
              size=0,color="black")
pdf(file = "H:/WRY/data/newRNA/SI/cluster/Vln_new_ratio_Late_Enterocyte.pdf",width = 5,height=4)
p1
dev.off()
##################
n <- readRDS("H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_cluster_E.rds")
info<-n@meta.data
write.csv(info, file="H:/WRY/data/newRNA/SI/cluster/cellinfo_SI_E.csv")
ratio<-read.csv(file = "H:/WRY/data/newRNA/SI/ntr/new_ratio_qc.csv",row.names = 1)
n$new_ratio<-ratio[rownames(n@meta.data),"new_ratio"]

data<-n@meta.data[,c("orig.ident","cell_type","new_ratio")]
data$type<-paste0(data$orig.ident,data$cell_type)
mean<-aggregate(data[,3],by=list(type=data$type),FUN=mean)
new_ratio<-data.frame(mean[c(1:8),c(1:2)],mean[c(9:16),2],mean[c(17:24),2],mean[c(25:32),2])
new_ratio$type<-gsub("s00h","",new_ratio$type)
rownames(new_ratio)<-new_ratio$type
new_ratio<-new_ratio[,c(2:5)]
levels <- c("ISC", 
            "TA",
            "Paneth", 
            "Goblet", 
            "EEC",
            "Tuft",
            "Early Enterocyte",
            "Late Enterocyte")
names(new_ratio)<-c("s00h","s02h","s06h","s72h")
new_ratio<-new_ratio[levels,]
write.csv(new_ratio,file="H:/WRY/data/newRNA/SI/cluster/newratio_mean_E.csv")

color<-c("#313695","#4575B4","#74ADD1","#ABD9E9","#FDAE61","#F46D43","#D73027","#A50026")
dffs <- cbind(n@reductions$umap@cell.embeddings,n@meta.data[,c("orig.ident","new_ratio")])
dffs_0h<-dffs[dffs$orig.ident == "s00h",]
dffs_0h1<-dffs_0h[order(dffs_0h[,"new_ratio"], decreasing = FALSE),]
dffs_0h1[dffs_0h1$new_ratio > 0.4,"new_ratio"] <- 0.4
pdf("H:/WRY/data/newRNA/SI/cluster/new_ratio_E_legend.pdf", width = 5, height = 5)
ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.1)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.4))
dev.off()
pdf("H:/WRY/data/newRNA/SI/cluster/new_ratio_E_0h.pdf", width = 5, height = 5)
ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.1)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.4))+NoLegend()
dev.off()
dffs_0h<-dffs[dffs$orig.ident == "s02h",]
dffs_0h1<-dffs_0h[order(dffs_0h[,"new_ratio"], decreasing = FALSE),]
dffs_0h1[dffs_0h1$new_ratio > 0.4,"new_ratio"] <- 0.4
pdf("H:/WRY/data/newRNA/SI/cluster/new_ratio_E_2h.pdf", width = 5, height = 5)
ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.1)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.4))+NoLegend()
dev.off()
dffs_0h<-dffs[dffs$orig.ident == "s06h",]
dffs_0h1<-dffs_0h[order(dffs_0h[,"new_ratio"], decreasing = FALSE),]
dffs_0h1[dffs_0h1$new_ratio > 0.4,"new_ratio"] <- 0.4
pdf("H:/WRY/data/newRNA/SI/cluster/new_ratio_E_6h.pdf", width = 5, height = 5)
ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.1)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.4))+NoLegend()
dev.off()
dffs_0h<-dffs[dffs$orig.ident == "s72h",]
dffs_0h1<-dffs_0h[order(dffs_0h[,"new_ratio"], decreasing = FALSE),]
dffs_0h1[dffs_0h1$new_ratio > 0.4,"new_ratio"] <- 0.4
pdf("H:/WRY/data/newRNA/SI/cluster/new_ratio_E_72h.pdf", width = 5, height = 5)
ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.1)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.4))+NoLegend()
dev.off()

pdf(file = "H:/WRY/data/newRNA/SI/cluster/Vln_newratio_E.pdf",width = 36,height = 5)
VlnPlot(n, features = c("new_ratio"),split.by = "orig.ident")
dev.off()
pdf("H:/WRY/data/newRNA/SI/cluster/new_ratio_E_nolegend.pdf", width = 5, height = 5)
FeaturePlot(n, features = "new_ratio", reduction = "umap", label = FALSE,order = TRUE,pt.size = 0.3)+NoLegend()+
  scale_color_gradientn(colors = brewer.pal(n = 9, name = "PuRd"))+labs(title = "")
dev.off()
pdf("H:/WRY/data/newRNA/SI/cluster/new_ratio_E_split.pdf", width = 21, height = 5)
FeaturePlot(n, features = "new_ratio", reduction = "umap", cols =c(colors = brewer.pal(n = 9, name = "PuRd")),label = FALSE,order = TRUE,split.by = "orig.ident",pt.size = 0.3)
dev.off()
pdf("H:/WRY/data/newRNA/SI/cluster/new_ratio_E_label.pdf", width = 6, height = 5)
FeaturePlot(n, features = "new_ratio", reduction = "umap", label = TRUE, order = TRUE,pt.size = 0.3)+
  scale_color_gradientn(colors = brewer.pal(n = 9, name = "PuRd"))+labs(title = "")
dev.off()

