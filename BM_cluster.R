library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)
library(harmony)
n00<-readRDS("H:/WRY/data/newRNA/BM/process/BM_clean/bm0h_clean.rds")
n02<-readRDS("H:/WRY/data/newRNA/BM/process/BM_clean/bm2h_clean.rds")
n06<-readRDS("H:/WRY/data/newRNA/BM/process/BM_clean/bm6h_clean.rds")
n12<-readRDS("H:/WRY/data/newRNA/BM/process/BM_clean/bm12h_clean.rds")
n24<-readRDS("H:/WRY/data/newRNA/BM/process/BM_clean/bm24h_clean.rds")
n48<-readRDS("H:/WRY/data/newRNA/BM/process/BM_clean/bm48h_clean.rds")
n72<-readRDS("H:/WRY/data/newRNA/BM/process/BM_clean/bm72h_clean.rds")
n <- merge(n00, c(n02,n06,n12,n24,n48,n72), add.cell.ids = c("n00h", "n02h", "n06h","n12h", "n24h", "n48h","n72h"), project = "bm")
saveRDS(n, file = "H:/WRY/data/newRNA/BM/cluster/merge_clean.rds")

cc<-read.csv("H:/WRY/code/ccgene/mouse_cc.csv",row.names = 1)
s.genes <- cc[cc$stage == "S","symbol"]
g2m.genes <- cc[cc$stage == "G2/M","symbol"]

n <- readRDS("H:/WRY/data/newRNA/BM/cluster/merge_clean.rds")
n@meta.data<-n@meta.data[rownames(n@meta.data),c("orig.ident","nCount_originalexp","nFeature_originalexp","percent.mt","percent.ribo","percent.dissociation")]
VlnPlot(n, features = c("nFeature_originalexp", "nCount_originalexp", "percent.mt"), ncol = 3)

n <- NormalizeData(n,normalization.method = "LogNormalize",scale.factor = 10000,assay="originalexp")
saveRDS(n, file = "H:/WRY/data/newRNA/BM/cluster/merge_clean_Norm.rds")
n <- FindVariableFeatures(n)
n <- CellCycleScoring(n, s.features = s.genes, g2m.features = g2m.genes)
n$CC.Difference <- n$S.Score - n$G2M.Score
n <- ScaleData(n,vars.to.regress = c('nCount_originalexp', 'percent.mt', 'percent.ribo', 'percent.dissociation', 'CC.Difference'))
n <- RunPCA(n)
n <- RunHarmony(n, group.by.vars = "orig.ident",assay.use = "originalexp",plot_convergence = TRUE)
n <- RunUMAP(n, reduction = "harmony", dims = 1:30)
n <- FindNeighbors(n, reduction = "harmony", dims = 1:30) 
n <- FindClusters(n,resolution=2)
pdf(file = "H:/WRY/data/newRNA/BM/cluster/umap_harmony.pdf",width = 6,height = 5)
DimPlot(n, reduction = "umap",label=TRUE)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/cluster/umap_ori_harmony.pdf",width = 36,height = 5)
DimPlot(n, reduction = "umap", split.by = "orig.ident")
dev.off()
saveRDS(n, file = "H:/WRY/data/newRNA/BM/cluster/combind_clean_harmony.rds")

FeaturePlot(n, features = c("Siglech","Itgax"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, pt.size = 0.3)
n0<-subset(n,subset = orig.ident == "n00h")
pdf(file = "H:/WRY/data/newRNA/BM/cluster/umap_harmony_n00h.pdf",width = 6,height = 5)
DimPlot(n0, reduction = "umap",label=TRUE,pt.size = 0.6)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/cluster/marker/cell_type/marker_all.pdf",width = 20,height = 20)
FeaturePlot(n0, features = c("Cd34","Mki67","Itgam","Ly6g","Csf1r","Adgre4","S100a8","Prss34","Gata2","Cd3e","Klrb1c","Ighd","Ighm","Jchain","Hba-a1"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, pt.size = 0.3)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/cluster/marker/cell_type/marker_HSPC.pdf",width = 15,height = 15)
FeaturePlot(n0, features = c("Cd34","H2afy","Mif","Npm1","Kit","Myc","Gata2","Lmo4"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, pt.size = 0.3)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/cluster/marker/cell_type/marker_T_NK.pdf",width = 10,height = 10)
FeaturePlot(n0, features = c("Cd3e","Trbc1","Klrb1c"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, pt.size = 0.3)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/cluster/marker/cell_type/marker_PreB.pdf",width = 10,height = 10)
FeaturePlot(n0, features = c("Vpreb3","Chchd10","Igll1","Ebf1"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, pt.size = 0.3)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/cluster/marker/cell_type/marker_B.pdf",width = 10,height = 15)
FeaturePlot(n0, features = c("Igkc","Igha","Jchain","Sdc1","Cd38"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, pt.size = 0.3)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/cluster/marker/cell_type/marker_B1.pdf",width = 10,height = 10)
FeaturePlot(n0, features = c("Ighm","Ighd","Cd79a","Cd79b"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, pt.size = 0.3)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/cluster/marker/cell_type/marker_Basophil.pdf",width = 10,height = 10)
FeaturePlot(n0, features = c("Prss34","Mcpt8","Cpa3","Fcer1a"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, pt.size = 0.3)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/cluster/marker/cell_type/marker_eosinophil.pdf",width = 10,height = 10)
FeaturePlot(n0, features = c("Gata1","Spi1","Id1"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, pt.size = 0.3)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/cluster/marker/cell_type/marker_Megakaryocyte.pdf",width = 10,height = 10)
FeaturePlot(n0, features = c("Pf4","Plek","Vwf","Gata2"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, pt.size = 0.3)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/cluster/marker/cell_type/marker_Erythroid.pdf",width = 10,height = 10)
FeaturePlot(n0, features = c("Hba-a1","Hbb-bs","Hbb-bt","Car2"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, pt.size = 0.3)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/cluster/marker/cell_type/marker_Macrophage.pdf",width = 10,height = 10)
FeaturePlot(n0, features = c("Ms4a6c","Ctss","Fcna"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, pt.size = 0.3)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/cluster/marker/cell_type/marker_Macrophage1.pdf",width = 10,height = 15)
FeaturePlot(n0, features = c("Adgre1","Adgre4","Apoe","Csf1r","Ly6g"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, pt.size = 0.3,order = TRUE)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/cluster/marker/cell_type/marker_Macrophage2.pdf",width = 36,height = 15)
FeaturePlot(n, features = c("Tnf","Il1b","Il6"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, pt.size = 0.3,split.by = "orig.ident")
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/cluster/marker/cell_type/marker_Neutrophil.pdf",width = 10,height = 15)
FeaturePlot(n0, features = c("Chil3","Elane","Fcnb","Ltf","Mmp8","Mpo"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, pt.size = 0.3)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/cluster/marker/cell_type/marker_Neutrophil1.pdf",width = 10,height = 15)
FeaturePlot(n0, features = c("Mmp9","Ctsg","Camp","Ltf","B2m"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, pt.size = 0.3)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/cluster/marker/cell_type/marker_pDC.pdf",width = 10,height = 10)
FeaturePlot(n0, features = c("Siglech","Cox6a2","Bst2","Ctsl"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, pt.size = 0.3)
dev.off()
VlnPlot(n, features = c("Mki67"))

n <- readRDS("H:/WRY/data/newRNA/BM/cluster/combind_clean_harmony.rds")
new.cluster.ids <- c("Monocyte",
                     "Neutrophil",
                     "Mki67 high Neutrophil",
                     "Monocyte",
                     "Neutrophil",
                     "Mki67 high Neutrophil",
                     "Neutrophil",
                     "Neutrophil",
                     "B",
                     "Mki67 high Monocyte",
                     "Neutrophil",
                     "Mki67 high Neutrophil",
                     "pDC",
                     "Erythroblast",
                     "Pre B",
                     "Neutrophil",
                     "Mki67 high Neutrophil",
                     "HSPC",
                     "T",
                     "Erythroblast",
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
            "Mki67 high Monocyte",
            "Monocyte", 
            "Macrophage", 
            "Mki67 high Neutrophil",
            "Neutrophil",
            "Basophil",
            "NK", 
            "pDC",
            "Pre B",
            "B",
            "T",
            "Megakaryocyte",
            "Erythroblast",
            "Erythroid")
Idents(n) <- factor(Idents(n), levels = levels)
pdf(file = "H:/WRY/data/newRNA/BM/cluster/umap_annotion.pdf",width = 7,height = 5)
DimPlot(n, reduction = "umap",label = TRUE, label.size = 4, pt.size = 0.15)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/cluster/umap_annotion_Nolegend.pdf",width = 4.8,height = 5)
DimPlot(n, reduction = "umap",label = FALSE, label.size = 4, pt.size = 0.15)+NoLegend()
dev.off()
saveRDS(n, file = "H:/WRY/data/newRNA/BM/cluster/combind_clean_harmony_cluster.rds")

n <- readRDS("H:/WRY/data/newRNA/BM/cluster/combind_clean_harmony_cluster.rds")
info<-cbind(n@reductions$umap@cell.embeddings,n@meta.data)
write.csv(info,file="H:/WRY/data/newRNA/BM/cluster/cellinfo_BM_all.csv")

pdf(file = "H:/WRY/data/newRNA/BM/cluster/umap_ori_annotion.pdf",width = 25,height = 5)
DimPlot(n, reduction = "umap", split.by = "orig.ident",label = FALSE,pt.size = 0.1)+NoLegend()
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/cluster/umap_ori_annotion_legend.pdf",width = 36,height = 5)
DimPlot(n, reduction = "umap", split.by = "orig.ident",label = FALSE,pt.size = 0.1)
dev.off()
ratio<-read.csv(file = "H:/WRY/data/newRNA/BM/ntr/new_ratio_qc.csv",row.names = 1)
n$new_ratio<-ratio[rownames(n@meta.data),"new_ratio"]
n$new_count<-ratio[rownames(n@meta.data),"new_count"]
n$total_count<-ratio[rownames(n@meta.data),"total_count"]
VlnPlot(n,features = c("new_count","total_count","new_ratio"),group.by="cell_type")
data<-n@meta.data[,c("orig.ident","cell_type","new_ratio")]
data$type<-paste0(data$orig.ident,data$cell_type)
mean<-aggregate(data[,3],by=list(type=data$type),FUN=mean)
new_ratio<-data.frame(mean[c(1:15),c(1:2)],mean[c(16:30),2],mean[c(31:45),2],mean[c(46:60),2],
                      mean[c(61:75),2],mean[c(76:90),2],mean[c(91:105),2])
new_ratio$type<-gsub("n00h","",new_ratio$type)
rownames(new_ratio)<-new_ratio$type
new_ratio<-new_ratio[,c(2:8)]
levels <- c("HSPC", 
            "Mki67 high Monocyte",
            "Monocyte", 
            "Macrophage", 
            "Mki67 high Neutrophil",
            "Neutrophil",
            "Basophil",
            "NK", 
            "pDC",
            "Pre B",
            "B",
            "T",
            "Megakaryocyte",
            "Erythroblast",
            "Erythroid")
names(new_ratio)<-c("n00h","n02h","n06h","n12h","n24h","n48h","n72h")
new_ratio<-new_ratio[levels,]
write.csv(new_ratio,file="H:/WRY/data/newRNA/BM/cluster/newratio_mean.csv")
# rev(brewer.pal(n = 9, name = "RdYlBu"))
color<-c("#313695","#4575B4","#74ADD1","#ABD9E9","#FDAE61","#F46D43","#D73027","#A50026")
dffs <- cbind(n@reductions$umap@cell.embeddings,n@meta.data[,c("orig.ident","new_ratio")])
dffs_0h<-dffs[dffs$orig.ident == "n00h",]
dffs_0h1<-dffs_0h[order(dffs_0h[,"new_ratio"], decreasing = FALSE),]
dffs_0h1[dffs_0h1$new_ratio > 0.3,"new_ratio"] <- 0.3
pdf("H:/WRY/data/newRNA/BM/cluster/new_ratio_legend.pdf", width = 5, height = 5)
ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.1)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.3))
dev.off()
pdf("H:/WRY/data/newRNA/BM/cluster/new_ratio_0h.pdf", width = 4.8, height = 5)
ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.1)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.3))+NoLegend()
dev.off()
dffs_0h<-dffs[dffs$orig.ident == "n12h",]
dffs_0h1<-dffs_0h[order(dffs_0h[,"new_ratio"], decreasing = FALSE),]
dffs_0h1[dffs_0h1$new_ratio > 0.3,"new_ratio"] <- 0.3
pdf("H:/WRY/data/newRNA/BM/cluster/new_ratio_12h.pdf", width = 4.8, height = 5)
ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.1)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.3))+NoLegend()
dev.off()
dffs_0h<-dffs[dffs$orig.ident == "n72h",]
dffs_0h1<-dffs_0h[order(dffs_0h[,"new_ratio"], decreasing = FALSE),]
dffs_0h1[dffs_0h1$new_ratio > 0.3,"new_ratio"] <- 0.3
pdf("H:/WRY/data/newRNA/BM/cluster/new_ratio_72h.pdf", width = 4.8, height = 5)
ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.1)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.3))+NoLegend()
dev.off()
dffs_0h<-dffs[dffs$orig.ident == "n06h",]
dffs_0h1<-dffs_0h[order(dffs_0h[,"new_ratio"], decreasing = FALSE),]
dffs_0h1[dffs_0h1$new_ratio > 0.3,"new_ratio"] <- 0.3
pdf("H:/WRY/data/newRNA/BM/cluster/new_ratio_6h.pdf", width = 4.8, height = 5)
ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.1)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.3))+NoLegend()
dev.off()
dffs_0h<-dffs[dffs$orig.ident == "n24h",]
dffs_0h1<-dffs_0h[order(dffs_0h[,"new_ratio"], decreasing = FALSE),]
dffs_0h1[dffs_0h1$new_ratio > 0.3,"new_ratio"] <- 0.3
pdf("H:/WRY/data/newRNA/BM/cluster/new_ratio_24h.pdf", width = 4.8, height = 5)
ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.1)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.3))+NoLegend()
dev.off()
dffs_0h<-dffs[dffs$orig.ident == "n48h",]
dffs_0h1<-dffs_0h[order(dffs_0h[,"new_ratio"], decreasing = FALSE),]
dffs_0h1[dffs_0h1$new_ratio > 0.3,"new_ratio"] <- 0.3
pdf("H:/WRY/data/newRNA/BM/cluster/new_ratio_48h.pdf", width = 4.8, height = 5)
ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.1)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.3))+NoLegend()
dev.off()
dffs_0h<-dffs[dffs$orig.ident == "n02h",]
dffs_0h1<-dffs_0h[order(dffs_0h[,"new_ratio"], decreasing = FALSE),]
dffs_0h1[dffs_0h1$new_ratio > 0.3,"new_ratio"] <- 0.3
pdf("H:/WRY/data/newRNA/BM/cluster/new_ratio_02h.pdf", width = 4.8, height = 5)
ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.1)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.3))+NoLegend()
dev.off()


n1<-subset(n, subset = orig.ident =="n00h" | orig.ident =="n12h" |orig.ident =="n72h")
FeaturePlot(n1, features = "new_ratio", reduction = "umap",label = FALSE,order = TRUE,split.by = "orig.ident",pt.size = 0.1)&
  scale_color_gradientn(colors = brewer.pal(n = 9, name = "PuRd"))&theme(legend.position = "right")

pdf("H:/WRY/data/newRNA/BM/cluster/new_ratio_split0.1.pdf", width = 36, height = 5)
FeaturePlot(n, features = "new_ratio", reduction = "umap", cols =c(colors = brewer.pal(n = 9, name = "PuRd")),label = FALSE,order = TRUE,split.by = "orig.ident",pt.size = 0.1)
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/cluster/Vln_newratio_BM.pdf",width = 36,height = 5)
VlnPlot(n, features = c("new_ratio"),split.by = "orig.ident")
dev.off()
pdf("H:/WRY/data/newRNA/BM/cluster/new_ratio_nolegend.pdf", width = 5, height = 5)
FeaturePlot(n, features = "new_ratio", reduction = "umap", label = FALSE,order = TRUE,pt.size = 0.3)+NoLegend()+
  scale_color_gradientn(colors = brewer.pal(n = 9, name = "PuRd"))+labs(title = "")
dev.off()
pdf("H:/WRY/data/newRNA/BM/cluster/new_ratio_label.pdf", width = 6, height = 5)
FeaturePlot(n, features = "new_ratio", reduction = "umap", label = TRUE, order = TRUE,pt.size = 0.3)+
  scale_color_gradientn(colors = brewer.pal(n = 9, name = "PuRd"))+labs(title = "")
dev.off()

data<-subset(n, subset = cell_type == "Neutrophil")
pdf(file = "H:/WRY/data/newRNA/BM/cluster/newratio/Vln_newratio_qc_Neutrophil.pdf",width = 12,height = 5)
VlnPlot(data, features = c("new_ratio"),group.by = "orig.ident")
dev.off()
data<-subset(n, subset = cell_type == "Monocyte")
pdf(file = "H:/WRY/data/newRNA/BM/cluster/newratio/Vln_newratio_qc_Monocyte.pdf",width = 12,height = 5)
VlnPlot(data, features = c("new_ratio"),group.by = "orig.ident")
dev.off()
data<-subset(n, subset = cell_type == "HSPC")
pdf(file = "H:/WRY/data/newRNA/BM/cluster/newratio/Vln_newratio_qc_HSPC.pdf",width = 12,height = 5)
VlnPlot(data, features = c("new_ratio"),group.by = "orig.ident")
dev.off()
data<-subset(n, subset = cell_type == "T")
pdf(file = "H:/WRY/data/newRNA/BM/cluster/newratio/Vln_newratio_qc_T.pdf",width = 12,height = 5)
VlnPlot(data, features = c("new_ratio"),group.by = "orig.ident")
dev.off()
data<-subset(n, subset = cell_type == "B")
pdf(file = "H:/WRY/data/newRNA/BM/cluster/newratio/Vln_newratio_qc_B.pdf",width = 12,height = 5)
VlnPlot(data, features = c("new_ratio"),group.by = "orig.ident")
dev.off()
data<-subset(n, subset = cell_type == "Macrophage")
pdf(file = "H:/WRY/data/newRNA/BM/cluster/newratio/Vln_newratio_qc_Macrophage.pdf",width = 12,height = 5)
VlnPlot(data, features = c("new_ratio"),group.by = "orig.ident")
dev.off()

library(org.Mm.eg.db)
library(ggplot2)
library(ggpubr)
library(cowplot)
info<-n@meta.data[,c("orig.ident","new_ratio","cell_type")]
datap<-info[info$cell_type == "Macrophage",]
mean<-aggregate(datap[,2],by=list(orig.ident=datap$orig.ident),FUN=mean)
names(mean)<-c("orig.ident","new_ratio")
p<-ggplot(mean,aes(orig.ident,new_ratio))
p2<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#070707")+
  geom_boxplot(data=datap,aes(group = orig.ident),size=0.5,color="#070707",width=0.5,outlier.alpha = 0)+
  geom_point(color="#f32f2f")+geom_line(aes(group=1),color="#f32f2f")+
  theme_bw()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  geom_signif(data=datap,comparisons = list(c("n00h","n02h"),c("n00h","n06h"),c("n00h","n12h"),c("n00h","n24h"),
                                            c("n00h","n48h"),c("n00h","n72h")),
              test = wilcox.test,
              map_signif_level = TRUE,textsize = 3,vjust = 0.2,tip_length = 0,
              y_position = c(max(datap$new_ratio)+0.01,max(datap$new_ratio)+0.02,max(datap$new_ratio)+0.03,
                             max(datap$new_ratio)+0.04,max(datap$new_ratio)+0.05,max(datap$new_ratio)+0.06),
              size=0,color="black")
pdf(file = "H:/WRY/data/newRNA/BM/cluster/newratio/boxplot_Macrophage_pvalue.pdf",width = 5,height = 5)
p2
dev.off()

p2<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#070707")+
  geom_boxplot(data=datap,aes(group = orig.ident),size=0.5,color="#070707",width=0.5,outlier.alpha = 0)+
  geom_point(color="#f32f2f")+geom_line(aes(group=1),color="#f32f2f")+
  theme_bw()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
pdf(file = "H:/WRY/data/newRNA/BM/cluster/newratio/boxplot_Macrophage.pdf",width = 5,height = 3)
p2
dev.off()

info<-n@meta.data[,c("orig.ident","new_ratio","cell_type")]
datap<-info[info$cell_type == "Neutrophil",]
mean<-aggregate(datap[,2],by=list(orig.ident=datap$orig.ident),FUN=mean)
names(mean)<-c("orig.ident","new_ratio")
p<-ggplot(mean,aes(orig.ident,new_ratio))
p2<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#070707")+
  geom_boxplot(data=datap,aes(group = orig.ident),size=0.5,color="#070707",width=0.5,outlier.alpha = 0)+
  geom_point(color="#f32f2f")+geom_line(aes(group=1),color="#f32f2f")+
  theme_bw()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  geom_signif(data=datap,comparisons = list(c("n00h","n02h"),c("n00h","n06h"),c("n00h","n12h"),c("n00h","n24h"),
                                            c("n00h","n48h"),c("n00h","n72h")),
              test = wilcox.test,
              map_signif_level = TRUE,textsize = 3,vjust = 0.2,tip_length = 0,
              y_position = c(max(datap$new_ratio)+0.01,max(datap$new_ratio)+0.02,max(datap$new_ratio)+0.03,
                             max(datap$new_ratio)+0.04,max(datap$new_ratio)+0.05,max(datap$new_ratio)+0.06),
              size=0,color="black")
pdf(file = "H:/WRY/data/newRNA/BM/cluster/newratio/boxplot_Neutrophil_pvalue.pdf",width = 5,height = 5)
p2
dev.off()

p2<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#070707")+
  geom_boxplot(data=datap,aes(group = orig.ident),size=0.5,color="#070707",width=0.5,outlier.alpha = 0)+
  geom_point(color="#f32f2f")+geom_line(aes(group=1),color="#f32f2f")+
  theme_bw()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
pdf(file = "H:/WRY/data/newRNA/BM/cluster/newratio/boxplot_Neutrophil.pdf",width = 5,height = 3)
p2
dev.off()

info<-n@meta.data[,c("orig.ident","new_ratio","cell_type")]
datap<-info[info$cell_type == "Monocyte",]
mean<-aggregate(datap[,2],by=list(orig.ident=datap$orig.ident),FUN=mean)
names(mean)<-c("orig.ident","new_ratio")
p<-ggplot(mean,aes(orig.ident,new_ratio))
p2<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#070707")+
  geom_boxplot(data=datap,aes(group = orig.ident),size=0.5,color="#070707",width=0.5,outlier.alpha = 0)+
  geom_point(color="#f32f2f")+geom_line(aes(group=1),color="#f32f2f")+
  theme_bw()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  geom_signif(data=datap,comparisons = list(c("n00h","n02h"),c("n00h","n06h"),c("n00h","n12h"),c("n00h","n24h"),
                                            c("n00h","n48h"),c("n00h","n72h")),
              test = wilcox.test,
              map_signif_level = TRUE,textsize = 3,vjust = 0.2,tip_length = 0,
              y_position = c(max(datap$new_ratio)+0.01,max(datap$new_ratio)+0.02,max(datap$new_ratio)+0.03,
                             max(datap$new_ratio)+0.04,max(datap$new_ratio)+0.05,max(datap$new_ratio)+0.06),
              size=0,color="black")
pdf(file = "H:/WRY/data/newRNA/BM/cluster/newratio/boxplot_Monocyte_pvalue.pdf",width = 5,height = 5)
p2
dev.off()

p2<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#070707")+
  geom_boxplot(data=datap,aes(group = orig.ident),size=0.5,color="#070707",width=0.5,outlier.alpha = 0)+
  geom_point(color="#f32f2f")+geom_line(aes(group=1),color="#f32f2f")+
  theme_bw()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
pdf(file = "H:/WRY/data/newRNA/BM/cluster/newratio/boxplot_Monocyte.pdf",width = 5,height = 3)
p2
dev.off()

info<-n@meta.data[,c("orig.ident","new_ratio","cell_type")]
datap<-info[info$cell_type == "B",]
mean<-aggregate(datap[,2],by=list(orig.ident=datap$orig.ident),FUN=mean)
names(mean)<-c("orig.ident","new_ratio")
p<-ggplot(mean,aes(orig.ident,new_ratio))
p2<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#070707")+
  geom_boxplot(data=datap,aes(group = orig.ident),size=0.5,color="#070707",width=0.5,outlier.alpha = 0)+
  geom_point(color="#f32f2f")+geom_line(aes(group=1),color="#f32f2f")+
  theme_bw()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  geom_signif(data=datap,comparisons = list(c("n00h","n02h"),c("n00h","n06h"),c("n00h","n12h"),c("n00h","n24h"),
                                            c("n00h","n48h"),c("n00h","n72h")),
              test = wilcox.test,
              map_signif_level = TRUE,textsize = 3,vjust = 0.2,tip_length = 0,
              y_position = c(max(datap$new_ratio)+0.01,max(datap$new_ratio)+0.02,max(datap$new_ratio)+0.03,
                             max(datap$new_ratio)+0.04,max(datap$new_ratio)+0.05,max(datap$new_ratio)+0.06),
              size=0,color="black")
pdf(file = "H:/WRY/data/newRNA/BM/cluster/newratio/boxplot_B_pvalue.pdf",width = 5,height = 5)
p2
dev.off()

p2<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#070707")+
  geom_boxplot(data=datap,aes(group = orig.ident),size=0.5,color="#070707",width=0.5,outlier.alpha = 0)+
  geom_point(color="#f32f2f")+geom_line(aes(group=1),color="#f32f2f")+
  theme_bw()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
pdf(file = "H:/WRY/data/newRNA/BM/cluster/newratio/boxplot_B.pdf",width = 5,height = 3)
p2
dev.off()

info<-n@meta.data[,c("orig.ident","new_ratio","cell_type")]
datap<-info[info$cell_type == "T",]
mean<-aggregate(datap[,2],by=list(orig.ident=datap$orig.ident),FUN=mean)
names(mean)<-c("orig.ident","new_ratio")
p<-ggplot(mean,aes(orig.ident,new_ratio))
p2<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#070707")+
  geom_boxplot(data=datap,aes(group = orig.ident),size=0.5,color="#070707",width=0.5,outlier.alpha = 0)+
  geom_point(color="#f32f2f")+geom_line(aes(group=1),color="#f32f2f")+
  theme_bw()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  geom_signif(data=datap,comparisons = list(c("n00h","n02h"),c("n00h","n06h"),c("n00h","n12h"),c("n00h","n24h"),
                                            c("n00h","n48h"),c("n00h","n72h")),
              test = wilcox.test,
              map_signif_level = TRUE,textsize = 3,vjust = 0.2,tip_length = 0,
              y_position = c(max(datap$new_ratio)+0.01,max(datap$new_ratio)+0.02,max(datap$new_ratio)+0.03,
                             max(datap$new_ratio)+0.04,max(datap$new_ratio)+0.05,max(datap$new_ratio)+0.06),
              size=0,color="black")
pdf(file = "H:/WRY/data/newRNA/BM/cluster/newratio/boxplot_T_pvalue.pdf",width = 5,height = 5)
p2
dev.off()

p2<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#070707")+
  geom_boxplot(data=datap,aes(group = orig.ident),size=0.5,color="#070707",width=0.5,outlier.alpha = 0)+
  geom_point(color="#f32f2f")+geom_line(aes(group=1),color="#f32f2f")+
  theme_bw()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
pdf(file = "H:/WRY/data/newRNA/BM/cluster/newratio/boxplot_T.pdf",width = 5,height = 3)
p2
dev.off()

info<-n@meta.data[,c("orig.ident","new_ratio","cell_type")]
datap<-info[info$cell_type == "HSPC",]
mean<-aggregate(datap[,2],by=list(orig.ident=datap$orig.ident),FUN=mean)
names(mean)<-c("orig.ident","new_ratio")
p<-ggplot(mean,aes(orig.ident,new_ratio))
p2<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#070707")+
  geom_boxplot(data=datap,aes(group = orig.ident),size=0.5,color="#070707",width=0.5,outlier.alpha = 0)+
  geom_point(color="#f32f2f")+geom_line(aes(group=1),color="#f32f2f")+
  theme_bw()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  geom_signif(data=datap,comparisons = list(c("n00h","n02h"),c("n00h","n06h"),c("n00h","n12h"),c("n00h","n24h"),
                                            c("n00h","n48h"),c("n00h","n72h")),
              test = wilcox.test,
              map_signif_level = TRUE,textsize = 3,vjust = 0.2,tip_length = 0,
              y_position = c(max(datap$new_ratio)+0.01,max(datap$new_ratio)+0.02,max(datap$new_ratio)+0.03,
                             max(datap$new_ratio)+0.04,max(datap$new_ratio)+0.05,max(datap$new_ratio)+0.06),
              size=0,color="black")
pdf(file = "H:/WRY/data/newRNA/BM/cluster/newratio/boxplot_HSPC_pvalue.pdf",width = 5,height = 5)
p2
dev.off()

p2<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#070707")+
  geom_boxplot(data=datap,aes(group = orig.ident),size=0.5,color="#070707",width=0.5,outlier.alpha = 0)+
  geom_point(color="#f32f2f")+geom_line(aes(group=1),color="#f32f2f")+
  theme_bw()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
pdf(file = "H:/WRY/data/newRNA/BM/cluster/newratio/boxplot_HSPC.pdf",width = 5,height = 3)
p2
dev.off()