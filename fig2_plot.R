library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)
library(harmony)
n <- readRDS("H:/WRY/data/newRNA/BM/cluster/combind_clean_harmony_cluster.rds")
####umap
pdf(file = "H:/WRY/data/newRNA/BM/fig2/BM_umap_ori_annotion.pdf",width = 36,height = 5)
DimPlot(n, reduction = "umap", split.by = "orig.ident",label = FALSE,pt.size = 0.1)
dev.off()
###new ratio
ratio<-read.csv(file = "H:/WRY/data/newRNA/BM/ntr/new_ratio_qc.csv",row.names = 1)
n$new_ratio<-ratio[rownames(n@meta.data),"new_ratio"]
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
write.csv(new_ratio,file="H:/WRY/data/newRNA/BM/fig2/newratio_mean.csv")
###
color<-c("#313695","#4575B4","#74ADD1","#ABD9E9","#FDAE61","#F46D43","#D73027","#A50026")
data <- cbind(n@reductions$umap@cell.embeddings,n@meta.data[,c("orig.ident","new_ratio")])

data_00h<-data[data$orig.ident == "n00h",]
data_00h_order<-data_00h[order(data_00h[,"new_ratio"], decreasing = FALSE),]
data_00h_order[data_00h_order$new_ratio > 0.3,"new_ratio"] <- 0.3
pdf("H:/WRY/data/newRNA/BM/fig2/BM_new_ratio_legend.pdf", width = 5.5, height = 5)
ggplot(data_00h_order,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.1)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.3))
dev.off()
pdf("H:/WRY/data/newRNA/BM/fig2/BM_new_ratio_0h.pdf", width = 6, height = 5)
ggplot(data_00h_order,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.1)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.3))+NoLegend()
dev.off()

data_06h<-data[data$orig.ident == "n06h",]
data_06h_order<-data_06h[order(data_06h[,"new_ratio"], decreasing = FALSE),]
data_06h_order[data_06h_order$new_ratio > 0.3,"new_ratio"] <- 0.3
pdf("H:/WRY/data/newRNA/BM/fig2/BM_new_ratio_6h.pdf", width = 6, height = 5)
ggplot(data_06h_order,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.1)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.3))+NoLegend()
dev.off()

data_72h<-data[data$orig.ident == "n72h",]
data_72h_order<-data_72h[order(data_72h[,"new_ratio"], decreasing = FALSE),]
data_72h_order[data_72h_order$new_ratio > 0.3,"new_ratio"] <- 0.3
pdf("H:/WRY/data/newRNA/BM/fig2/BM_new_ratio_72h.pdf", width = 6, height = 5)
ggplot(data_72h_order,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.1)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.3))+NoLegend()
dev.off()

#############################Heatmap
library(ComplexHeatmap)
library(circlize)
input<-read.csv(file="H:/WRY/data/newRNA/BM/fig2/gene_order.csv",header = FALSE)
new<-read.csv("H:/WRY/data/newRNA/BM/cluster/avg/avgnew.csv",row.names = 1)
stage<-c("n00h","n02h","n06h","n12h","n24h","n48h","n72h")
type   <- c("Monocyte", 
            "Macrophage", 
            "Neutrophil",
            "NK", 
            "pDC",
            "T",
            "B" )
cat<-c()
for(i in 1:length(type)){
  for (j in 1:length(stage)) {
    cat<-c(cat,paste(type[i],stage[j],sep=""))
  }
}
data <- new[input[,1],cat]
sc <- data.frame(t(scale(t(data))))

heatmap_show <- read.csv("H:/WRY/data/newRNA/BM/fig2/gene_select1.csv", header=FALSE)
cat<-factor(cat,levels = cat)
type<-factor(type,levels = type)
datah<-sc[,cat]
datah<-as.matrix(datah)
annotation_col <- data.frame(stage = factor(rep(stage, 7)), cell_type = factor(rep(type,each = 7)))
gene_select<-rownames(datah)[rownames(datah) %in% heatmap_show[,1]]
gene_pos<-which(rownames(datah) %in% gene_select)
row_anno<-rowAnnotation(gene=anno_mark(at=gene_pos,labels = gene_select,labels_gp = gpar(fontsize=15),side="left"))
col_fun <- colorRamp2(c(-1.5,0,3),c("#0000FF", "white", "#DD0000"))

pdf("H:/WRY/data/newRNA/BM/fig2/BM_heatmap_new_BM.pdf", width = 10, height = 10)
Heatmap(datah,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = annotation_col$cell_type,
        col = col_fun,
        column_title = NULL,
        left_annotation = row_anno)
dev.off()

###################new featureplot
n <- readRDS("H:/WRY/data/newRNA/BM/cluster/combind_clean_harmony_cluster.rds")
new <- readRDS("H:/WRY/data/newRNA/BM/new/merge_new.rds")
new <- new[,rownames(n@meta.data)]
new$cell_type<-n@meta.data[rownames(new@meta.data),"cell_type"]

n$Tnf_new<-new@assays$RNA@data["Tnf",rownames(n@meta.data)]
n$Il1b_new<-new@assays$RNA@data["Il1b",rownames(n@meta.data)]
n$Ccl4_new<-new@assays$RNA@data["Ccl4",rownames(n@meta.data)]
n$Ccr1_new<-new@assays$RNA@data["Ccr1",rownames(n@meta.data)]
n$Saa3_new<-new@assays$RNA@data["Saa3",rownames(n@meta.data)]
n$Camp_new<-new@assays$RNA@data["Camp",rownames(n@meta.data)]

n12<-subset(n, subset = orig.ident == "n12h")
pdf(file = "H:/WRY/data/newRNA/BM/fig2/BM_Tnf_new_12h.pdf",width = 5.5,height = 5)
FeaturePlot(n12, features = c("Tnf_new"), reduction = "umap", label = FALSE, pt.size = 0.15,order = T, max.cutoff = "q90")+
  scale_color_gradientn(colors = rev(brewer.pal(n = 9, name = "RdYlBu")))
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/fig2/BM_Il1b_new_12h.pdf",width = 5.5,height = 5)
FeaturePlot(n12, features = c("Il1b_new"), reduction = "umap", label = FALSE, pt.size = 0.15,order = T, max.cutoff = "q90")+
  scale_color_gradientn(colors = rev(brewer.pal(n = 9, name = "RdYlBu")))
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/fig2/BM_Ccl4_new_12h.pdf",width = 5.5,height = 5)
FeaturePlot(n12, features = c("Ccl4_new"), reduction = "umap", label = FALSE, pt.size = 0.15,order = T, max.cutoff = "q90")+
  scale_color_gradientn(colors = rev(brewer.pal(n = 9, name = "RdYlBu")))
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/fig2/BM_Ccr1_new_12h.pdf",width = 5.5,height = 5)
FeaturePlot(n12, features = c("Ccr1_new"), reduction = "umap", label = FALSE, pt.size = 0.15,order = T, max.cutoff = "q90")+
  scale_color_gradientn(colors = rev(brewer.pal(n = 9, name = "RdYlBu")))
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/fig2/BM_Camp_new_12h.pdf",width = 5.5,height = 5)
FeaturePlot(n12, features = c("Camp_new"), reduction = "umap", label = FALSE, pt.size = 0.15,order = T, max.cutoff = "q90")+
  scale_color_gradientn(colors = rev(brewer.pal(n = 9, name = "RdYlBu")))
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/fig2/BM_Saa3_new_12h.pdf",width = 5,height = 5)
FeaturePlot(n12, features = c("Saa3_new"), reduction = "umap", label = FALSE, pt.size = 0.15,order = T, max.cutoff = "q90")+
  scale_color_gradientn(colors = rev(brewer.pal(n = 9, name = "RdYlBu")))
dev.off()

##########GO score
library(stringr)
library(org.Mm.eg.db)
library(ggplot2)
library(ggpubr)
library(cowplot)
GO_term <- read.csv(file="H:/WRY/data/newRNA/BM/fig2/GO.csv",header = FALSE)
names(GO_term)<-c("ID","Description")
n <- readRDS("H:/WRY/data/newRNA/BM/cluster/combind_clean_harmony_cluster.rds")
new <- readRDS("H:/WRY/data/newRNA/BM/new/merge_new.rds")
new <- new[,rownames(n@meta.data)]
new$cell_type<-n@meta.data[rownames(new@meta.data),"cell_type"]
datan<-data.frame(new$cell_type,new$orig.ident)
for (i in 1:length(rownames(GO_term))) {
  DNA_geneID <- get(GO_term[i,"ID"], org.Mm.egGO2ALLEGS)
  DNA_geneSYMBOL <- mget(DNA_geneID, org.Mm.egSYMBOL) %>% unlist() 
  GO<-unique(DNA_geneSYMBOL)
  GO<-intersect(GO,rownames(new@assays$RNA@data))
  GO<-as.list(data.frame(GO))
  new<-AddModuleScore(new,features = GO,name = paste("GO_term",i,sep=""))
  id<-paste(paste("GO_term",i,sep=""),"1",sep="")
  datan<-cbind(datan,new@meta.data[,id])
}
names(datan)[1]<-"cell_type"
names(datan)[2]<-"stage"
names(datan)[c(3:12)]<-GO_term$Description

datat<-data.frame(n$cell_type,n$orig.ident)
for (i in 1:length(rownames(GO_term))) {
  DNA_geneID <- get(GO_term[i,"ID"], org.Mm.egGO2ALLEGS)
  DNA_geneSYMBOL <- mget(DNA_geneID, org.Mm.egSYMBOL) %>% unlist() 
  GO<-unique(DNA_geneSYMBOL)
  GO<-intersect(GO,rownames(n@assays$originalexp@data))
  GO<-as.list(data.frame(GO))
  n<-AddModuleScore(n,features = GO,name = paste("GO_term",i,sep=""))
  id<-paste(paste("GO_term",i,sep=""),"1",sep="")
  datat<-cbind(datat,n@meta.data[,id])
}
names(datat)[1]<-"cell_type"
names(datat)[2]<-"stage"
names(datat)[c(3:12)]<-GO_term$Description

datat$type<-paste(datat$cell_type,datat$stage,sep="_")
datan$type<-paste(datan$cell_type,datan$stage,sep="_")

#####################Inflammatory response Macrophage
k=1
datap<-datat[,c(1:2,k+2)]
names(datap)<-c("cell_type","stage","score")
datap$stage<-factor(datap$stage)
datap<-datap[datap$cell_type == "Macrophage",]
mean<-aggregate(datap[,3],by=list(stage=datap$stage),FUN=mean)
names(mean)<-c("stage","score")
p<-ggplot(mean,aes(stage,score))
p1<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#070707")+
  geom_boxplot(data=datap,aes(group = stage),size=0.5,color="#070707",width=0.5,outlier.alpha = 0)+
  geom_point(color="#070707")+geom_line(aes(group=1),color="#070707")+
  theme_bw()+NoLegend()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
p3<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#070707")+
  geom_boxplot(data=datap,aes(group = stage),size=0.5,color="#070707",width=0.5,outlier.alpha = 0)+
  geom_point(color="#070707")+geom_line(aes(group=1),color="#070707")+
  theme_bw()+NoLegend()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  geom_signif(data=datap,comparisons = list(c("n00h","n02h"),c("n00h","n06h"),c("n00h","n12h"),c("n00h","n24h"),
                                            c("n00h","n48h"),c("n00h","n72h")),
              test = wilcox.test,
              map_signif_level = TRUE,textsize = 3,vjust = 0.2,tip_length = 0,
              y_position = c(max(datap$score)+0.01,max(datap$score)+0.02,max(datap$score)+0.03,
                             max(datap$score)+0.04,max(datap$score)+0.05,max(datap$score)+0.06),
              size=0,color="black")

datap<-datan[,c(1:2,k+2)]
names(datap)<-c("cell_type","stage","score")
datap$stage<-factor(datap$stage)
datap<-datap[datap$cell_type == "Macrophage",]
mean<-aggregate(datap[,3],by=list(stage=datap$stage),FUN=mean)
names(mean)<-c("stage","score")
p<-ggplot(mean,aes(stage,score))
p2<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#f32f2f")+
  geom_boxplot(data=datap,aes(group = stage),size=0.5,color="#f32f2f",width=0.5,outlier.alpha = 0)+
  geom_point(color="#f32f2f")+geom_line(aes(group=1),color="#f32f2f")+
  theme_bw()+NoLegend()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
p4<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#f32f2f")+
  geom_boxplot(data=datap,aes(group = stage),size=0.5,color="#f32f2f",width=0.5,outlier.alpha = 0)+
  geom_point(color="#f32f2f")+geom_line(aes(group=1),color="#f32f2f")+
  theme_bw()+NoLegend()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  geom_signif(data=datap,comparisons = list(c("n00h","n02h"),c("n00h","n06h"),c("n00h","n12h"),c("n00h","n24h"),
                                            c("n00h","n48h"),c("n00h","n72h")),
              test = wilcox.test,
              map_signif_level = TRUE,textsize = 3,vjust = 0.2,tip_length = 0,
              y_position = c(max(datap$score)+0.01,max(datap$score)+0.02,max(datap$score)+0.03,
                             max(datap$score)+0.04,max(datap$score)+0.05,max(datap$score)+0.06),
              size=0,color="black")
pdf(file = "H:/WRY/data/newRNA/BM/fig2/boxplot_inflammatory_response_Macrophage.pdf",width = 7,height = 3)
p1+p2
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/fig2/boxplot_inflammatory_response_Macrophage_pvalue.pdf",width = 7,height = 3)
p3+p4
dev.off()

###########phagocytosis_Neutrophil
k=2
datap<-datat[,c(1:2,k+2)]
names(datap)<-c("cell_type","stage","score")
datap$stage<-factor(datap$stage)
datap<-datap[datap$cell_type == "Neutrophil",]
mean<-aggregate(datap[,3],by=list(stage=datap$stage),FUN=mean)
names(mean)<-c("stage","score")
p<-ggplot(mean,aes(stage,score))
p1<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#070707")+
  geom_boxplot(data=datap,aes(group = stage),size=0.5,color="#070707",width=0.5,outlier.alpha = 0)+
  geom_point(color="#070707")+geom_line(aes(group=1),color="#070707")+
  theme_bw()+NoLegend()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
p3<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#070707")+
  geom_boxplot(data=datap,aes(group = stage),size=0.5,color="#070707",width=0.5,outlier.alpha = 0)+
  geom_point(color="#070707")+geom_line(aes(group=1),color="#070707")+
  theme_bw()+NoLegend()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  geom_signif(data=datap,comparisons = list(c("n00h","n02h"),c("n00h","n06h"),c("n00h","n12h"),c("n00h","n24h"),
                                            c("n00h","n48h"),c("n00h","n72h")),
              test = wilcox.test,
              map_signif_level = TRUE,textsize = 3,vjust = 0.2,tip_length = 0,
              y_position = c(max(datap$score)+0.01,max(datap$score)+0.02,max(datap$score)+0.03,
                             max(datap$score)+0.04,max(datap$score)+0.05,max(datap$score)+0.06),
              size=0,color="black")

datap<-datan[,c(1:2,k+2)]
names(datap)<-c("cell_type","stage","score")
datap$stage<-factor(datap$stage)
datap<-datap[datap$cell_type == "Neutrophil",]
mean<-aggregate(datap[,3],by=list(stage=datap$stage),FUN=mean)
names(mean)<-c("stage","score")
p<-ggplot(mean,aes(stage,score))
p2<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#f32f2f")+
  geom_boxplot(data=datap,aes(group = stage),size=0.5,color="#f32f2f",width=0.5,outlier.alpha = 0)+
  geom_point(color="#f32f2f")+geom_line(aes(group=1),color="#f32f2f")+
  theme_bw()+NoLegend()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
p4<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#f32f2f")+
  geom_boxplot(data=datap,aes(group = stage),size=0.5,color="#f32f2f",width=0.5,outlier.alpha = 0)+
  geom_point(color="#f32f2f")+geom_line(aes(group=1),color="#f32f2f")+
  theme_bw()+NoLegend()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  geom_signif(data=datap,comparisons = list(c("n00h","n02h"),c("n00h","n06h"),c("n00h","n12h"),c("n00h","n24h"),
                                            c("n00h","n48h"),c("n00h","n72h")),
              test = wilcox.test,
              map_signif_level = TRUE,textsize = 3,vjust = 0.2,tip_length = 0,
              y_position = c(max(datap$score)+0.01,max(datap$score)+0.02,max(datap$score)+0.03,
                             max(datap$score)+0.04,max(datap$score)+0.05,max(datap$score)+0.06),
              size=0,color="black")

pdf(file = "H:/WRY/data/newRNA/BM/fig2/boxplot_phagocytosis_Neutrophil.pdf",width = 7,height = 3)
p1+p2
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/fig2/boxplot_phagocytosis_Neutrophil_pvalue.pdf",width = 7,height = 3)
p3+p4
dev.off()

##################################antigen processing and presentation Macrophage
k=7
datap<-datat[,c(1:2,k+2)]
names(datap)<-c("cell_type","stage","score")
datap$stage<-factor(datap$stage)
datap<-datap[datap$cell_type == "Macrophage",]
mean<-aggregate(datap[,3],by=list(stage=datap$stage),FUN=mean)
names(mean)<-c("stage","score")
p<-ggplot(mean,aes(stage,score))
p1<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#070707")+
  geom_boxplot(data=datap,aes(group = stage),size=0.5,color="#070707",width=0.5,outlier.alpha = 0)+
  geom_point(color="#070707")+geom_line(aes(group=1),color="#070707")+
  theme_bw()+NoLegend()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
p3<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#070707")+
  geom_boxplot(data=datap,aes(group = stage),size=0.5,color="#070707",width=0.5,outlier.alpha = 0)+
  geom_point(color="#070707")+geom_line(aes(group=1),color="#070707")+
  theme_bw()+NoLegend()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  geom_signif(data=datap,comparisons = list(c("n00h","n02h"),c("n00h","n06h"),c("n00h","n12h"),c("n00h","n24h"),
                                            c("n00h","n48h"),c("n00h","n72h")),
              test = wilcox.test,
              map_signif_level = TRUE,textsize = 3,vjust = 0.2,tip_length = 0,
              y_position = c(max(datap$score)+0.01,max(datap$score)+0.02,max(datap$score)+0.03,
                             max(datap$score)+0.04,max(datap$score)+0.05,max(datap$score)+0.06),
              size=0,color="black")

datap<-datan[,c(1:2,k+2)]
names(datap)<-c("cell_type","stage","score")
datap$stage<-factor(datap$stage)
datap<-datap[datap$cell_type == "Macrophage",]
mean<-aggregate(datap[,3],by=list(stage=datap$stage),FUN=mean)
names(mean)<-c("stage","score")
p<-ggplot(mean,aes(stage,score))
p2<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#f32f2f")+
  geom_boxplot(data=datap,aes(group = stage),size=0.5,color="#f32f2f",width=0.5,outlier.alpha = 0)+
  geom_point(color="#f32f2f")+geom_line(aes(group=1),color="#f32f2f")+
  theme_bw()+NoLegend()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))
p4<-p+stat_boxplot(data=datap,geom = 'errorbar',width=0.2,cex=1,color="#f32f2f")+
  geom_boxplot(data=datap,aes(group = stage),size=0.5,color="#f32f2f",width=0.5,outlier.alpha = 0)+
  geom_point(color="#f32f2f")+geom_line(aes(group=1),color="#f32f2f")+
  theme_bw()+NoLegend()+theme_classic()+xlab("")+ylab("")+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(axis.line = element_line(size=0.6, colour = 'black'))+
  theme(axis.ticks  = element_line(size=0.6, colour = 'black'))+
  geom_signif(data=datap,comparisons = list(c("n00h","n02h"),c("n00h","n06h"),c("n00h","n12h"),c("n00h","n24h"),
                                            c("n00h","n48h"),c("n00h","n72h")),
              test = wilcox.test,
              map_signif_level = TRUE,textsize = 3,vjust = 0.2,tip_length = 0,
              y_position = c(max(datap$score)+0.01,max(datap$score)+0.02,max(datap$score)+0.03,
                             max(datap$score)+0.04,max(datap$score)+0.05,max(datap$score)+0.06),
              size=0,color="black")
pdf(file = "H:/WRY/data/newRNA/BM/fig2/boxplot_antigen_processing_and_presentation_Macrophage.pdf",width = 7,height = 3)
p1+p2
dev.off()
pdf(file = "H:/WRY/data/newRNA/BM/fig2/boxplot_antigen_processing_and_presentation_Macrophage_pvalue.pdf",width = 7,height = 3)
p3+p4
dev.off()


