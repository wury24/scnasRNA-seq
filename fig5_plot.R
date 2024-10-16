library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)
library(harmony)
n <- readRDS("H:/WRY/data/newRNA/SI/cluster/B/combind_clean_harmony_cluster_B_new.rds")
#cols<-c("#F8766D","#D59100","#99A800","#00BC56","#06A4FF","#C77CFF","#FE6D8C")
cols<- c('#fe817d','#81b8df')
n0<-subset(n, subset = orig.ident == "s00h")
pdf(file = "H:/WRY/data/newRNA/SI/fig5/B_new_cluster_legend.pdf",width = 6,height = 5)
DimPlot(n0, reduction = "umap",label = TRUE,cols = cols,pt.size = 0.5)
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/fig5/B_new_cluster.pdf",width = 5,height = 5)
DimPlot(n0, reduction = "umap",label = FALSE,cols = cols,pt.size = 0.5)+NoLegend()
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/fig5/B_new_cluster_ori.pdf",width = 16,height = 5)
DimPlot(n, reduction = "umap",label = FALSE,cols = cols,pt.size = 0.5,split.by = "orig.ident")
dev.off()

color<-c("#313695","#4575B4","#74ADD1","#ABD9E9","#FDAE61","#F46D43","#D73027","#A50026")
ratio<-read.csv(file = "H:/WRY/data/newRNA/SI/ntr/new_ratio_qc.csv",row.names = 1)
n0$new_ratio<-ratio[rownames(n0@meta.data),"new_ratio"]
dffs <- cbind(n0@reductions$umap@cell.embeddings,n0@meta.data[,c("orig.ident","new_ratio")])
dffs_0h1<-dffs[order(dffs[,"new_ratio"], decreasing = FALSE),]
dffs_0h1[dffs_0h1$new_ratio > 0.4,"new_ratio"] <- 0.4
pdf(file = "H:/WRY/data/newRNA/SI/fig5/B_new_ratio_legend.pdf",width = 6,height = 5)
p1<-ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.5)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.4))
p1
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/fig5/B_new_ratio.pdf",width = 5,height = 5)
p1<-ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = new_ratio))+geom_point(size=0.5)+
  theme_classic()+scale_color_gradientn(colors = color,limits=c(0,0.4))+NoLegend()
p1
dev.off()
##########################
color<-c("#313695","#4575B4","#74ADD1","#ABD9E9","#FDAE61","#F46D43","#D73027","#A50026")
ratio<-read.csv(file = "H:/WRY/data/newRNA/SI/ntr/new_ratio_qc.csv",row.names = 1)
n$new_ratio<-ratio[rownames(n@meta.data),"new_ratio"]
dffs <- cbind(n@reductions$umap@cell.embeddings,n@meta.data[,c("orig.ident","new_ratio")])
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

pdf("H:/WRY/data/newRNA/SI/fig5/new_ratio_ori_B.pdf", width = 21, height = 4.5)
plot_grid(p1,p2,p3,p4,ncol = 4)
dev.off()
######################
mycolor<-cols
cell.prop<-data.frame(prop.table(table(Idents(n), n$orig.ident)))
colnames(cell.prop)<-c("cell_type","stage","proportion")
cell_count<-data.frame(table(Idents(n), n$orig.ident))
write.csv(cell_count, file = "H:/WRY/data/newRNA/SI/fig5/cell_count_B_new.csv")
pdf("H:/WRY/data/newRNA/SI/fig5/prop_orig_B_new.pdf", width = 5, height = 5)
ggplot(cell.prop,aes(stage,proportion,fill=cell_type))+
  geom_bar(stat="identity",position="fill")+scale_fill_manual(values=mycolor)+
  ggtitle("")+theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20))+
  theme_classic()+
  guides(fill=guide_legend(title=NULL))
dev.off()

#####################featureplot
b<-readRDS(file = "H:/WRY/data/newRNA/SI/cluster/B/B_total_data.rds")
n0$Igha_total<-b@assays$originalexp@data["Igha",rownames(n0@meta.data)]
n0$Jchain_total<-b@assays$originalexp@data["Jchain",rownames(n0@meta.data)]
n0$Ighm_total<-b@assays$originalexp@data["Ighm",rownames(n0@meta.data)]
n0$Ighd_total<-b@assays$originalexp@data["Ighd",rownames(n0@meta.data)]

pdf(file = "H:/WRY/data/newRNA/SI/fig5/B_marker_new.pdf",width = 11,height = 10)
FeaturePlot(n0, features = c("Igha","Jchain","Ighm","Ighd"), reduction = "umap", cols = c("grey","red"),label = FALSE, pt.size = 0.5,order = T, max.cutoff = "q90")
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/fig5/B_marker_total.pdf",width = 11,height = 10)
FeaturePlot(n0, features = c("Igha_total","Jchain_total","Ighm_total","Ighd_total"), reduction = "umap", cols = c("grey","red"),label = FALSE, pt.size = 0.5,order = T, max.cutoff = "q90")
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/fig5/B_marker_new_all.pdf",width =21,height = 20)
FeaturePlot(n, features = c("Igha","Jchain","Ighm","Ighd"), reduction = "umap", split.by = "orig.ident", cols = c("grey","red"),label = FALSE, pt.size = 0.5,order = T, max.cutoff = "q90")&theme(legend.position = "right")
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/fig5/B_marker_total_all.pdf",width =21,height = 20)
FeaturePlot(n, features = c("Igha_total","Jchain_total","Ighm_total","Ighd_total"), split.by = "orig.ident", reduction = "umap", cols = c("grey","red"),label = FALSE, pt.size = 0.5,order = T, max.cutoff = "q90")&theme(legend.position = "right")
dev.off()

n$Igha_total<-b@assays$originalexp@data["Igha",rownames(n@meta.data)]
n$Jchain_total<-b@assays$originalexp@data["Jchain",rownames(n@meta.data)]
n$Ighm_total<-b@assays$originalexp@data["Ighm",rownames(n@meta.data)]
n$Ighd_total<-b@assays$originalexp@data["Ighd",rownames(n@meta.data)]
n$Igha_new<-n@assays$RNA@data["Igha",rownames(n@meta.data)]
n$Jchain_new<-n@assays$RNA@data["Jchain",rownames(n@meta.data)]
n$Ighm_new<-n@assays$RNA@data["Ighm",rownames(n@meta.data)]
n$Ighd_new<-n@assays$RNA@data["Ighd",rownames(n@meta.data)]

color<-c("#ededed","#e20000")
gene_plot<-function(gene,min,max){
  dffs <- cbind(n@reductions$umap@cell.embeddings,n@meta.data[,c("orig.ident",gene)])
  dffs_0h<-dffs[dffs$orig.ident == "s00h",]
  dffs_0h1<-dffs_0h[order(dffs_0h[,gene], decreasing = FALSE),]
  names(dffs_0h1)[4]<-"exp"
  dffs_0h1[dffs_0h1$exp > max,"exp"] <- max
  dffs_0h1[dffs_0h1$exp < min,"exp"] <- min
  p1<-ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = exp))+geom_point(size=0.5)+
    theme_classic()+scale_color_gradient(low = "#ededed",high = "#e20000",limits=c(min,max))
  
  dffs_0h<-dffs[dffs$orig.ident == "s02h",]
  dffs_0h1<-dffs_0h[order(dffs_0h[,gene], decreasing = FALSE),]
  names(dffs_0h1)[4]<-"exp"
  dffs_0h1[dffs_0h1$exp > max,"exp"] <- max
  dffs_0h1[dffs_0h1$exp < min,"exp"] <- min
  p2<-ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = exp))+geom_point(size=0.5)+
    theme_classic()+scale_color_gradient(low = "#ededed",high = "#e20000",limits=c(min,max))
  
  dffs_0h<-dffs[dffs$orig.ident == "s06h",]
  dffs_0h1<-dffs_0h[order(dffs_0h[,gene], decreasing = FALSE),]
  names(dffs_0h1)[4]<-"exp"
  dffs_0h1[dffs_0h1$exp > max,"exp"] <- max
  dffs_0h1[dffs_0h1$exp < min,"exp"] <- min
  p3<-ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = exp))+geom_point(size=0.5)+
    theme_classic()+scale_color_gradient(low = "#ededed",high = "#e20000",limits=c(min,max))
  
  dffs_0h<-dffs[dffs$orig.ident == "s72h",]
  dffs_0h1<-dffs_0h[order(dffs_0h[,gene], decreasing = FALSE),]
  names(dffs_0h1)[4]<-"exp"
  dffs_0h1[dffs_0h1$exp > max,"exp"] <- max
  dffs_0h1[dffs_0h1$exp < min,"exp"] <- min
  p4<-ggplot(dffs_0h1,aes(x = UMAP_1, y = UMAP_2, color = exp))+geom_point(size=0.5)+
    theme_classic()+scale_color_gradient(low = "#ededed",high = "#e20000",limits=c(min,max))
  
  p<-plot_grid(p1,p2,p3,p4,ncol = 4)
  return(p)
}

pdf("H:/WRY/data/newRNA/SI/fig5/Igha_total_all.pdf", width = 21, height = 5)
gene_plot(gene="Igha_total",min=0,max=8)
dev.off()
pdf("H:/WRY/data/newRNA/SI/fig5/Jchain_total_all.pdf", width = 21, height = 5)
gene_plot(gene="Jchain_total",min=0,max=7)
dev.off()
pdf("H:/WRY/data/newRNA/SI/fig5/Ighm_total_all.pdf", width = 21, height = 5)
gene_plot(gene="Ighm_total",min=0,max=3)
dev.off()
pdf("H:/WRY/data/newRNA/SI/fig5/Ighd_total_all.pdf", width = 21, height = 5)
gene_plot(gene="Ighd_total",min=0,max=2.5)
dev.off()

pdf("H:/WRY/data/newRNA/SI/fig5/Igha_new_all.pdf", width = 21, height = 5)
gene_plot(gene="Igha_new",min=0,max=6)
dev.off()
pdf("H:/WRY/data/newRNA/SI/fig5/Jchain_new_all.pdf", width = 21, height = 5)
gene_plot(gene="Jchain_new",min=0,max=5)
dev.off()
pdf("H:/WRY/data/newRNA/SI/fig5/Ighm_new_all.pdf", width = 21, height = 5)
gene_plot(gene="Ighm_new",min=0,max=4)
dev.off()
pdf("H:/WRY/data/newRNA/SI/fig5/Ighd_new_all.pdf", width = 21, height = 5)
gene_plot(gene="Ighd_new",min=0,max=2.5)
dev.off()
###################GSEA
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(tidyverse)
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
write.csv(fgsea_cd8[,c(1:7)],file="H:/WRY/data/newRNA/SI/fig5/GSEA_M5_Plasma_72h.csv")

p1<-plotEnrichment(fgsea_sets1[["GOBP_IMMUNOGLOBULIN_PRODUCTION"]],
                   ranks_cd8) + labs(title="GOBP_IMMUNOGLOBULIN_PRODUCTION(Plasma)")
p2<-plotEnrichment(fgsea_sets1[["GOBP_HUMORAL_IMMUNE_RESPONSE"]],
                   ranks_cd8) + labs(title="GOBP_HUMORAL_IMMUNE_RESPONSE(Plasma)")
p3<-plotEnrichment(fgsea_sets1[["GOBP_POSITIVE_REGULATION_OF_B_CELL_ACTIVATION"]],
                   ranks_cd8) + labs(title="GOBP_POSITIVE_REGULATION_OF_B_CELL_ACTIVATION(Plasma)")
pdf(file="H:/WRY/data/newRNA/SI/fig5/GSEA_M5_Plasma_72h_plot.pdf",width = 10,height=6)
plot_grid(p1,p2,p3,ncol = 2)
dev.off()

p1<-plotEnrichment(fgsea_sets1[["GOCC_IMMUNOGLOBULIN_COMPLEX"]],
                   ranks_cd8) + labs(title="GOCC_IMMUNOGLOBULIN_COMPLEX(Plasma)")
p2<-plotEnrichment(fgsea_sets1[["GOBP_HUMORAL_IMMUNE_RESPONSE"]],
                   ranks_cd8) + labs(title="GOBP_HUMORAL_IMMUNE_RESPONSE(Plasma)")
p3<-plotEnrichment(fgsea_sets1[["GOBP_REGULATION_OF_B_CELL_ACTIVATION"]],
                   ranks_cd8) + labs(title="GOBP_REGULATION_OF_B_CELL_ACTIVATION(Plasma)")
pdf(file="H:/WRY/data/newRNA/SI/fig5/GSEA_M5_Plasma_72h_plot1.pdf",width = 10,height=6)
plot_grid(p1,p2,p3,ncol = 2)
dev.off()
#############################
new <- readRDS("H:/WRY/data/newRNA/SI/cluster/B/combind_clean_harmony_cluster_B_new.rds")
n<-readRDS(file = "H:/WRY/data/newRNA/SI/cluster/B/B_total_data.rds")
n$new_type<-new@meta.data[rownames(n@meta.data),"cell_type"]
Idents(new)<-new$cell_type
Idents(n)<-n$new_type
new$sample<-paste(new$cell_type,new$orig.ident,sep="_")
cat1<-c()
stage<-c("s00h","s02h","s06h","s72h")
type <- c("B",
          "Plasma")
for(i in 1:length(type)){
  for (j in 1:length(stage)) {
    cat1<-c(cat1,paste(type[i],stage[j],sep="_"))
  }
}

Idents(new)<-new$sample
Idents(new) <- factor(Idents(new), levels = cat1)
gene<-read.csv(file = "H:/WRY/data/newRNA/SI/cluster/sub/B_dotplot.csv",header = F)
pdf(file = "H:/WRY/data/newRNA/SI/fig5/B_dotplot_new_Nolegend.pdf",width =5,height = 3)
DotPlot(new,features = rev(gene[,1]),cols = c("lightgrey","#e20000"),idents=cat1) + RotatedAxis() + coord_flip()
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/fig5/B_dotplot_new.pdf",width =5,height = 5)
DotPlot(new,features = rev(gene[,1]),cols = c("lightgrey","#e20000"),idents=cat1) + RotatedAxis() + coord_flip()
dev.off()

n$sample<-paste(n$new_type,n$orig.ident,sep="_")
Idents(n)<-n$sample
Idents(n) <- factor(Idents(n), levels = cat1)
pdf(file = "H:/WRY/data/newRNA/SI/fig5/B_dotplot_total_Ighe_Nolegend.pdf",width =5,height = 3)
DotPlot(n,features = rev(gene[,1]),cols = c("lightgrey","#e20000"),idents=cat1) + RotatedAxis() + coord_flip()
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/fig5/B_dotplot_total_Ighe.pdf",width =5,height = 5)
DotPlot(n,features = rev(gene[,1]),cols = c("lightgrey","#e20000"),idents=cat1) + RotatedAxis() + coord_flip()
dev.off()
gene<-data.frame(gene[c(1:7),])
pdf(file = "H:/WRY/data/newRNA/SI/fig5/B_dotplot_total_Nolegend.pdf",width =5,height = 3)
DotPlot(n,features = rev(gene[,1]),cols = c("lightgrey","#e20000"),idents=cat1) + RotatedAxis() + coord_flip()
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/fig5/B_dotplot_total.pdf",width =5,height = 5)
DotPlot(n,features = rev(gene[,1]),cols = c("lightgrey","#e20000"),idents=cat1) + RotatedAxis() + coord_flip()
dev.off()
##############################
new <- readRDS("H:/WRY/data/newRNA/SI/cluster/B/combind_clean_harmony_cluster_B_new.rds")
b<-readRDS(file = "H:/WRY/data/newRNA/SI/cluster/B/B_total_data.rds")
b$new_type<-new@meta.data[rownames(b@meta.data),"cell_type"]
Idents(b)<-b$new_type
new$sample<-paste(new$cell_type,new$orig.ident,sep="_")
cat<-c()
stage<-c("s00h","s02h","s06h","s72h")
type <- c("B",
          "Plasma")
for(i in 1:length(type)){
  for (j in 1:length(stage)) {
    cat<-c(cat,paste(type[i],stage[j],sep="_"))
  }
}
Idents(new)<-new$sample
Idents(new) <- factor(Idents(new), levels = cat)

gene<-read.csv(file = "H:/WRY/data/newRNA/SI/cluster/sub/B_dotplot.csv",header = F)

data_Ig<-t(data.frame(new@assays$RNA@data[gene[c(1:7),1],rownames(new@meta.data)]))

data_Ig1<-data.frame(data_Ig)
data_Ig1[,"sample"]<-new@meta.data[rownames(data_Ig1),"cell_type"]
data_Ig1[,"sample1"]<-new@meta.data[rownames(data_Ig1),"orig.ident"]
data_Ig1$sample1<-factor(data_Ig1$sample1)
hist(data_Ig1[data_Ig1$sample == "B","Igha"])
summary(data_Ig1[data_Ig1$sample == "B","Igha"])

hist(data_Ig1[data_Ig1$sample == "Plasma","Igha"])
summary(data_Ig1[data_Ig1$sample == "Plasma","Igha"])
library(ggpubr)
hist_plot<-function(x,y){
  p1<-gghistogram(data_Ig1[data_Ig1$sample == "B" & data_Ig1$sample1 == y,x],fill = "grey80",xlab = paste("B",x,sep="_"))+yscale("log2")
  p2<-gghistogram(data_Ig1[data_Ig1$sample == "Plasma" & data_Ig1$sample1 == y,x],fill = "grey80",xlab = paste("Plasma",x,sep="_"))+yscale("log2")
  p<-plot_grid(p1,p2,nrow = 2)
  return(p)
}
p_1<-hist_plot("Ighd","s00h")
p_2<-hist_plot("Ighd","s02h")
p_3<-hist_plot("Ighd","s06h")
p_4<-hist_plot("Ighd","s72h")
p<-plot_grid(p_1,p_2,p_3,p_4,nrow = 1)
pdf(file = "H:/WRY/data/newRNA/SI/fig5/hist_Ighd_new.pdf",height=9,width=12)
p
dev.off()

p_1<-hist_plot("Ighm","s00h")
p_2<-hist_plot("Ighm","s02h")
p_3<-hist_plot("Ighm","s06h")
p_4<-hist_plot("Ighm","s72h")
p<-plot_grid(p_1,p_2,p_3,p_4,nrow = 1)
pdf(file = "H:/WRY/data/newRNA/SI/fig5/hist_Ighm_new.pdf",height=9,width=12)
p
dev.off()

p_1<-hist_plot("Ighg1","s00h")
p_2<-hist_plot("Ighg1","s02h")
p_3<-hist_plot("Ighg1","s06h")
p_4<-hist_plot("Ighg1","s72h")
p<-plot_grid(p_1,p_2,p_3,p_4,nrow = 1)
pdf(file = "H:/WRY/data/newRNA/SI/fig5/hist_Ighg1_new.pdf",height=9,width=12)
p
dev.off()

p_1<-hist_plot("Ighg2b","s00h")
p_2<-hist_plot("Ighg2b","s02h")
p_3<-hist_plot("Ighg2b","s06h")
p_4<-hist_plot("Ighg2b","s72h")
p<-plot_grid(p_1,p_2,p_3,p_4,nrow = 1)
pdf(file = "H:/WRY/data/newRNA/SI/fig5/hist_Ighg2b_new.pdf",height=9,width=12)
p
dev.off()

p_1<-hist_plot("Ighg2c","s00h")
p_2<-hist_plot("Ighg2c","s02h")
p_3<-hist_plot("Ighg2c","s06h")
p_4<-hist_plot("Ighg2c","s72h")
p<-plot_grid(p_1,p_2,p_3,p_4,nrow = 1)
pdf(file = "H:/WRY/data/newRNA/SI/fig5/hist_Ighg2c_new.pdf",height=9,width=12)
p
dev.off()

p_1<-hist_plot("Ighg3","s00h")
p_2<-hist_plot("Ighg3","s02h")
p_3<-hist_plot("Ighg3","s06h")
p_4<-hist_plot("Ighg3","s72h")
p<-plot_grid(p_1,p_2,p_3,p_4,nrow = 1)
pdf(file = "H:/WRY/data/newRNA/SI/fig5/hist_Ighg3_new.pdf",height=9,width=12)
p
dev.off()

p_1<-hist_plot("Igha","s00h")
p_2<-hist_plot("Igha","s02h")
p_3<-hist_plot("Igha","s06h")
p_4<-hist_plot("Igha","s72h")
p<-plot_grid(p_1,p_2,p_3,p_4,nrow = 1)
pdf(file = "H:/WRY/data/newRNA/SI/fig5/hist_Igha_new.pdf",height=9,width=12)
p
dev.off()

data_Ig<-t(data.frame(new@assays$RNA@data[gene[c(1:7),1],rownames(new@meta.data)]))
data_mul<-data.frame(data_Ig)
data_mul$Igha<-data_mul$Igha>3
data_mul[,c(2:7)]<-data_mul[,c(2:7)]>1

data_mul$Ig<-apply(data_mul[,c(1:7)],1,sum)

select<-function(x){
  gene<-c("NA")
  for(i in 1: 7){
    if(x[i]>0){
      gene<-paste(gene,names(data_mul)[i],sep = "_")
    }
  }
  return(gene)
}
data_mul$gene<-apply(data_mul[,c(1:7)],1,function(x) select(x))
data_mul[,"sample"]<-new@meta.data[rownames(data_mul),"sample"]
stage <- c("s00h","s02h","s06h","s72h")
type <- c("B",
          "Plasma")
cat<-c()
for(i in 1:length(type)){
  for (j in 1:length(stage)) {
    cat<-c(cat,paste(type[i],stage[j],sep="_"))
  }
}
percent1<-data.frame(table(data_mul[data_mul$sample == "B_s00h","gene"]))
names(percent1)<-c("gene","B_s00h")
percent<-percent1
for(i in 2:length(cat)){
  percent2<-data.frame(table(data_mul[data_mul$sample == cat[i],"gene"]))
  names(percent2)<-c("gene",cat[i])
  percent<-merge(percent,percent2,by="gene",all = TRUE)
}
percent_mul<-percent
percent_mul[is.na(percent_mul)]=0
count<-data.frame(table(data_mul$sample))
for(i in 1:length(cat)){
  percent_mul[,i+1]<-percent_mul[,i+1]*100/count[i,2]
}
percent_mul$gene<-gsub("NA_","",percent_mul$gene)

write.csv(percent_mul,file = "H:/WRY/data/newRNA/SI/fig5/Ig_percent_new.csv")

data<-read.csv(file = "H:/WRY/data/newRNA/SI/fig5/Ig_percent_new_plot.csv")
names(data)[1]<-"gene"
cell<-rep(names(data)[2:9],each = 10)
gene<-rep(data[,1],8)
data_plot<-data.frame(gene,c(data[,2],data[,3],data[,4],data[,5],data[,6],data[,7],
                             data[,8],data[,9]),cell)
names(data_plot)[2]<-"percent"
data_plot$gene<-factor(data_plot$gene,levels = data_plot$gene[1:10])
mycolor<-c("#E31A1C","#6A3D9A","#1F78B4","#33A02C","#FDBF6F","#FF7F00",
           "#F781BF","#B15928","#FFFF99","#999999")
pdf(file = "H:/WRY/data/newRNA/SI/fig5/Ig_percent_new.pdf",width = 6,height = 4)
ggplot(data=data_plot,aes(x="",y=percent,fill=gene))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = mycolor)+
  theme(axis.title =element_blank(),panel.background = element_blank(),
        axis.text =element_blank(),line = element_blank(),
        strip.text = element_text(size=10),strip.background = element_blank())+
  coord_polar(theta="y",start=0)+
  facet_wrap(~cell,nrow=2)
dev.off()

#################
##############################
new <- readRDS("H:/WRY/data/newRNA/SI/cluster/B/combind_clean_harmony_cluster_B_new.rds")
b<-readRDS(file = "H:/WRY/data/newRNA/SI/cluster/B/B_total_data.rds")
b$new_type<-new@meta.data[rownames(b@meta.data),"cell_type"]
Idents(b)<-b$new_type
b$sample<-paste(b$new_type,b$orig.ident,sep="_")
cat<-c()
stage<-c("s00h","s02h","s06h","s72h")
type <- c("B",
          "Plasma")
for(i in 1:length(type)){
  for (j in 1:length(stage)) {
    cat<-c(cat,paste(type[i],stage[j],sep="_"))
  }
}
Idents(b)<-b$sample
Idents(b) <- factor(Idents(b), levels = cat)

gene<-read.csv(file = "H:/WRY/data/newRNA/SI/cluster/sub/B_dotplot.csv",header = F)

data_Ig<-t(data.frame(b@assays$originalexp@data[gene[c(1:7),1],rownames(b@meta.data)]))
data_Ig1<-data.frame(data_Ig)
data_Ig1[,"sample"]<-b@meta.data[rownames(data_Ig1),"new_type"]
data_Ig1[,"sample1"]<-b@meta.data[rownames(data_Ig1),"orig.ident"]
data_Ig1$sample1<-factor(data_Ig1$sample1)
hist(data_Ig1[data_Ig1$sample == "B","Igha"])
summary(data_Ig1[data_Ig1$sample == "B","Igha"])

hist(data_Ig1[data_Ig1$sample == "Plasma","Igha"])
summary(data_Ig1[data_Ig1$sample == "Plasma","Igha"])
library(ggpubr)
hist_plot<-function(x,y){
  p1<-gghistogram(data_Ig1[data_Ig1$sample == "B" & data_Ig1$sample1 == y,x],fill = "grey80",xlab = paste("B",x,sep="_"))+yscale("log2")
  p2<-gghistogram(data_Ig1[data_Ig1$sample == "Plasma" & data_Ig1$sample1 == y,x],fill = "grey80",xlab = paste("Plasma",x,sep="_"))+yscale("log2")
  p<-plot_grid(p1,p2,nrow = 2)
  return(p)
}
p_1<-hist_plot("Ighd","s00h")
p_2<-hist_plot("Ighd","s02h")
p_3<-hist_plot("Ighd","s06h")
p_4<-hist_plot("Ighd","s72h")
p<-plot_grid(p_1,p_2,p_3,p_4,nrow = 1)
pdf(file = "H:/WRY/data/newRNA/SI/fig5/hist_Ighd_total.pdf",height=9,width=12)
p
dev.off()

p_1<-hist_plot("Ighm","s00h")
p_2<-hist_plot("Ighm","s02h")
p_3<-hist_plot("Ighm","s06h")
p_4<-hist_plot("Ighm","s72h")
p<-plot_grid(p_1,p_2,p_3,p_4,nrow = 1)
pdf(file = "H:/WRY/data/newRNA/SI/fig5/hist_Ighm_total.pdf",height=9,width=12)
p
dev.off()

p_1<-hist_plot("Ighg1","s00h")
p_2<-hist_plot("Ighg1","s02h")
p_3<-hist_plot("Ighg1","s06h")
p_4<-hist_plot("Ighg1","s72h")
p<-plot_grid(p_1,p_2,p_3,p_4,nrow = 1)
pdf(file = "H:/WRY/data/newRNA/SI/fig5/hist_Ighg1_total.pdf",height=9,width=12)
p
dev.off()

p_1<-hist_plot("Ighg2b","s00h")
p_2<-hist_plot("Ighg2b","s02h")
p_3<-hist_plot("Ighg2b","s06h")
p_4<-hist_plot("Ighg2b","s72h")
p<-plot_grid(p_1,p_2,p_3,p_4,nrow = 1)
pdf(file = "H:/WRY/data/newRNA/SI/fig5/hist_Ighg2b_total.pdf",height=9,width=12)
p
dev.off()

p_1<-hist_plot("Ighg2c","s00h")
p_2<-hist_plot("Ighg2c","s02h")
p_3<-hist_plot("Ighg2c","s06h")
p_4<-hist_plot("Ighg2c","s72h")
p<-plot_grid(p_1,p_2,p_3,p_4,nrow = 1)
pdf(file = "H:/WRY/data/newRNA/SI/fig5/hist_Ighg2c_total.pdf",height=9,width=12)
p
dev.off()

p_1<-hist_plot("Ighg3","s00h")
p_2<-hist_plot("Ighg3","s02h")
p_3<-hist_plot("Ighg3","s06h")
p_4<-hist_plot("Ighg3","s72h")
p<-plot_grid(p_1,p_2,p_3,p_4,nrow = 1)
pdf(file = "H:/WRY/data/newRNA/SI/fig5/hist_Ighg3_total.pdf",height=9,width=12)
p
dev.off()

p_1<-hist_plot("Igha","s00h")
p_2<-hist_plot("Igha","s02h")
p_3<-hist_plot("Igha","s06h")
p_4<-hist_plot("Igha","s72h")
p<-plot_grid(p_1,p_2,p_3,p_4,nrow = 1)
pdf(file = "H:/WRY/data/newRNA/SI/fig5/hist_Igha_total.pdf",height=9,width=12)
p
dev.off()

data_Ig<-t(data.frame(b@assays$originalexp@data[gene[c(1:7),1],rownames(b@meta.data)]))
data_mul<-data.frame(data_Ig)
data_mul$Igha<-data_mul$Igha>5
data_mul$Ighd<-data_mul$Ighd>1
data_mul$Ighm<-data_mul$Ighm>1
data_mul[,c(4:7)]<-data_mul[,c(4:7)]>2

data_mul$Ig<-apply(data_mul[,c(1:7)],1,sum)

select<-function(x){
  gene<-c("NA")
  for(i in 1: 7){
    if(x[i]>0){
      gene<-paste(gene,names(data_mul)[i],sep = "_")
    }
  }
  return(gene)
}
data_mul$gene<-apply(data_mul[,c(1:7)],1,function(x) select(x))
data_mul[,"sample"]<-b@meta.data[rownames(data_mul),"sample"]
stage <- c("s00h","s02h","s06h","s72h")
type <- c("B",
          "Plasma")
cat<-c()
for(i in 1:length(type)){
  for (j in 1:length(stage)) {
    cat<-c(cat,paste(type[i],stage[j],sep="_"))
  }
}
percent1<-data.frame(table(data_mul[data_mul$sample == "B_s00h","gene"]))
names(percent1)<-c("gene","B_s00h")
percent<-percent1
for(i in 2:length(cat)){
  percent2<-data.frame(table(data_mul[data_mul$sample == cat[i],"gene"]))
  names(percent2)<-c("gene",cat[i])
  percent<-merge(percent,percent2,by="gene",all = TRUE)
}
percent_mul<-percent
percent_mul[is.na(percent_mul)]=0
count<-data.frame(table(data_mul$sample))
for(i in 1:length(cat)){
  percent_mul[,i+1]<-percent_mul[,i+1]*100/count[i,2]
}
percent_mul$gene<-gsub("NA_","",percent_mul$gene)

write.csv(percent_mul,file = "H:/WRY/data/newRNA/SI/fig5/Ig_percent_total.csv")

data<-read.csv(file = "H:/WRY/data/newRNA/SI/fig5/Ig_percent_total_plot.csv")
names(data)[1]<-"gene"
cell<-rep(names(data)[2:9],each = 10)
gene<-rep(data[,1],8)
data_plot<-data.frame(gene,c(data[,2],data[,3],data[,4],data[,5],data[,6],data[,7],
                             data[,8],data[,9]),cell)
names(data_plot)[2]<-"percent"
data_plot$gene<-factor(data_plot$gene,levels = data_plot$gene[1:10])
mycolor<-c("#E31A1C","#6A3D9A","#1F78B4","#33A02C","#FDBF6F","#FF7F00",
           "#F781BF","#B15928","#FFFF99","#999999")
pdf(file = "H:/WRY/data/newRNA/SI/fig5/Ig_percent_total.pdf",width = 6,height = 4)
ggplot(data=data_plot,aes(x="",y=percent,fill=gene))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = mycolor)+
  theme(axis.title =element_blank(),panel.background = element_blank(),
        axis.text =element_blank(),line = element_blank(),
        strip.text = element_text(size=10),strip.background = element_blank())+
  coord_polar(theta="y",start=0)+
  facet_wrap(~cell,nrow=2)
dev.off()