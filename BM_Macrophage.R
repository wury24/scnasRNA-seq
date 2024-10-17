library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)
n <- readRDS("H:/WRY/data/newRNA/BM/cluster/combind_clean_harmony_cluster.rds")
levels <- c("Macrophage")
n<-subset(n,subset = cell_type == "Macrophage")
stage <- c("n00h","n02h","n06h","n12h","n24h","n48h","n72h")
main_dir <- "H:/WRY/data/newRNA/BM/fig3/stage"
count_total<-c()
for(j in 2:length(stage)){
    marker<-FindMarkers(n,ident.1 = stage[j],ident.2 = stage[1],min.pct = 0.25,group.by = "orig.ident",only.pos = TRUE)
    write.csv(marker,file =  paste("H:/WRY/data/newRNA/BM/fig3/stage/","Macrophage_",stage[j],"_",stage[1],"_total.csv",sep=""))
    marker<-marker[marker$p_val<0.05,]
    count_total<-c(count_total,length(rownames(marker)))
}

info<-n@meta.data
rm(n)
n <- readRDS("H:/WRY/data/newRNA/BM/new/merge_new.rds")
n<-n[,rownames(info)]
n$cell_type<-info[rownames(n@meta.data),"cell_type"]
Idents(n)<-n$cell_type
main_dir <- "H:/WRY/data/newRNA/BM/fig3/stage"
count_new<-c()
for(j in 2:length(stage)){
    marker<-FindMarkers(n,ident.1 = stage[j],ident.2 = stage[1],group.by = "orig.ident",only.pos = TRUE)
    write.csv(marker,file =  paste("H:/WRY/data/newRNA/BM/fig3/stage/","Macrophage_",stage[j],"_",stage[1],"_new.csv",sep=""))
    marker<-marker[marker$p_val<0.05,]
    count_new<-c(count_new, length(rownames(marker)))
}
data<-cbind(count_total,count_new)
rownames(data)<-stage[2:7]
write.csv(data,file="H:/WRY/data/newRNA/BM/fig3/diff_mac.csv")
###############calculate diff
gene<-data.frame()
stage<-c("n00h","n02h","n06h","n12h","n24h","n48h","n72h")
type   <- c( "Macrophage")
for (i in 1:length(type)){
  for(j in 2:length(stage)){
    marker1<-FindMarkers(n,subset.ident = type[i],ident.1 = stage[j],ident.2 = stage[j-1],group.by = "orig.ident",logfc.threshold = log(1.5,base = 2),min.pct = 0.25,only.pos = FALSE)
    marker1$gene<-rownames(marker1)
    rownames(marker1)<-NULL
    marker2<-FindMarkers(n,subset.ident = type[i],ident.1 = stage[j],ident.2 = stage[1],group.by = "orig.ident",logfc.threshold = log(1.5,base = 2),min.pct = 0.25,only.pos = FALSE)
    marker2$gene<-rownames(marker2)
    rownames(marker2)<-NULL
    marker<-rbind(marker1,marker2)
    if(length(gene) == 0){
      gene<-marker
    }else if(length(rownames(marker))==0){
      gene<-gene
    }else{
      gene<-rbind(gene,marker)
    }
  }
}
gene<-gene[gene$p_val_adj<0.05,]
input<-unique(gene$gene)
write.csv(input,file="H:/WRY/data/newRNA/BM/fig3/total_1.5_0.05_mac.csv")

#################################
library(ComplexHeatmap)
library(circlize)
input<-read.csv(file="H:/WRY/data/newRNA/BM/fig3/total_1.5_0.05_mac.csv",row.names = 1)
new<-read.csv("H:/WRY/data/newRNA/BM/cluster/avg/avgtotal.csv",row.names = 1)
stage<-c("n00h","n02h","n06h","n12h","n24h","n48h","n72h")
type   <- c("Macrophage")
cat<-c()
for(i in 1:length(type)){
  for (j in 1:length(stage)) {
    cat<-c(cat,paste(type[i],stage[j],sep=""))
  }
}
data<-new[input[,1],cat]
sc<-data.frame(t(scale(t(data))))
data$max<-0
for (i in 1: length(rownames(data))) {
  data[i,"max"]<-which.max(data[i,])
}
data1<-data[order(data$max, decreasing = TRUE),]
sc<-sc[rownames(data1),]
gene<-c()
for(i in 1:7){
  data2<-sc[data1$max == i,]
  data3<-data2[order(data2[,i], decreasing = TRUE),]
  gene<-c(gene,rownames(data3))
}
data<-new[input[,1],cat]
sc<-data.frame(t(scale(t(data))))
sc<-sc[gene,]
write.csv(gene,file="H:/WRY/data/newRNA/BM/fig3/Macrophage_order.csv")

heatmap_show <- read.csv("H:/WRY/data/newRNA/BM/fig3/Macrophage_order_show1.csv", header=FALSE)
cat<-factor(cat,levels = cat)
type<-factor(type,levels = type)
datah<-sc
datah<-as.matrix(datah)
annotation_col <- data.frame(stage = factor(rep(stage, 1)), cell_type = factor(rep(type,each = 7)))
gene_select<-rownames(datah)[rownames(datah) %in% heatmap_show[,1]]
gene_pos<-which(rownames(datah) %in% gene_select)
row_anno<-rowAnnotation(gene=anno_mark(at=gene_pos,labels = gene_select,labels_gp = gpar(fontsize=18),side="left"))
col_fun <- colorRamp2(c(-2,0,2),c("#0000FF", "white", "#DD0000"))

pdf("H:/WRY/data/newRNA/BM/fig3/heatmap_total_1.5_0.05_mac_label1.pdf", width = 7, height = 10)
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

total<-read.csv("H:/WRY/data/newRNA/BM/cluster/avg/avgnew.csv",row.names = 1)
cat<-c()
for(i in 1:length(type)){
  for (j in 1:length(stage)) {
    cat<-c(cat,paste(type[i],stage[j],sep=""))
  }
}
data<-total[gene,cat]
sc<-data.frame(t(scale(t(data))))
sc<-na.omit(sc)
sc<-as.matrix(sc)

heatmap_show <- read.csv("H:/WRY/data/newRNA/BM/cluster/heatmap/Macrophage_order_show1.csv", header=FALSE)
cat<-factor(cat,levels = cat)
type<-factor(type,levels = type)
datah<-sc
datah<-as.matrix(datah)
annotation_col <- data.frame(stage = factor(rep(stage, 1)), cell_type = factor(rep(type,each = 7)))
gene_select<-rownames(datah)[rownames(datah) %in% heatmap_show[,1]]
gene_pos<-which(rownames(datah) %in% gene_select)
row_anno<-rowAnnotation(gene=anno_mark(at=gene_pos,labels = gene_select,labels_gp = gpar(fontsize=18),side="left"))
col_fun <- colorRamp2(c(-2,0,2),c("#0000FF", "white", "#DD0000"))

pdf("H:/WRY/data/newRNA/BM/fig3/heatmap_new_1.5_0.05_mac_label1.pdf", width = 7, height = 10)
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
###############
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)
library(stringr)
library(org.Mm.eg.db)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggrepel)
library(clusterProfiler)
n <- readRDS("H:/WRY/data/newRNA/BM/cluster/combind_clean_harmony_cluster.rds")
n<-subset(n,subset = cell_type == "Macrophage")
n06h<-read.csv(file="H:/WRY/data/newRNA/BM/dynamo_new/final/dynamic/var_mac6h_2.csv",row.names = 1)
n06h<-n06h[,c("nCells","nCounts","ntr","alpha","gamma","half_life")]
n72h<-read.csv(file="H:/WRY/data/newRNA/BM/dynamo_new/final/dynamic/var_mac72h.csv",row.names = 1)
n72h<-n72h[,c("nCells","nCounts","ntr","alpha","gamma","half_life")]
marker_6h<-FindMarkers(n,ident.1 = "n06h",ident.2 = "n00h",min.pct = 0.25,group.by = "orig.ident",only.pos = TRUE)
marker_6h<-marker_6h[marker_6h$p_val_adj < 0.05,]
marker_72h<-FindMarkers(n,ident.1 = "n72h",ident.2 = "n00h",min.pct = 0.25,group.by = "orig.ident",only.pos = TRUE)
marker_72h<-marker_72h[marker_72h$p_val_adj < 0.05,]

GO_term <- read.csv(file="H:/WRY/data/newRNA/BM/fig3/GO_for_3G.csv")
GO_term$ratio<-GO_term$count
names(GO_term)[1]<-"Description"
GO_term[GO_term$stage == "6h","ratio"]<-GO_term[GO_term$stage == "6h","ratio"]/length(rownames(marker_6h))
GO_term[GO_term$stage == "72h","ratio"]<-GO_term[GO_term$stage == "72h","ratio"]/length(rownames(marker_72h))
color<- c(red='#fe817d', blue ='#81b8df')
pdf("H:/WRY/data/newRNA/BM/fig3/GO_dynamic_count1.pdf", width = 3, height =3)
ggplot(GO_term,aes(half_life,alpha,size=count,color=color,alpha=color))+geom_point()+
  scale_size(range = c(3,7))+
  scale_color_manual(values=color)+
  scale_alpha_manual(values=c(0.6,0.6))+
  theme_classic()+NoLegend()+
  theme(axis.text = element_text(color="black",size = 13))+
  theme(axis.title = element_text(color="black",size = 13))+
  theme(axis.line = element_line(size=0.3, colour = 'black'))+
  theme(axis.ticks = element_line(size=0.3, colour = 'black'))
dev.off()
pdf("H:/WRY/data/newRNA/BM/fig3/GO_dynamic_ratio_label.pdf", width = 10, height = 10)
ggplot(GO_term,aes(half_life,alpha,color=color,size=ratio))+geom_point()+
  scale_size(range = c(2,5))+
  scale_color_manual(values=color)+
  theme_classic()+NoLegend()+
  theme(axis.text = element_text(color="black",size = 13))+
  theme(axis.title = element_text(color="black",size = 13))+
  theme(axis.line = element_line(size=0.3, colour = 'black'))+
  theme(axis.ticks = element_line(size=0.3, colour = 'black'))+
  geom_text(aes(label = Description,fontface="plain"))
dev.off()
pdf("H:/WRY/data/newRNA/BM/fig3/GO_dynamic_count_legend.pdf", width = 3.5, height = 3)
ggplot(GO_term,aes(half_life,alpha,color=color,size=count,alpha=color))+geom_point()+
  scale_size(range = c(3,7))+
  scale_color_manual(values=color)+
  theme_classic()+
  scale_alpha_manual(values=c(0.6,0.6))+
  theme(axis.text = element_text(color="black",size = 13))+
  theme(axis.title = element_text(color="black",size = 13))+
  theme(axis.line = element_line(size=0.3, colour = 'black'))+
  theme(axis.ticks = element_line(size=0.3, colour = 'black'))
dev.off()
###########################################
n <- readRDS("H:/WRY/data/newRNA/BM/cluster/combind_clean_harmony_cluster.rds")
n<-subset(n,subset = cell_type == "Macrophage")
marker_6h<-FindMarkers(n,ident.1 = "n06h",ident.2 = "n00h",group.by = "orig.ident",only.pos = TRUE)
marker_6h<-marker_6h[marker_6h$p_val_adj < 0.05,]
type <-  c("Macrophage")
go <- read.csv(file="H:/WRY/data/newRNA/BM/fig3/GO_6h.csv",header = FALSE)
names(go)<-c("ID","Description")
GO0h<-data.frame()
input <- bitr(rownames(marker_6h), fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db",drop = TRUE)
ego <- enrichGO(gene = input$ENTREZID, 
                OrgDb = org.Mm.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                pvalueCutoff = 1,
                qvalueCutoff = 1,
                readable = TRUE)
GO <- ego@result
id<-intersect(GO$Description,go$Description)
rownames(GO)<-GO$Description
GO <- GO[id,c("ID","Description","GeneRatio","qvalue","Count")]
  

GO0h$cell_type <- factor(GO0h$cell_type,levels = type)
GO0h<-na.omit(GO0h)
order_GO<-go[c(18:1),2]
GO0h$Description <- factor(GO0h$Description,levels = order_GO)
write.csv(GO0h,file="H:/WRY/data/newRNA/BM/sup2/GO_new_0h_select_information.csv")

#######################
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)
library(stringr)
library(org.Mm.eg.db)
library(ggplot2)
library(ggpubr)
library(cowplot)
GO_term <- read.csv(file="H:/WRY/data/newRNA/BM/fig3/GO_mac.csv",header = FALSE)
names(GO_term)<-c("ID","Description")
n <- readRDS("H:/WRY/data/newRNA/BM/cluster/combind_clean_harmony_cluster.rds")
n<-subset(n, subset = cell_type == "Macrophage")
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
names(datan)[c(3:19)]<-GO_term$Description

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
names(datat)[c(3:19)]<-GO_term$Description

datatgo<-aggregate(datat[,3:19],by=list(stage=datat$stage),FUN=mean)
datango<-aggregate(datan[,3:19],by=list(stage=datan$stage),FUN=mean)

stage<-c("n00h","n02h","n06h","n12h","n24h","n48h","n72h")

rownames(datatgo)<-datatgo$stage
datatgo<-datatgo[,2:18]
datatgo<-data.frame(t(datatgo))

datatgo<-datatgo[,stage]
sc<-data.frame(t(scale(t(datatgo))))
library(pheatmap)
pdf(file = "H:/WRY/data/newRNA/BM/fig3/heatmap_GO_Macrophage_total.pdf",width = 15,height = 8)
p1<-pheatmap(sc, scale = "none",cluster_cols = FALSE, cluster_rows = FALSE, fontsize = 16,
             show_rownames = TRUE,show_colnames = FALSE,
             color = colorRampPalette(c("#19499c", "white", "#ce191b"))(199), 
             legend = TRUE)
#             breaks = unique(c(seq(-4,0,length=100),seq(0,4,length=100))))
p1
dev.off()


rownames(datango)<-datango$stage
datango<-datango[,2:18]
datango<-data.frame(t(datango))

datango<-datango[,stage]
sc<-data.frame(t(scale(t(datango))))
#sc<-sc[order[,1],]
library(pheatmap)
pdf(file = "H:/WRY/data/newRNA/BM/fig3/heatmap_GO_Macrophage_new.pdf",width = 15,height = 8)
p2<-pheatmap(sc, scale = "none",cluster_cols = FALSE, cluster_rows = FALSE, fontsize = 16,
             show_rownames = TRUE,show_colnames = FALSE,
             color = colorRampPalette(c("#19499c", "white", "#ce191b"))(199), 
             legend = TRUE)
p2
dev.off()