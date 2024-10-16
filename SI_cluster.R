library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)
library(harmony)
s00<-readRDS("H:/WRY/data/newRNA/SI/process/SI_clean/SI0h_clean.rds")
s02<-readRDS("H:/WRY/data/newRNA/SI/process/SI_clean/SI2h_clean.rds")
s06<-readRDS("H:/WRY/data/newRNA/SI/process/SI_clean/SI6h_clean.rds")
s72<-readRDS("H:/WRY/data/newRNA/SI/process/SI_clean/SI72h_clean.rds")
s <- merge(s00, c(s02,s06,s72), add.cell.ids = c("s00h", "s02h", "s06h","s72h"), project = "SI")
saveRDS(s, file = "H:/WRY/data/newRNA/SI/cluster/merge_clean.rds")

cc<-read.csv("H:/WRY/code/ccgene/mouse_cc.csv",row.names = 1)
s.genes <- cc[cc$stage == "S","symbol"]
g2m.genes <- cc[cc$stage == "G2/M","symbol"]

n <- readRDS("H:/WRY/data/newRNA/SI/cluster/merge_clean.rds")
n@meta.data<-n@meta.data[rownames(n@meta.data),c("orig.ident","nCount_originalexp","nFeature_originalexp","percent.mt","percent.ribo","percent.dissociation")]
VlnPlot(n, features = c("nFeature_originalexp", "nCount_originalexp", "percent.mt"), ncol = 3)
n <- NormalizeData(n,normalization.method = "LogNormalize",scale.factor = 10000,assay="originalexp")
saveRDS(n,"H:/WRY/data/newRNA/SI/dynamo/merge_clean_SI_norm.rds")
n <- FindVariableFeatures(n)
n <- CellCycleScoring(n, s.features = s.genes, g2m.features = g2m.genes)
n$CC.Difference <- n$S.Score - n$G2M.Score
n <- ScaleData(n,vars.to.regress = c('nCount_originalexp', 'percent.mt', 'percent.ribo', 'percent.dissociation', 'CC.Difference'))
n <- RunPCA(n)
n <- RunHarmony(n, group.by.vars = "orig.ident",assay.use = "originalexp",plot_convergence = TRUE)
n <- RunUMAP(n, reduction = "harmony", dims = 1:30)
n <- FindNeighbors(n, reduction = "harmony", dims = 1:30) 
n <- FindClusters(n,resolution=1)
pdf(file = "H:/WRY/data/newRNA/SI/cluster/umap_harmony.pdf",width = 6,height = 5)
DimPlot(n, reduction = "umap",label = TRUE)
dev.off()
pdf(file = "H:/WRY/data/newRNA/SI/cluster/umap_ori_harmony.pdf",width = 21,height = 5)
DimPlot(n, reduction = "umap", split.by = "orig.ident")
dev.off()
saveRDS(n, file = "H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony.rds")
FeaturePlot(n, features = c("Igkc"), reduction = "umap", cols = c("grey", "blue"), label = TRUE, label.size = 2,pt.size = 0.3)

n <- readRDS("H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony.rds")
new.cluster.ids <- c("E",
                     "E",
                     "E",
                     "Immune",
                     "E",
                     "E",
                     "E",
                     "E",
                     "Immune",
                     "E",
                     "E",
                     "Immune",
                     "E",
                     "Immune",
                     "Immune",
                     "Immune",
                     "N",
                     "Immune",
                     "E",
                     "E",
                     "E",
                     "E",
                     "N",
                     "Immune",
                     "N",
                     "N",
                     "N",
                     "N",
                     "Immune")
names(new.cluster.ids) <- levels(n)
n <- RenameIdents(n, new.cluster.ids)
DimPlot(n, reduction = "umap",label = TRUE)
n$cell_type <- Idents(n)
saveRDS(n, file = "H:/WRY/data/newRNA/SI/cluster/combind_clean_harmony_cluster.rds")