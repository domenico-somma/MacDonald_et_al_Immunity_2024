#!/usr/bin/env Rscript 
library(dplyr)
library(Seurat)
library(httr)
library(readr)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(patchwork)
library(unixtools)
library(ggrepel)
library(repr)
library(ggmin)
library(harmony)
library(SeuratWrappers)
library(Nebulosa)
library(ggthemes)
set_config(config(ssl_verifypeer = 0L))
ulimit::memory_limit(200000)
set.tempdir("/datastore/lucy/tmp/")
setwd("/datastore/lucy/SynovialAtlas")
options(warn = -1, verbose=FALSE)

setwd("/datastore/lucy/SynovialAtlas")
set.tempdir("/datastore/JackF/td")

tcells <- readRDS("./cache/tcells.rds")

tcells

DimPlot(tcells, label=TRUE)

setwd("/datastore/JackF/TCellPaperRevision")

tcells <- subset(tcells, subset = highDoublets == "Doublet", invert=TRUE)

unique(tcells@meta.data$Tissue)

#Continue with typical preprocessing pipeline
tcells <- tcells %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData() %>% RunPCA()

ElbowPlot(tcells, ndims=40)

#UMAP preintegration
tcells <- RunUMAP(object = tcells, reduction = "pca", dims = 1:40)

options(repr.plot.width=30, repr.plot.height=10)
DimPlot(tcells, raster=FALSE, group.by=c("Donor","Experiment","Tissue"))

#Change theta for tissue to 1 at this stage
tcells <- tcells %>% 
    RunHarmony(group.by.vars=c("Donor","Experiment","Tissue"), plot_convergence = TRUE, theta=c(0,2,1))
#Run UMAP and clustering 
tcells <- RunUMAP(object = tcells, reduction = "harmony", dims = 1:40)

options(repr.plot.width=30, repr.plot.height=10)
DimPlot(tcells, raster=FALSE, group.by=c("Donor","Experiment","Tissue"))

options(repr.plot.width=30, repr.plot.height=10)
DimPlot(tcells, raster=FALSE, split.by="Experiment", group.by="celltype")

tcells <- FindNeighbors(object = tcells, reduction = "harmony", dims = 1:40)
tcells <- FindClusters(tcells, res=0.6)

options(repr.plot.width=10, repr.plot.height=10)
DimPlot(tcells, label=TRUE)

options(repr.plot.width=20, repr.plot.height=10)
DimPlot(tcells, label=TRUE, split.by="Tissue")

genes=c("CD8A","GZMK","GZMB","GZMH","CCL5","CCL4","CCR7","TIGIT","DUSP4","CD52","MAF","IL7R", "FOXP3","SELL")
options(repr.plot.width=20, repr.plot.height=20)
FeaturePlot(tcells, features = genes,order=TRUE)

options(repr.plot.width=20, repr.plot.height=20)
plot_density(tcells, features=genes, reduction="umap")

cluster.markers <- FindAllMarkers(tcells, only.pos=TRUE, test.use="MAST", min.pct=0.4)
cluster.markers <- cluster.markers[which(cluster.markers$p_val_adj<0.05),]

top10 <- cluster.markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)

options(repr.plot.width=20, repr.plot.height=20)
DoHeatmap(subset(tcells, downsample = 1000), features=unique(top10$gene)) + NoLegend()

top10[top10$cluster=="14",]

#Changed this to 11 as Lucy had 12 but random seed will alter populations slightly
tcells <- subset(tcells, idents=c("6", "11", "12"), invert=TRUE)

#Continue with typical preprocessing pipeline
tcells <- FindVariableFeatures(tcells, selection.method = "vst", nfeatures = 2000)
tcells <- ScaleData(tcells)
tcells <- RunPCA(tcells)
tcells <- tcells %>% 
    RunHarmony(group.by.vars=c("Donor","Experiment","Tissue"), plot_convergence = TRUE, theta=c(0,2,1))

tcells <- RunUMAP(object = tcells, reduction = "harmony", dims = 1:40)
tcells <- FindNeighbors(object = tcells, reduction = "harmony", dims = 1:40)
tcells <- FindClusters(tcells, res=1.0)

options(repr.plot.width=10, repr.plot.height=10)
DimPlot(tcells, label=TRUE)

options(repr.plot.width=20, repr.plot.height=10)
DimPlot(tcells, label=TRUE, split.by="Tissue")

options(repr.plot.width=20, repr.plot.height=10)
DimPlot(tcells, label=TRUE, split.by="Experiment")

#Continue with typical preprocessing pipeline
tcells <- FindVariableFeatures(tcells, selection.method = "vst", nfeatures = 2000)
tcells <- ScaleData(tcells)
tcells <- RunPCA(tcells)
tcells <- tcells %>% 
    RunHarmony(group.by.vars=c("Donor","Experiment","Tissue"), plot_convergence = TRUE, theta=c(0,2,1))

tcells <- RunUMAP(object = tcells, reduction = "harmony", dims = 1:40)
tcells <- FindNeighbors(object = tcells, reduction = "harmony", dims = 1:40)
tcells <- FindClusters(tcells, res=1.3)

options(repr.plot.width=10, repr.plot.height=10)
DimPlot(tcells, label=TRUE)

options(repr.plot.width=20, repr.plot.height=10)
DimPlot(tcells, label=TRUE, split.by="Tissue")

options(repr.plot.width=30, repr.plot.height=10)
DimPlot(tcells, label=TRUE, split.by="Experiment")

genes=c("CD8A","GZMK","GZMB","GZMH","CCL5","CCL4","CCR7","TIGIT","DUSP4","CD52","MAF","IL7R", "FOXP3","SELL","LTB","GNLY", 'TNFSF13B', "ADGRG1")
options(repr.plot.width=20, repr.plot.height=20)
FeaturePlot(tcells, features = genes,order=TRUE)

options(repr.plot.width=25, repr.plot.height=20)
plot_density(tcells, features="CD4", reduction="umap")

options(repr.plot.width=25, repr.plot.height=20)
DefaultAssay(tcells)<-'ADT'
plot_density(tcells, features="Hu.GPR56", reduction="umap")

tcells.subset <- subset(tcells, downsample = 1000)

DefaultAssay(tcells.subset)<-'RNA'

cluster.markers <- FindAllMarkers(tcells.subset, only.pos=TRUE, test.use="MAST", min.pct=0.4)
cluster.markers <- cluster.markers[which(cluster.markers$p_val_adj<0.05),]

options(repr.plot.width=10, repr.plot.height=10)

BuildClusterTree(tcells) %>% PlotClusterTree()

top10 <- cluster.markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)

top10[top10$cluster=="12",]

DefaultAssay(tcells)<-"RNA"
options(repr.plot.width=30, repr.plot.height=8)
x <- list()
for(gene in c("KLRB1")){
  x[[gene]] <- VlnPlot(tcells, features = gene, pt.size=0) + theme_linedraw() + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.text.y= element_text(size = 15), axis.title.y=element_text(size = 20), plot.title=element_text(size=30))+ stat_summary(fun.y = median, geom='point', size = 2, colour = "black") 
}
plot_grid(plotlist=x, ncol=3)

DefaultAssay(tcells)<-"RNA"
options(repr.plot.width=30, repr.plot.height=8)
x <- list()
for(gene in c("CCR7", 'CD8B', 'CD8A')){
  x[[gene]] <- VlnPlot(tcells, features = gene, pt.size=0) + theme_linedraw() + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.text.y= element_text(size = 15), axis.title.y=element_text(size = 20), plot.title=element_text(size=30))+ stat_summary(fun.y = median, geom='point', size = 2, colour = "black") 
}
plot_grid(plotlist=x, ncol=3)

rownames(tcells@assays$ADT)

DefaultAssay(tcells)<-"RNA"
options(repr.plot.width=30, repr.plot.height=8)
x <- list()
for(gene in c("CD8A", "GZMK",'RORA', "GPR183", "SELL", "CCR7")){
  x[[gene]] <- VlnPlot(tcells, features = gene, pt.size=0) + theme_linedraw() + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.text.y= element_text(size = 15), axis.title.y=element_text(size = 20), plot.title=element_text(size=30))+ stat_summary(fun.y = median, geom='point', size = 2, colour = "black") 
}
plot_grid(plotlist=x, ncol=3)

options(repr.plot.width=20, repr.plot.height=20)
DoHeatmap(tcells.subset, features=unique(top10$gene)) + NoLegend()

#Jack T-Cell Cluster Annotation
tcells$deepclusters <- plyr::revalue(as.character(tcells$seurat_clusters),
                                        c("0"='CCR7+SELL+CD4+ Naive',
                                          "1"="CCR7+LTB+CD4+ TCM",
                                          "2"='GZMK+CCL5+CD8A+ TEM',
                                          "3"="CCR7+SELL+CD4+ Naive",
                                          "4"='KLRC2+ NK',
                                          "5"="MYOM2+ NK",
                                          "6"="GZMH+CCL5+CD8A+ TEM",
                                          "7"="ANXA1+CD4+ TEM",
                                          "8"="CCR7+AIF1+CD8B+ Naive",
                                          "9"="RORA+CD4+ TEM",
                                          "10"="KLRB1+ NKT",
                                          "11"="XCL1+ NK",
                                          "12"="CCL5+CD8B+CD8A+ Naive",
                                          "13"='MAF+ TPH',
                                         "14"='FOXP3+ Tregs',
                                         "15"="GZMH+CCL5+CD8A+ TEM",
                                         "16"="CCL5+CXCR6+MAF+ TPH",
                                         "17"='GZMK+ NKT',
                                         "18"='CCR7+SELL+CD4+ Naive',
                                         "19"='CCR7+SELL+CD4+ Naive',
                                         "20"='GPR183+ZFP36+ Naive',
                                         "21"='ISG15 Activated CD4'))


#Jack T-Cell Cluster Annotation
tcells$shallowclusters <- plyr::revalue(as.character(tcells$seurat_clusters),
                                        c("0"='CCR7+SELL+CD4+ Naive',
                                          "1"="CCR7+LTB+CD4+ TCM",
                                          "2"='CD8 T-Cell',
                                          "3"="CCR7+SELL+CD4+ Naive",
                                          "4"='NK',
                                          "5"="NK",
                                          "6"="CD8 T-Cell",
                                          "7"="ANXA1+CD4+ TEM",
                                          "8"="CD8 T-Cell",
                                          "9"="RORA+CD4+ TEM",
                                          "10"="NKT",
                                          "11"="NK",
                                          "12"="CD8 T-Cell",
                                          "13"='MAF+ TPH',
                                         "14"='FOXP3+ Tregs',
                                         "15"="CD8 T-Cell",
                                         "16"="CCL5+CXCR6+MAF+ TPH",
                                         "17"='NKT',
                                         "18"='CCR7+SELL+CD4+ Naive',
                                         "19"='CCR7+SELL+CD4+ Naive',
                                         "20"='GPR183+ZFP36+ Naive',
                                         "21"='ISG15 Activated CD4'))

options(repr.plot.width=10, repr.plot.height=10)

BuildClusterTree(tcells) %>% PlotClusterTree()

Idents(tcells)<-'deepclusters'
options(repr.plot.width=15, repr.plot.height=10)
DimPlot(tcells, label=TRUE)

Idents(tcells)<-'shallowclusters'
options(repr.plot.width=15, repr.plot.height=10)
DimPlot(tcells, label=TRUE)

options(repr.plot.width=15, repr.plot.height=10)
DimPlot(tcells, label=TRUE)

options(repr.plot.width=30, repr.plot.height=10)
DimPlot(tcells, label=TRUE, split.by='Tissue')

#Create data frame for information required to plot 
umap <- data.frame(rownames(tcells@meta.data),Embeddings(tcells, reduction = "umap")[,1], Embeddings(tcells, reduction = "umap")[,2],
                   tcells$Tissue,
                   tcells$SampleID, Idents(tcells))
colnames(umap) <- c("ID","UMAP1", "UMAP2", "Tissue", "Sample", "Celltype")

#Save
write.csv(umap, "single.cell-umap.plot.csv")

allcells <-  as.data.frame(cbind(umap$UMAP1,umap$UMAP2))
colnames(allcells) <- c("UMAP1", "UMAP2")

#Create an image showing split UMAPs of all conditions
options(repr.plot.width=14, repr.plot.height=8)
ggplot(umap, aes(x=UMAP1, y=UMAP2)) + 
  geom_point(data=allcells, colour="grey89", fill="grey89", size=3, shape=21) +
  geom_point(aes(fill=factor(Celltype)),color="black", size=2, shape=21, stroke = 1/4) + theme_powerpoint() + ggtitle("") + theme(axis.text=element_text(size=20),axis.title=element_text(size=20))

library(gridExtra)

allcells <-  as.data.frame(cbind(umap$UMAP1,umap$UMAP2))
colnames(allcells) <- c("UMAP1", "UMAP2")

#Create an image showing split UMAPs of all conditions
options(repr.plot.width=72, repr.plot.height=8)
ggplot(umap, aes(x=UMAP1, y=UMAP2)) + facet_grid(~Celltype, switch="x") + 
  geom_point(data=allcells, colour="grey89", fill="grey89", size=3, shape=21) +
  geom_point(aes(fill=factor(Celltype)),color="black", size=3, shape=21, stroke = 1/4) + theme_powerpoint() + ggtitle("") + guides(fill=FALSE) + theme(axis.text=element_text(size=20),axis.title=element_text(size=20))

allcells <-  as.data.frame(cbind(umap$UMAP1,umap$UMAP2))
colnames(allcells) <- c("UMAP1", "UMAP2")

#Create an image showing split UMAPs of all conditions
options(repr.plot.width=15, repr.plot.height=8)
ggplot(umap, aes(x=UMAP1, y=UMAP2)) + facet_grid(~Tissue, switch="x") + 
  geom_point(data=allcells, colour="grey89", fill="grey89", size=3, shape=21) +
  geom_point(aes(fill=factor(Celltype)),color="black", size=3, shape=21, stroke = 1/4) + theme_powerpoint() + ggtitle("") + guides(fill=FALSE) + theme(axis.text=element_text(size=20),axis.title=element_text(size=20))

options(repr.plot.width=20, repr.plot.height=10)
plot_density(tcells, features=c('CD3D', "CD8A", 'GNLY', 'CD4'), reduction="umap")

unique(tcells$Experiment)

Idents(tcells)<-'Experiment'
glas3<-subset(tcells, idents=c('Glasgow3'))

Idents(tcells)<-'clusters'

rownames(tcells@assays$ADT)

DefaultAssay(tcells)<-'ADT'
plot_density(tcells, features=c('Hu.CD4-RPA.T4', 'Hu.CD8'), reduction="umap")

DefaultAssay(glas3)<-'ADT'
options(repr.plot.width=15, )
plot_density(glas3, features=c('Hu.CD4-RPA.T4', 'Hu.CD8'), reduction="umap")

DefaultAssay(tcells)<-'RNA'
FindMarkers(tcells, ident.1='21', only.pos=TRUE)



top10[top10$cluster=="21",]

cluster.markers[cluster.markers$cluster=="20",]

options(repr.plot.width=20, repr.plot.height=20)

BuildClusterTree(tcells) %>% PlotClusterTree()

getwd()

saveRDS(tcells, 'jacktcells.rds')

tcells<-readRDS('jacktcells.rds')
tcells

options(repr.plot.width=15, repr.plot.height=11)
plot_density(tcells, features=c("IFNG", "IFNGR1"), reduction="umap")

Idents(tcells)<-'seurat_clusters'

Idents(tcells)<-'Tissue'

cd4t<-subset(tcells, idents=c('ST'))

unique(cd4t$group)

cd4t$subgroup <- cd4t$group

cd4t$group <- plyr::revalue(as.character(cd4t$group),
                                c('Naive RA'='Active RA',
                                 'Resistant to cDMARDs RA'='Active RA',
                                 'Resistant to bDMARDs RA'='Active RA',
                                 'Difficult to treat RA'='Active RA'))

cd4t$replicate <- cd4t$Unique_ID


unique(cd4t$shallowclusters)

Idents(cd4t)<-'shallowclusters'
cd4t<-subset(cd4t, idents=c('NK', 'CD8 T-Cell', 'NKT'), invert=TRUE)

options(repr.plot.width=32, repr.plot.height=20)
AbundancePlot(cd4t, group.by = "shallowclusters", split.by = "subgroup", replicate.by = "replicate", perform.stat.test=FALSE)

options(repr.plot.width=25, repr.plot.height=20)
AbundancePlot(cd4t, group.by = "shallowclusters", split.by = "group", replicate.by = "replicate", perform.stat.test=FALSE)

unique(cd4t$group)

unique(cd4t$condition)

cd4t$group <- as.ordered(factor(cd4t$group, levels=c("Healthy", "Active RA", "Remission RA")))

saveRDS(tcells, 'jacktcells.rds')

devtools::install_github("diegoalexespi/pochi")


library(pochi)

DefaultAssay(cd4t)<-"RNA"
Idents(cd4t)<-'shallowclusters'
options(repr.plot.width=30, repr.plot.height=8)
x <- list()
for(gene in c("CCL5")){
  x[[gene]] <- VlnPlot(cd4t, features = gene, pt.size=0) + theme_linedraw() + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.text.y= element_text(size = 15), axis.title.y=element_text(size = 20), plot.title=element_text(size=30))+ stat_summary(fun.y = median, geom='point', size = 2, colour = "black") 
}
plot_grid(plotlist=x, ncol=3)

treg<-subset(cd4t, idents=c('FOXP3+ Tregs'))

maf<-subset(cd4t, idents=c('MAF+ TPH'))

Idents(treg)<-'group'

DefaultAssay(treg)<-"RNA"
Idents(treg)<-'group'
options(repr.plot.width=30, repr.plot.height=8)
x <- list()
for(gene in c("CTLA4")){
  x[[gene]] <- VlnPlot(treg, features = gene, pt.size=1) + theme_linedraw() + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.text.y= element_text(size = 15), axis.title.y=element_text(size = 20), plot.title=element_text(size=30))+ stat_summary(fun.y = median, geom='point', size = 2, colour = "black") 
}
plot_grid(plotlist=x, ncol=3)

table(treg$group)

options(repr.plot.width=15, repr.plot.height=11)
plot_density(cd4t, features=c("CCL5"), reduction="umap")

tcells$subgroup <- tcells$group

tcells$group <- plyr::revalue(as.character(tcells$group),
                                c('Naive RA'='Active RA',
                                 'Resistant to cDMARDs RA'='Active RA',
                                 'Resistant to bDMARDs RA'='Active RA',
                                 'Difficult to treat RA'='Active RA'))



table(cd4t$group, cd4t$shallowclusters)
