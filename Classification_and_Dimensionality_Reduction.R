# Import packages
library(Seurat)
library(dplyr)
library(Matrix)
library(magrittr)
library(stringr)
#library(tidyverse)
library(ggplot2)
options(stringsAsFactors = F)
options(as.is = T)
# Color scheme
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175')
Mycol_2 <-c('#FF7F0EFF', '#1F77B4FF')
Mycol_3 <- c("#1B9E77","#D95F02","#7570B3")
#library("ggsci")
#Mycol_nejm= pal_nejm("default", alpha =0.7)(8)##Extract colors
Mycol_4 <- c("#A6CEE3" ,"#1F78B4", "#B2DF8A" ,"#33A02C" ,"#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00", "#CAB2D6" ,"#6A3D9A" ,"#FFFF99", "#B15928","#1B9E77","#D95F02" ,"#7570B3", "#E7298A" ,"#66A61E" ,"#E6AB02", "#A6761D","#666666")
Mycol_5<-c("#F8766D","#CD9600","#7CAE00","#0CB702","#00C19A","#00B8E7","#8494FF","#C77CFF","#FF61CC")
load('sce1.Rdata')
load('sce2.Rdata')
load('sce3.Rdata')
load('sce4.Rdata')
load('All_Data.Rdata')

Re_Pre_Data_1 <-sce1
Re_Post_Data_1  <- sce2
NoR_Pre_Data_1 <-sce3
NoR_Post_Data_1  <- sce4
Re_All_Data_1 <- Reduce(merge,list(sce1,sce2))
NoR_All_Data_1 <- Reduce(merge,list(sce3,sce4))

##### (1) Take log###################################
Re_All_Data_1 = NormalizeData(object = Re_All_Data_1,normalization.method =  "LogNormalize",  scale.factor = 1e6)
NoR_All_Data_1= NormalizeData(object = NoR_All_Data_1,normalization.method =  "LogNormalize",  scale.factor = 1e6)
Re_Pre_Data_1 = NormalizeData(object = Re_Pre_Data_1 ,normalization.method =  "LogNormalize",  scale.factor = 1e6)
Re_Post_Data_1 = NormalizeData(object = Re_Post_Data_1,normalization.method =  "LogNormalize",  scale.factor = 1e6)
NoR_Pre_Data_1 = NormalizeData(object = NoR_Pre_Data_1,normalization.method =  "LogNormalize",  scale.factor = 1e6)
NoR_Post_Data_1=NormalizeData(object=NoR_Post_Data_1,normalization.method =  "LogNormalize",  scale.factor = 1e6)
# NormalizeData() divides each gene's count by total, then multiplies by a scale.factor, then converts with natural logarithm
# Aims to eliminate the impact of different cell sequencing depths, taking log reduces data dispersion, +1 prevents read=0 case
# 2. Feature extraction
Re_All_Data_1 = FindVariableFeatures(object = Re_All_Data_1,selection.method = "vst", nfeatures = 4000)
NoR_All_Data_1 = FindVariableFeatures(object = NoR_All_Data_1,selection.method = "vst", nfeatures = 4000)
Re_Pre_Data_1 = FindVariableFeatures(object = Re_Pre_Data_1,selection.method = "vst", nfeatures = 4000)
Re_Post_Data_1 = FindVariableFeatures(object = Re_Post_Data_1,selection.method = "vst", nfeatures = 4000)
NoR_Pre_Data_1 = FindVariableFeatures(object = NoR_Pre_Data_1,selection.method = "vst", nfeatures = 4000)
NoR_Post_Data_1 = FindVariableFeatures(object = NoR_Post_Data_1,selection.method = "vst", nfeatures = 4000)
# FindVariableFeatures is hard filtering, based on some statistical metrics like sd, mad, vst etc. to judge from the 55,000+ genes in your input single-cell expression matrix, the most important 4000 genes, the remaining 52,000 genes are not considered in downstream analysis
# See what the most important 4000 genes are (top six)
head(VariableFeatures(Re_All_Data_1))  
head(VariableFeatures(NoR_All_Data))  
head(VariableFeatures(Re_Pre_Data))  
head(VariableFeatures(Re_Post_Data))  
head(VariableFeatures(NoR_Pre_Data))  
head(VariableFeatures(NoR_Post_Data))  
# 3. Normalization
Re_All_Data_1 = ScaleData(object = Re_All_Data_1)
NoR_All_Data_1 = ScaleData(object = NoR_All_Data_1)
Re_Pre_Data_1 = ScaleData(object = Re_Pre_Data_1)
Re_Post_Data_1 = ScaleData(object = Re_Post_Data_1)
NoR_Pre_Data_1 = ScaleData(object = NoR_Pre_Data_1)
NoR_Post_Data_1 = ScaleData(object = NoR_Post_Data_1)
# 4. PCA dimensionality reduction
Re_All_Data_1 =  RunPCA(object = Re_All_Data_1 ,do.print = FALSE)
NoR_All_Data_1 = RunPCA(object =NoR_All_Data_1 ,do.print = FALSE)
Re_Pre_Data_1 =  RunPCA(object =Re_Pre_Data_1 , do.print = FALSE)
Re_Post_Data_1 = RunPCA(object = Re_Post_Data_1,do.print = FALSE)
NoR_Pre_Data_1 = RunPCA(object =NoR_Pre_Data_1 , do.print = FALSE)
NoR_Post_Data_1= RunPCA(object =NoR_Post_Data_1 ,do.print = FALSE)
# 5000 genes will be transformed into 5000 dimensions, but we usually only look at the first dozen dimensions, so it's also a very efficient dimensionality reduction method
# Originally each cell has 5000 genes, each gene has expression levels in various cells, after doing PCA analysis, each cell has 50 PCs, each PC has embedding values in various cells. For ease of understanding, you can think of PCs as meta-genes that integrate multiple genes, think of cell embedding as the expression levels of these meta-genes
# RunPCA function calculates 50 PCs by default, if you want to see more PCs, you can modify the npcs parameter of RunPCA.

##TSNE clustering analysis
pcSelect=15
Re_All_Data_1     <-  FindNeighbors(object = Re_All_Data_1 ,dims = 1:pcSelect)    
NoR_All_Data_1  <- FindNeighbors(object =NoR_All_Data_1, dims = 1:pcSelect)    
Re_Pre_Data_1    <- FindNeighbors(object =Re_Pre_Data_1 , dims = 1:pcSelect)    
Re_Post_Data_1   <- FindNeighbors(object = Re_Post_Data_1 ,dims = 1:pcSelect)    
NoR_Pre_Data_1  <- FindNeighbors(object =NoR_Pre_Data_1 , dims = 1:pcSelect)    
NoR_Post_Data_1<- FindNeighbors(object =NoR_Post_Data_1, dims = 1:pcSelect)               
# Calculate adjacency distance
Re_All_Data_1     <- FindClusters(object = Re_All_Data_1  ,    resolution = 0.4)  
NoR_All_Data_1  <- FindClusters(object = NoR_All_Data_1,   resolution = 0.4)  
Re_Pre_Data_1    <- FindClusters(object = Re_Pre_Data_1   ,  resolution = 0.4)  
Re_Post_Data_1   <- FindClusters(object =Re_Post_Data_1   , resolution = 0.4)  
NoR_Pre_Data_1  <- FindClusters(object =NoR_Pre_Data_1  , resolution = 0.4)  
NoR_Post_Data_1<- FindClusters(object =NoR_Post_Data_1, resolution = 0.4)  
# Group cells, optimize standard modularity
Re_All_Data_1      <- RunTSNE(object = Re_All_Data_1  ,dims = 1:pcSelect)         
NoR_All_Data_1   <- RunTSNE(object =NoR_All_Data_1 ,dims = 1:pcSelect)                      
Re_Pre_Data_1     <- RunTSNE(object =Re_Pre_Data_1 , dims = 1:pcSelect)                      
Re_Post_Data_1    <- RunTSNE(object = Re_Post_Data_1 ,dims = 1:pcSelect)                      
NoR_Pre_Data_1   <- RunTSNE(object =NoR_Pre_Data_1, dims = 1:pcSelect)                      
NoR_Post_Data_1 <- RunTSNE(object= NoR_Post_Data_1 ,dims = 1:pcSelect)                                   
#TSNE clustering
Re_All_tSNE_image <- TSNEPlot(object = Re_All_Data_1, pt.size = 2, label = TRUE)    
NoR_All_tSNE_image <- TSNEPlot(object = NoR_All_Data_1, pt.size = 2, label = TRUE)    
Re_Pre_tSNE_image <- TSNEPlot(object = Re_Pre_Data_1, pt.size = 2, label = TRUE)    
Re_Post_tSNE_image <- TSNEPlot(object = Re_Post_Data_1, pt.size = 2, label = TRUE)    
NoR_Pre_tSNE_image <- TSNEPlot(object = NoR_Pre_Data_1 ,pt.size = 2, label = TRUE)    
NoR_Post_tSNE_image <- TSNEPlot(object =NoR_Post_Data_1, pt.size = 2, label = TRUE)    
#TSNE visualization
ggsave("177_Re_All_1_tSNE_image.png", plot = Re_All_tSNE_image, width = 10, height = 10)
ggsave("178_NoR_All_1_tSNE_image.png", plot = NoR_All_tSNE_image, width = 10, height = 10)
ggsave("179_Re_Pre_1_tSNE_image.png", plot = Re_Pre_tSNE_image, width = 10, height = 10)
ggsave("180_Re_Post_1_tSNE_image.png", plot = Re_Post_tSNE_image, width = 10, height = 10)
ggsave("181_NoR_Pre_1_tSNE_image.png", plot = NoR_Pre_tSNE_image, width = 10, height = 10)
ggsave("182_NoR_Post_1_tSNE_image.png", plot = NoR_Post_tSNE_image, width = 10, height = 10)
      
  
table(Re_All_Data_1     $RNA_snn_res.0.4)  #11
table(NoR_All_Data_1   $RNA_snn_res.0.4)  #10
table(Re_Pre_Data_1     $RNA_snn_res.0.4)  #8
table(Re_Post_Data_1    $RNA_snn_res.0.4)  #10
table(NoR_Pre_Data_1   $RNA_snn_res.0.4)  #9
table(NoR_Post_Data_1  $RNA_snn_res.0.4)  #11
   

cellTypes1 <- New_All_Data_2@meta.data$cellTypes[1: 2715]
cellTypes2 <- New_All_Data_2@meta.data$cellTypes[2716: 11653]
cellTypes3 <- New_All_Data_2@meta.data$cellTypes[1:1191]
cellTypes4 <- New_All_Data_2@meta.data$cellTypes[1192:2715]
cellTypes5 <- New_All_Data_2@meta.data$cellTypes[2716:5319]
cellTypes6 <- New_All_Data_2@meta.data$cellTypes[5320:11653]

Re_All_Data_1<- AddMetaData(Re_All_Data_1, cellTypes1,col.name = "cellTypes")
NoR_All_Data_1<- AddMetaData(NoR_All_Data_1, cellTypes2,col.name = "cellTypes")
Re_Pre_Data_1<- AddMetaData(Re_Pre_Data_1, cellTypes3,col.name = "cellTypes")
Re_Post_Data_1<- AddMetaData(Re_Post_Data_1, cellTypes4,col.name = "cellTypes")
NoR_Pre_Data_1<- AddMetaData(NoR_Pre_Data_1, cellTypes5,col.name = "cellTypes")
NoR_Post_Data_1<- AddMetaData(NoR_Post_Data_1, cellTypes6,col.name = "cellTypes")
Re_All_Data_1$cellTypes [Re_All_Data_1$cellTypes== "T cells"] = "γδT cells"
NoR_All_Data_1$cellTypes [NoR_All_Data_1$cellTypes== "T cells"] = "γδT cells"
Re_Pre_Data_1$cellTypes [Re_Pre_Data_1$cellTypes== "T cells"] = "γδT cells"
Re_Post_Data_1$cellTypes [Re_Post_Data_1$cellTypes== "T cells"] = "γδT cells"
NoR_Pre_Data_1$cellTypes [NoR_Pre_Data_1$cellTypes== "T cells"] = "γδT cells"
NoR_Post_Data_1$cellTypes [NoR_Post_Data_1$cellTypes== "T cells"] = "γδT cells"

Idents(Re_All_Data_1) <- Re_All_Data_1$cellTypes
Idents(NoR_All_Data_1) <- NoR_All_Data_1$cellTypes
Idents(Re_Pre_Data_1) <- Re_Pre_Data_1$cellTypes
Idents(Re_Post_Data_1) <- Re_Post_Data_1$cellTypes
Idents(NoR_Pre_Data_1) <- NoR_Pre_Data_1$cellTypes
Idents(NoR_Post_Data_1) <- NoR_Post_Data_1$cellTypes



###Cluster colors
CD8_T_cells_col= Mycol_5[3]
CD4_T_cells_col=Mycol_5[2]
Regulatory_T_cells_col= Mycol_5[6]
Monocyte_col= Mycol_5[4]
NKT_cells_col= Mycol_5[1]
B_cells_col= Mycol_5[7]
T_cells_col= Mycol_5[9]
pDCs_col= Mycol_5[5]
Plasma_cells_col= Mycol_5[8]



p1 <- DimPlot(Re_All_Data_1, reduction = "tsne", label = TRUE, pt.size = 1, cols= c(CD4_T_cells_col,Regulatory_T_cells_col , NKT_cells_col,CD8_T_cells_col, Monocyte_col, B_cells_col,pDCs_col, Plasma_cells_col, T_cells_col))
p2 <- DimPlot(NoR_All_Data_1, reduction = "tsne", label = TRUE, pt.size = 1, cols= c(NKT_cells_col,CD4_T_cells_col, CD8_T_cells_col, Monocyte_col, pDCs_col, Regulatory_T_cells_col , B_cells_col, Plasma_cells_col, T_cells_col))
p3 <- DimPlot(Re_Pre_Data_1, reduction = "tsne", label = TRUE, pt.size = 1, cols= c(CD4_T_cells_col,Regulatory_T_cells_col , NKT_cells_col,CD8_T_cells_col, Monocyte_col, B_cells_col,pDCs_col, Plasma_cells_col))
p4 <- DimPlot(Re_Post_Data_1, reduction = "tsne", label = TRUE, pt.size = 1, cols= c(CD8_T_cells_col, Regulatory_T_cells_col , B_cells_col,CD4_T_cells_col, NKT_cells_col, Monocyte_col, pDCs_col, T_cells_col ,Plasma_cells_col))
p5 <- DimPlot(NoR_Pre_Data_1, reduction = "tsne", label = TRUE, pt.size = 1, cols= c(NKT_cells_col,CD4_T_cells_col, CD8_T_cells_col, Monocyte_col, pDCs_col, Regulatory_T_cells_col , B_cells_col, Plasma_cells_col, T_cells_col))
p6 <- DimPlot(NoR_Post_Data_1 ,reduction = "tsne", label = TRUE, pt.size = 1, cols= c(CD8_T_cells_col, CD4_T_cells_col,Regulatory_T_cells_col , B_cells_col, Monocyte_col, NKT_cells_col, Plasma_cells_col, pDCs_col, T_cells_col))

ggsave("183_group1_subcluster_naming.png", p1, width=10 ,height=8)
ggsave("184_group2_subcluster_naming.png", p2, width=10 ,height=8)
ggsave("185_group3_subcluster_naming.png", p3, width=10 ,height=8)
ggsave("186_group4_subcluster_naming.png", p4, width=10 ,height=8)
ggsave("187_group5_subcluster_naming.png", p5, width=10 ,height=8)
ggsave("188_group6_subcluster_naming.png", p6, width=10 ,height=8)
  
  
  
# 9 9 8 9 9 9 
save(Re_All_Data_1      ,file='Re_All_Data_1.Rdata')
save(NoR_All_Data_1   ,file='NoR_All_Data_1.Rdata')
save(Re_Pre_Data_1     ,file='Re_Pre_Data_1.Rdata')
save(Re_Post_Data_1    ,file='Re_Post_Data_1.Rdata')
save(NoR_Pre_Data_1  ,file='NoR_Pre_Data_1.Rdata')
save(NoR_Post_Data_1,file='NoR_Post_Data_1.Rdata')

#Differential genes
# Find all differential genes (will take very long)
Re_All_Data_1.markers <- FindAllMarkers(object = Re_All_Data_1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
NoR_All_Data_1.markers <- FindAllMarkers(object = NoR_All_Data_1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
Re_Pre_Data_1.markers <- FindAllMarkers(object = Re_Pre_Data_1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
Re_Post_Data_1.markers <- FindAllMarkers(object = Re_Post_Data_1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
NoR_Pre_Data_1.markers <- FindAllMarkers(object = NoR_Pre_Data_1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
NoR_Post_Data_1.markers <- FindAllMarkers(object = NoR_Post_Data_1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(Re_All_Data_1.markers,file="_1_Re_All_markers.csv")
write.csv(NoR_All_Data_1.markers,file="_1_NoR_All_markers.csv")
write.csv(Re_Pre_Data_1.markers,file="_1_Re_Pre_markers.csv")
write.csv(Re_Post_Data_1.markers,file="_1_Re_Post_markers.csv")
write.csv(NoR_Pre_Data_1.markers,file="_1_NoR_Pre_markers.csv")
write.csv(NoR_Post_Data_1.markers,file="_1_NoR_Post_markers.csv")

#View top 20 genes
library(dplyr) 
Re_All_1_top40genes <- Re_All_Data_1.markers %>% group_by(cluster) %>% top_n(n=40, wt=avg_log2FC)
NoR_All_1_top20genes <- NoR_All_Data_1.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
#select_genes1 <- NoR_All_1_top20genes $gene [1:20]
#p1 <- VlnPlot (NoR_All_Data_1     , features =select_genes1  , pt.size=0, ncol=3)
#ggsave("196_test_VlnPlot.png", p1, width=30 ,height=20)
Re_Pre_1_top20genes <- Re_Pre_Data_1.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
Re_Post_1_top20genes <- Re_Post_Data_1.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
NoR_Pre_1_top20genes <- NoR_Pre_Data_1.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
NoR_Post_1_top20genes <- NoR_Post_Data_1.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
select_genes_1 <- Re_All_1_top40genes $gene [c(30,41,83,126,161,201,241,281,321)]
select_genes_2 <- NoR_All_1_top20genes $gene [c(1,23,41,61,81,101,121,141,161)]
select_genes_3 <- Re_Pre_1_top20genes $ gene[c(7,21,41,63,81,104,121,141)]
#select_genes_1 <- Re_Pre_1_top20genes $gene [21:40]
#p1 <- VlnPlot (Re_Pre_Data_1 , features =select_genes_1  , pt.size=0, ncol=3)
#ggsave("205_test_VlnPlot.png", p1, width=30 ,height=30)
select_genes_4 <- Re_Post_1_top20genes $gene[c(3,21,41,62,82,101,121,141,161)]
#select_genes_1 <- Re_Post_1_top20genes $gene [1:20]
#p1 <- VlnPlot (Re_Post_Data_1 , features =select_genes_1  , pt.size=0, ncol=3)
#ggsave("213_test_VlnPlot.png", p1, width=30 ,height=30)
select_genes_5 <- NoR_Pre_1_top20genes $gene[c(2,22,42,61,81,101,121,141,161)]
#select_genes_1 <- NoR_Pre_1_top20genes $gene [21:40]
#p1 <- VlnPlot (NoR_Pre_Data_1 , features =select_genes_1  , pt.size=0, ncol=3)
#ggsave("222_test_VlnPlot.png", p1, width=30 ,height=30)
select_genes_6 <- NoR_Post_1_top20genes $gene[c(1,22,41,61,81,102,121,141,161)]
#select_genes_1 <- NoR_Post_1_top20genes $gene [41:60]
#p1 <- VlnPlot (NoR_Post_Data_1 , features =select_genes_1  , pt.size=0, ncol=3)
#ggsave("231_test_VlnPlot.png", p1, width=30 ,height=30)


#featureplot display
p1 <- FeaturePlot(Re_All_Data_1     , features =select_genes_1  , reduction = "tsne", label=T, ncol=3)
p2 <- FeaturePlot(NoR_All_Data_1   , features =select_genes_2  , reduction = "tsne", label=T, ncol=3)
p3 <- FeaturePlot(Re_Pre_Data_1     , features =select_genes_3  , reduction = "tsne", label=T, ncol=3)
p4 <- FeaturePlot(Re_Post_Data_1   , features =select_genes_4  , reduction = "tsne", label=T, ncol=3)
p5 <- FeaturePlot(NoR_Pre_Data_1  , features =select_genes_5  , reduction = "tsne", label=T, ncol=3)
p6 <- FeaturePlot(NoR_Post_Data_1, features =select_genes_6  , reduction = "tsne", label=T, ncol=3)

ggsave("240_G1_FeaturePlot.png", p1, width=12 ,height=10)
ggsave("241_G2_FeaturePlot.png", p2, width=12 ,height=10)
ggsave("242_G3_FeaturePlot.png", p3, width=12 ,height=10)
ggsave("243_G4_FeaturePlot.png", p4, width=12 ,height=10)
ggsave("244_G5_FeaturePlot.png", p5, width=12 ,height=10)
ggsave("245_G6_FeaturePlot.png", p6, width=12 ,height=10)
  
  
  
p1 <- VlnPlot (Re_All_Data_1     , features =select_genes_1  , pt.size=0, ncol=3)
p2 <- VlnPlot (NoR_All_Data_1    , features =select_genes_2  , pt.size=0, ncol=3)
p3 <- VlnPlot (Re_Pre_Data_1     , features =select_genes_3 , pt.size=0, ncol=3)
p4 <- VlnPlot (Re_Post_Data_1     , features =select_genes_4  , pt.size=0, ncol=3)
p5 <- VlnPlot (NoR_Pre_Data_1     , features =select_genes_5  , pt.size=0, ncol=3)
p6 <- VlnPlot (NoR_Post_Data_1     , features =select_genes_6  , pt.size=0, ncol=3)
ggsave("246_G1_VlnPlot.png", p1, width=12 ,height=10)
ggsave("247_G2_VlnPlot.png", p2, width=12 ,height=10)
ggsave("248_G3_VlnPlot.png", p3, width=12 ,height=10)
ggsave("249_G4_VlnPlot.png", p4, width=12 ,height=10)
ggsave("250_G5_VlnPlot.png", p5, width=12 ,height=10)
ggsave("251_G6_VlnPlot.png", p6, width=12 ,height=10)
  
  
  
save(Re_All_Data_1      ,file='Re_All_Data_1.Rdata')
save(NoR_All_Data_1   ,file='NoR_All_Data_1.Rdata')
save(Re_Pre_Data_1     ,file='Re_Pre_Data_1.Rdata')
save(Re_Post_Data_1    ,file='Re_Post_Data_1.Rdata')
save(NoR_Pre_Data_1  ,file='NoR_Pre_Data_1.Rdata')
save(NoR_Post_Data_1,file='NoR_Post_Data_1.Rdata')

Step 4: Compare distribution differences of responder samples in each group
group_1<- Re_All_Data_1@meta.data
group_2<- NoR_All_Data_1@meta.data
group_3<- Re_Pre_Data_1@meta.data
group_4<- Re_Post_Data_1@meta.data
group_5<- NoR_Pre_Data_1@meta.data
group_6<- NoR_Post_Data_1@meta.data

Re_All_Sample_Cell_TotalNum <- as.data.frame(table(group_1$V5))
NoR_All_Sample_Cell_TotalNum <- as.data.frame(table(group_2$V5))
Re_Pre_Sample_Cell_TotalNum <- as.data.frame(table(group_3$V5))
Re_Post_Sample_Cell_TotalNum <- as.data.frame(table(group_4$V5))
NoR_Pre_Sample_Cell_TotalNum <- as.data.frame(table(group_5$V5))
NoR_Post_Sample_Cell_TotalNum <- as.data.frame(table(group_6$V5))
Re_All_Sample_Cell_TotalNum 
NoR_All_Sample_Cell_TotalNum
     

Re_Pre_Sample_Cell_TotalNum
Re_Post_Sample_Cell_TotalNum
   

NoR_Pre_Sample_Cell_TotalNum 
NoR_Post_Sample_Cell_TotalNum
     

PerSample_Cell_Proportion <- function(group,Sample_Cell_TotalNum){
	#group refers to the analyzed sample group (responder or non-responder...)
	#Sample_Cell_TotalNum refers to the dataframe of cell counts per sample in this group
	Cell_Proportion <- list()
	for(i in 1:dim(Sample_Cell_TotalNum)[1]){
		one_sample <- group[group$V5==Sample_Cell_TotalNum$Var1[i],]
		all_cellTypes <- as.data.frame(table(one_sample$cellTypes))
		rownames(all_cellTypes) <- all_cellTypes$Var1
		Cell_Proportion[[i]]<- (all_cellTypes)/Sample_Cell_TotalNum$Freq[i]
		Cell_Proportion[[i]]$Var1<- all_cellTypes$Var1
	}
	Cell_Proportion
}

Re_All_Cell_Proportion    <- PerSample_Cell_Proportion(group_1,Re_All_Sample_Cell_TotalNum)
NoR_All_Cell_Proportion <- PerSample_Cell_Proportion(group_2,NoR_All_Sample_Cell_TotalNum)
Re_Pre_Cell_Proportion   <- PerSample_Cell_Proportion(group_3, Re_Pre_Sample_Cell_TotalNum)
Re_Post_Cell_Proportion  <- PerSample_Cell_Proportion(group_4, Re_Post_Sample_Cell_TotalNum)
NoR_Pre_Cell_Proportion<-PerSample_Cell_Proportion(group_5, NoR_Pre_Sample_Cell_TotalNum)
NoR_Post_Cell_Proportion<-PerSample_Cell_Proportion(group_6,NoR_Post_Sample_Cell_TotalNum)
Re_All_Cell_Proportion   
   
NoR_All_Cell_Proportion
  
Re_Pre_Cell_Proportion  
   
Re_Post_Cell_Proportion 
   
NoR_Pre_Cell_Proportion
   
NoR_Post_Cell_Proportion
   


cell_Pro_dataframe <- function(Cell_Proportion){
	B_Cell_Pro<- c()
	CD4_T_Cell_Pro<- c()
	CD8_T_Cell_Pro <- c()
	Monocyte_Pro <- c()
	NKT_Cell_Pro <- c()
	pDCs_Pro <- c()
	Plasma_Cell_Pro <- c()
	Regulatory_T_Cell_Pro <- c()
	T_Cell_Pro <- c()

	for(i in 1:length(Cell_Proportion)){
		if(sum(Cell_Proportion[[i]]$Var1=='B cells')!=0){
			B_Cell_Pro[i] <- Cell_Proportion[[i]][Cell_Proportion[[i]]$Var1=='B cells',]$Freq
		}else B_Cell_Pro[i]=0
		
		if(sum(Cell_Proportion[[i]]$Var1=='CD4 T cells')!=0){
			CD4_T_Cell_Pro[i] <- Cell_Proportion[[i]][Cell_Proportion[[i]]$Var1=='CD4 T cells',]$Freq
		}else CD4_T_Cell_Pro[i]=0

		if(sum(Cell_Proportion[[i]]$Var1=='CD8 T cells')!=0){
			CD8_T_Cell_Pro[i] <- Cell_Proportion[[i]][Cell_Proportion[[i]]$Var1=='CD8 T cells',]$Freq
		}else CD8_T_Cell_Pro[i]=0
		
		if(sum(Cell_Proportion[[i]]$Var1=='Monocyte')!=0){
			Monocyte_Pro[i] <- Cell_Proportion[[i]][Cell_Proportion[[i]]$Var1=='Monocyte',]$Freq
		}else Monocyte_Pro[i]=0
		
		if(sum(Cell_Proportion[[i]]$Var1=='NKT cells')!=0){
			NKT_Cell_Pro[i] <- Cell_Proportion[[i]][Cell_Proportion[[i]]$Var1=='NKT cells',]$Freq
		}else NKT_Cell_Pro[i]=0
		
		if(sum(Cell_Proportion[[i]]$Var1=='pDCs')!=0){
			pDCs_Pro[i] <- Cell_Proportion[[i]][Cell_Proportion[[i]]$Var1=='pDCs',]$Freq
		}else pDCs_Pro[i]=0
		
		if(sum(Cell_Proportion[[i]]$Var1=='Plasma cells')!=0){
			Plasma_Cell_Pro[i] <- Cell_Proportion[[i]][Cell_Proportion[[i]]$Var1=='Plasma cells',]$Freq
		}else Plasma_Cell_Pro[i]=0
		
		if(sum(Cell_Proportion[[i]]$Var1=='Regulatory T cells')!=0){
			Regulatory_T_Cell_Pro[i] <- Cell_Proportion[[i]][Cell_Proportion[[i]]$Var1=='Regulatory T cells',]$Freq
		}else Regulatory_T_Cell_Pro[i]=0
		
		if(sum(Cell_Proportion[[i]]$Var1=='T cells')!=0){
			T_Cell_Pro[i] <- Cell_Proportion[[i]][Cell_Proportion[[i]]$Var1=='T cells',]$Freq
		}else T_Cell_Pro[i]=0
		
	}
	Proportion_List=list(B_Cell_Pro,CD4_T_Cell_Pro,CD8_T_Cell_Pro,Monocyte_Pro,NKT_Cell_Pro,pDCs_Pro,Plasma_Cell_Pro,Regulatory_T_Cell_Pro,T_Cell_Pro)
	names(Proportion_List) <-c('B cells','CD4 T cells','CD8 T cells','Monocyte','NKT cells','pDCs','Plasma cells','Regulatory T cells','T cells')
	Proportion_List
}
#B cells
#CD4 T cells
#CD8 T cells
#Monocyte
#NKT cells
#pDCs
#Plasma cells
#Regulatory T cells
#T cells

Re_All_cell_Prop_list<-cell_Pro_dataframe(Re_All_Cell_Proportion)
Re_All_cell_Prop_list
   
NoR_All_cell_Prop_list<-cell_Pro_dataframe(NoR_All_Cell_Proportion) 
NoR_All_cell_Prop_list
 

 
Re_Pre_cell_Prop_list<-cell_Pro_dataframe(Re_Pre_Cell_Proportion)
Re_Pre_cell_Prop_list
 
Re_Post_cell_Prop_list<-cell_Pro_dataframe(Re_Post_Cell_Proportion)
 
NoR_Pre_cell_Prop_list<-cell_Pro_dataframe(NoR_Pre_Cell_Proportion) 
NoR_Pre_cell_Prop_list
 

NoR_Post_cell_Prop_list<-cell_Pro_dataframe(NoR_Post_Cell_Proportion) 
NoR_Post_cell_Prop_list
 

 
 
Cell proportion situation				
## Calculate cell counts for each cluster
# Define cell types corresponding to each cluster, can change ident to corresponding cell types, calculate cell counts for each cell type
table(Idents(Re_All_Data_1))
 
# Calculate cell proportions
prop.table(table(Idents(Re_All_Data_1)))
 
## Draw stacked bar chart
library("ggplot2")
cell.prop_1<-data.frame(as.data.frame(prop.table(table(Idents(Re_All_Data_1)))),Tag='R')
#R(Responder)
cell.prop_2<-data.frame(as.data.frame(prop.table(table(Idents(NoR_All_Data_1)))),Tag='NR')
#NR(Non-Responder)
cell.prop_3<-data.frame(as.data.frame(prop.table(table(Idents(Re_Pre_Data_1)))),Tag='RB')
#RB(Responder-Pre)
cell.prop_4<-data.frame(as.data.frame(prop.table(table(Idents(Re_Post_Data_1)))),Tag='RP')
#RP(Responder-Post)
cell.prop_5<-data.frame(as.data.frame(prop.table(table(Idents(NoR_Pre_Data_1)))),Tag='NRB')
#NRB(Non-Responder-Pre)
cell.prop_6<-data.frame(as.data.frame(prop.table(table(Idents(NoR_Post_Data_1)))),Tag='NRP')
#NRP(Non-Responder-Post)
write.csv(cell.prop_1,file='_1Re_All_cell_proportion-stacked_chart.csv')
write.csv(cell.prop_2,file='_1NoR_All_cell_proportion-stacked_chart.csv')
write.csv(cell.prop_3,file='_1Re_Pre_cell_proportion-stacked_chart.csv')
write.csv(cell.prop_4,file='_1Re_Post_cell_proportion-stacked_chart.csv')
write.csv(cell.prop_5,file='_1NoR_Pre_cell_proportion-stacked_chart.csv')
write.csv(cell.prop_6,file='_1NoR_Post_cell_proportion-stacked_chart.csv')

colnames(cell.prop_1)[1:2]<-c("cluster","proportion")
colnames(cell.prop_2)[1:2]<-c("cluster","proportion")
colnames(cell.prop_3)[1:2]<-c("cluster","proportion")
colnames(cell.prop_4)[1:2]<-c("cluster","proportion")
colnames(cell.prop_5)[1:2]<-c("cluster","proportion")
colnames(cell.prop_6)[1:2]<-c("cluster","proportion")



#First group comparison shows R vs NR
cell.prop_compare1 <- rbind(cell.prop_1,cell.prop_2)
#Second group comparison shows R-Pre vs NR-Pre
cell.prop_compare2 <- rbind(cell.prop_3,cell.prop_5)
#Third group comparison shows R-Post vs NR-Post
cell.prop_compare3 <- rbind(cell.prop_4,cell.prop_6)
#Fourth group comparison shows R-Pre vs R-Post
cell.prop_compare4 <- rbind(cell.prop_3,cell.prop_4)
#Fifth group comparison shows NR-Pre vs NR-Post

#Stacked bar chart subpopulation fill colors
Bar_color <- function(cell.prop){
	col_res <-c()
	for(i in 1:length(cell.prop)){
		ifelse(cell.prop[i]=="CD8 T cells",col_res[i]<- Mycol_5[3],
		ifelse(cell.prop[i]=="CD5 T cells",
		col_res[i]<- Mycol_5[2],
		ifelse(cell.prop[i]=="Regulatory T cells",col_res[i]<- Mycol_5[6], 
		ifelse(cell.prop[i]=="Monocyte",col_res[i]<- Mycol_5[4], 
		ifelse(cell.prop[i]=="NKT cells",col_res[i]<- Mycol_5[1],
		ifelse(cell.prop[i]=="B cells", col_res[i]<- Mycol_5[7],
		ifelse(cell.prop[i]=="T cells", col_res[i]<- Mycol_5[9],
		ifelse(cell.prop[i]=="pDCs", col_res[i]<- Mycol_5[5],
		col_res[i]<- Mycol_5[8]))))))))
	}
	col_res
}
cell_prop_compare1 <- as.vector(unique(cell.prop_compare1$cluster))
cell_prop_compare2 <- as.vector(unique(cell.prop_compare2$cluster))
cell_prop_compare3 <- as.vector(unique(cell.prop_compare3$cluster))
cell_prop_compare4 <- as.vector(unique(cell.prop_compare4$cluster)) 
cell_prop_compare5 <- as.vector(unique(cell.prop_compare5$cluster)) 
bar_color1 <- Bar_color(cell_prop_compare1)
bar_color2 <- Bar_color(cell_prop_compare2)
bar_color3 <- Bar_color(cell_prop_compare3)
bar_color4 <- Bar_color(cell_prop_compare4)
bar_color5 <- Bar_color(cell_prop_compare5)



p_1 <- ggplot(cell.prop_compare1,aes(Tag,proportion,fill=cluster))+
geom_bar(stat="identity",position="fill")+labs(x='')+
ggtitle("")+
theme_bw()+
theme(axis.ticks.length=unit(0.1,'cm'))+
scale_fill_manual(values=bar_color1)+
guides(fill=guide_legend(title=NULL))

p_2 <- ggplot(cell.prop_compare2,aes(Tag,proportion,fill=cluster))+
geom_bar(stat="identity",position="fill")+labs(x='')+
ggtitle("")+
theme_bw()+
theme(axis.ticks.length=unit(0.1,'cm'))+
scale_fill_manual(values=bar_color2)+
guides(fill=guide_legend(title=NULL))

p_3 <- ggplot(cell.prop_compare3,aes(Tag,proportion,fill=cluster))+
geom_bar(stat="identity",position="fill")+labs(x='')+
ggtitle("")+
theme_bw()+
theme(axis.ticks.length=unit(0.1,'cm'))+
scale_fill_manual(values=bar_color3)+
guides(fill=guide_legend(title=NULL))

p_4 <- ggplot(cell.prop_compare4,aes(Tag,proportion,fill=cluster))+
geom_bar(stat="identity",position="fill")+labs(x='')+
ggtitle("")+
theme_bw()+
theme(axis.ticks.length=unit(0.1,'cm'))+
scale_fill_manual(values=bar_color4)+
guides(fill=guide_legend(title=NULL))

p_5 <- ggplot(cell.prop_compare5,aes(Tag,proportion,fill=cluster))+
geom_bar(stat="identity",position="fill")+labs(x='')+
ggtitle("")+
theme_bw()+
theme(axis.ticks.length=unit(0.1,'cm'))+
scale_fill_manual(values=bar_color5)+
guides(fill=guide_legend(title=NULL))
 
ggsave("252_RvsNR_cell_proportion_chart.png", p_1, width=4 ,height=6)
ggsave("253_RBvsNRB_cell_proportion_chart.png", p_2, width=4 ,height=6)
ggsave("254_RPvsNRP_cell_proportion_chart.png", p_3, width=4 ,height=6)
ggsave("255_RBvsRP_cell_proportion_chart.png", p_4, width=4 ,height=6)
ggsave("256_NRBvsNRP_cell_proportion_chart.png", p_5, width=4 ,height=6)
     

library(ggpubr)
change_to_dataframe <- function(cell_Prop_list,TagMode){
	data_frame <- data.frame()
	for (i in 1:9){ 
	#A total of 9 cell types summarized
		Prop <- cell_Prop_list[[i]]
		cellTp <- names(cell_Prop_list[i])
		cellTypes <-rep(cellTp,length(Prop))
		addData <- cbind(Prop,cellTypes)
		data_frame <- rbind(data_frame,addData)
	}
	Tag <- rep(TagMode,length(cell_Prop_list[[i]])*9)
	data_frame <- data.frame(data_frame,Tag)
}

Re_All_cell_Prop_df <-change_to_dataframe(Re_All_cell_Prop_list,"R")
NoR_All_cell_Prop_df <-change_to_dataframe(NoR_All_cell_Prop_list,"NR")
Re_Pre_cell_Prop_df <-change_to_dataframe(Re_Pre_cell_Prop_list,"RB")
Re_Post_cell_Prop_df <-change_to_dataframe(Re_Post_cell_Prop_list,"RP")
NoR_Pre_cell_Prop_df <-change_to_dataframe(NoR_Pre_cell_Prop_list,"NRB")
NoR_Post_cell_Prop_df <-change_to_dataframe(NoR_Post_cell_Prop_list,"NRP")
comp1 <- rbind(Re_All_cell_Prop_df,NoR_All_cell_Prop_df)
comp2 <- rbind(Re_Pre_cell_Prop_df,NoR_Pre_cell_Prop_df)
comp3 <- rbind(Re_Post_cell_Prop_df,NoR_Post_cell_Prop_df)
comp4 <- rbind(Re_Pre_cell_Prop_df,Re_Post_cell_Prop_df)
comp5 <- rbind(NoR_Pre_cell_Prop_df,NoR_Post_cell_Prop_df)

bar_Prop <- function(compData,compNum,col_order){
	compData <- data.frame(as.numeric(compData$Prop),compData$cellTypes,compData$Tag) 
	colnames(compData)<- c("Prop","cellTypes","Tag")
	p1 <- ggbarplot(compData, x = "cellTypes", y = "Prop", add = c("mean_se", "jitter"),
#shape= "Tag",
color = "Tag", palette = col_order, position = position_dodge(0.8))
	p1 <- p1 + stat_compare_means(aes(group = Tag), label = "p.format") 
	# label = "p.format" indicates using P values to indicate significance results, as shown below
	i0=256
	filename <- paste(i0+compNum,"Comp",compNum,"barplot.png",sep="_")
	ggsave(filename,width = 15,height = 8)
}
bar_Prop(comp1,1, rev(Mycol_2))
bar_Prop(comp2,2, rev(Mycol_2))
bar_Prop(comp3,3, rev(Mycol_2))
bar_Prop(comp4,4, Mycol_2)
bar_Prop(comp5,5, Mycol_2)
