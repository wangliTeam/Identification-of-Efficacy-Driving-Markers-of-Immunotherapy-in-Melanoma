save.image()

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
# Color palette
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175')
Mycol_2 <-c('#FF7F0EFF', '#1F77B4FF')
Mycol_3 <- c("#1B9E77","#D95F02","#7570B3")
#library("ggsci")
#Mycol_nejm= pal_nejm("default", alpha =0.7)(8)## Extract colors
Mycol_4 <- c("#A6CEE3" ,"#1F78B4", "#B2DF8A" ,"#33A02C" ,"#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00", "#CAB2D6" ,"#6A3D9A" ,"#FFFF99", "#B15928","#1B9E77", "#D95F02" ,"#7570B3", "#E7298A" ,"#66A61E" ,"#E6AB02", "#A6761D","#666666")

# Step 1: Data processing and clustering

# Import count matrix - sc_data
# Import metadata matrix - metadata
sc_data = read.table('GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt',header=T,row.names=1,skip = 1)
colnames(sc_data)[1:100] 
 
metadata = read.table("GSE120575_patient_ID_single_cells.txt",row.name=1, header=F,sep="\t")
metadata <- metadata[,1:6]
metadata <- metadata[-dim(metadata)[1],]
head(metadata)
 
dim(sc_data)
dim(metadata)
 
# Modify column names of sc_data
colnames(sc_data)=rownames(metadata)
# Confirm that row names of metadata and column names of sc_data are equal
colnames(sc_data)[1:100]
 
### save.image()

# Extract fifth and sixth columns of metadata
PP<- as.data.frame(str_split_fixed(metadata$V5, "_", 2))
Patient_information<-as.data.frame(cbind(PP , metadata$V6,metadata$V7))
# Label rows and columns
rownames(Patient_information)=rownames(metadata)
colnames(Patient_information) = c("Pre_or_Post","Patient","Re_or_NoR","Medicine")
head(Patient_information)
 
attach(Patient_information)
q <- aggregate(x= Patient, by = list(Medicine, Re_or_NoR),FUN=length)
  
C1<-table(Patient_information[(Medicine=="anti-CTLA4")&( Re_or_NoR=="Non-responder"),][,1:2])
C2<-table(Patient_information[(Medicine=="anti-CTLA4")&( Re_or_NoR=="Responder"),][,1:2])
C3<-table(Patient_information[(Medicine=="anti-CTLA4+PD1")&( Re_or_NoR=="Non-responder"),][,1:2])
C4<-table(Patient_information[(Medicine=="anti-CTLA4+PD1")&( Re_or_NoR=="Responder"),][,1:2])
C5<-table(Patient_information[(Medicine=="anti-PD1")&(Re_or_NoR=="Non-responder"),][,1:2])
C6<-table(Patient_information[(Medicine=="anti-PD1")&( Re_or_NoR=="Responder"),][,1:2])
### save.image()


# Classify cell samples based on above four categories - only take PD-1 inhibitor treated samples
attach( Patient_information)
ifPD_1 <- Medicine=="anti-PD1"
Patient_information <- Patient_information[ifPD_1,]
sc_data <- sc_data[,ifPD_1]
metadata <- metadata[ifPD_1,]
dim(sc_data)
dim(metadata)
 
############ Category 1: Anti-PD-1 treatment, responder patients, pre-treatment cells ####################
Re_Pre <- which( (Re_or_NoR=="Responder") & (Pre_or_Post=="Pre"))
sc_data1 <- sc_data[,Re_Pre]
metadata1 <- metadata[Re_Pre,]
metadata1$V5
dim(sc_data1)
dim(metadata1)
 
table(metadata1$V5)
############ Category 2: Anti-PD-1 treatment, responder patients, post-treatment cells ####################
Re_Post <- which( (Re_or_NoR=="Responder") & (Pre_or_Post=="Post"))
sc_data2 <- sc_data[,Re_Post]
metadata2 <- metadata[Re_Post,]
dim(sc_data2)
dim(metadata2)
 
table(metadata2$V5)
 
############ Category 3: Anti-PD-1 treatment, non-responder patients, pre-treatment cells ####################
NoR_Pre <- which( (Re_or_NoR=="Non-responder") & (Pre_or_Post=="Pre"))
sc_data3 <- sc_data[,NoR_Pre]
metadata3<- metadata[NoR_Pre,]
dim(sc_data3)
dim(metadata3)
 
table(metadata3$V5)

############ Category 4: Anti-PD-1 treatment, non-responder patients, post-treatment cells ####################
NoR_Post <- which( (Re_or_NoR=="Non-responder") & (Pre_or_Post=="Post"))
sc_data4 <- sc_data[,NoR_Post]
metadata4 <- metadata[NoR_Post,]
dim(sc_data4)
dim(metadata4)
 
table(metadata4$V5)
 
length(Re_Pre)
length(Re_Post)
length(NoR_Pre)
length(NoR_Post)
   # Results consistent with previous classification
# Need to filter ERCC (External RNA Control Consortium) here
# ERCC is exogenous RNA, mainly used for quality control
erccs = grep('^ERCC', x= rownames(x=sc_data),value = T) # value = T to get names
# Total of 10 (filter all four groups together)
# Calculate ERCC percentage
percent.ercc1 = Matrix::colSums(sc_data1[erccs, ])/Matrix::colSums(sc_data1)
percent.ercc2 = Matrix::colSums(sc_data2[erccs, ])/Matrix::colSums(sc_data2)
percent.ercc3 = Matrix::colSums(sc_data3[erccs, ])/Matrix::colSums(sc_data3)
percent.ercc4 = Matrix::colSums(sc_data4[erccs, ])/Matrix::colSums(sc_data4)

fivenum(percent.ercc1) 
fivenum(percent.ercc2)
fivenum(percent.ercc3)
fivenum(percent.ercc4) 
 

# Get indices
ercc.index = grep(pattern = "^ERCC", x = rownames(x = sc_data), value = FALSE) 
sc_data1 = sc_data1[-ercc.index,]
sc_data2 = sc_data2[-ercc.index,]
sc_data3 = sc_data3[-ercc.index,]
sc_data4 = sc_data4[-ercc.index,]

dim(sc_data1) 
dim(sc_data2) 
dim(sc_data3) 
dim(sc_data4)
 



# Create objects
# min.cells: minimum number of cells a gene must be expressed in
# min.features: minimum number of genes a cell must express
min.cells = 0  # Initial value 0
min.features = 1000   # Same as above

######## Group 1 ##########################
sce1 = CreateSeuratObject(counts =sc_data1,metadata = metadata1,min.cells =min.cells,min.features =min.features)
sce1
 
# meta.data in sce object is auto-generated when creating Seurat object, missing a lot of info, so need to add metadata
sce1 = AddMetaData(object = sce1, metadata = metadata1)
sce1 = AddMetaData(object = sce1, percent.ercc1, col.name = "percent.ercc")
# Check
head(sce1@meta.data)
 
######## Group 2  ##########################
sce2 = CreateSeuratObject(counts =sc_data2,metadata = metadata2,min.cells =min.cells,min.features =min.features)
sce2 = AddMetaData(object = sce2, metadata = metadata2)
sce2 = AddMetaData(object = sce2, percent.ercc2, col.name = "percent.ercc")
head(sce2@meta.data)


######## Group 3  ##########################
sce3 = CreateSeuratObject(counts =sc_data3,metadata = metadata3,min.cells =min.cells,min.features =min.features)
sce3 = AddMetaData(object = sce3, metadata = metadata3)
sce3 = AddMetaData(object = sce3, percent.ercc3, col.name = "percent.ercc")
head(sce3@meta.data)
 
######## Group 4  ##########################
sce4 = CreateSeuratObject(counts =sc_data4,metadata = metadata4,min.cells =min.cells,min.features =min.features)
sce4 = AddMetaData(object = sce4, metadata = metadata4)
sce4 = AddMetaData(object = sce4, percent.ercc4, col.name = "percent.ercc")
head(sce4@meta.data)
 



# Data cleaning
# Data cleaning principles: 1. Low quality; 2. Filter mitochondrial DNA; 3. Filter exogenous DNA;
# Calculate before filtering
######## Group 1  ##########################
table(grepl("^MT-",rownames(sce1)))   
table(grepl("^RP[SL][[:digit:]]",rownames(sce1)))
 

sce1[["percent.mt"]]  = PercentageFeatureSet(sce1, pattern = "^MT-")  # Add attribute column
sce1[["percent.rp"]]  = PercentageFeatureSet(sce1, pattern = "^RP[SL][[:digit:]]")
summary(sce1[["nCount_RNA"]])     # nCount_RNA is number of UMIs per cell
summary(sce1[["nFeature_RNA"]])     # nFeature_RNA is number of genes detected per cell
summary(sce1[["percent.mt"]]  )
summary(sce1[["percent.ercc"]] )
summary(sce1[["percent.rp"]]  )
 





######## Group 3  ##########################
sce2[["percent.mt"]]  = PercentageFeatureSet(sce2, pattern = "^MT-") 
sce2[["percent.rp"]]  = PercentageFeatureSet(sce2, pattern = "^RP[SL][[:digit:]]")
summary(sce2[["nCount_RNA"]])    
summary(sce2[["nFeature_RNA"]])    
summary(sce2[["percent.mt"]]  )
summary(sce2[["percent.ercc"]] )
summary(sce2[["percent.rp"]]  )
   

######## Group 3  ##########################
sce3[["percent.mt"]]  = PercentageFeatureSet(sce3, pattern = "^MT-") 
sce3[["percent.rp"]]  = PercentageFeatureSet(sce3, pattern = "^RP[SL][[:digit:]]")
summary(sce3[["nCount_RNA"]])    
summary(sce3[["nFeature_RNA"]])    
summary(sce3[["percent.mt"]]  )
summary(sce3[["percent.ercc"]] )
summary(sce3[["percent.rp"]]  )
 

######## Group 4  ##########################
sce4[["percent.mt"]]  = PercentageFeatureSet(sce4, pattern = "^MT-") 
sce4[["percent.rp"]]  = PercentageFeatureSet(sce4, pattern = "^RP[SL][[:digit:]]")
summary(sce4[["nCount_RNA"]])    
summary(sce4[["nFeature_RNA"]])    
summary(sce4[["percent.mt"]]  )
summary(sce4[["percent.ercc"]] )
summary(sce4[["percent.rp"]]  )
   



# Plot to check
violin1 <- VlnPlot(sce1,features = c("nCount_RNA","nFeature_RNA", "percent.mt","percent.rp"), cols =rainbow(4), pt.size = 0.01, ncol = 4)
violin2 <- VlnPlot(sce2,features = c("nCount_RNA","nFeature_RNA", "percent.mt","percent.rp"), cols =rainbow(1), pt.size = 0.01, ncol = 4)
violin3 <- VlnPlot(sce3,features = c("nCount_RNA","nFeature_RNA", "percent.mt","percent.rp"), cols =rainbow(4), pt.size = 0.01, ncol = 4)
violin4 <- VlnPlot(sce4,features = c("nCount_RNA","nFeature_RNA", "percent.mt","percent.rp"), cols =rainbow(4), pt.size = 0.01, ncol = 4)
ggsave("1_Re_Pre_Features.png", plot = violin1, width = 12, height = 6)
ggsave("2_Re_Post_Features.png", plot = violin2, width = 12, height = 6)
ggsave("3_NoR_Pre_Features.png", plot = violin3, width = 12, height = 6)
ggsave("4_NoR_Post_Features.png", plot = violin4, width = 12, height = 6)
  
  
All_Data <- Reduce(merge,list(sce1,sce2,sce3,sce4))  
# All_Data integrates previous four categories to form new Seurat object
All_Meta = All_Data@meta.data
Tag <- rep(c(1,2,3,4),times=c(1191,1524,2604,6334))
All_Data = AddMetaData(object = All_Data, metadata =All_Meta)
All_Data = AddMetaData(object = All_Data, Tag, col.name = "Tag")

Violin_Features_All <- VlnPlot(All_Data, features = c("nCount_RNA","nFeature_RNA","percent.mt", "percent.ercc","percent.rp") , group.by = "Tag",ncol=2,cols=my36colors)
ggsave("5_Violin_Features_All_2.png", plot = Violin_Features_All, width = 10, height = 10)
 
save(sce1,file='sce1.Rdata')
save(sce2,file='sce2.Rdata')
save(sce3,file='sce3.Rdata')
save(sce4,file='sce4.Rdata')
save(All_Data,file='All_Data.Rdata')

## First, perform clustering on the overall All_Data object
## 1. Start plotting
## Workflow: 1. log: NormalizeData              2. Find features: FindVariableFeatures 
##       3. Scale: ScaleData               4. Dimensionality reduction: ## RunTSNE
##       5. Build graph: FindNeighbors          6. Cluster: FindClusters 
##       7. tsne /umap: RunTSNE RunUMAP 8. ## Differential genes: FindAllMarkers / FindMarkers
##### (1) Take log ###################################
New_All_Data= All_Data
New_All_Data = NormalizeData(object = New_All_Data,normalization.method =  "LogNormalize",  scale.factor = 1e6)
# NormalizeData(), divides each gene's count by total, multiplies by scale.factor, then transforms with natural logarithm
# Aims to eliminate effects of different cell sequencing depths, log reduces data dispersion, +1 prevents read=0 case


# 2. Feature extraction
New_All_Data = FindVariableFeatures(object = New_All_Data,selection.method = "vst", nfeatures = 4000)
# FindVariableFeatures is hard filtering, based on statistical metrics like sd, mad, vst etc. to judge the most important 4000 genes among the 55k+ genes in input scRNA expression matrix, remaining 52k genes are not considered in downstream analysis
# Check what the top 4000 genes are (top six)
head(VariableFeatures(New_All_Data))
 

top10=head(VariableFeatures(New_All_Data),10)        # Select top 10 genes with highest expression dispersion
plot1= VariableFeaturePlot(New_All_Data)
plot2=LabelPoints(plot = plot1, points = top10 , repel = TRUE)+ FontSize(x.title = 20, y.title = 20) 
ggsave("6_VariableFeatures.png", plot = plot2, width = 15, height = 10)
 
# 3. Scaling
New_All_Data = ScaleData(object = New_All_Data)

# 4. PCA dimensionality reduction
New_All_Data = RunPCA(object = New_All_Data, do.print = FALSE)
# 5000 genes will be transformed into 5000 dimensions, but we usually just look at first dozen or so dimensions, so also very efficient dimensionality reduction method
# Originally each cell has 5000 genes, each gene has expression in each cell, after PCA analysis each cell has 50 PCs, each PC has embedding value in each cell. For easier understanding, you can think of PC as meta-gene integrating multiple genes, cell embedding as expression of these meta-genes
# RunPCA function calculates 50 PCs by default, if want to see more PCs can modify npcs parameter of RunPCA.

# Visualization
PCA_Elbow_Image <- ElbowPlot(New_All_Data, ndims=25)
# After ElbowPlot, observe with naked eye, when curve enters plateau phase with no obvious continued downward trend is okay. ElbowPlot shows 20 PCs by default, if find 20 PCs haven't reached plateau can modify ndims parameter of ElbowPlot function
ggsave("7_PCA_Elbow_Image.png", plot = PCA_Elbow_Image, width = 10, height = 6)





VizDim<-VizDimLoadings(New_All_Data, dims = 1:6, reduction = "pca", nfeatures = 10,ncol=2)
ggsave("8_VizDim_Image.png", plot = VizDim, width = 15, height = 15)
# Visualize top genes associated with reduction components
 
# Principal component analysis plot
Pca_Image <- DimPlot(object = New_All_Data, reduction = "pca")
ggsave("9_PCA_Image.png", plot = Pca_Image, width = 10, height = 10)
 

## TSNE clustering analysis
pcSelect=15
New_All_Data<- FindNeighbors(object = New_All_Data, dims = 1:pcSelect)               
# Calculate adjacency distance
save(New_All_Data,file='New_All_Data.Rdata')
###### Cluster choosing ############
# Check clustering stability at given resolution  
# Set different resolutions 
res.used <- seq(0.1,0.7,by=0.1)
res.used
# Loop over and perform clustering of different resolutions 
New_All_Data_2 <- New_All_Data
for(i in res.used){
  New_All_Data_2 <- FindClusters(object = New_All_Data_2, verbose = T, resolution = res.used)
}
# Make plot 
library(clustree)
clus.tree.out <- clustree(New_All_Data_2) +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
ggsave("10-_clustering_Image.png", plot = clus.tree.out, width = 14, height = 10)

 

# Maximum modularity in 10 random starts: 0.8746
# Number of communities: 11
New_All_Data_2<- FindClusters(object = New_All_Data_2, resolution = 0.4)  

# Group cells, optimize standard modularity
New_All_Data_2<- RunTSNE(object = New_All_Data_2, dims = 1:pcSelect)                      
# TSNE clustering
tSNE_image <- TSNEPlot(object = New_All_Data_2, pt.size = 2, label = TRUE)    
# TSNE visualization
ggsave("10_tSNE_Image.png", plot = tSNE_image, width = 10, height = 10)
 

table(New_All_Data_2$RNA_snn_res.0.4)
 

# 8. Differential genes
# Find differential genes between certain clusters
#cluster5.markers = FindMarkers(New_All_Data_3, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
# Find all differential genes (takes very long)
New_All_Data_2.markers <- FindAllMarkers(object = New_All_Data_2, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(New_All_Data_2.markers,file="New_All_Data_2_All.markers.csv")
 
library(dplyr) 
top20 <- New_All_Data_2.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
 


# This plot takes too long to draw, so only select 20
New_All_Data_2_top20 <- DoHeatmap(New_All_Data_2,top20$gene[c(1,21,41,61,81,101,121,141,161,181,201,221)],size=5) +scale_fill_gradientn(colors = c("#A6CEE3",'#FF7F00'))
#c(1,21,41,61,81,101,121,141,161,181,201,221)
ggsave(filename="111_New_All_Data_2_top20.markers_heatmap.png", plot = New_All_Data_2_top20, width = 10, height = 15)

write.csv(top20,file="Top20G_Each_cluster.csv")
  
New_All_Data_2.markers <- FindAllMarkers(object = New_All_Data_2, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
top20 <- New_All_Data_2.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
# Select some genes
select_genes <- top20$gene[c(1,21,41,61,81,101,121,141,161)]
# vlnplot display
p1 <- VlnPlot(New_All_Data_2, features = select_genes, pt.size=0, ncol=3)
ggsave("132_selectgenes_6_VlnPlot.png", p1, width=15 ,height=10)
 
# featureplot display
p2 <- FeaturePlot(New_All_Data_2, features = select_genes, reduction = "tsne", label=T, ncol=3)
ggsave("133_selectgenes_6_FeaturePlot.png", p2, width=12 ,height=8)
 




#Step3: SingleR cell type identification
save(New_All_Data_2,file='New_All_Data_2.Rdata')
genes_to_check = c("CD3E","CD8A","CD4","FOXP3","CD19","MZB1","CD138","FCGR3A","NCR1","NCAM1","MARCO","MERTK","CD14","FCER1A","MHCII","IL3RA","CLE4C","NRP1","LILRA4","FUT4","CEACAM1","CEACAM6","CEACAM3","ENPP3")

#CD8+ T cells (CD3E+CD8A+ CD4−); 
#CD4+ T cells (CD3E+CD4+);
#Regulatory T cells (Tregs) (CD3E+FOXP3+); 
#B cells (CD19+); 
#Plasma cells (MZB1+CD138+); 
#NK cells (FCGR3A+NCR1+NCAM1+ CD3E-); 
#NKT cells (FCGR3A+NCR1+NCAM1+CD3E+ FOXP3-);
#Macrophages (MARCO+ MERTK+);
#Monocyte (CD14+,FCER1A+);
#cDCs dendritic cells(MHCII+CD4+ CD3E-);
#pDCs(IL3RA+CLE4C+NRP1+LILRA4+MHCII+ CD3E-);
#Neutrophils(FUT4+ CD3E-)
#Granulocytes(CEACAM1+CEACAM6+CEACAM3+ENPP3+ CD3E-)

p = DotPlot(New_All_Data_2, features = genes_to_check) + coord_flip()
ggsave("15_Tcell_distribution.png", p, width=20 ,height=15)
 

library(GSEABase) # BiocManager:: install( 'GSEABase ')
library(clusterProfiler)
geneset <- read.gmt('Hallmark.gmt')
length( unique (geneset$term) )   #9
library(fgsea)
library(dplyr)
library(tibble)
library(Seurat)
GSEA_ana <- function(geneList,geneSet,groupnum,clusternum){
	cluster.genes<- geneList %>% arrange(desc(avg_log2FC)) %>% dplyr::select(gene,avg_log2FC) 
cluster.genes <- cluster.genes[,c("gene","avg_log2FC")]
# Genes sorted by logFC
	ranks<- deframe(cluster.genes)
	egmt <- GSEA( ranks, TERM2GENE=geneSet,minGSSize = 1,pvalueCutoff = 0.99 ,verbose=FALSE)
	head(egmt)
	gsea_results_df <- egmt@result
	name=paste('group',groupnum,'cluster',clusternum,'gsea_results_df.csv',sep='_')
	write.csv(gsea_results_df,file = name)
}

a <- as.data.frame(table(New_All_Data_2.markers $cluster))
startRow=1
for (i in 0:(length(a$Var1)-1)){
	nRow <- a[(i+1),]$Freq
	endRow <- startRow+nRow-1
	if (length(intersect(New_All_Data_2.markers [startRow:endRow,]$gene,geneset$gene))!=0){
		GSEA_ana(New_All_Data_2.markers [startRow:endRow,],geneset,0,i)
	}
	startRow=endRow+1
}


# Rename clusters
new.cluster.ids<- c("CD8 T cells","CD4 T cells","Regulatory T cells","CD8 T cells","CD8 T cells","Monocyte","NKT cells","B cells","T cells","pDCs","Plasma cells","B cells")

names(new.cluster.ids) <- levels(New_All_Data_2)
New_All_Data_2<- RenameIdents(New_All_Data_2, new.cluster.ids)
p <- DimPlot(New_All_Data_2, reduction = "tsne", label = TRUE, pt.size = 1)
ggsave("16_HER2_subcluster_naming.png", p, width=10 ,height=8)
 
cellTypes <- c()
attach(New_All_Data_2@meta.data)
for ( i in 1:11653){
	ifelse(seurat_clusters[i]==0, cellTypes[i]<- "CD8 T cells", ifelse(seurat_clusters[i]==1, cellTypes[i]<- "CD4 T cells", ifelse(seurat_clusters[i]==2, cellTypes[i]<- "Regulatory T cells", ifelse(seurat_clusters[i]==3, cellTypes[i]<- "CD8 T cells", ifelse(seurat_clusters[i]==4, cellTypes[i]<- "CD8 T cells", ifelse(seurat_clusters[i]==5, cellTypes[i]<- "Monocyte", ifelse(seurat_clusters[i]==6, cellTypes[i]<- "NKT cells", ifelse(seurat_clusters[i]==7, cellTypes[i]<- "B cells", ifelse(seurat_clusters[i]==8, cellTypes[i]<- "T cells", ifelse(seurat_clusters[i]==9, cellTypes[i]<-"pDCs",ifelse(seurat_clusters[i]==10, cellTypes[i]<- "Plasma cells",cellTypes[i]<- "B cells")))))))))))
}
New_All_Data_2<- AddMetaData(New_All_Data_2, cellTypes,col.name = "cellTypes")

#Step4: Compare distribution differences of responder samples in each component
Extract_metadata<- New_All_Data_2@meta.data
attach(Extract_metadata)
Extract_metadata[1:6,]

 
R_group <- Extract_metadata[V6=='Responder',]
NR_group <- Extract_metadata[V6!='Responder',]
dim(R_group)
dim(NR_group)
   
R_Sample_Cell_TotalNum <- as.data.frame(table(R_group$V5))
NR_Sample_Cell_TotalNum <- as.data.frame(table(NR_group$V5))
   
PerSample_Cell_Proportion <- function(group,Sample_Cell_TotalNum){
	# group refers to analyzed sample group (responder or non-responder...)
	# Sample_Cell_TotalNum refers to dataframe of cell counts per sample in this group
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

R_Cell_Proportion<- PerSample_Cell_Proportion(R_group,R_Sample_Cell_TotalNum)
NR_Cell_Proportion<- PerSample_Cell_Proportion(NR_group,NR_Sample_Cell_TotalNum)
 
 
table(Extract_metadata$cellTypes)  
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
R_cell_Prop_list<-cell_Pro_dataframe(R_Cell_Proportion)
NR_cell_Prop_list<-cell_Pro_dataframe(NR_Cell_Proportion)
 

 
library(ggpubr)
cell_proporation_box <- function(iNum,celltype, R_cell_Prop, NR_cell_Prop){
	# R_cell_Prop refers to proportion vector of this cellType in Re, NR_cell_Prop refers to NoR
	a_num <- length(R_cell_Prop)
	b_num <- length(NR_cell_Prop)
	data <- as.numeric(c(R_cell_Prop,NR_cell_Prop))*100
	tag <- rep(c('R','NR'),c(a_num,b_num))
	all_data <-data.frame(as.numeric(data),tag)
	colnames(all_data)[1]<- 'ProP'
    attach(all_data)
	p <- ggboxplot(all_data, x = 'tag', y = 'ProP',
          color ='tag',
		  ylab='% of CD45+',
		  xlab='',
          palette = Mycol_2,
          add = "jitter",
		  shape='tag'
		  )
	detach(all_data)
	p = p + stat_compare_means(aes(group = tag),label = "p.format")
	name <- paste(iNum,celltype,"R_NR_Box.pdf",sep = "_")
	ggexport(p,filename =name,width = 2, height =3.5)
}
cells <- names(R_cell_Prop_list)
for (i in 1:9){
	cell_proporation_box(i,cells[i],R_cell_Prop=R_cell_Prop_list[[i]],NR_cell_Prop=NR_cell_Prop_list[[i]])
}
                
# Order is B cells, CD4 T cells, CD8 T cells, Monocyte, NKT cells, pDCs, Plasma cells, Regulatory T cells, T cells
# B cells, CD4 T cells tend to be enriched in R group
# pDCs, Regulatory T cells, T cells tend to be enriched in NR group
# Step4: Compare distribution differences of responder samples in each cluster
PerSample_Cell_Proportion <- function(group,Sample_Cell_TotalNum){
	# group refers to analyzed sample group (responder or non-responder...)
	# Sample_Cell_TotalNum refers to dataframe of cell counts per sample in this group
	Cell_Proportion <- list()
	for(i in 1:dim(Sample_Cell_TotalNum)[1]){
		one_sample <- group[group$V5==Sample_Cell_TotalNum$Var1[i],]
		all_clusters <- as.data.frame(table(one_sample$seurat_clusters))
		rownames(all_clusters) <- all_clusters$Var1
		Cell_Proportion[[i]]<- (all_clusters)/Sample_Cell_TotalNum$Freq[i]
		Cell_Proportion[[i]]$Var1<- all_clusters$Var1
	}
	Cell_Proportion
}

R_Cell_Proportion<- PerSample_Cell_Proportion(R_group,R_Sample_Cell_TotalNum)
NR_Cell_Proportion<- PerSample_Cell_Proportion(NR_group,NR_Sample_Cell_TotalNum)
   
table(Extract_metadata$seurat_clusters)
 
cell_Pro_dataframe <- function(Cell_Proportion){
	c0_Pro <- c()
	c1_Pro <- c()
	c2_Pro <- c()
	c3_Pro <- c()
	c4_Pro <- c()
	c5_Pro <- c()
	c6_Pro <- c()
	c7_Pro <- c()
	c8_Pro <- c()
	c9_Pro <- c()
	c10_Pro <- c()
	c11_Pro <- c()

	for(i in 1:length(Cell_Proportion)){
		if(sum(Cell_Proportion[[i]]$Var1==0)!=0){
			c0_Pro[i] <- Cell_Proportion[[i]][Cell_Proportion[[i]]$Var1==0,]$Freq
		}else c0_Pro[i]=0
		if(sum(Cell_Proportion[[i]]$Var1==1)!=0){
			c1_Pro[i] <- Cell_Proportion[[i]][Cell_Proportion[[i]]$Var1==1,]$Freq
		}else c1_Pro[i]=0
		if(sum(Cell_Proportion[[i]]$Var1==2)!=0){
			c2_Pro[i] <- Cell_Proportion[[i]][Cell_Proportion[[i]]$Var1==2,]$Freq
		}else c2_Pro[i]=0
		if(sum(Cell_Proportion[[i]]$Var1==3)!=0){
			c3_Pro[i] <- Cell_Proportion[[i]][Cell_Proportion[[i]]$Var1==3,]$Freq
		}else c3_Pro[i]=0
		if(sum(Cell_Proportion[[i]]$Var1==4)!=0){
			c4_Pro[i] <- Cell_Proportion[[i]][Cell_Proportion[[i]]$Var1==4,]$Freq
		}else c4_Pro[i]=0
		if(sum(Cell_Proportion[[i]]$Var1==5)!=0){
			c5_Pro[i] <- Cell_Proportion[[i]][Cell_Proportion[[i]]$Var1==5,]$Freq
		}else c5_Pro[i]=0
		if(sum(Cell_Proportion[[i]]$Var1==6)!=0){
			c6_Pro[i] <- Cell_Proportion[[i]][Cell_Proportion[[i]]$Var1==6,]$Freq
		}else c6_Pro[i]=0
		if(sum(Cell_Proportion[[i]]$Var1==7)!=0){
			c7_Pro[i] <- Cell_Proportion[[i]][Cell_Proportion[[i]]$Var1==7,]$Freq
		}else c7_Pro[i]=0
		if(sum(Cell_Proportion[[i]]$Var1==8)!=0){
			c8_Pro[i] <- Cell_Proportion[[i]][Cell_Proportion[[i]]$Var1==8,]$Freq
		}else c8_Pro[i]=0
		if(sum(Cell_Proportion[[i]]$Var1==9)!=0){
			c9_Pro[i] <- Cell_Proportion[[i]][Cell_Proportion[[i]]$Var1==9,]$Freq
		}else c9_Pro[i]=0
		if(sum(Cell_Proportion[[i]]$Var1==10)!=0){
			c10_Pro[i] <- Cell_Proportion[[i]][Cell_Proportion[[i]]$Var1==10,]$Freq
		}else c10_Pro[i]=0
		if(sum(Cell_Proportion[[i]]$Var1==11)!=0){
			c11_Pro[i] <- Cell_Proportion[[i]][Cell_Proportion[[i]]$Var1==11,]$Freq
		}else c11_Pro[i]=0
	}
	list(c0_Pro,c1_Pro,c2_Pro,c3_Pro,c4_Pro,c5_Pro,c6_Pro,c7_Pro,c8_Pro,c9_Pro,c10_Pro,c11_Pro)
}
R_cell_Prop_list<-cell_Pro_dataframe(R_Cell_Proportion)
NR_cell_Prop_list<-cell_Pro_dataframe(NR_Cell_Proportion)
 

 
library(ggpubr)
cell_proporation_box <- function(celltype, R_cell_Prop, NR_cell_Prop){
	# R_cell_Prop refers to proportion vector of this cellType in Re, NR_cell_Prop refers to NoR
	a_num <- length(R_cell_Prop)
	b_num <- length(NR_cell_Prop)
	data <- as.numeric(c(R_cell_Prop,NR_cell_Prop))*100
	tag <- rep(c('R','NR'),c(a_num,b_num))
	all_data <-data.frame(as.numeric(data),tag)
	colnames(all_data)[1]<- 'ProP'
	attach(all_data)
	p <- ggboxplot(all_data, x = 'tag', y = 'ProP',
          color ='tag',
		  ylab='% of CD45+',
		  xlab='',
          palette = Mycol_2,
          add = "jitter",
		  shape='tag'
		  )
	detach(all_data)
	p = p + stat_compare_means(aes(group = tag),label = "p.format")
	name <- paste(celltype,"R_NR_Box.pdf",sep = "_")
	ggexport(p,filename =name,width = 2, height =3.5)
}
for (i in 1:12){
	cell_proporation_box(i,R_cell_Prop=R_cell_Prop_list[[i]],NR_cell_Prop=NR_cell_Prop_list[[i]])
}
       
       
       
# CD4 T cells(c1) tend to be in R
# B cells(c7), CD8 T cells(c0,c4), pDCs(c9), Regulatory T cells(c2), T cells(c8) tend to be in NR

#Step4.2: New differential genes
save(New_All_Data_2,file='New_All_Data_2.Rdata')
# Find differential genes for self-defined cell groups
Idents(New_All_Data_2)="cellTypes"
Nine_celltype_markers <- FindAllMarkers(object = New_All_Data_2, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
# Check number of marker genes per cluster
table(Nine_celltype_markers $cluster)
 
write.csv(Nine_celltype_markers,file="Nine_celltype_markers.csv")
library(dplyr) 
top20 <- Nine_celltype_markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
write.csv(top20 ,file="Top20G_Nine_cluster.csv")
 
# This plot takes too long to draw, so only select 20
New_All_Data_2@assays$RNA@scale.data <- scale(New_All_Data_2@assays$RNA@data, scale = TRUE)
Nine_celltype_markers_top20 <- DoHeatmap(New_All_Data_2,top20$gene[c(1,21,41,61,81,101,121,141,161)]) +scale_fill_gradientn(colors = c("#A6CEE3",'#FF7F00'))
ggsave(filename="17_Nine_celltype_markers_top10_heatmap.png", plot =Nine_celltype_markers_top20, width = 12, height = 15)
  
0

# featureplot display
p2 <- FeaturePlot(New_All_Data_2, features = select_genes, reduction = "tsne", label=T, ncol=3)
ggsave("13_selectgenes_9_FeaturePlot.png", p2, width=12 ,height=8)
 
#Step4.3: Differential gene enrichment analysis (GO & KEGG)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
difgene <- top20$gene 
# What is the gene order (corresponding to cell types)
GO_Enrichment_Analy <- function(Gene_List,cellName){
	go_id_trance <- bitr(Gene_List,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db",drop = T)
	Gene_List_entrezid <- go_id_trance$ ENTREZID
	#1. GO enrichment
	##CC represents cellular component, MF represents molecular function, BP represents biological process, ALL represents enriching all three processes simultaneously, choose what you need, I usually do BP, MF, CC these 3 groups then merge into one dataframe for easy extraction of partial pathways for plotting later.
	ego_ALL <- enrichGO(gene = Gene_List_entrezid,# We defined above
					   OrgDb=org.Hs.eg.db,
					   keyType = "ENTREZID",
					   ont = "ALL",# Type of GO to enrich
					   pAdjustMethod = "BH",# Don't worry about this, usually use BH
					   minGSSize = 1,
					   pvalueCutoff = 0.01,# P value can take 0.05
					   qvalueCutoff = 0.05,
					   readable = TRUE)
					   
	ego_CC <- enrichGO(gene = Gene_List_entrezid,
					   OrgDb=org.Hs.eg.db,
					   keyType = "ENTREZID",
					   ont = "CC",
					   pAdjustMethod = "BH",
					   minGSSize = 1,
					   pvalueCutoff = 0.01,
					   qvalueCutoff = 0.05,
					   readable = TRUE)

	ego_BP <- enrichGO(gene = Gene_List_entrezid,
					   OrgDb=org.Hs.eg.db,
					   keyType = "ENTREZID",
					   ont = "BP",
					   pAdjustMethod = "BH",
					   minGSSize = 1,
					   pvalueCutoff = 0.01,
					   qvalueCutoff = 0.05,
					   readable = TRUE)

	ego_MF <- enrichGO(gene = Gene_List_entrezid,
					   OrgDb=org.Hs.eg.db,
					   keyType = "ENTREZID",
					   ont = "MF",
					   pAdjustMethod = "BH",
					   minGSSize = 1,
					   pvalueCutoff = 0.01,
					   qvalueCutoff = 0.05,
					   readable = TRUE)
	# These three numbers represent number of BP, CC, MF pathways selected, set this yourself (before controlling number of GO terms, save)
	ego_result_ALL <- as.data.frame(ego_ALL)
	ego_result_BP <- as.data.frame(ego_BP)
	ego_result_CC <- as.data.frame(ego_CC)
	ego_result_MF <- as.data.frame(ego_MF)
	filename <- paste('cell', cellName,"GO_Analyse.csv",sep='_')
	write.csv(ego_result_ALL,file=filename)
####### Draw chord diagram ######################
	if(dim(ego_ALL)[1]!=0){
	go_enrich_df_3 <- as.data.frame(ego_ALL)
	if(dim(ego_ALL)[1]>20){
		go_enrich_df_3 <- as.data.frame(ego_ALL)[1:20,]
	}
		gene <- str_replace_all (go_enrich_df_3$geneID,"/" , "," )
		GO <- data.frame(ID=go_enrich_df_3$ID,Term=go_enrich_df_3$Description,Genes=gene,adj_pval=go_enrich_df_3$p.adjust,Category=go_enrich_df_3$ONTOLOGY)

		## Build gene matrix and randomly generate logFC
		genedata=data.frame(ID=go_id_trance$SYMBOL,logFC=rnorm(length(go_id_trance$SYMBOL),mean=0,sd=2))

		circ <- circle_dat(GO, genedata)
		###### Plotting #################
		chord <- chord_dat(data=circ, genes = genedata)
		# Generate matrix with selected gene list
		chord <- chord_dat(data=circ, process = GO$Term) # Generate matrix with selected GO term list
		chord <- chord_dat (data=circ, genes=genedata, process = GO$Term) # Build data
		p=GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5,ribbon.col=Mycol_4[1:length(GO$Term)])
		i0=17
		imageName=paste(i0+j,cellName,'Chord_image.png',sep = "_")
		ggsave(file=imageName, plot = p, width = 6, height = 8)
	}
	go_enrich_df <- data.frame(
	ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
	Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
	GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
	type=factor(c(rep("biological process", dim(ego_result_BP)[1]), 
				  rep("cellular component", dim(ego_result_CC)[1]),
				  rep("molecular function", dim(ego_result_MF)[1])), 
				  levels=c("biological process", "cellular component","molecular function" )))
	# After controlling number of displayed pathways
	if(dim(ego_result_BP)[1]>5){
	ego_result_BP =ego_result_BP [1:5,]}
	if(dim(ego_result_CC)[1]>5){
		ego_result_CC =ego_result_CC [1:5,]}
	if(dim(ego_result_MF)[1]>5){
		ego_result_MF =ego_result_MF [1:5,]}

		## Recombine above extracted partial pathways into dataframe
		go_enrich_df_2 <- data.frame(
		ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
		Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
		GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
		type=factor(c(rep("biological process", dim(ego_result_BP)[1]), 
					  rep("cellular component", dim(ego_result_CC)[1]),
					  rep("molecular function", dim(ego_result_MF)[1])), 
					  levels=c("biological process", "cellular component","molecular function" )))
		
		## Pathway names are too long, I selected first five words of pathway as pathway name
		if(nrow(go_enrich_df_2)!=0){
			for(i in 1:nrow(go_enrich_df_2)){
			  description_splite=strsplit(go_enrich_df_2$Description[i],split = " ")
			  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") # This 5 refers to 5 words, can change yourself
			  go_enrich_df_2$Description[i]=description_collapse
			 }
		}
		## Start drawing GO bar plot
		  ### Horizontal bar plot
		go_enrich_df_2$type_order=factor(rev(as.integer(rownames(go_enrich_df_2))),labels=rev(go_enrich_df_2$Description))# This step is necessary, to make bars display in order, not messy

		i0=25
		p <- ggplot(data= go_enrich_df_2, aes(x=type_order,y=GeneNumber, fill=type)) + # Horizontal and vertical axis values
		  geom_bar(stat="identity", width=0.8) + # Bar plot width, can set yourself
		  scale_fill_manual(values = Mycol_3) + ### Color
		  coord_flip() + ## This step makes bar plot horizontal, without this bar plot is vertical
		  xlab("GO term") + 
		  ylab("Gene_Number") + 
		  labs(title = "The Most Enriched GO Terms")+
		  theme_bw()
		imageName<-paste(i0+j,cellName,'GO_Enrichment.png',sep = "_") 
		ggsave(file=imageName, plot = p, width = 12, height = 10)

}

KEGG_Enrichment_Analy<- function(Gene_List,cellName){
	go_id_trance <- bitr(Gene_List,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db",drop = T)
	Gene_List_entrezid <- go_id_trance$ ENTREZID

	#1. KEGG enrichment
	kk <- enrichKEGG(gene =Gene_List_entrezid,keyType = "kegg",organism= "human", qvalueCutoff = 0.05, pvalueCutoff=0.05)
	
	#2. Visualization
	  ### Bar plot
	hh <- as.data.frame(kk)# Remember to save results yourself!
	if(dim(hh)[1]!=0){
		filename <- paste('cell',cellName,"KEGG_Analyse.csv",sep='_')
		write.csv(hh,file=filename)
		rownames(hh) <- 1:nrow(hh)
		hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
		p=ggplot(hh,aes(y=order,x=Count,fill=p.adjust))+
		  geom_bar(stat = "identity", width=0.5)+#### Bar width
		  #coord_flip()+## Swap horizontal and vertical axes
		  scale_fill_gradient(low = '#FFB90F',high ="#D95F02")+# Can change color yourself
		  labs(title = "KEGG Pathways Enrichment",
			   x = "Gene numbers", 
			   y = "Pathways")+
		  theme(axis.title.x = element_text(face = "bold",size = 16),
				axis.title.y = element_text(face = "bold",size = 16),
				legend.title = element_text(face = "bold",size = 16))+
		  theme_bw()
		i0=25
		imageName<-paste(i0,cellName,'cellType_KEGG_Enrichment.png',sep = "_")
		ggsave(file=imageName, plot = p, width = 8, height = 6)
		
		  ### Bubble plot
		rownames(hh) <- 1:nrow(hh)
		hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
		p=ggplot(hh,aes(y=order,x=Count))+
		geom_point(aes(size=Count,color=-1*p.adjust))+# Modify point size
		scale_color_gradient(low="green",high = "red")+
		labs(color=expression(p.adjust,size="Count"), 
			 x="Gene Number",y="Pathways",title="KEGG Pathway Enrichment")+
		theme_bw()
		i0=25
		imageName<-paste(i0,cellName,'KEGG_Enrichment.png',sep = "_")
		ggsave(file=imageName, plot = p, width = 8, height = 6)
	}
}

i=1
cells <- as.vector(unique(top20$cluster))
for (j in 1:9){
	one_gene_list <- difgene[i:(i+19)]
	GO_Enrichment_Analy(one_gene_list, cells[j])
	KEGG_Enrichment_Analy(one_gene_list,cells[j])
	i=i+20
}
#### 2. Compare differentially expressed genes between two groups
# Compare differentially expressed genes between cluster0 and cluster1
New_All_Data_3 <- New_All_Data_2
dge.cluster <- FindMarkers(New_All_Data_3,ident.1 = 0,ident.2 = 1)
sig_dge.cluster <- subset(dge.cluster, p_val_adj<0.01&abs(avg_logFC)>1)
# Compare differentially expressed genes between B_cell and T_cells
dge.celltype <- FindMarkers(scRNA, ident.1 = 'B_cell', ident.2 = 'T_cells', group.by = 'celltype')
sig_dge.celltype <- subset(dge.celltype, p_val_adj<0.01&abs(avg_logFC)>1)
# Compare differentially expressed genes between pseudotime State1 and State3
p_data <- subset(pData(mycds),select='State')
scRNAsub <- subset(scRNA, cells=row.names(p_data))
scRNAsub <- AddMetaData(scRNAsub,p_data,col.name = 'State')
dge.State <- FindMarkers(scRNAsub, ident.1 = 1, ident.2 = 3, group.by = 'State')
sig_dge.State <- subset(dge.State, p_val_adj<0.01&abs(avg_logFC)>1)



