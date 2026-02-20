save.image()

# Import packages
library(Seurat)
library(dplyr)
library(Matrix)
library(magrittr)
library(stringr)       #library(tidyverse)
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
#Mycol_nejm= pal_nejm("default", alpha =0.7)(8)##Extract color
Mycol_4 <- c("#A6CEE3" ,"#1F78B4", "#B2DF8A" ,"#33A02C" ,"#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00", "#CAB2D6" ,"#6A3D9A" ,"#FFFF99", "#B15928","#1B9E77", "#D95F02" ,"#7570B3", "#E7298A" ,"#66A61E" ,"#E6AB02", "#A6761D","#666666")

# Find all differentially expressed genes (takes a very long time)
#Re_All_Data.markers <- FindAllMarkers(object = Re_All_Data, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
#NoR_All_Data.markers <- FindAllMarkers(object = NoR_All_Data, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
#Re_Pre_Data.markers <- FindAllMarkers(object = Re_Pre_Data, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
#Re_Post_Data.markers <- FindAllMarkers(object = Re_Post_Data, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
#NoR_Pre_Data.markers <- FindAllMarkers(object = NoR_Pre_Data, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
#NoR_Post_Data.markers <- FindAllMarkers(object = NoR_Post_Data, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

library(GSEABase) # BiocManager:: install( 'GSEABase ')
library(clusterProfiler)
geneset_B <- read.gmt('Hallmark_Bcells.gmt')
geneset_pdcs <- read.gmt('Hallmark_pDCs.gmt')
B_cells<- c('CD19')	
pDCs<- c('IL3RA','CLE4C','NRP1','LILRA4','MHCII')	
# Extract significantly differentially expressed genes
Re_All_Data_newMarkers <- subset(Re_All_Data.markers , p_val_adj<0.05&abs(avg_log2FC)>1)
 
library(GSEABase) # BiocManager:: install( 'GSEABase ')
library(clusterProfiler)
library(fgsea)
library(dplyr)
library(tibble)
library(Seurat)
# Extract significantly differentially expressed genes
Re_All_Data_newMarkers <- subset(Re_All_Data.markers , p_val_adj<0.05&abs(avg_log2FC)>2)
Re_All_Data_Bcells_Markers <- subset(Re_All_Data_newMarkers,cluster=="B cells")$gene
#n=122
NoR_All_Data_newMarkers <- subset(NoR_All_Data.markers , p_val_adj<0.05&abs(avg_log2FC)>2)
NoR_All_Data_pDCs_Markers <- subset(NoR_All_Data_newMarkers,cluster=="pDCs")$gene
#n=350
############### Discussion of Re_All_Data B cells expression profile ##############
rowsChoose <- rownames(Re_All_Data@assays$RNA@data) %in% Re_All_Data_Bcells_Markers
colsChoose <- Re_All_Data@meta.data$cellTypes=="B cells"

# Re_All_Data B cells expression profile data
Re_All_Bcells_Data <- Re_All_Data@assays$RNA@data[rowsChoose,colsChoose]
#122 360
write.csv(Re_All_Bcells_Data, "Re_All_Bcells_Data.csv")

# Re_All_Data B cells selected gene names
rows_geneNames_Bcells <- rownames(Re_All_Bcells_Data)

G_List_Bcells<- data.frame(Re_All_Bcells_Data)



############### Discussion of NoR_All_Data pDCs cells expression profile ##############
rowsChoose <- rownames(NoR_All_Data@assays$RNA@data) %in% NoR_All_Data_pDCs_Markers
colsChoose <- NoR_All_Data@meta.data$cellTypes=="pDCs"

# NoR_All_Data pDCs cells expression profile data
NoR_All_pDCs_Data <- NoR_All_Data@assays$RNA@data[rowsChoose,colsChoose]
# 350 201
write.csv(NoR_All_pDCs_Data, "NoR_All_pDCs_Data.csv")

# NoR_All_Data pDCs cells selected gene names
rows_geneNames_pDCs <- rownames(NoR_All_pDCs_Data)

G_List_pDCs<- data.frame(Re_All_Bcells_Data)







CD8 T cell re-clustering
############### Discussion of All_Data CD8 T cells expression profile ##############
load('New_All_Data_2.Rdata')
cd8_sce1 = New_All_Data_2[,New_All_Data_2@meta.data$cellTypes=="CD8 T cells"]
# All_Data CD8 T cells expression profile data
CD8Tcells_Data <- cd8_sce1@assays$RNA@data
# 55727 genes   5240 samples
write.csv(CD8Tcells_Data, "CD8Tcells_Data.csv")

# 2. Gene filtering (by Median Absolute Deviation, MAD)
mads <- apply(CD8Tcells_Data,1,mad) # MAD measure
df <- CD8Tcells_Data[rev(order(mads))[1:5000],] # Extract top 5000 genes
# 3. Normalization
df <-  as.matrix(sweep(df,1, apply(df,1,median,na.rm=T))) # Subtract row median, default is subtraction

# 4. Run ConsensusClusterPlus
library(ConsensusClusterPlus)
maxK <-  6 # Try a K value
results <-  ConsensusClusterPlus(df,
                                 maxK = maxK,
                                 reps = 1000,              # Number of resampling iterations (usually 1000 or more)
                                 pItem = 0.8,              # Resampling proportion
                                 pFeature = 1,
                                 clusterAlg = "pam",       # Clustering algorithm
                                 distance="pearson",       # Distance calculation method
                                 title="CD8_cells_Concluster", # Result save path
                                 innerLinkage="complete",  # Default method "average" is not recommended here
                                 plot="png")               # Result save format
# 5. Determine optimal number of clusters using PAC method
#   Minimum area corresponds to optimal K
Kvec = 2:maxK
x1 = 0.1; x2 = 0.9        # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="")  # from 2 to maxK
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])          # M is the calculated consensus matrix
  PAC[i-1] = Fn(x2) - Fn(x1)
} 

optK = 5 # Optimal K value is 5
    
 
 

icl = calcICL(results,
              title="CD8_cells_Concluster",
              plot="pdf")   # icl contains two dataframes: "clusterConsensus" cluster consistency scores;
                            # "itemConsensus" cluster corresponding samples and sample consistency scores 
icl[["clusterConsensus"]][1:5,]
icl[["itemConsensus"]][1:5,]
 
itemCon_CD8 <- icl[["itemConsensus"]]
cluster_item <- subset(itemCon,k==5,select = c("cluster","item","itemConsensus"))
cluster_item$item[1:20]
 
cluster_item_sorted <- cluster_item[order(cluster_item$item, cluster_item$itemConsensus,decreasing=T),]

cluster_item_res=c()
for (i in 1:dim(cluster_item_sorted)[1]){
  if (i%%optK == 1){
  cluster_item_res=rbind(cluster_item_res,cluster_item_sorted[i,])
  }
}
head(cluster_item_res)
 
# Object is cd8_sce1
# CD8Tcells_Data <- cd8_sce1@assays$RNA@data
attach(cluster_item_res)
CD8_Cluster <- c()
for (i in 1:dim(cluster_item_res)[1]){
	CD8_Cluster[rownames(cd8_sce1@meta.data)==item[i]] <- cluster[i]
}
CD8_Cluster   # Cluster for cd8_sce1 sample order
detach(cluster_item_res)
 
cd8_sce1<- AddMetaData(cd8_sce1, CD8_Cluster,col.name = "CD8_Cluster")
save(cd8_sce1,file='cd8_sce1.Rdata')

cd8_sce1<- AddMetaData(cd8_sce1, CD8_Cluster,col.name = "CD8_Cluster")
save(cd8_sce1,file='cd8_sce1.Rdata')


####### Data construction ################
# For 5 clusters based on mean expression profile data
FiveCluster_meanEXPdata <- c()
CD8Meta <- cd8_sce1@meta.data
for (i in 1:optK){
	choose_i_sample <- rownames(CD8Meta[CD8Meta$CD8_Cluster==i,])
	chooseData <- CD8Tcells_Data[,colnames(CD8Tcells_Data)%in%choose_i_sample]
	add_row_mean = apply(chooseData,1,mean)
	FiveCluster_meanEXPdata=cbind(FiveCluster_meanEXPdata,add_row_mean)
}
colnames(FiveCluster_meanEXPdata) <- c("C1","C2","C3","C4","C5")
head(FiveCluster_meanEXPdata)

 
library(ggdendro)
library(ggplot2)

########## First perform hierarchical clustering, obtain clustering tree representing similarity between samples #######
# Hierarchical clustering
# Data is
# head(FiveCluster_meanEXPdata)
# write.csv(FiveCluster_meanEXPdata, "FiveCluster_meanEXPdata.csv")
# FiveCluster_meanEXPdata <- read.csv("FiveCluster_meanEXPdata.csv",row.names=1)
# install.packages ("vegan")
# library(vegan)

# Calculate distance between samples, using commonly used Bray-curtis distance in community analysis
dis_bray <- vegan::vegdist(t(FiveCluster_meanEXPdata), method = 'bray')
 
# Hierarchical clustering
tree <- hclust(dis_bray, method = 'average')
plot(tree)
############### Adjust the clustering tree ###############
dend1 <- as.dendrogram(tree)
p1 <- plot(dend1, 
     nodePar = list(pch = 17:16, cex = 1.2:0.8, col = 2:3),
     horiz = TRUE)# Place clustering tree horizontally
# pch: Specify the symbol used to draw points, value range [0,24], where 4 is "minus sign", 20 is "dot"
# cex: Specify the size of the symbol. cex is a numerical value, representing the multiple of pch, default is 1.5 times
# Clustering tree drawing, color branches by group

###### Bar chart part ##################

#### Data preparation ########
# Cluster and R_NR sample count matrix
CD8Meta <- cd8_sce1@meta.data[c("Tag","CD8_Cluster")]
CD8Meta_1 <- as.matrix(table(CD8Meta$Tag,CD8Meta$CD8_Cluster))
rownames(CD8Meta_1) <- c("RB","RP","NRB","NRP")
colnames(CD8Meta_1) <- c("C1","C2","C3","C4","C5")
 
  
# Here we want to evaluate the response tendency of each cell, and simply counting the number of cells of each CD8 T cell subtype in each component is inaccurate (since from the collected samples, the number of samples in NR group is significantly more, so its cell number is also significantly more).
# We calculated the proportion of each subtype cell in each component.
res <- apply(CD8Meta_1,1,function(x){
	x/sum(x)*100
})
 
CD8Meta_2<- as.data.frame(as.table(t(res)))
 
CD8Meta_2 <- data.frame(CD8Meta_2$Var2,CD8Meta_2$Var1,CD8Meta_2$Freq)
colnames(CD8Meta_2) <- c("cluster","Tag","SampleNum")
 
# write.csv(CD8Meta_2, "CD8Meta_2.csv")
# CD8Meta_2 <- read.csv("CD8Meta_2.csv",row.names=1)
 
CD8Meta_2$cluster = factor(CD8Meta_2$cluster, levels=c("C3","C1","C4","C2","C5")) ## Set bar order

aa1<-seq(0,1,by=.2)
p2<- ggplot(data=CD8Meta_2, aes(x=cluster, y=SampleNum, fill=Tag)) +
      geom_bar(stat="identity",width=1,color='black',
               position =position_fill(reverse=TRUE))+
      xlab('')+ylab('')+coord_flip()+
      scale_y_continuous(expand = c(0, 0),breaks=aa1)+ 
      theme(#axis.text.y = element_blank(),
            axis.text.x = element_text(size=10),
            axis.ticks.y=element_blank(),
            axis.line.x=element_line(colour="black"),
            panel.background=element_rect(fill="white")
            )
 
# For the various CD8 cell subtypes defined above, we need to screen their markers
# Find differentially expressed genes for self-defined clusters
Idents(cd8_sce1)="CD8_Cluster"
CD8_celltype_markers <- FindAllMarkers(object = cd8_sce1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
# Check number of marker genes per cluster
table(CD8_celltype_markers $cluster)
 

write.csv(CD8_celltype_markers,file="CD8_celltype_markers.csv")
library(dplyr) 
CD8_clusters_top60<- CD8_celltype_markers %>% group_by(cluster) %>% top_n(n=60, wt=avg_log2FC)
write.csv(CD8_clusters_top60 ,file="CD8_clusters_top60.csv")
# Select some genes
startNum <- 1
for (i in 1:5){
	
	endNum <- 60*i
	select_genes <- CD8_clusters_top60$gene[startNum:endNum]
	startNum=endNum+1
	p1 <- VlnPlot(cd8_sce1, features = select_genes, pt.size=0, ncol=4)
	Fname <- paste(380+i,i,"CD8subType_selectgenes_VlnPlot.png",sep="_")
	ggsave(Fname, p1, width=30,height=20)
}
###### Marker screening process ↑↑↑        ###############################

select_genes <- CD8_clusters_top60$gene[c(21,35,45,81,121,122,186,192,261)]
# vlnplot display
p1 <- VlnPlot(cd8_sce1, features = select_genes, pt.size=0, ncol=3)
ggsave("386_CD8subType_selectgenes_VlnPlot.png", p1, width=15 ,height=10)
 


# Total markers
CD8_C1_Marker <- CD8_clusters_top60$gene[c(21,35,45,49)]
# CD8_c1 markers: "ZBP1","LINC00861","FBXO44","RP11-761N21.2"
CD8_C2_Marker <- CD8_clusters_top60$gene[c(81,83,84,107,115)]
# CD8_c2 markers: "LMNA","KIAA1683","RP11-442H21.2","RALGAPA1","PPP1R10"
CD8_C3_Marker <- CD8_clusters_top60$gene[c(121:130)]
# CD8_c3 markers: "TYMS","ZWINT","TK1","CCNB2","KIAA0101","BIRC5","CENPW","ASF1B","SPC24","CDCA5"
CD8_C4_Marker <- CD8_clusters_top60$gene[c(186,192,198,199,211,212,228,231,232,233)]
# CD8_c4 markers: "IL7R","TCF7","FTH1P2","FTH1P23","BEST1","FTH1P20","SATB1","RPS4XP13","EEF1A1P16","MGAT4A"
CD8_C5_Marker <- CD8_clusters_top60$gene[c(261,276,284,288,289,295,297,298,299)]
# CD8_c5 markers: "MYO7A","PTMS","DUSP5","UBE2F","SEC14L1","KIAA0319L","ITPR1","SPEN","ITGA1"


# Find differentially expressed genes for self-defined clusters
cd8ClusterData <- cd8_sce1@meta.data$CD8_Cluster
 
##### Find differentially expressed genes for self-defined clusters  #########

####(1,2,5)vs(4,3)#######################
RvsNR_CD8Marker <- c()
for (i in 1:length(cd8ClusterData)){
	if(cd8ClusterData[i] %in% c(4,3)){
		RvsNR_CD8Marker[i]="NR"
	}else RvsNR_CD8Marker[i]="R"
}
cd8_sce1<- AddMetaData(cd8_sce1, RvsNR_CD8Marker,col.name = "RvsNR_CD8Marker")

Idents(cd8_sce1)="RvsNR_CD8Marker"
CD8_RvsNR_markers <- FindAllMarkers(object = cd8_sce1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
# Check number of marker genes per cluster
table(CD8_RvsNR_markers $cluster)
 
top100 <- CD8_RvsNR_markers %>% group_by(cluster) %>% top_n(n=100, wt=avg_log2FC)
startNum <- 1
for (i in 1:4){
	
	endNum <- 50*i
	select_genes <- top100$gene[startNum:endNum]
	startNum=endNum+1
	p1 <- VlnPlot(cd8_sce1, features = select_genes, pt.size=0, ncol=5)
	Fname <- paste(389+i,i,"CD8_RvsNR_markers_VlnPlot.png",sep="_")
	ggsave(Fname, p1, width=40,height=30)
}

select_genes <- top100$gene[c(34,73,77,101,102,103)]
# vlnplot display
p1 <- VlnPlot(cd8_sce1, features = select_genes, pt.size=0, ncol=3)
ggsave("394_CD8_RvsNR_markers_VlnPlot.png", p1, width=10 ,height=5)
R_CD8_Marker <- top100$gene[c(34,73,77,85,91)]
NR_CD8_Marker <- top100$gene[c(101:110)]
 
# R group markers: "AC069363.1","ITGA1","SYNC","PTPN4","CRIM1"
# NR group markers: "TYMS","CDCA7","KIAA0101","NME1","STMN1","WDR76","MRPL37","TMEM106C","DUT","CKS1B" 



####(1)vs(2,5)vs(4)vs(3)#####################
RBvsRPvsNRBvsNRP_CD8Marker <- c()
for (i in 1:length(cd8ClusterData)){
	if(cd8ClusterData[i] %in% c(2,5)){
		RBvsRPvsNRBvsNRP_CD8Marker[i]="RP"
	}
	else if (cd8ClusterData[i] == 4){
		RBvsRPvsNRBvsNRP_CD8Marker[i]="NRB"
	}
	else if (cd8ClusterData[i] == 3){
		RBvsRPvsNRBvsNRP_CD8Marker[i]="NRP"
	}else RBvsRPvsNRBvsNRP_CD8Marker[i]="RB"
}

cd8_sce1<- AddMetaData(cd8_sce1, RBvsRPvsNRBvsNRP_CD8Marker,col.name = "RBvsRPvsNRBvsNRP_CD8Marker")

Idents(cd8_sce1)="RBvsRPvsNRBvsNRP_CD8Marker"
CD8_RBvsRPvsNRBvsNRP_markers <- FindAllMarkers(object = cd8_sce1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
# Check number of marker genes per cluster
table(CD8_RBvsRPvsNRBvsNRP_markers $cluster)
  
top50 <- CD8_RBvsRPvsNRBvsNRP_markers %>% group_by(cluster) %>% top_n(n=50, wt=avg_log2FC)

startNum <- 1
for (i in 1:4){
	
	endNum <- 50*i
	select_genes <- top50$gene[startNum:endNum]
	startNum=endNum+1
	p1 <- VlnPlot(cd8_sce1, features = select_genes, pt.size=0, ncol=5)
	Fname <- paste(394+i,i,"CD8_RvsNR_markers_VlnPlot.png",sep="_")
	ggsave(Fname, p1, width=40,height=30)
}

select_genes <- top50$gene[c(17,30,38,76,101,102,155,160,166)]

# vlnplot display
p1 <- VlnPlot(cd8_sce1, features = select_genes, pt.size=0, ncol=3)
ggsave("399_RBvsRPvsNRBvsNRP_CD8Marker_VlnPlot.png", p1, width=15 ,height=10)
 
 


RB_CD8_Marker <- top50$gene[c(17,30,38,42)]
RP_CD8_Marker <- top50$gene[c(73,74,75,76,91)]
NRP_CD8_Marker <- top50$gene[c(101:110)]
NRB_CD8_Marker <- top50$gene[c(155,160,166,167,178,179,185,192,195,196,197)]

# First is RB group: "ZBP1","LINC00861","FBXO44","RP11-761N21.2"
## Second is RP group: "PDE3B","TGIF1","VCAM1","DFNB31","RALGAPA1"
## Third is NRP group: "ZWINT","TK1","CCNB2","KIAA0101","BIRC5","CENPW","ASF1B","SPC24","CDCA5","DHFRP1"
# Fourth is NRB group: "IL7R","TCF7","FTH1P2","FTH1P23","BEST1","FTH1P20","LMNA","SATB1","RPS4XP13","EEF1A1P16","MGAT4A"
 
Note: Dendrogram only displays top 7 genes at most









GSEA_ana <- function(geneList,geneSet){
	cluster.genes<- geneList %>% arrange(desc(Data)) %>% dplyr::select(Gene,Data) # Genes sorted by geneEXPdata
	ranks<- deframe(cluster.genes)
	egmt <- GSEA( ranks, TERM2GENE=geneSet,minGSSize = 1,pvalueCutoff = 1 ,verbose=FALSE)
	head(egmt)
	gsea_results_df <- egmt@result
	rbind(enrich_re,gsea_results_df)
}

enrich_re <- data.frame()
for (i in 1:dim(G_List)[2]){
	GL <- data.frame(G_List[,i],rownames(G_List))
	names(GL)<- c("Data","Gene")
	if(sum(GL[GL$Gene %in% B_cells,]$Data)!=0){
		enrich_re =GSEA_ana(GL,geneset_B)
	}else enrich_re=rbind(enrich_re,0)
}
rownames(enrich_re)<- colnames(Re_All_Bcells_Data)
### Re group B cells each cell contribution score is enrich_re$enrichmentScore
