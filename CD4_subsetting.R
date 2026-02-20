# Extract CD4 cell expression profile and metadata information
CD4_sce1 = New_All_Data_2[,New_All_Data_2@meta.data$cellTypes=="CD4 T cells"]
# All_Data CD4 cells expression profile data
CD4_cells_Data <- CD4_sce1@assays$RNA@data
# 55727 genes, 1951 samples
write.csv(CD4_cells_Data, "CD4_cells_Data.csv")


# Extract CD4 cell marker genes related to R
CD4cells_Marker <- Re_All_Data_1.markers[Re_All_Data_1.markers$cluster == "CD4 T cells",]

######## Select top 20%  #69 genes  ############
CD4cells_Marker_top_genes <- CD4cells_Marker %>% top_frac(.2,wt=avg_log2FC)  
dim(CD4cells_Marker_top_genes)

select_CD4_top_genes <- CD4cells_Marker_top_genes$gene
p1 <- VlnPlot (Re_All_Data_1 , features =select_CD4_top_genes[1:9]  , pt.size=0, ncol=3)
ggsave("400_select_CD4_top_genes_VlnPlot.png", p1, width=10 ,height=8)
select_CD4_top_genes[1:9]
 


# Calculate mean expression of CD4 cell specific genes in Re_All_Data_1
Re_All_ExpData <- Re_All_Data_1@assays$RNA@data
select_CD4_top_genes_ExpData <- Re_All_ExpData[rownames(Re_All_ExpData) %in% select_CD4_top_genes,Re_All_Data_1@meta.data$cellTypes=="CD4 T cells"]
dim(select_CD4_top_genes_ExpData)
# 69 genes, 844 samples

select_CD4_top_genes_meanExp<- apply(select_CD4_top_genes_ExpData,1,mean)
select_CD4_top_genes_meanExp
write.table(select_CD4_top_genes_meanExp,"select_CD4_top_genes_meanExp.txt")
 





############### Analysis of All_Data CD4 cell expression profile ##############
# 2. Gene filtering (by Median Absolute Deviation, MAD)
mads <- apply(CD4_cells_Data,1,mad) # MAD measure
df <- CD4_cells_Data[rev(order(mads))[1:5000],] # Extract top 5000 genes
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
                                 title="CD4_cells_Concluster", # Result save path
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

optK = 5  # Optimal K value is 5

icl = calcICL(results,
              title="CD4_cells_Concluster",
              plot="pdf")   # icl contains two dataframes: "clusterConsensus" cluster consistency scores;
                            # "itemConsensus" cluster corresponding samples and sample consistency scores 
icl[["clusterConsensus"]][1:5,]
icl[["itemConsensus"]][1:5,]
 

itemCon_CD4_cells <- icl[["itemConsensus"]]
cluster_item <- subset(itemCon_CD4_cells,k==5,select = c("cluster","item","itemConsensus"))
cluster_item$item[1:20]
 

cluster_item_sorted <- cluster_item[order(cluster_item$item, cluster_item$itemConsensus,decreasing=T),]

cluster_item_res=c()
for (i in 1:dim(cluster_item_sorted)[1]){
  if (i%%optK == 1){
  cluster_item_res=rbind(cluster_item_res,cluster_item_sorted[i,])
  }
}
head(cluster_item_res)
 


# Object is CD4_sce1
# CD4_Data <- CD4_sce1@assays$RNA@data
attach(cluster_item_res)
CD4_Cluster <- c()
for (i in 1:dim(cluster_item_res)[1]){
	CD4_Cluster[rownames(CD4_sce1@meta.data)==item[i]] <- cluster[i]
}
CD4_Cluster   # Cluster for CD4_sce1 sample order
detach(cluster_item_res)
 



CD4_sce1<- AddMetaData(CD4_sce1, CD4_Cluster,col.name = "CD4_Cluster")
save(CD4_sce1,file='CD4_sce1.Rdata')

CD4_C1_Sample <- colnames(CD4_cells_Data[,CD4_sce1@meta.data$CD4_Cluster == 1])
CD4_C2_Sample <- colnames(CD4_cells_Data[,CD4_sce1@meta.data$CD4_Cluster == 2])
CD4_C3_Sample <- colnames(CD4_cells_Data[,CD4_sce1@meta.data$CD4_Cluster == 3])
CD4_C4_Sample <- colnames(CD4_cells_Data[,CD4_sce1@meta.data$CD4_Cluster == 4])
CD4_C5_Sample <- colnames(CD4_cells_Data[,CD4_sce1@meta.data$CD4_Cluster == 5])

length(CD4_C1_Sample)   #605
length(CD4_C2_Sample)   #450
length(CD4_C3_Sample)   #345
length(CD4_C4_Sample)   #313
length(CD4_C5_Sample)   #238


All_ExpData <- New_All_Data_2@assays$RNA@data
CD4_C1_select_genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% select_CD4_top_genes,colnames(All_ExpData) %in% CD4_C1_Sample]
dim(CD4_C1_select_genes_ExpData)
# 69 genes, 605 samples
CD4_C2_select_genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% select_CD4_top_genes,colnames(All_ExpData) %in% CD4_C2_Sample]
dim(CD4_C2_select_genes_ExpData)
# 69 genes, 450 samples
CD4_C3_select_genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% select_CD4_top_genes,colnames(All_ExpData) %in% CD4_C3_Sample]
dim(CD4_C3_select_genes_ExpData)
# 69 genes, 345 samples
CD4_C4_select_genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% select_CD4_top_genes,colnames(All_ExpData) %in% CD4_C4_Sample]
dim(CD4_C4_select_genes_ExpData)
# 69 genes, 313 samples
CD4_C5_select_genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% select_CD4_top_genes,colnames(All_ExpData) %in% CD4_C5_Sample]
dim(CD4_C5_select_genes_ExpData)
# 69 genes, 238 samples

CD4_C1_select_genes_meanExp<- apply(CD4_C1_select_genes_ExpData,1,mean)
CD4_C2_select_genes_meanExp<- apply(CD4_C2_select_genes_ExpData,1,mean)
CD4_C3_select_genes_meanExp<- apply(CD4_C3_select_genes_ExpData,1,mean)
CD4_C4_select_genes_meanExp<- apply(CD4_C4_select_genes_ExpData,1,mean)
CD4_C5_select_genes_meanExp<- apply(CD4_C5_select_genes_ExpData,1,mean)


######### Correlation matrix calculation ##################
CD4subType_cor <- t(rbind(select_CD4_top_genes_meanExp,CD4_C1_select_genes_meanExp,CD4_C2_select_genes_meanExp,CD4_C3_select_genes_meanExp,CD4_C4_select_genes_meanExp,CD4_C5_select_genes_meanExp))
colnames(CD4subType_cor)=c("Re_CD4Tcells","CD4_c1","CD4_c2","CD4_c3","CD4_c4","CD4_c5")

res <- cor(CD4subType_cor)
res2 <- signif(res, digits = 3)
res2
library(Hmisc)# Load package
 

Hmisc::rcorr(as.matrix(res2), type = "pearson") -> corrlist
# Note: Need to install Hmisc package first, otherwise it will error, same below, too lazy to write automated code, if many people ask me to write, I might勉强也行, put the matrix or dataframe that needs correlation calculation inside as.matrix()

library(tidyr)
# Correlation coefficient matrix
corrlist$r %>%  # Extract r matrix
  as_tibble() %>%  # Set to tibble format
  mutate(v = colnames(.)) %>%  # Take r matrix column name vector as a new vector in column
  select(v, everything()) %>%  # Select all tables
  pivot_longer(2:7) -> corrdf # Convert from wide to long data

# p value matrix
corrlist$P %>% 
  as_tibble() %>% 
  mutate(v = colnames(.)) %>% 
  select(v, everything()) %>% 
  pivot_longer(2:7) %>%  # 2:7 means I have 7 variables
  mutate(label = case_when(  # Set label, add judgment, when P value meets specific conditions display "\n" plus specific number of *
    is.na(value) ~ " ", # NA value assigned to space
    value <= 0.001 ~ "\n***", # P<0.001 display carriage return plus three asterisks
    between(value, 0.001, 0.01) ~ "\n**", # P 0.001-0.01 display carriage return plus two *
    between(value, 0.01, 0.05) ~ "\n*", # P 0.01-0.05 display carriage return plus one *
    T ~ ""
  )) -> pdf
corrdf %>% 
  left_join(pdf, by = c("v", "name")) %>% # Merge r value and p matrix by v and name fields respectively
  rename(corr = value.x, p = value.y) %>% # Rename value.x in merged dataframe to corr, value.y to p
  mutate(corr = signif(corr, digits = 3)) -> corrdf # Change decimal to two digits
 

library(grDevices)
# windowsFonts("Arial" = windowsFont("Arial")) # Set font to prevent error in code below
corrdf %>% 
  mutate(v = forcats::fct_reorder(v, corr), # Reorder v and name
         name = forcats::fct_reorder(name, corr)) ->corrdf2
corrdf2		 
		
 

write.csv(corrdf2,"corrdf2_CD4.csv")
corrdf2 <- read.csv("corrdf2_CD4.csv",row.names=1)
# library(dplyr)  pipe %>% package

# Red and blue ↓
corrdf2 %>%  
  ggplot(aes(x = v, y = name)) + 
  geom_tile(aes(fill = corr)) + # Header fill mapping based on corr
  geom_text(aes(label = paste(corr, label, sep = "")),size=6) + ggthemes::scale_fill_gradient2_tableau("Red-Blue Diverging")+ theme(axis.text=element_text(size=14),axis.title.x=element_blank(),axis.title.y=element_blank())
# Blue and orange ↓
corrdf2 %>%  
  ggplot(aes(x = v, y = name)) + 
  geom_tile(aes(fill = corr)) + # Header fill mapping based on corr
  geom_text(aes(label = paste(corr, label, sep = "")),size=6)+scale_fill_gradient(low='#1F77B4FF',high='#FF7F0EFF')+ theme(axis.text=element_text(size=14),axis.title.x=element_blank(),axis.title.y=element_blank())
 


# For the various CD4 cell subtypes defined above, we need to screen their markers
# Find differentially expressed genes for self-defined clusters
Idents(CD4_sce1)="CD4_Cluster"
CD4_celltype_markers <- FindAllMarkers(object = CD4_sce1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
# Check number of marker genes per cluster
table(CD4_celltype_markers $cluster)
# write.csv(CD4_celltype_markers,file="CD4_celltype_markers.csv")
CD4_celltype_markers <- read.csv("CD4_celltype_markers.csv")
library(dplyr)
 

### Select major cell type specific genes ##############
CD4cells_Marker <- Re_All_Data_1.markers[Re_All_Data_1.markers$cluster == "CD4 T cells",]
######## Select top 20%  #69 genes  ############
CD4cells_Marker_top_genes <- CD4cells_Marker %>% top_frac(.2,wt=avg_log2FC)  
dim(CD4cells_Marker_top_genes)
write.csv(CD4cells_Marker_top_genes, "CD4cells_Marker_top20%_genes.csv")

########### Gene extraction ########
select_CD4_top_genes <- CD4cells_Marker_top_genes $gene

# Select subtype specific genes ######
select_genes_CD4_c1 <- (CD4_celltype_markers[CD4_celltype_markers$cluster==1,] %>% top_n(n=30,wt=avg_log2FC))$gene     #top30
write.csv(CD4cells_Marker_top_genes, "CD4cells_Marker_top20%_genes.csv")

res <- c()
CD4_inter_c1_genes <- intersect(select_genes_CD4_c1,select_CD4_top_genes)
res[1]=length(CD4_inter_c1_genes)
res

# top20-c1:"LINC00861"
# top30-c1:"LINC00861" "TMEM63A" 
# Total 2 genes: "LINC00861" "TMEM63A"

CD4Tcells-top20%-DiffGenes=69, CD4-c1-DiffGenes=30, CD4Tcells-top20%-DiffGenes & CD4-c1-DiffGenes =2
# Visualization website  https://asntech.shinyapps.io/intervene/
# Parameters: color: #E09A4A  #317EAB
#         Plot width: 440, Plot height=300, ratio=0.5 , font all 1.5
# Connecting point size=4, Connecting point size=1


############### Analysis of All_Data CD8 T cell expression profile ##############
load('New_All_Data_2.Rdata')
cells <- c("B cells","CD4 T cells","pDCs","Regulatory T cells","T cells")

CD4_sce1 = New_All_Data_2[,New_All_Data_2@meta.data$cellTypes=="CD4 T cells"]
# All_Data CD4 T cell expression profile data
CD4cells_Data <- CD4_sce1@assays$RNA@data
# 55727 genes, 1951 samples
write.csv(CD4cells_Data, "CD4cells_Data.csv")

# 2. Gene filtering (by Median Absolute Deviation, MAD)
mads <- apply(CD4cells_Data,1,mad) # MAD measure
df <- CD4cells_Data[rev(order(mads))[1:5000],] # Extract top 5000 genes
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
                                 title="CD4cell_Concluster", # Result save path
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

optK = 5  # Optimal K value is 5

icl = calcICL(results,
              title="CD4cell_Concluster2",
              plot="pdf")   # icl contains two dataframes: "clusterConsensus" cluster consistency scores;
                            # "itemConsensus" cluster corresponding samples and sample consistency scores 
icl[["clusterConsensus"]][1:5,]
icl[["itemConsensus"]][1:5,]
   
   
 

itemCon_CD4cells <- icl[["itemConsensus"]]
cluster_item_CD4cells <- subset(itemCon_CD4cells,k==5,select = c("cluster","item","itemConsensus"))
cluster_item_CD4cells$item[1:20]
 

cluster_item_CD4cells_sorted <- cluster_item_CD4cells[order(cluster_item_CD4cells$item, cluster_item_CD4cells$itemConsensus,decreasing=T),]

cluster_item_CD4cells_res=c()
for (i in 1:dim(cluster_item_CD4cells_sorted)[1]){
  if (i%%optK == 1){
  cluster_item_CD4cells_res=rbind(cluster_item_CD4cells_res,cluster_item_CD4cells_sorted[i,])
  }
}
head(cluster_item_CD4cells_res)
 


# Object is CD4_sce1
# CD4cells_Data <- CD4_sce1@assays$RNA@data
attach(cluster_item_CD4cells_res)
CD4_Cluster <- c()
for (i in 1:dim(cluster_item_CD4cells_res)[1]){
	CD4_Cluster[rownames(CD4_sce1@meta.data)==item[i]] <- cluster[i]
}
CD4_Cluster   # Cluster for CD4_sce1 sample order
detach(cluster_item_CD4cells_res)
 

CD4_sce1<- AddMetaData(CD4_sce1, CD4_Cluster,col.name = "CD4_Cluster")
head(CD4_sce1@meta.data)
 

 
####### Data construction ################
# For 5 clusters based on mean expression profile data

Five_CD4_Clusters_meanEXPdata <- c()
CD4Meta <- CD4_sce1@meta.data
for (i in 1:optK){
	choose_i_sample <- rownames(CD4Meta[CD4Meta$CD4_Cluster==i,])
	chooseData <- CD4cells_Data[,colnames(CD4cells_Data)%in%choose_i_sample]
	add_row_mean = apply(chooseData,1,mean)
	Five_CD4_Clusters_meanEXPdata=cbind(Five_CD4_Clusters_meanEXPdata,add_row_mean)
}
colnames(Five_CD4_Clusters_meanEXPdata) <- c("C1","C2","C3","C4","C5")
head(Five_CD4_Clusters_meanEXPdata)

 

library(ggdendro)
library(ggplot2)

########## First perform hierarchical clustering, obtain clustering tree representing similarity between samples #######
# Hierarchical clustering
# Data is
# head(Five_CD4_Clusters_meanEXPdata)
# write.csv(Five_CD4_Clusters_meanEXPdata, "Five_CD4_Clusters_meanEXPdata.csv")
# Five_CD4_Clusters_meanEXPdata <- read.csv("Five_CD4_Clusters_meanEXPdata.csv",row.names=1)
# install.packages ("vegan")
library(vegan)

# Calculate distance between samples, using commonly used Bray-curtis distance in community analysis
dis_bray_CD4cells <- vegan::vegdist(t(Five_CD4_Clusters_meanEXPdata), method = 'bray')
 

# Hierarchical clustering
tree_Bcells <- hclust(dis_bray_CD4cells, method = 'average')
plot(tree_Bcells)


############### Adjust the clustering tree ###############
dend_B <- as.dendrogram(tree_Bcells)
p1 <- plot(dend_B, 
     nodePar = list(pch = 17:16, cex = 1.2:0.8, col = 2:3),
     horiz = TRUE)# Place clustering tree horizontally
# pch: Specify the symbol used to draw points, value range [0,24], where 4 is "minus sign", 20 is "dot"
# cex: Specify the size of the symbol. cex is a numerical value, representing the multiple of pch, default is 1.5 times
# Clustering tree drawing, color branches by group

###### Bar chart part ##################

#### Data preparation ########
# Cluster and R_NR sample count matrix
CD4Meta <- CD4_sce1@meta.data[c("Tag","CD4_Cluster")]
CD4Meta_1 <- as.matrix(table(CD4Meta$Tag,CD4Meta$CD4_Cluster))
rownames(CD4Meta_1) <- c("RB","RP","NRB","NRP")
colnames(CD4Meta_1) <- c("C1","C2","C3","C4","C5")
 

# Here we want to evaluate the response tendency of each cell, and simply counting the number of cells of each CD4 T cell subtype in each component is inaccurate (since from the collected samples, the number of samples in NR group is significantly more, so its cell number is also significantly more).
# We calculated the proportion of each subtype cell in each component.

res <- apply(CD4Meta_1,1,function(x){
	x/sum(x)*100
})
 

CD4Meta_2<- as.data.frame(as.table(t(res)))
 

CD4Meta_2 <- data.frame(CD4Meta_2$Var2,CD4Meta_2$Var1,CD4Meta_2$Freq)
colnames(CD4Meta_2) <- c("cluster","Tag","SampleNum")
 

# write.csv(CD4Meta_2, "CD4Meta_2.csv")
# CD4Meta_2 <- read.csv("CD4Meta_2.csv",row.names=1)
 

CD4Meta_2$cluster = factor(CD4Meta_2$cluster, levels=c("C5","C1","C2","C3","C4")) ## Set bar order


aa1<-seq(0,1,by=.2)
p2<- ggplot(data=CD4Meta_2, aes(x=cluster, y=SampleNum, fill=Tag)) +
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
Idents(CD4_sce1)="CD4_Cluster"
CD4_celltype_markers <- FindAllMarkers(object = CD4_sce1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
# Check number of marker genes per cluster
table(CD4_celltype_markers $cluster)
 


write.csv(CD4_celltype_markers,file="CD4_celltype_markers.csv")
library(dplyr) 
CD4_clusters_top50<- CD4_celltype_markers %>% group_by(cluster) %>% top_n(n=50, wt=avg_log2FC)
write.csv(CD4_clusters_top50 ,file="CD4_clusters_top50.csv")
# Select some genes
startNum <- 1
for (i in 1:optK){
	
	endNum <- 50*i
	select_genes <- CD4_clusters_top50$gene[startNum:endNum]
	startNum=endNum+1
	p1 <- VlnPlot(CD4_sce1, features = select_genes, pt.size=0, ncol=4)
	Fname <- paste(349+i,i,"CD4subType_selectgenes_VlnPlot.png",sep="_")
	ggsave(Fname, p1, width=30,height=20)
}
###### Marker screening process ↑↑↑        ###############################

select_genes <- CD4_clusters_top50$gene[c(5,51,117,172,205)]
# vlnplot display
p1 <- VlnPlot(CD4_sce1, features = select_genes, pt.size=0, ncol=3)
ggsave("355_CD8subType_selectgenes_VlnPlot.png", p1, width=10 ,height=4)

 

# Total markers
CD4_M1 <- CD4_clusters_top50$gene[c(5,8,14,19,25,33,34,39,46,47,48)]
# "CD8B"      "LINC00623" "INO80D"    "CYTH4"     "YIPF4"     "CTC1"      "A2M"       "DBNL"      "C19orf66"  "SPN"       "UBA7"

CD4_M2 <- CD4_clusters_top50$gene[c(51,54,55,56,57,59,63,66,77,80,82,86,88,90,91,93,96)]
# "CTLA4"       "TIGIT"       "HLA-DRB1"    "TOX"         "PDCD1"      "HLA-DRB5"    "MIR4435-1HG" "CORO1B"      "TBC1D4"      "EIF4G3"   "BATF"        "SIRPG"       "FLOT2"       "CERKL"       "OAS2"    "SEC16A"      "MSI2" 
CD4_M3 <- CD4_clusters_top50$gene[c(117,139,145,149)]
# "RP11-464D20.2" "RP11-277L2.3"  "SNRPB2"        "AC120194.1" 
CD4_M4 <- CD4_clusters_top50$gene[c(172,173,174,181,184,187,189,190,192,194,195:200)]
# "NME1-NME2"     "RP11-771F20.1" "EEF1A1P7"      "RPS7P10"      "RPL7AP11"      "LEF1"          "NOSIP"         "AC010468.1"   "SBDSP1"        "TMEM230"       "C1QBP"         "AD000092.3" "RP11-252A24.7" "RPL29P11"      "RP11-282K24.1" "ATP6AP2"
CD4_M5 <- CD4_clusters_top50$gene[c(205,212,213,223,224,225,227,228,229,231,232,234,236,237,239,241:244,246:250)]
# "ELL2"          "NR4A1"         "NR4A3"         "CTB-119C2.1"   "DCTN6"         "SERTAD1"       "DNAJA1P3"      "DDIT3"     "HSPA8P5"       "RP11-290L1.3"  "PRMT10"        "NAMPTL"  "RP11-1100L3.8" "PBX4"          "HSPA8P8"       "DUSP10"       "AC002306.1"    "BCL2L11"       "SLC2A14"       "FTH1P23"    "MTHFD2"        "NEU1"          "ETV3"          "JMJD6" 
p1 <- VlnPlot(CD4_sce1, features = CD4_M1, pt.size=0, ncol=3)
ggsave("356_CD8subType1_selectgenes_VlnPlot.png", p1, width=20 ,height=10)
p2 <- VlnPlot(CD4_sce1, features = CD4_M2, pt.size=0, ncol=3)
ggsave("357_CD8subType2_selectgenes_VlnPlot.png", p2, width=20 ,height=10)
p3 <- VlnPlot(CD4_sce1, features = CD4_M3, pt.size=0, ncol=3)
ggsave("358_CD8subType3_selectgenes_VlnPlot.png", p3, width=20 ,height=10)
p4 <- VlnPlot(CD4_sce1, features = CD4_M4, pt.size=0, ncol=3)
ggsave("359_CD8subType4_selectgenes_VlnPlot.png", p4, width=20 ,height=10)
p5 <- VlnPlot(CD4_sce1, features = CD4_M5, pt.size=0, ncol=3)
ggsave("360_CD8subType5_selectgenes_VlnPlot.png", p5, width=20 ,height=10)

#### Find differentially expressed genes for self-defined clusters  #########

#### (1,2)vs(3,4,5)#######################
CD4ClusterData <- CD4_sce1@meta.data$CD4_Cluster
RvsNR_CD4Marker <- c()
for (i in 1:length(CD4ClusterData)){
	if(CD4ClusterData[i] %in% c(1,2)){
		RvsNR_CD4Marker[i]="NR"
	}else RvsNR_CD4Marker[i]="R"
}

CD4_sce1<- AddMetaData(CD4_sce1, RvsNR_CD4Marker,col.name = "RvsNR_CD4Marker")

Idents(CD4_sce1)="RvsNR_CD4Marker"
CD8_RvsNR_markers <- FindAllMarkers(object = CD4_sce1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
# Check number of marker genes per cluster
table(CD8_RvsNR_markers $cluster)
 

top20 <- CD8_RvsNR_markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)

startNum <- 1
for (i in 1:2){
	
	endNum <- 20*i
	select_genes <- top20$gene[startNum:endNum]
	startNum=endNum+1
	p1 <- VlnPlot(CD4_sce1, features = select_genes, pt.size=0, ncol=4)
	Fname <- paste(360+i,i,"CD4subType_selectgenes_VlnPlot.png",sep="_")
	ggsave(Fname, p1, width=30,height=20)
}
select_genes <- top20$gene[c(1,3,4,6,22,26,27,28)]
# vlnplot display
p1 <- VlnPlot(CD4_sce1, features = select_genes, pt.size=0, ncol=3)
ggsave("363_CD4_RvsNR_markers_VlnPlot.png", p1, width=15 ,height=10)

CD4_NR_Marker <- top20$gene[c(1,3,4,6,9:12,14,15,17:19)]
CD4_R_Marker <- top20$gene[c(22,26,27,28,29,30,32:39)]
# First 13 R group: "NKG7","SLFN12L","CD84","KLRK1","HERC2P2","STAG3L2","TIGIT","APOBEC3G","RBCK1","LRBA","GRAP2","DOCK11","TTN"
# Last 14 NR group: "CDKN1A","BEST1","FAM177A1","LMNA","RASGEF1B","FTH1P23","FTH1P20","RGCC","KDM6B","CSRNP1","CHMP1B","PRMT10","NR4A1","IDI1"


#### (5),(4,3)vs(1)vs(2)#####################
RBvsRPvsNRBvsNRP_CD4Marker <- c()
for (i in 1:length(CD4ClusterData)){
	if(CD4ClusterData[i] %in% c(4,3)){
		RBvsRPvsNRBvsNRP_CD4Marker[i]="RB"
	}
	else if (CD4ClusterData[i] == 1){
		RBvsRPvsNRBvsNRP_CD4Marker[i]="NRB"
	}
	else if (CD4ClusterData[i] == 2){
		RBvsRPvsNRBvsNRP_CD4Marker[i]="NRP"
	}else RBvsRPvsNRBvsNRP_CD4Marker[i]="RP"
}

CD4_sce1<- AddMetaData(CD4_sce1, RBvsRPvsNRBvsNRP_CD4Marker,col.name = "RBvsRPvsNRBvsNRP_CD4Marker")

Idents(CD4_sce1)="RBvsRPvsNRBvsNRP_CD4Marker"
CD4_RBvsRPvsNRBvsNRP_markers <- FindAllMarkers(object = CD4_sce1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
# Check number of marker genes per cluster
table(CD4_RBvsRPvsNRBvsNRP_markers $cluster)
 


top20 <- CD4_RBvsRPvsNRBvsNRP_markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)

startNum <- 1
for (i in 1:4){
	
	endNum <- 20*i
	select_genes <- top20$gene[startNum:endNum]
	startNum=endNum+1
	p1 <- VlnPlot(CD4_sce1, features = select_genes, pt.size=0, ncol=4)
	Fname <- paste(363+i,i,"CD4subType_selectgenes_VlnPlot.png",sep="_")
	ggsave(Fname, p1, width=30,height=20)
}

select_genes <- top20$gene[c(4,6,9,21,22,23,52,53,56,61,64,67)]

# vlnplot display
p1 <- VlnPlot(CD4_sce1, features = select_genes, pt.size=0, ncol=3)
ggsave("368_RBvsRPvsNRBvsNRP_CD4Marker_VlnPlot.png", p1, width=15 ,height=10)

CD4_NRB_Marker <- top20$gene[c(4,6,9,13,15,20)]
CD4_NRP_Marker <- top20$gene[c(21,22,23,24,25,26,29,31,36,38,39)]
CD4_RB_Marker <- top20$gene[c(52,53,56,57,59)]
CD4_RP_Marker <- top20$gene[c(61,64,67,68,69,71,72,73:80)]
# First is NRB group: "CD8B","LINC00623","INO80D","CYTH4","YIPF4","NBPF20"
# Second is NRP group: "CTLA4","TIGIT","HLA-DRB1","TOX","PDCD1","HLA-DRB5","MIR4435-1HG" ,"CORO1B","TBC1D4","EIF4G3","BATF" 
# Third is RB group: "RPL24P4","AC079250.1","RPL21P28","RPS4XP13","RP11-771F20.1"
# Fourth is RP group: "CDKN1A","ELL2","NR4A1","NR4A3","RGCC","CTB-119C2.1","DCTN6","SERTAD1","DNAJA1P3","DDIT3","HSPA8P5","RP11-290L1.3","NAMPTL","RP11-1100L3.8" ,"HSPA8P8"
