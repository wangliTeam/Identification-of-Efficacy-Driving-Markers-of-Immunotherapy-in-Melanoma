# Extract T cell expression profile and metadata information
T_sce1 = New_All_Data_2[,New_All_Data_2@meta.data$cellTypes=="T cells"]
# All_Data T cells expression profile data
T_cells_Data <- T_sce1@assays$RNA@data
# 55727 genes   362 samples
write.csv(T_cells_Data, "T_cells_Data.csv")


# Extract marker genes related to R for T cells
Tcells_Marker <- NoR_All_Data_1.markers[NoR_All_Data_1.markers$cluster == "T cells",]

Tcells_Marker_top20genes <- Tcells_Marker %>% top_n(n=20,wt=avg_log2FC)

select_T_top20genes <- Tcells_Marker_top20genes$gene
p1 <- VlnPlot (NoR_All_Data_1 , features =select_T_top20genes  , pt.size=0, ncol=3)
ggsave("407_select_T_top20genes_VlnPlot.png", p1, width=30 ,height=30)
 

# Calculate mean expression of T cell-specific genes (in T cells) in NoR_All_Data_1
NoR_All_ExpData <- NoR_All_Data_1@assays$RNA@data
select_T_top20genes_ExpData <- NoR_All_ExpData[rownames(NoR_All_ExpData) %in% c("TRDV1","TRDC","AE000661.37","TRGV8","RTKN2","SPC24","GTSE1","HIST1H3G","AURKB","CDCA8","ASPM","UBE2C","CDCA5","KIFC1","BIRC5","CCNB2","KIF18B","POLQ","CDC20","NCAPG"),NoR_All_Data_1@meta.data$cellTypes=="T cells"]
# 20 genes, 359 samples

select_T_top20genes_meanExp<- apply(select_T_top20genes_ExpData,1,mean)
select_T_top20genes_meanExp
write.table(select_T_top20genes_meanExp,"select_T_top20genes_meanExp.txt")
 

############### Discussion of All_Data T cells expression profile ##############
# 2. Filter genes (using median absolute deviation, MAD)
mads <- apply(T_cells_Data,1,mad) # MAD measure
df <- T_cells_Data[rev(order(mads))[1:5000],] # Extract top 5000 genes
# 3. Normalization
df <-  as.matrix(sweep(df,1, apply(df,1,median,na.rm=T))) # Subtract row median, default is subtraction

# 4. Run ConsensusClusterPlus
library(ConsensusClusterPlus)
maxK <-  6 # Try a K value
results <-  ConsensusClusterPlus(df,
                                 maxK = maxK,
                                 reps = 1000,              # Number of subsamples (usually 1000 or more)
                                 pItem = 0.8,              # Subsample ratio
                                 pFeature = 1,
                                 clusterAlg = "pam",       # Clustering algorithm
                                 distance="pearson",       # Distance calculation method
                                 title="T_Concluster", # Result save path
                                 innerLinkage="complete",  # Do not recommend using default "average" method here
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

optK = 5  # Ideal K value is 5

icl = calcICL(results,
              title="T_Concluster",
              plot="pdf")   # icl contains two dataframes: "clusterConsensus" cluster consensus scores;
                            #"itemConsensus" cluster corresponding samples and sample consensus scores 
icl[["clusterConsensus"]][1:5,]
icl[["itemConsensus"]][1:5,]
 


itemCon_T_cells <- icl[["itemConsensus"]]
cluster_item <- subset(itemCon_T_cells,k==5,select = c("cluster","item","itemConsensus"))
cluster_item$item[1:20]
 


cluster_item_sorted <- cluster_item[order(cluster_item$item, cluster_item$itemConsensus,decreasing=T),]

cluster_item_res=c()
for (i in 1:dim(cluster_item_sorted)[1]){
  if (i%%optK == 1){
  cluster_item_res=rbind(cluster_item_res,cluster_item_sorted[i,])
  }
}
head(cluster_item_res)
 


# Object is T_sce1
# T_Data <- T_sce1@assays$RNA@data
attach(cluster_item_res)
T_Cluster <- c()
for (i in 1:dim(cluster_item_res)[1]){
	T_Cluster[rownames(T_sce1@meta.data)==item[i]] <- cluster[i]
}
T_Cluster   # Cluster corresponding to T_sce1 sample order
detach(cluster_item_res)
 


T_sce1<- AddMetaData(T_sce1, T_Cluster,col.name = "T_Cluster")
save(T_sce1,file='T_sce1.Rdata')

T_C1_Sample <- colnames(T_cells_Data[,T_sce1@meta.data$T_Cluster == 1])
T_C2_Sample <- colnames(T_cells_Data[,T_sce1@meta.data$T_Cluster == 2])
T_C3_Sample <- colnames(T_cells_Data[,T_sce1@meta.data$T_Cluster == 3])
T_C4_Sample <- colnames(T_cells_Data[,T_sce1@meta.data$T_Cluster == 4])
T_C5_Sample <- colnames(T_cells_Data[,T_sce1@meta.data$T_Cluster == 5])
length(T_C1_Sample)   #84
length(T_C2_Sample)   #71
length(T_C3_Sample)   #55
length(T_C4_Sample)   #88
length(T_C5_Sample)   #64

All_ExpData <- New_All_Data_2@assays$RNA@data
T_C1_select20genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% c("TRDV1","TRDC","AE000661.37","TRGV8","RTKN2","SPC24","GTSE1","HIST1H3G","AURKB","CDCA8","ASPM","UBE2C","CDCA5","KIFC1","BIRC5","CCNB2","KIF18B","POLQ","CDC20","NCAPG"),colnames(All_ExpData) %in% T_C1_Sample]
dim(T_C1_select20genes_ExpData)
#20 genes   84 samples
T_C2_select20genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% c("TRDV1","TRDC","AE000661.37","TRGV8","RTKN2","SPC24","GTSE1","HIST1H3G","AURKB","CDCA8","ASPM","UBE2C","CDCA5","KIFC1","BIRC5","CCNB2","KIF18B","POLQ","CDC20","NCAPG"),colnames(All_ExpData) %in% T_C2_Sample]
dim(T_C2_select20genes_ExpData)
#20 genes   71 samples
T_C3_select20genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% c("TRDV1","TRDC","AE000661.37","TRGV8","RTKN2","SPC24","GTSE1","HIST1H3G","AURKB","CDCA8","ASPM","UBE2C","CDCA5","KIFC1","BIRC5","CCNB2","KIF18B","POLQ","CDC20","NCAPG"),colnames(All_ExpData) %in% T_C3_Sample]
dim(T_C3_select20genes_ExpData)
#20 genes   55 samples
T_C4_select20genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% c("TRDV1","TRDC","AE000661.37","TRGV8","RTKN2","SPC24","GTSE1","HIST1H3G","AURKB","CDCA8","ASPM","UBE2C","CDCA5","KIFC1","BIRC5","CCNB2","KIF18B","POLQ","CDC20","NCAPG"),colnames(All_ExpData) %in% T_C4_Sample]
dim(T_C4_select20genes_ExpData)
#20 genes   88 samples
T_C5_select20genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% c("TRDV1","TRDC","AE000661.37","TRGV8","RTKN2","SPC24","GTSE1","HIST1H3G","AURKB","CDCA8","ASPM","UBE2C","CDCA5","KIFC1","BIRC5","CCNB2","KIF18B","POLQ","CDC20","NCAPG"),colnames(All_ExpData) %in% T_C5_Sample]
dim(T_C5_select20genes_ExpData)
#20 genes   64 samples


T_C1_select20genes_meanExp<- apply(T_C1_select20genes_ExpData,1,mean)
T_C2_select20genes_meanExp<- apply(T_C2_select20genes_ExpData,1,mean)
T_C3_select20genes_meanExp<- apply(T_C3_select20genes_ExpData,1,mean)
T_C4_select20genes_meanExp<- apply(T_C4_select20genes_ExpData,1,mean)
T_C5_select20genes_meanExp<- apply(T_C5_select20genes_ExpData,1,mean)

######### Calculate correlation matrix ##################
TsubType_cor <- t(rbind(select_T_top20genes_meanExp,T_C1_select20genes_meanExp,T_C2_select20genes_meanExp,T_C3_select20genes_meanExp,T_C4_select20genes_meanExp,T_C5_select20genes_meanExp))
colnames(TsubType_cor)=c("NoR_Tcells","T_c1","T_c2","T_c3","T_c4","T_c5")
TsubType_cor
 


res <- cor(TsubType_cor)
res2 <- signif(res, digits = 3)
res2
library(Hmisc)# Load package
 



Hmisc::rcorr(as.matrix(res2), type = "pearson") -> corrlist
# Note: Need to install Hmisc package first, otherwise it will throw an error, same below. Because I'm lazy, I won't write automated code. If many people ask me to write, maybe I will. Put the matrix or dataframe that needs correlation calculation inside as.matrix()

library(tidyr)
# Correlation coefficient matrix
corrlist$r %>%  # Extract r matrix
  as_tibble() %>%  # Set to tibble format
  mutate(v = colnames(.)) %>%  # Add column names vector as a new vector to the column
  select(v, everything()) %>%  # Select all tables
  pivot_longer(2:7) -> corrdf # Convert data from wide to long format

# p-value matrix
corrlist$P %>% 
  as_tibble() %>% 
  mutate(v = colnames(.)) %>% 
  select(v, everything()) %>% 
  pivot_longer(2:7) %>%  # 2:7 means I have 7 variables
  mutate(label = case_when(  # Set label, and add judgment, when P value meets specific conditions display "\n" plus specific number of asterisks
    is.na(value) ~ " ", # NA values assigned to space
    value <= 0.001 ~ "\n***", # P<0.001 display newline plus three asterisks
    between(value, 0.001, 0.01) ~ "\n**", # P 0.001-0.01 display newline plus two asterisks
    between(value, 0.01, 0.05) ~ "\n*", # P 0.01-0.05 display newline plus one asterisk
    T ~ ""
  )) -> pdf
corrdf %>% 
  left_join(pdf, by = c("v", "name")) %>% # Merge r value and p matrix by v and name fields
  rename(corr = value.x, p = value.y) %>% # Rename value.x to corr, value.y to p in merged dataframe
  mutate(corr = signif(corr, digits = 3)) -> corrdf # Change decimal places to two
corrdf
 


library(grDevices)
#windowsFonts("Arial" = windowsFont("Arial")) # Set font to prevent error in the following code
corrdf %>% 
  mutate(v = forcats::fct_reorder(v, corr), # Reorder v and name
         name = forcats::fct_reorder(name, corr)) ->corrdf2
corrdf2
 


write.csv(corrdf2,"corrdf2_T.csv")
corrdf2 <- read.csv("corrdf2_T.csv",row.names=1)
#library(dplyr)  Pipe %>% package

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
 



# For the various γδT cell subtypes divided above, we need to screen their markers
# Find differentially expressed genes for self-defined cell groups
Idents(T_sce1)="T_Cluster"
T_celltype_markers <- FindAllMarkers(object = T_sce1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
# Check number of marker genes per cluster
table(T_celltype_markers $cluster)
#write.csv(T_celltype_markers,file="T_celltype_markers.csv")
T_celltype_markers<- read.csv("T_celltype_markers.csv")
library(dplyr)
 


### Select major cell type-specific genes ##############
Tcells_Marker <- NoR_All_Data_1.markers[NoR_All_Data_1.markers$cluster == "γδT cells",]
######## Select top 20%  #280 genes  ############
Tcells_Marker_top_genes <- Tcells_Marker %>% top_frac(.2,wt=avg_log2FC)  
write.csv(Tcells_Marker_top_genes, "Tcells_Marker_top_genes.csv")
########### Extract genes ########
select_T_top_genes <- Tcells_Marker_top_genes$gene

# Select subtype-specific genes ######
select_genes_T_c2 <- (T_celltype_markers[T_celltype_markers$cluster==2,] %>% top_n(n=30,wt=avg_log2FC))$gene  #top30
select_genes_T_c3 <- (T_celltype_markers[T_celltype_markers$cluster==3,] %>% top_n(n=30,wt=avg_log2FC))$gene    #top30
select_genes_T_c4 <- (T_celltype_markers[T_celltype_markers$cluster==4,] %>% top_n(n=30,wt=avg_log2FC))$gene
select_genes_T_c5 <- (T_celltype_markers[T_celltype_markers$cluster==5,] %>% top_n(n=30,wt=avg_log2FC))$gene     #top30

res <- c()
T_inter_c2_genes <- intersect(select_genes_T_c2,select_T_top_genes)
res[1]=length(T_inter_c2_genes)
T_inter_c3_genes <- intersect(select_genes_T_c3,select_T_top_genes)
res[2]=length(T_inter_c3_genes)
T_inter_c4_genes <- intersect(select_genes_T_c4,select_T_top_genes)
res[3]=length(T_inter_c4_genes)
T_inter_c5_genes <- intersect(select_genes_T_c5,select_T_top_genes)
res[4]=length(T_inter_c5_genes)

res
#top30-c2:"LYAR"  "ACOT7"
#top30-c3:"KIF20B"    "ARHGAP11A"
#top30-c4:"EPAS1" "TRGV5" "RGPD2"
#top30-c5:
#Total 7 genes: "LYAR" , "ACOT7","KIF20B"   , "ARHGAP11A","EPAS1", "TRGV5" ,"RGPD2"
#γδTcells-top20%-DiffGenes=280, γδTcells-c2-DiffGenes=30, γδTcells-c3-DiffGenes=30, γδTcells-c4-DiffGenes=30, γδTcells-c5-DiffGenes=30, γδTcells-top20%-DiffGenes&γδTcells-c2-DiffGenes =2, γδTcells-top20%-DiffGenes&γδTcells-c3-DiffGenes =2, γδTcells-top20%-DiffGenes&γδTcells-c4-DiffGenes =3, γδTcells-top20%-DiffGenes&γδTcells-c2-DiffGenes &γδTcells-c3-DiffGenes &γδTcells-c4-DiffGenes&γδTcells-c5-DiffGenes =7
# Visualization website: https://asntech.shinyapps.io/intervene/
# Parameters: color: #E09A4A  #317EAB
#           Plot width: 440, Plot height=300, ratio=0.5, all text 1.5
# Connecting point size=4, Connecting point size=1
