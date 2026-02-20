# Extract B cell expression profile and metadata information
cells <- c("B cells","CD4 T cells","pDCs","Regulatory T cells","T cells")
B_sce1 = New_All_Data_2[,New_All_Data_2@meta.data$cellTypes=="B cells"]
# All_Data B cells expression profile data
B_cells_Data <- B_sce1@assays$RNA@data
# 55727 genes, 696 samples
write.csv(B_cells_Data, "B_cells_Data.csv")


# Extract B cell marker genes related to R
Bcells_Marker <- Re_All_Data_1.markers[Re_All_Data_1.markers$cluster == "B cells",]

######## Select top 20%  #109 genes  ############
Bcells_Marker_top_genes <- Bcells_Marker %>% top_frac(.2,wt=avg_log2FC)  
########### Gene extraction ########
select_B_top_genes <- Bcells_Marker_top_genes$gene

p1 <- VlnPlot (Re_All_Data_1 , features =select_B_top_genes[1:9]  , pt.size=0, ncol=3)
ggsave("387_select_B_top__precent_genes_VlnPlot.png", p1, width=10 ,height=8)
select_B_top_genes[1:9]
 

# Calculate mean expression of B cell specific genes in Re_All_Data_1
Re_All_ExpData <- Re_All_Data_1@assays$RNA@data
select_B_top100genes_ExpData <- Re_All_ExpData[rownames(Re_All_ExpData) %in% select_B_top100genes,Re_All_Data_1@meta.data$cellTypes=="B cells"]
# 100 genes, 384 samples

select_B_top_genes_meanExp<- apply(select_B_top_genes_ExpData,1,mean)
select_B_top_genes_meanExp
write.table(select_B_top_genes_meanExp,"select_B_top_genes_meanExp.txt")
 





############### Analysis of All_Data B cell expression profile ##############
# 2. Gene filtering (by Median Absolute Deviation, MAD)
mads <- apply(B_cells_Data,1,mad) # MAD measure
df <- B_cells_Data[rev(order(mads))[1:5000],] # Extract top 5000 genes
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
                                 title="B_cells_Concluster", # Result save path
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

optK = 4  # Optimal K value is 4

icl = calcICL(results,
              title="B_cells_Concluster",
              plot="pdf")   # icl contains two dataframes: "clusterConsensus" cluster consistency scores;
                            # "itemConsensus" cluster corresponding samples and sample consistency scores 
icl[["clusterConsensus"]][1:5,]
icl[["itemConsensus"]][1:5,]

itemCon_B_cells <- icl[["itemConsensus"]]
cluster_item <- subset(itemCon_B_cells,k==4,select = c("cluster","item","itemConsensus"))
cluster_item$item[1:20]
 
cluster_item_Bcells_sorted <- cluster_item_Bcells[order(cluster_item_Bcells$item, cluster_item_Bcells$itemConsensus,decreasing=T),]

cluster_item_res=c()
for (i in 1:dim(cluster_item_sorted)[1]){
  if (i%%optK == 1){
  cluster_item_res=rbind(cluster_item_res,cluster_item_sorted[i,])
  }
}
head(cluster_item_res)
 



# Object is B_sce1
# B_Data <- B_sce1@assays$RNA@data
attach(cluster_item_res)
B_Cluster <- c()
for (i in 1:dim(cluster_item_res)[1]){
	B_Cluster[rownames(B_sce1@meta.data)==item[i]] <- cluster[i]
}
B_Cluster   # Cluster for B_sce1 sample order
detach(cluster_item_res)
 

B_sce1<- AddMetaData(B_sce1, B_Cluster,col.name = "B_Cluster")
 

B_C1_Sample <- colnames(B_cells_Data[,B_sce1@meta.data$B_Cluster == 1])
B_C2_Sample <- colnames(B_cells_Data[,B_sce1@meta.data$B_Cluster == 2])
B_C3_Sample <- colnames(B_cells_Data[,B_sce1@meta.data$B_Cluster == 3])
B_C4_Sample <- colnames(B_cells_Data[,B_sce1@meta.data$B_Cluster == 4])
length(B_C1_Sample)   #378
length(B_C2_Sample)   #167
length(B_C3_Sample)   #73
length(B_C4_Sample)   #78


All_ExpData <- New_All_Data_2@assays$RNA@data
B_C1_select_genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% select_B_top_genes,colnames(All_ExpData) %in% B_C1_Sample]
dim(B_C1_select_genes_ExpData)
# 109 genes, 378 samples
B_C2_select_genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% select_B_top_genes,colnames(All_ExpData) %in% B_C2_Sample]
dim(B_C2_select_genes_ExpData)
# 109 genes, 167 samples
B_C3_select_genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% select_B_top_genes,colnames(All_ExpData) %in% B_C3_Sample]
dim(B_C3_select_genes_ExpData)
# 109 genes, 73 samples
B_C4_select_genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% select_B_top_genes,colnames(All_ExpData) %in% B_C4_Sample]
dim(B_C4_select_genes_ExpData)
# 109 genes, 78 samples


B_C1_select_genes_meanExp<- apply(B_C1_select_genes_ExpData,1,mean)
B_C2_select_genes_meanExp<- apply(B_C2_select_genes_ExpData,1,mean)
B_C3_select_genes_meanExp<- apply(B_C3_select_genes_ExpData,1,mean)
B_C4_select_genes_meanExp<- apply(B_C4_select_genes_ExpData,1,mean)


######### Correlation matrix calculation ##################
BsubType_cor <- t(rbind(select_B_top_genes_meanExp,B_C1_select_genes_meanExp,B_C2_select_genes_meanExp,B_C3_select_genes_meanExp,B_C4_select_genes_meanExp))
colnames(BsubType_cor)=c("Re_Bcells","B_c1","B_c2","B_c3","B_c4")

res <- cor(BsubType_cor)
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
  pivot_longer(2:6) -> corrdf # Convert from wide to long data

# p value matrix
corrlist$P %>% 
  as_tibble() %>% 
  mutate(v = colnames(.)) %>% 
  select(v, everything()) %>% 
  pivot_longer(2:6) %>%  # 2:6 means I have 6 variables
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
write.csv(corrdf2,"corrdf2.csv")
corrdf2 <- read.csv("corrdf2.csv",row.names=1)
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
 


# For the various B cell subtypes defined above, we need to screen their markers
# Find differentially expressed genes for self-defined clusters
Idents(B_sce1)="B_Cluster"
B_celltype_markers <- FindAllMarkers(object = B_sce1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
# Check number of marker genes per cluster
table(B_celltype_markers $cluster)
write.csv(B_celltype_markers,file="B_celltype_markers.csv")
library(dplyr)
 

######################### top30 ############################
### Select major cell type specific genes ##############
Bcells_Marker <- Re_All_Data_1.markers[Re_All_Data_1.markers$cluster == "B cells",]
######## Select top 20%  #109 genes  ############
Bcells_Marker_top_genes <- Bcells_Marker %>% top_frac(.2,wt=avg_log2FC)  
########### Gene extraction ########
select_B_top_genes <- Bcells_Marker_top_genes$gene
write.csv(Bcells_Marker_top_genes, "Bcells_Marker_top20%_genes.csv")

# Select subtype c1/c2 specific genes ######
select_genes_B_c1 <- (B_celltype_markers[B_celltype_markers$cluster==1,] %>% top_n(n=30,wt=avg_log2FC))$gene     #top30
select_genes_B_c2 <- (B_celltype_markers[B_celltype_markers$cluster==2,] %>% top_n(n=30,wt=avg_log2FC))$gene  #top30

# Take intersection with B gene top 20% specific genes
res <- c()
B_inter_c1_genes <- intersect(select_genes_B_c1,select_B_top_genes) #4
res[1]=length(B_inter_c1_genes)
B_inter_c2_genes <- intersect(select_genes_B_c2,select_B_top_genes) #1
res[2]=length(B_inter_c2_genes)
res
# top30-c1:"IGHG1" "IGHGP" "IGHG2" "IGHG4" "SCIMP"
# top30-c2:"IGHV3-30"
B_inter_c1_2_genes <- union(B_inter_c1_genes,B_inter_c2_genes)
# Total 5 genes: "IGHG1"    "IGHGP"    "IGHG2"    "IGHG4"    "SCIMP"    "IGHV3-30"

# Visualization website  https://asntech.shinyapps.io/intervene/
# Parameters: color: #E09A4A  #317EAB
#         Plot width: 440, Plot height=300, ratio=0.5 , font all 1.5
