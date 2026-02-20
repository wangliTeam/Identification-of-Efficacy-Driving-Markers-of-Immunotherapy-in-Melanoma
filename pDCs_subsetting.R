# Extract pDCs cell expression profile and metadata information
pDCs_sce1 = New_All_Data_2[,New_All_Data_2@meta.data$cellTypes=="pDCs"]
# All_Data pDCs cells expression profile data
pDCs_cells_Data <- pDCs_sce1@assays$RNA@data
# 55727 genes   225 samples
write.csv(pDCs_cells_Data, "pDCs_cells_Data.csv")



# Extract marker genes related to R for pDCs cells
pDCscells_Marker <- NoR_All_Data_1.markers[NoR_All_Data_1.markers$cluster == "pDCs",]

######## Select top 20%  #357 genes  ############
pDCscells_Marker_top_genes <- pDCscells_Marker %>% top_frac(.2,wt=avg_log2FC)  
########### Extract genes ########
select_pDCs_top_genes <- pDCscells_Marker_top_genes$gene

p1 <- VlnPlot (NoR_All_Data_1 , features =select_pDCs_top_genes[1:9]  , pt.size=0, ncol=3)
ggsave("401_select_pDCs_top_genes_VlnPlot.png", p1, width=10 ,height=8)
select_pDCs_top_genes[1:9]
 



# Calculate mean expression of pDCs cell-specific genes (in pDCs cells) in NoR_All_Data_1
NoR_All_ExpData <- NoR_All_Data_1@assays$RNA@data
select_pDCs_top_genes_ExpData <- NoR_All_ExpData[rownames(NoR_All_ExpData) %in% select_pDCs_top_genes,NoR_All_Data_1@meta.data$cellTypes=="pDCs"]
#385 genes, 199 samples

select_pDCs_top_genes_meanExp<- apply(select_pDCs_top_genes_ExpData,1,mean)
select_pDCs_top_genes_meanExp
write.table(select_pDCs_top_genes_meanExp,"select_pDCs_top_genes_meanExp.txt") 

############### Discussion of All_Data pDCs cells expression profile ##############
# 2. Filter genes (using median absolute deviation, MAD)
mads <- apply(pDCs_cells_Data,1,mad) # MAD measure
df <- pDCs_cells_Data[rev(order(mads))[1:5000],] # Extract top 5000 genes
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
                                 title="pDCs_Concluster", # Result save path
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

optK = 4  # Ideal K value is 4

icl = calcICL(results,
              title="pDCs_Concluster",
              plot="pdf")   # icl contains two dataframes: "clusterConsensus" cluster consensus scores;
                            #"itemConsensus" cluster corresponding samples and sample consensus scores 
icl[["clusterConsensus"]][1:5,]
icl[["itemConsensus"]][1:5,]
 



itemCon_pDCs_cells <- icl[["itemConsensus"]]
cluster_item <- subset(itemCon_pDCs_cells,k==4,select = c("cluster","item","itemConsensus"))
cluster_item$item[1:20]
 



cluster_item_sorted <- cluster_item[order(cluster_item$item, cluster_item$itemConsensus,decreasing=T),]

cluster_item_res=c()
for (i in 1:dim(cluster_item_sorted)[1]){
  if (i%%optK == 1){
  cluster_item_res=rbind(cluster_item_res,cluster_item_sorted[i,])
  }
}
head(cluster_item_res)
 



# Object is pDCs_sce1
# pDCs_Data <- pDCs_sce1@assays$RNA@data
attach(cluster_item_res)
pDCs_Cluster <- c()
for (i in 1:dim(cluster_item_res)[1]){
	pDCs_Cluster[rownames(pDCs_sce1@meta.data)==item[i]] <- cluster[i]
}
pDCs_Cluster   # Cluster corresponding to pDCs_sce1 sample order
detach(cluster_item_res)
 



pDCs_sce1<- AddMetaData(pDCs_sce1, pDCs_Cluster,col.name = "pDCs_Cluster")

pDCs_C1_Sample <- colnames(pDCs_cells_Data[,pDCs_sce1@meta.data$pDCs_Cluster == 1])
pDCs_C2_Sample <- colnames(pDCs_cells_Data[,pDCs_sce1@meta.data$pDCs_Cluster == 2])
pDCs_C3_Sample <- colnames(pDCs_cells_Data[,pDCs_sce1@meta.data$pDCs_Cluster == 3])
pDCs_C4_Sample <- colnames(pDCs_cells_Data[,pDCs_sce1@meta.data$pDCs_Cluster == 4])
length(pDCs_C1_Sample)   #127
length(pDCs_C2_Sample)   #82
length(pDCs_C3_Sample)   #11
length(pDCs_C4_Sample)   #5

All_ExpData <- New_All_Data_2@assays$RNA@data
pDCs_C1_select_genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% select_pDCs_top_genes,colnames(All_ExpData) %in% pDCs_C1_Sample]
dim(pDCs_C1_select_genes_ExpData)
#385 genes   127 samples
pDCs_C2_select_genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% select_pDCs_top_genes,colnames(All_ExpData) %in% pDCs_C2_Sample]
dim(pDCs_C2_select_genes_ExpData)
#385 genes   82 samples
pDCs_C3_select_genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% select_pDCs_top_genes,colnames(All_ExpData) %in% pDCs_C3_Sample]
dim(pDCs_C3_select_genes_ExpData)
#385 genes   11 samples
pDCs_C4_select_genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% select_pDCs_top_genes,colnames(All_ExpData) %in% pDCs_C4_Sample]
dim(pDCs_C4_select_genes_ExpData)
#385 genes   5 samples



pDCs_C1_select_genes_meanExp<- apply(pDCs_C1_select_genes_ExpData,1,mean)
pDCs_C2_select_genes_meanExp<- apply(pDCs_C2_select_genes_ExpData,1,mean)
pDCs_C3_select_genes_meanExp<- apply(pDCs_C3_select_genes_ExpData,1,mean)
pDCs_C4_select_genes_meanExp<- apply(pDCs_C4_select_genes_ExpData,1,mean)


######### Calculate correlation matrix ##################
pDCssubType_cor <- t(rbind(select_pDCs_top_genes_meanExp,pDCs_C1_select_genes_meanExp,pDCs_C2_select_genes_meanExp,pDCs_C3_select_genes_meanExp,pDCs_C4_select_genes_meanExp))
colnames(pDCssubType_cor)=c("NoR_pDCs","pDCs_c1","pDCs_c2","pDCs_c3","pDCs_c4")


res <- cor(pDCssubType_cor)
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
  pivot_longer(2:6) -> corrdf # Convert data from wide to long format

# p-value matrix
corrlist$P %>% 
  as_tibble() %>% 
  mutate(v = colnames(.)) %>% 
  select(v, everything()) %>% 
  pivot_longer(2:6) %>%  # 2:6 means I have 37 variables
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
 


# For the various pDCs cell subtypes divided above, we need to screen their markers
# Find differentially expressed genes for self-defined cell groups
Idents(pDCs_sce1)="pDCs_Cluster"
pDCs_celltype_markers <- FindAllMarkers(object = pDCs_sce1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
# Check number of marker genes per cluster
table(pDCs_celltype_markers $cluster)
write.csv(pDCs_celltype_markers,file="pDCs_celltype_markers.csv")
pDCs_celltype_markers<- read.csv("pDCs_celltype_markers.csv")
library(dplyr) 
 




### Select major cell type-specific genes ##############
pDCscells_Marker <- NoR_All_Data_1.markers[NoR_All_Data_1.markers$cluster == "pDCs",]
######## Select top 20%  #385 genes  ############
pDCscells_Marker_top_genes <- pDCscells_Marker %>% top_frac(.2,wt=avg_log2FC)  
write.csv(pDCscells_Marker_top_genes ,"pDCscells_Marker_top_genes.csv")
########### Extract genes ########
select_pDCs_top_genes <- pDCscells_Marker_top_genes$gene

# Select subtype specific genes ######
select_genes_pDCs_c1 <- (pDCs_celltype_markers[pDCs_celltype_markers$cluster==1,] %>% top_n(n=30,wt=avg_log2FC))$gene     #top30
select_genes_pDCs_c2 <- (pDCs_celltype_markers[pDCs_celltype_markers$cluster==2,] %>% top_n(n=30,wt=avg_log2FC))$gene  #top30

res <- c()
pDCs_inter_c1_genes <- intersect(select_genes_pDCs_c1,select_pDCs_top_genes)
res[1]=length(pDCs_inter_c1_genes)
pDCs_inter_c2_genes <- intersect(select_genes_pDCs_c2,select_pDCs_top_genes)
res[2]=length(pDCs_inter_c2_genes)
res
#top20-c1:
#top20-c2:"RP1-313I6.12"

#top30-c1:"BCL7A","RUFY4","AC006978.6"
#top30-c2:"SLC32A1","RP1-313I6.12"
#Total 5 genes: "BCL7A","RUFY4","AC006978.6","SLC32A1","RP1-313I6.12"
#pDCs-top20%-DiffGenes=385, pDCs-c1-DiffGenes=30, pDCs-c2-DiffGenes=30, pDCs-top20%-DiffGenes & pDCs-c1-DiffGenes =3, pDCs-top20%-DiffGenes & pDCs-c2-DiffGenes =2, pDCs-top20%-DiffGenes & pDCs-c1-DiffGenes & pDCs-c1-DiffGenes =5
# Visualization website: https://asntech.shinyapps.io/intervene/
# Parameters: color: #E09A4A  #317EAB
#           Plot width: 440, Plot height=300, ratio=0.5, all text 1.5
# Connecting point size=4, Connecting point size=1
