# Extract Treg cell expression profile and metadata information
Treg_sce1 = New_All_Data_2[,New_All_Data_2@meta.data$cellTypes=="Regulatory T cells"]
# All_Data Regulatory T cells expression profile data
Regulatory_T_cells_Data <- Treg_sce1@assays$RNA@data
# 55727 genes   696 samples
write.csv(Regulatory_T_cells_Data, "Regulatory_T_cells_Data.csv")



# Extract marker genes related to R for Treg cells
######## Select top 20%  #136 genes  ############
Tregcells_Marker_top_genes <- Tregcells_Marker %>% top_frac(.2,wt=avg_log2FC)  
########### Extract genes ########
select_Treg_top_genes <- Tregcells_Marker_top_genes$gene

p1 <- VlnPlot (NoR_All_Data_1 , features =select_Treg_top_genes[1:9]  , pt.size=0, ncol=3)
ggsave("406_select_Treg_top_genes_VlnPlot.png", p1, width=10 ,height=8)
select_Treg_top_genes[1:9]
 



# Calculate mean expression of Treg cell-specific genes (in Treg cells) in NoR_All_Data_1
NoR_All_ExpData <- NoR_All_Data_1@assays$RNA@data
select_Treg_top_genes_ExpData <- NoR_All_ExpData[rownames(NoR_All_ExpData) %in% select_Treg_top_genes,NoR_All_Data_1@meta.data$cellTypes=="Regulatory T cells"]
#136 genes, 1258 samples

select_Treg_top_genes_meanExp<- apply(select_Treg_top_genes_ExpData,1,mean)
write.table(select_Treg_top_genes_meanExp,"select_Treg_top_genes_meanExp.txt")
select_Treg_top_genes_meanExp
 





############### Discussion of All_Data Treg cells expression profile ##############
# 2. Filter genes (using median absolute deviation, MAD)
mads <- apply(Regulatory_T_cells_Data,1,mad) # MAD measure
df <- Regulatory_T_cells_Data[rev(order(mads))[1:5000],] # Extract top 5000 genes
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
                                 title="Regulatory_T_cells_Concluster", # Result save path
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
              title="Regulatory_T_cells_Concluster",
              plot="pdf")   # icl contains two dataframes: "clusterConsensus" cluster consensus scores;
                            #"itemConsensus" cluster corresponding samples and sample consensus scores 
icl[["clusterConsensus"]][1:5,]
icl[["itemConsensus"]][1:5,]
 



itemCon_Regulatory_T_cells <- icl[["itemConsensus"]]
cluster_item <- subset(itemCon_Regulatory_T_cells,k==5,select = c("cluster","item","itemConsensus"))
cluster_item$item[1:20]
 



cluster_item_sorted <- cluster_item[order(cluster_item$item, cluster_item$itemConsensus,decreasing=T),]

cluster_item_res=c()
for (i in 1:dim(cluster_item_sorted)[1]){
  if (i%%optK == 1){
  cluster_item_res=rbind(cluster_item_res,cluster_item_sorted[i,])
  }
}
head(cluster_item_res)
 



# Object is Treg_sce1
# Treg_Data <- Treg_sce1@assays$RNA@data
attach(cluster_item_res)
Treg_Cluster <- c()
for (i in 1:dim(cluster_item_res)[1]){
	Treg_Cluster[rownames(Treg_sce1@meta.data)==item[i]] <- cluster[i]
}
Treg_Cluster   # Cluster corresponding to Treg_sce1 sample order
detach(cluster_item_res)
 




Treg_sce1<- AddMetaData(Treg_sce1, Treg_Cluster,col.name = "Treg_Cluster")
save(Treg_sce1,file='Treg_sce1.Rdata')



Treg_C1_Sample <- colnames(Regulatory_T_cells_Data[,Treg_sce1@meta.data$Treg_Cluster == 1])
Treg_C2_Sample <- colnames(Regulatory_T_cells_Data[,Treg_sce1@meta.data$Treg_Cluster == 2])
Treg_C3_Sample <- colnames(Regulatory_T_cells_Data[,Treg_sce1@meta.data$Treg_Cluster == 3])
Treg_C4_Sample <- colnames(Regulatory_T_cells_Data[,Treg_sce1@meta.data$Treg_Cluster == 4])
Treg_C5_Sample <- colnames(Regulatory_T_cells_Data[,Treg_sce1@meta.data$Treg_Cluster == 5])

length(Treg_C1_Sample)   #266
length(Treg_C2_Sample)   #465
length(Treg_C3_Sample)   #220
length(Treg_C4_Sample)   #269
length(Treg_C5_Sample)   #229

All_ExpData <- New_All_Data_2@assays$RNA@data
Treg_C1_select_genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% select_Treg_top_genes,colnames(All_ExpData) %in% Treg_C1_Sample]
dim(Treg_C1_select_genes_ExpData)
#136 genes   266 samples
Treg_C2_select_genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% select_Treg_top_genes,colnames(All_ExpData) %in% Treg_C2_Sample]
dim(Treg_C2_select_genes_ExpData)
#136 genes   465 samples
Treg_C3_select_genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% select_Treg_top_genes,colnames(All_ExpData) %in% Treg_C3_Sample]
dim(Treg_C3_select_genes_ExpData)
#136 genes   220 samples
Treg_C4_select_genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% select_Treg_top_genes,colnames(All_ExpData) %in% Treg_C4_Sample]
dim(Treg_C4_select_genes_ExpData)
#136 genes   269 samples
Treg_C5_select_genes_ExpData <- All_ExpData[rownames(All_ExpData) %in% select_Treg_top_genes,colnames(All_ExpData) %in% Treg_C5_Sample]
dim(Treg_C5_select_genes_ExpData)
#136 genes   229 samples

Treg_C1_select_genes_meanExp<- apply(Treg_C1_select_genes_ExpData,1,mean)
Treg_C2_select_genes_meanExp<- apply(Treg_C2_select_genes_ExpData,1,mean)
Treg_C3_select_genes_meanExp<- apply(Treg_C3_select_genes_ExpData,1,mean)
Treg_C4_select_genes_meanExp<- apply(Treg_C4_select_genes_ExpData,1,mean)
Treg_C5_select_genes_meanExp<- apply(Treg_C5_select_genes_ExpData,1,mean)

######### Calculate correlation matrix ##################
TregsubType_cor <- t(rbind(select_Treg_top_genes_meanExp,Treg_C1_select_genes_meanExp,Treg_C2_select_genes_meanExp,Treg_C3_select_genes_meanExp,Treg_C4_select_genes_meanExp,Treg_C5_select_genes_meanExp))
colnames(TregsubType_cor)=c("NoR_Tregs","Treg_c1","Treg_c2","Treg_c3","Treg_c4","Treg_c5")

res <- cor(TregsubType_cor)
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
 



write.csv(corrdf2,"corrdf2_Treg.csv")
corrdf2 <- read.csv("corrdf2_Treg.csv",row.names=1)
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
 


# For the various Treg cell subtypes divided above, we need to screen their markers
# Find differentially expressed genes for self-defined cell groups
Idents(Treg_sce1)="Treg_Cluster"
Treg_celltype_markers <- FindAllMarkers(object = Treg_sce1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
# Check number of marker genes per cluster
table(Treg_celltype_markers $cluster)
write.csv(Treg_celltype_markers,file="Treg_celltype_markers.csv")
Treg_celltype_markers<- read.csv("Treg_celltype_markers.csv")
library(dplyr)
 



### Select major cell type-specific genes ##############
Tregcells_Marker <- NoR_All_Data_1.markers[NoR_All_Data_1.markers$cluster == "Regulatory T cells",]
######## Select top 20%  #136 genes  ############
Tregcells_Marker_top_genes <- Tregcells_Marker %>% top_frac(.2,wt=avg_log2FC)  
write.csv(Tregcells_Marker_top_genes, "Tregcells_Marker_top_genes.csv")
########### Extract genes ########
select_Treg_top_genes <- Tregcells_Marker_top_genes$gene

# Select subtype-specific genes ######
select_genes_Treg_c2 <- (Treg_celltype_markers[Treg_celltype_markers$cluster==2,] %>% top_n(n=30,wt=avg_log2FC))$gene  #top30
select_genes_Treg_c4 <- (Treg_celltype_markers[Treg_celltype_markers$cluster==4,] %>% top_n(n=30,wt=avg_log2FC))$gene  #top30

res <- c()
Treg_inter_c2_genes <- intersect(select_genes_Treg_c2,select_Treg_top_genes)
Treg_inter_c4_genes <- intersect(select_genes_Treg_c4,select_Treg_top_genes)
res[1]=length(Treg_inter_c2_genes)
res[2]=length(Treg_inter_c4_genes)
res
#top20-c2:"FOXP3" "IL2RA"
#top30-c2:"FOXP3","RAB11FIP1","IL2RA"
#top30-c4: "SUSD3"
#Total 4 genes: "FOXP3","RAB11FIP1","IL2RA","SUSD3"
#Tregs-top20%-DEGs=136, Tregs-c2-DEGs=30, Tregs-c4-DEGs=30,Tregs-top20%-DEGs & Tregs-c2-DEGs =3, Tregs-top20%-DEGs & Tregs-c4-DEGs =1, Tregs-top20%-DEGs & Tregs-c2-DEGs& Tregs-c4-DEGs =4
# Visualization website: https://asntech.shinyapps.io/intervene/
# Parameters: color: #E09A4A  #317EAB
#           Plot width: 440, Plot height=300, ratio=0.5, all text 1.5
# Connecting point size=4, Connecting point size=1
