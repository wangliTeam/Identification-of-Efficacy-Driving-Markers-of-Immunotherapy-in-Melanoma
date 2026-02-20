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
#Mycol_nejm= pal_nejm("default", alpha =0.7)(8)## Extract colors
Mycol_4 <- c("#A6CEE3" ,"#1F78B4", "#B2DF8A" ,"#33A02C" ,"#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00", "#CAB2D6" ,"#6A3D9A" ,"#FFFF99", "#B15928","#1B9E77", "#D95F02" ,"#7570B3", "#E7298A" ,"#66A61E" ,"#E6AB02", "#A6761D","#666666")
Mycol_5<-c("#F8766D","#CD9600","#7CAE00","#0CB702","#00C19A","#00B8E7","#8494FF","#C77CFF","#FF61CC")
# For the various B cell subtypes divided above, we need to screen their markers
# Find differentially expressed genes for self-defined cell groups
Idents(B_sce1)="B_Cluster"
B_celltype_markers <- FindAllMarkers(object = B_sce1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
# Check number of marker genes per cluster
table(B_celltype_markers $cluster)
#write.csv(B_celltype_markers,file="B_celltype_markers.csv")
B_celltype_markers <- read.csv("B_celltype_markers.csv")
library(dplyr)
 


######################### top30 ############################
### Select major cell type-specific genes ##############
Bcells_Marker <- Re_All_Data_1.markers[Re_All_Data_1.markers$cluster == "B cells",]
######## Select top 20%  #109 genes  ############
Bcells_Marker_top_genes <- Bcells_Marker %>% top_frac(.2,wt=avg_log2FC)  
########### Extract genes ########
select_B_top_genes <- Bcells_Marker_top_genes$gene
write.csv(Bcells_Marker_top_genes, "Bcells_Marker_top20%_genes.csv")

# Select subtype c1/c2 specific genes ######
select_genes_B_c1 <- (B_celltype_markers[B_celltype_markers$cluster==1,] %>% top_n(n=30,wt=avg_log2FC))$gene     #top30
select_genes_B_c2 <- (B_celltype_markers[B_celltype_markers$cluster==2,] %>% top_n(n=30,wt=avg_log2FC))$gene  #top30

# Take intersection with top 20% B cell specific genes
res <- c()
B_inter_c1_genes <- intersect(select_genes_B_c1,select_B_top_genes) #4
res[1]=length(B_inter_c1_genes)
B_inter_c2_genes <- intersect(select_genes_B_c2,select_B_top_genes) #1
res[2]=length(B_inter_c2_genes)
res
#top30-c1:"IGHG1" "IGHGP" "IGHG2" "IGHG4" "SCIMP"
#top30-c2:"IGHV3-30"
B_inter_c1_2_genes <- union(B_inter_c1_genes,B_inter_c2_genes)
#Total 5 genes: "IGHG1"    "IGHGP"    "IGHG2"    "IGHG4"    "SCIMP"    "IGHV3-30"

#B-top20%-DEGs=109,B-c1-DEGs=30, B-c2-DEGs=30, B-top20%-DEGs&B-c1-DEGs=4, B-top20%-DEGs&B-c2-DEGs=1, B-top20%-DEGs& B-c1-DEGs& B-c2-DEGs=5
# Visualization website: https://asntech.shinyapps.io/intervene/
# Parameters: color: #E09A4A  #317EAB
#           Plot width: 440, Plot height=300, ratio=0.5, all text 1.5
# Connecting point size=4, Connecting point size=1
   


# For the various CD4 cell subtypes divided above, we need to screen their markers
# Find differentially expressed genes for self-defined cell groups
Idents(CD4_sce1)="CD4_Cluster"
CD4_celltype_markers <- FindAllMarkers(object = CD4_sce1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
# Check number of marker genes per cluster
table(CD4_celltype_markers $cluster)
#write.csv(CD4_celltype_markers,file="CD4_celltype_markers.csv")
CD4_celltype_markers <- read.csv("CD4_celltype_markers.csv")
library(dplyr)
 


### Select major cell type-specific genes ##############
CD4cells_Marker <- Re_All_Data_1.markers[Re_All_Data_1.markers$cluster == "CD4 T cells",]
######## Select top 20%  #69 genes  ############
CD4cells_Marker_top_genes <- CD4cells_Marker %>% top_frac(.2,wt=avg_log2FC)  
dim(CD4cells_Marker_top_genes)
write.csv(CD4cells_Marker_top_genes, "CD4cells_Marker_top20%_genes.csv")

########### Extract genes ########
select_CD4_top_genes <- CD4cells_Marker_top_genes $gene

# Select subtype specific genes ######
select_genes_CD4_c1 <- (CD4_celltype_markers[CD4_celltype_markers$cluster==1,] %>% top_n(n=30,wt=avg_log2FC))$gene     #top30
write.csv(CD4cells_Marker_top_genes, "CD4cells_Marker_top20%_genes.csv")

res <- c()
CD4_inter_c1_genes <- intersect(select_genes_CD4_c1,select_CD4_top_genes)
res[1]=length(CD4_inter_c1_genes)
res

#top20-c1:"LINC00861"
#top30-c1:"LINC00861" "TMEM63A" 
#Total 2 genes: "LINC00861" "TMEM63A"

#CD4Tcells-top20%-DEGs=69, CD4-c1-DEGs=30, CD4Tcells-top20%-DEGs & CD4-c1-DEGs =2
# Visualization website: https://asntech.shinyapps.io/intervene/
# Parameters: color: #E09A4A  #317EAB
#           Plot width: 440, Plot height=300, ratio=0.5, all text 1.5
# Connecting point size=4, Connecting point size=1
   


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
#pDCs-top20%-DEGs=385, pDCs-c1-DEGs=30, pDCs-c2-DEGs=30, pDCs-top20%-DEGs & pDCs-c1-DEGs =3, pDCs-top20%-DEGs & pDCs-c2-DEGs =2, pDCs-top20%-DEGs & pDCs-c1-DEGs & pDCs-c1-DEGs =5
# Visualization website: https://asntech.shinyapps.io/intervene/
# Parameters: color: #E09A4A  #317EAB
#           Plot width: 440, Plot height=300, ratio=0.5, all text 1.5
# Connecting point size=4, Connecting point size=1




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

# Select subtype specific genes ######
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




# For the various γδT cells cell subtypes divided above, we need to screen their markers
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

# Select subtype specific genes ######
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
#γδTcells-top20%-DEGs=280, γδTcells-c2-DEGs=30, γδTcells-c3-DEGs=30, γδTcells-c4-DEGs=30, γδTcells-c5-DEGs=30, γδTcells-top20%-DEGs&γδTcells-c2-DEGs =2, γδTcells-top20%-DEGs&γδTcells-c3-DEGs =2, γδTcells-top20%-DEGs&γδTcells-c4-DEGs =3, γδTcells-top20%-DEGs&γδTcells-c2-DEGs &γδTcells-c3-DEGs &γδTcells-c4-DEGs&γδTcells-c5-DEGs =7
# Visualization website: https://asntech.shinyapps.io/intervene/
# Parameters: color: #E09A4A  #317EAB
#           Plot width: 440, Plot height=300, ratio=0.5, all text 1.5
# Connecting point size=4, Connecting point size=1
	



####### PPI analysis of marker gene sets between responder and non-responder groups –(STRING)############
####### Screened responder and non-responder markers ##############################
###↑↑↑ Could not find good PPI interactions in string, so proceed directly with subsequent analysis ######

############ Non-responder group PPI analysis  ########################
######NoR_Marker_top30###########
NoR_Marker_top30 <- c("BCL7A","RUFY4","AC006978.6","SLC32A1","RP1-313I6.12","FOXP3","RAB11FIP1","IL2RA","LYAR","ACOT7","KIF20B","ARHGAP11A","EPAS1","TRGV5","RGPD2","IL7R","TCF7","FTH1P2","FTH1P23","BEST1","FTH1P20","LMNA","SATB1","RPS4XP13","EEF1A1P16","MGAT4A")
  
################ Visualize in cytoscape and analyze important modules via Mcode ####
 
#Module1: IL7R,FOXP3,IL2RA
#Module2: EPAS1
#Mediator TCF7: TCF7,SATB1
######## Also perform KEGG enrichment analysis on non-responder genes - (DAVID) #########
#Pathways with P<0.05: hsa05200:Pathways in cancer, we display KEGG pathway for these genes
 


###### We found 4 non-responder marker genes enriched in the risk pathway: TCF7,IL2RA,IL7R,EPAS1
###### We want to further explore whether they are related to patient survival (subsequent)  #########
   


#Respectively IL2RA,TRGV5,BEST1
 
#Responder group also has one: LINC00861

########### Responder group gene co-expression analysis ###############
######### 
####Correlation matrix construction ###################
##All_ExpData is the full gene-sample expression value matrix
Re_Markers_Exp <- All_ExpData[rownames(All_ExpData) %in% Re_Marker_top30,]
dim(Re_Markers_Exp)  #8，11653
Re_Markers_Exp_matrix <- as.matrix(t(Re_Markers_Exp))
res <- cor(Re_Markers_Exp_matrix)
res2 <- signif(res, digits = 3)
library(Hmisc)# Load package



Hmisc::rcorr(as.matrix(res2), type = "pearson") -> corrlist
# Note: Need to install Hmisc package first, otherwise it will throw an error, same below. Because I'm lazy, I won't write automated code. If many people ask me to write, maybe I will. Put the matrix or dataframe that needs correlation calculation inside as.matrix()

library(tidyr)
# Correlation coefficient matrix
corrlist$r %>%  # Extract r matrix
  as_tibble() %>%  # Set to tibble format
  mutate(v = colnames(.)) %>%  # Add column names vector as a new vector to the column
  select(v, everything()) %>%  # Select all tables
  pivot_longer(2:9) -> corrdf # Convert data from wide to long format

# p-value matrix
corrlist$P %>% 
  as_tibble() %>% 
  mutate(v = colnames(.)) %>% 
  select(v, everything()) %>% 
  pivot_longer(2:9) %>%  # 2:9 means I have 8 variables
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

library(grDevices)
#windowsFonts("Arial" = windowsFont("Arial")) # Set font to prevent error in the following code
corrdf %>% 
  mutate(v = forcats::fct_reorder(v, corr), # Reorder v and name
         name = forcats::fct_reorder(name, corr)) ->corrdf2
write.csv(corrdf2,"corrdf2_Re.test.csv")
corrdf2 <- read.csv("corrdf2_Re.test.csv",row.names=1)
#library(dplyr)  Pipe %>% package

# Red and blue ↓
corrdf2 %>%  
  ggplot(aes(x = v, y = name)) + 
  geom_tile(aes(fill = corr)) + # Header fill mapping based on corr
  geom_text(aes(label = paste(corr, label, sep = "")),size=8) + ggthemes::scale_fill_gradient2_tableau("Red-Blue Diverging")+ theme(axis.text=element_text(size=14),axis.title.x=element_blank(),axis.title.y=element_blank())
 
### Extract genes with co-expression to construct interaction network (Re_Cor_relation.xlsx) ###
 
##### Construct network in Cytoscape  ####
 
#### Here explain how these proteins affect immune response ###
