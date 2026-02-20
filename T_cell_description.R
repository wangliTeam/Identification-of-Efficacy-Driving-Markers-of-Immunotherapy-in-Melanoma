# Since in the comparison between R and NR groups, T cells showed a specific phenomenon in the NR group. This paper considers this cell to be a specific cell subpopulation in the non-responder patient group.
# We perform functional enrichment analysis on this cell
NoR_All_1_Data.markers <- read.csv("_1分好_NoR_All_Data_markers.csv",row.names=1)
NoR_All_1_top10genes <- NoR_All_1_Data.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
T_top9genes <- NoR_All_1_top10genes[NoR_All_1_top10genes$cluster=="T cells",]$gene[1:9]
# Still using those two functions GO_Enrichment_Analy, KEGG_Enrichment_Analy

# Select some genes
select_genes <- T_top9genes
# vlnplot display
p1 <- VlnPlot(New_All_Data_2, features = select_genes, pt.size=0, ncol=3)
ggsave("460_selectTcells_genes_9_VlnPlot.png", p1, width=15 ,height=10)
 
cells <- c("pDCs","Regulatory T cells","T cells")
for (j in 1:length(cells)){
	difgene<- NoR_All_1_top20genes[NoR_All_1_top20genes$cluster==cells[j],]$gene
	one_gene_list <- difgene
	GO_Enrichment_Analy(one_gene_list, cells[j],294+j,"NRP")
	KEGG_Enrichment_Analy(one_gene_list,cells[j],290+j,"NRP")
}

Re_All_1_Data.markers <- read.csv("Re_All_Data_markers.csv",row.names=1)
Re_All_1_top20genes <- Re_All_1_Data.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)

cells <- c("B cells","CD4 T cells")
for (j in 1:length(cells)){
	difgene<- Re_All_1_top20genes[Re_All_1_top20genes$cluster==cells[j],]$gene
	one_gene_list <- difgene
	GO_Enrichment_Analy(one_gene_list, cells[j],299+j,"NRP")
	KEGG_Enrichment_Analy(one_gene_list,cells[j],232+j,"NRP")
}

# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("GOplot")
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)

# What is the gene order (corresponding to cell types)
GO_Enrichment_Analy <- function(Gene_List,cellName,iNum,dataType){
	go_id_trance <- bitr(Gene_List,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db",drop = T)
	Gene_List_entrezid <- go_id_trance$ ENTREZID
	# 1. GO enrichment
	## CC represents cellular component, MF represents molecular function, BP represents biological process, ALL represents enrichment of all three processes at the same time, choose what you need, I usually do BP, MF, CC these 3 groups and then merge into a dataframe, which is convenient for extracting some pathways for plotting later.
	ego_ALL <- enrichGO(gene = Gene_List_entrezid,# We defined above
					   OrgDb=org.Hs.eg.db,
					   keyType = "ENTREZID",
					   ont = "ALL",# Enriched GO type
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
	# These three numbers represent the number of selected BP, CC, MF pathways, just set this yourself (before controlling the number of GO terms, save)
	ego_result_ALL <- as.data.frame(ego_ALL)
	ego_result_BP <- as.data.frame(ego_BP)
	ego_result_CC <- as.data.frame(ego_CC)
	ego_result_MF <- as.data.frame(ego_MF)
filename <- paste(dataType, cellName,"GO_Analyse.csv",sep='_')
write.csv(ego_result_ALL,file=filename)
####### Draw chord diagram ######################
if(dim(ego_ALL)[1]!=0){
go_enrich_df_3 <- as.data.frame(ego_ALL)
if(dim(ego_ALL)[1]>20){
	go_enrich_df_3 <- as.data.frame(ego_ALL)[1:20,]
}
	gene <- str_replace_all (go_enrich_df_3$geneID,"/" , "," )
	GO <- data.frame(ID=go_enrich_df_3$ID,Term=go_enrich_df_3$Description,Genes=gene,adj_pval=go_enrich_df_3$p.adjust,Category=go_enrich_df_3$ONTOLOGY)

	## Construct gene matrix and randomly generate logFC
	genedata=data.frame(ID=go_id_trance$SYMBOL,logFC=rnorm(length(go_id_trance$SYMBOL),mean=0,sd=2))

	circ <- circle_dat(GO, genedata)
	###### Plot #################
	chord <- chord_dat(data=circ, genes = genedata)
	# Generate matrix with selected gene list
	chord <- chord_dat(data=circ, process = GO$Term) # Generate matrix with selected GO term list
	chord <- chord_dat (data=circ, genes=genedata, process = GO$Term) # Construct data
	p=GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5,ribbon.col=Mycol_4[1:length(GO$Term)])
	i0=iNum
	imageName=paste(i0,cellName,'Chord_image.png',sep = "_")
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
# After controlling the number of displayed pathways
	if(dim(ego_result_BP)[1]>5){
		ego_result_BP =ego_result_BP [1:5,]}
if(dim(ego_result_CC)[1]>5){
	ego_result_CC =ego_result_CC [1:5,]}
if(dim(ego_result_MF)[1]>5){
	ego_result_MF =ego_result_MF [1:5,]}

	## Recombine the partial pathways we extracted above into a dataframe
	go_enrich_df_2 <- data.frame(
	ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
	Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
	GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
	type=factor(c(rep("biological process", dim(ego_result_BP)[1]), 
				  rep("cellular component", dim(ego_result_CC)[1]),
				  rep("molecular function", dim(ego_result_MF)[1])), 
				  levels=c("biological process", "cellular component","molecular function" )))
	
	## The pathway name is too long, I selected the first five words of the pathway as the pathway name
	if(nrow(go_enrich_df_2)!=0){
		for(i in 1:nrow(go_enrich_df_2)){
		  description_splite=strsplit(go_enrich_df_2$Description[i],split = " ")
		  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") # Here 5 means 5 words, can be changed by yourself
		  go_enrich_df_2$Description[i]=description_collapse
		 }
	}
	## Start drawing GO bar chart
	  ### Horizontal bar chart
	go_enrich_df_2$type_order=factor(rev(as.integer(rownames(go_enrich_df_2))),labels=rev(go_enrich_df_2$Description))# This step is necessary, to make the bars display in order, not messy

	i0=iNum+1
	p <- ggplot(data= go_enrich_df_2, aes(x=type_order,y=GeneNumber, fill=type)) + # Horizontal and vertical axis values
	  geom_bar(stat="identity", width=0.8) + # Bar chart width, can be set by yourself
	  scale_fill_manual(values = Mycol_3) + ### Color
	  coord_flip() + ## This step makes the bar chart horizontal, without this the bar chart is vertical
	  xlab("GO term") + 
	  ylab("Gene_Number") + 
	  labs(title = "The Most Enriched GO Terms")+
	  theme_bw()
	imageName<-paste(i0,cellName,'GO_Enrichment.png',sep = "_") 
	ggsave(file=imageName, plot = p, width = 12, height = 10)

}

KEGG_Enrichment_Analy<- function(Gene_List,cellName,iNum,dataType){
	go_id_trance <- bitr(Gene_List,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db",drop = T)
	Gene_List_entrezid <- go_id_trance$ ENTREZID

	# 1. KEGG enrichment
	kk <- enrichKEGG(gene =Gene_List_entrezid,keyType = "kegg",organism= "human", qvalueCutoff = 0.05, pvalueCutoff=0.05)
	
	# 2. Visualization
	  ### Bar chart
	hh <- as.data.frame(kk)# Remember to save the results yourself!
	if(dim(hh)[1]!=0){
		filename <- paste(dataType,cellName,"KEGG_Analyse.csv",sep='_')
		write.csv(hh,file=filename)
		rownames(hh) <- 1:nrow(hh)
		hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
		p=ggplot(hh,aes(y=order,x=Count,fill=p.adjust))+
		  geom_bar(stat = "identity", width=0.5)+#### Bar width
		  #coord_flip()+## Swap horizontal and vertical axes
		  scale_fill_gradient(low = '#FFB90F',high ="#D95F02")+# Color can be changed by yourself
		  labs(title = "KEGG Pathways Enrichment",
			   x = "Gene numbers", 
			   y = "Pathways")+
		  theme(axis.title.x = element_text(face = "bold",size = 16),
				axis.title.y = element_text(face = "bold",size = 16),
				legend.title = element_text(face = "bold",size = 16))+
		  theme_bw()
		i0=iNum
		imageName<-paste(i0,cellName,'cellType_KEGG_Enrichment.png',sep = "_")
		ggsave(file=imageName, plot = p, width = 8, height = 6)
		 ### Bubble chart
		rownames(hh) <- 1:nrow(hh)
		hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
		p=ggplot(hh,aes(y=order,x=Count))+
		geom_point(aes(size=Count,color=-1*p.adjust))+# Modify point size
		scale_color_gradient(low="green",high = "red")+
		labs(color=expression(p.adjust,size="Count"), 
			 x="Gene Number",y="Pathways",title="KEGG Pathway Enrichment")+
		theme_bw()
		i0=iNum+1
		imageName<-paste(i0,cellName,'KEGG_Enrichment.png',sep = "_")
		ggsave(file=imageName, plot = p, width = 8, height = 6)
	}
}
