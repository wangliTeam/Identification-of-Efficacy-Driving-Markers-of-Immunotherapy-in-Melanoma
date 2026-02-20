library(stringr)

exprSet_Unamand  <- New_All_Data_2@assays$RNA@data
metarSet_Unamand <- New_All_Data_2@meta.data

meta_before <- read.csv("Survival_Condation.csv",header=T,row.names=1)
V5_before <- strsplit(metarSet_Unamand$V5,split="_")
for (i in 1:length(V5_before)){
   if (length(V5_before[[i]])>2){
       V5_before[[i]] = V5_before[[i]][1:2]
   }
}
V5_split <- t(as.data.frame(V5_before))
rownames(V5_split) <- rownames(metarSet_Unamand)

metarSet_amand1 <- metarSet_Unamand[,c('V6','Tag','cellTypes')]
metarSet_amand2 <- cbind(V5_split,metarSet_amand1)
colnames(metarSet_amand2) <- c("Pre_or_Post","Patient","Re_or_Nor","Tag","cellTypes")

dim(metarSet_amand2)

Add_Col <- data.frame()
for (i in 1:dim(metarSet_amand2)[1]){
   ad_data <- subset(meta_before,Patient.ID==metarSet_amand2[i,"Patient"],c(Overall.survival..days.,Status..Alive.0..Dead.1.))
   if(length(row.names(ad_data))==0){
       ad_data[1,]<- c(NA,NA)
   }
   Add_Col = rbind( Add_Col,ad_data)
}

meta <- cbind(metarSet_amand2,Add_Col)
colnames(meta)[6:7]<- c("Overall_survival_days","Status")

k = !(is.na(meta$"Overall_survival_days")|is.na(meta$"Status"))
exprSet = exprSet_Unamand[,k]
meta = meta[k,]
save(meta,exprSet,file = "survival_model.Rdata")

rm(list = ls())
load("survival_model.Rdata")

library(survival)
library(survminer)

Re_Marker_top30 <- c("IGHG1","IGHGP","IGHG2","IGHG4","SCIMP","IGHV3-30","LINC00861","TMEM63A")
NoR_Marker_top30 <- c("BCL7A","RUFY4","AC006978.6","SLC32A1","RP1-313I6.12","FOXP3","RAB11FIP1","IL2RA","LYAR","ACOT7","KIF20B","ARHGAP11A","EPAS1","TRGV5","RGPD2","IL7R","TCF7","FTH1P2","FTH1P23","BEST1","FTH1P20","LMNA","SATB1","RPS4XP13","EEF1A1P16","MGAT4A")

gs=NoR_Marker_top30
splots <- lapply(gs, function(g){
  meta$gene=ifelse(exprSet[g,] > median(exprSet[g,]),'high','low')
  sfit1=survfit(Surv(Overall_survival_days, Status)~gene, data=meta)
  ggsurvplot(sfit1,pval =TRUE, data = meta, risk.table = TRUE)
}) 
res <- arrange_ggsurvplots(splots, print = TRUE,  
                    ncol = 2, nrow = 1,conf.int = TRUE,
                    risk.table.height = 0.25,font.title=gene)

ggsave("NoR_Survival.pdf", res,width = 10,height = 6)

NoR_Marker_top30[c(1,5,8,9,10,14,16,17,20,22,23)]

gs=Re_Marker_top30
splots <- lapply(gs, function(g){
  meta$gene=ifelse(exprSet[g,] > median(exprSet[g,]),'high','low')
  sfit1=survfit(Surv(Overall_survival_days, Status)~gene, data=meta)
  ggsurvplot(sfit1,pval =TRUE, data = meta, risk.table = TRUE)
}) 
res <- arrange_ggsurvplots(splots, print = TRUE,  
                    ncol = 2, nrow = 1,conf.int = TRUE,
                    risk.table.height = 0.25,font.title=gene)

ggsave("Re_Survival.pdf", res,width = 10,height = 6)

Re_Marker_top30[7]
