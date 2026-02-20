rm(list = ls())

RandomWalk2igraph<-function(igraphM,VertexWeight,EdgeWeight=TRUE,gamma=0.7){
  adjM<-igraphM
  res<-rw(adjM,VertexWeight,gamma)
  return(drop(res))
}

rw<-function(W,p0,gamma) {
  p0<-t(p0)
  p0 <- p0/sum(p0)
  PT <- p0
  k <- 0
  delta <- 1
  Ng <- dim(W)[2]
  for (i in 1:Ng) {
    sumr<-sum(W[i,])
    if(sumr==0){
      W[i,] <-numeric(length=length(W[i,]))
    }
    if(sumr>0){
      W[i,] <- W[i,]/sumr
    }
  }
  W<-as.matrix(W)
  W <- t(W)
  while(delta>1e-10) {
    PT1 <- (1-gamma)*W
    PT2 <- PT1 %*% t(PT)
    PT3 <- (gamma*p0)
    PT4 <- t(PT2) + PT3
    delta <- sum(abs(PT4 - PT))
    PT <- PT4
    k <- k + 1
  }
  PT<-t(PT)
  rownames(PT)<-NULL
  return(PT)
}

Network2AdjacentMatrix<- function(network_PPI){
  geneid1<-as.character(network_PPI[,1])
  geneid2<-as.character(network_PPI[,2])
  Allgeneid<-unique(c(geneid1,geneid2))
  matrixPPI<-matrix(0,length(Allgeneid),length(Allgeneid))
  colnames(matrixPPI)<-Allgeneid;rownames(matrixPPI)<-Allgeneid
  for(i in 1:dim(network_PPI)[1]){
    n1<-which(Allgeneid==network_PPI[i,1])
    n2<-which(Allgeneid==network_PPI[i,2])
    matrixPPI[n1,n2]<-1
    matrixPPI[n2,n1]<-1
  }
  return(matrixPPI)
}

PPI <- read.table("NoR_PPI.txt",sep = "\t",header = T)
PPI_matrix <- Network2AdjacentMatrix(PPI)
Allgene <- colnames(PPI_matrix)
Allgene

NoR_Marker <- c("BCL7A","RUFY4","AC006978.6","SLC32A1","RP1-313I6.12","FOXP3","RAB11FIP1","IL2RA","LYAR","ACOT7","KIF20B","ARHGAP11A","EPAS1","TRGV5","RGPD2","IL7R","TCF7","FTH1P2","FTH1P23","BEST1","FTH1P20","LMNA","SATB1","RPS4XP13","EEF1A1P16","MGAT4A")

seedgene1<-NoR_Marker

seedgene_in_PPI1 <- intersect(seedgene1,Allgene)

Vector1 <- rep(0,length=length(Allgene))
Vector1[match(seedgene_in_PPI1,Allgene)] <- 1

visProbs1 <- RandomWalk2igraph(PPI_matrix,Vector1,EdgeWeight=TRUE,gamma=0.7)

result1 <- as.data.frame(cbind(Allgene,visProbs1))
colnames(result1) <- c("geneId","globalScore")

result11 <- result1[order(result1$globalScore,decreasing=T),]

write.table(result1,"random_score1.txt",sep = "\t",row.names = F,col.names = T)
write.table(result11,"random_score_sort1.txt",sep = "\t",row.names = F,col.names = T)

spyScore <- function(randomwalkScore){
    geneSet <- NoR_Marker
    index <- randomwalkScore$"geneId" %in% geneSet
    geneSetScore <- mean(as.numeric(randomwalkScore[index, "globalScore"]))
    return(geneSetScore)
}
subpathwayScore <-spyScore(result11)

randomspyScore <- function(randomwalkScore)
{
  index <- sample(nrow(randomwalkScore))
  randomwalkScore <- data.frame(geneId = randomwalkScore$geneId, globalScore = randomwalkScore[index,2])
  spy_Score <- spyScore(randomwalkScore = randomwalkScore)
  return(spy_Score)
}

spyScore_random1000 <- c()
for(i in 1:1000)
{
  ransubpathscore <- randomspyScore(result11)
  spyScore_random1000[i] <- ransubpathscore
}
