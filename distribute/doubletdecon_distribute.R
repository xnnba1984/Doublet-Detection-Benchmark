library(devtools)
library(Matrix)
library(Seurat)
library(parallel)
library(DoubletDecon)
library(dplyr)
library(pbapply)
setwd('/media/nxi/nxi/doublet')
source('utility.R')
getwd()

# read data; parameter setting
system.time(sim.data <- readRDS("real_data/pbmc-ch.rds"))
precision.list <- list()
recall.list <- list()
tnr.list <- list()

#####################################
# batch=2
#####################################
data <- sim.data[[1]]; dim(data)
labels <- sim.data[[2]]; table(labels)
set.seed(10)
index <- sample(1:dim(data)[2], round(dim(data)[2]/2))
data1 <- data[,index]; dim(data1)
data2 <- data[,-index]; dim(data2)
label1 <- labels[index]; table(label1)
label2 <- labels[-index]; table(label2)
label <- c(label1, label2)
label <- ifelse(label=='doublet', 1, 0); table(label)

result2.list <- pblapply(list(data1, data2), function(x){
  doublet.seurat <- CreateSeuratObject(counts = x, project = "doublet",
                                       min.cells = 1 , min.features = 1); doublet.seurat
  doublet.seurat <- NormalizeData(doublet.seurat)
  doublet.seurat <- ScaleData(doublet.seurat)
  doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
  doublet.seurat <- RunPCA(doublet.seurat)
  doublet.seurat <- FindNeighbors(doublet.seurat, dims = 1:10)
  doublet.seurat <- FindClusters(doublet.seurat, resolution = 0.5)
  newFiles <- Improved_Seurat_Pre_Process(doublet.seurat, num_genes=50, write_files=FALSE)
  result <- Main_Doublet_Decon(rawDataFile=newFiles$newExpressionFile, 
                                groupsFile=newFiles$newGroupsFile, 
                                filename='DoubletDecon', 
                                location='/media/nxi/nxi/doublet/paper_result',
                                fullDataFile=NULL, 
                                removeCC=FALSE, 
                                species="hsa", 
                                rhop=1.1, 
                                write=TRUE, 
                                PMF=TRUE, 
                                useFull=FALSE, 
                                heatmap=FALSE,
                                centroids=TRUE,
                                num_doubs=100, 
                                only50=FALSE,
                                min_uniq=4,
                                nCores=-1)
  return(result)
})

# save prediction
preds <- c()
for(i in 1:length(result2.list)){
  results <- result2.list[[i]]
  pred <- results$DRS_doublet_table$isADoublet; table(pred)
  pred <- ifelse(pred==T,1,0); table(pred)
  preds <- c(preds, pred)
}; table(preds)

# result
index.doublet <- which(label==1)
index.singlet <- which(label==0)
tp <- sum(preds[which(label==1)]==1); tp
fp <- sum(preds[which(label==0)]==1); fp
fn <- sum(preds[which(label==1)]==0); fn
tn <- sum(preds[which(label==0)]==0); tn

precision <- tp/(tp + fp); precision
recall <- tp/(tp + fn); recall
tnr <- tn/(tn + fp); tnr
f1 <- 2 * precision * recall / (precision + recall); f1

#####################################
# batch=4
#####################################
data <- sim.data[[1]]; dim(data)
labels <- sim.data[[2]]; table(labels)
set.seed(10)
index <- sample(1:dim(data)[2], dim(data)[2])
chunk <- split(index, cut(seq_along(index), 4, labels = FALSE)) 
data.chunk <- pblapply(chunk, function(x){
  return(data[,x])
})
labels4 <- labels[unlist(chunk,use.names = F)]

result4.list <- pblapply(data.chunk, function(x){
  doublet.seurat <- CreateSeuratObject(counts = x, project = "doublet",
                                       min.cells = 1 , min.features = 1); doublet.seurat
  doublet.seurat <- NormalizeData(doublet.seurat)
  doublet.seurat <- ScaleData(doublet.seurat)
  doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
  doublet.seurat <- RunPCA(doublet.seurat)
  doublet.seurat <- FindNeighbors(doublet.seurat, dims = 1:10)
  doublet.seurat <- FindClusters(doublet.seurat, resolution = 0.5)
  newFiles <- Improved_Seurat_Pre_Process(doublet.seurat, num_genes=50, write_files=FALSE)
  result <- Main_Doublet_Decon(rawDataFile=newFiles$newExpressionFile, 
                               groupsFile=newFiles$newGroupsFile, 
                               filename='DoubletDecon', 
                               location='/media/nxi/nxi/doublet/paper_result',
                               fullDataFile=NULL, 
                               removeCC=FALSE, 
                               species="hsa", 
                               rhop=1.1, 
                               write=TRUE, 
                               PMF=TRUE, 
                               useFull=FALSE, 
                               heatmap=FALSE,
                               centroids=TRUE,
                               num_doubs=100, 
                               only50=FALSE,
                               min_uniq=4,
                               nCores=-1)
  return(result)
})
