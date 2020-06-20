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

decon_reslut <- list()

# read data; parameter setting
#loc <- 'real_data/pbmc-ch.rds'
#loc <- 'real_data/cline-ch.rds'
#loc <- 'real_data/mkidney-ch.rds'
#loc <- 'real_data/hm-12k.rds'
#loc <- 'real_data/hm-6k.rds'
#loc <- 'real_data/pbmc-1A-dm.rds'
#loc <- 'real_data/pbmc-1B-dm.rds'
#loc <- 'real_data/pbmc-1C-dm.rds'
#loc <- 'real_data/pbmc-2ctrl-dm.rds'
#loc <- 'real_data/pbmc-2stim-dm.rds'
#loc <- 'real_data/J293t-dm.rds'
#loc <- 'real_data/pdx-MULTI.rds'
#loc <- 'real_data/HMEC-orig-MULTI.rds'
#loc <- 'real_data/HMEC-rep-MULTI.rds'
#loc <- 'real_data/HEK-HMEC-MULTI.rds'
loc <- 'real_data/nuc-MULTI.rds'

data <- readRDS(loc)
count <- data[[1]]; dim(count)
label <- data[[2]]; table(label)
label <- ifelse(label == 'doublet', 1, 0); table(label)

doublet.seurat <- CreateSeuratObject(counts = count, project = "doublet",
                                     min.cells = 1 , min.features = 1); doublet.seurat
doublet.seurat <- NormalizeData(doublet.seurat)
doublet.seurat <- ScaleData(doublet.seurat)
doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
doublet.seurat <- RunPCA(doublet.seurat)
doublet.seurat <- FindNeighbors(doublet.seurat, dims = 1:10)
doublet.seurat <- FindClusters(doublet.seurat, resolution = 0.5)

system.time({
newFiles <- Improved_Seurat_Pre_Process(doublet.seurat, num_genes=50, write_files=FALSE)
results <- Main_Doublet_Decon(rawDataFile=newFiles$newExpressionFile, 
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
})
# save prediction
pred <- results$DRS_doublet_table$isADoublet; table(pred)
pred <- ifelse(pred==T,1,0); table(pred)
# result
index.doublet <- which(label==1)
index.singlet <- which(label==0)
tp <- sum(pred[which(label==1)]==1); tp
fp <- sum(pred[which(label==0)]==1); fp
fn <- sum(pred[which(label==1)]==0); fn
tn <- sum(pred[which(label==0)]==0); tn

precision <- tp/(tp + fp); precision
recall <- tp/(tp + fn); recall
tnr <- tn/(tn + fp); tnr
f1 <- 2 * precision * recall / (precision + recall); f1

decon_reslut <- append(decon_reslut, list(pred))

#####################################################################################
# running time
#####################################################################################
# read data; parameter setting
datasets <- c(
  'real_data/pbmc-ch.rds',
  'real_data/cline-ch.rds',
  'real_data/mkidney-ch.rds',
  #'real_data/hm-12k.rds',
  'real_data/hm-6k.rds',
  'real_data/pbmc-1A-dm.rds',
  'real_data/pbmc-1B-dm.rds',
  'real_data/pbmc-1C-dm.rds',
  #'real_data/pbmc-2ctrl-dm.rds',
  'real_data/pbmc-2stim-dm.rds',
  #'real_data/J293t-dm.rds',
  'real_data/pdx-MULTI.rds',
  'real_data/HMEC-orig-MULTI.rds',
  'real_data/HMEC-rep-MULTI.rds',
  'real_data/HEK-HMEC-MULTI.rds'
  #'real_data/nuc-MULTI.rds'
)
times <- pbsapply(datasets, function(loc){
  data <- readRDS(loc)
  count <- data[[1]]; dim(count)
  label <- data[[2]]; table(label)
  label <- ifelse(label == 'doublet', 1, 0); table(label)
  doublet.rate <- sum(label==1) / length(label); doublet.rate
  time <- system.time({
    doublet.seurat <- CreateSeuratObject(counts = count, project = "doublet",
                                         min.cells = 1 , min.features = 1); doublet.seurat
    doublet.seurat <- NormalizeData(doublet.seurat)
    doublet.seurat <- ScaleData(doublet.seurat)
    doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
    doublet.seurat <- RunPCA(doublet.seurat)
    doublet.seurat <- FindNeighbors(doublet.seurat, dims = 1:10)
    doublet.seurat <- FindClusters(doublet.seurat, resolution = 0.5)
    newFiles <- Improved_Seurat_Pre_Process(doublet.seurat, num_genes=50, write_files=FALSE)
    results <- Main_Doublet_Decon(rawDataFile=newFiles$newExpressionFile, 
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
  })
  return(time[3])
}, simplify = T); times




