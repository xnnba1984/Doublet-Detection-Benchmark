library(devtools)
library(scDesign)
library(Matrix)
library(Seurat)
library(simulator)
library(pbapply)
setwd('/media/nxi/nxi/doublet')
getwd()

system.time(sim.data <- readRDS("paper_sim/sim_type.rds"))
sim.doublet <- sim.data[["8"]][[1]][[4]]; dim(sim.doublet)
sim.label <- sim.data[["8"]][[2]][[4]]; table(sim.label)
sim.label <- ifelse(sim.label=='doublet', 1, 0); table(sim.label)
sim.type <- sim.data[["8"]][[2]][[4]]; table(sim.type)
sim.type <- ifelse(sim.type=='doublet', length(table(sim.type))-1, sim.type); table(sim.type)

###########################################################
# clusters identified on different threshold 
###########################################################
score.list <- readRDS("paper_result/hybrid_type.rds")
score <- score.list[[3]][[7]][[4]]
pred.num <- seq(from = 0, by = 0.005, to = 0.25) * 4000; pred.num
set.seed(10)
cluster <- pbsapply(pred.num, function(x){
  threshhold <- sort(score, decreasing = TRUE)[x]
  pred <- as.numeric(score > threshhold)
  pred.index <- which(pred == 1)
  if(length(pred.index) > 0){
    doublet.clean <- sim.doublet[,-pred.index]
  }else{
    doublet.clean <- sim.doublet
  }
  # cluster cleaned data
  clean.seurat <- CreateSeuratObject(counts = doublet.clean, project = "doublet", min.cells = 1, min.features = 1); clean.seurat
  clean.seurat <- NormalizeData(clean.seurat)
  clean.seurat <- ScaleData(clean.seurat)
  clean.seurat <- FindVariableFeatures(clean.seurat, selection.method = "vst", nfeatures = 2000)
  clean.seurat <- RunPCA(clean.seurat)
  clean.seurat <- FindNeighbors(clean.seurat, dims = 1:10)
  clean.seurat <- FindClusters(clean.seurat, resolution = 0.5)
  return(nlevels(clean.seurat@meta.data[["seurat_clusters"]]))
}); cluster

result <- as.data.frame(cbind(pred.num, cluster)); head(result)
library(ggplot2)
ggplot(result, aes(x=pred.num, y=cluster)) + geom_line() + xlab('Doublets Prediction Number') + ylab('Cluster Number') + 
  ggtitle('Cluster Number over Doublets Prediction Number - bcds')
saveRDS(result, "paper_result/hybrid_cluster_8.rds")

#############################################################
# cluster quality
#############################################################
result <- readRDS("paper_result/hybrid_cluster_8.rds")
cluster <- 8
correct <- result[result$cluster==cluster,]; correct$pred.num; length(correct$pred.num)
set.seed(10)

error <- pbsapply(correct$pred.num, function(x){
  pred.num <- x
  threshhold <- sort(score, decreasing = TRUE)[pred.num]
  pred <- as.numeric(score > threshhold)
  pred.index <- which(pred == 1)
  doublet.clean <- sim.doublet[,-pred.index]; dim(doublet.clean)
  type.clean <- sim.type[-pred.index]; length(type.clean); table(type.clean); table(sim.type)
  
  # cluster cleaned data
  clean.seurat <- CreateSeuratObject(counts = doublet.clean, project = "doublet", min.cells = 1, min.features = 1); clean.seurat
  clean.seurat <- NormalizeData(clean.seurat)
  clean.seurat <- ScaleData(clean.seurat)
  clean.seurat <- FindVariableFeatures(clean.seurat, selection.method = "vst", nfeatures = 2000)
  clean.seurat <- RunPCA(clean.seurat)
  clean.seurat <- FindNeighbors(clean.seurat, dims = 1:10)
  clean.seurat <- FindClusters(clean.seurat, resolution = 0.5)
  cluster.label <- Idents(clean.seurat); table(cluster.label)
  clean.seurat <- RunTSNE(clean.seurat, dims = 1:10)
  
  tsne <- as.data.frame(clean.seurat@reductions[["tsne"]]@cell.embeddings)
  tsne$type <- as.numeric(type.clean)
  tsne$cluster <- cluster.label
  error <- pbsapply(0:(cluster-1), function(x){
    sub <- tsne[tsne$cluster==x,]
    return(mean(sub$type == cluster))
  })
  return(error)
}); error
colnames(error) <- correct$pred.num; error
saveRDS(error, "paper_result/hybrid_purity_8.rds")
hist(as.vector(error), breaks = 50)
mean(as.vector(error)); sd(as.vector(error))



























