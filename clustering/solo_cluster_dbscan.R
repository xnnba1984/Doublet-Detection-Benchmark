library(reticulate)
library(devtools)
library(Matrix)
library(Seurat)
library(simulator)
library(pbapply)
library(dbscan)
np <- import("numpy")

# read simulation datasets with 4, 6, 8 cell types
# change cell type and location accordingly
system.time(sim.data <- readRDS("paper_sim/sim_type.rds"))
# 5th simulation in each cell type group
sim.doublet <- sim.data[["8"]][[1]][[5]]; dim(sim.doublet)
# read cell type and doublet labels
sim.type <- sim.data[["8"]][[2]][[5]]; table(sim.type)
sim.type <- ifelse(sim.type=='doublet', length(table(sim.type))-1, sim.type); table(sim.type)

###########################################################
# clusters identified on different threshold 
###########################################################
# read doublet score of simulation data with 4, 6, 8 cell types
score <- as.numeric(np$load('paper_sim/loom/output_type/type_7_5.loom/softmax_scores.npy')); length(score); hist(score)
# 19 simulation datasets with cell type increased from 2 to 20
# each cell type contains 5 simulations
# we use 4, 6, 8 cell types, the fifth simulation dataset
score <- score.list[[3]][[7]][[5]]
# identification rate increased from 0 to 25% with a step size 0.5%
pred.num <- seq(from = 0, by = 0.005, to = 0.25) * 4000; pred.num
# loop over each identification rate
cluster <- pbsapply(pred.num, function(x){
  # select threhold based on identification rate
  threshhold <- sort(score, decreasing = TRUE)[x]
  pred <- as.numeric(score > threshhold)
  pred.index <- which(pred == 1)
  # remove predicted doublets
  if(length(pred.index) > 0){
    doublet.clean <- sim.doublet[,-pred.index]
  }else{
    doublet.clean <- sim.doublet
  }
  # cluster cleaned data by louvain clustering
  clean.seurat <- CreateSeuratObject(counts = doublet.clean, project = "doublet", 
                                     min.cells = 1, min.features = 1); clean.seurat
  clean.seurat <- NormalizeData(clean.seurat)
  clean.seurat <- ScaleData(clean.seurat)
  clean.seurat <- FindVariableFeatures(clean.seurat, selection.method = "vst", nfeatures = 2000)
  clean.seurat <- RunPCA(clean.seurat, npcs = 30)
  clean.seurat <- RunTSNE(clean.seurat, perplexity=30)
  x <- clean.seurat@reductions[["tsne"]]@cell.embeddings; dim(x); class(x)
  res <- dbscan(x, eps = 3, minPts = 5)
  print(table(res$cluster))
  cluster <- names(table(res$cluster)); cluster
  return(ifelse('0'%in%cluster, length(cluster)-1, length(cluster)))
}); cluster
# save identification rate and corresponding clustering result
# change the name based on cell type and location
result <- as.data.frame(cbind(pred.num, cluster)); head(result)
saveRDS(result, "paper_result/solo_cluster_dbscan_8.rds")

##########################################################################################################################
# proportion of singlets in correctly identified clusters (cluster quality)
##########################################################################################################################
# 4, 6, 8 clusters
result <- readRDS("paper_result/solo_cluster_dbscan_8.rds")
cluster <- 8
# select identification rates that have correct clustering result
correct <- result[result$cluster==cluster,]; correct$pred.num; length(correct$pred.num)
# loop over identification rates
error <- pbsapply(correct$pred.num, function(x){
  pred.num <- x
  # select threhold based on identification rate
  threshhold <- sort(score, decreasing = TRUE)[pred.num]
  pred <- as.numeric(score > threshhold)
  pred.index <- which(pred == 1)
  # remove predicted doublets
  doublet.clean <- sim.doublet[,-pred.index]; dim(doublet.clean)
  type.clean <- sim.type[-pred.index]; length(type.clean); table(type.clean); table(sim.type)
  
  # cluster cleaned data by louvain clustering
  clean.seurat <- CreateSeuratObject(counts = doublet.clean, project = "doublet", 
                                     min.cells = 1, min.features = 1); clean.seurat
  clean.seurat <- NormalizeData(clean.seurat)
  clean.seurat <- ScaleData(clean.seurat)
  clean.seurat <- FindVariableFeatures(clean.seurat, selection.method = "vst", nfeatures = 2000)
  clean.seurat <- RunPCA(clean.seurat, npcs = 30)
  clean.seurat <- RunTSNE(clean.seurat, perplexity=30)
  x <- clean.seurat@reductions[["tsne"]]@cell.embeddings; dim(x); class(x)
  res <- dbscan(x, eps = 3, minPts = 5)
  cluster.label <- res$cluster; table(cluster.label)
  noise <- which(cluster.label==0)
  if(length(noise) != 0){
    x <- x[-noise,]
    cluster.label <- cluster.label[-noise]
    type.clean <- type.clean[-noise]
  }
  # low dimension (tsne) embeddings
  tsne <- as.data.frame(x)
  # assign cell type and doubelt labels
  tsne$type <- as.numeric(type.clean)
  # assign clustring labels
  tsne$cluster <- cluster.label
  # calculate doublet rate in each cluster (1-purity)
  error <- pbsapply(1:cluster, function(x){
    sub <- tsne[tsne$cluster==x,]
    return(mean(sub$type == cluster))
  })
  return(error)
}); error
# save result, change cell type accordingly
colnames(error) <- correct$pred.num; error
saveRDS(error, "paper_result/solo_purity_dbscan_8.rds")




