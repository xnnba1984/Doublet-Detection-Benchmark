library(devtools)
library(scDesign)
library(Matrix)
library(Seurat)
library(simulator)
library(pbapply)
library(ggplot2)
setwd('/media/nxi/nxi/doublet')
getwd()

system.time(sim.data <- readRDS("paper_sim/sim_type.rds"))
sim.doublet <- sim.data[["8"]][[1]][[5]]; dim(sim.doublet)
sim.label <- sim.data[["8"]][[2]][[5]]; table(sim.label)
sim.label <- ifelse(sim.label=='doublet', 1, 0); table(sim.label)
sim.type <- sim.data[["8"]][[2]][[5]]; table(sim.type)
sim.type <- ifelse(sim.type=='doublet', length(table(sim.type))-1, sim.type); table(sim.type)
set.seed(10)
system.time({
# cluster original doublets data
doublet.seurat <- CreateSeuratObject(counts = sim.doublet, project = "doublet", 
                                     min.cells = 1, min.features = 1); doublet.seurat
doublet.seurat <- NormalizeData(doublet.seurat)
doublet.seurat <- ScaleData(doublet.seurat)
doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
doublet.seurat <- RunPCA(doublet.seurat)
doublet.seurat <- FindNeighbors(doublet.seurat, dims = 1:10)
doublet.seurat <- FindClusters(doublet.seurat, resolution = 0.5)
})
cluster.label <- Idents(doublet.seurat); table(cluster.label)
doublet.seurat <- RunTSNE(doublet.seurat, dims = 1:10)

tsne <- as.data.frame(doublet.seurat@reductions[["tsne"]]@cell.embeddings)
tsne$type <- sim.label
tsne$type <- ifelse(tsne$type==0, 'Singlet', 'Doublet'); table(tsne$type)

############################################################################################################
# before clean
############################################################################################################
pdf('paper_figure/cluster_4.pdf')
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(aes(color = type, size = type)) +
  xlab('tSNE 1') + ylab('tSNE 2') + labs(title = 'Clustering with Doublets') + 
  scale_color_manual(breaks = c("Doublet", "Singlet"), values=c("red", "gray60")) + 
  scale_size_manual(values=c(0.7,0.3)) + theme(plot.title = element_text(hjust = 0.5))
dev.off()

############################################################################################################
# after clean 
############################################################################################################
score.list <- readRDS("paper_result/scrublet_type.rds")
score <- score.list[[3]][[3]][[5]]
threshhold <- sort(score, decreasing = TRUE)[400]
pred <- as.numeric(score > threshhold)
pred.index <- which(pred == 1)
tsne.clean <- tsne[-pred.index,]

pdf('paper_figure/cluster_clean_4.pdf')
ggplot(tsne.clean, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(aes(color = type, size = type)) +
  xlab('tSNE 1') + ylab('tSNE 2') + labs(title = 'Clustering after Dropping Detected Doublets') + 
  scale_color_manual(breaks = c("Doublet", "Singlet"), values=c("red", "gray60")) + 
  scale_size_manual(values=c(0.7,0.3)) + theme(plot.title = element_text(hjust = 0.5))
dev.off()

###########################################################
# clusters identified on different threshold 
###########################################################
score.list <- readRDS("paper_result/scrublet_type.rds")
score <- score.list[[3]][[7]][[5]]
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
  ggtitle('Cluster Number over Doublets Prediction Number - Scrublet')
saveRDS(result, "paper_result/scrublet_cluster_8.rds")

#############################################################
# cluster quality
#############################################################
result <- readRDS("paper_result/scrublet_cluster_8.rds")
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
saveRDS(error, "paper_result/scrublet_purity_8.rds")
hist(as.vector(error), breaks = 50)
mean(as.vector(error)); sd(as.vector(error))

#########################################
# plot cluster result
#########################################
# tsne <- as.data.frame(clean.seurat@reductions[["tsne"]]@cell.embeddings)
# tsne$type <- as.numeric(type.clean)
# tsne$type <- sapply(tsne$type, function(x){
#   return(ifelse(x%in%c(0,1,2,3), paste('Singlet'), paste('Doublet')))
# })
# 
# ggplot(tsne, aes(x = tSNE_1, y = tSNE_2)) +
#   geom_point(aes(color = type, size = type)) +
#   xlab('Dim 1') + ylab('Dim 2') + labs(title = 'tSNE Plot of Clustering after Cleaning - Scrublet') +
#   scale_color_manual(breaks = c("Doublet", "Singlet"), values=c("red", "gray60")) + 
#   scale_size_manual(values=c(0.7,0.3))
# ggsave("figure/scrublet_clean_cluster_4.png")






















