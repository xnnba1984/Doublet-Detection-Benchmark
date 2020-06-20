library(Matrix)
library(Seurat)
library(scds)
library(SingleCellExperiment)
library(PRROC)
setwd('/media/nxi/nxi/doublet')
source('utility.R')
getwd()

# read data; parameter setting
sim.data <- readRDS("paper_sim/sim_hetero.rds")
prauc.list <- list()
rocauc.list <- list()
score.list <- list()
set.seed(10)
system.time({
  for(i in 1:length(sim.data)){
    print(i)
    print('====================')
    counts <- sim.data[[i]][[1]]
    labels <- sim.data[[i]][[2]]
    prauc <- c()
    rocauc <- c()
    scores <- list()
    for(j in 1:length(counts)){
      print(j)
      sim.doublets <- counts[[j]]; dim(sim.doublets)
      sim.labels <- labels[[j]]; table(sim.labels)
      sim.labels <- ifelse(sim.labels=='doublet', 1, 0); table(sim.labels)
      sim.doublets <- sim.doublets[which(rowSums(sim.doublets) != 0),]; dim(sim.doublets)

      sce <- SingleCellExperiment(assays = list(counts = sim.doublets))
      # co-expression, boost, hybrid
      system.time({
        sce <- cxds_bcds_hybrid(sce)
        #sce <- cxds(sce)
        #sce <- bcds(sce)
      })
      CD <- colData(sce)
      #score <- CD$cxds_score
      #score <- CD$bcds_score
      score <- CD$hybrid_score
      fg <- score[sim.labels == 1]
      bg <- score[sim.labels == 0]
      pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
      roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
      prauc[j] <- pr$auc.integral
      rocauc[j] <- roc$auc
      scores[[j]] <- score
    }
    prauc.list[[i]] <- prauc
    rocauc.list[[i]] <- rocauc
    score.list[[i]] <- scores
  }
})
saveRDS(list(prauc.list, rocauc.list, score.list), 'paper_result/hybrid_hetero.rds')
pr <- unlist(lapply(prauc.list, mean)); pr
plot(pr); lines(pr)
roc <- unlist(lapply(rocauc.list, mean)); roc
plot(roc); lines(roc)


















