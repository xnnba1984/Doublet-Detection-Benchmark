library(Matrix)
library(Seurat)
library(dplyr)
library(PRROC)
library(DoubletFinder)
setwd('/media/nxi/nxi/doublet')
source('utility.R')
getwd()

# read data; parameter setting
sim.data <- readRDS("paper_sim/sim_hetero.rds")
prauc.list <- list()
rocauc.list <- list()
score.list <- list()
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
      
      ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
      doublet.seurat <- CreateSeuratObject(counts = sim.doublets, project = "doublet", min.cells = 1); doublet.seurat
      doublet.seurat <- NormalizeData(doublet.seurat)
      doublet.seurat <- ScaleData(doublet.seurat)
      doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
      doublet.seurat <- RunPCA(doublet.seurat)
      
      ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
      sweep.res.doublet <- paramSweep_v3(doublet.seurat, PCs = 1:10, sct = FALSE, num.cores = 44)
      sweep.stats.doublet <- summarizeSweep(sweep.res.doublet, GT = FALSE)
      bcmvn.doublet <- find.pK(sweep.stats.doublet)
      pK <- bcmvn.doublet$pK[which.max(bcmvn.doublet$BCmetric)]; pK <- as.numeric(levels(pK))[pK]; pK
      doublet.seurat <- doubletFinder_v3(doublet.seurat, PCs = 1:10, pN = 0.25, pK = pK, nExp = sum(sim.labels==1))
      attribute <- paste('pANN', 0.25, pK, sum(sim.labels==1), sep = '_'); attribute
      score <- doublet.seurat@meta.data[[attribute]]
      
      # prediction result
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
saveRDS(list(prauc.list, rocauc.list, score.list), 'paper_result/doubletfinder_hetero.rds')
pr <- unlist(lapply(prauc.list, mean)); pr
plot(pr); lines(pr)
roc <- unlist(lapply(rocauc.list, mean)); roc
plot(roc); lines(roc)





















