library(Matrix)
library(Seurat)
library(dplyr)
library(PRROC)
library(DoubletFinder)
source('utility.R')
set.seed(2020)

###############################################################################################
# calculate auprc and auroc on simulation datasets under different doublet rates, cell types, 
# sequencing depth, and heterogeneity between cell types
###############################################################################################
# read simulation data under different doublet rates, cell types, 
# sequencing depth, and heterogeneity between cell types
# comment out unused data to choose one out of four
sim.data <- readRDS("paper_sim/sim_rate.rds")
sim.data <- readRDS("paper_sim/sim_type.rds")
sim.data <- readRDS("paper_sim/sim_depth.rds")
sim.data <- readRDS("paper_sim/sim_hetero.rds")
prauc.list <- list()
rocauc.list <- list()
score.list <- list()
# loop over each simulation dataset
system.time({
  for(i in 1:length(sim.data)){
    print(i)
    print('====================')
    counts <- sim.data[[i]][[1]]
    labels <- sim.data[[i]][[2]]
    prauc <- c()
    rocauc <- c()
    scores <- list()
    # each experimental condition incldues five simulation datasets
    # calculate on all of them, we use the first one in the paper
    for(j in 1:length(counts)){
      print(j)
      # read data and labels
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
# save the result, change the file name based on which experimental condition we use
saveRDS(list(prauc.list, rocauc.list, score.list), 'paper_result/doubletfinder_hetero.rds')

##############################################################################################################################
# Impact of doublet detection on the identification of highly variable genes (HVGs)
##############################################################################################################################
# simulation data with different doublet rates
sim.data <- readRDS("paper_sim/sim_rate.rds")
# we use three doublet rates: 0.1, 0.2, 0.4
rate <- 0.4
data <- sim.data[[as.character(rate)]]
# use the first simulation dataset out of five
data.matrix <- data[[1]][[1]]; dim(data.matrix)
data.type <- data[[2]][[1]]; table(data.type)
data.type <- ifelse(data.type=='doublet', 1, 0); table(data.type)
index.doublet <- which(data.type==1)
data.clean <- data.matrix[,-index.doublet]; dim(data.clean)
# read calculated doublet scores, the third element in the list
score.list <- readRDS('paper_result/doubletfinder_rate.rds')[[3]]
# use the score of the first simulation dataset
score <- score.list[[20]][[1]]
# choose threshold 
threshold <- sort(score, decreasing = T)[floor(length(score)*rate)]; threshold
index.pred <- which(score > threshold)
data.pred <- data.matrix[,-index.pred]; dim(data.pred)

# clean variable genes
doublet.seurat <- CreateSeuratObject(counts = data.clean, project = "doublet", min.cells = 0); doublet.seurat
doublet.seurat <- NormalizeData(doublet.seurat)
doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
clean.variable <- VariableFeatures(doublet.seurat)

# contiminated variable genes
doublet.seurat <- CreateSeuratObject(counts = data.matrix, project = "doublet", min.cells = 0); doublet.seurat
doublet.seurat <- NormalizeData(doublet.seurat)
doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
doublet.variable <- VariableFeatures(doublet.seurat)

# post-doublet-detection variable genes
doublet.seurat <- CreateSeuratObject(counts = data.pred, project = "doublet", min.cells = 0); doublet.seurat
doublet.seurat <- NormalizeData(doublet.seurat)
doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
pred.variable <- VariableFeatures(doublet.seurat)

# Jaccard index: negative control
J.clean.doublet <- 
  length(intersect(clean.variable, doublet.variable))/length(union(clean.variable, doublet.variable)); J.clean.doublet
# # Jaccard index: post-doublet-detection
J.clean.pred <- 
  length(intersect(clean.variable, pred.variable))/length(union(clean.variable, pred.variable)); J.clean.pred



