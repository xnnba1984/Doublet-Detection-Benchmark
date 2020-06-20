library(reticulate)
use_python('/home/nxi/anaconda3/bin', required = T)
py_config()
library(Matrix)
library(Seurat)
library(dplyr)
library(PRROC)
setwd('/media/nxi/nxi/doublet')
source('utility.R')
getwd()

# read python module
scr <- import('scrublet')
scipy.io <- import('scipy.io')
np <- import('numpy')
os <- import('os')

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
    
    result <- scr$Scrublet(counts_matrix = t(sim.doublets), expected_doublet_rate = 0.2, random_state = 10L)
    results <- result$scrub_doublets(min_counts=2, 
                                     min_cells=3, 
                                     min_gene_variability_pctl=85, 
                                     n_prin_comps=30L)
    # prediction result
    score <- as.vector(results[[1]])
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
saveRDS(list(prauc.list, rocauc.list, score.list), 'paper_result/scrublet_hetero.rds')
pr <- unlist(lapply(prauc.list, mean)); pr
plot(pr); lines(pr)
roc <- unlist(lapply(rocauc.list, mean)); roc
plot(roc); lines(roc)

















