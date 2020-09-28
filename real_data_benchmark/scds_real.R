library(Matrix)
library(Seurat)
library(scds)
library(SingleCellExperiment)
library(PRROC)
source('utility.R')
set.seed(2020)
##############################################################
# calculate doublet score on 16 benchmark datasets
##############################################################
# list to save doublet scores
score.list.cxds <- list()
score.list.bcds <- list()
score.list.hybrid <- list()
# read real data from local
# change the location accordingly
locs <- c('real_data/pbmc-ch.rds', 'real_data/cline-ch.rds', 'real_data/mkidney-ch.rds', 'real_data/hm-12k.rds', 
          'real_data/hm-6k.rds', 'real_data/pbmc-1A-dm.rds', 'real_data/pbmc-1B-dm.rds', 'real_data/pbmc-1C-dm.rds',
          'real_data/pbmc-2ctrl-dm.rds', 'real_data/pbmc-2stim-dm.rds', 'real_data/J293t-dm.rds', 'real_data/pdx-MULTI.rds',
          'real_data/HMEC-orig-MULTI.rds', 'real_data/HMEC-rep-MULTI.rds', 'real_data/HEK-HMEC-MULTI.rds', 
          'real_data/nuc-MULTI.rds')

# loop over each dataset
for(loc in locs){
  data <- readRDS(loc)
  count <- data[[1]]; dim(count)
  label <- data[[2]]; table(label)
  label <- ifelse(label == 'doublet', 1, 0); table(label)
  doublet.rate <- sum(label==1) / length(label); doublet.rate
  sce <- SingleCellExperiment(assays = list(counts = count))
  # cxds, bcds, hybrid
  system.time({
    sce <- cxds_bcds_hybrid(sce)
  })
  CD <- colData(sce)
  score.cxds <- CD$cxds_score
  score.bcds <- CD$bcds_score
  score <- CD$hybrid_score
  # save scores
  score.list.cxds <- append(score.list.cxds, list(score.cxds))
  score.list.bcds <- append(score.list.bcds, list(score.bcds))
  score.list.hybrid <- append(score.list.hybrid, list(score.hybrid))
}
# save results, change the location accordingly
saveRDS(score.list.cxds, 'paper_result/cxds_real_score.rds')
saveRDS(score.list.bcds, 'paper_result/bcds_real_score.rds')
saveRDS(score.list.hybrid, 'paper_result/hybrid_real_score.rds')

##########################################################################################################################
# pr, recall, and tnr under 10%, 20%, and 40% identification rates
##########################################################################################################################
# read score list
# change to cxds or hybrid accordingly
score.list <- readRDS('paper_result/bcds_real_score.rds')
# 16 benchmark datasets locations
locs <- c('real_data/pbmc-ch.rds', 'real_data/cline-ch.rds', 'real_data/mkidney-ch.rds', 'real_data/hm-12k.rds', 
          'real_data/hm-6k.rds', 'real_data/pbmc-1A-dm.rds', 'real_data/pbmc-1B-dm.rds', 'real_data/pbmc-1C-dm.rds',
          'real_data/pbmc-2ctrl-dm.rds', 'real_data/pbmc-2stim-dm.rds', 'real_data/J293t-dm.rds', 'real_data/pdx-MULTI.rds',
          'real_data/HMEC-orig-MULTI.rds', 'real_data/HMEC-rep-MULTI.rds', 'real_data/HEK-HMEC-MULTI.rds', 
          'real_data/nuc-MULTI.rds')
# identification rates
rs <- c(0.1, 0.2, 0.4)
# result matrix; 16 rows: each dataset per row; 9 cols: precision, recall, tnr per method
results <- matrix(, nrow = length(locs), ncol = 0)

# loop over identification rates
for(r in rs){
  print('====================')
  print(r)
  precisions <- c()
  recalls <- c()
  tnrs <- c()
  result <- matrix(data = 0, nrow = length(locs), ncol=3)
  for(i in 1:length(locs)){
    print(locs[i])
    data <- readRDS(locs[i])
    # obtain the doublet labels
    label <- data[[2]]; table(label)
    label <- ifelse(label == 'doublet', 1, 0); table(label)
    # calculate threshold based on identification rate
    score <- score.list[[i]]
    d <- floor(length(label) * r); d
    thresh <- sort(score, decreasing = T)[d]; thresh
    # predict doublet based on threshold
    pred <- score > thresh; table(pred)
    # result
    tp <- sum(pred[which(label==1)]==1); tp
    fp <- sum(pred[which(label==0)]==1); fp
    fn <- sum(pred[which(label==1)]==0); fn
    tn <- sum(pred[which(label==0)]==0); tn
    
    precision <- tp/(tp + fp); precision
    recall <- tp/(tp + fn); recall
    tnr <- tn/(tn + fp); tnr
    
    precisions[i] <- precision
    recalls[i] <- recall
    tnrs[i] <- tnr
  }
  result <- cbind(precisions, recalls, tnrs)
  colnames(result) <- paste(colnames(result), r, sep = '_')
  results <- cbind(results, result)
}
# changel the location and name accordingly
write.table(round(results,3), 'threshold.txt', row.names = F)

##########################################################################################################################
# pr, recall, and tnr under the thresholds determined by doubletdecon
##########################################################################################################################
# read doublet scores
# change to cxds or bcds accordingly
score.list <- readRDS('paper_result/hybrid_real_score.rds')
# doublet # selected by doubletdecon
d <- c(7872,2822,8417,2813,1223,1493,1961,4077,3479,18007,8448,3124)
precisions <- c()
recalls <- c()
tnrs <- c()

# loop over 16 datasets
for(i in 1:length(locs)){
  # obtain doublet labels
  data <- readRDS(locs[i])
  label <- data[[2]]; table(label)
  label <- ifelse(label == 'doublet', 1, 0); table(label)
  score <- score.list[[i]]
  # calculate threshold based on doublet number 
  thresh <- sort(score, decreasing = T)[d[i]]
  # predict doublets
  pred <- score > thresh; table(pred)
  # result
  tp <- sum(pred[which(label==1)]==1); tp
  fp <- sum(pred[which(label==0)]==1); fp
  fn <- sum(pred[which(label==1)]==0); fn
  tn <- sum(pred[which(label==0)]==0); tn
  
  precision <- tp/(tp + fp); precision
  recall <- tp/(tp + fn); recall
  tnr <- tn/(tn + fp); tnr
  
  precisions[i] <- precision
  recalls[i] <- recall
  tnrs[i] <- tnr
}
# save the result accordingly
names(precisions) <- locs; precisions
names(recalls) <- locs; recalls
names(tnrs) <- locs; tnrs

#####################################################################################
# running time on 16 benchmark datasets
#####################################################################################
# location of 16 datasets
# change the location accordingly
locs <- c('real_data/pbmc-ch.rds', 'real_data/cline-ch.rds', 'real_data/mkidney-ch.rds', 'real_data/hm-12k.rds', 
          'real_data/hm-6k.rds', 'real_data/pbmc-1A-dm.rds', 'real_data/pbmc-1B-dm.rds', 'real_data/pbmc-1C-dm.rds',
          'real_data/pbmc-2ctrl-dm.rds', 'real_data/pbmc-2stim-dm.rds', 'real_data/J293t-dm.rds', 'real_data/pdx-MULTI.rds',
          'real_data/HMEC-orig-MULTI.rds', 'real_data/HMEC-rep-MULTI.rds', 'real_data/HEK-HMEC-MULTI.rds', 
          'real_data/nuc-MULTI.rds')

# loop over each dataset
times <- pbsapply(datasets, function(loc){
  data <- readRDS(loc)
  count <- data[[1]]; dim(count)
  label <- data[[2]]; table(label)
  label <- ifelse(label == 'doublet', 1, 0); table(label)
  doublet.rate <- sum(label==1) / length(label); doublet.rate
  time <- system.time({
    sce <- SingleCellExperiment(assays = list(counts = count))
    # change to other method accordingly
    sce <- cxds_bcds_hybrid(sce)
  })
  return(time[3])
}, simplify = T); times

##################################################################################
# stability on subsamples of real datasets
##################################################################################
# read subsamples
loc <- 'real_data/pbmc-2ctrl-dm_sub.rds'
data <- readRDS(loc)
counts <- data[[1]]
labels <- data[[2]]
prauc <- c()
rocauc <- c()
scores <- list()

# loop over 20 subsamples with 90% droplets and 90% genes
system.time({
  for(i in 1:length(counts)){
    print(i)
    count <- counts[[i]]
    label <- labels[[i]]
    label <- ifelse(label == 'doublet', 1, 0); table(label)
    sce <- SingleCellExperiment(assays = list(counts = count))
    # change to other method accordingly
    sce <- cxds_bcds_hybrid(sce)
    CD <- colData(sce)
    score <- CD$hybrid_score
    scores <- append(scores, list(score))
    fg <- score[label==1]
    bg <- score[label==0]
    pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    prauc[i] <- pr$auc.integral
    roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    rocauc[i] <- roc$auc
  }
})
# save the result accordingly
saveRDS(list(prauc, rocauc,scores), 'paper_result/hybrid_stability_score_2.rds')


