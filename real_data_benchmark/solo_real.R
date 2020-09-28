library(reticulate)
library(Matrix)
library(PRROC)
library(RcppCNPy)
np <- import("numpy")
set.seed(2020)
############################################################################################################################
# calculate auprc and auroc by doublet scores of 16 benchmark datasets, calculated by solo command
############################################################################################################################
# read calculated doublet scores from local
# change the location accordingly
files <- c('real_data/solo_out/pbmc-ch/softmax_scores.npy', 'real_data/solo_out/cline-ch/softmax_scores.npy',
           'real_data/solo_out/mkidney-ch/softmax_scores.npy', 'real_data/solo_out/hm-12k/softmax_scores.npy',
           'real_data/solo_out/hm-6k/softmax_scores.npy', 'real_data/solo_out/pbmc-1A-dm/softmax_scores.npy',
           'real_data/solo_out/pbmc-1B-dm/softmax_scores.npy', 'real_data/solo_out/pbmc-1C-dm/softmax_scores.npy',
           'real_data/solo_out/pbmc-2ctrl-dm/softmax_scores.npy', 'real_data/solo_out/pbmc-2stim-dm/softmax_scores.npy',
           'real_data/solo_out/J293t-dm/softmax_scores.npy', 'real_data/solo_out/pdx-MULTI/softmax_scores.npy',
           'real_data/solo_out/HMEC-orig-MULTI/softmax_scores.npy', 'real_data/solo_out/HMEC-rep-MULTI/softmax_scores.npy',
           'real_data/solo_out/HEK-HMEC-MULTI/softmax_scores.npy', 'real_data/solo_out/nuc-MULTI/softmax_scores.npy')

# read real data from local
# change the location accordingly
locs <- c('real_data/pbmc-ch.rds', 'real_data/cline-ch.rds', 'real_data/mkidney-ch.rds', 'real_data/hm-12k.rds', 
          'real_data/hm-6k.rds', 'real_data/pbmc-1A-dm.rds', 'real_data/pbmc-1B-dm.rds', 'real_data/pbmc-1C-dm.rds',
          'real_data/pbmc-2ctrl-dm.rds', 'real_data/pbmc-2stim-dm.rds', 'real_data/J293t-dm.rds', 'real_data/pdx-MULTI.rds',
          'real_data/HMEC-orig-MULTI.rds', 'real_data/HMEC-rep-MULTI.rds', 'real_data/HEK-HMEC-MULTI.rds', 
          'real_data/nuc-MULTI.rds')
prs <- c()
rocs <- c()
# loop over each dataset and its scores
for(i in length(locs)){
  score <- as.numeric(np$load(files[i])); length(score); hist(score)
  data <- readRDS(locs[i])
  label <- data[[2]]; table(label)
  label <- ifelse(label == 'doublet', 1, 0); table(label)
  # auc
  fg <- score[label==1]
  bg <- score[label==0]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T); pr$auc.integral
  roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc
  prs[i] <- pr
  rocs[i] <- roc[i]
}
# save the result accordingly

##########################################################################################################################
# pr, recall, and tnr under 10%, 20%, and 40% identification rates
##########################################################################################################################
files <- c('real_data/solo_out/pbmc-ch/softmax_scores.npy', 'real_data/solo_out/cline-ch/softmax_scores.npy',
          'real_data/solo_out/mkidney-ch/softmax_scores.npy', 'real_data/solo_out/hm-12k/softmax_scores.npy',
          'real_data/solo_out/hm-6k/softmax_scores.npy', 'real_data/solo_out/pbmc-1A-dm/softmax_scores.npy',
          'real_data/solo_out/pbmc-1B-dm/softmax_scores.npy', 'real_data/solo_out/pbmc-1C-dm/softmax_scores.npy',
          'real_data/solo_out/pbmc-2ctrl-dm/softmax_scores.npy', 'real_data/solo_out/pbmc-2stim-dm/softmax_scores.npy',
          'real_data/solo_out/J293t-dm/softmax_scores.npy', 'real_data/solo_out/pdx-MULTI/softmax_scores.npy',
          'real_data/solo_out/HMEC-orig-MULTI/softmax_scores.npy', 'real_data/solo_out/HMEC-rep-MULTI/softmax_scores.npy',
          'real_data/solo_out/HEK-HMEC-MULTI/softmax_scores.npy', 'real_data/solo_out/nuc-MULTI/softmax_scores.npy')
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
    file <- files[i]
    # calculate threshold based on identification rate
    score <- as.numeric(np$load(file))
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
files <- c('real_data/solo_out/pbmc-ch/softmax_scores.npy', 'real_data/solo_out/cline-ch/softmax_scores.npy',
           'real_data/solo_out/mkidney-ch/softmax_scores.npy', 'real_data/solo_out/hm-12k/softmax_scores.npy',
           'real_data/solo_out/hm-6k/softmax_scores.npy', 'real_data/solo_out/pbmc-1A-dm/softmax_scores.npy',
           'real_data/solo_out/pbmc-1B-dm/softmax_scores.npy', 'real_data/solo_out/pbmc-1C-dm/softmax_scores.npy',
           'real_data/solo_out/pbmc-2ctrl-dm/softmax_scores.npy', 'real_data/solo_out/pbmc-2stim-dm/softmax_scores.npy',
           'real_data/solo_out/J293t-dm/softmax_scores.npy', 'real_data/solo_out/pdx-MULTI/softmax_scores.npy',
           'real_data/solo_out/HMEC-orig-MULTI/softmax_scores.npy', 'real_data/solo_out/HMEC-rep-MULTI/softmax_scores.npy',
           'real_data/solo_out/HEK-HMEC-MULTI/softmax_scores.npy', 'real_data/solo_out/nuc-MULTI/softmax_scores.npy')
score.list <- list()
for(file in files){
  score <- as.numeric(np$load(file))
  score.list <- append(score.list, list(score))
}
# 16 benchmark datasets locations
locs <- c('real_data/pbmc-ch.rds', 'real_data/cline-ch.rds', 'real_data/mkidney-ch.rds', 'real_data/hm-12k.rds', 
          'real_data/hm-6k.rds', 'real_data/pbmc-1A-dm.rds', 'real_data/pbmc-1B-dm.rds', 'real_data/pbmc-1C-dm.rds',
          'real_data/pbmc-2ctrl-dm.rds', 'real_data/pbmc-2stim-dm.rds', 'real_data/J293t-dm.rds', 'real_data/pdx-MULTI.rds',
          'real_data/HMEC-orig-MULTI.rds', 'real_data/HMEC-rep-MULTI.rds', 'real_data/HEK-HMEC-MULTI.rds', 
          'real_data/nuc-MULTI.rds')
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

#####################################################################################################################################
# running time on 16 benchmark datasets
# stability on subsamples of real datasets
# The previous studies have been done by solo command
####################################################################################################################################
