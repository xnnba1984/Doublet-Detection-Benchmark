#####################################################################################################################################
# time for f in *.loom; do solo solo_params_example.json $f -o 'output_psudotime_DE/'$f; done;
####################################################################################################################################
library(reticulate)
use_python('/home/nxi/anaconda3/bin', required = T)
library(Matrix)
library(PRROC)
library(RcppCNPy)
setwd('/media/nxi/nxi/doublet')
getwd()
np <- import("numpy")

#file <- 'real_data/solo_out/pbmc-ch/softmax_scores.npy'
#file <- 'real_data/solo_out/cline-ch/softmax_scores.npy'
#file <- 'real_data/solo_out/mkidney-ch/softmax_scores.npy'
#file <- 'real_data/solo_out/hm-12k/softmax_scores.npy'
#file <- 'real_data/solo_out/hm-6k/softmax_scores.npy'
#file <- 'real_data/solo_out/pbmc-1A-dm/softmax_scores.npy'
#file <- 'real_data/solo_out/pbmc-1B-dm/softmax_scores.npy'
#file <- 'real_data/solo_out/pbmc-1C-dm/softmax_scores.npy'
#file <- 'real_data/solo_out/pbmc-2ctrl-dm/softmax_scores.npy'
#file <- 'real_data/solo_out/pbmc-2stim-dm/softmax_scores.npy'
#file <- 'real_data/solo_out/J293t-dm/softmax_scores.npy'
#file <- 'real_data/solo_out/pdx-MULTI/softmax_scores.npy'
#file <- 'real_data/solo_out/HMEC-orig-MULTI/softmax_scores.npy'
#file <- 'real_data/solo_out/HMEC-rep-MULTI/softmax_scores.npy'
#file <- 'real_data/solo_out/HEK-HMEC-MULTI/softmax_scores.npy'
file <- 'real_data/solo_out/nuc-MULTI/softmax_scores.npy'


# data reading
score <- as.numeric(np$load(file)); length(score); hist(score)

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
label <- data[[2]]; table(label)
label <- ifelse(label == 'doublet', 1, 0); table(label)
# auc
fg <- score[label==1]
bg <- score[label==0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T); pr$auc.integral
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc

###################################################################################
files <- c('real_data/solo_out/pbmc-ch/softmax_scores.npy',
 'real_data/solo_out/cline-ch/softmax_scores.npy',
 'real_data/solo_out/mkidney-ch/softmax_scores.npy',
 'real_data/solo_out/hm-12k/softmax_scores.npy',
 'real_data/solo_out/hm-6k/softmax_scores.npy',
 'real_data/solo_out/pbmc-1A-dm/softmax_scores.npy',
 'real_data/solo_out/pbmc-1B-dm/softmax_scores.npy',
 'real_data/solo_out/pbmc-1C-dm/softmax_scores.npy',
 'real_data/solo_out/pbmc-2ctrl-dm/softmax_scores.npy',
 'real_data/solo_out/pbmc-2stim-dm/softmax_scores.npy',
 'real_data/solo_out/J293t-dm/softmax_scores.npy',
 'real_data/solo_out/pdx-MULTI/softmax_scores.npy',
 'real_data/solo_out/HMEC-orig-MULTI/softmax_scores.npy',
 'real_data/solo_out/HMEC-rep-MULTI/softmax_scores.npy',
 'real_data/solo_out/HEK-HMEC-MULTI/softmax_scores.npy',
 'real_data/solo_out/nuc-MULTI/softmax_scores.npy')
score.list <- list()
for(file in files){
  score <- as.numeric(np$load(file))
  score.list <- append(score.list, list(score))
}
locs <- c(
  loc <- 'real_data/pbmc-ch.rds',
  loc <- 'real_data/cline-ch.rds',
  loc <- 'real_data/mkidney-ch.rds',
  #loc <- 'real_data/hm-12k.rds'
  loc <- 'real_data/hm-6k.rds',
  loc <- 'real_data/pbmc-1A-dm.rds',
  loc <- 'real_data/pbmc-1B-dm.rds',
  loc <- 'real_data/pbmc-1C-dm.rds',
  #loc <- 'real_data/pbmc-2ctrl-dm.rds'
  loc <- 'real_data/pbmc-2stim-dm.rds',
  #loc <- 'real_data/J293t-dm.rds'
  loc <- 'real_data/pdx-MULTI.rds',
  loc <- 'real_data/HMEC-orig-MULTI.rds',
  loc <- 'real_data/HMEC-rep-MULTI.rds',
  loc <- 'real_data/HEK-HMEC-MULTI.rds'
  #loc <- 'real_data/nuc-MULTI.rds'
)
n <- c(1,2,3,5,6,7,8,10,12,13,14,15)
d <- c(7872,2822,8417,2813,1223,1493,1961,4077,3479,18007,8448,3124)
precisions <- c()
recalls <- c()
tnrs <- c()
f1s <- c()

for(i in 1:length(locs)){
  print(i)
  data <- readRDS(locs[i])
  label <- data[[2]]; table(label)
  label <- ifelse(label == 'doublet', 1, 0); table(label)
  score <- score.list[[n[i]]]
  thresh <- sort(score, decreasing = T)[d[i]]
  pred <- score > thresh; table(pred)
  # result
  tp <- sum(pred[which(label==1)]==1); tp
  fp <- sum(pred[which(label==0)]==1); fp
  fn <- sum(pred[which(label==1)]==0); fn
  tn <- sum(pred[which(label==0)]==0); tn
  
  precision <- tp/(tp + fp); precision
  recall <- tp/(tp + fn); recall
  tnr <- tn/(tn + fp); tnr
  f1 <- 2 * precision * recall / (precision + recall); f1
  
  precisions[i] <- precision
  recalls[i] <- recall
  tnrs[i] <- tnr
  f1s[i] <- f1
}
names(precisions) <- locs; precisions
names(recalls) <- locs; recalls
names(tnrs) <- locs; tnrs
names(f1s) <- locs; f1s




