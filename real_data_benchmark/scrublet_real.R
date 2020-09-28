library(reticulate)
use_python('/home/nxi/anaconda3/bin', required = T)
py_config()
library(Matrix)
library(Seurat)
library(dplyr)
library(PRROC)
library(pbapply)
setwd('/media/nxi/nxi/doublet')
source('utility.R')
getwd()

# read python module
scr <- import('scrublet')
scipy.io <- import('scipy.io')
np <- import('numpy')
os <- import('os')
score.list <- list()

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
count <- data[[1]]; dim(count)
label <- data[[2]]; table(label)
label <- ifelse(label == 'doublet', 1, 0); table(label)
doublet.rate <- sum(label==1) / length(label); doublet.rate

system.time({
result <- scr$Scrublet(counts_matrix = t(count), expected_doublet_rate = doublet.rate, random_state = 10L)
results <- result$scrub_doublets(min_counts=2, 
                                 min_cells=3, 
                                 min_gene_variability_pctl=85, 
                                 n_prin_comps=30L)
})
score <- as.vector(results[[1]])
score.list <- append(score.list, list(score))

# pr-auc
fg <- score[label==1]
bg <- score[label==0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T); pr$auc.integral
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc

saveRDS(score.list, 'paper_result/scrublet_real_score.rds')


# pr
score.list <- readRDS('paper_result/scrublet_real_score.rds')
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

#####################################################################################
# running time
#####################################################################################
# read data; parameter setting
datasets <- c(
'real_data/pbmc-ch.rds',
'real_data/cline-ch.rds',
'real_data/mkidney-ch.rds',
'real_data/hm-12k.rds',
'real_data/hm-6k.rds',
'real_data/pbmc-1A-dm.rds',
'real_data/pbmc-1B-dm.rds',
'real_data/pbmc-1C-dm.rds',
'real_data/pbmc-2ctrl-dm.rds',
'real_data/pbmc-2stim-dm.rds',
'real_data/J293t-dm.rds',
'real_data/pdx-MULTI.rds',
'real_data/HMEC-orig-MULTI.rds',
'real_data/HMEC-rep-MULTI.rds',
'real_data/HEK-HMEC-MULTI.rds',
'real_data/nuc-MULTI.rds'
)
times <- pbsapply(datasets, function(loc){
  data <- readRDS(loc)
  count <- data[[1]]; dim(count)
  label <- data[[2]]; table(label)
  label <- ifelse(label == 'doublet', 1, 0); table(label)
  doublet.rate <- sum(label==1) / length(label); doublet.rate
  time <- system.time({
    result <- scr$Scrublet(counts_matrix = t(count), expected_doublet_rate = doublet.rate, random_state = 10L)
    results <- result$scrub_doublets(min_counts=2, 
                                     min_cells=3, 
                                     min_gene_variability_pctl=85, 
                                     n_prin_comps=30L)
  })
  return(time[3])
}, simplify = T); times

##################################################################################
# stability
##################################################################################
loc <- 'real_data/pbmc-2ctrl-dm_sub.rds'
data <- readRDS(loc)
counts <- data[[1]]
labels <- data[[2]]
prauc <- c()
rocauc <- c()
scores <- list()
system.time({
for(i in 1:length(counts)){
  print(i)
  count <- counts[[i]]
  label <- labels[[i]]
  label <- ifelse(label == 'doublet', 1, 0); table(label)
  doublet.rate <- sum(label==1) / length(label); doublet.rate
  result <- scr$Scrublet(counts_matrix = t(count), expected_doublet_rate = doublet.rate, random_state = 10L)
  results <- result$scrub_doublets(min_counts=2, 
                                   min_cells=3, 
                                   min_gene_variability_pctl=85, 
                                   n_prin_comps=30L)
  score <- as.vector(results[[1]])
  scores <- append(scores, list(score))
  fg <- score[label==1]
  bg <- score[label==0]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  prauc[i] <- pr$auc.integral
  roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  rocauc[i] <- roc$auc
}
})
boxplot(prauc)
boxplot(rocauc)
sd(prauc)
sd(rocauc)

saveRDS(list(prauc, rocauc,scores), 'paper_result/scrublet_stability_score_2.rds')









