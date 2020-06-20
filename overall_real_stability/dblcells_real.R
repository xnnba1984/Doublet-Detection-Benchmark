library(Matrix)
library(Seurat)
library(SingleCellExperiment)
library(scran)
library(PRROC)
library(pbapply)
setwd('/media/nxi/nxi/doublet')
getwd()

score.list <- list()

# read data; parameter setting
loc <- 'real_data/pbmc-ch.rds'
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
#loc <- 'real_data/nuc-MULTI.rds'

data <- readRDS(loc)
count <- data[[1]]; dim(count)
label <- data[[2]]; table(label)
label <- ifelse(label == 'doublet', 1, 0); table(label)
doublet.rate <- sum(label==1) / length(label); doublet.rate

set.seed(10)
sce <- SingleCellExperiment(list(counts=count))
clusters <- quickCluster(sce, method="igraph", min.mean=0.1)
table(clusters)
sce <- computeSumFactors(sce, clusters=clusters, min.mean=0.1)
summary(sizeFactors(sce))
sce <- normalize(sce)
assayNames(sce)
tech.trend <- makeTechTrend(x=sce)
fit <- trendVar(sce, use.spikes=FALSE)
plot(fit$mean, fit$var, pch=16, 
     xlab="Mean log-expression",
     ylab="Variance of log-expression")
curve(tech.trend(x), add=TRUE, col="red")
sce <- denoisePCA(sce, technical=tech.trend)
ncol(reducedDim(sce))
score <- doubletCells(sce)
hist(score)

system.time(score <- doubletCells(count))
score.list <- append(score.list, list(score))

# pr-auc
fg <- score[label==1]
bg <- score[label==0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T); pr$auc.integral
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc

saveRDS(score.list, 'paper_result/dblcells_real_score.rds')


# pr
score.list <- readRDS('paper_result/dblcells_real_score.rds')
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
    sce <- SingleCellExperiment(list(counts=count))
    score <- doubletCells(sce)
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
    score <- doubletCells(count)
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
saveRDS(list(prauc, rocauc,scores), 'paper_result/dblcells_stability_score_2.rds')




