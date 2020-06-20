library(Matrix)
library(Seurat)
library(scds)
library(SingleCellExperiment)
library(PRROC)
library(pbapply)
setwd('/media/nxi/nxi/doublet')
source('utility.R')
getwd()# read data; parameter setting

system.time(sim.data <- readRDS("real_data/pbmc-2ctrl-dm.rds"))
score.cxds.list <- list()
score.bcds.list <- list()
score.hybrid.list <- list()
#####################################
# batch=1
#####################################
data <- sim.data[[1]]; dim(data)
sce <- SingleCellExperiment(assays = list(counts = data))

# co-expression, boost, hybrid
system.time({sce <- cxds(sce); sce <- bcds(sce); sce <- cxds_bcds_hybrid(sce)})
CD <- colData(sce)

# save scores
score1.cxds <- CD$cxds_score
score1.bcds <- CD$bcds_score
score1.hybrid <- CD$hybrid_score

# analyze result
labels <- sim.data[[2]]; table(labels)
labels <- ifelse(labels == 'doublet', 1, 0); table(labels)
fg <- score1.hybrid[labels == 1]
bg <- score1.hybrid[labels == 0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr1 <- pr$auc.integral; pr1

# save merge prediction result
score.cxds.list <- append(score.cxds.list, list(score1.cxds))
score.bcds.list <- append(score.bcds.list, list(score1.bcds))
score.hybrid.list <- append(score.hybrid.list, list(score1.hybrid))

#####################################
# batch=2
#####################################
data <- sim.data[[1]]; dim(data)
labels <- sim.data[[2]]; table(labels)
set.seed(10)
index <- sample(1:dim(data)[2], round(dim(data)[2]/2))
data1 <- data[,index]; dim(data1)
data2 <- data[,-index]; dim(data2)
label1 <- labels[index]; table(label1)
label2 <- labels[-index]; table(label2)

score2.list <- pblapply(list(data1, data2), function(x){
  sce <- SingleCellExperiment(assays = list(counts = x))
  sce <- cxds(sce); sce <- bcds(sce); sce <- cxds_bcds_hybrid(sce)
  CD <- colData(sce)
  return(list(CD$cxds_score,CD$bcds_score,CD$hybrid_score))
})

score2.cxds <- c(score2.list[[1]][[1]], score2.list[[2]][[1]])
score2.bcds <- c(score2.list[[1]][[2]], score2.list[[2]][[2]])
score2.hybrid <- c(score2.list[[1]][[3]], score2.list[[2]][[3]])

labels2 <- c(label1, label2)

# analyze result
labels2 <- ifelse(labels2 == 'doublet', 1, 0); table(labels2)
fg <- score2.cxds[labels2 == 1]
bg <- score2.cxds[labels2 == 0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr2 <- pr$auc.integral; pr2
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc

# save merge prediction result
score.bcds.list <- append(score.bcds.list, list(score2.bcds))
score.cxds.list <- append(score.cxds.list, list(score2.cxds))
score.hybrid.list <- append(score.hybrid.list, list(score2.hybrid))

#####################################
# batch=4
#####################################
data <- sim.data[[1]]; dim(data)
labels <- sim.data[[2]]; table(labels)
set.seed(10)
index <- sample(1:dim(data)[2], dim(data)[2])
chunk <- split(index, cut(seq_along(index), 4, labels = FALSE)) 
data.chunk <- pblapply(chunk, function(x){
  return(data[,x])
})
labels4 <- labels[unlist(chunk,use.names = F)]

score4.list <- pblapply(data.chunk, function(x){
  sce <- SingleCellExperiment(assays = list(counts = x))
  sce <- cxds(sce); sce <- bcds(sce); sce <- cxds_bcds_hybrid(sce)
  CD <- colData(sce)
  return(list(CD$cxds_score,CD$bcds_score,CD$hybrid_score))
})
score4.cxds <- c(score4.list[[1]][[1]], score4.list[[2]][[1]], score4.list[[3]][[1]], score4.list[[4]][[1]])
score4.bcds <- c(score4.list[[1]][[2]], score4.list[[2]][[2]], score4.list[[3]][[2]], score4.list[[4]][[2]])
score4.hybrid <- c(score4.list[[1]][[3]], score4.list[[2]][[3]], score4.list[[3]][[3]], score4.list[[4]][[3]])

# analyze result
labels4 <- ifelse(labels4 == 'doublet', 1, 0); table(labels4)
fg <- score4.cxds[labels4 == 1]
bg <- score4.cxds[labels4 == 0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr4 <- pr$auc.integral; pr4
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc

# save merge prediction result
score.bcds.list <- append(score.bcds.list, list(score4.bcds))
score.cxds.list <- append(score.cxds.list, list(score4.cxds))
score.hybrid.list <- append(score.hybrid.list, list(score4.hybrid))

#####################################
# batch=6
#####################################
data <- sim.data[[1]]; dim(data)
labels <- sim.data[[2]]; table(labels)
set.seed(10)
index <- sample(1:dim(data)[2], dim(data)[2])
chunk <- split(index, cut(seq_along(index), 6, labels = FALSE)) 
data.chunk <- pblapply(chunk, function(x){
  return(data[,x])
})
labels6 <- labels[unlist(chunk,use.names = F)]

score6.list <- pblapply(data.chunk, function(x){
  sce <- SingleCellExperiment(assays = list(counts = x))
  sce <- cxds(sce); sce <- bcds(sce); sce <- cxds_bcds_hybrid(sce)
  CD <- colData(sce)
  return(list(CD$cxds_score,CD$bcds_score,CD$hybrid_score))
})
score6.cxds <- c(score6.list[[1]][[1]], score6.list[[2]][[1]], score6.list[[3]][[1]], score6.list[[4]][[1]],
                 score6.list[[5]][[1]], score6.list[[6]][[1]])
score6.bcds <- c(score6.list[[1]][[2]], score6.list[[2]][[2]], score6.list[[3]][[2]], score6.list[[4]][[2]],
                 score6.list[[5]][[2]], score6.list[[6]][[2]])
score6.hybrid <- c(score6.list[[1]][[3]], score6.list[[2]][[3]], score6.list[[3]][[3]], score6.list[[4]][[3]],
                   score6.list[[5]][[3]], score6.list[[6]][[3]])

# analyze result
labels6 <- ifelse(labels6 == 'doublet', 1, 0); table(labels6)
fg <- score6.cxds[labels6 == 1]
bg <- score6.cxds[labels6 == 0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr6 <- pr$auc.integral; pr6
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc

# save merge prediction result
score.bcds.list <- append(score.bcds.list, list(score6.bcds))
score.cxds.list <- append(score.cxds.list, list(score6.cxds))
score.hybrid.list <- append(score.hybrid.list, list(score6.hybrid))

#####################################
# batch=8
#####################################
data <- sim.data[[1]]; dim(data)
labels <- sim.data[[2]]; table(labels)
set.seed(10)
index <- sample(1:dim(data)[2], dim(data)[2])
chunk <- split(index, cut(seq_along(index), 8, labels = FALSE)) 
data.chunk <- pblapply(chunk, function(x){
  return(data[,x])
})
labels8 <- labels[unlist(chunk,use.names = F)]

score8.list <- pblapply(data.chunk, function(x){
  sce <- SingleCellExperiment(assays = list(counts = x))
  sce <- cxds(sce); sce <- bcds(sce); sce <- cxds_bcds_hybrid(sce)
  CD <- colData(sce)
  return(list(CD$cxds_score,CD$bcds_score,CD$hybrid_score))
})
score8.cxds <- c(score8.list[[1]][[1]], score8.list[[2]][[1]], score8.list[[3]][[1]], score8.list[[4]][[1]],
                 score8.list[[5]][[1]], score8.list[[6]][[1]],score8.list[[7]][[1]], score8.list[[8]][[1]])
score8.bcds <- c(score8.list[[1]][[2]], score8.list[[2]][[2]], score8.list[[3]][[2]], score8.list[[4]][[2]],
                 score8.list[[5]][[2]], score8.list[[6]][[2]],score8.list[[7]][[2]], score8.list[[8]][[2]])
score8.hybrid <- c(score8.list[[1]][[3]], score8.list[[2]][[3]], score8.list[[3]][[3]], score8.list[[4]][[3]],
                   score8.list[[5]][[3]], score8.list[[6]][[3]],score8.list[[7]][[3]], score8.list[[8]][[3]])

# analyze result
labels8 <- ifelse(labels8 == 'doublet', 1, 0); table(labels8)
fg <- score8.cxds[labels8 == 1]
bg <- score8.cxds[labels8 == 0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr8 <- pr$auc.integral; pr8
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc

# save merge prediction result
score.bcds.list <- append(score.bcds.list, list(score8.bcds))
score.cxds.list <- append(score.cxds.list, list(score8.cxds))
score.hybrid.list <- append(score.hybrid.list, list(score8.hybrid))

#####################################
# batch=10
#####################################
data <- sim.data[[1]]; dim(data)
labels <- sim.data[[2]]; table(labels)
set.seed(10)
index <- sample(1:dim(data)[2], dim(data)[2])
chunk <- split(index, cut(seq_along(index), 10, labels = FALSE)) 
data.chunk <- pblapply(chunk, function(x){
  return(data[,x])
})
labels10 <- labels[unlist(chunk,use.names = F)]

score10.list <- pblapply(data.chunk, function(x){
  sce <- SingleCellExperiment(assays = list(counts = x))
  sce <- cxds(sce); sce <- bcds(sce); sce <- cxds_bcds_hybrid(sce)
  CD <- colData(sce)
  return(list(CD$cxds_score,CD$bcds_score,CD$hybrid_score))
})
score10.cxds <- c(score10.list[[1]][[1]], score10.list[[2]][[1]], score10.list[[3]][[1]], score10.list[[4]][[1]],
                 score10.list[[5]][[1]], score10.list[[6]][[1]],score10.list[[7]][[1]], score10.list[[8]][[1]],
                 score10.list[[9]][[1]], score10.list[[10]][[1]])
score10.bcds <- c(score10.list[[1]][[2]], score10.list[[2]][[2]], score10.list[[3]][[2]], score10.list[[4]][[2]],
                 score10.list[[5]][[2]], score10.list[[6]][[2]],score10.list[[7]][[2]], score10.list[[8]][[2]],
                 score10.list[[9]][[2]], score10.list[[10]][[2]])
score10.hybrid <- c(score10.list[[1]][[3]], score10.list[[2]][[3]], score10.list[[3]][[3]], score10.list[[4]][[3]],
                   score10.list[[5]][[3]], score10.list[[6]][[3]],score10.list[[7]][[3]], score10.list[[8]][[3]],
                   score10.list[[9]][[3]], score10.list[[10]][[3]])

# analyze result
labels10 <- ifelse(labels10 == 'doublet', 1, 0); table(labels10)
fg <- score10.hybrid[labels10 == 1]
bg <- score10.hybrid[labels10 == 0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr10 <- pr$auc.integral; pr10
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc

# save merge prediction result
score.bcds.list <- append(score.bcds.list, list(score10.bcds))
score.cxds.list <- append(score.cxds.list, list(score10.cxds))
score.hybrid.list <- append(score.hybrid.list, list(score10.hybrid))

score.list <- list(score.cxds.list, score.bcds.list, score.hybrid.list)
# save all scores
system.time(saveRDS(score.list, "sim_result/scds_distribute_pbmc-2ctrl-dm.rds"))









