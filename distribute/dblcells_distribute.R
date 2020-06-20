library(Matrix)
library(Seurat)
library(SingleCellExperiment)
library(scran)
library(PRROC)
library(pbapply)
setwd('/media/nxi/nxi/doublet')
getwd()

# read data; parameter setting
system.time(sim.data <- readRDS("real_data/pbmc-2ctrl-dm.rds"))
score.list <- list()
#####################################
# batch=1
#####################################
data <- sim.data[[1]]
system.time({
  score1 <- doubletCells(data)
})
# save merge prediction result
score.list <- append(score.list, list(score1))

# analyze result
labels <- sim.data[[2]]; table(labels)
labels <- ifelse(labels == 'doublet', 1, 0); table(labels)
fg <- score1[labels == 1]
bg <- score1[labels == 0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr1 <- pr$auc.integral; pr1


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
  score1 <- doubletCells(x)
  return(score1)
})
score2 <- unlist(score2.list)
labels2 <- c(label1, label2)

# analyze result
labels2 <- ifelse(labels2 == 'doublet', 1, 0); table(labels2)
fg <- score2[labels2 == 1]
bg <- score2[labels2 == 0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr2 <- pr$auc.integral; pr2
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc

score.list <- append(score.list, list(score2.list))
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
  score1 <- doubletCells(x)
  return(score1)
})
score4 <- unlist(score4.list, use.names = F)

# analyze result
labels4 <- ifelse(labels4 == 'doublet', 1, 0); table(labels4)
fg <- score4[labels4 == 1]
bg <- score4[labels4 == 0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr4 <- pr$auc.integral; pr4
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc

score.list <- append(score.list, list(score4.list))
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
  score1 <- doubletCells(x)
  return(score1)
})
score6 <- unlist(score6.list, use.names = F)

# analyze result
labels6 <- ifelse(labels6 == 'doublet', 1, 0); table(labels6)
fg <- score6[labels6 == 1]
bg <- score6[labels6 == 0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr6 <- pr$auc.integral; pr6
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc

score.list <- append(score.list, list(score6.list))
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
  score1 <- doubletCells(x)
  return(score1)
})
score8 <- unlist(score8.list, use.names = F)

# analyze result
labels8 <- ifelse(labels8 == 'doublet', 1, 0); table(labels8)
fg <- score8[labels8 == 1]
bg <- score8[labels8 == 0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr8 <- pr$auc.integral; pr8
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc

score.list <- append(score.list, list(score8.list))
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
  score1 <- doubletCells(x)
  return(score1)
})
score10 <- unlist(score10.list, use.names = F)

# analyze result
labels10 <- ifelse(labels10 == 'doublet', 1, 0); table(labels10)
fg <- score10[labels10 == 1]
bg <- score10[labels10 == 0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr10 <- pr$auc.integral; pr10
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc

score.list <- append(score.list, list(score10.list))

# save all scores
system.time(saveRDS(score.list, "paper_result/dblcells_distribute_pbmc-2ctrl-dm.rds"))







