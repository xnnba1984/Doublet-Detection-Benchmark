library(reticulate)
use_python('/home/nxi/anaconda3/bin', required = T)
library(Matrix)
library(Seurat)
library(PRROC)
library(pbapply)
setwd('/media/nxi/nxi/doublet')
getwd()
np <- import("numpy")

# read data; parameter setting
system.time(sim.data <- readRDS("real_data/pbmc-2ctrl-dm.rds"))
score.list <- list()

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

# save sim data to loom
library(loomR)
system.time({
  create('paper_sim/loom/distribute_2_1_.loom', data1, transpose = T)
  create('paper_sim/loom/distribute_2_2_.loom', data2, transpose = T)
})

file1 <- 'paper_sim/loom/output_distribute_/distribute_2_1_.loom/softmax_scores.npy'
score1 <- as.numeric(np$load(file1)); length(score1); hist(score1)
file2 <- 'paper_sim/loom/output_distribute_/distribute_2_2_.loom/softmax_scores.npy'
score2 <- as.numeric(np$load(file2)); length(score2); hist(score2)
score <- c(score1, score2)
label <- c(label1, label2)
label <- ifelse(label == 'doublet', 1, 0); table(label)
fg <- score[label == 1]
bg <- score[label == 0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr2 <- pr$auc.integral; pr2
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc

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
labels4 <- labels[unlist(chunk,use.names = F)]; table(labels4)

system.time({
  for(i in 1:length(data.chunk)){
    print(paste('===================', i))
    counts <- data.chunk[[i]]; dim(counts)
    counts <- counts[which(rowSums(counts) != 0),]; dim(counts)
    print(dim(counts))
    create(paste(paste('paper_sim/loom/distribute_4_', i, sep='_'), 'loom', sep='.'), counts, transpose = T)
  }
})

file1 <- 'paper_sim/loom/output_distribute_/distribute_4__1.loom/softmax_scores.npy'
score1 <- as.numeric(np$load(file1)); length(score1); hist(score1)
file2 <- 'paper_sim/loom/output_distribute_/distribute_4__2.loom/softmax_scores.npy'
score2 <- as.numeric(np$load(file2)); length(score2); hist(score2)
file3 <- 'paper_sim/loom/output_distribute_/distribute_4__3.loom/softmax_scores.npy'
score3 <- as.numeric(np$load(file3)); length(score3); hist(score3)
file4 <- 'paper_sim/loom/output_distribute_/distribute_4__4.loom/softmax_scores.npy'
score4 <- as.numeric(np$load(file4)); length(score4); hist(score4)

score <- c(score1, score2, score3, score4)
label <- ifelse(labels4 == 'doublet', 1, 0); table(label)
fg <- score[label == 1]
bg <- score[label == 0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr2 <- pr$auc.integral; pr2
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc

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

system.time({
  for(i in 1:length(data.chunk)){
    print(paste('===================', i))
    counts <- data.chunk[[i]]; dim(counts)
    counts <- counts[which(rowSums(counts) != 0),]; dim(counts)
    print(dim(counts))
    create(paste(paste('paper_sim/loom/distribute_6_', i, sep='_'), 'loom', sep='.'), counts, transpose = T)
  }
})

file1 <- 'paper_sim/loom/output_distribute_/distribute_6__1.loom/softmax_scores.npy'
score1 <- as.numeric(np$load(file1)); length(score1); hist(score1)
file2 <- 'paper_sim/loom/output_distribute_/distribute_6__2.loom/softmax_scores.npy'
score2 <- as.numeric(np$load(file2)); length(score2); hist(score2)
file3 <- 'paper_sim/loom/output_distribute_/distribute_6__3.loom/softmax_scores.npy'
score3 <- as.numeric(np$load(file3)); length(score3); hist(score3)
file4 <- 'paper_sim/loom/output_distribute_/distribute_6__4.loom/softmax_scores.npy'
score4 <- as.numeric(np$load(file4)); length(score4); hist(score4)
file5 <- 'paper_sim/loom/output_distribute_/distribute_6__5.loom/softmax_scores.npy'
score5 <- as.numeric(np$load(file5)); length(score5); hist(score5)
file6 <- 'paper_sim/loom/output_distribute_/distribute_6__6.loom/softmax_scores.npy'
score6 <- as.numeric(np$load(file6)); length(score6); hist(score6)

score <- c(score1, score2, score3, score4, score5, score6)
label <- ifelse(labels6 == 'doublet', 1, 0); table(label)
fg <- score[label == 1]
bg <- score[label == 0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr2 <- pr$auc.integral; pr2
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc

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

system.time({
  for(i in 1:length(data.chunk)){
    print(paste('===================', i))
    counts <- data.chunk[[i]]; dim(counts)
    counts <- counts[which(rowSums(counts) != 0),]; dim(counts)
    print(dim(counts))
    create(paste(paste('paper_sim/loom/distribute_8_', i, sep='_'), 'loom', sep='.'), counts, transpose = T)
  }
})

file1 <- 'paper_sim/loom/output_distribute_/distribute_8__1.loom/softmax_scores.npy'
score1 <- as.numeric(np$load(file1)); length(score1); hist(score1)
file2 <- 'paper_sim/loom/output_distribute_/distribute_8__2.loom/softmax_scores.npy'
score2 <- as.numeric(np$load(file2)); length(score2); hist(score2)
file3 <- 'paper_sim/loom/output_distribute_/distribute_8__3.loom/softmax_scores.npy'
score3 <- as.numeric(np$load(file3)); length(score3); hist(score3)
file4 <- 'paper_sim/loom/output_distribute_/distribute_8__4.loom/softmax_scores.npy'
score4 <- as.numeric(np$load(file4)); length(score4); hist(score4)
file5 <- 'paper_sim/loom/output_distribute_/distribute_8__5.loom/softmax_scores.npy'
score5 <- as.numeric(np$load(file5)); length(score5); hist(score5)
file6 <- 'paper_sim/loom/output_distribute_/distribute_8__6.loom/softmax_scores.npy'
score6 <- as.numeric(np$load(file6)); length(score6); hist(score6)
file7 <- 'paper_sim/loom/output_distribute_/distribute_8__7.loom/softmax_scores.npy'
score7 <- as.numeric(np$load(file7)); length(score7); hist(score7)
file8 <- 'paper_sim/loom/output_distribute_/distribute_8__8.loom/softmax_scores.npy'
score8 <- as.numeric(np$load(file8)); length(score8); hist(score8)

score <- c(score1, score2, score3, score4, score5, score6, score7, score8)
label <- ifelse(labels8 == 'doublet', 1, 0); table(label)
fg <- score[label == 1]
bg <- score[label == 0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr2 <- pr$auc.integral; pr2
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc

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

system.time({
  for(i in 1:length(data.chunk)){
    print(paste('===================', i))
    counts <- data.chunk[[i]]; dim(counts)
    counts <- counts[which(rowSums(counts) != 0),]; dim(counts)
    print(dim(counts))
    create(paste(paste('paper_sim/loom/distribute_10_', i, sep='_'), 'loom', sep='.'), counts, transpose = T)
  }
})

file1 <- 'paper_sim/loom/output_distribute_/distribute_10__1.loom/softmax_scores.npy'
score1 <- as.numeric(np$load(file1)); length(score1); hist(score1)
file2 <- 'paper_sim/loom/output_distribute_/distribute_10__2.loom/softmax_scores.npy'
score2 <- as.numeric(np$load(file2)); length(score2); hist(score2)
file3 <- 'paper_sim/loom/output_distribute_/distribute_10__3.loom/softmax_scores.npy'
score3 <- as.numeric(np$load(file3)); length(score3); hist(score3)
file4 <- 'paper_sim/loom/output_distribute_/distribute_10__4.loom/softmax_scores.npy'
score4 <- as.numeric(np$load(file4)); length(score4); hist(score4)
file5 <- 'paper_sim/loom/output_distribute_/distribute_10__5.loom/softmax_scores.npy'
score5 <- as.numeric(np$load(file5)); length(score5); hist(score5)
file6 <- 'paper_sim/loom/output_distribute_/distribute_10__6.loom/softmax_scores.npy'
score6 <- as.numeric(np$load(file6)); length(score6); hist(score6)
file7 <- 'paper_sim/loom/output_distribute_/distribute_10__7.loom/softmax_scores.npy'
score7 <- as.numeric(np$load(file7)); length(score7); hist(score7)
file8 <- 'paper_sim/loom/output_distribute_/distribute_10__8.loom/softmax_scores.npy'
score8 <- as.numeric(np$load(file8)); length(score8); hist(score8)
file9 <- 'paper_sim/loom/output_distribute_/distribute_10__9.loom/softmax_scores.npy'
score9 <- as.numeric(np$load(file9)); length(score9); hist(score9)
file10 <- 'paper_sim/loom/output_distribute_/distribute_10__10.loom/softmax_scores.npy'
score10 <- as.numeric(np$load(file10)); length(score10); hist(score10)

score <- c(score1, score2, score3, score4, score5, score6, score7, score8, score9, score10)
label <- ifelse(labels8 == 'doublet', 1, 0); table(label)
fg <- score[label == 1]
bg <- score[label == 0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr2 <- pr$auc.integral; pr2
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc

















