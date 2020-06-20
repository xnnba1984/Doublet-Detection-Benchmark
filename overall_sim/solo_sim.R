library(reticulate)
use_python('/home/nxi/anaconda3/bin', required = T)
library(Matrix)
library(PRROC)
library(RcppCNPy)
setwd('/media/nxi/nxi/doublet')
getwd()
np <- import("numpy")

sim.data <- readRDS("paper_sim/sim_hetero.rds")
prauc <- matrix(nrow = 21, ncol = 5)
rocauc <- matrix(nrow = 21, ncol = 5)
for(i in 1:21){
  for(j in 1:5){
    file <- paste('paper_sim/loom/output_hetero/hetero_', i, '_', j, '.loom/softmax_scores.npy', sep = '')
    # data reading
    score <- as.numeric(np$load(file)); length(score); hist(score)
    label <- sim.data[[i]][[2]][[j]]; table(label)
    label <- ifelse(label=='doublet', 1, 0); table(label)
    # auc
    fg <- score[label==1]
    bg <- score[label==0]
    pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T); pr$auc.integral
    roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc
    prauc[i, j] <- pr$auc.integral
    rocauc[i, j] <- roc$auc
  }
}
saveRDS(list(prauc, rocauc), 'paper_result/solo_hetero.rds')
pr <- apply(prauc, 1, mean); pr
plot(pr); lines(pr)
roc <- apply(rocauc, 1, mean); roc
plot(roc); lines(roc)






