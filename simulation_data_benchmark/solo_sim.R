library(reticulate)
library(Matrix)
library(PRROC)
library(RcppCNPy)
np <- import("numpy")
set.seed(2020)

###############################################################################################
# calculate auprc and auroc on simulation datasets under different doublet rates, cell types, 
# sequencing depth, and heterogeneity between cell types
###############################################################################################
# read simulation data under different doublet rates, cell types, 
# sequencing depth, and heterogeneity between cell types
# comment out unused data to choose one out of four
sim.data <- readRDS("paper_sim/sim_rate.rds")
sim.data <- readRDS("paper_sim/sim_type.rds")
sim.data <- readRDS("paper_sim/sim_depth.rds")
sim.data <- readRDS("paper_sim/sim_hetero.rds")
# matrix to save auprc and auroc, col: five simulation datasets per experimental condition
prauc <- matrix(nrow = 21, ncol = 5)
rocauc <- matrix(nrow = 21, ncol = 5)
# loop over each simulation dataset
for(i in 1:21){
  # each experimental condition incldues five simulation datasets
  # calculate on all of them, we use the first one in the paper
  for(j in 1:5){
    # read doublet score already calculated by Solo command
    file <- paste('paper_sim/loom/output_hetero/hetero_', i, '_', j, '.loom/softmax_scores.npy', sep = '')
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
# save the result, change the file name based on which experimental condition we use
saveRDS(list(prauc, rocauc), 'paper_result/solo_hetero.rds')

##############################################################################################################################
# Impact of doublet detection on the identification of highly variable genes (HVGs)
##############################################################################################################################
# simulation data with different doublet rates
sim.data <- readRDS("paper_sim/sim_rate.rds")
# we use three doublet rates: 0.1, 0.2, 0.4
rate <- 0.4
data <- sim.data[[as.character(rate)]]
# use the first simulation dataset out of five
data.matrix <- data[[1]][[1]]; dim(data.matrix)
data.type <- data[[2]][[1]]; table(data.type)
data.type <- ifelse(data.type=='doublet', 1, 0); table(data.type)
index.doublet <- which(data.type==1)
data.clean <- data.matrix[,-index.doublet]; dim(data.clean)
# read calculated doublet scores
file <- 'paper_sim/loom/output_rate/rate_20_1.loom/softmax_scores.npy'
score <- as.numeric(np$load(file)); length(score); hist(score)
# choose threshold 
threshold <- sort(score, decreasing = T)[floor(length(score)*rate)]; threshold
index.pred <- which(score > threshold)
data.pred <- data.matrix[,-index.pred]; dim(data.pred)

# clean variable genes
doublet.seurat <- CreateSeuratObject(counts = data.clean, project = "doublet", min.cells = 0); doublet.seurat
doublet.seurat <- NormalizeData(doublet.seurat)
doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
clean.variable <- VariableFeatures(doublet.seurat)

# contiminated variable genes
doublet.seurat <- CreateSeuratObject(counts = data.matrix, project = "doublet", min.cells = 0); doublet.seurat
doublet.seurat <- NormalizeData(doublet.seurat)
doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
doublet.variable <- VariableFeatures(doublet.seurat)

# post-doublet-detection variable genes
doublet.seurat <- CreateSeuratObject(counts = data.pred, project = "doublet", min.cells = 0); doublet.seurat
doublet.seurat <- NormalizeData(doublet.seurat)
doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
pred.variable <- VariableFeatures(doublet.seurat)

# Jaccard index: negative control
J.clean.doublet <- 
  length(intersect(clean.variable, doublet.variable))/length(union(clean.variable, doublet.variable)); J.clean.doublet
# # Jaccard index: post-doublet-detection
J.clean.pred <- 
  length(intersect(clean.variable, pred.variable))/length(union(clean.variable, pred.variable)); J.clean.pred





