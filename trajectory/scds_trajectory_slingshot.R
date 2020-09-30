library(Matrix)
library(Seurat)
library(dplyr)
library(reticulate)
library(PRROC)
library(scds)
source('utility.R')
library(SingleCellExperiment)
library(slingshot, quietly = TRUE)
set.seed(2020)

############################################################################################
# Effects of doublet detection on cell trajectory inference
############################################################################################
# read simulation data with bifurcating trajectories
system.time(sim.data <- readRDS('sim_data/psudotime_500_100.rds'))
# count matrix
counts <- sim.data[[1]]; dim(counts)
# doublet labels
sim.type <- sim.data[[2]]; table(sim.type)

# cxds, bcds, or hybrid to remove doublets
sce <- SingleCellExperiment(assays = list(counts = counts))
sce <- cxds(sce); sce <- bcds(sce); sce <- cxds_bcds_hybrid(sce)
CD <- colData(sce)

# change the scores accordingly
score <- CD$hybrid_score
# select threshold to predict doublets
# identification rate is the true doublet rate
threshhold <- sort(score, decreasing = TRUE)[100]
pred <- as.numeric(score > threshhold)
pred.index <- which(pred == 1)

# lineage construction by slingshot after removing doublets
counts <- data[,-pred.index]; dim(counts)
sim.type <- sim.type[-pred.index]; table(sim.type)
sim <- SingleCellExperiment(assays = List(counts = counts))

# quantile normalization based on the recommendation of slingshot
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min'); dim(rk)
  counts.sort <- apply(counts,2,sort); dim(counts.sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sim)$norm <- FQnorm(assays(sim)$counts)

# pca dimension reduction
pca <- prcomp(t(log1p(assays(sim)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]
reducedDims(sim) <- SimpleList(PCA = rd1)

# GMM clustering
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification; table(cl1)
colData(sim)$GMM <- cl1

# lineage reconstruction
sim <- slingshot(sim, clusterLabels = 'GMM', reducedDim = 'PCA')
pca <- as.data.frame(reducedDims(sim)$PCA)
pca$type <- as.numeric(sim.type)
pca$type <- ifelse(pca$type==0, 'singlet', 'doublet')
pca$type <- as.factor(pca$type)

# plot constructed trajectories, singlets, and doublets
palette(c("red","grey"))
plot(pca$PC1, pca$PC2, col = pca$type, pch=16, asp = 0)
lines(SlingshotDataSet(sim), lwd=2, col='black')

####################################################################################
# temporally expressed genes analysis on single lineage 
####################################################################################
# read simulation data with single trajectory
system.time(sim.data <- readRDS('sim_data/psudotime_20%_1875_60.rds'))
counts <- sim.data[[1]]; dim(counts)
sim.type <- sim.data[[2]]; table(sim.type)
sim <- SingleCellExperiment(assays = List(counts = counts))

# cxds, bcds, or hybrid to remove doublets
sce <- SingleCellExperiment(assays = list(counts = counts))
sce <- cxds(sce); sce <- bcds(sce); sce <- cxds_bcds_hybrid(sce)
CD <- colData(sce)
# change the score accordingly
score <- CD$hybrid_score
threshhold <- sort(score, decreasing = TRUE)[60]
pred <- as.numeric(score > threshhold)
pred.index <- which(pred == 1)

# temporally expressed genes analysis after removing doublets
# remove detected doublets
counts <- data[,-pred.index]; dim(counts)
sim.type <- sim.type[-pred.index]; table(sim.type)
sim <- SingleCellExperiment(assays = List(counts = counts))

# quantile normalization recommended by slingshot
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min'); dim(rk)
  counts.sort <- apply(counts,2,sort); dim(counts.sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sim)$norm <- FQnorm(assays(sim)$counts)

# pca dimension reduction
pca <- prcomp(t(log1p(assays(sim)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]
reducedDims(sim) <- SimpleList(PCA = rd1)

# GMM clustering requested by slingshot
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification; table(cl1)
colData(sim)$GMM <- cl1

# lineage reconstruction
sim <- slingshot(sim, clusterLabels = 'GMM', reducedDim = 'PCA')
sim$slingPseudotime_1[which(sim.type==0)]
cor(1:length(sim$slingPseudotime_1[which(sim.type==0)]), sim$slingPseudotime_1[which(sim.type==0)])

# temporally expressed genes analysis
require(gam)
t <- sim$slingPseudotime_1
Y <- log1p(assays(sim)$norm); dim(Y)

# fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  suppressWarnings({
    tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
  })
  p <- summary(tmp)[3][[1]][2,3]
  p
})
# p-value threshold
p <- .05
# find temporally expressed genes
# calculated precision, recall, true negative rate (TNR)
findDE <- names(gam.pval[gam.pval<=p]); length(findDE)
findnonDE <- names(gam.pval[gam.pval>p]); length(findnonDE)
nonDE <- rownames(counts)[1:500]
DE <- rownames(counts)[501:750]
tp <- length(intersect(findDE, DE)); tp
fp <- length(setdiff(findDE, DE)); fp
fn <- length(setdiff(findnonDE, nonDE)); fn
tn <- length(intersect(findnonDE, nonDE)); tn

precision <- tp/(tp+fp); precision
recall <- tp/(tp+fn); recall
tnr <- tn/(tn+fp); tnr
