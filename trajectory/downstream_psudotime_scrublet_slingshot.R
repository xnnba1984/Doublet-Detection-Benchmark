library(reticulate)
use_python('/home/nxi/anaconda3/bin', required = T)
py_config()
library(Matrix)
library(Seurat)
library(dplyr)
library(PRROC)
source('utility.R')
library(SingleCellExperiment)
library(slingshot, quietly = TRUE)
setwd('/media/nxi/nxi/doublet')
getwd()

scr <- import('scrublet')
scipy.io <- import('scipy.io')
np <- import('numpy')
os <- import('os')

##############################################
# two branches
##############################################
system.time(sim.data <- readRDS('sim_data/psudotime_500_100.rds'))
counts <- sim.data[[1]]; dim(counts)
sim.type <- sim.data[[2]]; table(sim.type)

data <- as(counts, "sparseMatrix"); dim(data)
result <- scr$Scrublet(counts_matrix = t(data), expected_doublet_rate = 0.17, random_state = 1L)
results <- result$scrub_doublets(min_counts=2, 
                                 min_cells=3, 
                                 min_gene_variability_pctl=85, 
                                 n_prin_comps=30L)
# prediction result
score <- as.vector(results[[1]])
threshhold <- sort(score, decreasing = TRUE)[100]
pred <- as.numeric(score > threshhold)
pred.index <- which(pred == 1)

##########################################
# redo lineage
##########################################
counts <- data[,-pred.index]; dim(counts)
sim.type <- sim.type[-pred.index]; table(sim.type)
sim <- SingleCellExperiment(assays = List(counts = counts))
# quantile normalization
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min'); dim(rk)
  counts.sort <- apply(counts,2,sort); dim(counts.sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sim)$norm <- FQnorm(assays(sim)$counts)
# pca
pca <- prcomp(t(log1p(assays(sim)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]
plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 0)
reducedDims(sim) <- SimpleList(PCA = rd1)

# GMM clustering
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification; table(cl1)
colData(sim)$GMM <- cl1
library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 0)

# lineage reconstruction
sim <- slingshot(sim, clusterLabels = 'GMM', reducedDim = 'PCA')
pca <- as.data.frame(reducedDims(sim)$PCA)
pca$type <- as.numeric(sim.type)
pca$type <- ifelse(pca$type==0, 'singlet', 'doublet')
pca$type <- as.factor(pca$type)
palette(c("red","grey"))
#png(filename="figure/psudaotime_scrublet_clean.png")
plot(pca$PC1, pca$PC2, col = pca$type, pch=16, asp = 0)
lines(SlingshotDataSet(sim), lwd=2, col='black')
#dev.off()

############################################
# plot in original space
############################################
counts <- sim.data[[1]]; dim(counts)
sim.type <- sim.data[[2]]; table(sim.type)

sim.original <- SingleCellExperiment(assays = List(counts = counts)); dim(sim.original)
# quantile normalization
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min'); dim(rk)
  counts.sort <- apply(counts,2,sort); dim(counts.sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sim.original)$norm <- FQnorm(assays(sim.original)$counts)
# pca
pca <- prcomp(t(log1p(assays(sim.original)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]
plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 0)
reducedDims(sim.original) <- SimpleList(PCA = rd1)

# GMM clustering
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification; table(cl1)
colData(sim.original)$GMM <- cl1
library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 0)

# lineage reconstruction
sim.original <- slingshot(sim.original, clusterLabels = 'GMM', reducedDim = 'PCA')
pca <- as.data.frame(reducedDims(sim.original)$PCA)
pca$type <- as.numeric(sim.type)
pca$type <- ifelse(pca$type==0, 'singlet', 'doublet')
pca$type <- as.factor(pca$type)
pca <- pca[-pred.index,]; dim(pca)
palette(c("red","grey"))
#png(filename="figure/psudaotime_scrublet_clean.png")
pdf('paper_figure/psudotime_scrublet_splingshot.pdf')
plot(pca$PC1, pca$PC2, col = pca$type, pch=16, asp = 0, main='Scrublet', xlab='PC1', ylab='PC2')
lines(SlingshotDataSet(sim), lwd=2, col='black')
dev.off()


##########################################
# single lineage DE
##########################################
system.time(sim.data <- readRDS('sim_data/psudotime_20%_1875_60.rds'))
counts <- sim.data[[1]]; dim(counts)
sim.type <- sim.data[[2]]; table(sim.type)

set.seed(10)
data <- as(counts, "sparseMatrix"); dim(data)
result <- scr$Scrublet(counts_matrix = t(data), expected_doublet_rate = 0.2, random_state = 1L)
results <- result$scrub_doublets(min_counts=2, 
                                 min_cells=3, 
                                 min_gene_variability_pctl=85, 
                                 n_prin_comps=30L)
# prediction result
score <- as.vector(results[[1]])
threshhold <- sort(score, decreasing = TRUE)[60]
pred <- as.numeric(score > threshhold)
pred.index <- which(pred == 1)

###################################################
# downstream
###################################################
counts <- data[,-pred.index]; dim(counts)
sim.type <- sim.type[-pred.index]; table(sim.type)
sim <- SingleCellExperiment(assays = List(counts = counts))
# quantile normalization
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min'); dim(rk)
  counts.sort <- apply(counts,2,sort); dim(counts.sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sim)$norm <- FQnorm(assays(sim)$counts)

# pca
pca <- prcomp(t(log1p(assays(sim)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]
plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 0)
reducedDims(sim) <- SimpleList(PCA = rd1)

# GMM clustering
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification; table(cl1)
colData(sim)$GMM <- cl1
library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 0)

# lineage reconstruction
sim <- slingshot(sim, clusterLabels = 'GMM', reducedDim = 'PCA')
sim$slingPseudotime_1[which(sim.type==0)]
cor(1:length(sim$slingPseudotime_1[which(sim.type==0)]), sim$slingPseudotime_1[which(sim.type==0)])

# plot trajectory
#colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
#plotcol <- colors[cut(sim$slingPseudotime_1, breaks=100)]
#png(filename="figure/psudaotime_nodoublet.png")
plot(reducedDims(sim)$PCA, col = 'gray', pch=16, asp = 0)
lines(SlingshotDataSet(sim), lwd=2, col='black')
#dev.off()

#plot(reducedDims(sim)$PCA, pch=16, asp = 0)


# MST
#plot(reducedDims(sim)$PCA, col = brewer.pal(9,'Set1')[sim$GMM], pch=16, asp = 0)
#lines(SlingshotDataSet(sim), lwd=2, type = 'lineages', col = 'black')

# DE gene
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
p <- .05
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
f1 <- 2 * precision * recall / (precision + recall); f1











































