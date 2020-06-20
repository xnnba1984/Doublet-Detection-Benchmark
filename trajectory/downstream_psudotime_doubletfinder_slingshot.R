library(Matrix)
library(Seurat)
library(dplyr)
library(DoubletFinder)
library(PRROC)
source('utility.R')
library(SingleCellExperiment)
library(slingshot, quietly = TRUE)
setwd('/media/nxi/nxi/doublet')
getwd()

system.time(sim.data <- readRDS('sim_data/psudotime_20%_1875_60.rds'))
counts <- sim.data[[1]]; dim(counts)
sim.type <- sim.data[[2]]; table(sim.type)

set.seed(10)
## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
doublet.seurat <- CreateSeuratObject(counts = counts, project = "doublet", min.cells = 1, min.features = 1); doublet.seurat
doublet.seurat <- NormalizeData(doublet.seurat)
doublet.seurat <- ScaleData(doublet.seurat)
doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
doublet.seurat <- RunPCA(doublet.seurat)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.doublet <- paramSweep_v3(doublet.seurat, PCs = 1:10, sct = FALSE, num.cores = 44)
sweep.stats.doublet <- summarizeSweep(sweep.res.doublet, GT = FALSE)
bcmvn.doublet <- find.pK(sweep.stats.doublet)
pK <- bcmvn.doublet$pK[which.max(bcmvn.doublet$BCmetric)]; pK <- as.numeric(levels(pK))[pK]; pK
doublet.seurat <- doubletFinder_v3(doublet.seurat, PCs = 1:10, pN = 0.25, pK = pK, nExp = 60)
attribute <- paste('pANN', 0.25, pK, 60, sep = '_'); attribute
score <- doublet.seurat@meta.data[[attribute]]; score

# prediction result
threshhold <- sort(score, decreasing = TRUE)[60]
pred <- as.numeric(score > threshhold)
pred.index <- which(pred == 1)

###################################################
# downstream
###################################################
counts <- counts[,-pred.index]; dim(counts)
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











































