library(reticulate)
library(scater)
library(SingleCellExperiment)
library(slingshot, quietly = TRUE)
set.seed(2020)

#################################################################################
# plot cell trajectory by MST with 20% doublet (negative control)
#################################################################################
# read simulation data with 3 sequential trajectories
system.time(sim.data <- readRDS('sim_data/psudotime_sequential.rds'))
# count matrix
sim.doublet <- sim.data[[1]]; dim(sim.doublet)
# doublet labels
sim.types <- sim.data[[2]]; table(sim.types)
doublet.o <- SingleCellExperiment(assays = List(counts = sim.doublet)); dim(doublet.o)
doublet.o <- normalize(doublet.o)

# quantile normalization
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(doublet.o)$norm <- FQnorm(assays(doublet.o)$counts)

# pca dimension reduction
pca <- prcomp(t(log1p(assays(doublet.o)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]
reducedDims(doublet.o) <- SimpleList(PCA = rd1)

# GMM clustering
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification
colData(doublet.o)$GMM <- cl1
doublet.o <- slingshot(doublet.o, clusterLabels = 'GMM', reducedDim = 'PCA')
pca.doublets <- as.data.frame(reducedDims(doublet.o)$PCA)
pca.doublets$type <- as.numeric(sim.types)
pca.doublets$type <- ifelse(pca.doublets$type==0, 'singlet', 'doublet')
pca.doublets$type <- as.factor(pca.doublets$type)
palette(c("red","grey"))

# dataset without doublets
pca.singlet <- pca.doublets[pca.doublets$type=='singlet',]; dim(pca.singlet)

# plot inferred trajectories with doublets
plot(pca.doublets$PC1, pca.doublets$PC2, col = pca.doublets$type, pch=16, asp = 0, 
     main='Original Doublets', xlim=c(min(pca.doublets$PC1),max(pca.doublets$PC1)),
     xlab=NA, ylab=NA,xaxt='n', yaxt='n', cex.main=2, cex=1.3,font.main = 1)
lines(SlingshotDataSet(doublet.o), lwd=7, type = 'lineages', col = 'black')

# plot inferred trajectories without doublets (positive control)
plot(pca.singlet$PC1, pca.singlet$PC2, col = 'grey', pch=16, asp = 0, main='No Doublets',
     xlim=c(min(pca.singlet$PC1), max(pca.singlet$PC1)),
     xlab=NA, ylab=NA,xaxt='n', yaxt='n', cex.main=2, cex=1.3,font.main = 1)
lines(SlingshotDataSet(sim), lwd=7, type = 'lineages', col = 'black')

##################################################################################################################
# apply scrublet to remove doublet and contruct the trajectories
##################################################################################################################
py_config()
library(Matrix)
library(Seurat)
library(dplyr)
library(PRROC)
source('utility.R')

# import python modules
scr <- import('scrublet')
scipy.io <- import('scipy.io')
np <- import('numpy')
os <- import('os')

# scrublet
data <- as(sim.doublet, "sparseMatrix"); dim(data)
result <- scr$Scrublet(counts_matrix = t(data), expected_doublet_rate = .2, random_state = 10L)
results <- result$scrub_doublets(min_counts=2, 
                                 min_cells=3, 
                                 min_gene_variability_pctl=85, 
                                 n_prin_comps=30L)
# prediction result
score <- as.vector(results[[1]])
threshhold <- sort(score, decreasing = TRUE)[250]
pred <- as.numeric(score > threshhold); table(pred)
pred.index <- which(pred == 1)

# redo lineage on the after-removing-doublets dataset
counts <- data[,-pred.index]; dim(counts)
type <- sim.types[-pred.index]; table(type)
sim.scrublet <- SingleCellExperiment(assays = List(counts = counts))

# quantile normalization
assays(sim.scrublet)$norm <- FQnorm(assays(sim.scrublet)$counts)

# pca
pca.scrublet <- prcomp(t(log1p(assays(sim.scrublet)$norm)), scale. = FALSE)
rd1 <- pca.scrublet$x[,1:2]; dim(rd1)
reducedDims(sim.scrublet) <- SimpleList(PCA = rd1)

# GMM clustering
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification; table(cl1)
colData(sim.scrublet)$GMM <- cl1

# lineage reconstruction
sim.scrublet <- slingshot(sim.scrublet, clusterLabels = 'GMM', reducedDim = 'PCA')
pca.scrublet <- as.data.frame(reducedDims(sim.scrublet)$PCA)
pca.scrublet$type <- as.numeric(type)
pca.scrublet$type <- ifelse(pca.scrublet$type==0, 'singlet', 'doublet')
pca.scrublet$type <- as.factor(pca.scrublet$type)
palette(c("red","grey"))

# plot trajectory after scrublet
pca.scrublet <- pca.doublets[-pred.index,]; dim(pca.scrublet)
palette(c("red","grey"))
plot(pca.scrublet$PC1, pca.scrublet$PC2, col = pca.scrublet$type, pch=16, asp = 0, main='Scrublet',
     xlab=NA, ylab=NA,xaxt='n', yaxt='n', cex.main=2, cex=1.3,font.main = 1)
lines(SlingshotDataSet(sim.scrublet), lwd=7, type = 'lineages', col = 'black')

##################################################################################################################
# apply doubletcells to remove doublet and contruct the trajectories
##################################################################################################################
library(scran)
score <- doubletCells(data)
# prediction result
threshhold <- sort(score, decreasing = TRUE)[250]
pred <- as.numeric(score > threshhold)
pred.index <- which(pred == 1)

# redo lineage on after-doublet-removing datasets
counts <- data[,-pred.index]; dim(counts)
type <- sim.types[-pred.index]; table(type)
sim.dblcell <- SingleCellExperiment(assays = List(counts = counts))

# quantile normalization
assays(sim.dblcell)$norm <- FQnorm(assays(sim.dblcell)$counts)

# pca
pca.dblcell <- prcomp(t(log1p(assays(sim.dblcell)$norm)), scale. = FALSE)
rd1 <- pca.dblcell$x[,1:2]
reducedDims(sim.dblcell) <- SimpleList(PCA = rd1)

# GMM clustering
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification; table(cl1)
colData(sim.dblcell)$GMM <- cl1

# lineage reconstruction
sim.dblcell <- slingshot(sim.dblcell, clusterLabels = 'GMM', reducedDim = 'PCA')
pca.dblcell <- as.data.frame(reducedDims(sim.dblcell)$PCA)
pca.dblcell$type <- as.numeric(type)
pca.dblcell$type <- ifelse(pca.dblcell$type==0, 'singlet', 'doublet')
pca.dblcell$type <- as.factor(pca.dblcell$type)
palette(c("red","grey"))

# plot trajectory after applying doubletcells
pca.dblcell <- pca.doublets[-pred.index,]; dim(pca.dblcell)
palette(c("red","grey"))
plot(pca.dblcell$PC1, pca.dblcell$PC2, col = pca.dblcell$type, pch=16, asp = 0, main='doubletCells',
     xlab=NA, ylab=NA,xaxt='n', yaxt='n', cex.main=2, cex=1.3,font.main = 1)
line.dblcell <- SlingshotDataSet(sim.dblcell)
lines(SlingshotDataSet(sim.dblcell), lwd=7, type = 'lineages', col = 'black')

##################################################################################################################
# apply cxds, bcds, or hybrid to remove doublet and contruct the trajectories
##################################################################################################################
library(scds)
sce <- SingleCellExperiment(assays = list(counts = data))
sce <- cxds(sce); sce <- bcds(sce); sce <- cxds_bcds_hybrid(sce)
CD <- colData(sce)
# change scores based on different methods
score <- CD$cxds_score
# prediction result
threshhold <- sort(score, decreasing = TRUE)[250]
pred <- as.numeric(score > threshhold); table(pred)
pred.index <- which(pred == 1)

# redo lineage after doublet detection
counts <- data[,-pred.index]; dim(counts)
type <- sim.types[-pred.index]; table(type)
sim.cxds <- SingleCellExperiment(assays = List(counts = counts))

# quantile normalization
assays(sim.cxds)$norm <- FQnorm(assays(sim.cxds)$counts)

# pca
pca.cxds <- prcomp(t(log1p(assays(sim.cxds)$norm)), scale. = FALSE)
rd1 <- pca.cxds$x[,1:2]
reducedDims(sim.cxds) <- SimpleList(PCA = rd1)

# GMM clustering
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification; table(cl1)
colData(sim.cxds)$GMM <- cl1

# lineage reconstruction
sim.cxds <- slingshot(sim.cxds, clusterLabels = 'GMM', reducedDim = 'PCA')
pca.cxds <- as.data.frame(reducedDims(sim.cxds)$PCA)
pca.cxds$type <- as.numeric(type)
pca.cxds$type <- ifelse(pca.cxds$type==0, 'singlet', 'doublet')
pca.cxds$type <- as.factor(pca.cxds$type)

# plot trajectory after cxds, bcds, or hybrid
pca.cxds <- pca.doublets[-pred.index,]; dim(pca.cxds)
palette(c("red","grey"))
plot(pca.cxds$PC1, pca.cxds$PC2, col = pca.cxds$type, pch=16, asp = 0, main='cxds', 
     xlab=NA, ylab=NA,xaxt='n', yaxt='n', cex.main=2, cex=1.3,font.main = 1)
lines(SlingshotDataSet(sim.cxds), lwd=7, type = 'lineages', col = 'black')

##################################################################################################################
# apply cxds, bcds, or hybrid to remove doublet and contruct the trajectories
##################################################################################################################
# read python module
doubletdetection <- import('doubletdetection')
np <- import('numpy')
# use doubletdetection to remove doublets
drop <- which(apply(counts, 1, sum)==0)
counts <- data[-drop,];dim(data)
clf <- doubletdetection$BoostClassifier(n_iters=5L, use_phenograph=FALSE, standard_scaling=TRUE, random_state = 10L)
fit <- clf$fit(t(counts))
p.values <- abs(fit$all_log_p_values_)
score <- abs(apply(p.values, 2, mean))

# prediction result
threshhold <- sort(score, decreasing = TRUE)[250]
pred <- as.numeric(score > threshhold)
pred.index <- which(pred == 1)

# construct lineage after applying doubletdetection
counts <- data[,-pred.index]; dim(counts)
type <- sim.types[-pred.index]; table(type)
sim.doubletdetection <- SingleCellExperiment(assays = List(counts = counts))

# quantile normalization
assays(sim.doubletdetection)$norm <- FQnorm(assays(sim.doubletdetection)$counts)

# pca
pca.doubletdetection <- prcomp(t(log1p(assays(sim.doubletdetection)$norm)), scale. = FALSE)
rd1 <- pca.doubletdetection$x[,1:2]
reducedDims(sim.doubletdetection) <- SimpleList(PCA = rd1)

# GMM clustering
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification; table(cl1)
colData(sim.doubletdetection)$GMM <- cl1

# lineage reconstruction
sim.doubletdetection <- slingshot(sim.doubletdetection, clusterLabels = 'GMM', reducedDim = 'PCA')
pca.doubletdetection <- as.data.frame(reducedDims(sim.doubletdetection)$PCA)
pca.doubletdetection$type <- as.numeric(type)
pca.doubletdetection$type <- ifelse(pca.doubletdetection$type==0, 'singlet', 'doublet')
pca.doubletdetection$type <- as.factor(pca.doubletdetection$type)

# plot trajectory after doubletdetection
pca.doubletdetection <- pca.doublets[-pred.index,]; dim(pca.doubletdetection)
palette(c("red","grey"))
plot(pca.doubletdetection$PC1, pca.doubletdetection$PC2, col = pca.doubletdetection$type, pch=16, asp = 0, 
     main='DoubletDetection', xlab=NA, ylab=NA,xaxt='n', yaxt='n', cex.main=2, cex=1.3,font.main = 1)
lines(SlingshotDataSet(sim), lwd=7, type = 'lineages', col = 'black')

##################################################################################################################
# apply doubletfinder to remove doublet and contruct the trajectories
##################################################################################################################
library(DoubletFinder)
## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
doublet.seurat <- CreateSeuratObject(counts = data, project = "doublet", min.cells = 1, min.features = 1); doublet.seurat
doublet.seurat <- NormalizeData(doublet.seurat)
doublet.seurat <- ScaleData(doublet.seurat)
doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
doublet.seurat <- RunPCA(doublet.seurat)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.doublet <- paramSweep_v3(doublet.seurat, PCs = 1:10, sct = FALSE, num.cores = 44)
sweep.stats.doublet <- summarizeSweep(sweep.res.doublet, GT = FALSE)
bcmvn.doublet <- find.pK(sweep.stats.doublet)
pK <- bcmvn.doublet$pK[which.max(bcmvn.doublet$BCmetric)]; pK <- as.numeric(levels(pK))[pK]; pK
doublet.seurat <- doubletFinder_v3(doublet.seurat, PCs = 1:10, pN = 0.25, pK = pK, nExp = 250)
attribute <- paste('pANN', 0.25, pK, 250, sep = '_'); attribute
score <- doublet.seurat@meta.data[[attribute]]; score

# prediction result
threshhold <- sort(score, decreasing = TRUE)[250]
pred <- as.numeric(score > threshhold); table(pred)
pred.index <- which(pred == 1)

# redo lineage after removing doublets by doubletfinder
counts <- data[,-pred.index]; dim(counts)
type <- sim.types[-pred.index]; table(type)
sim.doubletfinder <- SingleCellExperiment(assays = List(counts = counts))

# quantile normalization
assays(sim.doubletfinder)$norm <- FQnorm(assays(sim.doubletfinder)$counts)

# pca
pca.doubletfinder <- prcomp(t(log1p(assays(sim.doubletfinder)$norm)), scale. = FALSE)
rd1 <- pca.doubletfinder$x[,1:2]
reducedDims(sim.doubletfinder) <- SimpleList(PCA = rd1)

# GMM clustering
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification; table(cl1)
colData(sim.doubletfinder)$GMM <- cl1

# lineage reconstruction
sim.doubletfinder <- slingshot(sim.doubletfinder, clusterLabels = 'GMM', reducedDim = 'PCA')
pca.doubletfinder <- as.data.frame(reducedDims(sim.doubletfinder)$PCA)
pca.doubletfinder$type <- as.numeric(type)
pca.doubletfinder$type <- ifelse(pca.doubletfinder$type==0, 'singlet', 'doublet')
pca.doubletfinder$type <- as.factor(pca.doubletfinder$type)

# plot trajectory after doubletfinder
pca.doubletfinder <- pca.doublets[-pred.index,]; dim(pca.doubletfinder)
palette(c("red","grey"))
plot(pca.doubletfinder$PC1, pca.doubletfinder$PC2, col = pca.doubletfinder$type, pch=16, asp = 0,
     main='DoubletFinder', xlab=NA, ylab=NA,xaxt='n', yaxt='n', cex.main=2, cex=1.3,font.main = 1)
lines(SlingshotDataSet(sim.doubletfinder), lwd=7, type = 'lineages', col = 'black')

##################################################################################################################
# apply solo to remove doublet and contruct the trajectories
# solo was running by commands
##################################################################################################################
np <- import("numpy")
file <- 'paper_sim/loom/output_psudotime_mst/psudotime_mst.loom/softmax_scores.npy'
score <- as.numeric(np$load(file)); length(score); hist(score)
threshhold <- sort(score, decreasing = TRUE)[250]
pred <- as.numeric(score > threshhold); table(pred)
pred.index <- which(pred == 1)

# redo lineage after 
counts <- data[,-pred.index]; dim(counts)
type <- sim.types[-pred.index]; table(type)
sim.scrublet <- SingleCellExperiment(assays = List(counts = counts))

# quantile normalization
assays(sim.scrublet)$norm <- FQnorm(assays(sim.scrublet)$counts)

# pca
pca.scrublet <- prcomp(t(log1p(assays(sim.scrublet)$norm)), scale. = FALSE)
rd1 <- pca.scrublet$x[,1:2]; dim(rd1)
reducedDims(sim.scrublet) <- SimpleList(PCA = rd1)

# GMM clustering
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification; table(cl1)
colData(sim.scrublet)$GMM <- cl1

# lineage reconstruction
sim.scrublet <- slingshot(sim.scrublet, clusterLabels = 'GMM', reducedDim = 'PCA')
pca.scrublet <- as.data.frame(reducedDims(sim.scrublet)$PCA)
pca.scrublet$type <- as.numeric(type)
pca.scrublet$type <- ifelse(pca.scrublet$type==0, 'singlet', 'doublet')
pca.scrublet$type <- as.factor(pca.scrublet$type)

# plot trajectory after solo
pca.scrublet <- pca.doublets[-pred.index,]; dim(pca.scrublet)
palette(c("red","grey"))
plot(pca.scrublet$PC1, pca.scrublet$PC2, col = pca.scrublet$type, pch=16, asp = 0, main='solo',
     xlab=NA, ylab=NA,xaxt='n', yaxt='n', cex.main=2, cex=1.3,font.main = 1)
lines(SlingshotDataSet(sim.scrublet), lwd=7, type = 'lineages', col = 'black')
