library(reticulate)
use_python('/home/nxi/anaconda3/bin', required = T)
library(SingleCellExperiment)
library(slingshot, quietly = TRUE)
library(TSCAN)
set.seed(2020)

# read simulation data with single trajectory and temporally expressed genes
data <- readRDS('sim_data/psudotime_tscan.rds')
sim.doublet <- data[[1]]; dim(sim.doublet)
sim.types <- data[[2]]; table(sim.types)
nonDE <- data[[3]]
DE <- data[[4]]
counts <- sim.doublet[,sim.types==0]; dim(counts)

####################################################################################################
# use tscan to identify temporally expressed genes without doublets (postive control)
####################################################################################################
# log transform
procdata <- log2(counts + 1)
lpsmclust <- exprmclust(procdata)
lpsorder <- TSCANorder(lpsmclust)
diffval <- difftest(procdata,lpsorder)

# select temporally expressed genes under qvlue cutoff of 0.05
findDE <- row.names(diffval)[diffval$qval <= 0.05]; length(findDE)
findnonDE <- row.names(diffval)[diffval$qval > 0.05]; length(findnonDE)

# calculate precision, recall, and true negative rate
tp <- length(intersect(findDE, DE)); tp
fp <- length(setdiff(findDE, DE)); fp
fn <- length(setdiff(findnonDE, nonDE)); fn
tn <- length(intersect(findnonDE, nonDE)); tn

precision <- tp/(tp+fp); precision
recall <- tp/(tp+fn); recall
tnr <- tn/(tn+fp); tnr

####################################################################################################
# use tscan to identify temporally expressed genes with doublets (negative control)
####################################################################################################
# log transform
procdata <- log2(sim.doublet + 1)
lpsmclust <- exprmclust(procdata)
lpsorder <- TSCANorder(lpsmclust)
diffval <- difftest(sim.doublet,lpsorder)
# selected temporally expressed genes under qvlue cutoff of 0.05
findDE <- row.names(diffval)[diffval$qval <= 0.05]; length(findDE)
findnonDE <- row.names(diffval)[diffval$qval > 0.05]; length(findnonDE)

# calculate precision, recall, and true negative rate
tp <- length(intersect(findDE, DE)); tp
fp <- length(setdiff(findDE, DE)); fp
fn <- length(setdiff(findnonDE, nonDE)); fn
tn <- length(intersect(findnonDE, nonDE)); tn

precision <- tp/(tp+fp); precision
recall <- tp/(tp+fn); recall
tnr <- tn/(tn+fp); tnr

####################################################################################################
# use tscan to identify temporally expressed genes after removing doublets by scrublet
####################################################################################################
# read python module
scr <- import('scrublet')
scipy.io <- import('scipy.io')
np <- import('numpy')
os <- import('os')

# scrublet
data <- as(sim.doublet, "sparseMatrix"); dim(data)
result <- scr$Scrublet(counts_matrix = t(data), expected_doublet_rate = 0.2, random_state = 10L)
results <- result$scrub_doublets(min_counts=2, 
                                 min_cells=3, 
                                 min_gene_variability_pctl=85, 
                                 n_prin_comps=30L)
# prediction result
score <- as.vector(results[[1]])
threshhold <- sort(score, decreasing = TRUE)[60]
pred <- as.numeric(score > threshhold)
pred.index <- which(pred == 1)

# identify temporally expressed genes
counts <- data[,-pred.index]; dim(counts)
procdata <- log2(counts + 1)
lpsmclust <- exprmclust(as.matrix(procdata))
lpsorder <- TSCANorder(lpsmclust)
diffval <- difftest(sim.doublet,lpsorder)

# selected temporally expressed genes under qvlue cutoff of 0.05
findDE <- row.names(diffval)[diffval$qval <= 0.05]; length(findDE)
findnonDE <- row.names(diffval)[diffval$qval > 0.05]; length(findnonDE)

# calculate precision, recall, and true negative rate
tp <- length(intersect(findDE, DE)); tp
fp <- length(setdiff(findDE, DE)); fp
fn <- length(setdiff(findnonDE, nonDE)); fn
tn <- length(intersect(findnonDE, nonDE)); tn

precision <- tp/(tp+fp); precision
recall <- tp/(tp+fn); recall
tnr <- tn/(tn+fp); tnr

####################################################################################################
# use tscan to identify temporally expressed genes after removing doublets by doubletcells
####################################################################################################
library(scran)
# prediction result
score <- doubletCells(sim.doublet)
# prediction result
threshhold <- sort(score, decreasing = TRUE)[60]
pred <- as.numeric(score > threshhold)
pred.index <- which(pred == 1)

# identify temporally expressed genes
counts <- data[,-pred.index]; dim(counts)
procdata <- log2(counts + 1)
lpsmclust <- exprmclust(as.matrix(procdata))
lpsorder <- TSCANorder(lpsmclust)
diffval <- difftest(sim.doublet,lpsorder)

# selected temporally expressed genes under qvlue cutoff of 0.05
findDE <- row.names(diffval)[diffval$qval <= 0.05]; length(findDE)
findnonDE <- row.names(diffval)[diffval$qval > 0.05]; length(findnonDE)

# calculate precision, recall, and true negative rate
tp <- length(intersect(findDE, DE)); tp
fp <- length(setdiff(findDE, DE)); fp
fn <- length(setdiff(findnonDE, nonDE)); fn
tn <- length(intersect(findnonDE, nonDE)); tn

precision <- tp/(tp+fp); precision
recall <- tp/(tp+fn); recall
tnr <- tn/(tn+fp); tnr

####################################################################################################
# use tscan to identify temporally expressed genes after removing doublets by cxds, bcds, or hybrid
####################################################################################################
library(scds)
# prediction result
sce <- SingleCellExperiment(assays = list(counts = sim.doublet))
sce <- cxds(sce); sce <- bcds(sce); sce <- cxds_bcds_hybrid(sce)
CD <- colData(sce)
# save scores, change the mothod accordingly
score <- CD$cxds_score
threshhold <- sort(score, decreasing = TRUE)[60]
pred <- as.numeric(score > threshhold)
pred.index <- which(pred == 1)

# identify temporally expressed genes
counts <- data[,-pred.index]; dim(counts)
procdata <- log2(counts + 1)
lpsmclust <- exprmclust(as.matrix(procdata))
lpsorder <- TSCANorder(lpsmclust)
diffval <- difftest(sim.doublet,lpsorder)

# selected temporally expressed genes under qvlue cutoff of 0.05
findDE <- row.names(diffval)[diffval$qval <= 0.05]; length(findDE)
findnonDE <- row.names(diffval)[diffval$qval > 0.05]; length(findnonDE)

# calculate precision, recall, and true negative rate
tp <- length(intersect(findDE, DE)); tp
fp <- length(setdiff(findDE, DE)); fp
fn <- length(setdiff(findnonDE, nonDE)); fn
tn <- length(intersect(findnonDE, nonDE)); tn

precision <- tp/(tp+fp); precision
recall <- tp/(tp+fn); recall
tnr <- tn/(tn+fp); tnr

####################################################################################################
# use tscan to identify temporally expressed genes after removing doublets by doubletdetection
####################################################################################################
# read python module
doubletdetection <- import('doubletdetection')
np <- import('numpy')

# doubletdetection
clf <- doubletdetection$BoostClassifier(n_iters=5L, use_phenograph=FALSE, standard_scaling=TRUE, random_state = 10L)
fit <- clf$fit(t(sim.doublet))
p.values <- abs(fit$all_log_p_values_)
score <- abs(apply(p.values, 2, mean))
threshhold <- sort(score, decreasing = TRUE)[60]
pred <- as.numeric(score > threshhold)
pred.index <- which(pred == 1)

# identify temporally expressed genes
counts <- data[,-pred.index]; dim(counts)
procdata <- log2(counts + 1)
lpsmclust <- exprmclust(as.matrix(procdata))
lpsorder <- TSCANorder(lpsmclust)
diffval <- difftest(sim.doublet,lpsorder)

# selected temporally expressed genes under qvlue cutoff of 0.05
findDE <- row.names(diffval)[diffval$qval <= 0.05]; length(findDE)
findnonDE <- row.names(diffval)[diffval$qval > 0.05]; length(findnonDE)

# calculate precision, recall, and true negative rate
tp <- length(intersect(findDE, DE)); tp
fp <- length(setdiff(findDE, DE)); fp
fn <- length(setdiff(findnonDE, nonDE)); fn
tn <- length(intersect(findnonDE, nonDE)); tn

precision <- tp/(tp+fp); precision
recall <- tp/(tp+fn); recall
tnr <- tn/(tn+fp); tnr

####################################################################################################
# use tscan to identify temporally expressed genes after removing doublets by doubletfinder
####################################################################################################
library(DoubletFinder)
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
threshhold <- sort(score, decreasing = TRUE)[60]
pred <- as.numeric(score > threshhold)
pred.index <- which(pred == 1)

# identify temporally expressed genes
counts <- data[,-pred.index]; dim(counts)
procdata <- log2(counts + 1)
lpsmclust <- exprmclust(as.matrix(procdata))
lpsorder <- TSCANorder(lpsmclust)
diffval <- difftest(sim.doublet,lpsorder)

# selected temporally expressed genes under qvlue cutoff of 0.05
findDE <- row.names(diffval)[diffval$qval <= 0.05]; length(findDE)
findnonDE <- row.names(diffval)[diffval$qval > 0.05]; length(findnonDE)

# calculate precision, recall, and true negative rate
tp <- length(intersect(findDE, DE)); tp
fp <- length(setdiff(findDE, DE)); fp
fn <- length(setdiff(findnonDE, nonDE)); fn
tn <- length(intersect(findnonDE, nonDE)); tn

precision <- tp/(tp+fp); precision
recall <- tp/(tp+fn); recall
tnr <- tn/(tn+fp); tnr

####################################################################################################
# use tscan to identify temporally expressed genes after removing doublets by solo
# solo was run by linux command, we read the calculated score directly
####################################################################################################
# read score
file <- 'paper_sim/loom/output_tscan/tscan.loom/softmax_scores.npy'
score <- as.numeric(np$load(file)); length(score); hist(score)
threshhold <- sort(score, decreasing = TRUE)[60]
pred <- as.numeric(score > threshhold)
pred.index <- which(pred == 1)

# identify temporally expressed genes
counts <- data[,-pred.index]; dim(counts)
procdata <- log2(counts + 1)
lpsmclust <- exprmclust(as.matrix(procdata))
lpsorder <- TSCANorder(lpsmclust)
diffval <- difftest(sim.doublet,lpsorder)

# selected temporally expressed genes under qvlue cutoff of 0.05
findDE <- row.names(diffval)[diffval$qval <= 0.05]; length(findDE)
findnonDE <- row.names(diffval)[diffval$qval > 0.05]; length(findnonDE)

# calculate precision, recall, and true negative rate
tp <- length(intersect(findDE, DE)); tp
fp <- length(setdiff(findDE, DE)); fp
fn <- length(setdiff(findnonDE, nonDE)); fn
tn <- length(intersect(findnonDE, nonDE)); tn

precision <- tp/(tp+fp); precision
recall <- tp/(tp+fn); recall
tnr <- tn/(tn+fp); tnr

