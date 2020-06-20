library(reticulate)
use_python('/home/nxi/anaconda3/bin', required = T)
library(SingleCellExperiment)
library(slingshot, quietly = TRUE)
library(TSCAN)
seed <- 10

###########################################
# add doublets
###########################################
set.seed(seed)
# generate synthetic count data representing a single lineage
means <- rbind(
  # non-DE genes 500
  matrix(rep(rep(c(0.1,0.5,1,2,3), each = 300),100),
         ncol = 300, byrow = TRUE),
  # early deactivation 50
  matrix(rep(exp(atan( ((300:1)-200)/50 )),50), ncol = 300, byrow = TRUE),
  # late deactivation 50
  matrix(rep(exp(atan( ((300:1)-100)/50 )),50), ncol = 300, byrow = TRUE),
  # early activation 50
  matrix(rep(exp(atan( ((1:300)-100)/50 )),50), ncol = 300, byrow = TRUE),
  # late activation 50
  matrix(rep(exp(atan( ((1:300)-200)/50 )),50), ncol = 300, byrow = TRUE),
  # transient 50
  matrix(rep(exp(atan( c((1:100)/33, rep(3,100), (100:1)/33) )),50), 
         ncol = 300, byrow = TRUE)
); dim(means)

counts <- apply(means,2,function(cell_means){
  total <- rnbinom(1, mu = round(7500/4), size = 4)
  rmultinom(1, total, cell_means)
}); dim(counts)
rownames(counts) <- paste0('G',1:750)
colnames(counts) <- paste0('c',1:300)
nonDE <- rownames(counts)[1:500]
DE <- rownames(counts)[501:750]

##################################################################
# DE without doublets
##################################################################
procdata <- log2(counts + 1)
lpsmclust <- exprmclust(procdata)
lpsorder <- TSCANorder(lpsmclust)
diffval <- difftest(procdata,lpsorder)
#Selected differentially expressed genes under qvlue cutoff of 0.05
findDE <- row.names(diffval)[diffval$qval <= 0.05]; length(findDE)
findnonDE <- row.names(diffval)[diffval$qval > 0.05]; length(findnonDE)
nonDE <- rownames(counts)[1:500]
DE <- rownames(counts)[501:750]
tp <- length(intersect(findDE, DE)); tp
fp <- length(setdiff(findDE, DE)); fp
fn <- length(setdiff(findnonDE, nonDE)); fn
tn <- length(intersect(findnonDE, nonDE)); tn

precision <- tp/(tp+fp); precision
recall <- tp/(tp+fn); recall
tnr <- tn/(tn+fp); tnr

#############################################################
# add doublets
#############################################################
doublet.number <- 60
# cell merage ratio: for each cell c = d*(a*c1 + b*c2), a+1=2, a~N(1, 0.1), d~N(1, 0.1)
a <- rnorm(doublet.number, 1, 0.1)
b <- 2 - a
c <- rnorm(doublet.number, 1, 0.1)

# sample doublets
doublet.pair <- combn(ncol(counts), 2); dim(doublet.pair)
doublet.pair <- doublet.pair[,sample(ncol(doublet.pair), doublet.number)]; dim(doublet.pair)
doublets <- sapply(1:ncol(doublet.pair), function(pair) {
  pair.1 <- doublet.pair[,pair][1]
  pair.2 <- doublet.pair[,pair][2]
  doublet.cell <- round(c[pair]*(a[pair]*counts[,pair.1] + b[pair]*counts[,pair.2]))
  return(doublet.cell)
}); class(doublets); dim(doublets)

# drop doublets ingradients; combine doublets and non-doublets
#drop.index <- (unique(as.vector(doublet.pair))); length(drop.index)
#keep.index <- (1:300)[-drop.index]; length(keep.index)
#length.singlet <- ncol(counts)-length(drop.index); length.singlet
#sim.doublet <- counts[,-drop.index]; dim(sim.doublet)
sim.doublet <- cbind(counts, doublets); dim(sim.doublet)
colnames(sim.doublet) <- as.character(1:ncol(sim.doublet))
#sim.types <- c(rep(0, length.singlet), rep(1, doublet.number)); table(sim.types)

##################################################################
# DE with doublets
##################################################################
procdata <- log2(sim.doublet + 1)
lpsmclust <- exprmclust(procdata)
lpsorder <- TSCANorder(lpsmclust)
diffval <- difftest(sim.doublet,lpsorder)
#Selected differentially expressed genes under qvlue cutoff of 0.05
findDE <- row.names(diffval)[diffval$qval <= 0.05]; length(findDE)
findnonDE <- row.names(diffval)[diffval$qval > 0.05]; length(findnonDE)
nonDE <- rownames(counts)[1:500]
DE <- rownames(counts)[501:750]
tp <- length(intersect(findDE, DE)); tp
fp <- length(setdiff(findDE, DE)); fp
fn <- length(setdiff(findnonDE, nonDE)); fn
tn <- length(intersect(findnonDE, nonDE)); tn

precision <- tp/(tp+fp); precision
recall <- tp/(tp+fn); recall
tnr <- tn/(tn+fp); tnr

##########################################
# scrublet DE
##########################################
set.seed(10)
scr <- import('scrublet')
scipy.io <- import('scipy.io')
np <- import('numpy')
os <- import('os')
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

###################################################
# downstream
###################################################
counts <- data[,-pred.index]; dim(counts)
procdata <- log2(counts + 1)
lpsmclust <- exprmclust(as.matrix(procdata))
lpsorder <- TSCANorder(lpsmclust)
diffval <- difftest(sim.doublet,lpsorder)
#Selected differentially expressed genes under qvlue cutoff of 0.05
findDE <- row.names(diffval)[diffval$qval <= 0.05]; length(findDE)
findnonDE <- row.names(diffval)[diffval$qval > 0.05]; length(findnonDE)

tp <- length(intersect(findDE, DE)); tp
fp <- length(setdiff(findDE, DE)); fp
fn <- length(setdiff(findnonDE, nonDE)); fn
tn <- length(intersect(findnonDE, nonDE)); tn

precision <- tp/(tp+fp); precision
recall <- tp/(tp+fn); recall
tnr <- tn/(tn+fp); tnr

##########################################
# dblcells DE
##########################################
library(scran)
set.seed(10)
# prediction result
score <- doubletCells(sim.doublet)
threshhold <- sort(score, decreasing = TRUE)[60]
pred <- as.numeric(score > threshhold)
pred.index <- which(pred == 1)

###################################################
# downstream
###################################################
counts <- sim.doublet[,-pred.index]; dim(counts)
procdata <- log2(counts + 1)
lpsmclust <- exprmclust(as.matrix(procdata))
lpsorder <- TSCANorder(lpsmclust)
diffval <- difftest(sim.doublet,lpsorder)
#Selected differentially expressed genes under qvlue cutoff of 0.05
findDE <- row.names(diffval)[diffval$qval <= 0.05]; length(findDE)
findnonDE <- row.names(diffval)[diffval$qval > 0.05]; length(findnonDE)

tp <- length(intersect(findDE, DE)); tp
fp <- length(setdiff(findDE, DE)); fp
fn <- length(setdiff(findnonDE, nonDE)); fn
tn <- length(intersect(findnonDE, nonDE)); tn

precision <- tp/(tp+fp); precision
recall <- tp/(tp+fn); recall
tnr <- tn/(tn+fp); tnr

##########################################
# scds DE
##########################################
library(scds)
set.seed(10)
# prediction result
sce <- SingleCellExperiment(assays = list(counts = sim.doublet))
sce <- cxds(sce); sce <- bcds(sce); sce <- cxds_bcds_hybrid(sce)
CD <- colData(sce)
# save scores
score <- CD$cxds_score
threshhold <- sort(score, decreasing = TRUE)[60]
pred <- as.numeric(score > threshhold)
pred.index <- which(pred == 1)

###################################################
# downstream
###################################################
counts <- sim.doublet[,-pred.index]; dim(counts)
procdata <- log2(counts + 1)
lpsmclust <- exprmclust(as.matrix(procdata))
lpsorder <- TSCANorder(lpsmclust)
diffval <- difftest(sim.doublet,lpsorder)
#Selected differentially expressed genes under qvlue cutoff of 0.05
findDE <- row.names(diffval)[diffval$qval <= 0.05]; length(findDE)
findnonDE <- row.names(diffval)[diffval$qval > 0.05]; length(findnonDE)

tp <- length(intersect(findDE, DE)); tp
fp <- length(setdiff(findDE, DE)); fp
fn <- length(setdiff(findnonDE, nonDE)); fn
tn <- length(intersect(findnonDE, nonDE)); tn

precision <- tp/(tp+fp); precision
recall <- tp/(tp+fn); recall
tnr <- tn/(tn+fp); tnr

##########################################
# doubletdetection DE
##########################################
doubletdetection <- import('doubletdetection')
np <- import('numpy')
set.seed(10)
clf <- doubletdetection$BoostClassifier(n_iters=5L, use_phenograph=FALSE, standard_scaling=TRUE, random_state = 10L)
fit <- clf$fit(t(sim.doublet))
p.values <- abs(fit$all_log_p_values_)
score <- abs(apply(p.values, 2, mean))
threshhold <- sort(score, decreasing = TRUE)[60]
pred <- as.numeric(score > threshhold)
pred.index <- which(pred == 1)

###################################################
# downstream
###################################################
counts <- sim.doublet[,-pred.index]; dim(counts)
procdata <- log2(counts + 1)
lpsmclust <- exprmclust(as.matrix(procdata))
lpsorder <- TSCANorder(lpsmclust)
diffval <- difftest(sim.doublet,lpsorder)
#Selected differentially expressed genes under qvlue cutoff of 0.05
findDE <- row.names(diffval)[diffval$qval <= 0.05]; length(findDE)
findnonDE <- row.names(diffval)[diffval$qval > 0.05]; length(findnonDE)
nonDE <- rownames(counts)[1:500]
DE <- rownames(counts)[501:750]
tp <- length(intersect(findDE, DE)); tp
fp <- length(setdiff(findDE, DE)); fp
fn <- length(setdiff(findnonDE, nonDE)); fn
tn <- length(intersect(findnonDE, nonDE)); tn

precision <- tp/(tp+fp); precision
recall <- tp/(tp+fn); recall
tnr <- tn/(tn+fp); tnr

##########################################
# doubletfinder DE
##########################################
library(DoubletFinder)
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
threshhold <- sort(score, decreasing = TRUE)[60]
pred <- as.numeric(score > threshhold)
pred.index <- which(pred == 1)

###################################################
# downstream
###################################################
counts <- sim.doublet[,-pred.index]; dim(counts)
procdata <- log2(counts + 1)
lpsmclust <- exprmclust(as.matrix(procdata))
lpsorder <- TSCANorder(lpsmclust)
diffval <- difftest(sim.doublet,lpsorder)
#Selected differentially expressed genes under qvlue cutoff of 0.05
findDE <- row.names(diffval)[diffval$qval <= 0.05]; length(findDE)
findnonDE <- row.names(diffval)[diffval$qval > 0.05]; length(findnonDE)
nonDE <- rownames(counts)[1:500]
DE <- rownames(counts)[501:750]
tp <- length(intersect(findDE, DE)); tp
fp <- length(setdiff(findDE, DE)); fp
fn <- length(setdiff(findnonDE, nonDE)); fn
tn <- length(intersect(findnonDE, nonDE)); tn

precision <- tp/(tp+fp); precision
recall <- tp/(tp+fn); recall
tnr <- tn/(tn+fp); tnr

###################################################
# loom
###################################################
library(loomR)
create('paper_sim/loom/tscan.loom', sim.doublet, transpose = T)

file <- 'paper_sim/loom/output_tscan/tscan.loom/softmax_scores.npy'
score <- as.numeric(np$load(file)); length(score); hist(score)
threshhold <- sort(score, decreasing = TRUE)[60]
pred <- as.numeric(score > threshhold)
pred.index <- which(pred == 1)

###################################################
# downstream
###################################################
counts <- sim.doublet[,-pred.index]; dim(counts)
procdata <- log2(counts + 1)
lpsmclust <- exprmclust(as.matrix(procdata))
lpsorder <- TSCANorder(lpsmclust)
diffval <- difftest(sim.doublet,lpsorder)
#Selected differentially expressed genes under qvlue cutoff of 0.05
findDE <- row.names(diffval)[diffval$qval <= 0.05]; length(findDE)
findnonDE <- row.names(diffval)[diffval$qval > 0.05]; length(findnonDE)
nonDE <- rownames(counts)[1:500]
DE <- rownames(counts)[501:750]
tp <- length(intersect(findDE, DE)); tp
fp <- length(setdiff(findDE, DE)); fp
fn <- length(setdiff(findnonDE, nonDE)); fn
tn <- length(intersect(findnonDE, nonDE)); tn

precision <- tp/(tp+fp); precision
recall <- tp/(tp+fn); recall
tnr <- tn/(tn+fp); tnr







