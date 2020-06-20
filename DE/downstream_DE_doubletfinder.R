library(Matrix)
library(Seurat)
library(SingleCellExperiment)
library(DoubletFinder)
library(PRROC)
library(pbapply)
setwd('/media/nxi/nxi/doublet')
getwd()

sim.data <- readRDS('paper_sim/sim_DE.rds')
sim.doublet <- sim.data[[1]]; dim(sim.doublet)
cluster <- sim.data[[2]]; table(cluster)
de.up <- sim.data[[3]]; length(de.up)
de.down <- sim.data[[4]]; length(de.down)
de.truth <- c(de.up, de.down); length(de.truth)

###########################################################
# after cleaning
## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
doublet.seurat <- CreateSeuratObject(counts = sim.doublet, project = "doublet", 
                                     min.cells = 1, min.features = 1); doublet.seurat
doublet.seurat <- NormalizeData(doublet.seurat)
doublet.seurat <- ScaleData(doublet.seurat)
doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
doublet.seurat <- RunPCA(doublet.seurat)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
system.time(sweep.res.doublet <- paramSweep_v3(doublet.seurat, PCs = 1:10, sct = FALSE, num.cores = 44))
sweep.stats.doublet <- summarizeSweep(sweep.res.doublet, GT = FALSE)
bcmvn.doublet <- find.pK(sweep.stats.doublet)
pK <- bcmvn.doublet$pK[which.max(bcmvn.doublet$BCmetric)]; pK <- as.numeric(levels(pK))[pK]; pK
doublet.seurat <- doubletFinder_v3(doublet.seurat, PCs = 1:10, pN = 0.25, pK = pK, nExp = 667)
attribute <- paste('pANN', 0.25, pK, 667, sep = '_'); attribute
score <- doublet.seurat@meta.data[[attribute]]; score

###########################################################
# best threshold
t <- sort(score, decreasing = TRUE)[667]
pred.index <- which(score > t); length(pred.index)
cluster.clean <- cluster[-pred.index]; table(cluster.clean)
###########################################################
doublet.clean <- sim.doublet[,-pred.index]; dim(doublet.clean)
clean.seurat <- CreateSeuratObject(counts = doublet.clean, project = "doublet",
                                   min.cells = 1, min.features = 1); clean.seurat
clean.seurat <- NormalizeData(clean.seurat)
clean.seurat <- ScaleData(clean.seurat)
clean.seurat <- FindVariableFeatures(clean.seurat, selection.method = "vst", nfeatures = 2000)
clean.seurat <- RunPCA(clean.seurat)
Idents(clean.seurat) <- as.factor(cluster.clean); levels(clean.seurat)
clean.seurat[["RNA"]]@counts <- as.matrix(clean.seurat[["RNA"]]@counts) + 1
mode(clean.seurat[["RNA"]]@counts) <- 'integer'
system.time(marker.clean <- FindMarkers(clean.seurat, ident.1 = '0', ident.2 = "1", test.use = 'DESeq2', 
                                        slot="counts", max.cells.per.ident = 500, random.seed = 10))
#system.time(marker.clean <- FindMarkers(clean.seurat, ident.1 = '0', ident.2 = "1", test.use = 'wilcox'))
de.doublet.clean <- marker.clean[marker.clean$p_val_adj <= 0.05,]; dim(de.doublet.clean)
#de.doublet.clean <- marker.clean[marker.clean$power >= 0,]; dim(de.doublet.clean)

de.doublet.clean <- rownames(de.doublet.clean); length(de.doublet.clean)
tp.clean <- length(intersect(de.doublet.clean, de.truth)); tp.clean
fp.clean <- length(setdiff(de.doublet.clean, de.truth)); fp.clean
fn.clean <- length(setdiff(de.truth, de.doublet.clean)); fn.clean
nde.truth <- setdiff(rownames(sim.doublet), de.truth); length(nde.truth)
nde.clean <- setdiff(rownames(sim.doublet), de.doublet.clean); length(nde.clean)
tn.clean <- length(intersect(nde.truth, nde.clean)); tn.clean

precision.clean <- tp.clean / (tp.clean + fp.clean); precision.clean
recall.clean <- tp.clean / (tp.clean + fn.clean); recall.clean
tnr.clean <- tn.clean / (tn.clean + fn.clean); tnr.clean
f1.clean <- 2 * precision.clean * recall.clean / (precision.clean + recall.clean); f1.clean










