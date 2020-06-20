library(reticulate)
use_python('/home/nxi/anaconda3/bin', required = T)
library(Matrix)
library(Seurat)
library(pbapply)
library(ggplot2)
library(parallel)
setwd('/media/nxi/nxi/doublet')
getwd()

# read python module
scr <- import('scrublet')
scipy.io <- import('scipy.io')
np <- import('numpy')
os <- import('os')

sim.data <- readRDS('paper_sim/sim_DE.rds')
sim.doublet <- sim.data[[1]]; dim(sim.doublet)
cluster <- sim.data[[2]]; table(cluster)
de.up <- sim.data[[3]]; length(de.up)
de.down <- sim.data[[4]]; length(de.down)
de.truth <- c(de.up, de.down); length(de.truth)

##########################################################
# no doublet
##########################################################
nodoublet <- sim.doublet[,c(1:500, 835:1334)]; dim(nodoublet)
#nodoublet <- sim.doublet[,c(1:500, 626:1125)]; dim(nodoublet)
doublet.seurat <- CreateSeuratObject(counts = nodoublet, project = "doublet", min.cells = 1, 
                                     min.features = 1); doublet.seurat
nodoublet <- doublet.seurat[["RNA"]]@counts; dim(nodoublet)
doublet.seurat <- NormalizeData(doublet.seurat)
doublet.seurat <- ScaleData(doublet.seurat)
doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
doublet.seurat <- RunPCA(doublet.seurat)
Idents(doublet.seurat) <- as.factor(c(rep(0,500),rep(1,500))); levels(doublet.seurat)
#doublet.seurat[["RNA"]]@counts <- as.matrix(doublet.seurat[["RNA"]]@counts) + 1
#mode(doublet.seurat[["RNA"]]@counts) <- 'integer'
#system.time(marker.doublet <- FindMarkers(doublet.seurat, ident.1 = '0', ident.2 = "1", test.use = 'DESeq2', 
#                                          slot="counts", max.cells.per.ident = 500, random.seed = 10))
system.time(marker.doublet <- FindMarkers(doublet.seurat, ident.1 = '0', ident.2 = "1", test.use = 'wilcox'))
de.doublet <- marker.doublet[marker.doublet$p_val_adj <= 0.05,]; dim(de.doublet)
#de.doublet <- marker.doublet[marker.doublet$power > 0.2,]; dim(de.doublet)
de.doublet <- rownames(de.doublet)

tp <- length(intersect(de.doublet, de.truth)); tp
fp <- length(setdiff(de.doublet, de.truth)); fp
fn <- length(setdiff(de.truth, de.doublet)); fn
nde.truth <- setdiff(rownames(nodoublet), de.truth); length(nde.truth)
nde.doublet <- setdiff(rownames(nodoublet), de.doublet); length(nde.doublet)
tn <- length(intersect(nde.truth, nde.doublet)); tn

precision <- tp / (tp + fp); precision
recall <- tp / (tp + fn); recall
tnr <- tn / (tn + fn); tnr
f1 <- 2 * precision * recall / (precision + recall); f1

##########################################################
# before cleaning
##########################################################
doublet.seurat <- CreateSeuratObject(counts = sim.doublet, project = "doublet", min.cells = 1, 
                                     min.features = 1); doublet.seurat
doublet.seurat <- NormalizeData(doublet.seurat)
doublet.seurat <- ScaleData(doublet.seurat)
doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
doublet.seurat <- RunPCA(doublet.seurat)
Idents(doublet.seurat) <- as.factor(cluster); levels(doublet.seurat)
#doublet.seurat[["RNA"]]@counts <- as.matrix(doublet.seurat[["RNA"]]@counts) + 1
#mode(doublet.seurat[["RNA"]]@counts) <- 'integer'
#system.time(marker.doublet <- FindMarkers(doublet.seurat, ident.1 = '0', ident.2 = "1", test.use = 'DESeq2', 
#                                          slot="counts", max.cells.per.ident = 500, random.seed = 10))
system.time(marker.doublet <- FindMarkers(doublet.seurat, ident.1 = '0', ident.2 = "1", test.use = 'wilcox'))
de.doublet <- marker.doublet[marker.doublet$p_val_adj <= 0.05,]; dim(de.doublet)
#de.doublet <- marker.doublet[marker.doublet$power > 0.2,]; dim(de.doublet)
de.doublet <- rownames(de.doublet)

tp <- length(intersect(de.doublet, de.truth)); tp
fp <- length(setdiff(de.doublet, de.truth)); fp
fn <- length(setdiff(de.truth, de.doublet)); fn
nde.truth <- setdiff(rownames(nodoublet), de.truth); length(nde.truth)
nde.doublet <- setdiff(rownames(nodoublet), de.doublet); length(nde.doublet)
tn <- length(intersect(nde.truth, nde.doublet)); tn

precision <- tp / (tp + fp); precision
recall <- tp / (tp + fn); recall
tnr <- tn / (tn + fn); tnr
f1 <- 2 * precision * recall / (precision + recall); f1

###########################################################
# after cleaning
###########################################################
data <- as(sim.doublet, "sparseMatrix"); dim(data)
result <- scr$Scrublet(counts_matrix = t(data), expected_doublet_rate = 0.2, random_state = 10L)
results <- result$scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, 
                                 n_prin_comps=30L)

#pred <- as.vector(results[[2]]); table(pred)
#pred.index <- which(pred==TRUE); length(pred.index)
#cluster.clean <- cluster[-pred.index]; table(cluster.clean)
###########################################################
# best threshold
score <- results[[1]]
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
#clean.seurat[["RNA"]]@counts <- as.matrix(clean.seurat[["RNA"]]@counts) + 1
#mode(clean.seurat[["RNA"]]@counts) <- 'integer'
#system.time(marker.clean <- FindMarkers(clean.seurat, ident.1 = '0', ident.2 = "1", test.use = 'DESeq2', 
#                                          slot="counts", max.cells.per.ident = 500, random.seed = 10))
system.time(marker.clean <- FindMarkers(clean.seurat, ident.1 = '0', ident.2 = "1", test.use = 'wilcox'))
de.doublet.clean <- marker.clean[marker.clean$p_val_adj <= 0.05,]; dim(de.doublet.clean)
#de.doublet.clean <- marker.clean[marker.clean$power > 0.4,]; dim(de.doublet.clean)

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










