library(Matrix)
library(Seurat)
library(SingleCellExperiment)
library(scran)
library(PRROC)
library(reticulate)
library(pbapply)
set.seed(2020)

# read python module
doubletdetection <- import('doubletdetection')
np <- import('numpy')

# read simulation data with groud truth DE genes
sim.data <- readRDS('paper_sim/sim_DE.rds')
# count matrix
sim.doublet <- sim.data[[1]]; dim(sim.doublet)
# two cell types
cluster <- sim.data[[2]]; table(cluster)
# up and down regulated genes
de.up <- sim.data[[3]]; length(de.up)
de.down <- sim.data[[4]]; length(de.down)
# merge to get ground truth DE genes
de.truth <- c(de.up, de.down); length(de.truth)

####################################################################################################################
# DE analysis on the dataset without doublets (postive control)
# DE analysis on the dataset before doublet detection (negative control)
# The previous two functionalities were implemented in scrublet_DE.R
####################################################################################################################

######################################################################################################################
# DE analysis on the dataset after doublet detection
######################################################################################################################
# doubletdetection
doublet.seurat <- CreateSeuratObject(counts = sim.doublet, project = "doublet", 
                                     min.cells = 1, min.features = 1); doublet.seurat
raw_counts <- t(doublet.seurat[['RNA']]@counts)
clf <- doubletdetection$BoostClassifier(n_iters=5L, use_phenograph=FALSE, 
                                        standard_scaling=TRUE, random_state = 10L)
system.time(fit <- clf$fit(raw_counts))
p.values <- abs(fit$all_log_p_values_)
score <- apply(p.values, 2, mean)
# select threshold based on identification rate equal to the true doublet rate
score <- results[[1]]
t <- sort(score, decreasing = TRUE)[667]
# remove predicted doublet
pred.index <- which(score > t); length(pred.index)
cluster.clean <- cluster[-pred.index]; table(cluster.clean)
doublet.clean <- sim.doublet[,-pred.index]; dim(doublet.clean)
clean.seurat <- CreateSeuratObject(counts = doublet.clean, project = "doublet",
                                   min.cells = 1, min.features = 1); clean.seurat
clean.seurat <- NormalizeData(clean.seurat)
clean.seurat <- ScaleData(clean.seurat)
clean.seurat <- FindVariableFeatures(clean.seurat, selection.method = "vst", nfeatures = 2000)
clean.seurat <- RunPCA(clean.seurat)
Idents(clean.seurat) <- as.factor(cluster.clean); levels(clean.seurat)
# three DE methods: Wlicoxon, MAST, DESeq2; controled by test.use
system.time(marker.clean <- FindMarkers(clean.seurat, ident.1 = '0', ident.2 = "1", test.use = 'wilcox'))
de.doublet.clean <- marker.clean[marker.clean$p_val_adj <= 0.05,]; dim(de.doublet.clean)
de.doublet.clean <- rownames(de.doublet.clean); length(de.doublet.clean)

# calculate precision, recall, and true negative rate (TNR)
tp.clean <- length(intersect(de.doublet.clean, de.truth)); tp.clean
fp.clean <- length(setdiff(de.doublet.clean, de.truth)); fp.clean
fn.clean <- length(setdiff(de.truth, de.doublet.clean)); fn.clean
nde.truth <- setdiff(rownames(sim.doublet), de.truth); length(nde.truth)
nde.clean <- setdiff(rownames(sim.doublet), de.doublet.clean); length(nde.clean)
tn.clean <- length(intersect(nde.truth, nde.clean)); tn.clean

precision.clean <- tp.clean / (tp.clean + fp.clean); precision.clean
recall.clean <- tp.clean / (tp.clean + fn.clean); recall.clean
tnr.clean <- tn.clean / (tn.clean + fn.clean); tnr.clean
