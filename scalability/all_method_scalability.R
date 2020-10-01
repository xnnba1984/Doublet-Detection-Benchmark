library(reticulate)
library(Matrix)
library(Seurat)
library(dplyr)
library(parallel)
library(PRROC)
library(pbapply)
library(scds)
library(scran)
library(DoubletFinder)
source('utility.R')
set.seed(2020)

# read python module
scr <- import('scrublet')
scipy.io <- import('scipy.io')
np <- import('numpy')
os <- import('os')
doubletdetection <- import('doubletdetection')

# read simulation data with 10000 cells
system.time(sim.data <- readRDS("sim_data/sim_293t_5000_2_distribute_1.rds"))
data <- sim.data[[1]][[1]]
# 5000 genes
index.gene <- sample(1:dim(data)[1], 5000)
# cell number from 400 to 10000, by 400
number.cell <- seq(400, 10000, 400)

##########################################################################
# scrublet
##########################################################################
# loop over each cell number
time.scrublet <- pbsapply(number.cell, function(cell){
  count <- data[index.gene, sample(1:dim(data)[2], cell)]; dim(count)
  time <- system.time({
    count <- as(count, "sparseMatrix"); dim(count)
    result <- scr$Scrublet(counts_matrix = t(count), expected_doublet_rate = 0.2, random_state = 10L)
    results <- result$scrub_doublets(min_counts=2, 
                                       min_cells=3, 
                                       min_gene_variability_pctl=85, 
                                       n_prin_comps=30L)
  })
  return(time[3])
}, simplify = T)
# save running time
saveRDS(time.scrublet,'sim_result/scrublet_time.rds')

##########################################################################
# doubletcells
##########################################################################
# loop over each cell number
time.dblcell <- pbsapply(number.cell, function(cell){
  count <- data[index.gene, sample(1:dim(data)[2], cell)]; dim(count)
  time <- system.time({
    score <- doubletCells(count)
  })
  return(time[3])
}, simplify = T); time.dblcell

plot(time.dblcell)
model <- lm(time.dblcell~number.cell)
summary(model)
saveRDS(time.dblcell,'sim_result/dblcell_time.rds')

##########################################################################
# cxds
##########################################################################
# loop over each cell number
time.cxds <- pbsapply(number.cell, function(cell){
  count <- data[index.gene, sample(1:dim(data)[2], cell)]; dim(count)
  time <- system.time({
    sce <- SingleCellExperiment(assays = list(counts = count))
    sce <- cxds(sce)
  })
  return(time[3])
}, simplify = T); time.cxds
# save running time
saveRDS(time.cxds,'sim_result/cxds_time.rds')

##########################################################################
# bcds
##########################################################################
# loop over each cell number
time.bcds <- pbsapply(number.cell, function(cell){
  count <- data[index.gene, sample(1:dim(data)[2], cell)]; dim(count)
  time <- system.time({
    sce <- SingleCellExperiment(assays = list(counts = count))
    sce <- bcds(sce)
  })
  return(time[3])
}, simplify = T); time.bcds
# save running time
saveRDS(time.bcds,'sim_result/bcds_time.rds')

##########################################################################
# hybrid
##########################################################################
# loop over each cell number
time.hybrid <- pbsapply(number.cell, function(cell){
  count <- data[index.gene, sample(1:dim(data)[2], cell)]; dim(count)
  time <- system.time({
    sce <- SingleCellExperiment(assays = list(counts = count))
    sce <- cxds_bcds_hybrid(sce)
  })
  return(time[3])
}, simplify = T); time.hybrid
# save running time
saveRDS(time.hybrid,'sim_result/hybrid_time.rds')

##########################################################################
# doubletdetection
##########################################################################
# loop over each cell number
times <- vector()
system.time({
  for(cell in number.cell){
    print(cell)
    count <- data[index.gene, sample(1:dim(data)[2], cell)]; dim(count)
    time <- system.time({
      doublet.seurat <- CreateSeuratObject(counts = count, project = "doublet", min.cells = 1, min.features = 1)
      count <- as.matrix(doublet.seurat[['RNA']]@counts); dim(count)
      clf <- doubletdetection$BoostClassifier(n_iters=1L, use_phenograph=FALSE, standard_scaling=F, random_state = 3L)
      fit <- clf$fit(t(count))
    })
    times <- append(times, time[3])
  }
});times
# the repeat parameter was set to 1 in doubletdetection to accelerate the running
# times 5 to compenstate
times <- times * 5
saveRDS(times,'sim_result/doubletdetection_time.rds')

##########################################################################
# doubletfinder
##########################################################################
# loop over each cell number
time.doubletfinder <- pbsapply(number.cell, function(cell){
  count <- data[index.gene, sample(1:dim(data)[2], cell)]; dim(count)
  time <- system.time({
    ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
    doublet.seurat <- CreateSeuratObject(counts = count, project = "doublet", min.cells = 1, min.features = 1); doublet.seurat
    doublet.seurat <- NormalizeData(doublet.seurat)
    doublet.seurat <- ScaleData(doublet.seurat)
    doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
    doublet.seurat <- RunPCA(doublet.seurat)
    
    ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
    sweep.res.doublet <- paramSweep_v3(doublet.seurat, PCs = 1:10, sct = FALSE, num.cores = 44)
    sweep.stats.doublet <- summarizeSweep(sweep.res.doublet, GT = FALSE)
    bcmvn.doublet <- find.pK(sweep.stats.doublet)
    pK <- bcmvn.doublet$pK[which.max(bcmvn.doublet$BCmetric)]; pK <- as.numeric(levels(pK))[pK]; pK
    doublet.seurat <- doubletFinder_v3(doublet.seurat, PCs = 1:10, pN = 0.25, pK = pK, nExp = 200)
    attribute <- paste('pANN', 0.25, pK, 200, sep = '_'); attribute
    score <- doublet.seurat@meta.data[[attribute]]; score
  })
  return(time[3])
}, simplify = T); time.doubletfinder
saveRDS(time.doubletfinder, "sim_result/doubletfinder_time.rds")

##########################################################################
# The running time of Solo was recorded in linux command
##########################################################################

