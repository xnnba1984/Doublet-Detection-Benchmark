library(Matrix)
library(Seurat)
library(dplyr)
library(parallel)
library(reticulate)
library(PRROC)
library(pbapply)
library(scds)
library(scran)
library(DoubletFinder)
setwd('/media/nxi/nxi/doublet')
source('utility.R')
getwd()

# read python module
doubletdetection <- import('doubletdetection')
np <- import('numpy')

# read data; parameter setting
system.time(sim.data <- readRDS("sim_data/sim_293t_5000_2_distribute_1.rds"))
data <- sim.data[[1]][[1]]
set.seed(1)
index.gene <- sample(1:dim(data)[1], 5000)
number.cell <- seq(200, 10000, 200)

##########################################################################
# dblcell
##########################################################################
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

plot(time.doubletfinder)
model <- lm(time.doubletfinder~I(number.cell^1.5))
summary(model)
saveRDS(time.doubletfinder, "sim_result/doubletfinder_time.rds")

# 100k cells
index.gene.100k <- sample(1:dim(data)[1], 10000)
dim(data)
data.100k <- data[index.gene.100k,]; dim(data.100k)
data.100k <- cbind(data.100k,data.100k,data.100k,data.100k,data.100k,
                   data.100k,data.100k,data.100k,data.100k,data.100k); dim(data.100k)

time <- system.time({
  ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
  doublet.seurat <- CreateSeuratObject(counts = data.100k, project = "doublet", min.cells = 1, min.features = 1); doublet.seurat
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



















