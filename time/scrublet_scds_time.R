library(reticulate)
use_python('/home/nxi/anaconda3/bin', required = T)
library(Matrix)
library(Seurat)
library(dplyr)
library(parallel)
library(PRROC)
library(pbapply)
library(scds)
library(scran)
library(DoubletFinder)
setwd('/media/nxi/nxi/doublet')
source('utility.R')
getwd()

# read python module
scr <- import('scrublet')
scipy.io <- import('scipy.io')
np <- import('numpy')
os <- import('os')

# read data; parameter setting
system.time(sim.data <- readRDS("sim_data/sim_293t_5000_2_distribute_1.rds"))
data <- sim.data[[1]][[1]]
set.seed(1)
#index.gene <- sample(1:dim(data)[1], 5000)
#number.cell <- seq(200, 10000, 200)
number.cell <- sample(1:dim(data)[2], 1000)
number.gene <- seq(500, 15000, 1000)

##########################################################################
# scrublet
##########################################################################
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

plot(time.scrublet)
model <- lm(time.scrublet~I(number.cell^1.2));summary(model)
saveRDS(time.scrublet,'sim_result/scrublet_time.rds')

##########################################################################
# cxds
##########################################################################
time.cxds <- pbsapply(number.cell, function(cell){
  count <- data[index.gene, sample(1:dim(data)[2], cell)]; dim(count)
  time <- system.time({
    sce <- SingleCellExperiment(assays = list(counts = count))
    sce <- cxds(sce)
  })
  return(time[3])
}, simplify = T); time.cxds

plot(time.cxds)
model <- lm(time.cxds~I(number.cell^1.3));summary(model)
saveRDS(time.cxds,'sim_result/cxds_time.rds')

##########################################################################
# bcds
##########################################################################
time.bcds <- pbsapply(number.cell, function(cell){
  count <- data[index.gene, sample(1:dim(data)[2], cell)]; dim(count)
  time <- system.time({
    sce <- SingleCellExperiment(assays = list(counts = count))
    sce <- bcds(sce)
  })
  return(time[3])
}, simplify = T); time.bcds

plot(time.bcds)
model <- lm(time.bcds~I(number.cell^.5));summary(model)
saveRDS(time.bcds,'sim_result/bcds_time.rds')
##########################################################################
# hybrid
##########################################################################
time.hybrid <- pbsapply(number.cell, function(cell){
  count <- data[index.gene, sample(1:dim(data)[2], cell)]; dim(count)
  time <- system.time({
    sce <- SingleCellExperiment(assays = list(counts = count))
    sce <- cxds_bcds_hybrid(sce)
  })
  return(time[3])
}, simplify = T); time.hybrid

plot(time.hybrid)
model <- lm(time.hybrid~I(number.cell^.5));summary(model)
saveRDS(time.hybrid,'sim_result/hybrid_time.rds')

time.dblcell <- readRDS('sim_result/dblcell_time.rds')
time.doubletfinder <- readRDS('sim_result/doubletfinder_time.rds')
time.doubletdetection <- readRDS('sim_result/doubletdetection_time.rds')
plot(time.doubletdetection)

# search best complexity
exp <- seq(.1,2,.1)
adj.r2 <- vector()
for(e in exp){
  model <- lm(time.doubletdetection~I(number.cell^e))
  adj.r2 <- append(adj.r2, summary(model)$adj.r.squared)
}
max(adj.r2); exp[which(adj.r2==max(adj.r2))]

# 100k cells
index.gene.100k <- sample(1:dim(data)[1], 10000)
dim(data)
data.100k <- data[index.gene.100k,]; dim(data.100k)
data.100k <- cbind(data.100k,data.100k,data.100k,data.100k,data.100k,
                   data.100k,data.100k,data.100k,data.100k,data.100k); dim(data.100k)

# scrublet
count <- as(data.100k, "sparseMatrix"); dim(data.100k)
system.time({
  result <- scr$Scrublet(counts_matrix = t(count), expected_doublet_rate = 0.2, random_state = 10L)
  results <- result$scrub_doublets(min_counts=2, 
                                   min_cells=3, 
                                   min_gene_variability_pctl=85, 
                                   n_prin_comps=30L)
})

# cxds
system.time({
  sce <- SingleCellExperiment(assays = list(counts = data.100k))
  sce <- cxds_bcds_hybrid(sce)
})




