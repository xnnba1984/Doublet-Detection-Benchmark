library(reticulate)
use_python('/home/nxi/anaconda3/bin', required = T)
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

times <- times * 5
saveRDS(times,'sim_result/doubletdetection_time.rds')

# 100k cells
index.gene.100k <- sample(1:dim(data)[1], 10000)
dim(data)
data.100k <- data[index.gene.100k,]; dim(data.100k)
data.100k <- cbind(data.100k,data.100k,data.100k,data.100k,data.100k,
                   data.100k,data.100k,data.100k,data.100k,data.100k); dim(data.100k)

system.time({
  #doublet.seurat <- CreateSeuratObject(counts = data.100k, project = "doublet", min.cells = 1, min.features = 1)
  #count <- as.matrix(doublet.seurat[['RNA']]@counts); dim(count)
  clf <- doubletdetection$BoostClassifier(n_iters=1L, use_phenograph=FALSE, standard_scaling=F, random_state = 3L)
  fit <- clf$fit(t(data.100k))
})
4818.881 * 5 / 3600

#####################################################################################
# running time
#####################################################################################
# read data; parameter setting
datasets <- c(
  #'real_data/pbmc-ch.rds',
  'real_data/cline-ch.rds',
  'real_data/mkidney-ch.rds',
  'real_data/hm-12k.rds',
  'real_data/hm-6k.rds',
  'real_data/pbmc-1A-dm.rds',
  'real_data/pbmc-1B-dm.rds',
  'real_data/pbmc-1C-dm.rds',
  'real_data/pbmc-2ctrl-dm.rds',
  'real_data/pbmc-2stim-dm.rds',
  'real_data/J293t-dm.rds',
  #'real_data/pdx-MULTI.rds',
  'real_data/HMEC-orig-MULTI.rds',
  'real_data/HMEC-rep-MULTI.rds',
  'real_data/HEK-HMEC-MULTI.rds',
  'real_data/nuc-MULTI.rds'
)
times <- pbsapply(datasets, function(loc){
  data <- readRDS(loc)
  count <- data[[1]]; dim(count)
  label <- data[[2]]; table(label)
  label <- ifelse(label == 'doublet', 1, 0); table(label)
  doublet.rate <- sum(label==1) / length(label); doublet.rate
  time <- system.time({
    doublet.seurat <- CreateSeuratObject(counts = count, project = "doublet", min.cells = 1, min.features = 1)
    count <- as.matrix(doublet.seurat[['RNA']]@counts); dim(count)
    clf <- doubletdetection$BoostClassifier(n_iters=1L, use_phenograph=FALSE, standard_scaling=F, random_state = 3L)
    fit <- clf$fit(t(count))
  })
  return(time[3])
}, simplify = T); times




















