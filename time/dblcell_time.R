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
scr <- import('scrublet')
scipy.io <- import('scipy.io')
np <- import('numpy')
os <- import('os')

# read data; parameter setting
system.time(sim.data <- readRDS("sim_data/sim_293t_5000_2_distribute_1.rds"))
data <- sim.data[[1]][[1]]
set.seed(1)
index.gene <- sample(1:dim(data)[1], 5000)
number.cell <- seq(200, 10000, 200)

##########################################################################
# dblcell
##########################################################################
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

# 100k cells
index.gene.100k <- sample(1:dim(data)[1], 10000)
dim(data)
data.100k <- data[index.gene.100k,]; dim(data.100k)
data.100k <- cbind(data.100k,data.100k,data.100k,data.100k,data.100k,
                   data.100k,data.100k,data.100k,data.100k,data.100k); dim(data.100k)
system.time({
  score <- doubletCells(data.100k)
})




















