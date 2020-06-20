library(reticulate)
use_python('/home/nxi/anaconda3/bin', required = T)
library(Matrix)
library(PRROC)
library(RcppCNPy)
setwd('/media/nxi/nxi/doublet')
getwd()
np <- import("numpy")

# read data; parameter setting
system.time(sim.data <- readRDS("sim_data/sim_293t_5000_2_distribute_1.rds"))
data <- sim.data[[1]][[1]]; dim(data)
data <- as(data, "sparseMatrix"); dim(data)
set.seed(1)
index.gene <- sample(1:dim(data)[1], 5000)
number.cell <- seq(400, 10000, 400)

# save sim data to loom
library(loomR)
system.time({
  for(i in 1:length(number.cell)){
      print(paste('===================', i))
      count <- data[index.gene, sample(1:dim(data)[2], number.cell[i])]; dim(count)
      create(paste(paste('paper_sim/loom/time', number.cell[i], sep = '_'), 'loom', sep='.'), 
             count, transpose = T)
    }
})

times <- c(15, 21, 39, 43, 61, 93, 101, 88, 77, 98, 123, 162, 145, 176, 164, 184, 229, 161, 214, 218, 205, 424, 
           408, 341, 340)
saveRDS(times, 'sim_result/solo_time.rds')
readRDS('sim_result/solo_time.rds')

# scalibility
for(i in seq(from = 400, to = 10000, by = 400)){
  file <- paste('paper_sim/loom/output_time/time_', i, '.loom', '/softmax_scores.npy', sep = '')
  print(file)
}







