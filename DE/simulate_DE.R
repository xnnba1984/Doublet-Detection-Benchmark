library(devtools)
library(scDesign)
library(Matrix)
library(Seurat)
library(simulator)
library(pbapply)
setwd('/media/nxi/nxi/doublet')
getwd()

# get real data summary
real <- readRDS('paper_sim/real_summary.rds')
fit <- density(real$rsize)
bw.size <- fit$bw
fit <- density(real$rgene)
bw.gene <- fit$bw
mean(real$rsize); mean(real$rgene)

# read raw counte matrix
raw.matrix <- as.matrix(readMM('data/293t_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/matrix.mtx'))
dim(raw.matrix)

cell.number <- 500
doublet.rate <- .4
doublet.number <- round(1000 * doublet.rate / (1 - doublet.rate))
lsize <- 10000
squence.depth <- lsize * cell.number
cell.type <- 2
pUp <- .03
pDown <- .03
fU <- 3
fL <- 1.5 
ncores <- 44

set.seed(10)
# drop non-expressed genes
raw.matrix <- raw.matrix[apply(raw.matrix, 1, function(x) !all(x==0)), ]; dim(raw.matrix)
system.time(simdata <- design_data(realcount = raw.matrix, S = rep(squence.depth, cell.type), 
                                  ncell = rep(cell.number, cell.type), ngroup = cell.type, 
                                  pUp = pUp, pDown = pDown, fU = fU, fL = fL, ncores = ncores))
sim1 <- simdata$count[[1]]; dim(sim1)
sim2 <- simdata$count[[2]]; dim(sim2)
sim <- cbind(sim1, sim2); dim(sim)
sim <- as(sim, "sparseMatrix") 

# sample doublets index
doublet.pair <- combn(ncol(sim), 2); dim(doublet.pair)
doublet.pair <- doublet.pair[,sample(ncol(doublet.pair), doublet.number)]; dim(doublet.pair)

# sample doublets lsize ratio
means <- sample(real$rsize, doublet.number, replace = TRUE)
r.size <- rnorm(doublet.number, mean = means, sd = bw.size)

# sample doublets ngene ratio
ngenes <- sample(real$rgene, doublet.number, replace = TRUE)
r.genes <- rnorm(doublet.number, mean = means, sd = bw.gene)

# sample doublets
doublets <- sapply(1:ncol(doublet.pair), function(pair) {
  pair.1 <- doublet.pair[,pair][1]
  pair.2 <- doublet.pair[,pair][2]
  # similate lsize ratio
  doublet.cell <- round((sim[,pair.1] + sim[,pair.2]) * r.size[pair] / 2, digits = 5)
  # simulate ngene ratio
  gap <- round(sum(doublet.cell>0) - mean(c(sum(sim[,pair.1]>0), sum(sim[,pair.2]>0))) * r.genes[pair]); gap
  if(gap > 0 ){
    gene.express <- which(doublet.cell > 0)
    gene.drop <- sample(gene.express, gap, replace = F)
    doublet.cell[gene.drop] = 0
  }
  #doublet.cell <- sim[,pair.1] + sim[,pair.2]
  return(doublet.cell)
}); doublets <- as(doublets, 'sparseMatrix'); dim(doublets)

index1 <- sample(1:ncol(doublets), round(ncol(doublets)/2)); length(index1)
dim(doublets[,index1])
sim1 <- cbind(sim1, doublets[,index1]); dim(sim1)
sim2 <- cbind(sim2, doublets[,-index1]); dim(sim2)
sim.doublet <- cbind(sim1, sim2); dim(sim.doublet)
colnames(sim.doublet) <- as.character(1:ncol(sim.doublet))
cluster <- c(rep(0, cell.number + dim(doublets[,index1])[2]), rep(1, cell.number + dim(doublets[,-index1])[2])); table(cluster)

# save result to lists
sim.data <- list(sim.doublet, cluster, simdata[["genesUp"]][[2]], simdata[["genesDown"]][[2]])
saveRDS(sim.data, 'paper_sim/sim_DE.rds')

# save sim data to loom
library(loomR)
create('paper_sim/loom/DE.loom', sim.doublet, transpose = T)














