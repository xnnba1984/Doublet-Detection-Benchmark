library(devtools)
library(scDesign)
library(Matrix)
library(Seurat)
library(simulator)
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
raw.matrix <- as.matrix(readMM('data/293t_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/matrix.mtx')); dim(raw.matrix)

# drop non-expressed genes
raw.matrix <- raw.matrix[apply(raw.matrix, 1, function(x) !all(x==0)), ]; dim(raw.matrix)

# estimate model parameters
system.time(estpa <- estimate_pa(raw.matrix, ncores = 44))

# parameter setting
cell.number <- 500
#cell.types <- seq(2, 20, 1)
cell.type <- 2
lsize <- 2000
#lsizes <- seq(500, 10000, 500)
squence.depth <- lsize * cell.number
#doublet.rates <- seq(0.02, 0.4, 0.02)
doublet.rate <- .2
cores <- 44
seed <- 10
sim.num <- 5
pUp <- seq(from = .01, by = .002, to = 0.05)
pDown <- seq(from = .01, by = .002, to = 0.05)
fU <- seq(from = 1, by = .2, to = 5)
fL <- seq(from = .5, by = .1, to = 2.5)

sim.list <- list()
system.time({
for(i in 1:length(pUp)){
  print(paste('==========', pUp[i]))
  sim.doublet.list <- list()
  sim.label.list <- list()
  doublet.number <- cell.number * cell.type * doublet.rate
  for(j in 1:sim.num){
    print(j)
    # simulate count matrix
    system.time(sim <- simulate_de_mfo(raw.matrix, estpa, S = rep(squence.depth, cell.type), 
                                       Js = rep(cell.number, cell.type),
                                       ngroup = cell.type, pUp = pUp[i], pDown = pDown[i], fU = fU[i], fL = fL[i]))
    sim <- do.call(cbind, sim$count); dim(sim); sum(sim) / ncol(sim)
    sim <- as(sim, "sparseMatrix") 
    
    # keep cell type labels
    type.label <- as.vector(sapply(0:(cell.type-1), function(x){
      return(rep(x, cell.number))
    }, simplify = TRUE)); table(type.label)
    
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
      return(doublet.cell)
    }); doublets <- as(doublets, 'sparseMatrix'); dim(doublets)
    
    mean(colSums(doublets>0)) / mean(colSums(sim>0)[1:200])
    mean(colSums(doublets)) / mean(colSums(sim)[1:200])
    
    # drop doublets ingradients; combine doublets and non-doublets
    drop.index <- (unique(as.vector(doublet.pair))); length(drop.index)
    sim.singlets <- sim[,-drop.index]; dim(sim.singlets)
    type.label <- type.label[-drop.index]; table(type.label)
    type.label <- c(type.label, rep('doublet', ncol(doublets))); table(type.label)
    sim.data <- cbind(sim.singlets, doublets); dim(sim.data)
    colnames(sim.data) <- as.character(1:ncol(sim.data))
    
    # save result to lists
    sim.doublet.list <- append(sim.doublet.list, list(sim.data))
    sim.label.list <- append(sim.label.list, list(type.label))
  }
  sim.list[[as.character(i)]] <- list(sim.doublet.list, sim.label.list)
}
})
saveRDS(sim.list, 'paper_sim/sim_hetero.rds')

# save sim data to loom
library(loomR)
system.time({
for(i in 1:length(sim.list)){
  print(paste('===================', i))
  counts <- sim.list[[i]][[1]]
  for(j in 1:length(counts)){
    print(j)
    count <- counts[[j]]
    count <- count[which(rowSums(count) != 0),]
    print(dim(count))
    create(paste(paste('paper_sim/loom/hetero', i, j, sep='_'), 'loom', sep='.'), count, transpose = T)
  }
}
})







