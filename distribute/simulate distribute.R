library(devtools)
library(scDesign)
library(Matrix)
library(Seurat)
library(simulator)
setwd('/media/nxi/nxi/doublet')
getwd()

# read raw counte matrix
raw.matrix <- as.matrix(readMM('data/293t_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/matrix.mtx'))
dim(raw.matrix)

# drop non-expressed genes
raw.matrix <- raw.matrix[apply(raw.matrix, 1, function(x) !all(x==0)), ]; dim(raw.matrix)
# sequencing depth
depth <- sum(raw.matrix) / ncol(raw.matrix)
cores <- 44
seed <- 7
# estimate model parameters
set.seed(seed)
system.time(estpa <- estimate_pa(raw.matrix, ncores = cores))

# parameter setting
cell.number <- 6000
cell.type <- 2
squence.depth <- round(depth*cell.number)
doublet.rate <- 0.2
doublet.number <- round(cell.number * cell.type * doublet.rate)
sim.number <- 1

# simulate count matrix
sim.doublet.list <- list()
sim.label.list <- list()
sim.type.list <- list()
set.seed(seed)
system.time(
  for(i in 1:sim.number){
    print(paste("create simualtion data", i))
    system.time(sim <- simulate_de_mfo(raw.matrix, estpa, S = rep(squence.depth, cell.type), 
                                       Js = rep(cell.number, cell.type),
                                       ngroup = cell.type, pUp = 0.03, pDown = 0.03, fU = 3, fL = 1.5))
    sim <- do.call(cbind, sim$count); dim(sim); sum(sim) / ncol(sim)
    
    # drop non-expressed genes
    doublet.seurat <- CreateSeuratObject(counts = sim, project = "doublet", min.cells = 1, min.features = 1)
    sim <- as.matrix(doublet.seurat[['RNA']]@counts); dim(sim)
    
    # keep cell type labels
    type.label <- as.vector(sapply(0:(cell.type-1), function(x){
      return(rep(x, cell.number))
    }, simplify = TRUE))
    
    # cell merage ratio: for each cell c = d*(a*c1 + b*c2), a+1=2, a~N(1, 0.1), d~N(1, 0.1)
    a <- rnorm(doublet.number, 1, 0.1)
    b <- 2 - a
    c <- rnorm(doublet.number, 1, 0.1)
    
    # sample doublets
    doublet.pair <- combn(ncol(sim), 2); dim(doublet.pair)
    doublet.pair <- doublet.pair[,sample(ncol(doublet.pair), doublet.number)]; dim(doublet.pair)
    doublets <- sapply(1:ncol(doublet.pair), function(pair) {
      pair.1 <- doublet.pair[,pair][1]
      pair.2 <- doublet.pair[,pair][2]
      doublet.cell <- round(c[pair]*(a[pair]*sim[,pair.1] + b[pair]*sim[,pair.2]))
      return(doublet.cell)
    }); class(doublets); dim(doublets)
    
    # label homo and hetero doublets; 1: homo, 2: hetero
    doublet.type <- sapply(1:ncol(doublet.pair), function(pair) {
      pair.1 <- doublet.pair[,pair][1]
      pair.2 <- doublet.pair[,pair][2]
      type <- ifelse(ceiling(pair.1 / cell.number) !=  ceiling(pair.2 / cell.number), 2, 1)
      return(type)
    }); sum(doublet.type==2); sum(doublet.type==1); doublet.type
    
    # drop doublets ingradients; combine doublets and non-doublets
    drop.index <- (unique(as.vector(doublet.pair))); length(drop.index)
    length.singlet <- ncol(sim)-length(drop.index); length.singlet
    sim.doublet <- sim[,-drop.index]; dim(sim.doublet)
    sim.doublet <- cbind(sim.doublet, doublets); dim(sim.doublet)
    colnames(sim.doublet) <- as.character(1:ncol(sim.doublet))
    
    # add cell type labels to doublets
    type.label <- type.label[-drop.index]
    type.label <- c(type.label, rep(cell.type, doublet.number))
    
    # doublet true label
    doublet.label <- c(rep(0, length.singlet), doublet.type); length(doublet.label)
    
    # save result to lists
    sim.doublet.list <- append(sim.doublet.list, list(sim.doublet))
    sim.label.list <- append(sim.label.list, list(doublet.label))
    sim.type.list <- append(sim.type.list, list(type.label))
  })

# sim.doublet.list: expression matrix; sim.label.list: singlet, homo doublet, heter doublet label; 
# sim.type.list: singlet doublet label
sim.data <- list(sim.doublet.list, sim.label.list, sim.type.list)

system.time(saveRDS(sim.data, "sim_data/sim_293t_5000_2_distribute_1.rds"))








