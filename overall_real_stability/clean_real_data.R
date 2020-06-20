library(Matrix)
library(Seurat)
library(SingleCellExperiment)
library(scran)
library(PRROC)
library(pbapply)
library(deMULTIplex)
library(ggplot2)
library(biomaRt)
setwd('/media/nxi/nxi/doublet')
getwd()

#################################################
# ch_pbmc
#################################################
# Load in the UMI matrix
pbmc.umis <- readRDS("real_data/pbmc_umi_mtx.rds"); dim(pbmc.umis)
summary(colSums(pbmc.umis))

# Load in the HTO count matrix
pbmc.htos <- readRDS("pbmc_hto_mtx.rds"); dim(pbmc.htos)
#rownames(pbmc.htos)
# Select cell barcodes detected by both RNA and HTO In the example datasets 
joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos)); length(joint.bcs)

# Subset RNA and HTO counts by joint cell barcodes
pbmc.umis <- pbmc.umis[, joint.bcs]; dim(pbmc.umis)
pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs]); dim(pbmc.htos)
pbmc.htos[,1:10]
# Confirm that the HTO have the correct names
rownames(pbmc.htos)

# Setup Seurat object
pbmc.hashtag <- CreateSeuratObject(counts = pbmc.umis, min.cells = 1); pbmc.hashtag

# Normalize RNA data with log normalization
pbmc.hashtag <- NormalizeData(pbmc.hashtag)
# Find and scale variable features
pbmc.hashtag <- FindVariableFeatures(pbmc.hashtag, selection.method = "mean.var.plot")
pbmc.hashtag <- ScaleData(pbmc.hashtag, features = VariableFeatures(pbmc.hashtag)); dim(pbmc.hashtag)

# Add HTO data as a new assay independent from RNA
pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = pbmc.htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR"); dim(pbmc.hashtag)

pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)
pbmc.hashtag <- subset(pbmc.hashtag, idents = "Negative", invert = TRUE)
label <- pbmc.hashtag$HTO_classification.global; length(label); table(label)
counts <- pbmc.hashtag[['RNA']]@counts; dim(counts)
counts <- CreateSeuratObject(counts = counts, min.cells = 1); counts
counts <- counts[['RNA']]@counts; dim(counts)
sum(colnames(counts) != names(label))
summary(colSums(counts))
summary(colSums(counts != 0))
label <- ifelse(label=='Doublet', 'doublet', 'singlet')

# save to file
saveRDS(list(counts, label), 'real_data/pbmc-ch.rds')

#################################################
# ch_cell_line
#################################################
#Read in UMI count matrix for RNA
hto12.umis <- readRDS("real_data/hto12_umi_mtx.rds"); dim(hto12.umis)
summary(colSums(hto12.umis))
#hist(colSums(hto12.umis), breaks = 100)
# Read in HTO count matrix
hto12.htos <- readRDS("real_data/hto12_hto_mtx.rds"); dim(hto12.htos)

# Select cell barcodes detected in both RNA and HTO
cells.use <- intersect(rownames(hto12.htos), colnames(hto12.umis)); length(cells.use)

# Create Seurat object and add HTO data
hto12 <- CreateSeuratObject(counts = hto12.umis[, cells.use], min.features = 300, min.cells = 1); hto12
hto12[["HTO"]] <- CreateAssayObject(counts = t(x = hto12.htos[colnames(hto12), 1:12])); dim(hto12[["HTO"]])

# Normalize data
hto12 <- NormalizeData(hto12)
hto12 <- NormalizeData(hto12, assay = "HTO", normalization.method = "CLR")
hto12 <- HTODemux(hto12, assay = "HTO", positive.quantile = 0.99)
hto12 <- subset(hto12, idents = "Negative", invert = TRUE)
label <- hto12$HTO_classification.global; length(label); table(label)
counts <- hto12[['RNA']]@counts; dim(counts)

hto12 <- CreateSeuratObject(counts = counts, min.cells = 1); hto12
counts <- hto12[['RNA']]@counts; dim(counts)
sum(colnames(counts) != names(label))
summary(colSums(counts))
summary(colSums(counts != 0))
label <- ifelse(label=='Doublet', 'doublet', 'singlet')

# save to file
saveRDS(list(counts, label), 'real_data/cline-ch.rds')

#####################################################################
# hg_mm_12k / 6k
#####################################################################
library(biomaRt)

hg.barcode <- read.table('real_data/hgmm_12k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/barcodes.tsv',
                         stringsAsFactors = F)
hg.barcode <- hg.barcode$V1; length(hg.barcode)
hg.counts <- readMM('real_data/hgmm_12k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/matrix.mtx'); dim(hg.counts)
hg.genes <- read.table('real_data/hgmm_12k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/genes.tsv', 
                       stringsAsFactors = F)
hg.genes <- hg.genes$V1; length(hg.genes)
hg.genes <- substring(hg.genes, 6); length(hg.genes)
rownames(hg.counts) <- hg.genes
colnames(hg.counts) <- hg.barcode

mm.barcode <- read.table('real_data/hgmm_12k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/mm10/barcodes.tsv',
                         stringsAsFactors = F)
mm.barcode <- mm.barcode$V1; length(mm.barcode) 
mm.counts <- readMM('real_data/hgmm_12k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/mm10/matrix.mtx'); dim(mm.counts)
mm.genes <- read.table('real_data/hgmm_12k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/mm10/genes.tsv', 
                       stringsAsFactors = F)
mm.genes <- mm.genes$V1; length(mm.genes)
mm.genes <- substring(mm.genes, 6); length(mm.genes)
rownames(mm.counts) <- mm.genes
colnames(mm.counts) <- mm.barcode

# doublet cells
doublet.barcode <- intersect(hg.barcode, mm.barcode); length(doublet.barcode)

# orthologs genes
mart1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mart2 = useMart("ensembl", dataset="mmusculus_gene_ensembl") 

# human / mouse
og <- getLDS(attributes=c("ensembl_gene_id"),
       filters="ensembl_gene_id", values=mm.genes, mart=mart2,
       attributesL=c("ensembl_gene_id"), martL=mart1); dim(og)
index <- og$Gene.stable.ID.1%in%hg.genes; sum(index)
og <- og[index,]; dim(og)
og <- og[!duplicated(og$Gene.stable.ID),]; dim(og)
og <- og[!duplicated(og$Gene.stable.ID.1),]; dim(og)
mm.og <- og$Gene.stable.ID; length(mm.og)
hg.og <- og$Gene.stable.ID.1; length(hg.og)

hg.counts <- hg.counts[rownames(hg.counts)%in%hg.og,]; dim(hg.counts)
mm.counts <- mm.counts[rownames(mm.counts)%in%mm.og,]; dim(mm.counts)

rownames(hg.counts)[1:10]
rownames(mm.counts)[1:10]
# match gene names from mm to hg
mm_hg <- sapply(rownames(mm.counts), function(x){
  return(og[which(x == og$Gene.stable.ID),]$Gene.stable.ID.1)
}); length(mm_hg)
names(mm_hg) <- NULL
rownames(mm.counts) <- mm_hg

# order both by gene names
hg.counts <- hg.counts[order(rownames(hg.counts)),]
rownames(hg.counts)[1:10]
mm.counts <- mm.counts[order(rownames(mm.counts)),]
rownames(mm.counts)[1:10]
sum(rownames(hg.counts)!=rownames(mm.counts))
#View(as.matrix(hg.counts))
#View(as.matrix(mm.counts))
hg.doublet <- hg.counts[, doublet.barcode]; dim(hg.doublet)
mm.doublet <- mm.counts[, doublet.barcode]; dim(mm.doublet)
doublet.counts <- hg.doublet + mm.doublet; dim(doublet.counts)
#View(as.matrix(doublet.counts))
hg.singlet <- hg.counts[,!colnames(hg.counts)%in%doublet.barcode]; dim(hg.singlet)
mm.singlet <- mm.counts[,!colnames(mm.counts)%in%doublet.barcode]; dim(mm.singlet)
hg_mm.counts <- cbind(hg.singlet, mm.singlet, doublet.counts); dim(hg_mm.counts)
hg_mm.seurat<- CreateSeuratObject(counts = hg_mm.counts, min.cell = 1); hg_mm.seurat
hg_mm.counts <- hg_mm.seurat[["RNA"]]@counts; dim(hg_mm.counts)
label <- c(rep('singlet', dim(hg.singlet)[2]+dim(mm.singlet)[2]), rep('doublet', dim(doublet.counts)[2])); table(label)
summary(colSums(hg_mm.counts))
summary(colSums(hg_mm.counts != 0))

saveRDS(list(hg_mm.counts, label), 'real_data/hg_mm_6k.rds')

############################################################################
# multiplex batch2 control
############################################################################
multiplex.barcode <- read.table('real_data/GSE96583_RAW/GSM2560247_barcodes.tsv', stringsAsFactors = F)
multiplex.barcode <- multiplex.barcode$V1; length(multiplex.barcode)
multiplex.gene <- read.table('real_data/GSE96583_RAW/GSE96583_batch1.genes.tsv', stringsAsFactors = F)
multiplex.gene <- multiplex.gene$V1; length(multiplex.gene)
multiplex.counts <- readMM('real_data/GSE96583_RAW/GSM2560249_2.2.mtx'); dim(multiplex.counts)
colnames(multiplex.counts) <- multiplex.barcode
rownames(multiplex.counts) <- multiplex.gene

label.table <- read.table('real_data/GSE96583_RAW/GSE96583_batch1.total.tsne.df.tsv',header = T, stringsAsFactors = F, 
                          sep = '\t'); dim(label.table)
label.table <- label.table[label.table$barcode%in%multiplex.barcode,]; dim(label.table)
label.table <- label.table[label.table$multiplets!='ambs',]; dim(label.table)
table(label.table$multiplets)
multiplex.counts <- multiplex.counts[,colnames(multiplex.counts)%in%label.table$barcode]; dim(multiplex.counts)
multiplex.counts <- multiplex.counts[,order(colnames(multiplex.counts))]; dim(multiplex.counts)
label.table <- label.table[order(label.table$barcode),]; dim(label.table)
sum(colnames(multiplex.counts) != label.table$barcode)

multiplex.seurat<- CreateSeuratObject(counts = multiplex.counts, min.cell = 1); multiplex.seurat
multiplex.counts <- multiplex.seurat[["RNA"]]@counts; dim(multiplex.counts); class(multiplex.counts)
summary(colSums(multiplex.counts))

saveRDS(list(multiplex.counts, label.table$multiplets), 'real_data/demuxlet_2_stimulated.rds')

############################################################################
# multiplex batch1
############################################################################
# transfer triplet to sparse matrix
triplet <- read.table('real_data/GSE96583_RAW/GSM2560247_C.txt', stringsAsFactors = F)
i <- as.numeric(triplet$V1[-1])
j <- as.numeric(triplet$V2[-1])
x <- as.numeric(triplet$V3[-1])
multiplex.counts <- sparseMatrix(i=i, j=j, x=x, dims = c(triplet$V1[1], triplet$V2[1])); dim(multiplex.counts)

multiplex.barcode <- read.table('real_data/GSE96583_RAW/GSM2560247_barcodes.tsv', stringsAsFactors = F)
multiplex.barcode <- multiplex.barcode$V1; length(multiplex.barcode)
multiplex.gene <- read.table('real_data/GSE96583_RAW/GSE96583_batch1.genes.tsv', stringsAsFactors = F)
multiplex.gene <- multiplex.gene$V1; length(multiplex.gene)
colnames(multiplex.counts) <- multiplex.barcode
rownames(multiplex.counts) <- multiplex.gene

label.table <- read.table('real_data/GSE96583_RAW/GSE96583_batch1.total.tsne.df.tsv',header = T, stringsAsFactors = F, 
                          sep = '\t', row.names = NULL); dim(label.table)
colnames(label.table)[1] <- 'barcode'
label.table <- label.table[label.table$barcode%in%multiplex.barcode,]; dim(label.table)
label.table <- label.table[label.table$multiplets!='ambs',]; dim(label.table)
table(label.table$multiplets)
multiplex.counts <- multiplex.counts[,colnames(multiplex.counts)%in%label.table$barcode]; dim(multiplex.counts)
multiplex.counts <- multiplex.counts[,order(colnames(multiplex.counts))]; dim(multiplex.counts)
label.table <- label.table[order(label.table$barcode),]; dim(label.table)
sum(colnames(multiplex.counts) != label.table$barcode)

multiplex.seurat<- CreateSeuratObject(counts = multiplex.counts, min.cell = 1); multiplex.seurat
multiplex.counts <- multiplex.seurat[["RNA"]]@counts; dim(multiplex.counts); class(multiplex.counts)
summary(colSums(multiplex.counts))

saveRDS(list(multiplex.counts, label.table$multiplets), 'real_data/demuxlet_1_C.rds')

############################################################################
# multiplex jurkat 293t
############################################################################
label.table <- read.table('real_data/jurkat_293t_demuxlet.txt',header = T, stringsAsFactors = F, 
                                     sep = '\t', row.names = NULL); dim(label.table)
label.table$BEST <- sapply(strsplit(label.table$BEST, '-'), function(x){
  return(x[[1]])
})
label.table$BEST <- ifelse(label.table$BEST == 'DBL', 'doublet', 'singlet'); table(label.table$BEST)
multiplex.barcode <- read.table('real_data/jurkat_293t_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/barcodes.tsv', 
                                stringsAsFactors = F)
multiplex.barcode <- multiplex.barcode$V1; length(multiplex.barcode)
#multiplex.barcode <- intersect(label.table$BARCODE, multiplex.barcode); length(multiplex.barcode)
multiplex.gene <- read.table('real_data/jurkat_293t_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/genes.tsv', 
                                stringsAsFactors = F)
multiplex.gene <- multiplex.gene$V1; length(multiplex.gene)
multiplex.counts <- readMM('real_data/jurkat_293t_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/matrix.mtx'); dim(multiplex.counts)
rownames(multiplex.counts) <- multiplex.gene
colnames(multiplex.counts) <- multiplex.barcode

multiplex.counts <- multiplex.counts[,colnames(multiplex.counts)%in%label.table$BARCODE]; dim(multiplex.counts)
multiplex.counts <- multiplex.counts[,order(colnames(multiplex.counts))]; dim(multiplex.counts)
label.table <- label.table[order(label.table$BARCODE),]; dim(label.table)
sum(colnames(multiplex.counts) != label.table$BARCODE)

multiplex.seurat<- CreateSeuratObject(counts = multiplex.counts, min.cell = 1); multiplex.seurat
multiplex.counts <- multiplex.seurat[["RNA"]]@counts; dim(multiplex.counts); class(multiplex.counts)
summary(colSums(multiplex.counts))
table(label.table$BEST)

saveRDS(list(multiplex.counts, label.table$BEST), 'real_data/demuxlet_jurkat_293t.rds')

############################################################################
# kidney
############################################################################
code <- readMM('real_data/solo/GSE140262_RAW/GSM4158566_kidney2_CMO.mtx'); dim(code)


counts <- readMM('real_data/solo/GSE140262_RAW/GSM4158565_kidney2_RNA.mtx'); dim(counts)
counts <- t(counts); dim(counts)
label.table <- read.table('real_data/solo/GSE140262_kidney2metadata.csv',header = T, stringsAsFactors = F, 
                          sep = '\t', row.names = NULL); dim(label.table)
label.table$Category_with_clustering <- ifelse(label.table$Category_with_clustering != 'Doublet', 'singlet', 'doublet')
table(label.table$Category_with_clustering)
genes <- read.table('real_data/solo/GSE140262_kidney2_gene_data.csv',header = T, stringsAsFactors = F, 
                    sep = '\t', row.names = NULL); dim(genes)
rownames(counts) <- genes$gene_id
colnames(counts) <- label.table$index

counts <- counts[,order(colnames(counts))]; dim(counts)
label.table <- label.table[order(label.table$index),]; dim(label.table)
sum(colnames(counts) != label.table$index)

kidney.seurat<- CreateSeuratObject(counts = counts, min.cell = 1); kidney.seurat
counts <- kidney.seurat[["RNA"]]@counts; dim(counts); class(counts)
summary(colSums(counts))
table(label.table$Category_with_clustering)
summary(colSums(counts != 0))

saveRDS(list(counts, label.table$Category_with_clustering), 'real_data/mouse_kidney_2.rds')

############################################################################
# pdx
############################################################################
bar.table <- read.table('real_data/pdx/GSE129578_processed_data_files.csv/PDX_MULTI_matrix.csv',header = T, 
                            stringsAsFactors = F, sep = '\t', row.names = NULL); dim(bar.table)
rownames(bar.table) <- bar.table$CellID
bar.table <- bar.table[,-1]; dim(bar.table)
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

## Identify ideal inter-maxima quantile to set barcode-specific thresholds
threshold.results1 <- findThresh(call.list=bar.table_sweep.list)

## Finalize round 1 classifications, remove negative cells
round1.calls <- classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema)); length(round1.calls); table(round1.calls)
neg.cells1 <- names(round1.calls)[which(round1.calls == "Negative")]; length(neg.cells1)
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells1), ]; dim(bar.table)

## Round 2 -----------------------------------------------------------------------------------------------------
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results2 <- findThresh(call.list=bar.table_sweep.list)
round2.calls <- classifyCells(bar.table, q=findQ(threshold.results2$res, threshold.results2$extrema)); length(round2.calls); table(round2.calls)
neg.cells2 <- names(round2.calls)[which(round2.calls == "Negative")]; length(neg.cells2)
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells2), ]; dim(bar.table)

## Round 3 -----------------------------------------------------------------------------------------------------
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results3 <- findThresh(call.list=bar.table_sweep.list)
round3.calls <- classifyCells(bar.table, q=findQ(threshold.results3$res, threshold.results3$extrema)); length(round3.calls); table(round3.calls)
neg.cells3 <- names(round3.calls)[which(round3.calls == "Negative")]; length(neg.cells3)
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells3), ]; dim(bar.table)

## Round 4 -----------------------------------------------------------------------------------------------------
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results4 <- findThresh(call.list=bar.table_sweep.list)
round4.calls <- classifyCells(bar.table, q=findQ(threshold.results4$res, threshold.results4$extrema)); length(round4.calls); table(round4.calls)
label <- ifelse(round4.calls != 'Doublet', 'singlet', 'doublet'); table(label); length(label)

# select human cells
barcode.hg <- read.table('real_data/pdx/GSE129578_processed_data_files.tsv/PDX_hs_barcodes.tsv',header = F,
                        stringsAsFactors = F, sep = '\t', row.names = NULL); dim(barcode.hg)
barcode.hg <- barcode.hg$V1; length(barcode.hg)
barcode.hg <- sapply(strsplit(barcode.hg, '-'), function(x){
  return(x[[1]])
}); length(barcode.hg)
hg.genes <- read.table('real_data/pdx/GSE129578_processed_data_files.tsv/PDX_hs_genes.tsv',header = F,
                        stringsAsFactors = F, sep = '\t', row.names = NULL); dim(hg.genes)
hg.genes <- hg.genes$V1; length(hg.genes)
hg.genes <- substring(hg.genes, 6); length(hg.genes)
counts.hg <- readMM('real_data/pdx/GSE129578_PDX_hs_matrix.mtx'); dim(counts.hg)
rownames(counts.hg) <- hg.genes
colnames(counts.hg) <- barcode.hg
barcode.hg <- intersect(barcode.hg, names(label)); length(barcode.hg)
counts.hg <- counts.hg[,which(colnames(counts.hg)%in%barcode.hg)]; dim(counts.hg)
#summary(colSums(counts.hg))

# select mouse cells
barcode.mm <- read.table('real_data/pdx/GSE129578_processed_data_files.tsv/PDX_mm_barcodes.tsv',header = F, 
                         stringsAsFactors = F, sep = '\t', row.names = NULL); dim(barcode.mm)
barcode.mm <- barcode.mm$V1; length(barcode.mm)
barcode.mm <- sapply(strsplit(barcode.mm, '-'), function(x){
  return(x[[1]])
}); length(barcode.mm)
mm.genes <- read.table('real_data/pdx/GSE129578_processed_data_files.tsv/PDX_mm_genes.tsv',header = F, 
                       stringsAsFactors = F, sep = '\t', row.names = NULL); dim(mm.genes)
mm.genes <- mm.genes$V1; length(mm.genes)
mm.genes <- substring(mm.genes, 6); length(mm.genes)
counts.mm <- readMM('real_data/pdx/GSE129578_PDX_mm_matrix.mtx'); dim(counts.mm)
rownames(counts.mm) <- mm.genes
colnames(counts.mm) <- barcode.mm
barcode.mm <- intersect(barcode.mm, names(label)); length(barcode.mm)
counts.mm <- counts.mm[,which(colnames(counts.mm)%in%barcode.mm)]; dim(counts.mm)
counts.mm <- counts.mm[,order(colnames(counts.mm))]; dim(counts.mm)
label <- label[order(names(label))]; length(label)
sum(colnames(counts.mm) != names(label))

# transfer human genes to mouse genes 
mart1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mart2 = useMart("ensembl", dataset="mmusculus_gene_ensembl") 

# human / mouse
og <- getLDS(attributes=c("ensembl_gene_id"),
             filters="ensembl_gene_id", values=mm.genes, mart=mart2,
             attributesL=c("ensembl_gene_id"), martL=mart1); dim(og)
index <- og$Gene.stable.ID.1%in%hg.genes; sum(index)
og <- og[index,]; dim(og)
og <- og[!duplicated(og$Gene.stable.ID),]; dim(og)
og <- og[!duplicated(og$Gene.stable.ID.1),]; dim(og)
mm.og <- og$Gene.stable.ID; length(mm.og)
hg.og <- og$Gene.stable.ID.1; length(hg.og)
counts.hg <- counts.hg[rownames(counts.hg)%in%hg.og,]; dim(counts.hg)
counts.mm <- counts.mm[rownames(counts.mm)%in%mm.og,]; dim(counts.mm)

rownames(counts.hg)[1:10]
rownames(counts.mm)[1:10]
# match gene names from mm to hg
mm_hg <- sapply(rownames(counts.mm), function(x){
  return(og[which(x == og$Gene.stable.ID),]$Gene.stable.ID.1)
}); length(mm_hg)
names(mm_hg) <- NULL
rownames(counts.mm) <- mm_hg

# order both by gene names
counts.hg <- counts.hg[order(rownames(counts.hg)),]
rownames(counts.hg)[1:10]
counts.mm <- counts.mm[order(rownames(counts.mm)),]
rownames(counts.mm)[1:10]
sum(rownames(counts.hg)!=rownames(counts.mm))
sum(colnames(counts.hg)!=colnames(counts.mm))
hg_mm.counts <- counts.hg + counts.mm; dim(hg_mm.counts)
sum(colnames(hg_mm.counts) != names(label))

hg_mm.seurat<- CreateSeuratObject(counts = hg_mm.counts, min.cell = 1); hg_mm.seurat
hg_mm.counts <- hg_mm.seurat[["RNA"]]@counts; dim(hg_mm.counts); class(hg_mm.counts)
summary(colSums(hg_mm.counts))
summary(colSums(hg_mm.counts != 0))
table(label)

# only save mouse data
saveRDS(list(hg_mm.counts, label), 'real_data/pdx.rds')

############################################################################
# HMEC orig
############################################################################
bar.table <- read.table('real_data/pdx/GSE129578_processed_data_files.csv/HMEC_orig_MULTI_matrix.csv',header = T, 
                        stringsAsFactors = F, sep = '\t', row.names = NULL); dim(bar.table)
rownames(bar.table) <- bar.table$CellID
bar.table <- bar.table[,-1]; dim(bar.table)
# remove summary column
bar.table <- bar.table[,-dim(bar.table)[2]]; dim(bar.table)
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
## Identify ideal inter-maxima quantile to set barcode-specific thresholds
threshold.results1 <- findThresh(call.list=bar.table_sweep.list)
## Finalize round 1 classifications, remove negative cells
round1.calls <- classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema)); length(round1.calls); table(round1.calls)
neg.cells1 <- names(round1.calls)[which(round1.calls == "Negative")]; length(neg.cells1)
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells1), ]; dim(bar.table)

## Round 2 -----------------------------------------------------------------------------------------------------
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results2 <- findThresh(call.list=bar.table_sweep.list)
round2.calls <- classifyCells(bar.table, q=findQ(threshold.results2$res, threshold.results2$extrema)); length(round2.calls); table(round2.calls)
neg.cells2 <- names(round2.calls)[which(round2.calls == "Negative")]; length(neg.cells2)
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells2), ]; dim(bar.table)

## Round 3 -----------------------------------------------------------------------------------------------------
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results3 <- findThresh(call.list=bar.table_sweep.list)
round3.calls <- classifyCells(bar.table, q=findQ(threshold.results3$res, threshold.results3$extrema)); length(round3.calls); table(round3.calls)
neg.cells3 <- names(round3.calls)[which(round3.calls == "Negative")]; length(neg.cells3)
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells3), ]; dim(bar.table)

## Round 4 -----------------------------------------------------------------------------------------------------
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results4 <- findThresh(call.list=bar.table_sweep.list)
round4.calls <- classifyCells(bar.table, q=findQ(threshold.results4$res, threshold.results4$extrema)); length(round4.calls); table(round4.calls)
label <- ifelse(round4.calls != 'Doublet', 'singlet', 'doublet'); table(label); length(label)

# select cells
barcode.HMEC.orig <- read.table('real_data/pdx/GSE129578_processed_data_files.tsv/HMEC_orig_barcodes.tsv',header = F, 
                         stringsAsFactors = F, sep = '\t', row.names = NULL); dim(barcode.HMEC.orig)
barcode.HMEC.orig <- barcode.HMEC.orig$V1; length(barcode.HMEC.orig)
HMEC.orig.genes <- read.table('real_data/pdx/GSE129578_processed_data_files.tsv/HMEC_orig_genes.tsv',header = F, 
                       stringsAsFactors = F, sep = '\t', row.names = NULL); dim(HMEC.orig.genes)
HMEC.orig.genes <- HMEC.orig.genes$V1; length(HMEC.orig.genes)
counts.HMEC.orig <- readMM('real_data/pdx/GSE129578_HMEC_orig_matrix.mtx'); dim(counts.HMEC.orig)
rownames(counts.HMEC.orig) <- HMEC.orig.genes
colnames(counts.HMEC.orig) <- barcode.HMEC.orig
barcode.HMEC.orig <- intersect(barcode.HMEC.orig, names(label)); length(barcode.HMEC.orig)
counts.HMEC.orig <- counts.HMEC.orig[,which(colnames(counts.HMEC.orig)%in%barcode.HMEC.orig)]; dim(counts.HMEC.orig)
counts.HMEC.orig <- counts.HMEC.orig[,order(colnames(counts.HMEC.orig))]; dim(counts.HMEC.orig)
label <- label[order(names(label))]; length(label)
sum(colnames(counts.HMEC.orig) != names(label))
HMEC.orig.seurat<- CreateSeuratObject(counts = counts.HMEC.orig, min.cell = 1); HMEC.orig.seurat
counts.HMEC.orig <- HMEC.orig.seurat[["RNA"]]@counts; dim(counts.HMEC.orig); class(counts.HMEC.orig)
summary(colSums(counts.HMEC.orig))
table(label)

# only save mouse data
saveRDS(list(counts.HMEC.orig, label), 'real_data/HMEC_orig.rds')

############################################################################
# HMEC rep
############################################################################
bar.table <- read.table('real_data/pdx/GSE129578_processed_data_files.csv/HMEC_techrep_MULTI_matrix.csv',header = T, 
                        stringsAsFactors = F, sep = '\t', row.names = NULL); dim(bar.table)
rownames(bar.table) <- bar.table$CellID
bar.table <- bar.table[,-1]; dim(bar.table)
# remove summary column
bar.table <- bar.table[,-dim(bar.table)[2]]; dim(bar.table)
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
## Identify ideal inter-maxima quantile to set barcode-specific thresholds
threshold.results1 <- findThresh(call.list=bar.table_sweep.list)
## Finalize round 1 classifications, remove negative cells
round1.calls <- classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema)); length(round1.calls); table(round1.calls)
neg.cells1 <- names(round1.calls)[which(round1.calls == "Negative")]; length(neg.cells1)
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells1), ]; dim(bar.table)

## Round 2 -----------------------------------------------------------------------------------------------------
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results2 <- findThresh(call.list=bar.table_sweep.list)
round2.calls <- classifyCells(bar.table, q=findQ(threshold.results2$res, threshold.results2$extrema)); length(round2.calls); table(round2.calls)
neg.cells2 <- names(round2.calls)[which(round2.calls == "Negative")]; length(neg.cells2)
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells2), ]; dim(bar.table)

## Round 3 -----------------------------------------------------------------------------------------------------
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results3 <- findThresh(call.list=bar.table_sweep.list)
round3.calls <- classifyCells(bar.table, q=findQ(threshold.results3$res, threshold.results3$extrema)); length(round3.calls); table(round3.calls)
neg.cells3 <- names(round3.calls)[which(round3.calls == "Negative")]; length(neg.cells3)
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells3), ]; dim(bar.table)

## Round 4 -----------------------------------------------------------------------------------------------------
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results4 <- findThresh(call.list=bar.table_sweep.list)
round4.calls <- classifyCells(bar.table, q=findQ(threshold.results4$res, threshold.results4$extrema)); length(round4.calls); table(round4.calls)
label <- ifelse(round4.calls != 'Doublet', 'singlet', 'doublet'); table(label); length(label)

# select cells
barcode.HMEC.rep <- read.table('real_data/pdx/GSE129578_processed_data_files.tsv/HMEC1_techrep_barcodes.tsv',header = F, 
                                stringsAsFactors = F, sep = '\t', row.names = NULL); dim(barcode.HMEC.rep)
barcode.HMEC.rep <- barcode.HMEC.rep$V1; length(barcode.HMEC.rep)
HMEC.rep.genes <- read.table('real_data/pdx/GSE129578_processed_data_files.tsv/HMEC1_techrep_genes.tsv',header = F, 
                              stringsAsFactors = F, sep = '\t', row.names = NULL); dim(HMEC.rep.genes)
HMEC.rep.genes <- HMEC.rep.genes$V1; length(HMEC.rep.genes)
counts.HMEC.rep1 <- readMM('real_data/pdx/GSE129578_HMEC1_techrep_matrix.mtx'); dim(counts.HMEC.rep1)
counts.HMEC.rep2 <- readMM('real_data/pdx/GSE129578_HMEC2_techrep_matrix.mtx'); dim(counts.HMEC.rep2)
counts.HMEC.rep3 <- readMM('real_data/pdx/GSE129578_HMEC3_techrep_matrix.mtx'); dim(counts.HMEC.rep3)
counts.HMEC.rep4 <- readMM('real_data/pdx/GSE129578_HMEC4_techrep_matrix.mtx'); dim(counts.HMEC.rep4)
counts.HMEC.rep <- counts.HMEC.rep1 + counts.HMEC.rep2 + counts.HMEC.rep3 + counts.HMEC.rep4; dim(counts.HMEC.rep)

rownames(counts.HMEC.rep) <- HMEC.rep.genes
colnames(counts.HMEC.rep) <- barcode.HMEC.rep
barcode.HMEC.rep <- intersect(barcode.HMEC.rep, names(label)); length(barcode.HMEC.rep)
counts.HMEC.rep <- counts.HMEC.rep[,which(colnames(counts.HMEC.rep)%in%barcode.HMEC.rep)]; dim(counts.HMEC.rep)
counts.HMEC.rep <- counts.HMEC.rep[,order(colnames(counts.HMEC.rep))]; dim(counts.HMEC.rep)
label <- label[names(label)%in%barcode.HMEC.rep]; length(label); table(label)
label <- label[order(names(label))]; length(label)
sum(colnames(counts.HMEC.rep) != names(label))
HMEC.rep.seurat<- CreateSeuratObject(counts = counts.HMEC.rep, min.cell = 1); HMEC.rep.seurat
counts.HMEC.rep <- HMEC.rep.seurat[["RNA"]]@counts; dim(counts.HMEC.rep); class(counts.HMEC.rep)
summary(colSums(counts.HMEC.rep))
table(label)

# only save mouse data
saveRDS(list(counts.HMEC.rep, label), 'real_data/HMEC_rep.rds')

############################################################################
# HMEC poc
############################################################################
bar.table <- read.table('real_data/pdx/GSE129578_processed_data_files.csv/POC_MULTI_matrix.csv',header = T, 
                        stringsAsFactors = F, sep = '\t', row.names = NULL); dim(bar.table)
bar.table <- bar.table[bar.table$LaneID != 'UN',]; dim(bar.table)
rownames(bar.table) <- bar.table$CellID
bar.table <- bar.table[,-c(1,2,6)]; dim(bar.table)
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
## Identify ideal inter-maxima quantile to set barcode-specific thresholds
threshold.results1 <- findThresh(call.list=bar.table_sweep.list)
## Finalize round 1 classifications, remove negative cells
round1.calls <- classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema)); length(round1.calls); table(round1.calls)
neg.cells1 <- names(round1.calls)[which(round1.calls == "Negative")]; length(neg.cells1)
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells1), ]; dim(bar.table)

## Round 2 -----------------------------------------------------------------------------------------------------
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results2 <- findThresh(call.list=bar.table_sweep.list)
round2.calls <- classifyCells(bar.table, q=findQ(threshold.results2$res, threshold.results2$extrema)); length(round2.calls); table(round2.calls)
label <- ifelse(round2.calls != 'Doublet', 'singlet', 'doublet'); table(label); length(label)

# select cells
barcode.poc <- read.table('real_data/pdx/GSE129578_processed_data_files.tsv/POC_barcodes.tsv',header = F, 
                               stringsAsFactors = F, sep = '\t', row.names = NULL); dim(barcode.poc)
barcode.poc <- barcode.poc$V1; length(barcode.poc)
poc.genes <- read.table('real_data/pdx/GSE129578_processed_data_files.tsv/POC_genes.tsv',header = F, 
                             stringsAsFactors = F, sep = '\t', row.names = NULL); dim(poc.genes)
poc.genes <- poc.genes$V1; length(poc.genes)
counts.poc <- readMM('real_data/pdx/GSE129578_POC_matrix.mtx'); dim(counts.poc)

rownames(counts.poc) <- poc.genes
colnames(counts.poc) <- barcode.poc
barcode.poc <- intersect(barcode.poc, names(label)); length(barcode.poc)
counts.poc <- counts.poc[,which(colnames(counts.poc)%in%barcode.poc)]; dim(counts.poc)
counts.poc <- counts.poc[,order(colnames(counts.poc))]; dim(counts.poc)
label <- label[order(names(label))]; length(label)
sum(colnames(counts.poc) != names(label))
poc.seurat<- CreateSeuratObject(counts = counts.poc, min.cell = 1); poc.seurat
counts.poc <- poc.seurat[["RNA"]]@counts; dim(counts.poc); class(counts.poc)
summary(colSums(counts.poc))
table(label)

# only save mouse data
saveRDS(list(counts.poc, label), 'real_data/poc.rds')

############################################################################
# HMEC poc nuc
############################################################################
bar.table <- read.table('real_data/pdx/GSE129578_processed_data_files.csv/POC_nuc_MULTI_matrix.csv',header = T, 
                        stringsAsFactors = F, sep = '\t', row.names = NULL); dim(bar.table)
bar.table <- bar.table[bar.table$nUMI_Bar != 0,]; dim(bar.table)
bar.table <- bar.table[!duplicated(bar.table$CellID),]; dim(bar.table)
rownames(bar.table) <- bar.table$CellID
bar.table <- bar.table[,-c(1,14)]; dim(bar.table)
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
## Identify ideal inter-maxima quantile to set barcode-specific thresholds
threshold.results1 <- findThresh(call.list=bar.table_sweep.list)
## Finalize round 1 classifications, remove negative cells
round1.calls <- classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema)); length(round1.calls); table(round1.calls)
neg.cells1 <- names(round1.calls)[which(round1.calls == "Negative")]; length(neg.cells1)
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells1), ]; dim(bar.table)

## Round 2 -----------------------------------------------------------------------------------------------------
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results2 <- findThresh(call.list=bar.table_sweep.list)
round2.calls <- classifyCells(bar.table, q=findQ(threshold.results2$res, threshold.results2$extrema)); length(round2.calls); table(round2.calls)
neg.cells2 <- names(round2.calls)[which(round2.calls == "Negative")]; length(neg.cells2)
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells2), ]; dim(bar.table)

## Round 3 -----------------------------------------------------------------------------------------------------
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results3 <- findThresh(call.list=bar.table_sweep.list)
round3.calls <- classifyCells(bar.table, q=findQ(threshold.results3$res, threshold.results3$extrema)); length(round3.calls); table(round3.calls)
label <- ifelse(round3.calls != 'Doublet', 'singlet', 'doublet'); table(label); length(label)

# select cells
barcode.poc.nuc <- read.table('real_data/pdx/GSE129578_processed_data_files.tsv/POC_nuc_barcodes.tsv',header = F, 
                          stringsAsFactors = F, sep = '\t', row.names = NULL); dim(barcode.poc.nuc)
barcode.poc.nuc <- barcode.poc.nuc$V1; length(barcode.poc.nuc)
poc.nuc.genes <- read.table('real_data/pdx/GSE129578_processed_data_files.tsv/POC_nuc_genes.tsv',header = F, 
                        stringsAsFactors = F, sep = '\t', row.names = NULL); dim(poc.nuc.genes)
poc.nuc.genes <- poc.nuc.genes$V1; length(poc.nuc.genes)
counts.poc.nuc <- readMM('real_data/pdx/GSE129578_POC_nuc_matrix.mtx'); dim(counts.poc.nuc)

rownames(counts.poc.nuc) <- poc.nuc.genes
colnames(counts.poc.nuc) <- barcode.poc.nuc
barcode.poc.nuc <- intersect(barcode.poc.nuc, names(label)); length(barcode.poc.nuc)
counts.poc.nuc <- counts.poc.nuc[,which(colnames(counts.poc.nuc)%in%barcode.poc.nuc)]; dim(counts.poc.nuc)
counts.poc.nuc <- counts.poc.nuc[,order(colnames(counts.poc.nuc))]; dim(counts.poc.nuc)
label <- label[order(names(label))]; length(label)
sum(colnames(counts.poc.nuc) != names(label))
poc.nuc.seurat<- CreateSeuratObject(counts = counts.poc.nuc, min.cell = 1); poc.nuc.seurat
counts.poc.nuc <- poc.nuc.seurat[["RNA"]]@counts; dim(counts.poc.nuc); class(counts.poc.nuc)

# unify genes
index <- sapply(rownames(counts.poc.nuc), function(x){
  return(grepl('hg19', x, fixed = T))
},simplify = T); table(index)
counts.hg <- counts.poc.nuc[index,]; dim(counts.hg)
counts.mm <- counts.poc.nuc[!index,]; dim(counts.mm)
summary(colSums(counts.hg))
summary(colSums(counts.mm))
hg.genes <- substring(rownames(counts.hg), 6); length(hg.genes)
mm.genes <- substring(rownames(counts.mm), 6); length(mm.genes)
rownames(counts.hg) <- hg.genes
rownames(counts.mm) <- mm.genes

# orthologs genes
mart1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mart2 = useMart("ensembl", dataset="mmusculus_gene_ensembl") 

# human / mouse
og <- getLDS(attributes=c("ensembl_gene_id"),
             filters="ensembl_gene_id", values=mm.genes, mart=mart2,
             attributesL=c("ensembl_gene_id"), martL=mart1); dim(og)
index <- og$Gene.stable.ID.1%in%hg.genes; sum(index)
og <- og[index,]; dim(og)
og <- og[!duplicated(og$Gene.stable.ID),]; dim(og)
og <- og[!duplicated(og$Gene.stable.ID.1),]; dim(og)
og.mm <- og$Gene.stable.ID
og.hg <- og$Gene.stable.ID.1
counts.mm <- counts.mm[rownames(counts.mm)%in%og.mm,]; dim(counts.mm)
mm_hg <- sapply(rownames(counts.mm), function(x){
  return(og[which(x == og$Gene.stable.ID),]$Gene.stable.ID.1)
}); length(mm_hg)
names(mm_hg) <- NULL
rownames(counts.mm) <- mm_hg
mrow <- match(rownames(counts.mm),rownames(counts.hg))
m <- counts.hg[mrow,] + counts.mm; dim(m)
m_ <- counts.hg[-mrow,]; dim(m_)
counts <- rbind(m, m_); dim(counts)
counts <- counts[order(rownames(counts)),]
sum(colnames(counts) != names(label))
table(label)
summary(colSums(counts))
summary(colSums(counts!=0))

# only save mouse data
saveRDS(list(counts, label), 'real_data/nuc-MULTI.rds')


######################################################################
# reread everything
######################################################################
# read data; parameter setting
#loc <- 'real_data/pbmc-ch.rds'
#loc <- 'real_data/cline-ch.rds'
#loc <- 'real_data/mkidney-ch.rds'
#loc <- 'real_data/hm-12k.rds'
#loc <- 'real_data/hm-6k.rds'
#loc <- 'real_data/pbmc-1A-dm.rds'
#loc <- 'real_data/pbmc-1B-dm.rds'
#loc <- 'real_data/pbmc-1C-dm.rds'
#loc <- 'real_data/pbmc-2ctrl-dm.rds'
#loc <- 'real_data/pbmc-2stim-dm.rds'
#loc <- 'real_data/J293t-dm.rds'
#loc <- 'real_data/pdx-MULTI.rds'
#loc <- 'real_data/HMEC-orig-MULTI.rds'
#loc <- 'real_data/HMEC-rep-MULTI.rds'
#loc <- 'real_data/HEK-HMEC-MULTI.rds'
loc <- 'real_data/nuc-MULTI.rds'

data <- readRDS(loc)
counts <- data[[1]]; dim(counts)
label <- data[[2]]; table(label); length(label)
#write.csv(label, file ="label.csv", row.names=FALSE,)
#label <- ifelse(label == 'doublet', 1, 0); table(label)
# save to loom
library(loomR)
create(paste(paste('real_data',strsplit(strsplit(loc,'/')[[1]][2], '\\.')[[1]][1], sep='/'), 'loom', sep='.'), 
       counts, transpose = T)







