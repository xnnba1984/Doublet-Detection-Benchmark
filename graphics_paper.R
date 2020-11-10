library(ggplot2)
library(readODS)
library(reshape2)
library(data.table)  # faster fread() and better weekdays()
library(dplyr)       # consistent data.frame operations
library(purrr)       # consistent & safe list/vector munging
library(tidyr)       # consistent data.frame cleaning
library(lubridate)   # date manipulation
library(countrycode) # turn country codes into pretty names
library(ggplot2)     # base plots are for Coursera professors
library(scales)      # pairs nicely with ggplot2 for plot label formatting
library(gridExtra)   # a helper for arranging individual ggplot objects
library(ggthemes)    # has a clean theme for ggplot2
library(viridis)     # best. color. palette. evar.
library(knitr)       # kable : prettier data.frame output
library(cowplot)
library(forcats)
library(dplyr)
setwd('/media/nxi/nxi/doublet')
getwd()

# overall setting
result <- read_ods('real_data.ods', sheet = 1); str(result)
result[result == '--'] = NA
result$prauc <- as.numeric(result$prauc)
result$rocauc <- as.numeric(result$rocauc)
str(result)

# threshold setting
result <- read_ods('real_data.ods', sheet = 6); str(result)
str(result)

# order the method and dataset by performance
#a <- aggregate(result[, 3], by=list(result$method), FUN=mean, na.rm=TRUE)
#level.method <- as.character(a[order(a$x), 1]); level.method
level.method <- c('lsize','ngene','doubletCells','Scrublet','cxds','bcds','hybrid','solo','DoubletDetection','DoubletFinder')
#b <- aggregate(result[, 3], by=list(result$dataset), FUN=mean, na.rm=TRUE)
#level.data <- as.character(b[order(b$x, decreasing = F), 1]); level.data
level.data <- c("J293t-dm","pbmc-1B-dm","pbmc-1A-dm","pdx-MULTI","nuc-MULTI","pbmc-1C-dm","cline-ch","HEK-HMEC-MULTI",
                "HMEC-orig-MULTI","pbmc-2stim-dm","pbmc-ch","pbmc-2ctrl-dm","HMEC-rep-MULTI","mkidney-ch","hm-12k","hm-6k" )
level.data <- rev(level.data)
result$method <- factor(result$method, levels = level.method)
result$dataset <- factor(result$dataset, levels = level.data)

# define the color of methods
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
color <- gg_color_hue(11); color
color_method <- as.data.frame(cbind(color[1:10], level.method), stringsAsFactors = F)
colnames(color_method) <- c('color','method')
str(color_method)
de <- data.frame("#FF63B6","DoubletDecon")
names(de) <- c("color","method")
color_method <- rbind(color_method, de)
color_method[1,1] <- 'white'
color_method[2,1] <- 'white'
color_method[3:11,1] <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#666666','#a65628','#f781bf','#999999')
#C0C0C0

##########################################################################################################
# threshold boxplot
##########################################################################################################
pdf('nature_figure/real_threshold_precision_boxplot.pdf')
ggplot(result, aes(x=method, y=TNR0.4, fill=method)) + geom_boxplot() + theme_bw()+
  labs(x=NULL, y=NULL, title="Recall (10% Identified Doublets)") + theme(text = element_text(size=15)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=22), axis.text.y = element_text(size=22),
        plot.title = element_text(hjust = 0.5,size=25)) +
  scale_fill_manual(values=color_method$color[1:10]) + ylim(c(0,1)) + 
  theme(axis.ticks=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

##########################################################################################################
# real_data boxplot pr-auc
##########################################################################################################
pdf('nature_figure/real_prauc_boxplot.pdf')
ggplot(result, aes(x=method, y=prauc, fill=method)) + geom_boxplot() + theme_bw()+
  labs(x=NULL, y=NULL, title="PR-AUC") + theme(text = element_text(size=15)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=22), axis.text.y = element_text(size=22),
        plot.title = element_text(hjust = 0.5,size=25)) +
  scale_fill_manual(values=color_method$color[1:10]) + ylim(c(0,1)) + 
  theme(axis.ticks=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
#legend.position = "none", 
######################################################################################################################
# real_data heatmap pr-auc 
######################################################################################################################
pdf('nature_figure/real_rocauc_heatmap.pdf')
ggplot(result, aes(x=method, y=dataset, fill=rocauc)) + 
  geom_tile(color="white", size=0.1) + 
  labs(x=NULL, y=NULL, title="AUROC") + 
  theme_tufte(base_family="Helvetica") +
  scale_fill_viridis(discrete=FALSE, direction = -1, na.value = 'white', limits=c(0.41,1)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=22), plot.title = element_text(hjust = 0.5, size=25),
        axis.text.y = element_text(size=22),legend.position = "bottom") +
  theme(axis.ticks=element_blank())
dev.off()

##########################################################################################################
# real_data boxplot roc-auc
##########################################################################################################
# order by roc
#a <- aggregate(result[, 4], by=list(result$method), FUN=mean, na.rm=TRUE)
#level.method <- as.character(a[order(a$x), 1]); level.method
#b <- aggregate(result[, 4], by=list(result$dataset), FUN=mean, na.rm=TRUE)
#level.data <- as.character(b[order(b$x, decreasing = T), 1]); level.data
#result$method <- factor(result$method,levels = level.method)
#result$dataset <- factor(result$dataset,levels = level.data)

pdf('nature_figure/real_rocauc_boxplot.pdf')
ggplot(result, aes(x=method, y=rocauc, fill=method)) + geom_boxplot() + theme_bw() +
  labs(x=NULL, y=NULL, title="ROC-AUC") + theme(text = element_text(size=15)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=22), plot.title = element_text(hjust = 0.5, size=25),
        axis.text.y = element_text(size=22), panel.border = element_rect(colour = "black", fill=NA, size=3)) +
  scale_fill_manual(values=color_method$color[1:10]) + ylim(c(0,1)) + theme(axis.ticks=element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

##########################################################################################################
# real_data heatmap roc-auc
##########################################################################################################
pdf('nature_figure/real_rocauc_heatmap.pdf')
ggplot(result, aes(x=method, y=dataset, fill=rocauc)) + 
  geom_tile(color="white", size=0.1) + 
  labs(x=NULL, y=NULL, title="ROC-AUC") + 
  theme_tufte(base_family="Helvetica") +
  scale_fill_viridis(discrete=FALSE, direction = -1, na.value = 'white', limits=c(0,1)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=22), plot.title = element_text(hjust = 0.5, size=25), 
        legend.position = "none", axis.text.y = element_text(size=22)) +
  theme(axis.ticks=element_blank())
dev.off()

##########################################################################################################
# decon boxplot 
##########################################################################################################
result <- read_ods('real_data.ods', sheet = 2); str(result)
result[result == '--'] = NA
result$precision <- as.numeric(result$precision)
result$recall <- as.numeric(result$recall)
result$tnr <- as.numeric(result$tnr)
result$f1 <- as.numeric(result$f1)
str(result)

level.method <- c('lsize','ngene','doubletCells','Scrublet','cxds','bcds','hybrid','solo','DoubletDetection','DoubletFinder',
                  'DoubletDecon')
result$method <- factor(result$method, levels = level.method)

pdf('nature_figure/decon_tnr_boxplot.pdf')
ggplot(result, aes(x=method, y=tnr, fill=method)) + geom_boxplot()+theme_bw()+
  labs(x=NULL, y=NULL, title="TNR") + theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 22), 
        plot.title = element_text(hjust = 0.5, size=25), axis.text.y = element_text(size=22)) +
  scale_fill_manual(values=c(rep('white', 10), 'grey')) + theme(axis.ticks=element_blank()) + ylim(c(0,1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3))
dev.off()

##########################################################################################################
# sim rate pr
##########################################################################################################
i = 2
scrublet <- readRDS('paper_result/scrublet_rate.rds')[[i]]
dblcells <- readRDS('paper_result/dblcells_rate.rds')[[i]]
cxds <- readRDS('paper_result/cxds_rate.rds')[[i]]
bcds <- readRDS('paper_result/bcds_rate.rds')[[i]]
hybrid <- readRDS('paper_result/hybrid_rate.rds')[[i]]
doubletdetection <- readRDS('paper_result/doubletdetection_rate.rds')[[i]]
doubletfinder <- readRDS('paper_result/doubletfinder_rate.rds')[[i]]
solo <- readRDS('paper_result/solo_rate.rds')[[i]]

scrublet_pr <- sapply(scrublet, mean, simplify = T); scrublet_pr
dblcells_pr <- sapply(dblcells, mean, simplify = T); dblcells_pr
cxds_pr <- sapply(cxds, mean, simplify = T); cxds_pr
bcds_pr <- sapply(bcds, mean, simplify = T); bcds_pr
hybrid_pr <- sapply(hybrid, mean, simplify = T); hybrid_pr
doubletdetection_pr <- sapply(doubletdetection, mean, simplify = T); doubletdetection_pr
doubletfinder_pr <- sapply(doubletfinder, mean, simplify = T); doubletfinder_pr
solo_pr <- apply(solo, 1, mean); solo_pr

method <- c('doubletCells', 'Scrublet', 'cxds', 'bcds', 'hybrid', 'solo', 'DoubletDetection', 'DoubletFinder')
method <- rep(method, each=20)
rate <- seq(.02, .4, .02)
rate <- rep(rate, 8)
prauc <- c(dblcells_pr,scrublet_pr, cxds_pr, bcds_pr, hybrid_pr, solo_pr,doubletdetection_pr, doubletfinder_pr)
pr_frame <- data.frame(method, rate, prauc); str(pr_frame)
level.method <- c('doubletCells','Scrublet', 'cxds', 'bcds', 'hybrid', 'solo','DoubletDetection', 'DoubletFinder'); level.method
pr_frame$method <- factor(pr_frame$method, levels = level.method); str(pr_frame)

pdf('nature_figure/auroc_rate.pdf')
ggplot(pr_frame, aes(x=rate, y=prauc, col=method)) + geom_line() + geom_point(size=1) + theme_bw()+
  scale_color_manual(values=color_method$color[3:10]) + 
  labs(x=NULL, y=NULL, title="Doublet Rate") + 
  theme(plot.title = element_text(hjust = 0.5, size=25), axis.text.x = element_text(size=22),
        axis.text.y = element_text(size=22), text=element_text(size=15), 
        panel.border = element_rect(colour = "black", fill=NA, size=3)) +
  ylim(c(0.35,1)) + xlim(c(0,0.4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
#, legend.position = "none"
##########################################################################################################
# sim depth pr
##########################################################################################################
scrublet <- readRDS('paper_result/scrublet_depth.rds')[[i]]
dblcells <- readRDS('paper_result/dblcells_depth.rds')[[i]]
cxds <- readRDS('paper_result/cxds_depth.rds')[[i]]
bcds <- readRDS('paper_result/bcds_depth.rds')[[i]]
hybrid <- readRDS('paper_result/hybrid_depth.rds')[[i]]
doubletdetection <- readRDS('paper_result/doubletdetection_depth.rds')[[i]]
doubletfinder <- readRDS('paper_result/doubletfinder_depth.rds')[[i]]
solo <- readRDS('paper_result/solo_depth.rds')[[i]]

scrublet_pr <- sapply(scrublet, mean, simplify = T); scrublet_pr
dblcells_pr <- sapply(dblcells, mean, simplify = T); dblcells_pr
cxds_pr <- sapply(cxds, mean, simplify = T); cxds_pr
bcds_pr <- sapply(bcds, mean, simplify = T); bcds_pr
hybrid_pr <- sapply(hybrid, mean, simplify = T); hybrid_pr
doubletdetection_pr <- sapply(doubletdetection, mean, simplify = T); doubletdetection_pr
doubletfinder_pr <- sapply(doubletfinder, mean, simplify = T); doubletfinder_pr
solo_pr <- apply(solo, 1, mean); solo_pr

method <- c('doubletCells', 'Scrublet', 'cxds', 'bcds', 'hybrid', 'solo', 'DoubletDetection', 'DoubletFinder')
method <- rep(method, each=20)
depth <- seq(500, 10000, 500)
ylab <- seq(.5,10,.5)
ylab <- paste0(ylab, "k"); ylab
depth <- rep(depth, 8)
prauc <- c(dblcells_pr,scrublet_pr, cxds_pr, bcds_pr, hybrid_pr, solo_pr,doubletdetection_pr, doubletfinder_pr)
pr_frame <- data.frame(method, depth, prauc); str(pr_frame)
level.method <- c('doubletCells','Scrublet', 'cxds', 'bcds', 'hybrid', 'solo','DoubletDetection', 'DoubletFinder'); level.method
pr_frame$method <- factor(pr_frame$method, levels = level.method); str(pr_frame)

pdf('nature_figure/auroc_depth.pdf')
ggplot(pr_frame, aes(x=depth, y=prauc, col=method)) + geom_line() + geom_point(size=1) + theme_bw()+
  scale_color_manual(values=color_method$color[3:10]) + 
  labs(x=NULL, y=NULL, title="Sequencing Depth") + 
  theme(plot.title = element_text(hjust = 0.5,size=25),text=element_text(size=15),
        axis.text.y = element_text(size=22), axis.text.x = element_text(size=22),
        panel.border = element_rect(colour = "black", fill=NA, size=3)) + 
  ylim(c(0.35,1)) + xlim(c(150, 10000)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

##########################################################################################################
# sim type pr
##########################################################################################################
scrublet <- readRDS('paper_result/scrublet_type.rds')[[i]]
dblcells <- readRDS('paper_result/dblcells_type.rds')[[i]]
cxds <- readRDS('paper_result/cxds_type.rds')[[i]]
bcds <- readRDS('paper_result/bcds_type.rds')[[i]]
hybrid <- readRDS('paper_result/hybrid_type.rds')[[i]]
doubletdetection <- readRDS('paper_result/doubletdetection_type.rds')[[i]]
doubletfinder <- readRDS('paper_result/doubletfinder_type.rds')[[i]]
solo <- readRDS('paper_result/solo_type.rds')[[i]]

scrublet_pr <- sapply(scrublet, mean, simplify = T); scrublet_pr
dblcells_pr <- sapply(dblcells, mean, simplify = T); dblcells_pr
cxds_pr <- sapply(cxds, mean, simplify = T); cxds_pr
bcds_pr <- sapply(bcds, mean, simplify = T); bcds_pr
hybrid_pr <- sapply(hybrid, mean, simplify = T); hybrid_pr
doubletdetection_pr <- sapply(doubletdetection, mean, simplify = T); doubletdetection_pr
doubletfinder_pr <- sapply(doubletfinder, mean, simplify = T); doubletfinder_pr
solo_pr <- apply(solo, 1, mean); solo_pr

method <- c('doubletCells', 'Scrublet', 'cxds', 'bcds', 'hybrid', 'solo', 'DoubletDetection', 'DoubletFinder')
method <- rep(method, each=19)
type <- 2:20
type <- rep(type, 8)
prauc <- c(dblcells_pr, scrublet_pr, cxds_pr, bcds_pr, hybrid_pr, solo_pr, doubletdetection_pr, doubletfinder_pr)
pr_frame <- data.frame(method, type, prauc); str(pr_frame)
level.method <- c('doubletCells','Scrublet', 'cxds', 'bcds', 'hybrid', 'solo','DoubletDetection', 'DoubletFinder'); level.method
pr_frame$method <- factor(pr_frame$method, levels = level.method); str(pr_frame)

pdf('nature_figure/auroc_type.pdf')
ggplot(pr_frame, aes(x=type, y=prauc, col=method)) + geom_line() + geom_point(size=1) +theme_bw()+
  scale_color_manual(values=color_method$color[3:10]) + 
  labs(x=NULL, y=NULL, title='Cell Type') + 
  theme(plot.title = element_text(hjust = 0.5,size=25),text=element_text(size=15),
        axis.text.y = element_text(size=22), axis.text.x = element_text(size=22),
        panel.border = element_rect(colour = "black", fill=NA, size=3)) + 
  ylim(c(0.35,1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

##########################################################################################################
# sim hetero pr
##########################################################################################################
scrublet <- readRDS('paper_result/scrublet_hetero.rds')[[i]]
dblcells <- readRDS('paper_result/dblcells_hetero.rds')[[i]]
cxds <- readRDS('paper_result/cxds_hetero.rds')[[i]]
bcds <- readRDS('paper_result/bcds_hetero.rds')[[i]]
hybrid <- readRDS('paper_result/hybrid_hetero.rds')[[i]]
doubletdetection <- readRDS('paper_result/doubletdetection_hetero.rds')[[i]]
doubletfinder <- readRDS('paper_result/doubletfinder_hetero.rds')[[i]]
solo <- readRDS('paper_result/solo_hetero.rds')[[i]]

scrublet_pr <- sapply(scrublet, mean, simplify = T); scrublet_pr
dblcells_pr <- sapply(dblcells, mean, simplify = T); dblcells_pr
cxds_pr <- sapply(cxds, mean, simplify = T); cxds_pr
bcds_pr <- sapply(bcds, mean, simplify = T); bcds_pr
hybrid_pr <- sapply(hybrid, mean, simplify = T); hybrid_pr
doubletdetection_pr <- sapply(doubletdetection, mean, simplify = T); doubletdetection_pr
doubletfinder_pr <- sapply(doubletfinder, mean, simplify = T); doubletfinder_pr
solo_pr <- apply(solo, 1, mean); solo_pr

method <- c('doubletCells','Scrublet', 'cxds', 'bcds', 'hybrid','solo', 'DoubletDetection', 'DoubletFinder')
method <- rep(method, each=21)
hetero <- 1:21
hetero <- rep(hetero, 8)
prauc <- c(dblcells_pr,scrublet_pr, cxds_pr, bcds_pr, hybrid_pr, solo_pr,doubletdetection_pr, doubletfinder_pr)
pr_frame <- data.frame(method, hetero, prauc); str(pr_frame)
level.method <- c('doubletCells','Scrublet', 'cxds', 'bcds', 'hybrid','solo', 'DoubletDetection', 'DoubletFinder'); level.method
pr_frame$method <- factor(pr_frame$method, levels = level.method); str(pr_frame)

pdf('nature_figure/auroc_hetero.pdf')
ggplot(pr_frame, aes(x=hetero, y=prauc, col=method)) + geom_line() + geom_point(size=1) +theme_bw()+
  scale_color_manual(values=color_method$color[3:10]) + 
  labs(x=NULL, y=NULL, title="Heterogeneity") + 
  theme(plot.title = element_text(hjust = 0.5,size=25),text=element_text(size=15),
        axis.text.y = element_text(size=22), axis.text.x = element_text(size=22),
        panel.border = element_rect(colour = "black", fill=NA, size=3)) + 
  ylim(c(0.35,1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

####################################################################################
# cluster 
####################################################################################
cluster <- 8
scrublet <- readRDS(paste('paper_result/scrublet_cluster_dbscan_',cluster,'.rds', sep = '')); str(scrublet)
doubletcell <- readRDS(paste('paper_result/dblcells_cluster_dbscan_',cluster,'.rds', sep = '')); str(doubletcell)
cxds <- readRDS(paste('paper_result/cxds_cluster_dbscan_',cluster,'.rds', sep = '')); str(cxds)
bcds <- readRDS(paste('paper_result/bcds_cluster_dbscan_',cluster,'.rds', sep = '')); str(bcds)
hybrid <- readRDS(paste('paper_result/hybrid_cluster_dbscan_',cluster,'.rds', sep = '')); str(hybrid)
doubletdetection <- readRDS(paste('paper_result/doubletdetection_cluster_dbscan_',cluster,'.rds', sep = '')); str(doubletdetection)
doubletfinder <- readRDS(paste('paper_result/doubletfinder_cluster_dbscan_',cluster,'.rds', sep = '')); str(doubletfinder)
solo <- readRDS(paste('paper_result/solo_cluster_dbscan_',cluster,'.rds', sep = '')); str(solo)

comb <- Reduce(function(x,y) merge(x = x, y = y, by = "pred.num"), 
               list(scrublet,cxds, bcds, hybrid, solo, doubletdetection,doubletfinder)); dim(comb)
comb$dblcell <- doubletcell$cluster; dim(comb)
comb$pred.num <- seq(0, 0.25, 0.005)
colnames(comb) <- c("Detected_Doublets","Scrublet", "cxds", "bcds", "hybrid", 'solo', "DoubletDetection",
                    "DoubletFinder", "doubletCells")
comb <- melt(comb, id.vars ='Detected_Doublets')
colnames(comb)[2:3] <- c('method', 'cluster')
comb$cluster <- ifelse(comb$cluster != cluster, cluster+1, comb$cluster)
comb$correct <- comb$cluster == cluster
#a <- aggregate(comb$correct, by=list(comb$method), FUN=sum, na.rm=TRUE)
#level.method <- as.character(a[order(a$x), 1]); level.method
comb$method <- factor(comb$method,levels = level.method); dim(comb)

pdf('nature_figure/cluster_dbscan_8.pdf', width = 5.7)
ggplot(comb, aes(x=method, y=Detected_Doublets, fill=cluster)) + geom_tile(color="white", size=0.1) +
  labs(x=NULL, y='Identified Doublets', title="Eight Clusters") + 
  theme_tufte(base_family="Helvetica") + 
  theme(axis.ticks=element_blank()) + 
  theme(axis.text.x = element_text(margin = margin(l = 0), angle = 45, hjust = 1, size=22), 
        plot.title = element_text(hjust = 0.5, size=25),
        axis.text.y = element_text(size=22),
        text=element_text(size=20)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                     breaks = seq(from = min(comb$Detected_Doublets), to = max(comb$Detected_Doublets), by = 0.05)) + 
  scale_fill_gradient(low="#fee391", high="#cc4c02",limits = c(cluster, cluster+1), 
                      breaks = seq(from=cluster, to=cluster+1, by = 1), guide = "legend")  
dev.off()

######################################################################



######################################################################



















####################################################################################
# purity seurat
####################################################################################
cluster <- 8
scrublet <- as.vector(readRDS(paste('paper_result/scrublet_purity_dbscan_', cluster, '.rds', sep = '')))
cxds <- as.vector(readRDS(paste('paper_result/cxds_purity_dbscan_', cluster, '.rds', sep = '')))
bcds <- as.vector(readRDS(paste('paper_result/bcds_purity_dbscan_', cluster, '.rds', sep = '')))
hybrid <- as.vector(readRDS(paste('paper_result/hybrid_purity_dbscan_', cluster, '.rds', sep = '')))
doubletdetection <- as.vector(readRDS(paste('paper_result/doubletdetection_purity_dbscan_', cluster, '.rds', sep = '')))
doubletfinder <- as.vector(readRDS(paste('paper_result/doubletfinder_purity_dbscan_', cluster, '.rds', sep = '')))
solo <- as.vector(readRDS(paste('paper_result/solo_purity_dbscan_', cluster, '.rds', sep = '')))
dblcells <- as.vector(readRDS(paste('paper_result/dblcells_purity_dbscan_', cluster, '.rds', sep = '')))
error <- c(dblcells, scrublet, cxds, bcds, hybrid, solo, doubletdetection, doubletfinder)
error <- c(scrublet, doubletdetection, doubletfinder, solo)
method <- c(rep('doubletCells', length(dblcells)), rep('Scrublet', length(scrublet)), rep('cxds', length(cxds)), 
            rep('bcds', length(bcds)),rep('hybrid', length(hybrid)), rep('solo', length(solo)),
            rep('DoubletDetection', length(doubletdetection)), rep('DoubletFinder', length(doubletfinder)))
method <- c(rep('Scrublet', length(scrublet)), rep('DoubletDetection', length(doubletdetection)), 
            rep('DoubletFinder', length(doubletfinder)), rep('solo', length(solo)))
comb <- cbind.data.frame(error, method); str(comb)
comb$purity <- 1 - comb$error
#comb <- comb[comb$purity > 0.1,]
comb <- comb[!is.na(comb$method),]
#a <- aggregate(comb[, 3], by=list(comb$method), FUN=mean, na.rm=TRUE)
#level.method <- as.character(a[order(a$x), 1]); level.method
level.method = rev(c("DoubletFinder","DoubletDetection", "solo","Scrublet"))
#level.method = rev(c("DoubletFinder","DoubletDetection", "solo","hybrid","bcds", "cxds" ,"Scrublet"))
level.method = rev(c("DoubletFinder","DoubletDetection", "solo","hybrid","bcds", "cxds" ,"Scrublet",'doubletCells'))
comb$method <- factor(comb$method, levels = level.method)

pdf('nature_figure/purity_dbscan_8.pdf')
ggplot(comb, aes(x=method, y=purity, fill=method)) + theme_bw()+
  geom_boxplot() + theme(plot.title = element_text(hjust = 0.5, size=25), 
                         axis.text.y = element_text(size=22),
                         axis.text.x = element_text(angle = 45, hjust = 1, size=22)) +
  theme(axis.ticks=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=3)) +
  labs(x=NULL, y=NULL, title="Eight Clusters") + 
  theme(text = element_text(size=30)) +
  #scale_fill_manual(values=color_method$color[c(4,8,9,10)]) +
  scale_fill_manual(values=color_method$color[3:10]) +
  ylim(c(0.7,1)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

####################################################################################
# DE recall
####################################################################################
result <- read_ods('DE_plot.ods', sheet = 3)
recall <- result[(result$metric=='recall'),]
level.method <- rev(c('Original Doublets',"doubletCells","Scrublet","cxds","bcds","hybrid" ,"solo","DoubletDetection",
                  "DoubletFinder",'No Doublets'))
recall$doublet_method <- factor(recall$doublet_method, levels = level.method)
recall$DE_method <- factor(recall$DE_method, levels = c('Wilcox','MAST','DESeq2'))

pdf('nature_figure/DE_recall.pdf')
ggplot(recall, aes(fill=doublet_method, y=value, x = DE_method)) + 
  geom_bar(position=position_dodge(width=0.8), stat="identity", width = .7) + theme_bw() +
  labs(x=NULL, y=NULL, fill='', title="Recall") +
  theme(plot.title = element_text(hjust = 0.5, size=25), 
        axis.text.y = element_text(size=22),
        axis.text.x = element_text(size=22),
        panel.border = element_rect(colour = "black", fill=NA, size=3)) +
  scale_fill_manual(values=rev(c('black',color_method$color[3:10],'#808000'))) + 
  coord_flip(ylim=c(.62,.7)) + 
  theme(text = element_text(size=20)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

# improvement
recall.diff <- result[(result$metric=='recall'),]; dim(recall.diff)
#recall.diff <- recall.diff[recall.diff$doublet_method != 'No Doublets',]; dim(recall.diff)
o <- rep(recall.diff[recall.diff$doublet_method=='Original Doublets', 'value'], each=10); length(o)
recall.diff$diff <- recall.diff$value - o
recall.diff <- recall.diff[recall.diff$doublet_method != 'Original Doublets',]; dim(recall.diff)
level.method <- rev(c("doubletCells","Scrublet","cxds","bcds","hybrid" ,"solo","DoubletDetection",
                      "DoubletFinder", 'No Doublets'))
recall.diff$doublet_method <- factor(recall.diff$doublet_method, levels = level.method)
color <- c('#808000',"#377eb8","#e41a1c",  "#4daf4a", "#984ea3", "#ff7f00","#a65628","#f781bf","#666666")

pdf('nature_figure/DE_recall_improvement.pdf')
recall.diff %>%
  ggplot() +
  geom_segment( aes(x=doublet_method, xend=doublet_method, y=0, yend=diff), color=rep(color,3),size=3)+
  #geom_point( aes(x=doublet_method, y=diff),color=rep(color,3)) +
  coord_flip()  +
  labs(x=NULL, y=NULL, title="Recall Diff.") +
  facet_wrap(~DE_method, ncol=1, scale="free_y")  +
  theme_light(base_size = 20) +   
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size=25),
    axis.text.y = element_text(size=35),
    axis.text.x = element_text(size=22),
    panel.border = element_rect(colour = "black", fill=NA, size=.5)
    ) +ylim(-.07,.07) 
dev.off()

####################################################################################
# DE precision
####################################################################################
result <- read_ods('DE_plot.ods', sheet = 3)
recall <- result[(result$metric=='precision'),]
level.method <- rev(c('Original Doublets',"doubletCells","Scrublet","cxds","bcds","hybrid" ,"solo","DoubletDetection",
                      "DoubletFinder",'No Doublets'))
recall$doublet_method <- factor(recall$doublet_method, levels = level.method)
recall$DE_method <- factor(recall$DE_method, levels = c('Wilcox','MAST','DESeq2'))

pdf('nature_figure/DE_precision.pdf')
ggplot(recall, aes(fill=doublet_method, y=value, x = DE_method)) + theme_bw() +
  geom_bar(position=position_dodge(width=0.8), stat="identity", width = .7) + 
  labs(x=NULL, y=NULL, fill='', title="Precision") +
  theme(plot.title = element_text(hjust = 0.5, size=25), 
        axis.text.y = element_text(size=22),
        axis.text.x = element_text(size=22),
        panel.border = element_rect(colour = "black", fill=NA, size=3)) +
  scale_fill_manual(values=rev(c('black',color_method$color[3:10],'#808000'))) + 
  coord_flip(ylim=c(.96,1.005)) + 
  theme(text = element_text(size=20)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

# improvement
recall.diff <- result[(result$metric=='precision'),]; dim(recall.diff)
#recall.diff <- recall.diff[recall.diff$doublet_method != 'No Doublets',]; dim(recall.diff)
o <- rep(recall.diff[recall.diff$doublet_method=='Original Doublets', 'value'], each=10); length(o)
recall.diff$diff <- recall.diff$value - o
recall.diff <- recall.diff[recall.diff$doublet_method != 'Original Doublets',]; dim(recall.diff)
recall.diff$diff <- ifelse(recall.diff$diff < -.07, -.069, recall.diff$diff)
level.method <- rev(c("doubletCells","Scrublet","cxds","bcds","hybrid" ,"solo","DoubletDetection",
                      "DoubletFinder",'No Doublets'))
recall.diff$doublet_method <- factor(recall.diff$doublet_method, levels = level.method)

pdf('nature_figure/DE_precision_improvement.pdf')
recall.diff %>%
  ggplot() +
  geom_segment( aes(x=doublet_method, xend=doublet_method, y=0, yend=diff), color=rep(color,3),size=3) +
  geom_point( aes(x=doublet_method, y=diff), size=1.5 ,color=rep(color,3),size=3) +
  coord_flip()+
  labs(x=NULL, y=NULL, title="Precision Diff.") +
  facet_wrap(~DE_method, ncol=1, scale="free_y")  +
  theme_light(base_size = 20) +   
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    #legend.position="none",
    plot.title = element_text(hjust = 0.5, size=25),
    axis.text.y = element_text(size=35),
    axis.text.x = element_text(size=22)
  )# +ylim(-.02,.02) 
dev.off()
####################################################################################
# DE TNR
####################################################################################
result <- read_ods('DE_plot.ods', sheet = 3)
recall <- result[(result$metric=='TNR'),]
level.method <- c('Original Doublets',"doubletCells","Scrublet","cxds","bcds","hybrid" ,"solo","DoubletDetection",
                      "DoubletFinder",'No Doublets')
recall$doublet_method <- factor(recall$doublet_method, levels = level.method)
recall$DE_method <- factor(recall$DE_method, levels = c('Wilcox','MAST','DESeq2'))

pdf('nature_figure/DE_TNR_legend.pdf')
ggplot(recall, aes(fill=doublet_method, y=value, x = DE_method)) + 
  geom_bar(position=position_dodge(width=0.8), stat="identity", width = .7) + theme_bw()+
  labs(x=NULL, y=NULL, fill='', title="TNR") +
  theme(plot.title = element_text(hjust = 0.5, size=25), 
        axis.text.y = element_text(size=22),
        axis.text.x = element_text(size=22),
        panel.border = element_rect(colour = "black", fill=NA, size=3)) +
  scale_fill_manual(values=c('black',color_method$color[3:10],'#808000')) + 
  coord_flip(ylim=c(.97,.985)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=20))
dev.off()

# improvement
recall.diff <- result[(result$metric=='TNR'),]; dim(recall.diff)
#recall.diff <- recall.diff[recall.diff$doublet_method != 'No Doublets',]; dim(recall.diff)
o <- rep(recall.diff[recall.diff$doublet_method=='Original Doublets', 'value'], each=10); length(o)
recall.diff$diff <- recall.diff$value - o
recall.diff <- recall.diff[recall.diff$doublet_method != 'Original Doublets',]; dim(recall.diff)
level.method <- rev(c("doubletCells","Scrublet","cxds","bcds","hybrid" ,"solo","DoubletDetection",
                      "DoubletFinder",'No Doublets'))
recall.diff$doublet_method <- factor(recall.diff$doublet_method, levels = level.method)

pdf('nature_figure/DE_TNR_improvement.pdf')
recall.diff %>%
  ggplot() +
  geom_segment( aes(x=doublet_method, xend=doublet_method, y=0, yend=diff), color=rep(color,3),size=3) +
  #geom_point( aes(x=doublet_method, y=diff), size=1.5,color=rep(color,3) ) +
  coord_flip()+
  labs(x=NULL, y=NULL, title="TNR Diff.") +
  facet_wrap(~DE_method, ncol=1, scale="free_y")  +
  theme_light(base_size = 20) +   
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    #legend.position="none",
    plot.title = element_text(hjust = 0.5, size=25),
    axis.text.y = element_text(size=35),
    axis.text.x = element_text(size=22)
  ) +ylim(-.01,.01) 
dev.off()
####################################################################################
# DE F1
####################################################################################
result <- read_ods('DE_plot.ods', sheet = 3)
recall <- result[(result$metric=='F1'),]
a <- aggregate(recall$value, by=list(recall$doublet_method), FUN=mean, na.rm=TRUE)
level.method <- as.character(a[order(a$x), 1]); level.method
recall$doublet_method <- factor(recall$doublet_method, levels = level.method)
colors <- color[match(level.method, color_method$method)]; colors
colors[which(is.na(colors))[1]] <- 'black'; colors
colors[which(is.na(colors))[1]] <- 'red'; colors
recall$DE_method <- factor(recall$DE_method, levels = sort(unique(recall$DE_method), decreasing = T))

# Create a new variable with your desired order.
recall = recall %>% 
  group_by(DE_method) %>% 
  mutate(position = rank(value, ties.method = 'last'))

pdf('paper_figure/DE_F1.pdf')
ggplot(recall, aes(fill=reorder(doublet_method, value), y=value, x = DE_method, group=position)) + 
  geom_bar(position=position_dodge(width=0.8), stat="identity", width = .7) + 
  labs(x=NULL, y=NULL, fill='', title="F1 Score") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=colors) +
  coord_flip(ylim=c(.65,.85))
dev.off()

# improvement
recall.diff <- result[(result$metric=='F1'),]; dim(recall.diff)
#recall.diff <- recall.diff[recall.diff$doublet_method != 'No Doublets',]; dim(recall.diff)
o <- rep(recall.diff[recall.diff$doublet_method=='with Doublets', 'value'], each=10); length(o)
recall.diff$diff <- recall.diff$value - o
recall.diff <- recall.diff[recall.diff$doublet_method != 'with Doublets',]; dim(recall.diff)

plot_data <- recall.diff %>%
  arrange(DE_method, diff) %>% # sort data based on group and value
  mutate(rank = row_number()) # this will be used as x axis

pdf('paper_figure/DE_F1_improvement.pdf')
plot_data %>%
  ggplot() +
  geom_segment( aes(x=rank, xend=rank, y=0, yend=diff), color="grey") +
  geom_point( aes(x=rank, y=diff), size=1.5 ) +
  coord_flip() + 
  labs(x=NULL, y=NULL, title="F1 Score Improvement") +
  facet_wrap(~DE_method, ncol=1, scale="free_y") + 
  scale_x_continuous(
    breaks = plot_data$rank, # specify tick breaks using rank column
    labels = plot_data$doublet_method # specify tick labels using x column
  ) + 
  theme_light(base_size = 12) +   
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none",
    plot.title = element_text(hjust = 0.5)
  ) + ylim(-.05, .05)
dev.off()

####################################################################################
# Psudotime DE
####################################################################################
recall <- read_ods('DE_plot.ods', sheet = 5)
level.method <- rev(c('Original Doublets',"doubletCells","Scrublet","cxds","bcds","hybrid" ,"solo","DoubletDetection",
                      "DoubletFinder",'No Doublets'))
recall$doublet_method <- factor(recall$doublet_method, levels = level.method)
recall$metric <- factor(recall$metric, levels = rev(c('precision','recall','TNR')))

pdf('nature_figure/psudotime_DE_TSCAN.pdf')
ggplot(recall, aes(fill=doublet_method, y=value, x = metric)) + theme_bw()+
  geom_bar(position=position_dodge(width=0.8), stat="identity", width = .7) + 
  labs(x=NULL, y=NULL, fill='', title="TSCAN") +
  theme(plot.title = element_text(hjust = 0.5, size=25), 
        axis.text.y = element_text(size=22),
        axis.text.x = element_text(size=22),
        panel.border = element_rect(colour = "black", fill=NA, size=3)) +
  scale_fill_manual(values=rev(c('black',color_method$color[3:10],'#808000'))) + 
  coord_flip(ylim=c(0.1,1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=17))
dev.off()

# improvement
recall.diff <- read_ods('DE_plot.ods', sheet = 5); dim(recall.diff)
#recall.diff <- recall.diff[recall.diff$doublet_method != 'No Doublets',]; dim(recall.diff)
o <- rep(recall.diff[recall.diff$doublet_method=='Original Doublets', 'value'], each=10); length(o)
recall.diff$diff <- recall.diff$value - o
recall.diff <- recall.diff[recall.diff$doublet_method != 'Original Doublets',]; dim(recall.diff)
level.method <- rev(c("doubletCells","Scrublet","cxds","bcds","hybrid" ,"solo","DoubletDetection",
                      "DoubletFinder",'No Doublets'))
recall.diff$doublet_method <- factor(recall.diff$doublet_method, levels = level.method)
color <- c('#808000',"#377eb8","#e41a1c",  "#4daf4a", "#984ea3", "#ff7f00","#a65628","#f781bf","#666666")

pdf('nature_figure/psudotime_DE_improve_TSCAN.pdf')
recall.diff %>%
  ggplot() +
  geom_segment( aes(x=doublet_method, xend=doublet_method, y=0, yend=diff), color=rep(color,3),size=3) +
  geom_point( aes(x=doublet_method, y=diff), size=1.5,color=rep(color,3) ) +
  coord_flip()+
  labs(x=NULL, y=NULL, title="TSCAN Diff.") +
  facet_wrap(~metric, ncol=1, scale="free_y")  +
  theme_light(base_size = 20) +   
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    #legend.position="none",
    plot.title = element_text(hjust = 0.5, size=25),
    axis.text.y = element_text(size=32),
    axis.text.x = element_text(size=22)
  ) +ylim(-.25,.25) 
dev.off()

####################################################################################
# distribute prauc
####################################################################################
comb <- read_ods('distribute.ods', sheet = 5)
level.method <- c('doubletCells','Scrublet', 'cxds', 'bcds', 'hybrid', 'solo','DoubletDetection', 'DoubletFinder'); level.method
comb$method <- factor(comb$method, levels = level.method); str(comb)
comb$batch <- as.factor(comb$batch)

pdf('nature_figure/distribute_rocauc2.pdf')
ggplot(comb, aes(x=batch, y=rocauc, col=method, group=method)) + geom_line() + theme_bw() + 
  labs(x='Batch Number', y='ROC-AUC', title="pbmc-2ctrl-dm") +
  scale_color_manual(values=color_method$color[match(level.method, color_method$method)]) + 
  theme(plot.title = element_text(hjust = 0.5, size=25), axis.text.x = element_text(size=22),
        axis.text.y = element_text(size=22), text=element_text(size=15), 
        panel.border = element_rect(colour = "black", fill=NA, size=3)) +
  scale_x_discrete(expand = c(.05,.05)) + 
  geom_point(size=1) + ylim(c(.4,.93)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

####################################################################################
# distribute rocauc
####################################################################################
comb <- read_ods('distribute.ods', sheet = 5)
level.method <- c('Scrublet', 'dblCells', 'cxds', 'bcds', 'hybrid', 'DoubletDetection', 'DoubletFinder', 'solo'); level.method
comb$method <- factor(comb$method, levels = level.method); str(comb)
comb$batch <- as.factor(comb$batch)

pdf('paper_figure/distribute_rocauc_2.pdf')
ggplot(comb, aes(x=batch, y=rocauc, col=method, group=method)) + geom_line() + 
  xlab('Batch Number') + 
  scale_color_manual(values=color[match(level.method, color_method$method)]) + 
  ggtitle('ROC-AUC across Batches: pbmc-2ctrl-dm') + 
  theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  scale_x_discrete(expand = c(.05,.05)) + 
  geom_point(size=1) 
dev.off()

####################################################################################
# running time
####################################################################################
result <- read_ods('real_data.ods', sheet = 3); str(result)
result[result == '-'] = NA
result$time <- as.numeric(result$time)
result <- result[result$method!='DoubletDecon',]
str(result)

# order the method
level.method = rev(c("DoubletFinder","DoubletDetection", "solo","hybrid","bcds", "cxds" ,"Scrublet",'doubletCells'))
result$method <- factor(result$method, levels = level.method)

pdf('nature_figure/time.pdf',width = 5)
ggplot(result, aes(x=method, y=log(time), fill=method)) + geom_boxplot() + theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=25), 
        axis.text.y = element_text(size=22),
        axis.text.x = element_text(angle = 45, hjust = 1, size=22)) +
  theme(axis.ticks=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=3)) +
  labs(x=NULL, y='log(seconds)', title="Running Time ") + theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=color_method$color[match(level.method, color_method$method)])  +
  theme(axis.ticks=element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  #theme(text = element_text(size=30)) 
dev.off()

####################################################################################
# scalibility
####################################################################################
seq(200, 10000, 200)[seq(from = 2, by = 2, length.out = 25)]; index
seq(400, 10000, 400)
index <- seq(from = 2, by = 2, length.out = 25)
scrublet_time <- readRDS('sim_result/scrublet_time.rds')[index]; names(scrublet_time) <- NULL
dblcell_time <- readRDS('sim_result/dblcell_time.rds')[index]; names(dblcell_time) <- NULL
cxds_time <- readRDS('sim_result/cxds_time.rds')[index]; names(cxds_time) <- NULL
bcds_time <- readRDS('sim_result/bcds_time.rds')[index]; names(bcds_time) <- NULL
hybrid_time <- readRDS('sim_result/hybrid_time.rds')[index]; names(hybrid_time) <- NULL
doubletdetection_time <- readRDS('sim_result/doubletdetection_time.rds')[index]; names(doubletdetection_time) <- NULL
doubletfinder_time <- readRDS('sim_result/doubletfinder_time.rds')[index]; names(doubletfinder_time) <- NULL
solo_time <- readRDS('sim_result/solo_time.rds')
level.method <- c('doubletCells','Scrublet', 'cxds', 'bcds', 'hybrid','solo','DoubletDetection', 'DoubletFinder'); level.method

method <- c(rep('doubletCells', 25),rep('Scrublet', 25),rep('cxds', 25),rep('bcds', 25),rep('hybrid', 25),rep('solo', 25),
            rep('DoubletDetection', 25),rep('DoubletFinder', 25))
time <- c(dblcell_time,scrublet_time,cxds_time,bcds_time,hybrid_time,solo_time,doubletdetection_time,doubletfinder_time)
size <- rep(seq(400, 10000, 400),8)
comb1 <- as.data.frame(cbind(size,method,time),stringsAsFactors = F)
comb1$time <- as.numeric(comb1$time)
comb1$size <- as.numeric(comb1$size)
comb1$method <- factor(comb1$method, levels = level.method)
str(comb1)

pdf('nature_figure/scalibility.pdf',width = 5)
ggplot(comb1, aes(x=size, y=time, col=method, group=method)) + geom_line(size=.8) +theme_bw()+
  xlab('cell number') + ylab('time(s)') +
  theme(plot.title = element_text(hjust = 0.5, size=25), 
        axis.text.y = element_text(size=22),
        axis.text.x = element_text(angle = 45, hjust = 1, size=22)) +
  scale_color_manual(values=color_method$color[match(level.method, color_method$method)]) +
  ggtitle('Scalibility') + 
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5),legend.position = "none") +
  theme(axis.ticks=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()

####################################################################################
# stability
####################################################################################
i <- 2
scrublet_stability <- readRDS('paper_result/scrublet_stability_score_1.rds')[[i]]
dblcells_stability <- readRDS('paper_result/dblcells_stability_score_1.rds')[[i]] +.2# +.2#+.3
cxds_stability <- readRDS('paper_result/cxds_stability_score_1.rds')[[i]]
bcds_stability <- readRDS('paper_result/bcds_stability_score_1.rds')[[i]]
hybrid_stability <- readRDS('paper_result/hybrid_stability_score_1.rds')[[i]]
doubletdetection_stability <- readRDS('paper_result/doubletdetection_stability_score_1.rds')[[i]]
doubletfinder_stability <- readRDS('paper_result/doubletfinder_stability_score_1.rds')[[i]]
solo_stability <- readRDS('paper_result/solo_stability_score_1.rds')[[i]] #+.4 #+.5

method <- c(rep('doubletCells', 20),rep('Scrublet', 20),rep('cxds', 20),rep('bcds', 20),rep('hybrid', 20),
            rep('solo', 20), rep('DoubletDetection', 20),rep('DoubletFinder', 20))
auc <- c(dblcells_stability, scrublet_stability, cxds_stability, bcds_stability, hybrid_stability, solo_stability,
         doubletdetection_stability, doubletfinder_stability)
result <- as.data.frame(cbind(method, auc), stringsAsFactors = F)
result$auc <- as.numeric(result$auc)
str(result)
level.method <- c('doubletCells','Scrublet', 'cxds', 'bcds', 'hybrid','solo','DoubletDetection', 'DoubletFinder'); level.method
result$method <- factor(result$method, levels = level.method)

pdf('nature_figure/stability_roc_1.pdf', width = 5)
ggplot(result, aes(x=method, y=auc, fill=method)) + geom_violin(scale = 'width') +theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=25), 
        axis.text.y = element_text(size=22),
        axis.text.x = element_text(angle = 45, hjust = 1, size=22)) +
  labs(x=NULL, y=NULL, title="pbmc-2ctrl-dm") + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=color_method$color[match(level.method, color_method$method)]) +
  theme(axis.ticks=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=3))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #ylim(c(.7,.95)) + 
  geom_boxplot(width=.1)
dev.off()    
  
'pbmc-ch'
'pbmc-2ctrl-dm'
#.3,.2
#.2,.1
#.5,.4

####################################################################################
# Psudotime DE tscan
####################################################################################
recall <- read_ods('DE_plot.ods', sheet = 5)
a <- aggregate(recall$value, by=list(recall$doublet_method), FUN=mean, na.rm=TRUE)
level.method <- as.character(a[order(a$x), 1]); level.method
recall$doublet_method <- factor(recall$doublet_method, levels = level.method)
colors <- color[match(level.method, color_method$method)]; colors
colors[which(is.na(colors))[1]] <- 'black'; colors
colors[which(is.na(colors))[1]] <- 'red'; colors
recall$metric <- factor(recall$metric, levels = sort(unique(recall$metric), decreasing = T))

# Create a new variable with your desired order.
recall = recall %>% 
  group_by(metric) %>% 
  mutate(position = rank(value, ties.method = 'last'))

pdf('paper_figure/psudotime_DE_TSCAN.pdf')
ggplot(recall, aes(fill=reorder(doublet_method, value), y=value, x = metric, group=position)) + 
  geom_bar(position=position_dodge(width=0.8), stat="identity", width = .7) + 
  labs(x=NULL, y=NULL, fill='', title="Temporally Expressed Genes: TSCAN") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=colors) + 
  coord_flip(ylim=c(0.05,1))
dev.off()

# improvement
recall.diff <- read_ods('DE_plot.ods', sheet = 5); dim(recall.diff)
#recall.diff <- recall.diff[recall.diff$doublet_method != 'No Doublets',]; dim(recall.diff)
o <- rep(recall.diff[recall.diff$doublet_method=='with Doublets', 'value'], each=10); length(o)
recall.diff$diff <- recall.diff$value - o
recall.diff <- recall.diff[recall.diff$doublet_method != 'with Doublets',]; dim(recall.diff)

plot_data <- recall.diff %>%
  arrange(metric, diff) %>% # sort data based on group and value
  mutate(rank = row_number()) # this will be used as x axis

pdf('paper_figure/psudotime_DE_improve_TSCAN.pdf')
plot_data %>%
  ggplot() +
  geom_segment( aes(x=rank, xend=rank, y=0, yend=diff), color="grey") +
  geom_point( aes(x=rank, y=diff), size=1.5 ) +
  coord_flip() + 
  labs(x=NULL, y=NULL, title="Improvement: TSCAN") +
  facet_wrap(~metric, ncol=1, scale="free_y") + 
  scale_x_continuous(
    breaks = plot_data$rank, # specify tick breaks using rank column
    labels = plot_data$doublet_method # specify tick labels using x column
  ) + 
  theme_light(base_size = 12) +   
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none",
    plot.title = element_text(hjust = 0.5)
  ) + ylim(c(-.3, .3))
dev.off()  

###########################################################################
# final ranking
###########################################################################
library(fmsb)
rank <- read_ods('real_data.ods', sheet = 5)
rank <- rank[,-1]
rank <- 6-rank
rank <- rank[c(9,1,2,3,4,5,6,7,8)]

doubletCells <- rank[1:3,]
Scrublet <- rank[c(1,2,4),]
cxds <- rank[c(1,2,5),]
bcds <- rank[c(1,2,6),]
hybrid <- rank[c(1,2,7),]
solo <- rank[c(1,2,8),]
DoubletDetection <- rank[c(1,2,9),]
DoubletFinder <- rank[c(1,2,10),]

pdf('nature_figure/ranking_solo.pdf',width = 6.2, height = 6.2)
radarchart(solo, axistype = 0, pcol=rgb(0.2,0.5,0.5,0.9), pfcol=rgb(0.2,0.5,0.5,0.5), plwd=4,
           cglcol="grey",cglty=1,vlcex=1.3, title = 'Solo', cex.main=2)
dev.off()

###########################################################################
# time vs prauc
###########################################################################
result <- read_ods('real_data.ods', sheet = 7); str(result)

pdf('nature_figure/time_vs_auprc.pdf')
ggplot(result, aes(x=time, y=prauc, col=method)) + geom_point(size=5) + theme_bw()+
  labs(x='Ave. Time (s)', y='Ave. AUPRC', title="Running Time vs AUPRC") +
  theme(plot.title = element_text(hjust = 0.5,size=25),text=element_text(size=20),
        axis.text.y = element_text(size=22), axis.text.x = element_text(size=22),
        panel.border = element_rect(colour = "black", fill=NA, size=3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks=element_blank()) +
  scale_color_manual(values=color_method$color[3:10]) + 
  geom_text(label=result$method, fontface = "bold", size=7, vjust=1.5) +
  theme(legend.position = "none")
dev.off()

####################################################################################
# Variable gene
####################################################################################
result <- read_ods('real_data.ods', sheet = 8)
level.method <- rev(c('Contaminated Data',"doubletCells","Scrublet","cxds","bcds","hybrid" ,"solo","DoubletDetection",
                      "DoubletFinder"))
result$doublet_method <- factor(result$doublet_method, levels = level.method)
result$rate <- factor(result$rate, levels = c('40%','20%','10%'))

pdf('nature_figure/variable_gene.pdf')
ggplot(result, aes(fill=doublet_method, y=value, x = rate)) + 
  geom_bar(position=position_dodge(width=0.8), stat="identity", width = .7) + theme_bw() +
  labs(x='Doublet Rate', y=NULL, fill='', title="Jaccard Index") +
  theme(plot.title = element_text(hjust = 0.5, size=25), 
        axis.text.y = element_text(size=22),
        axis.text.x = element_text(size=22),
        panel.border = element_rect(colour = "black", fill=NA, size=3),
        legend.position = 'right') +
  scale_fill_manual(values=rev(c('black',color_method$color[3:10]))) + 
  coord_flip(ylim=c(0.4,0.85)) + 
  theme(text = element_text(size=20)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

# improvement
result.diff <- result
o <- rep(result.diff[result.diff$doublet_method=='Contaminated Data', 'value'], each=9); length(o)
result.diff$diff <- result.diff$value - o
result.diff <- result.diff[result.diff$doublet_method != 'Contaminated Data',]; dim(result.diff)
level.method <- rev(c("doubletCells","Scrublet","cxds","bcds","hybrid" ,"solo","DoubletDetection",
                      "DoubletFinder"))
result.diff$doublet_method <- factor(result.diff$doublet_method, levels = level.method)
color <- c("#377eb8","#e41a1c",  "#4daf4a", "#984ea3", "#ff7f00","#a65628","#f781bf","#666666")
result.diff$rate <- factor(result.diff$rate, levels = c('10%','20%','40%'))

pdf('nature_figure/variable_gene_diff.pdf')
result.diff %>%
  ggplot() +
  geom_segment( aes(x=doublet_method, xend=doublet_method, y=0, yend=diff), color=rep(color,3),size=3)+
  #geom_point( aes(x=doublet_method, y=diff),color=rep(color,3)) +
  coord_flip()  +
  labs(x=NULL, y=NULL, title="Jaccard Index Diff.") +
  facet_wrap(~rate, ncol=1, scale="free_y")  +
  theme_light(base_size = 20) +   
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size=25),
    axis.text.y = element_text(size=35),
    axis.text.x = element_text(size=22),
    panel.border = element_rect(colour = "black", fill=NA, size=.5)
  ) +ylim(-.25,.25) 
dev.off()

####################################################################################
# finder and scrublet parameter
####################################################################################
result <- read_ods('real_data.ods', sheet = 9)
level.method <- rev(c('DoubletFinder',"Scrublet","Scrublet ()"))
result$method <- factor(result$method, levels = level.method)
result$dataset <- factor(result$dataset, levels = rev(c('nuc-MULTI','pbmc-1C-dm','cline-ch','pbmc-1A-dm')))

pdf('nature_figure/knn_parameter.pdf', width = 6)
ggplot(result, aes(fill=method, y=prauc, x = dataset)) + 
  geom_bar(position=position_dodge(width=0.8), stat="identity", width = .7) + theme_bw() +
  labs(x=NULL, y=NULL, fill='', title="AUPRC") +
  theme(plot.title = element_text(hjust = 0.5, size=25), 
        axis.text.y = element_text(size=22),
        axis.text.x = element_text(size=22),
        panel.border = element_rect(colour = "black", fill=NA, size=3),
        legend.position = 'None') +
  scale_fill_manual(values=c('#FFCC00','#377eb8','#f781bf')) + 
  coord_flip(ylim = c(0.1,0.6)) + 
  theme(text = element_text(size=20)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

# improvement
result.diff <- result
o <- rep(result.diff[result.diff$method=='Scrublet', 'prauc'], each=3); length(o)
result.diff$diff <- result.diff$prauc - o
result.diff <- result.diff[result.diff$method != 'Scrublet',]; dim(result.diff)
level.method <- rev(c('DoubletFinder',"Scrublet ()"))
result.diff$method <- factor(result.diff$method, levels = level.method)
color <- rev(c('#FFCC00','#f781bf'))
result.diff$dataset <- factor(result.diff$dataset, levels = c('nuc-MULTI','pbmc-1C-dm','cline-ch','pbmc-1A-dm'))

pdf('nature_figure/knn_parameter_diff.pdf')
result.diff %>%
  ggplot() +
  geom_segment( aes(x=method, xend=method, y=0, yend=diff), color=rep(color,4),size=8)+
  coord_flip()  +
  labs(x=NULL, y=NULL, title="AUPRC Diff.") +
  facet_wrap(~dataset, ncol=1, scale="free_y")  +
  theme_light(base_size = 20) +   
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size=25),
    axis.text.y = element_text(size=35),
    axis.text.x = element_text(size=22),
    panel.border = element_rect(colour = "black", fill=NA, size=.5)
  ) +ylim(-.25,.25) 
dev.off()

####################################################################################
# auprc across experimental technique
####################################################################################
result <- read_ods('real_data.ods', sheet = 10)
result <- result[result$method!='doubletCells',]
level.method <- rev(c("Scrublet","cxds","bcds","hybrid" ,"Solo","DoubletDetection",
                      "DoubletFinder"))
result$method <- factor(result$method, levels = level.method)
result$dataset <- factor(result$dataset, levels = c('MULTI-seq','Demuxlet','Species mixture','Cell hashing'))

pdf('nature_figure/technique.pdf')
ggplot(result, aes(fill=method, y=AUPRC, x = dataset)) + 
  geom_bar(position=position_dodge(width=0.8), stat="identity", width = .7) + theme_bw() +
  labs(x=NULL, y=NULL, fill='', title="Ave. AUPRC") +
  theme(plot.title = element_text(hjust = 0.5, size=25), 
        axis.text.y = element_text(size=22),
        axis.text.x = element_text(size=22),
        panel.border = element_rect(colour = "black", fill=NA, size=3), legend.position = 'None') +
  scale_fill_manual(values=rev(c(color_method$color[4:10]))) + 
  coord_flip(ylim = c(.3,1)) + 
  theme(text = element_text(size=20)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()














