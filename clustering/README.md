# Effects of doublet detection on cell clustering
## The name of the file indicates the doublet detection methods and clustering methods

doubletcells_cluster_louvain.R/doubletcells_cluster_dbscan.R: doubletCells with louvain and dbscan clustering

doubletdetection_cluster_louvain.R/doubletdetection_cluster_dbscan.R: DoubletDetection with louvain and dbscan clustering

doubletfinder_cluster_louvain.R/doubletfinder_cluster_dbscan.R: DoubletFinder with louvain and dbscan clustering

scds_cluster_louvain.R/scds_cluster_dbscan.R: cxds, scds, and hybrid with louvain and dbscan clustering

scrublet_cluster_louvain.R/scrublet_cluster_dbscan.R: Scrublet with louvain and dbscan clustering

solo_cluster_louvain.R/solo_cluster_dbscan.R: Solo with louvain and dbscan clustering

## The functionality of each file

Each file includes the following code blocks that implement different functionality. They are indicated by detailed comments within each file.

1. Clusters identified based on different identification rates. The identification rates increased from 0 to 0.25, increased by 0.01. The true cluster number is 4, 6, or 8.

2. The proportion of singlets in correctly identified clusters from step one (cluster quality or purity)



