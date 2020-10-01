# Performance of doublet-detection methods under distributed computing
This folder contains the benchmark of eight doublet-detection methods in distributed computing. Two real datasets were randomly split into 2, 4, 6, 8, or 10 batches. Each method was executed on each batch, and the corresponding batch scores were combined together. The AUPRC and AUROC calculated on merged scores were compared with the original metrics. 

## The name of the file indicates the doublet detection methods

1. doubletcells_distribute.R: doubletCells

2. doubletdetection_distribute.R: DoubletDetection

3. doubletfinder_distribute.R: DoubletFinder

4. scds_distribute.R: cxds, bcds, and hybrid

5. scrublet_distribute.R: Scrublet

6. solo_distribute.R: Solo

## The functionality of each file

1. Each file contains the code block to randomly split the dataset into 2, 4, 6, 8, and 10 batches. 

2. Each doublet-detection method was applied on each post-doublet-detection batch to obtain doublet scores. 

3. The AUPRC and AUROC were calculated on the doublet scores merged by all batch scores.

## Side notes

The doublet scores of Solo were calculated by its Linux command at https://github.com/calico/solo. The code in this folder read and merged the calculated scores to compute AUPRC and AUROC.
