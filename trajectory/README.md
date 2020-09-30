# Effects of doublet detection on cell trajectory inference

This folder contains the benchmark of doublet-detection methods in cell trajectory inference. Two simulation datasets, one bifurcating and one sequential, were used. Slingshot and minimum spanning tree (MST) were used to construct cell trajectories. Slingshot, TSCAN, and general additive model were used to infer the temporally expressed genes.

## Cell trajectory inference and temporally expressed genes analysis conducted by Slingshot

doubletcells_trajectory_slingshot.R: doubletCells

doubletdetection_trajectory_slingshot.R: DoubletDetection

doubletfinder_trajectory_slingshot.R: DoubletFinder

scds_trajectory_slingshot.R: cxds, bcds, and hybrid

scrublet_trajectory_slingshot.R: Scrublet

solo_trajectory_slingshot.R: Solo

## Cell trajectory inference conducted by MST and GAM

MST.R: all eight doublet detection methods are included in this file, including negative and positive control

## Temporally expressed genes analysis conducted by TSCAN and GAM

TSCAN.R: all eight doublet detection methods are included in this file, including negative and positive control

## The functionality of each file
### Slingshot

1. Cell trajectory inference result by Slingshot, after applying each doublet-detection method to remove doublets.

2. Temporally expressed gene analysis by Slingshot and GAM, after applying each doublet-detection method to remove doublets.

### MST

Cell trajectory inference result by MST on clean data without doublets, contaminated data with doublets, and the data after applying each doublet-detection method to remove doublets.

### TSCAN

Temporally expressed gene analysis by TSCAN and GAM, on clean data without doublets, contaminated data with doublets, and the data after applying each doublet-detection method to remove doublets.


