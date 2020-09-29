# Effects of doublet detection on DE gene analysis
## The name of the file indicates the doublet detection methods 

doubletcells_DE.R: doubletCells

doubletdetection_DE.R: DoubletDetection

doubletfinder_DE.R: DoubletFinder

scds_DE.R: cxds, scds, and hybrid

scrublet_DE.R: Scrublet

solo_DE.R: Solo

## The functionality of each file

Each file includes the following code blocks that implement different functionality. They are indicated by detailed comments within each file.

1. DE analysis on the dataset without doublets (positive control).

2. DE analysis on the dataset before doublet detection (negative control).

3. DE analysis on the dataset after doublet detection (the result to compare among different doublet-detection methods)

## Some notes

Only scrublet_DE.R includes all three functionality since positive control and negative control need only be calculated once. Other code files only include DE analysis after doublet detection by each method.

