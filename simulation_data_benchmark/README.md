# Benchmark on simulation datasets under different doublet rates, cell types, sequencing depth, and heterogeneity between cell types
## The name of the file indicates the doublet detection methods

doubletFinder_sim.R: DoubletFinder

doubletcells_sim.R: doubletCells

doubletdetection_sim.R: DoubletDetection

scds_sim.R: cxds, bcds, and hybrid

scrublet_sim.R: Scrublet

solo_sim.R: Solo

## The functionality of each file

Each file includes the following code blocks that implement different functionality. They are indicated by comments within each file.

1. Calculate doublet scores, AUPRC, and AUROC on simulation datasets under different doublet rates, cell types, sequencing depth, and heterogeneity between cell types.

2. Impact of doublet detection on the identification of highly variable genes (HVGs). The Jaccard index between clean HVGs and contaminated HGVs, between clean HVGs and post-doublet-detection HVGs, were calculated to measure the impact.

## Some notes
### Solo
The double scores, AUPRC, and AUROC of Solo were calculated by Linux commands, as shown at https://github.com/calico/solo.
