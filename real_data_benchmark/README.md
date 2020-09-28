# Benchmark on 16 real datasets with experimentally annotated doublets
## The name of the file indicates the doublet detection methods

doubeltdecon_real.R: DoubletDecon

doubletFinder_real.R: DoubletFinder

doubletcells_real.R: doubletCells

doubletdetection_real.R: DoubletDetection

lsize_ngene_real.R: two baseline methods lsize and ngene

scds_real.R: cxds, bcds, and hybrid

scrublet_real.R: Scrublet

solo_real.R: Solo

## The functionality of each file

Each file includes the following code blocks that implement different functionality. They are indicated by comments within each file.

1. Calculate doublet scores, auprc, and auroc on 16 benchmark datasets.

2. Precision, recall, and true negative rate (TNR) under 10%, 20%, and 40% identification rates.

3. Precision, recall, and true TNR under the thresholds determined by DoubletDecon.

4. Running time on 16 benchmark datasets

5. Stability measured by running on subsamples of 2 real datasets.

## Some notes
### lsize and ngene
The code of these two baseline methods does not contain the functionality of running time or stability. These two items are not compared in the manuscript.
### DoubetDecon
The code of DoubletDecon only calculated precision, recall, TNR, and running time because it provided doublet prediction without doublet scores.
### Solo
The double scores, running time, and stability of Solo were calculated by Linux commands, as shown at https://github.com/calico/solo.
