# Doublet Detection Benchmark

This repository contains the code to generate the result in the paper 'Benchmarking computational doublet-detection methods for single-cell RNA sequencing data'.

## Code Description

1. overall_real_stability: This folder contains the code of detecting doublets in 16 real scRNA-seq datasets by 11 doublet-detection methods (including two baseline methods). In addition, within each file there are functions to reproduce stability evaluation.
2. overall_sim: This folder contains the code of detecting doublets in synthetic datasets by 8 doublet-detection methods under various experimental conditions.
3. cluster: This folder contains the code to reproduce the effect of doublet-detection methods on downstream cell clustering analysis.
4. DE: This folder contains the code to reproduce the effect of doublet-detection methods on downstream DE gene analysis.
5. trajectory: This folder contains the code to reproduce the effect of doublet-detection methods on cell trajectory and pseudo time inference analysis.
6. distribute: This folder contains the code to measure the performance of doublet-detection methods under distributed computing.
7. time: This folder contains the code to calcuate the running time of each doublet-detection method on 16 real scRNA-seq datasets.
