# Benchmarking computational doublet-detection methods for single-cell RNA sequencing data

This repository contains the code and datasets used in the paper 'Benchmarking computational doublet-detection methods for single-cell RNA sequencing data'. The preprint can be found at https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3646565.

## Code Description

1. overall_real_stability: the code of detecting doublets in 16 real scRNA-seq datasets by 11 doublet-detection methods (including two baseline methods). Within each file, there are additional functions to reproduce stability evaluation.
2. overall_sim: the code of detecting doublets in synthetic datasets by 8 doublet-detection methods under various experimental conditions.
3. cluster: the code of measuring the effect of doublet-detection methods on downstream cell clustering analysis.
4. DE: the code of measuring the effect of doublet-detection methods on downstream DE gene analysis.
5. trajectory: the code of measuring the effect of doublet-detection methods on cell trajectory and pseudo time inference analysis.
6. distribute: the code of measuring the performance of doublet-detection methods under distributed computing.
7. time: the code of calcuating the running time of each doublet-detection method on 16 real scRNA-seq datasets.

## Dataset

1. 16 real scRNA-seq datasets with experimentally validated doublet annotations used in the benchmark are available at 
https://drive.google.com/drive/folders/1QEdDC_nxn9BvxypvwQV1lioE0d_aeAGI?usp=sharing.
2. Synthetic datasets used in the benchmark, including varying doublet rates, cell types, sequencing depth, heterogeneity, cell clustering, and DE gene ananlysis, are available at 
https://drive.google.com/drive/folders/19Yk35U1Onj6Y4-AuYuvwFYUrGL7ysiqJ?usp=sharing.
3. The datasets of cell trajectory inference are available in the code file to conduct the analysis.
