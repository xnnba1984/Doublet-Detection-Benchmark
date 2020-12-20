# Benchmarking computational doublet-detection methods for single-cell RNA sequencing data

This repository contains the code and datasets used in the paper 'Benchmarking computational doublet-detection methods for single-cell RNA sequencing data'. The full text published on *Cell Systems* is available at https://www.sciencedirect.com/science/article/pii/S2405471220304592?dgcid=author. The accepted preprint is available at http://jsb.ucla.edu/sites/default/files/publications/final_author_manuscript_compressed.pdf.

## Code Structure

Each folder contains different benchmark studies described in the paper. There are detailed documents within each fold to further explain their functionality. 

1. real_data_benchmark: the benchmark of 11 doublet-detection methods on 16 real scRNA-seq datasets, including two baseline methods. In addition, accuracy under different identification rates, running time, and stability are also included.

2. simulation_data_benchmark: the benchmark of 8 doublet-detection methods on simulated datasets under different experimental conditions, including doublet rates, cell types, sequencing depth, and heterogeneity between cell types. The impact of doublet detection on the identification of highly expressed genes (HVGs) is also included.

3. clustering: effects of doublet-detection methods on downstream cell clustering analysis.

4. DE: effects of doublet-detection methods on downstream DE genes analysis.

5. trajectory: effects of doublet-detection methods on cell trajectory and temporally expressed gene analysis.

6. distribute computing: performance of doublet-detection methods under distributed computing.

7. time: scalability of doublet-detection methods.

8. graphics_paper.R: code scripts used to generate the figures reported in the paper.

## Dataset

1. 16 real scRNA-seq datasets with experimentally validated doublet annotations used in the benchmark are available at https://drive.google.com/drive/folders/1QEdDC_nxn9BvxypvwQV1lioE0d_aeAGI?usp=sharing.

2. Simulation datasets with ground truth doublets used in the benchmark are available at https://drive.google.com/drive/folders/19Yk35U1Onj6Y4-AuYuvwFYUrGL7ysiqJ?usp=sharing. It includes datasets under varying doublet rates, cell types, sequencing depth, heterogeneity, and datasets used in the downstream analysis of cell clustering, DE gene, and cell trajectory inference.

3. All the datasets are also available at Zenodo https://zenodo.org/record/4062232#.X3YR9Hn0kuU.
