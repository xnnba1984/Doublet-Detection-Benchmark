#  The scalability of doublet-detection methods 

This folder contains the code that measures how fast the running time of doublet-detection increases as the number of droplets grows. 25 synthetic scRNA-seq datasets were generated with the number of droplets ranging from 400 to 10,000. Then we applied each doublet-detection method to these datasets and recorded its running time.

## The functionality of each file

All doublet-detection methods' running time, except Solo, was calculated in all_method_scalability.R. The running time of Solo was record by its Linux command provided at https://github.com/calico/solo. 
