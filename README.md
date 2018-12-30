Matlab code for
"Nonparametric Bayes Models of Fiber Curves Connecting Brain Regions"
==============
This Matlab code is to reproduce some key results presented in the paper "Nonparametric Bayes Models of Fiber Curves Connecting Brain Regions".

The data used in the paper are from two datasets: (a) a test-retest dataset collected at the University of Sherbrooke; (b) the Human Connectome Project Dataset. Since we mainly focus our analysis on brain structural connectivity, only dMRI data in these datasets are used. The whole processed dMRI (after fiber reconstruction) take a lot of disk space. Therefore, we are not able to upload all the original data. After preprocessing using the method described in the paper, i.e., variance decomposition, we dramatically reduced the data size. We release some of the processed data (after variance decomposition) in .mat format, which can be loaded to MATLAB directly.

If you use this code, please cite:

```
Z. Zhang, M. Descoteaux, D. Dunson,  "Nonparametric Bayes Models of Fiber Curves Connecting Brain Regions",arXiv:1612.01014, 2019
```

Folders and files description
-----
**curve_processing_functions**
This folder contains fiber curve preprocessing functions - i.e. to decompose the fibers into the rotation, translation, scaling and reparameterization parts.

**data**
This folder contains most of the data used in the paper

**distributions**
This folder contains all functions to manipulate different distributions, e.g, randomly sample data from different distributions. 

**mcmc_sampler**
This folder contains different MCMC samplers developed by this paper

**MCMCDiag**
This contains functions to perform MCMC diagnosis 

**utilities**
This folder contains some functions (e.g., for plotting, writing a logger file and so on) that are useful for the main Matlab files. 

*main_simulation.m* This Matlab code is for the simulation study (variation decomposition with simulated data) presented in "Nonparametric Bayes Models of Fiber Curves Connecting Brain Regions". It will generate subplots in Figure 3 in the paper.

*main_fibers_in_one_subject.m* This main file demonstrates the process of modeling fiber curves between two regions of interest in one subject. It can be used to reproduce results presented in Figure 7 in the main paper (with some modification, it can also help to produce Figure 2, 3 and 4 in the Supplement), and Table 1.

*main_fibers_across_3subject.m* This main file demonstrates the process of modeling fiber curves between two regions of interest in multiple subjects. It is an implementation of MCMC sampling algorithm in Section 4.2 of the main paper. We use 3 subjects' data from the Sherbrooke test-retest dataset. Results of Figure 8 and 9, and Table 1 in Supplement can be reproduced with this file.

*main_fibers_across_5subjects.m* This main file demonstrates the process of modeling fiber curves between two regions of interest in multiple subjects. It is an implementation of MCMC sampling algorithm in Section 4.2 of the main paper. We use 5 subjects' data from the Sherbrooke test-retest dataset. Results of Table 2 in the main paper can be reproduced here. 

*main_fibers_across_20hcp_subjects.m* Similar to *main_fibers_across_3subject.m* and *main_fibers_across_5subjects.m* for modeling data from multiple subjects, here we use 20 subjects' data from the HCP dataset. Among the 20 subjects, 10 have very low reading scores with a mean of 67.4 (sd = 3.3) and 10 have very high scores with a mean of 134.9 (sd = 2.6 ) (Refer to Supplementary
Section 10 for more details on how these subjects were selected).



Usage
-----
Use Matlab to run *main_simulation.m*, *main_fibers_in_one_subject.m*, *main_fibers_across_3subject.m*, *main_fibers_across_5subjects.m* and *main_fibers_across_20hcp_subjects.m*