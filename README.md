# SM_CCA
Repository accompanying the paper "Canonical Correlation Analysis on Data With Structured Missingness" by Radosavljević et. al. (2025)

The file ```SM_CCA.R``` contains all R functions necessary for replicating all CCA estimations in the paper. The following are the most important functions:

* ```mat_summary.R```, which summarises missingness patterns in the rows of input data. The summary object is used for fast implementation of covariance estimation from missing data usig the EM-algorithm.

* ```Block_EM.R```, which performs covariance estimation from missing data using the EM-algorithm.


