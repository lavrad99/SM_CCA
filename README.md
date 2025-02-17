# SM_CCA
Repository accompanying the paper "Canonical Correlation Analysis on Data With Structured Missingness" by Radosavljević et. al. (2025)

The folder ```R``` contains all R functions necessary for replicating all CCA estimations in the paper. The following are the most important functions:

* ```geigen_CCA.R```, which performs CCA via the generalised eigenvalue problem.

* ```mat_summary.R```, which summarises missingness patterns in the rows of input data. The summary object is used for fast implementation of covariance estimation from missing data usig the EM-algorithm.

* ```Block_EM.R```, which performs covariance estimation from missing data using the EM-algorithm.

* ```SM_CCA.R```, which performs the proposed method for CCA on data with structured missingness.

* ```MICE_CCA.R```, which performs CCA using Multiple Imputation by Chained Equations (MICE). This is the main competitor to our proposed method.


