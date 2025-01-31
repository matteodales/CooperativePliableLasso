# CooperativePliableLasso
Code relative to the article "Integrating Multiple Data Sources with Interactions in Multi-Omics Using Cooperative Learning".

- File **functions.r** contains the implementation of the algorithms in the article.
- Folder **simulation_studies** contains the code used to generate the results in the simulation studies section.
- Folder **real_multiomics_studies** contains the code use to analyse the real data.

### Installation Instructions

Since the <tt>pliable</tt> package is no longer available on CRAN, a workaround is needed for its installation.

#### Installing the Pliable Package via Workaround

You can install the <tt>pliable</tt> package for your work using the following workaround:

- **Step 1:** Install an older version of R (R version 4.2 or earlier, can be found [here](https://cran.r-project.org/bin/windows/base/old/))and switch to it on RStudio from Global Options.
- **Step 2:** Install an older version of GCC (GCC 8.4 for Windows and gfortran 8.2.0 for macOS) which can be found [here](https://ftp.gnu.org/gnu/gcc/). This is necessary because the <tt>pliable</tt> package relies on Fortran subroutines that are compatible with these versions of GCC.
- **Step 3:** Download pliable version 1.1.1 from the [CRAN archive](https://cran.r-project.org/src/contrib/Archive/pliable/)

Once you've set up the environment with the older versions of R and GCC, you can install the <tt>pliable</tt> package using the standard R installation command:

```R
install.packages("pliable_1.1.1.tar.gz", repos = NULL, type="source")
```

#### Using the svreg Package (Alternative)
Otherwise, the alternative implementation of pliable lasso in the <tt>svreg</tt> package by Kim et al. can be used. This package is fully implemented in R, and although it may increase computational time due to the lack of Fortran subroutines, it provides a valid alternative. The package can be installed from the available GitHub through the command

```R
devtools::install_github("Tanya-Garcia-Lab/svreg")
```
