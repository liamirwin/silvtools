
<!-- README.md is generated from README.Rmd. Please edit that file -->

# silvtools

<!-- badges: start -->

# silvtools <img src="man/figures/logo.png" align="right" height="277" />

<!-- badges: end -->

## Overview

The goal of silvtools is to provide a reproducible workflow for
tree-level analyses of fine-scale remote sensing datasets.

The package contains functions for a wide variety of tasks including
segmentation, tree matching, competition index calculation, crown
structural metric creation with alpha shapes, and solar simulation with
rayshader

Please note that silvtools is a relatively new package and its functions
are still in the development stage. As a result, issues may be common
and some functions may not yet have the level of stability and
functionality desired. Despite this, silvtools builds upon the work of
other excellent packages such as lidR and rayshader, and has the goal of
providing a reproducible workflow for tree-level analyses of fine-scale
remote sensing datasets. Nevertheless, caution should be taken when
using the package’s functions as they are not thoroughly tested and may
have bugs or limitations that need to be addressed.

## Installation :computer:

You can install the development version of silvtools from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("liamirwin/silvtools")
library(silvtools)
```
