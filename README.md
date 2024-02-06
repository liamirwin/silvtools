
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

## Example - Crown segmentation and metric attribution

``` r

las <- readLAS(system.file("extdata", "uls.laz", package = "silvtools"))

tree_las <- segment_trees(las, mcwatershed(chm = rasterize_canopy(las, res = 1, p2r()),
treetops = locate_trees(las, lmf(ws = 5, hmin = 5))))

# Get alphashape metrics produces convex hulls for each treeID and generates crown volumes as well as other structural metrics
ashape_df <- get_alphashape_metrics(tree_las, prog_bar = TRUE)
```

## Acknowledgements

Thank you to extensive help and code contributions from many members of
the UBC Integrated Remote Sensing Studio (IRSS) including but not
limited to; Martin Queinnec, Samuel Grubinger, Rik Nuijten. See [IRSS
LAB](https://irsslab.forestry.ubc.ca) and [IRSS
Github](https://github.com/IRSS-UBC) for more.
