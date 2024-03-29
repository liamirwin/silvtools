% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_pixel_metrics.R
\name{generate_pixel_metrics}
\alias{generate_pixel_metrics}
\title{Generate lidar pixel metrics}
\usage{
generate_pixel_metrics(
  proj_dir,
  res = 1,
  metrics = c("basic"),
  zmin = NA,
  num_cores = 1L,
  chunk_buf = NULL,
  acq = NULL
)
}
\arguments{
\item{proj_dir}{Project directory containing the normalized LAS files.}

\item{res}{CHM resolution in meters. Default is 1.}

\item{metrics}{Sets of metrics to generate. Default is 'basic'. 'percentiles'}

\item{zmin}{Minimum Z value of points to consider in metrics. Default NA but EFI's often use 2m}

\item{num_cores}{Number of cores to use for parallel processing. Default is 1.}

\item{chunk_buf}{Chunk buffer size in meters. Default is 5\% of tile size.}

\item{acq}{Acquisition name. Default is NULL.}
}
\value{
Saves lidar pixel metrics locally within proj_dir
}
\description{
This function generates a pixel metrics from tiled LAS files.
}
