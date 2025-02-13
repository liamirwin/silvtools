% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_results_polygons.R
\name{get_results_polygons}
\alias{get_results_polygons}
\title{Retrieve RESULTS Polygons Intersecting a Given Extent}
\usage{
get_results_polygons(extent, plot = TRUE, symbology_field = NULL)
}
\arguments{
\item{extent}{An \code{sf} object representing the spatial extent for the query.
It must be a simple feature or a simple feature collection, and will be
converted to a simple bounding box polygon. If the CRS is not BC Albers
(EPSG:3005), it will be transformed.}

\item{plot}{Logical, whether to plot the resulting polygons. Default is
\code{TRUE}.}

\item{symbology_field}{Character, the field to use for symbology in the plot.
If \code{NULL} (default), the function will attempt to use "FOREST_COVER_ID" if
it exists in the data. If not provided and "FOREST_COVER_ID" is not found,
an error is thrown. This field is converted to a factor for plotting.}
}
\value{
An \code{sf} object containing the RESULTS polygons that intersect the
given extent. If no polygons are found, an empty \code{sf} object is returned.
If \code{plot} is \code{TRUE}, a \code{ggplot} object is also displayed, showing the
polygons filled according to the \code{symbology_field}.
}
\description{
This function queries the BC Data Catalogue for RESULTS polygons that
intersect a specified spatial extent. It performs CRS validation and
transformation, retrieves data using \code{bcdata}, intersects the data with
the extent, and optionally plots the resulting polygons using \code{ggplot2}.
}
\details{
The function performs several steps:
\enumerate{
\item Converts the input \code{extent} to its bounding box.
\item Validates and transforms the CRS of the extent to BC Albers (EPSG:3005).
\item Queries the BC Data Catalogue for RESULTS polygons intersecting the extent.
\item Performs a spatial intersection between the retrieved polygons and the extent.
\item Optionally plots the resulting polygons using \code{ggplot2}, with symbology
based on the specified field.
}

The function uses \code{cli} for console output, providing informative messages,
warnings, and errors during execution.
}
\examples{
\dontrun{
library(sf)

# Example extent (replace with your actual extent)
extent <- st_read("path/to/your/extent.gpkg")

# Retrieve and plot RESULTS polygons
results <- get_results_polygons(extent)

# Retrieve without plotting, using a custom symbology field
results_no_plot <- get_results_polygons(extent, plot = FALSE, symbology_field = "MY_FIELD")
}

}
