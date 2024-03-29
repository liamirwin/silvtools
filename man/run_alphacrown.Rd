% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_alphacrown.R
\name{run_alphacrown}
\alias{run_alphacrown}
\title{Run the Alphacrown Shiny App}
\usage{
run_alphacrown()
}
\value{
No return value. This function runs a Shiny app.
}
\description{
This function starts the Alphacrown Shiny App, which allows users to analyze tree crowns using the specified folder path and alpha value.
Simply input a folder containing individual tree las files and the app will allow you to visualize different trees with different alpha values and provide a set of summary metrics.
}
\examples{
\dontrun{

library(silvtools)
# Get the example dataset path
example_data_path <- system.file("extdata", "individual_trees", package = "silvtools")

# Launch the Shiny app with the example dataset
run_alphacrown(example_data_path)
}

}
