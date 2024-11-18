#' Setup standard directories for an LAS project
#'
#' This function creates a set of standardized directories for managing and organizing
#' Airborne LiDAR (ALS) or other LAS data and associated spatial data, including input, output,
#' documentation, and script folders.
#'
#' @param proj_dir A character string representing the path to the root directory of the ALS project.
#'
#' @return A list of created directories, invisibly.
#'
#' @examples
#' # Set up directories for an example ALS project
#' proj_path <- "path/to/your/project"
#' setup_als_dirs(proj_path)
#'
#' @export
setup_als_dirs <- function(proj_dir){

  # Generate standardized list of inputs and outputs for associated spatial data
  # Multilevel, /input, /output, /doc, /script
  # Doc intended to store relevant literature, site info, data collection info
  # Script intended to house R code and other code used to process input dataset

  targ_dirs <- file.path(proj_dir, c('input','output',
                                     'input/las',
                                     'input/las/norm',
                                     'input/las/class',
                                     'input/las/tile',
                                     'input/las/raw',
                                     'input/vector',
                                     'input/raster',
                                     'output/vector',
                                     'output/raster'))

  # Apply dir.create to generate standard folders

  lapply(targ_dirs, FUN = dir.create , showWarnings = FALSE)

}
