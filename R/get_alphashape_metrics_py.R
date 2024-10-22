#' Calculate Crown Metrics for a Set of Tree Points in a LAS Object Using Python
#'
#' This function processes a LAS object or file path, filters out trees with fewer than 4 points, and computes crown attributes for the remaining tree points using Python's `scipy.spatial` module for convex hull computation. The function supports RGB-based chromatic metrics if the input data contains color channels.
#'
#' The function primarily uses the `get_crown_attributes_py()` function to compute crown height, size, and complexity metrics for each tree. If the LAS object contains fewer than 4 points for a tree, that tree is discarded from the analysis.
#'
#' @param chunk A LAS object or a file path to a LAS file containing tree point cloud data. The data must include tree IDs and (optionally) RGB color channels.
#' @param RGB Logical; if \code{TRUE}, the function computes chromatic coordinates (RCC, GCC, BCC) using the provided color channels (R, G, B). Default is \code{FALSE}.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item If \code{chunk} is a file path, it reads the LAS file.
#'   \item Filters out trees with fewer than 4 points, as convex hulls cannot be computed with fewer points.
#'   \item Computes crown metrics using \code{get_crown_attributes_py()}, including convex hull volume, crown height, size, complexity metrics, and (optionally) chromatic coordinates if RGB data is available.
#' }
#'
#' The Python environment for convex hull computation must be set up as described in the \code{get_crown_attributes_py} function documentation.
#'
#' @section Python Environment Setup:
#' This function requires Python with `numpy` and `scipy` installed in a conda environment. Set up a conda environment using the following commands:
#'
#' \code{conda create -p C:/crown_vol_env python=3.9 numpy scipy}
#'
#' Then, specify the path to the Python executable:
#'
#' \code{use_python("G:/calc_convex_volume/crown_vol_conda/python.exe", required = TRUE)}
#'
#' Ensure that the required Python packages (`numpy`, `scipy`) are installed in this environment.
#'
#' @return A data frame containing crown metrics for each tree, including:
#' \itemize{
#'   \item \strong{Zmax}: Maximum height of the tree points.
#'   \item \strong{Zq999}: 99.9th percentile of the height values.
#'   \item \strong{Z_mean}: Mean height of the tree points.
#'   \item \strong{n_points}: Number of points in the crown.
#'   \item \strong{vol_convex}: Volume of the convex hull around the points.
#'   \item \strong{CV_Z}: Coefficient of variation of the height (Z) values.
#'   \item \strong{CRR}: Crown ratio, defined as the ratio of the difference between mean and minimum height to the difference between maximum and minimum height.
#'   \item \strong{RCC_mean, GCC_mean, BCC_mean}: Mean relative chromatic coordinates (if \code{RGB = TRUE}).
#'   \item \strong{RCC_median, GCC_median, BCC_median}: Median relative chromatic coordinates (if \code{RGB = TRUE}).
#' }
#'
#' @note
#' Trees with fewer than 4 points are discarded since a 3D convex hull cannot be computed for such cases.
#'
#' @examples
#' \dontrun{
#' library(lidR)
#'
#' # Load and decimate a LAS file
#' las <- readLAS('G:/Ontario_2023/Block_18/LAS/DLS23/output/tree_las/B18_DLS23_472499.9_5341401_norm_knnidw_treeid.laz')
#' las_dec <- decimate_points(las, random(5))
#'
#' # Compute crown metrics using Python (without RGB channels)
#' tic('Python')
#' metrics_py <- get_alphashape_metrics_py(las_dec, RGB = FALSE)
#' toc()
#'
#' # Compute crown metrics with RGB channels
#' metrics_py_rgb <- get_alphashape_metrics_py(las_dec, RGB = TRUE)
#' print(metrics_py_rgb)
#' }
#'
#' @export
get_alphashape_metrics_py <- function(chunk, RGB = FALSE, np = np, scipy_spatial = scipy_spatial, QhullError = QhullError) {
  if ("LAS" %in% class(chunk)) {
    tree_las <- chunk
    message('Individual LAS object input into function')
  } else {
    tree_las <- lidR::readLAS(chunk)
  }

  if (lidR::is.empty(tree_las)) return(NULL)

  message('Beginning crown metric generation for chunk')

  print('Using python...')

  tree_las <- lidR::filter_duplicates(tree_las)

  if (RGB) {
    obs <- tree_las@data %>%
      dplyr::filter(!is.na(treeID)) %>%
      dplyr::group_by(treeID) %>%
      dplyr::summarise(n = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(n < 4)
  } else {
    obs <- tree_las@data %>%
      dplyr::filter(!is.na(treeID)) %>%
      dplyr::group_by(treeID) %>%
      dplyr::summarise(n = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(n < 4)
  }

  if (nrow(obs) > 0) {
    message(glue::glue('{nrow(obs)} treeIDs had fewer than 4 points and were discarded'))
  }

  # Without RGB
  mets <- tree_las@data %>%
    dplyr::filter(!is.na(treeID)) %>%
    dplyr::filter(!treeID %in% obs$treeID) %>%
    dplyr::group_by(treeID) %>%
    dplyr::do(get_crown_attributes_py(
      X = .$X, Y = .$Y, Z = .$Z,
      np = np, scipy_spatial = scipy_spatial, QhullError = QhullError
    ))

  # With RGB
  mets <- tree_las@data %>%
    dplyr::filter(!is.na(treeID)) %>%
    dplyr::filter(!treeID %in% obs$treeID) %>%
    dplyr::group_by(treeID) %>%
    dplyr::do(get_crown_attributes_py(
      X = .$X, Y = .$Y, Z = .$Z,
      R = .$R, G = .$G, B = .$B,
      np = np, scipy_spatial = scipy_spatial, QhullError = QhullError
    ))






  return(mets)
}













