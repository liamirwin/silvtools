#' Calculate Crown Attributes for a Set of 3D Tree Points Using Python
#'
#' This function calculates various crown attributes for a set of 3D tree points by leveraging Python's `scipy.spatial` module for convex hull computation.
#' The attributes computed include crown height statistics (e.g., Zmax, Zq999, Z_mean), crown size (number of points, convex hull volume), crown complexity metrics (e.g., coefficient of variation of Z, crown ratio), and chromatic coordinates if color channels are provided.
#'
#' The convex hull computation is based on the `ConvexHull` function from the `scipy.spatial` module, which may raise an error if the input points are insufficient to define a 3D shape (fewer than 4 points).
#'
#' @param X Numeric vector representing the x-coordinates of the tree points.
#' @param Y Numeric vector representing the y-coordinates of the tree points.
#' @param Z Numeric vector representing the z-coordinates of the tree points.
#' @param R Numeric vector (optional) representing the red channel of the tree points.
#' @param G Numeric vector (optional) representing the green channel of the tree points.
#' @param B Numeric vector (optional) representing the blue channel of the tree points.
#'
#' @details
#' This function computes the following metrics:
#' \itemize{
#'   \item \strong{Zmax}: Maximum height of the tree points.
#'   \item \strong{Zq999}: 99.9th percentile of the height values.
#'   \item \strong{Zq99}: 99th percentile of the height values.
#'   \item \strong{Z_mean}: Mean height of the tree points.
#'   \item \strong{n_points}: Number of points in the crown.
#'   \item \strong{vol_convex}: Volume of the convex hull around the points.
#'   \item \strong{CV_Z}: Coefficient of variation of the height (Z) values.
#'   \item \strong{CRR}: Crown ratio, defined as the ratio of the difference between mean and minimum height to the difference between maximum and minimum height.
#'   \item \strong{RCC_mean, GCC_mean, BCC_mean}: Mean relative chromatic coordinates if RGB channels are provided.
#'   \item \strong{RCC_median, GCC_median, BCC_median}: Median relative chromatic coordinates if RGB channels are provided.
#' }
#'
#' @section Python Environment Setup:
#' This function uses Python to compute the convex hull volume. To run this function, ensure that a Python environment is set up with the necessary packages (`numpy`, `scipy`). You can set up a conda environment with these dependencies by running the following commands:
#'
#' \code{conda create -p C:/crown_vol_env python=3.9 numpy scipy}
#'
#' After creating the environment, specify the path to the Python executable:
#'
#' \code{use_python("G:/calc_convex_volume/crown_vol_conda/python.exe", required = TRUE)}
#'
#' Ensure that the Python packages `numpy` and `scipy` are installed in this environment.
#'
#' @return A data frame with the calculated crown attributes. If RGB channels are provided, the data frame will also include chromatic coordinates.
#'
#' If fewer than 4 points are provided for the convex hull computation, a warning will be issued, and \code{NULL} will be returned.
#'
#' @note
#' The Python-based convex hull computation can fail if the input points are insufficient to create a valid convex shape. In such cases, a warning is issued, and the function returns \code{NA} for the convex hull volume.
#'
#' @examples
#' \dontrun{
#' # Miniconda path with numpy and scipy G:/calc_convex_volume/crown_vol_conda
#' conda_path <- "G:/calc_convex_volume/crown_vol_conda/python.exe"
#' # Specify the path to your Python executable if necessary
#' use_python(conda_path, required = TRUE)
#' np <- import("numpy")
#' scipy_spatial <- import("scipy.spatial")
#' QhullError <- import("scipy.spatial.qhull", convert = FALSE)$QhullError
#' }
#'
#' @export
get_crown_attributes_py <- function(X, Y, Z, R = NULL, G = NULL, B = NULL, np, scipy_spatial, QhullError) {
  if (length(X) < 4) {
    warning('Cannot compute a 3D convex hull from fewer than 4 points')
    return(NULL)
  }

  # Prepare the points as a numpy array
  points <- np$array(np$c_[X, Y, Z])


  # Try computing the convex hull
  vol_convex <- tryCatch({
    hull <- scipy_spatial$ConvexHull(points)
    hull$volume
  }, error = function(e) {
    if (inherits(e, "python.builtin.Exception") && e$`__class__`$`__name__` == "QhullError") {
      warning('Convex hull computation failed for this tree')
      return(NA)
    } else {
      stop(e)
    }
  })


  # Compute the centroid of the crown
  centroid_x <- mean(X)
  centroid_y <- mean(Y)
  centroid_z <- mean(Z)

  # Update the data frame with centroid coordinates
  df <- data.frame(
    # Crown height
    Zmax = max(Z),
    Zq999 = quantile(Z, 0.999),
    Zq99 = quantile(Z, 0.990),
    Z_mean = mean(Z),
    # Crown size
    n_points = length(Z),
    vol_convex = vol_convex,
    # Crown complexity
    CV_Z = sd(Z) / mean(Z),
    CRR = (mean(Z) - min(Z)) / (max(Z) - min(Z)),
    # Use centroid coordinates for consistent positioning
    X = centroid_x,
    Y = centroid_y,
    Z = centroid_z
  )


  # If RGB values are provided, compute chromatic coordinates
  if (!is.null(R) && !is.null(G) && !is.null(B)) {
    # Avoid division by zero
    sum_rgb <- R + G + B
    sum_rgb[sum_rgb == 0] <- NA
    RCC <- R / sum_rgb
    GCC <- G / sum_rgb
    BCC <- B / sum_rgb

    # Add chromatic coordinates stats to the data frame
    df$RCC_mean <- mean(RCC, na.rm = TRUE)
    df$GCC_mean <- mean(GCC, na.rm = TRUE)
    df$BCC_mean <- mean(BCC, na.rm = TRUE)
    df$RCC_median <- median(RCC, na.rm = TRUE)
    df$GCC_median <- median(GCC, na.rm = TRUE)
    df$BCC_median <- median(BCC, na.rm = TRUE)
  }

  return(df)
}
