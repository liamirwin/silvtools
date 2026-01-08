#' Calculate crown attributes for a set of tree points
#'
#' This function calculates various crown attributes for a set of 3D tree points including:
#' crown height (Zmax, Zq999, Zq99, Z_mean), point structure metrics (n_points, CV_Z, CRR,
#' additional vertical distribution percentiles, skewness, and entropy), and a pzabove2 metric
#' (proportion of returns above 2m). If color channels (R, G, B) are provided, it computes the
#' relative chromatic coordinates (RCC, GCC, BCC), basic RGB statistics, and HSV metrics.
#'
#' @param X Numeric vector, x-coordinates of tree points.
#' @param Y Numeric vector, y-coordinates of tree points.
#' @param Z Numeric vector, z-coordinates of tree points.
#' @param R Numeric vector, red channel of tree points (optional).
#' @param G Numeric vector, green channel of tree points (optional).
#' @param B Numeric vector, blue channel of tree points (optional).
#' @param ReturnNumber Numeric vector, return number for each tree point (optional).
#' @param filter_return1 Logical, if TRUE and ReturnNumber is provided, only use points where ReturnNumber == 1 for spectral metrics. Defaults to FALSE.
#' @param filter_above_zmean Logical, if TRUE, only use points with Z greater than Z_mean for spectral metrics. Defaults to FALSE.
#'
#' @return A data frame with columns for each crown attribute.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' attributes <- get_crown_attributes_2025(X, Y, Z, R, G, B, ReturnNumber,
#'                                          filter_return1 = TRUE, filter_above_zmean = TRUE)
#' }
#'
#' @export
get_crown_attributes_2025 <- function(X, Y, Z, R = NULL, G = NULL, B = NULL,
                                      ReturnNumber = NULL, filter_return1 = FALSE,
                                      filter_above_zmean = FALSE) {

  # Helper function to compute standard deviation safely:
  safe_sd <- function(x, na.rm = FALSE) {
    s <- sd(x, na.rm = na.rm)
    if (is.na(s)) return(0) else return(s)
  }

  # Crown top (highest point)
  top_idx <- which.max(Z)
  top_x <- X[top_idx]
  top_y <- Y[top_idx]

  # Crown height metrics
  Zmax <- max(Z)
  Zq999 <- as.numeric(quantile(Z, 0.999, names = FALSE))
  Zq99  <- as.numeric(quantile(Z, 0.990, names = FALSE))
  Z_mean <- mean(Z)

  # Additional vertical distribution percentiles [Brandtberg (2007); DOI: 10.1016/j.rse.2006.08.004]
  Zq25 <- as.numeric(quantile(Z, 0.25, names = FALSE))
  Zq50 <- as.numeric(quantile(Z, 0.50, names = FALSE))  # median
  Zq75 <- as.numeric(quantile(Z, 0.75, names = FALSE))

  # Vertical distribution skewness (using third moment) [Li et al. 2013; DOI: 10.1016/j.rse.2012.11.016]
  if (length(Z) <= 1) {
    skewness_Z <- 0
  } else {
    skewness_Z <- mean((Z - Z_mean)^3) / (safe_sd(Z)^3)
  }

  # Shannon entropy of Z distribution (vertical complexity) [Li et al. 2013; DOI: 10.1016/j.rse.2012.11.016]
  shannon_entropy <- function(z, bins = 10) {
    hist_vals <- hist(z, breaks = bins, plot = FALSE)
    p <- hist_vals$counts / sum(hist_vals$counts)
    p <- p[p > 0]  # remove zeros
    entropy <- -sum(p * log(p))
    return(entropy)
  }
  Z_entropy <- shannon_entropy(Z, bins = 10)

  # Crown structure metrics using points only
  n_points <- length(Z)
  CV_Z <- if(mean(Z) == 0) 0 else safe_sd(Z) / mean(Z)
  CRR <- (Z_mean - min(Z)) / (Zmax - min(Z))  # Crown Relief Ratio

  # pzabove2 metric (proportion of returns above 2m)
  pzabove2 <- mean(Z > 2)

  df <- data.frame(
    Zmax       = Zmax,
    Zq999      = Zq999,
    Zq99       = Zq99,
    Z_mean     = Z_mean,
    Zq25       = Zq25,       # 25th percentile [Brandtberg (2007)]
    Zq50       = Zq50,       # 50th percentile (median) [Brandtberg (2007)]
    Zq75       = Zq75,       # 75th percentile [Brandtberg (2007)]
    skewness_Z = skewness_Z,  # Skewness [Li et al. 2013]
    Z_entropy  = Z_entropy,   # Shannon entropy [Li et al. 2013]
    n_points   = n_points,
    CV_Z       = CV_Z,
    CRR        = CRR,
    pzabove2   = pzabove2,    # Proportion of returns above 2m
    X          = top_x,
    Y          = top_y
  )

  # Optional color metrics
  if (!is.null(R) && !is.null(G) && !is.null(B)) {
    if (filter_return1 && !is.null(ReturnNumber)) {
      idx1 <- (ReturnNumber == 1)
      if (sum(idx1) == 0) {
        idx1 <- rep(TRUE, length(Z))
        warning("No points with ReturnNumber==1; using all points for color metrics.")
      }
    } else {
      idx1 <- rep(TRUE, length(Z))
    }

    # Additional filter: only use points above Z_mean
    if (filter_above_zmean) {
      idx_z <- (Z > Z_mean)
      if (sum(idx1 & idx_z) == 0) {
        warning("No points above Z_mean after applying filters; using all points for color metrics.")
      }
      idx1 <- idx1 & idx_z
    }

    R_sub <- R[idx1]
    G_sub <- G[idx1]
    B_sub <- B[idx1]

    total <- R_sub + G_sub + B_sub
    total[total == 0] <- NA

    RCC <- R_sub / total
    GCC <- G_sub / total
    BCC <- B_sub / total

    df$RCC_mean   <- mean(RCC, na.rm = TRUE)
    df$GCC_mean   <- mean(GCC, na.rm = TRUE)
    df$BCC_mean   <- mean(BCC, na.rm = TRUE)
    df$RCC_median <- median(RCC, na.rm = TRUE)
    df$GCC_median <- median(GCC, na.rm = TRUE)
    df$BCC_median <- median(BCC, na.rm = TRUE)

    # Additional color metrics in original RGB space [Onishi and Ise 2021; DOI: 10.1016/j.rse.2020.111000]
    df$R_mean <- mean(R_sub, na.rm = TRUE)
    df$G_mean <- mean(G_sub, na.rm = TRUE)
    df$B_mean <- mean(B_sub, na.rm = TRUE)
    df$R_sd   <- safe_sd(R_sub, na.rm = TRUE)
    df$G_sd   <- safe_sd(G_sub, na.rm = TRUE)
    df$B_sd   <- safe_sd(B_sub, na.rm = TRUE)

    # Determine the maxColorValue based on the input range for proper HSV conversion
    max_val <- max(c(R_sub, G_sub, B_sub), na.rm = TRUE)
    if (max_val <= 1) {
      max_color_value <- 1
    } else if (max_val <= 255) {
      max_color_value <- 255
    } else {
      max_color_value <- max_val  # fallback if values exceed 255
    }

    # Convert RGB to HSV using grDevices::rgb2hsv.
    hsv_mat <- grDevices::rgb2hsv(r = R_sub, g = G_sub, b = B_sub, maxColorValue = max_color_value)
    df$hue_mean        <- mean(hsv_mat["h", ], na.rm = TRUE)        # Mean hue [Onishi and Ise 2021]
    df$hue_sd          <- safe_sd(hsv_mat["h", ], na.rm = TRUE)       # Hue standard deviation [Onishi and Ise 2021]
    df$saturation_mean <- mean(hsv_mat["s", ], na.rm = TRUE)         # Mean saturation [Onishi and Ise 2021]
    df$saturation_sd   <- safe_sd(hsv_mat["s", ], na.rm = TRUE)       # Saturation standard deviation [Onishi and Ise 2021]
    df$brightness_mean <- mean(hsv_mat["v", ], na.rm = TRUE)         # Mean brightness (value) [Onishi and Ise 2021]
    df$brightness_sd   <- safe_sd(hsv_mat["v", ], na.rm = TRUE)       # Brightness standard deviation [Onishi and Ise 2021]
  }

  return(df)
}
