#' Calculate crown attributes for a set of tree points
#'
#' This function calculates various crown attributes for a set of 3D tree points including:
#' crown height (Zmax, Zq999, Zq99, Z_mean), crown size (n_points, vol_convex, vol_concave, vol_a05), and crown complexity (CV_Z, CRR)
#'
#' @param X Numeric vector, x-coordinates of tree points
#' @param Y Numeric vector, y-coordinates of tree points
#' @param Z Numeric vector, z-coordinates of tree points
#' @param pb progress bar object, used for progress tracking (optional)
#'
#' @return A data frame with columns for each crown attribute
#'
#' @examples
#' \dontrun{
#' # For example usage, see the parent function get_alphashape_metrics()
#' }
#'
#' @export
#'
get_crown_attributes <- function(X,Y,Z,pb = NULL){
  if (length(X) <= 3 || length(Y) <= 3 || length(Z) <= 3) {
    print('Cannot compute a 3D hull from 3 or fewer points')
    return(NULL)
  }
  # alphashadep3d
  a3d = cbind(X, Y, Z)
  # get treetop location
  top_x <- a3d[which.max(a3d[,3]),1]
  top_y <- a3d[which.max(a3d[,3]),2]
  # normalize X and Y
  a3d[,1] = a3d[,1] - mean(a3d[,1]) #center points around 0,0,0
  a3d[,2] = a3d[,2] - mean(a3d[,2]) #center points around 0,0,0
  alpha <- c(Inf, 1)
  ashape <- alphashape3d::ashape3d(x = a3d, alpha = alpha, pert = TRUE)
  # Return Number Distribution Metrics (function of crown transparency)
  # if(length(Z) == length(ReturnNumber)){
  # ZRN = cbind(Z, ReturnNumber)
  # Z1 <- subset.matrix(ZRN, ReturnNumber == 1)
  # Z2 <- subset.matrix(ZRN, ReturnNumber == 2)
  # Z3 <- subset.matrix(ZRN, ReturnNumber == 3)
  # } else{
  #   Z1 <- NA
  #   Z2 <- NA
  #   Z3 <- NA
  #   print('Return Number and Z length not equal...')
  # }
  # calculate crown metrics
  df <- data.frame(
    # Crown height
    Zmax = max(Z),
    Zq999 = as.numeric(quantile(Z, 0.999)),
    Zq99 = as.numeric(quantile(Z, 0.990)),
    Z_mean = mean(Z),
    # Crown size
    n_points = length(Z),
    #volume = raster::cellStats(chm_vol, sum),
    #n_green_points = length(las_green@data$Z),
    vol_concave = alphashape3d::volume_ashape3d(ashape, indexAlpha = 1),
    vol_convex = alphashape3d::volume_ashape3d(ashape, indexAlpha = 2),
    # Crown complexity
    CV_Z = sd(Z) / mean(Z),
    CRR = (mean(Z) - min(Z)) / (max(Z) - min(Z)),
    # Crown Transparency
    #Rn_1_Z_mean = ifelse(!is.na(Z1[1]), mean(Z1[,1]), NA),
    #Rn_2_Z_mean = ifelse(!is.na(Z2[1]), mean(Z2[,1]), NA),
    #Rn_3_Z_mean = ifelse(!is.na(Z3[1]), mean(Z3[,1]), NA),
    #nmbr_rn_1 = nrow(Z1),
    #nmbr_rn_2 = nrow(Z2),
    #nmbr_rn_3 = nrow(Z3),
    # Append georeferenced tree top coordinates
    X = top_x,
    Y = top_y)
  if(!is.null(pb)){
  pb$tick()
  }
  return(df)
}
