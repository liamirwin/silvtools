#' Generate Crown Structural Metrics
#'
#' @param las single tree las object
#'
#' @return returns a dataframe with 3D convex hull metrics
#' @export
#'
#' @examples
#' \dontrun{
#' las <- lidR::filter_poi(tree_las, treeID == 17)
#' ashape <- alphashape_metrics(las)
#' }
alphashape_metrics <- function(las){
  Z = las@data$Z

  chm = lidR::pixel_metrics(las, func = ~max(Z), res = 0.05)
  chm_mean = lidR::pixel_metrics(las, func = ~mean(Z), res = 0.1)
  chm_mean[chm_mean < (as.numeric(terra::global(chm_mean, fun = 'max', na.rm = T)))] = NA

  chm_mean_trim = terra::trim(chm_mean)

  # alphashadep3d
  a3d = cbind(las@data$X, las@data$Y, las@data$Z)

  a3d[,1] = a3d[,1] - mean(a3d[,1]) #center points around 0,0,0
  a3d[,2] = a3d[,2] - mean(a3d[,2]) #center points around 0,0,0

  #shape = alphashape3d::ashape3d(x = a3d, alpha = 1)
  #plot(shape)

  df <- data.frame(

    # Crown height
    Zq999 = as.numeric(quantile(Z, 0.999)),
    Zq99 = as.numeric(quantile(Z, 0.990)),
    Z_mean = mean(Z),

    # Crown size
    n_points = length(las@data$Z),
    area = lidR::area(las),
    #volume = raster::cellStats(chm_vol, sum),
    #n_green_points = length(las_green@data$Z),
    vol_convex = alphashape3d::volume_ashape3d(alphashape3d::ashape3d(x = a3d, alpha = Inf)),
    vol_concave = alphashape3d::volume_ashape3d(alphashape3d::ashape3d(x = a3d, alpha = 1)),
    vol_a05= alphashape3d::volume_ashape3d(alphashape3d::ashape3d(x = a3d, alpha = .5)),

    # Crown complexity
    CV_Z = sd(Z) / mean(Z),
    rumple = lidR::rumple_index(chm),
    CRR = (mean(Z) - min(Z)) / (max(Z) - min(Z)))

  return(df)

}
