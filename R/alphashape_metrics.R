#' Generate Crown Structural Metrics
#'
#' @param las single tree las object
#' @param snapshot Logical indicating whether to generate snapshots of the 3D alphashape (default: `FALSE`)
#' @param image_output Path to directory where the snapshots will be saved (default: `NULL`)
#' @param acq Identifier to include in the filename of the snapshot (default: `NULL`)
#'
#' @return returns an sf object with 3D convex hull metrics
#' @export
#'
#' @examples
#' \dontrun{
#' las <- lidR::filter_poi(tree_las, treeID == 17)
#' ashape <- alphashape_metrics(las)
#' }
alphashape_metrics <- function(las, snapshot = FALSE, image_output = NULL, acq = NULL){
  Z = las@data$Z

  chm = lidR::pixel_metrics(las, func = ~max(Z), res = 0.05)
  #chm_mean = lidR::pixel_metrics(las, func = ~mean(Z), res = 0.1)
  #chm_mean[chm_mean < (as.numeric(terra::global(chm_mean, fun = 'max', na.rm = T)))] = NA
  #chm_mean_trim = terra::trim(chm_mean)

  # alphashadep3d
  a3d = cbind(las@data$X, las@data$Y, las@data$Z)

  a3d[,1] = a3d[,1] - mean(a3d[,1]) #center points around 0,0,0
  a3d[,2] = a3d[,2] - mean(a3d[,2]) #center points around 0,0,0

  if(snapshot == TRUE){
  shape = alphashape3d::ashape3d(x = a3d, alpha = 1)
  plot(shape)
  rgl::rgl.snapshot(filename = glue::glue("{image_output}/{acq}_{unique(las@data$treeID)}_alpha1_ashape.png"))
  Sys.sleep(1)
  rgl::rgl.close()
  }

  # Remove massive crown volume changes when ground points exist?
  # Filter only points in top 50% of tree
  las_filt <- as.data.frame(las@data) %>% filter(Z > max(Z) * 0.5)

  a3d_filt = cbind(las_filt$X, las_filt$Y, las_filt$Z)

  a3d_filt[,1] = a3d_filt[,1] - mean(a3d_filt[,1]) #center points around 0,0,0
  a3d_filt[,2] = a3d_filt[,2] - mean(a3d_filt[,2]) #center points around 0,0,0

  if(snapshot == TRUE){
    shape = alphashape3d::ashape3d(x = a3d_filt, alpha = 1)
    plot(shape)
    rgl::rgl.snapshot(filename = glue::glue("{image_output}/{acq}_{unique(las@data$treeID)}_filt_alpha1_ashape.png"))
    Sys.sleep(1)
    rgl::rgl.close()
  }

  top_x <- las_filt %>% slice_max(Z, n = 1) %>% select(X)
  top_y <- las_filt %>% slice_max(Z, n = 1) %>% select(Y)

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
    # Filtered crown metrics
    vol_convex_filt = alphashape3d::volume_ashape3d(alphashape3d::ashape3d(x = a3d_filt, alpha = Inf)),
    vol_concave_filt = alphashape3d::volume_ashape3d(alphashape3d::ashape3d(x = a3d_filt, alpha = 1)),
    vol_a05_filt= alphashape3d::volume_ashape3d(alphashape3d::ashape3d(x = a3d_filt, alpha = .5)),
    # Crown complexity
    CV_Z = sd(Z) / mean(Z),
    rumple = lidR::rumple_index(chm),
    CRR = (mean(Z) - min(Z)) / (max(Z) - min(Z)),
    X = top_x,
    Y = top_y)

  df_sf <- sf::st_as_sf(df, coords = c('X','Y'), crs = sf::st_crs(las))

  return(df_sf)

}
