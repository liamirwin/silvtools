

library(lidR)
las = readLAS(system.file("extdata", "MixedConifer.laz", package="lidR"))

# Current process involves loading in individual las files and converting to alphashape object
# to compute convex_hull volume, etc...
# What if instead we just grouped a chunk (or even block) las dataframe by treeID
# and created ashapes for each value; re-appending to the treetop location (max(z)[x,y])

convert_las_to_ashape_df <- function(las){





}

tree_las <- segment_trees(las, dalponte2016(chm = rasterize_canopy(las, res = 1, p2r()),
                                            treetops = locate_trees(las, lmf(ws = 5, hmin = 5))))

a_list <- list()

pb <- progress_bar$new(total = length(unique(tree_las@data$treeID)))

X <- x$X
Y <- x$Y
Z <- x$Z


shape = alphashape3d::ashape3d(x = a3d, alpha = 1)

get_crown_attributes <- function(X,Y,Z){
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
    vol_convex = alphashape3d::volume_ashape3d(alphashape3d::ashape3d(x = a3d, alpha = Inf)),
    vol_concave = alphashape3d::volume_ashape3d(alphashape3d::ashape3d(x = a3d, alpha = 1)),
    vol_a05= alphashape3d::volume_ashape3d(alphashape3d::ashape3d(x = a3d, alpha = .5)),
    # Crown complexity
    CV_Z = sd(Z) / mean(Z),
    CRR = (mean(Z) - min(Z)) / (max(Z) - min(Z)),
    X = top_x,
    Y = top_y)
  return(df)
}

trees <- las@data %>% as.data.frame() %>% dplyr::select(X,Y,Z,treeID) %>%
  dplyr::filter(!is.na(treeID)) %>% dplyr::group_by(treeID) %>%
  summarise(ashape_metrics = get_crown_attributes(X,Y,Z))


for(i in 1:length(na.omit(unique(tree_las@data$treeID)))){
  id <- na.omit(unique(tree_las@data$treeID))[i]
  tree <- filter_poi(tree_las, treeID == id)
  a_list[[i]] <- amet(tree)
  pb$tick()$print()
}

amet <- function(las, snapshot = FALSE, image_output = NULL, acq = NULL){
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

  top_x <- las@data %>% slice_max(Z, n = 1) %>% dplyr::select(X)
  top_y <- las@data %>% slice_max(Z, n = 1) %>% dplyr::select(Y)

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
    CRR = (mean(Z) - min(Z)) / (max(Z) - min(Z)),
    X = top_x,
    Y = top_y)

  df_sf <- sf::st_as_sf(df, coords = c('X','Y'), crs = sf::st_crs(las))

  return(df_sf)

}
