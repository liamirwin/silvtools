#
# February 15, 2023
# Liam Irwin
# Experimenting with alphashapes; novel metrics? speed increase?
#
# Crown transparency?
# Radius at different points along hull?

library(lidR)

las_list <- list.files('H:/Quesnel_2022/process/CT1/output/individual_trees', pattern = '.laz', full.names = T)
las_list <- las_list[10000:10100]



tree_df <- data.frame(treeID = c(), quality = c())

for(i in 1:length(las_list)){
  las <- readLAS(las_list[i])
  tree_df[i, 1] <- na.omit(unique(las@data$treeID))
  plot(las, color = 'RGB', bg = 'white')
  rgl::par3d(windowRect = c(10, 5, 1000, 1000))
  Sys.sleep(2)
  rgl::rgl.close()
  # 1 = bad, 2 = fine, 3 = best
  tree_df[i, 2] <- readline(prompt = 'segmentation quality? (1-3)')
  print(glue::glue('Tree {i} out of {1:length(las_list)}'))
}


readLAS("-thin_with_voxel 0.02")


gt <- tree_df %>% dplyr::filter(V2 == 3)
las <- readLAS(str_subset(las_list, pattern = '19119'))
las <- readLAS(str_subset(las_list, pattern = '19119'), filter = "-thin_with_voxel 0.02")
las <- readLAS(str_subset(las_list, pattern = '19120'), filter = "-thin_with_voxel 0.02")
las <- readLAS(str_subset(las_list, pattern = '19122'), filter = "-thin_with_voxel 0.02")
X <- las@data$X
Y <- las@data$Y
Z <- las@data$Z

get_crown_attributes
function(X,Y,Z,pb = NULL){
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

  #if(snapshot == TRUE){
    shape_concave = alphashape3d::ashape3d(x = a3d, alpha = Inf)
    plot(shape_concave)
    shape_convex = alphashape3d::ashape3d(x = a3d, alpha = 5)
    plot(shape_convex)
    shape_alf05 = alphashape3d::ashape3d(x = a3d, alpha = 0.5)
    plot(shape_alf05)

    rgl::rgl.snapshot(filename = glue::glue("{image_output}/{acq}_{unique(las@data$treeID)}_filt_alpha1_ashape.png"))
    Sys.sleep(1)
    rgl::rgl.close()
  #}



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
  if(!is.null(pb)){
    pb$tick()
  }
  return(df)
}


# For new values of alpha, we can use ashape3d.obj as input (faster)
alpha_test <- function(las){

  a3d <- las@data %>%
    dplyr::select(X,Y,Z) %>%
    distinct(X, Y, .keep_all = TRUE) %>%
    as.matrix()

  # get treetop location
  top_x <- a3d[which.max(a3d[,3]),1]
  top_y <- a3d[which.max(a3d[,3]),2]
  # normalize X and Y
  a3d[,1] = a3d[,1] - mean(a3d[,1]) #center points around 0,0,0
  a3d[,2] = a3d[,2] - mean(a3d[,2]) #center points around 0,0,0

  tictoc::tic()
  alpha <- c(0, 1, Inf)
  shape <- ashape3d(a3d, alpha = alpha)
  tictoc::toc()

  in3d <- inashape3d(shape, points = a3d)
  plot(shape, indexAlpha = 2, transparency = 0.2)
  colors <- ifelse(in3d, "green", "blue")
  points3d(a3d, col = colors)

return(shape)

}

alpha_test(las)

