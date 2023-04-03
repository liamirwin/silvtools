

#' @param X Numeric vector, x-coordinates of tree points
#' @param Y Numeric vector, y-coordinates of tree points
#' @param Z Numeric vector, z-coordinates of tree points
#' @param pb progress bar object, used for progress tracking (optional)
library(lidR)
library(alphashape3d)

tree_dir <- 'H:/Quesnel_2022/process/CT1/output/individual_trees'

tree_list <- list.files(tree_dir, pattern = '.laz', full.names = T)

snapshot_dir <- glue::glue('{tree_dir}/snapshots')

if(!dir.exists(snapshot_dir)){
  dir.create(snapshot_dir)
}

analyze_ashape <- function(las_file, decimate = FALSE){

  las <- lidR::readLAS(las_file)

  if(decimate == TRUE){
    n1 <- length(las@data$Z)
    las <- lidR::decimate_points(las, random_per_voxel(res = 0.25, n = 1))
    n2 <- length(las@data$Z)
    #print(glue::glue('LAS decimated from {n1} to {n2} returns - kept {round(((n2/n1) * 100), 2)}%'))
  }

  X <- las@data$X
  Y <- las@data$Y
  Z <- las@data$Z
  ID <- unique(las@data$treeID)


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
  ashape <- alphashape3d::ashape3d(x = a3d, alpha = alpha)


  # calculate crown metrics
  df <- data.frame(
    # Crown height
    Zmax = max(Z),
    Zq999 = as.numeric(quantile(Z, 0.999)),
    Zq99 = as.numeric(quantile(Z, 0.990)),
    Z_mean = mean(Z),
    # Crown size
    n_points = length(Z),
    vol_concave = alphashape3d::volume_ashape3d(ashape, indexAlpha = 1),
    vol_convex = alphashape3d::volume_ashape3d(ashape, indexAlpha = 2),
    # Crown complexity
    CV_Z = sd(Z) / mean(Z),
    CRR = (mean(Z) - min(Z)) / (max(Z) - min(Z)),
    # Append georeferenced tree top coordinates
    X = top_x,
    Y = top_y)

  # Plot Concave
  plot(ashape, indexAlpha = 1, transparency = 0.4, axes = FALSE)
  # Add Axes
  axes3d(
    edges=c('x--', 'y+-', 'z--'),
    labels=T
  )
  # Add Points
  points3d(a3d, color = 'black')
  in3d <- inashape3d(ashape, points = a3d)
  # Rotate Plot
  delta_angle <- 2
  angle <- rep(delta_angle * pi / 8, 8 / delta_angle)[1]
  i = 1
  while(i < 8) {
    view3d(userMatrix = rotate3d(par3d("userMatrix"), angle, 0, 0, 1))
    i = i + 1
  }
  # Add Title
  title3d(main = glue::glue('ID: {ID} - v_conc: {round(df$vol_concave)} v_conv: {round(df$vol_convex)} n: {df$n_points} CRR: {round(df$CRR, 3)}'))
  par3d(windowRect = c(0, 0, 600, 600))
  rgl.snapshot(filename = glue::glue('{snapshot_dir}/ashape_concave_{ID}.png'))
  close3d()
  #print(glue::glue('Took snapshot of {ID}'))

  return(df)
}

ashape_list <- list()

for(k in 1:length(tree_list)){
  las_file <- tree_list[k]
  ashape_list[[k]] <- analyze_ashape(las_file, decimate = FALSE)
}



# Extract Individual Tree Las Files

tree_las_list <- list.files('H:/Quesnel_2022/process/CT1/output/tree_las')
