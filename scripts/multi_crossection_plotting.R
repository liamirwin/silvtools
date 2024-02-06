# Clip DLS23, SPL20, and SPL18 to same extent, visualize cross sections


library(lidR)

# Read in the data

DLS23 <- catalog('G:/Ontario_2023/Block_18/NW/input/las/norm')
SPL20 <- catalog('G:/Ontario_2023/Block_18/LAS/SPL18/input/las/norm')
SPL18 <- catalog('G:/Ontario_2023/Block_18/LAS/SPL20/input/las/norm')

roi <- st_read('G:/Ontario_2023/Block_18/scratch/figure_extent.shp')

roi23 <- clip_roi(DLS23, roi)
roi20 <- clip_roi(SPL20, roi)
roi18 <- clip_roi(SPL18, roi)

writeLAS(roi23, 'G:/Ontario_2023/Block_18/scratch/roi23.laz')
writeLAS(roi20, 'G:/Ontario_2023/Block_18/scratch/roi20.laz')
writeLAS(roi18, 'G:/Ontario_2023/Block_18/scratch/roi18.laz')

roi23 <- readLAS('G:/Ontario_2023/Block_18/scratch/roi23.laz')
roi20 <- readLAS('G:/Ontario_2023/Block_18/scratch/roi20.laz')
roi18 <- readLAS('G:/Ontario_2023/Block_18/scratch/roi18.laz')

# Plot cross section

plot_crossection <- function(las,
                             p1 = c(min(las@data$X), mean(las@data$Y)),
                             p2 = c(max(las@data$X), mean(las@data$Y)),
                             width = 4, colour_by = NULL)
{
  colour_by <- rlang::enquo(colour_by)
  data_clip <- clip_transect(las, p1, p2, width)
  p <- ggplot(data_clip@data, aes(X,Z)) + geom_point(size = 0.5) + coord_equal() + theme_minimal()

  if (!is.null(colour_by))
    p <- p + aes(color = !!colour_by) + labs(color = "")

  return(p)
}


# Revise plot cross seciton to plot all three lidar acquisitons simultaneously

plot_multi_crossection <- function(las1, las2, las3,
                                   p1 = c(min(las1@data$X), mean(las1@data$Y)),
                                   p2 = c(max(las1@data$X), mean(las1@data$Y)),
                                   width = 4, point_size = 0.5, colour_by = NULL, transect_length = NULL,
                                   push = 0)
{

  if(!is.null(transect_length)){
    p1 <- c(p1[1], p1[2] - transect_length/2)
    p2 <- c(p2[1], p2[2] - transect_length/2)
  }

  if(push != 0){
    p1 <- c(p1[1] - push, p1[2])
    p2 <- c(p2[1] - push, p2[2])
  }


  data_clip1 <- clip_transect(las1, p1, p2, width)
  data_clip2 <- clip_transect(las2, p1, p2, width)
  data_clip3 <- clip_transect(las3, p1, p2, width)

  df1 <- data_clip1@data %>% select(X,Y,Z) %>% mutate(lidar = 'DLS23')
  df2 <- data_clip2@data %>% select(X,Y,Z) %>% mutate(lidar = 'SPL20')
  df3 <- data_clip3@data %>% select(X,Y,Z) %>% mutate(lidar = 'SPL18')

  df <- rbind(df1, df2, df3)

  p <- ggplot(df, aes(X,Z, color = lidar)) + geom_point(size = point_size) + coord_equal() + theme_minimal() +
    theme(legend.position = 'bottom')

  return(p)
}

plot_multi_crossection(roi23, roi20, roi18, width = 0.5, point_size = 1, colour_by = NULL, transect_length = 2, push = 5)

