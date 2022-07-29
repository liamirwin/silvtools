dir_setup_als <- function(proj_dir){

  # Generate standardized list of inputs and outputs for assoicated spatial data
  # Multilevel, /input, /output, /doc, /script
  # Doc intended to store relevant literature, site info, data collection info
  # Script intended to house R code and other code used to process input dataset

  targ_dirs <- file.path(proj_dir, c('input','output','doc',
                                     'script','input/lidar','input/lidar/norm','input/lidar/class','input/vector',
                                     'input/raster', 'output/lidar',
                                     'output/vector','output/raster',
                                     'output/vector/treetops', 'output/vector/crowns',
                                     'output/raster/crowns'))

  # Apply dir.create to generate standard folders

  lapply(targ_dirs, FUN = dir.create , showWarnings = FALSE)

}
