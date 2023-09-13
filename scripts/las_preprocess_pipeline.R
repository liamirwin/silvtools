# las pipeline

library(silvtools)

proj_dir <- 'D:/scratch/test_proj'

setup_als_dirs(proj_dir)

perform_tiling(proj_dir = proj_dir, tile_size = 500, num_cores = 4L)

perform_ground_classification(proj_dir = proj_dir, num_cores = 4L)

perform_height_normalization(proj_dir = proj_dir, num_cores = 4L)

x <- generate_dtm(proj_dir, dtm_res = 5, num_cores = 4L)
