# las pipeline
load_all('F:/R/silvtools')
proj_dir <- 'G:/MKRF_2023/DLS_2024/A_POST'
proj_dir <- 'G:/MKRF_2023/DLS_2024/C_POST'
proj_dir <- 'G:/MKRF_2023/DLS_2024/PILE'
proj_dir <- 'G:/MKRF_2023/ALS_2022'
proj_dir <- 'G:/Misc/still_creek_als'

proj_dir <- 'G:/Ontario_2023/Block_18/LAS/DLS23'
proj_dir <- 'G:/Ontario_2023/Block_18/LAS/SPL18'
proj_dir <- 'G:/Ontario_2023/Block_18/LAS/SPL20'

setup_als_dirs(proj_dir)
num_cores <- 6L
silvtools::perform_tiling(proj_dir = proj_dir, tile_size = 250, num_cores = num_cores, tile_dir = glue::glue("{proj_dir}/input/las/class"))

# perform_ground_classification(proj_dir = proj_dir, num_cores = 4L)

perform_height_normalization(proj_dir = proj_dir, num_cores = num_cores)

generate_dtm(proj_dir, res = 0.25, num_cores = num_cores, dtm_algorithm = 'tin')

generate_chm(proj_dir, res = 0.25, num_cores = num_cores, algorithm = 'p2r', clamp = T, min_clamp = 0, max_clamp = 65)

approximate_chm_treetops(proj_dir, fixed_window = TRUE, auto_window = FALSE, variable_window = FALSE ,
                         hmin = 10, fix_ws = 3, chm_ext = "smooth_p2r")

chm <- rast('G:/MKRF_2023/ALS_2022/output/raster/chm/ALS_2022_chm_fill_p2r_sc0_0.25m.tif')
chm <- chm * 1
ttops <- locate_trees(chm, algorithm = lmf(ws = 3, hmin = 10))


generate_chm(proj_dir, res = 0.05, num_cores = num_cores, algorithm = 'p2r', clamp = F, dsm = T)
