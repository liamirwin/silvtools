
library(lasR)
proj_dir <- 'G:/MKRF_2023/DLS_2024/PILE'
proj_dir <- 'G:/Ontario_2023/Block_18/LAS/DLS23'
silvtools::setup_als_dirs(proj_dir)

f <- list.files(glue::glue('{proj_dir}/input/las/class'), pattern = '.laz', full.names = T)
norm_dir <- glue::glue('{proj_dir}/input/las/norm')
dtm_dir <- glue::glue('{proj_dir}/output/raster/dtm')
chm_dir <- glue::glue('{proj_dir}/output/raster/chm')
ttops_dir <- glue::glue('{proj_dir}/output/vector/treetops')



# Create all dirs
dir.create(norm_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(dtm_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(chm_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(ttops_dir, showWarnings = FALSE, recursive = TRUE)

dtm_res <- 0.25
chm_res <- 0.10
metric_res <- 1


tri <- triangulate(filter = keep_ground())
dtm <- rasterize(res = dtm_res, tri, ofile = paste0(dtm_dir, "/*_dtm_", dtm_res, "m.tif"))
norm <- transform_with(tri, "-")
chm <- rasterize(res = chm_res, "max", ofile = paste0(chm_dir, "/*_chm_", chm_res, "m.tif"))

ws <- 3
min_height <- 10

lmf <- local_maximum_raster(ws = ws, min_height = min_height, ofile = paste0(ttops_dir, "/*_lmf_ws", ws, "_hmin", min_height, ".gpkg"), raster = chm)

# Reads tiled, classified LAZ files, triangulates, generates a DTM, normalizes heights, and generates a CHM
pipeline <- reader_las() + tri + dtm + norm + write_las(paste0(norm_dir, "/*.laz")) + write_lax() + chm +
  lmf




# Pipeline without writing normalized LAZ files
pipeline <- reader_las() + tri + dtm + norm + chm



# Execute the pipeline
tictoc::tic()
exec(pipeline, on = f, ncores = nested(4, 4), with = list(progress = TRUE))
tictoc::toc()
