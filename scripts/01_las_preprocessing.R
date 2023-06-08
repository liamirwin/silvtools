# Script to perform pre-processing to create products from a las file
# Can tile, classify ground, normalize, generate CHM, DSM, DTM, and metrics
# ----- Load Packages ----- #
library(lidR)
library(sf)
library(future)
library(glue)
library(lidRmetrics)
library(terra)
library(terra)
library(stringr)
library(geometry) # Required for rumple metrics
library(silvtools)
# ---- Processing switches ----

# Set Switches via ShinyApp

configure_las_process()

# ULS or DAP?
is_dap <- FALSE
# Run in parallel?
run_parallel <- T
num_cores <- 2L
# Tile area?
make_tile <- T
# Tile size (m)
tile_size <- 10
chunk_buf <- 10
# Classify ground points?
ground_classify <- T
# Normalize points?
normalize <- T
# Filter out outlier normalized returns?
filter_normalize <- F
# Create DSM?
make_dsm <- T
dsm_res <- 0.05
# Create CHM?
make_chm <- T
chm_res <- 0.05
# Create DTM?
make_dtm <- T
dtm_res <- 0.05
# Calculate Metrics?
make_mets <- T
met_res <- 1
# Is ALS?
is_als <- F
is_mls <- T
# List directories (each is one acquisiton of ULS/DAP)
blocks_dir <- list.dirs('H:/Quesnel_2022/process', recursive = FALSE)
# Omit these blocks from processing stream
processed <- c('CT1','CT2','CT3','CT4','CT5', 'CT1-T-DAP', 'CT1-DAP')
# target <- c('CT1')
blocks_dir <- blocks_dir[!basename(blocks_dir) %in% processed]
# blocks_dir <- blocks_dir[basename(blocks_dir) %in% target]
blocks_dir <- 'G:/Block_18/blocks/N'
blocks_dir <- 'I:/NZ_2023/Cass/ULS'
blocks_dir <- 'G:/Scantiques_Roadshow/Sites/BC/Vaseux Lake/Vaseux_2017'
blocks_dir <- list.dirs('I:/NZ_2023/Cass/MLS/plots', recursive = FALSE)
processed <- c('PLOT_1','PLOT_2','PLOT_3', 'PLOT_4')
blocks_dir <- blocks_dir[!basename(blocks_dir) %in% processed]

################################################################################
# START BUTTON
################################################################################

for(i in 1:length(blocks_dir)){

tictoc::tic()

if(length(blocks_dir) == 1){
    i <- 1
}

proj_dir <- blocks_dir[i]

# ----- Output directories -----

raster_output <- glue::glue('{proj_dir}/output/raster')
vector_output <- glue::glue('{proj_dir}/output/vector')


if(stringr::str_detect(basename(proj_dir), pattern = 'DAP') | is_dap){
  is_dap = TRUE
  # Set acquisition name (DAPYY_blockname)
  acq <- paste0('DAP22_', stringr::str_replace(basename(proj_dir), pattern = "-DAP", replacement = ""))
  print(paste0('Set acqusition type as DAP named ', acq))
} else if(stringr::str_detect(basename(proj_dir), pattern = 'MLS')){
  is_mls = TRUE
  # Set acquisition name (DAPYY_blockname)
  acq <- paste0('MLS22_', stringr::str_replace(basename(proj_dir), pattern = "-MLS", replacement = ""))
  print(paste0('Set acqusition type as MLS named ', acq))
} else{
  is_dap = FALSE
  # Set acquisition name (ULSYY_blockname)
  acq <- paste0('ULS23_',basename(proj_dir))
  print(paste0('Set acqusition type as lidar (ULS) named ', acq))
}
if(is_als){
  acq <- paste0('ALS_', stringr::str_replace(basename(proj_dir), pattern = "-DAP", replacement = ""))
  print(paste0('Set acqusition type as lidar (ALS) named ', acq))
}
if(is_mls){
  acq <- paste0('MLS_', stringr::str_replace(basename(proj_dir), pattern = "-DAP", replacement = ""))
  print(paste0('Set acqusition type as lidar (MLS) named ', acq))
}

# ---- Initialize Parallel Processing ----

if(run_parallel == TRUE){
  future::plan(future::multisession, workers = num_cores)
  print(glue::glue('Parallel processing initiated using {num_cores} cores'))
}

# ---- Tiling ----

if(make_tile == TRUE){
  tictoc::tic()
  # Create single las file catalog with raw file (careful with file size)
  ctg_tile <- catalog(glue::glue('{proj_dir}/input/las/raw'))
  tile_size <- tile_size
  print(glue::glue('Raw {acq} point cloud contains {density(ctg_tile) * sf::st_area(ctg_tile)} points at {round(density(ctg_tile))} pts/m2'))
  # Processing options
  opt_chunk_size(ctg_tile) <- tile_size
  opt_chunk_buffer(ctg_tile) <- 0
  opt_progress(ctg_tile) <- T
  opt_laz_compression(ctg_tile) <- T
  tile_dir <- glue::glue('{proj_dir}/input/las/tile')
  opt_output_files(ctg_tile) <- '{tile_dir}/{acq}_{XLEFT}_{YBOTTOM}'
  ctg_tile@output_options$drivers$LAS$param = list(index = TRUE)
  # Create directory for tiles if needed
  if(!dir.exists(tile_dir)){
    dir.create(tile_dir, recursive = T)
    print(glue::glue('Created a directory for tiled laz files at {tile_dir}'))
  }
  # Perform tiling of raw las file
  print(glue::glue('Beginning tiling process for {acq} at {tile_size}m'))
  ctg_tile <- catalog_retile(ctg_tile)
  print(glue::glue('Tiling process complete for {acq} {nrow(ctg_tile)} {tile_size}m tiles created'))
  # Index Tiles
  lidR:::catalog_laxindex(ctg_tile)
  print(glue::glue('Indexing of tiles for {acq} {nrow(ctg_tile)} {tile_size}m complete'))
  tictoc::toc()
}

# ---- Classify ground points ----

if (ground_classify == TRUE) {
  tile_dir <- glue::glue('{proj_dir}/input/las/tile')
  ctg_tile <- catalog(tile_dir)
  # Processing options
  opt_progress(ctg_tile) <- T
  opt_laz_compression(ctg_tile) <- T
  opt_chunk_size(ctg_tile) <- tile_size
  opt_chunk_buffer <- chunk_buf
  class_dir <- glue::glue('{proj_dir}/input/las/class')
  opt_output_files(ctg_tile) <- '{class_dir}/{acq}_{XLEFT}_{YBOTTOM}_class'
  print(glue::glue('Beginning ground classification process for {acq}'))
  # Create directory for classified point clouds if needed
  if(!dir.exists(class_dir)){
    dir.create(class_dir, recursive = T)
    print(glue::glue('Created a directory for classified laz files at {class_dir}'))
  }
  if(is_dap){
    las <- readLAS(ctg_tile[1,])
    if('confidence' %in% names(las)){
    opt_filter(ctg_tile) = "-keep_random_fraction 0.01 -keep_attribute_above 0 2"
    print(glue::glue('Confidence attribute detected in DAP point cloud; filtering high confidence points for ground classification'))
    } else{
    opt_filter(ctg_tile) = "-keep_random_fraction 0.01"
    print(glue::glue('No confidence attribute detected in DAP point cloud'))
    }
    ctg_class <- classify_ground(ctg_tile, algorithm = csf(
      sloop_smooth = FALSE,
      class_threshold = 0.07,
      # larger resolution = coarser DTM
      cloth_resolution = .7,
      rigidness = 2L,
      # even for non flat sites, this should be 3L in my experience
      iterations = 500L,
      time_step = 0.65))
  } else if(is_mls == TRUE){
    opt_filter(ctg_tile) = "-keep_random_fraction 0.01"
    ctg_class <- classify_ground(ctg_tile, algorithm = csf(
      sloop_smooth = FALSE,
      class_threshold = 0.07,
      # larger resolution = coarser DTM
      cloth_resolution = .7,
      rigidness = 2L,
      # even for non flat sites, this should be 3L in my experience
      iterations = 500L,
      time_step = 0.65))
    }else{
    ctg_class <-

    }
  lidR:::catalog_laxindex(ctg_class)
  print(glue::glue('Ground classification process for {acq} complete'))
}

# ----- DTM -----
if(make_dtm == TRUE){
  ctg_class <- catalog(glue::glue('{proj_dir}/input/las/class'))
  opt_progress(ctg_class) <- T
  opt_chunk_buffer(ctg_class) <- chunk_buf
  opt_stop_early(ctg_class) <- FALSE
  if (!dir.exists(glue::glue('{raster_output}/dtm'))) {
    dir.create(glue::glue('{raster_output}/dtm'), showWarnings = FALSE, recursive = T)
    print(glue::glue('Created a directory for DTM {raster_output}/dtm'))
  }
  dtm_res <- dtm_res
  print(glue::glue('Generating terrain model for {acq} at {dtm_res}m'))
  if(!dir.exists(glue::glue('{raster_output}/dtm/tiles'))){
    dir.create(glue::glue('{raster_output}/dtm/tiles'), showWarnings = FALSE)
    print(glue::glue('Created a directory for DTM tiles {raster_output}/dtm/tiles'))
  }
  opt_output_files(ctg_class) <-
    '{raster_output}/dtm/tiles/{acq}_dtm_{dtm_res}_{XLEFT}_{YBOTTOM}'
  ctg_class@output_options$drivers$SpatRaster$param$overwrite <-
    TRUE

  if(is_mls){
  dtm_algorithm <- knnidw()
  } else{
   dtm_algorithm <- tin()
  }

  rasterize_terrain(ctg_class, res = dtm_res, algorithm = dtm_algorithm)

  #--- Load DTM Tiles as virtual raster dataset ---
  dtm_tiles <-
    list.files(glue::glue('{raster_output}/dtm/tiles/'),
               pattern = '.tif$',
               full.names = T)
  dtm <-
    terra::vrt(dtm_tiles,
        glue::glue('{raster_output}/dtm/tiles/{acq}_dtm.vrt'),
        overwrite = T)

  terra::writeRaster(dtm, filename = glue::glue('{raster_output}/dtm/{acq}_dtm_tin_{dtm_res}m.tif'), overwrite = T)

  dtm_smooth <- dtm %>%
    focal(w = matrix(1, 25, 25),
          fun = mean,
          na.rm = TRUE,
          pad = TRUE)

  terra::writeRaster(dtm_smooth, filename = glue::glue('{raster_output}/dtm/{acq}_dtm_tin_smooth_{dtm_res}m.tif'), overwrite = T)
  print(glue::glue('DTM generation process complete for {acq}'))
}


# ---- Normalize point cloud ----

if(normalize == TRUE){
  class_dir <- glue::glue('{proj_dir}/input/las/class')
  ctg_class <- catalog(class_dir)
  # Processing options
  opt_progress(ctg_class) <- T
  opt_laz_compression(ctg_class) <- T
  opt_chunk_buffer(ctg_class) <- chunk_buf
  opt_filter(ctg_class) <- '-drop_as_witheld'
  norm_dir <- glue::glue('{proj_dir}/input/las/norm')
  opt_output_files(ctg_class) <- '{norm_dir}/{acq}_{XLEFT}_{YBOTTOM}_norm'
  print(glue::glue('Beginnng height normalization process for {acq}'))
  # Create directory for normalized point clouds if needed
  if(!dir.exists(norm_dir)){
    dir.create(norm_dir, recursive = T)
    print(glue::glue('Created a directory for normalized laz files at {norm_dir}'))
  }
  if(is_dap | is_mls){
    tile_dir <- glue::glue('{proj_dir}/input/las/tile')
    ctg_tile <- catalog(tile_dir)
    opt_progress(ctg_tile) <- T
    opt_laz_compression(ctg_tile) <- T
    opt_chunk_buffer(ctg_tile) <- chunk_buf
    norm_dir <- glue::glue('{proj_dir}/input/las/norm')
    opt_output_files(ctg_tile) <- '{norm_dir}/{acq}_{XLEFT}_{YBOTTOM}_norm'
    dtm <- rast(glue::glue('{raster_output}/dtm/{acq}_dtm_tin_smooth_{dtm_res}m.tif'))
    ctg_norm <- lidR::normalize_height(ctg_tile, dtm, na.rm = TRUE)
  } else{
  # Normalize point cloud
  ctg_norm <- lidR::normalize_height(ctg_class, tin())
  }
  if(filter_normalize == F){
  lidR:::catalog_laxindex(ctg_norm)
  print('Indexed normalized tiles...')
  }
}


if(filter_normalize == T){
  print(glue::glue('Filtering potential outliers from normalized tiles'))
  norm_dir <- glue::glue('{proj_dir}/input/las/norm')
  ctg_norm <- lidR::catalog(norm_dir)
  opt_output_files(ctg_norm) <- '{norm_dir}/{acq}_{XLEFT}_{YBOTTOM}_norm'
  opt_filter(ctg_norm)       <- "-drop_z_below 0"
  opt_chunk_buffer(ctg_norm) <- chunk_buf
  opt_laz_compression(ctg_norm) <- T
  opt_progress(ctg_norm) <- T
  # Removes outliers in by calculating the 95th percentile
  # of height in 10x10m pixels, and filtering out points above the 95th percentile plus 20%
  filter_noise = function(las, sensitivity)
  {
    if (is(las, "LAS"))
    {
      p95 <- pixel_metrics(las, ~quantile(Z, probs = 0.95), 10)
      las <- merge_spatial(las, p95, "p95")
      las <- filter_poi(las, Z < p95*sensitivity)
      las$p95 <- NULL
      return(las)
    }

    if (is(las, "LAScatalog"))
    {
      res <- catalog_map(las, filter_noise, sensitivity = sensitivity)
      return(res)
    }
  }
  filter_noise(ctg_norm, sensitivity = 1.2)
  ctg_norm <- lidR::catalog(norm_dir)
  lidR:::catalog_laxindex(ctg_norm)
  print('Indexed normalized tiles...')
}

# ---- Metrics -----

if(make_mets == TRUE){
  if(!dir.exists(glue::glue('{raster_output}/metrics'))){
    dir.create(glue::glue('{raster_output}/metrics'), recursive = T)
    print(glue::glue('Created a directory for metrics at {raster_output}/metrics'))
  }
  met_res <- met_res
  ctg_norm <- catalog(glue::glue('{proj_dir}/input/las/norm'))
  opt_progress(ctg_norm) <- T
  opt_chunk_buffer(ctg_norm) <- chunk_buf
  opt_select(ctg_norm) <- "xyz"
  opt_filter(ctg_norm) <- "-drop_withheld -drop_z_below 0"
  opt_progress(ctg_norm) <- T


  # Basic metrics suite
  # -------------------

  print(glue::glue('Generating basic metrics for {acq} at {met_res}m'))
  basic <- pixel_metrics(ctg_norm, ~lidRmetrics::metrics_basic(Z), res = met_res)
  terra::writeRaster(x = basic, filename = glue::glue('{raster_output}/metrics/{acq}_basic_{met_res}m.tif'), overwrite = TRUE)


}

# ----- Canopy Height Model -----
if (make_chm == TRUE) {
  if (!dir.exists(glue::glue('{raster_output}/chm'))) {
    dir.create(glue::glue('{raster_output}/chm'), showWarnings = FALSE)
    print(glue::glue('Created a directory for CHM {raster_output}/chm'))
  }
  # Processing options
  ctg_norm <- catalog(glue::glue('{proj_dir}/input/las/norm'))
  opt_progress(ctg_norm) <- T
  opt_chunk_buffer(ctg_norm) <- chunk_buf
  chm_res <- chm_res
  print(glue::glue('Generating canopy height model for {acq} at {chm_res}m'))
  # Write as tiles for  for later mosaic
  dir.create(glue::glue('{raster_output}/chm'), showWarnings = FALSE)
  dir.create(glue::glue('{raster_output}/chm/tiles'), showWarnings = FALSE)
  opt_output_files(ctg_norm) <-
    '{raster_output}/chm/tiles/{acq}_chm_{chm_res}_{XLEFT}_{YBOTTOM}'
  ctg_norm@output_options$drivers$SpatRaster$param$overwrite <-
    TRUE
  rasterize_canopy(ctg_norm,
                   res = chm_res,
                   algorithm = p2r(na.fill = knnidw()))
  #--- Load CHM Tiles as virtual raster dataset ---
  chm_tiles <-
    list.files(glue::glue('{raster_output}/chm/tiles/'),
               pattern = '.tif',
               full.names = T)
  chm <-
    vrt(chm_tiles,
        glue::glue('{raster_output}/chm/tiles/{acq}_chm.vrt'),
        overwrite = T)
  #--- Remove extreme CHM values. Extreme points get 0 or 30 ---
  #chm <- terra::clamp(chm, 0, 30, values = TRUE)
  #--- Set layer name to Z ---
  names(chm) <- 'Z'
  # ----- Fill CHM -----
  print(glue::glue('Filling canopy height model for {acq} at {chm_res}m'))
  chm_filled <- terra::focal(
    chm,
    w = 3,
    fun = "mean",
    na.policy = "only",
    na.rm = TRUE
  )
  names(chm_filled) <- 'Z'
  # ----- Smooth CHM -----
  fgauss <- function(sigma, n = ws) {
    m <- matrix(ncol = n, nrow = n)
    col <- rep(1:n, n)
    row <- rep(1:n, each = n)
    x <- col - ceiling(n/2)
    y <- row - ceiling(n/2)
    m[cbind(row, col)] <- 1/(2 * pi * sigma^2) * exp(-(x^2 +
                                                         y^2)/(2 * sigma^2))
    m/sum(m)
  }
  print(glue::glue('Smoothing canopy height model for {acq} at {chm_res}m'))
  chm_smooth <- terra::focal(chm_filled,
                             w = fgauss(1, n = 5))
  names(chm_smooth) <- 'Z'
  # ---- Write CHM ----
  print(glue::glue('Writing canopy height models for {acq} at {chm_res}m'))
  terra::writeRaster(chm,
                     glue::glue('{raster_output}/chm/{acq}_chm_{chm_res}m_sub0_p2r.tif'),
                     overwrite = T)
  terra::writeRaster(
    chm_filled,
    filename = glue::glue('{raster_output}/chm/{acq}_chm_fill_p2r_{chm_res}m.tif'),
    overwrite = T
  )
  terra::writeRaster(
    chm_smooth,
    filename = glue::glue(
      '{raster_output}/chm/{acq}_chm_smooth_p2r_{chm_res}m.tif'
    ),
    overwrite = T
  )
  # ---- Delete intermediate raster tiles ----
  if(dir.exists(glue::glue('{raster_output}/chm/tiles'))){
    unlink(glue::glue('{raster_output}/chm/tiles'), recursive = TRUE)
    print('Deleted intermediate CHM tiles')
  }
}

# ----- Digital Surface Model -----

if (make_dsm == TRUE) {
  if (!dir.exists(glue::glue('{raster_output}/dsm'))) {
    dir.create(glue::glue('{raster_output}/dsm'), showWarnings = FALSE)
    print(glue::glue('Created a directory for DSM {raster_output}/dsm'))
  }
  ctg_class <- catalog(glue::glue('{proj_dir}/input/las/tile'))
  opt_progress(ctg_class) <- T
  opt_chunk_buffer(ctg_class) <- chunk_buf
  dsm_res <- dsm_res
  print(glue::glue('Generating digital surface model for {acq} at {dsm_res}m'))

  #--- Process keeps failing as we cant load the raster into memory ---
  #--- Write as tiles for each tile instead for later mosaic ---
  dir.create(glue::glue('{raster_output}/dsm'), showWarnings = FALSE)
  dir.create(glue::glue('{raster_output}/dsm/tiles'), showWarnings = FALSE)
  opt_output_files(ctg_class) <-
    '{raster_output}/dsm/tiles/{acq}_dsm_{dsm_res}_{XLEFT}_{YBOTTOM}'
  ctg_class@output_options$drivers$SpatRaster$param$overwrite <-
    TRUE
  opt_progress(ctg_class) <- T
  opt_chunk_buffer(ctg_class) <- chunk_buf
  rasterize_canopy(ctg_class,
                   res = dsm_res,
                   algorithm = p2r(na.fill = knnidw()))
  #--- Load DSM Tiles as virtual raster dataset ---
  dsm_tiles <-
    list.files(glue::glue('{raster_output}/dsm/tiles/'),
               pattern = '.tif',
               full.names = T)
  dsm <-
    vrt(dsm_tiles,
        glue::glue('{raster_output}/dsm/tiles/{acq}_dsm.vrt'),
        overwrite = T)
  #--- Set layer name to Z ---
  names(dsm) <- 'Z'
  #--- Write merged raster to disk ---
  terra::writeRaster(dsm,
                     glue::glue('{raster_output}/dsm/{acq}_dsm_{dsm_res}m_sub0_p2r.tif'),
                     overwrite = T)
  # ----- Fill dsm -----
  print(glue::glue('Filling digital surface model for {acq} at {dsm_res}m'))
  dsm_filled <- terra::focal(
    dsm,
    w = 3,
    fun = "mean",
    na.policy = "only",
    na.rm = TRUE
  )
  names(dsm_filled) <- 'Z'
  # ----- Smooth dsm -----
  print(glue::glue('Smoothing dsm for {acq} at {dsm_res}m'))
  dsm_smooth <- terra::focal(dsm_filled,
                             w = fgauss(1, n = 5))
  names(dsm_smooth) <- 'Z'
  # ---- Write DSM ----
  print(glue::glue('Writing dsm for {acq} at {dsm_res}m'))
  terra::writeRaster(
    dsm_filled,
    filename = glue::glue('{raster_output}/dsm/{acq}_dsm_fill_p2r_{dsm_res}m.tif'),
    overwrite = T
  )
  terra::writeRaster(
    dsm_smooth,
    filename = glue::glue(
      '{raster_output}/dsm/{acq}_dsm__smooth_p2r_{dsm_res}m.tif'
    ),
    overwrite = T
  )
  # ---- Delete intermediate raster tiles ----
  if(dir.exists(glue::glue('{raster_output}/dsm/tiles'))){
    unlink(glue::glue('{raster_output}/dsm/tiles'), recursive = TRUE)
    print('Deleted intermediate DSM tiles')
  }
}

print(glue::glue('Finished point cloud processing for {acq}'))
tictoc::toc()

}
