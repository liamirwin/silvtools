# ----- Load Packages ----- #

library(lidR)
library(sf)
library(future)
library(glue)
library(lidRmetrics)
library(terra)
library(mapview)
library(terra)
library(stringr)
library(geometry) # Required for rumple metrics


# List directories (each is one acquisiton of ULS/DAP)
blocks_dir <- list.dirs('H:/Quesnel_2022/process', recursive = FALSE)
# Omit these already processed blocks from processing stream
processed <- c('CT1','CT2','CT3','CT4','CT5')
blocks_dir <- blocks_dir[!basename(blocks_dir) %in% processed]
# ULS or DAP?
is_dap <- FALSE

# ---- Treetop Loop (CHM) ----

for(i in 1:length(blocks_dir)){

  tictoc::tic()

  proj_dir <- blocks_dir[i]

  if(stringr::str_detect(basename(proj_dir), pattern = 'DAP')){
    is_dap == TRUE
    # Set acquisition name (DAPYY_blockname)
    acq <- paste0('DAP22_', stringr::str_replace(basename(proj_dir), pattern = "-DAP", replacement = ""))
    print(paste0('Set acqusition type as DAP named ', acq))
  } else{
    is_dap == FALSE
    # Set acquisition name (ULSYY_blockname)
    acq <- paste0('ULS22_',basename(proj_dir))
    print(paste0('Set acqusition type as lidar (ULS) named ', acq))
  }

  # ----- Output directories -----

  raster_output <- glue::glue('{proj_dir}/output/raster')
  vector_output <- glue::glue('{proj_dir}/output/vector')

  # ---- Find Treetops ----

  chm_files <- list.files(glue::glue('{proj_dir}/output/raster/chm'), pattern = '.tif$', full.names = T)
  chm <- terra::rast(str_subset(chm_files, pattern = 'smooth'))
  # Fixed window size (2m)
  lmf_ws2 <- locate_trees(chm, lmf(ws = 2, hmin = 5), uniqueness = 'incremental')
  # Lmf auto
  lmf_auto <- locate_trees(chm, lidRplugins::lmfauto() , uniqueness = 'incremental')
  # Variable window size
  f <- function(x) {x * 0.07 + 1}
  lmf_v <- locate_trees(chm, lmf(f, hmin = 5), uniqueness = 'incremental')

  if(!dir.exists(glue::glue('{vector_output}/treetops'))){
    dir.create(glue::glue('{vector_output}/treetops'), recursive = T)
    print(glue::glue('Created tree top output directory for {acq}'))
  }

  sf::st_write(lmf_ws2, glue::glue('{vector_output}/treetops/{acq}_lmf_ws2_ttops.gpkg'))
  sf::st_write(lmf_auto, glue::glue('{vector_output}/treetops/{acq}_lmf_auto_ttops.gpkg'))
  sf::st_write(lmf_v, glue::glue('{vector_output}/treetops/{acq}_lmf_v_ttops.gpkg'))

  print(glue::glue('Found {nrow(lmf_ws2)} ws2 treetops, {nrow(lmf_auto)} auto treetops, and {nrow(lmf_v)} variable ws treetops in {acq}'))
  tictoc::toc()
}

# ---- Crown Segmentation ----


for(i in 1:length(blocks_dir))

tictoc::tic()

proj_dir <- blocks_dir[i]

if(stringr::str_detect(basename(proj_dir), pattern = 'DAP')){
  is_dap == TRUE
  # Set acquisition name (DAPYY_blockname)
  acq <- paste0('DAP22_', stringr::str_replace(basename(proj_dir), pattern = "-DAP", replacement = ""))
  print(paste0('Set acqusition type as DAP named ', acq))
} else{
  is_dap == FALSE
  # Set acquisition name (ULSYY_blockname)
  acq <- paste0('ULS22_',basename(proj_dir))
  print(paste0('Set acqusition type as lidar (ULS) named ', acq))
}

# ----- Output directories -----

raster_output <- glue::glue('{proj_dir}/output/raster')
vector_output <- glue::glue('{proj_dir}/output/vector')

# ---- Raster Watershed Segmentation ----

chm <- terra::rast(stringr::str_subset(list.files(paste0(raster_output, '/chm'), pattern = '.tif$',full.names = T), pattern = 'smooth'))
ttops_files <- list.files(glue::glue('{vector_output}/treetops'), pattern = '.gpkg', full.names = T)
ttops <- st_read(str_subset(ttops_files, pattern = 'ws2'))
ttops <- st_read(glue('{vector_output}/treetops/{acq}_lmf_ws2_ttops.gpkg$'))

crowns <- silvtools::crown_mask(chunk = chm, ttops = ttops, chm_cutoff_percent = 0.25, vis = FALSE)
crowns_auto <- silvtools::crown_mask(chunk = chm, ttops = ttops, chm_cutoff_percent = 0.25, vis = FALSE)
crowns_p <- sf::st_as_sf(terra::as.polygons(crowns)) %>% convert_multi_to_single_polygons(polygons = ., fill_holes = TRUE)

# crowns_auto <- silvtools::mcwatershed(chm = chm, treetops = ttops_auto, th_tree = 2, tree_id = 'treeID')()

dir.create(glue('{vector_output}/crowns'), showWarnings = FALSE, recursive = T)
dir.create(glue('{raster_output}/crowns'), showWarnings = FALSE, recursive = T)

terra::writeRaster(crowns, glue('{raster_output}/crowns/{acq}_lmf_ws2_watershed_crowns.tif'), overwrite = T)
sf::st_write(crowns_p, glue('{vector_output}/crowns/{acq}_lmf_ws2_watershed_crowns.shp'), append = FALSE)
# writeRaster(crowns_auto, glue('{vector_output}/crowns/{acq}_lmf_auto_watershed_crowns.tif'))

print(paste0('Wrote crown rasters for ', acq))

}

# ---- Treetop Loop ----

for(i in 1:length(blocks_dir)){

  tictoc::tic()

  proj_dir <- blocks_dir[i]

  if(stringr::str_detect(basename(proj_dir), pattern = 'DAP')){
    is_dap == TRUE
    # Set acquisition name (DAPYY_blockname)
    acq <- paste0('DAP22_', stringr::str_replace(basename(proj_dir), pattern = "-DAP", replacement = ""))
    print(paste0('Set acqusition type as DAP named ', acq))
  } else{
    is_dap == FALSE
    # Set acquisition name (ULSYY_blockname)
    acq <- paste0('ULS22_',basename(proj_dir))
    print(paste0('Set acqusition type as lidar (ULS) named ', acq))
  }

  ctg_class <- catalog(glue('{proj_dir}/input/class'))
  opt_progress(ctg_class) <- T
  ctg_norm <- catalog(glue('{proj_dir}/input/norm'))
  opt_progress(ctg_norm) <- T

  print(glue('Catalogs setup for {acq}'))

  # ---- Initialize Parallel Processing ----

  plan(multisession, workers = 4L)

  # ----- Output directories -----

  raster_output <- glue::glue('{proj_dir}/output/raster')
  vector_output <- glue::glue('{proj_dir}/output/vector')

  # ---- Find Treetops ----

  opt_chunk_buffer(ctg_norm) <- 5
  lmf_ws2 <- locate_trees(ctg_norm, lmf(ws = 2, hmin = 2), uniqueness = 'incremental')
  # lmf_auto <- locate_trees(ctg_norm, lidRplugins::lmfauto() , uniqueness = 'incremental')
  lmf_ws2[sapply(lmf_ws2, function(x) dim(x)[1]) > 0]
  lmf_ws2_df <- do.call(rbind, lmf_ws2)
  lmf_ws2_df <- lmf_ws2_df %>% dplyr::mutate(treeID = 1:nrow(lmf_ws2_df))


  dir.create(glue('{vector_output}/treetops'), showWarnings = FALSE, recursive = T)

  st_write(lmf_ws2_df, glue::glue('{vector_output}/treetops/{acq}_lmf_ws2_ttops.shp'))
  # st_write(lmf_auto, glue('{vector_output}/treetops/{acq}_lmf_auto_ttops.shp'))

  print(glue('Found {nrow(lmf_ws2_df)} treetops in {acq}'))
  tictoc::toc()

}
