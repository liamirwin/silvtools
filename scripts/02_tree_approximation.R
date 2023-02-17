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
library(lidRplugins)

# List directories (each is one acquisiton of ULS/DAP)
blocks_dir <- list.dirs('E:/Quesnel_2022/process', recursive = FALSE)
# Omit these already processed blocks from processing stream
processed <- c('CT1','CT1-T-DAP','CT3','CT4','CT5')
blocks_dir <- blocks_dir[!basename(blocks_dir) %in% processed]
# ULS or DAP?
is_dap <- FALSE

# ---- Treetop Loop (CHM) ----

for(i in 1:length(blocks_dir)){

  tictoc::tic()

  # ---- Project setup ----

  proj_dir <- blocks_dir[i]

  if(stringr::str_detect(basename(proj_dir), pattern = 'DAP')){
    is_dap == TRUE
    # Set acquisition name (DAPYY_blockname)
    acq <- glue::glue('DAP22_{stringr::str_replace(basename(proj_dir), pattern = "-DAP", replacement = "")}')
    print(glue::glue('Set acqusition type as DAP named {acq}'))
  } else{
    is_dap == FALSE
    # Set acquisition name (ULSYY_blockname)
    acq <- glue::glue('ULS22_{basename(proj_dir)}')
    print(glue::glue('Set acqusition type as lidar (ULS) named {acq}'))
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


for(i in 1:length(blocks_dir)){

tictoc::tic()

# ---- Project setup ----

proj_dir <- blocks_dir[i]

if(stringr::str_detect(basename(proj_dir), pattern = 'DAP')){
  is_dap == TRUE
  # Set acquisition name (DAPYY_blockname)
  acq <- glue::glue('DAP22_{stringr::str_replace(basename(proj_dir), pattern = "-DAP", replacement = "")}')
  print(glue::glue('Set acqusition type as DAP named {acq}'))
} else{
  is_dap == FALSE
  # Set acquisition name (ULSYY_blockname)
  acq <- glue::glue('ULS22_{basename(proj_dir)}')
  print(glue::glue('Set acqusition type as lidar (ULS) named {acq}'))
}

# ----- Output directories -----

raster_output <- glue::glue('{proj_dir}/output/raster')
vector_output <- glue::glue('{proj_dir}/output/vector')

# ---- Raster Watershed Segmentation ----

chm <- terra::rast(stringr::str_subset(list.files(paste0(raster_output, '/chm'), pattern = '.tif$',full.names = T), pattern = 'smooth'))
ttops_files <- list.files(glue::glue('{vector_output}/treetops'), pattern = '.gpkg', full.names = T)

# ---- Load Tree tops ----

# Three different lmf methods lmf(ws = 2) lmfauto( ), variable ws lmf
ttops_ws2 <- st_read(str_subset(ttops_files, pattern = 'ws2'), quiet = TRUE)
ttops_auto <- st_read(str_subset(ttops_files, pattern = 'auto'), quiet = TRUE)
ttops_v <- st_read(str_subset(ttops_files, pattern = '_v_'), quiet = TRUE)

print(glue::glue('Successfully loaded all sets of treetops for {acq}
                 lmfws2 = {nrow(ttops_ws2)} trees
                 lmfauto = {nrow(ttops_auto)} trees
                 lmfv = {nrow(ttops_v)}'))

# ---- Segment Crowns ----

print(glue::glue('Begnning Crown Segmentation for {acq}'))

# lmf(ws = 2)
crowns <- silvtools::crown_mask(chunk = chm, ttops = ttops_ws2, crown_height_threshold = 0.25, vis = FALSE)
print('lmf(ws = 2) crown segmentation complete...')
# lmfauto( )
crowns_auto <- silvtools::crown_mask(chunk = chm, ttops = ttops_auto, crown_height_threshold = 0.25, vis = FALSE)
print('lmfauto( ) crown segmentation complete...')
# lmfv()
crowns_v <- silvtools::crown_mask(chunk = chm, ttops = ttops_v, crown_height_threshold = 0.25, vis = FALSE)
print('lmfv( ) crown segmentation complete...')

# ---- Clean Crowns ----

# Clean crowns, take largest polygon for each treeID; fill holes within crowns

print('Cleaning lmf(ws = 2) crowns...')
# lmf(ws = 2)
crowns_p <- sf::st_as_sf(terra::as.polygons(crowns)) %>%
  convert_multi_to_single_polygons(polygons = ., fill_holes = TRUE)

print('Cleaning lmfauto( ) crowns...')
# lmfauto( )
crowns_auto_p <- sf::st_as_sf(terra::as.polygons(crowns_auto)) %>%
  convert_multi_to_single_polygons(polygons = ., fill_holes = TRUE)

print('Cleaning lmfv( ) crowns...')
# lmfv()
crowns_v_p <- sf::st_as_sf(terra::as.polygons(crowns_v)) %>%
  convert_multi_to_single_polygons(polygons = ., fill_holes = TRUE)s

dir.create(glue::glue('{vector_output}/crowns'), showWarnings = FALSE, recursive = T)
dir.create(glue::glue('{raster_output}/crowns'), showWarnings = FALSE, recursive = T)

# ---- Write Crowns ----

print(glue::glue('Segmentation process complete; attempting to write raster and vector crowns for {acq}'))

# Write lmf(ws = 2) crowns
terra::writeRaster(crowns, glue::glue('{raster_output}/crowns/{acq}_lmf_ws2_watershed_crowns.tif'), overwrite = T)
sf::st_write(crowns_p, glue::glue('{vector_output}/crowns/{acq}_lmf_ws2_watershed_crowns.shp'), append = FALSE)
# Write lmfauto( ) crowns
terra::writeRaster(crowns_auto, glue::glue('{raster_output}/crowns/{acq}_lmf_auto_watershed_crowns.tif'), overwrite = T)
sf::st_write(crowns_auto_p, glue::glue('{vector_output}/crowns/{acq}_lmf_auto_watershed_crowns.shp'), append = FALSE)
# Write lmfv( ) crowns
terra::writeRaster(crowns_v, glue::glue('{raster_output}/crowns/{acq}_lmf_v_watershed_crowns.tif'), overwrite = T)
sf::st_write(crowns_v_p, glue::glue('{vector_output}/crowns/{acq}_lmf_v_watershed_crowns.shp'), append = FALSE)

print(glue::glue('Wrote crown rasters for {acq}'))

tictoc::toc()

}















# ---- Treetop Loop LAS ----

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

  ctg_class <- catalog(glue::glue('{proj_dir}/input/class'))
  opt_progress(ctg_class) <- T
  ctg_norm <- catalog(glue::glue('{proj_dir}/input/norm'))
  opt_progress(ctg_norm) <- T

  print(glue::glue('Catalogs setup for {acq}'))

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


  dir.create(glue::glue('{vector_output}/treetops'), showWarnings = FALSE, recursive = T)

  st_write(lmf_ws2_df, glue::glue('{vector_output}/treetops/{acq}_lmf_ws2_ttops.shp'))
  # st_write(lmf_auto, glue('{vector_output}/treetops/{acq}_lmf_auto_ttops.shp'))

  print(glue::glue('Found {nrow(lmf_ws2_df)} treetops in {acq}'))
  tictoc::toc()

}
