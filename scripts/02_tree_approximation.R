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
library(furrr)
library(exactextractr)



# List directories (each is one acquisiton of ULS/DAP)
blocks_dir <- list.dirs('H:/Quesnel_2022/process', recursive = FALSE)
# Omit these already processed blocks from processing stream
processed <- c("CT1",'CT2','CT3','CT4','CT1-T-DAP', 'CT2-T-DAP','CT2-DAP','CT1-DAP')
blocks_dir <- blocks_dir[!basename(blocks_dir) %in% processed]
# blocks_dir <- 'F:/Quesnel_2022/GeoSLAM/plot_las/CT1P1'
# blocks_dir <- 'Y:/Irwin/NZ_2023/Campbell/Campbell_ULS'
blocks_dir <- 'G:/Block_18/blocks/N'


# Processing Switches

# ULS or DAP?
is_dap <- TRUE
# Run in parallel?
run_parallel <- T
num_cores <- 3L
make_tiles <- TRUE
tile_size <- 100
# ---- Treetop Loop (CHM) ----

for(i in 1:length(blocks_dir)){

  tictoc::tic()

  # -------------------------
  # Section 1: Project setup
  # -------------------------


  if(length(blocks_dir) == 1){
    i <- 1
  }

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

  # -----------------------------
  # Section 2: Output directories
  # -----------------------------

  raster_output <- glue::glue('{proj_dir}/output/raster')
  vector_output <- glue::glue('{proj_dir}/output/vector')

  # -------------------
  # Section 3: Treetops
  # -------------------

  # Make sure we load correct CHM

  chm_file <- stringr::str_subset(list.files(paste0(raster_output, '/chm'), pattern = '.tif$',full.names = T), pattern = 'smooth')

  if(length(chm_file) != 1){

    print(glue::glue("{length(chm_file)} CHM files matched expected input in {raster_output}/chm"))

    for(k in 1:length(chm_file)){
      print(glue::glue('{k}.  {chm_file[k]}'))
    }

    repeat {
      user_choice <- as.integer(readline("Enter the number corresponding to the correct CHM file: "))
      if (user_choice >= 1 && user_choice <= length(chm_file)) {
        break
      } else {
        cat("Invalid choice. Please enter a number between 1 and", length(chm_file), ".\n")
      }
    }

    chm_file <- chm_file[user_choice]
    cat(glue("\nSelected CHM file: {chm_file}\n"))
  }

  chm <- rast(chm_file)

  print(glue::glue('{chm_file} loaded as CHM, beginning local maxima finding...'))

  # Fixed window size (2m)
  tictoc::tic()
  lmf_ws2 <- locate_trees(chm, lmf(ws = 2, hmin = 5), uniqueness = 'incremental')
  elapsed_time_ws2 <- tictoc::toc(log = TRUE)
  print(glue::glue("Processed fixed window size (2m) treetops for {acq} in {elapsed_time_ws2$toc - elapsed_time_ws2$tic} seconds"))

  # Lmf auto
  tictoc::tic()
  lmf_auto <- locate_trees(chm, lidRplugins::lmfauto() , uniqueness = 'incremental')
  elapsed_time_auto <- tictoc::toc(log = TRUE)
  print(glue::glue("Processed lmf auto treetops for {acq} in {elapsed_time_auto$toc - elapsed_time_auto$tic} seconds"))

  # Variable window size
  tictoc::tic()
  f <- function(x) {x * 0.07 + 1}
  lmf_v <- locate_trees(chm, lmf(f, hmin = 5), uniqueness = 'incremental')
  elapsed_time_v <- tictoc::toc(log = TRUE)
  print(glue::glue("Processed variable window size treetops for {acq} in {elapsed_time_v$toc - elapsed_time_v$tic} seconds"))

  print(glue::glue("Summary for {acq} treetop detection:"))
  print(glue::glue("Fixed window size (2m) treetops: {nrow(lmf_ws2)} detected in {elapsed_time_ws2$toc - elapsed_time_ws2$tic} seconds"))
  print(glue::glue("LMF auto treetops: {nrow(lmf_auto)} detected in {elapsed_time_auto$toc - elapsed_time_auto$tic} seconds"))
  print(glue::glue("Variable window size treetops: {nrow(lmf_v)} detected in {elapsed_time_v$toc - elapsed_time_v$tic} seconds"))

  if(!dir.exists(glue::glue('{vector_output}/treetops'))){
    dir.create(glue::glue('{vector_output}/treetops'), recursive = T)
    print(glue::glue('Created tree top output directory for {acq}'))
  }

  sf::st_write(lmf_ws2, glue::glue('{vector_output}/treetops/{acq}_lmf_ws2_ttops.gpkg'), append = FALSE)
  sf::st_write(lmf_auto, glue::glue('{vector_output}/treetops/{acq}_lmf_auto_ttops.gpkg'), append = FALSE)
  sf::st_write(lmf_v, glue::glue('{vector_output}/treetops/{acq}_lmf_v_ttops.gpkg'), append = FALSE)

  print(glue::glue('Found and wrote {nrow(lmf_ws2)} ws2 treetops, {nrow(lmf_auto)} auto treetops, and {nrow(lmf_v)} variable ws treetops in {acq}'))
  tictoc::toc()
}

# ---- Crown Segmentation ----


for(i in 1:length(blocks_dir)){

tictoc::tic()

# ---- Project setup ----

proj_dir <- blocks_dir[i]

if(stringr::str_detect(basename(proj_dir), pattern = 'DAP') | is_dap){
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

chm_file <- stringr::str_subset(list.files(paste0(raster_output, '/chm'), pattern = '.tif$',full.names = T), pattern = 'smooth')


# Make sure we load correct CHM

if(length(chm_file) != 1){

  print(glue::glue("{length(chm_file)} CHM files matched expected input in {raster_output}/chm"))

  for(k in 1:length(chm_file)){
    print(glue::glue('{k}.  {chm_file[k]}'))
  }

  repeat {
    user_choice <- as.integer(readline("Enter the number corresponding to the correct CHM file: "))
    if (user_choice >= 1 && user_choice <= length(chm_file)) {
      break
    } else {
      cat("Invalid choice. Please enter a number between 1 and", length(chm_file), ".\n")
    }
  }

  chm_file <- chm_file[user_choice]
  cat(glue("\nSelected CHM file: {chm_file}\n"))
}

chm <- terra::rast(chm_file)
# chm_matrix <- as.matrix(chm)
# chm_r <- rast(chm_matrix)
# crs(chm_r) <- crs(chm)
# chm <- chm_r
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



#
# if(make_tiles == T){
#
# # Tile CHM
#
#   make.boxes = function(rstr, # Input raster
#                         size, # Box size (in m)
#                         overlap = 30, # How much do you want the boxes to overlap (in m)?
#                         progress.bar = T, # Progress bar? Need pbapply to be installed
#                         keep.all = F){# Keep all of the boxes or just the ones that contain valid pixels?
#
#     # Install packages if they aren't already installed
#     if(length(find.package("pbapply", quiet = T)) == 0){cat("\nInstalling pbapply"); install.packages("pbapply")}
#     if(length(find.package("exactextractr", quiet = T)) == 0){cat("\nInstalling exactextractr"); install.packages("exactextractr")}
#     if(length(find.package("terra", quiet = T)) == 0){cat("\nInstalling terra"); install.packages("terra")}
#     if(length(find.package("stringr", quiet = T)) == 0){cat("\nInstalling stringr"); install.packages("stringr")}
#     if(length(find.package("sf", quiet = T)) == 0){cat("\nInstalling sf"); install.packages("sf")}
#     if(length(find.package("purrr", quiet = T)) == 0){cat("\nInstalling purrr"); install.packages("purrr")}
#
#     if(progress.bar == F){pbapply::pboptions(type = "none")} # For the losers who don't want a progress bar
#
#     # Make a matrix of all possible boxes in the spatial extent of the raster
#     xv = seq(xmin(rstr), xmax(rstr), by = size)
#     yv = seq(ymin(rstr), ymax(rstr), by = size)
#
#     box.ctrs = expand.grid(xv, yv) |> as.matrix()
#
#     # Pad the box names so that they are the same length
#     n.pad = nrow(box.ctrs) |> stringr::str_length()
#
#     # Make all of the boxes - obviously if you have many boxes to make, this might take a bit
#     boxes = pbapply::pbapply(box.ctrs, 1, FUN = function(rw){
#       x = rw[1]
#       y = rw[2]
#
#       n.name = which(box.ctrs[,1] == x & box.ctrs[,2] == y) |>
#         stringr::str_pad(width = n.pad, side = "left", pad = "0")
#
#       # Make a box
#       box.sub <- cbind(c(x - overlap, x - overlap, x + size + overlap, x + size + overlap),
#                        c(y - overlap, y + size + overlap, y + size + overlap, y - overlap)) |>
#         terra::vect(type="polygons", crs = crs(rstr))
#
#       box.sub$name = n.name
#
#       # If you only want to keep boxes with valid pixels (helpful for non-rectangular rasters)
#       if(keep.all == F){
#         # See if the box overlaps the raster
#         test1 = sf::st_as_sf(box.sub) %>% exact_extract(rstr, .) |> dplyr::bind_rows() |> dplyr::pull(value) |> na.omit()
#
#         if(length(test1) > 0){return(box.sub)}
#         else(return())
#       }
#
#       # If you want to return all boxes, including potentially blank ones
#       return(box.sub)
#
#     })
#     return(vect(purrr::compact(boxes)))
#   }
#
#   print(glue:glue('Tiling CHM for {acq}'))
#
#   bxs <- make.boxes(chm, size = tile_size)
#
#   # Save CHM Tiles
#
#   tile_dir <- glue::glue('{raster_output}/chm/tile')
#
#   if(!dir.exists(tile_dir)){
#     dir.create(tile_dir)
#     print(glue::glue('Created directory for {size}m CHM tiles: {tile_dir}'))
#     }
#
#   pblapply(1:length(bxs), FUN = function(b){
#     nm <- bxs$name[b]
#
#     # What will be the output file path?
#     out.fp = glue::glue("{tile_dir}/{basename(tools::file_path_sans_ext(chm_file))}_{nm}.tif")
#
#     # Skip it if it already exists
#     if(file.exists(out.fp)){return(NA)}
#
#     r.sub <- terra::crop(chm, bxs[b]) # This can also be a stack if needed, or whatever you need to process something on
#
#     v.sub <- st_intersection(st_as_sf(ttops), st_as_sf(bxs[b,]))
#
#     r.out = silvtools::crown_mask(chunk = r.sub, ttops = v.sub,
#                                   crown_height_threshold = 0.25, vis = FALSE)
#
#     terra::writeRaster(r.sub, out.fp, overwrite = T)
#     return(NA)
#   })
#
#   plot(chm, ext = ext(bxs))
#   plot(st_geometry(st_as_sf(bxs)), add = T)
#   text(bxs, bxs$name)
#
# }


# ---- Segment Crowns ----

# if(run_parallel == T){
#   plan(multiprocess, workers = num_cores)
# }
#
# print(glue::glue('Begnning Crown Segmentation for {acq}'))
#
# chm_tiles <- list.files(tile_dir, full.names = T, pattern = '.tif$')
# crown_tiles <- glue::glue('{raster_output}/crowns/tile')
#
# if(!dir.exists(crown_tiles)){
#   dir.create(crown_tiles)
#   print(glue::glue('Created directory for Crown tiles: {crown_tiles}'))
# }
#
# b <- 100
#
# pblapply(1:length(chm_tiles), FUN = function(b){
#
#   chm_tile <- terra::rast(chm_tiles[b])
#
#   # What will be the output file path?
#   out.fp = glue::glue("{crown_tiles}/{basename(tools::file_path_sans_ext(chm_tiles[b]))}_{b}.tif")
#
#   # Skip it if it already exists
#   if(file.exists(out.fp)){return(NA)}
#
#   # Run some function on the subset
#   r.out = silvtools::crown_mask(chunk = chm, ttops = ttops_ws2, crown_height_threshold = 0.25, vis = FALSE)
#
#   terra::writeRaster(r.sub, out.fp, overwrite = T)
#   return(NA)
# })
#

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
  silvtools::convert_multi_to_single_polygons(polygons = ., fill_holes = TRUE)

print('Cleaning lmfauto( ) crowns...')
# lmfauto( )
crowns_auto_p <- sf::st_as_sf(terra::as.polygons(crowns_auto)) %>%
  silvtools::convert_multi_to_single_polygons(polygons = ., fill_holes = TRUE)

print('Cleaning lmfv( ) crowns...')
# lmfv()
crowns_v_p <- sf::st_as_sf(terra::as.polygons(crowns_v)) %>%
  silvtools::convert_multi_to_single_polygons(polygons = ., fill_holes = TRUE)

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
