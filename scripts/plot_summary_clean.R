# Match Postex and Lidar Stems - Assess Crown Relationships
# May 2023
# Liam Irwin

# ---- Load Packages ----

library(sf)
library(dplyr)
library(tidyverse)
library(terra)
library(tmap)
library(plotly)
library(silvtools)
library(lidR)
# ---- Processing Setup ----
plots_dir <- 'G:/Quesnel_2022/plot_summary/plots'
process_dir <- 'G:/Quesnel_2022/process'
blocks_dir <- list.dirs(plots_dir, recursive = FALSE)

# Area around plot stems to crop CHM to
chm_radius <- 10
coords <- 26910
# Output lists for matching and plotting
matches <- list()
tmap_plots <- list()
# Generate tmap validation plots? - Takes some time
plot_result <- FALSE
# Should irradiance rasters be generated for each plot? - Takes some time
generate_insol = FALSE
# Extract and include solar irradiance?
extract_insol = T
# Extract and include topographic wetness?
extract_twi = T
# Generate Pairwise competition index?
generate_heygi <- T
# Generate Area Potentially Avaliable Index
generate_apa <- T
# Use manually adjusted maxima
use_manual_maxima <- T


comp_input <- 'vol_concave'
maxR <- 6
# Match Reference Stems with Closest Detected Treetops

for(i in 1:length(blocks_dir)){

  # ---------- Set up directories and output paths ----------
  proj_dir <- blocks_dir[i]
  plot_name <- basename(proj_dir)
  block_name <- str_extract(plot_name, "CT\\d+")
  block_dir <- glue::glue('{process_dir}/{block_name}')
  raster_output <- glue::glue('{block_dir}/output/raster')
  vector_output <- glue::glue('{block_dir}/output/vector')

  print(glue::glue('Loaded plot directory for {plot_name}'))

  # ---------- Setup directories ----------
  silvtools::setup_als_dirs(proj_dir)

  # ---------- Load and filter reference trees ----------
  area_stems <- st_read(glue::glue('{dirname(file.path(plots_dir))}/postex_w_cores_2022.gpkg'), quiet = TRUE)
  plot_stems <- area_stems %>% filter(PlotID == plot_name)
  print(glue::glue('Filtered {nrow(plot_stems)} reference trees for {plot_name} from overall {nrow(area_stems)} in dataset...'))

  # ---------- Filter tree cores ----------
  core_stems <- plot_stems %>% filter(!is.na(mean_bai_5))
  print(glue::glue('Of the {nrow(plot_stems)} reference trees for {plot_name}, {nrow(core_stems)} had tree cores taken...'))

  # ---------- Load CHM, treetops and crowns files ----------
  chm_file <- stringr::str_subset(list.files(glue::glue('{block_dir}/output/raster/chm'), pattern = '.tif$',full.names = T), pattern = 'fill')

  # Choose if manual maxima or automatically detected are used
  if (use_manual_maxima == TRUE) {
    ttops_file <- 'G:/Quesnel_2022/plot_summary/ct_ttops_ws1.gpkg'
    print(glue::glue('Loaded manually adjusted tree tops'))
  } else{
    ttops_file <-
      stringr::str_subset(list.files(
        glue::glue('{block_dir}/output/vector/treetops'),
        pattern = '.gpkg',
        full.names = T
      ),
      pattern = 'ashape')
  }

  crown_file <- stringr::str_subset(list.files(glue::glue('{block_dir}/output/vector/crowns'), pattern = '.shp$', full.names = T), pattern = 'ws2')

  # If multiple CHM files found, allow user to select correct one
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
    cat(glue::glue("\nSelected CHM file: {chm_file}\n"))
  }

  chm <- rast(chm_file)
  ttops <- st_read(ttops_file, quiet = TRUE)
  crowns <- st_read(crown_file, quiet = TRUE)


  # Define a bounding box based on the plot_stems data, adjusting it to include a radius around the plot
  bbox <- plot_stems %>% st_bbox() %>% st_as_sfc() %>% terra::vect() %>% terra::ext()
  # Here we are expanding the bounding box by the value of chm_radius in each direction
  bbox[1,] <- bbox[1,] - chm_radius
  bbox[2,] <- bbox[2,] + chm_radius
  bbox[3,] <- bbox[3,] - chm_radius
  bbox[4,] <- bbox[4,] + chm_radius
  # The bounding box is converted to an sf object
  bbox_sf <- st_as_sf(vect(bbox))
  # Setting the CRS of the bounding box to match that of the treetop data (ttops)
  st_crs(bbox_sf) <- st_crs(ttops)
  # Crop the canopy height model (chm) to the extent of the bounding box
  plot_chm <- terra::crop(chm, bbox)
  # Intersect the treetops and crowns data with the bounding box to get only data within the bounding box
  plot_ttops <- ttops %>% st_intersection(bbox_sf)
  plot_crowns <- crowns %>% st_intersection(bbox_sf)
  # Apply the crown_mask function from the silvtools package, convert the output to polygons and then to an sf object.
  # Then convert multipolygons to single polygons.
  plot_crowns <- silvtools::crown_mask(chunk = plot_chm, ttops = ttops, crown_height_threshold = 0.25, vis = TRUE) %>% terra::as.polygons(.) %>% sf::st_as_sf(.) %>% silvtools::convert_multi_to_single_polygons(polygons = ., fill_holes = TRUE)
  # Set up a catalog of normalized lidar data
  ctg_norm <- catalog(glue::glue('{block_dir}/input/las/norm'))
  # Clip lidar data to the extent of the bounding box
  plot_las <- clip_roi(ctg_norm, bbox_sf)
  print('Clipped normalized las to plot extent')

  apply_treeid_to_las <- function(chunk, crowns){
    tictoc::tic()
    if('LAS' %in% class(chunk)){
      las <- chunk
    } else{
      las <- readLAS(chunk)
    }
    box <- st_as_sfc(st_bbox(chunk))
    if (is.empty(las)) return(NULL)
    # if (!'treeID' %in% names(las@data)) return(NULL)
    box_buf <- st_buffer(box, dist = 5)
    # Select crowns relevant to chunk of interest
    chunk_crowns <- crowns[sf::st_intersects(crowns, box_buf, sparse = FALSE),]
    if(nrow(chunk_crowns) == 0) return(NULL)
    if("Z" %in% names(chunk_crowns)){
      names(chunk_crowns) <- c('treeID','geometry')
    }
    # Merge tree IDs with las
    tree_las <- merge_spatial(las, chunk_crowns, attribute = 'treeID')
    tree_las = add_lasattribute(tree_las, name="treeID", desc="ID of a tree")
    tree_las <- filter_poi(tree_las, !is.na(treeID))
    glue::glue('Merged {nrow(chunk_crowns)} crowns with {length(las@data$Z)} points in chunk')
    tictoc::toc()
    return(tree_las)
  }

  tree_las <- apply_treeid_to_las(plot_las, plot_crowns)
  print('Applied crown tree IDs to plot las')
  plot_ashapes <- get_alphashape_metrics(tree_las, RGB = TRUE)
  plot_ttops <- merge(plot_ttops, plot_ashapes, by = 'treeID')

  if(generate_insol == TRUE){
    # ---------- Setup Solar Position Dataframe ----------
    start_date <- as.POSIXct("2022-05-15 00:00:00", tz = "America/Los_Angeles")
    end_date <- as.POSIXct("2022-09-15 00:00:00", tz = "America/Los_Angeles")
    interval <- '1 hour'
    lat <- 53.371759251761766
    lon <- -122.76808019195623
    time_df <- get_solar_pos(start_date, end_date, interval, lat, lon)
    # Filter positions of interest
    solar_pos <- time_df %>% filter(wday == 2 & as.numeric(hour) %in% c(9,10,11,12,13,14,15,16))
    # Find the required file
    dsm_file <- list.files(glue::glue('{block_dir}/output/raster/dsm'), pattern = 'fill_p2r_0.1m.tif$', full.names = T)
    dsm <- rast(dsm_file) %>% terra::crop(., bbox)


    # Call solar simulator
    solar_simulator(dsm_file = dsm,
                    proj_dir = proj_dir,
                    site_name = plot_name,
                    sun_pos = solar_pos,
                    num_cores = 1
    )

    create_animation(glue::glue('{proj_dir}/output/raster/irradiance/{plot_name}'), sitename = plot_name)

    # Define function to calculate mean irradiance
    mean_irradiance <- function(proj_dir, sitename){
      irradiance_dir <- paste0(proj_dir, '/output/raster/irradiance/', sitename)
      irradiance_files <- list.files(irradiance_dir, pattern = '.tif$', full.names = T)
      print(glue::glue('Found {length(irradiance_files)} irradiance rasters in folder for {sitename}'))
      r <- terra::rast(irradiance_files)
      print(glue::glue('Calculating mean value across {length(irradiance_files)} rasters'))
      r_mean <- terra::app(r, mean)
      terra::writeRaster(r_mean, paste0(proj_dir, '/output/raster/irradiance/', sitename, '_mean_irradiance.tif'), overwrite = T)
    }

    # Call the mean_irradiance function
    mean_irradiance(proj_dir = proj_dir, sitename = plot_name)

  }

  # ---------- Extract Solar Irradiance ----------
  insol_file <- stringr::str_subset(list.files(glue::glue('{proj_dir}/output/raster/irradiance'), pattern = '.tif', full.names = T), pattern = 'mean')
  # If mean irradiance raster found and extraction requested, perform extraction
  if(length(insol_file) == 0){

    print(glue::glue('No mean irradiance raster was found in {proj_dir}/output/raster/irradiance, consider rerunning with generate_insol = TRUE'))

  }


  if(extract_insol == TRUE & length(insol_file) == 1){

    print(glue::glue('Found a mean irradiance raster for {basename(proj_dir)} clipping to crown extents'))

    insol <- rast(insol_file) %>% terra::crop(., bbox)

    insol_crowns <- exactextractr::exact_extract(insol, plot_crowns, fun = 'mean', append_cols = TRUE)

    names(insol_crowns) <- c('treeID','irr_mean')

    plot_ttops <- merge(plot_ttops, insol_crowns, by = 'treeID')

  }

  # Extract Topographic Wetness
  # If TWI raster found and extraction requested, perform extraction
  twi_file <- stringr::str_subset(list.files(glue::glue('{block_dir}/output/raster/topography'), pattern = '.tif$', full.names = T), pattern = 'twi')

  if(length(twi_file) == 0){

    print(glue::glue('No topographic wetness index raster was found in {block_dir}/output/raster/topography'))

  }

  if(extract_twi == TRUE & length(twi_file) == 1){

    print(glue::glue('Found a topographic wetness index raster for {basename(block_dir)} clipping to crown extents'))

    twi <- rast(twi_file) %>% terra::crop(., bbox)

    twi_crowns <- exactextractr::exact_extract(twi, plot_crowns, fun = 'mean', append_cols = TRUE)

    names(twi_crowns) <- c('treeID','twi_mean')

    plot_ttops <- merge(plot_ttops, twi_crowns, by = 'treeID')

  }


  if(generate_heygi == TRUE){

    heygi_ttops <- heygi_cindex(plot_ttops, comp_input = comp_input, maxR = maxR) %>%
      st_drop_geometry() %>% select(cindex, treeID)

    heygi_ttops$treeID <- plot_ttops$treeID
    plot_ttops <- merge(plot_ttops, heygi_ttops, by = 'treeID')

  }

  if(generate_apa == TRUE){
    apa_ttops <- apa_cindex(plot_ttops, comp_input = comp_input) %>%
      st_drop_geometry() %>%
      select(aindex, afree, treeID)

    apa_ttops$treeID <- plot_ttops$treeID

    plot_ttops <- merge(plot_ttops, apa_ttops, by = 'treeID')

  }


  # Perform tree matching and score results
  matches[[i]] <- tree_matching(plot_stems, plot_ttops, plot_name)
  score <- tree_matching_scores(matches[[i]]) %>% mutate(PlotID = plot_name)

  # Return the score
  score

}


matches_df <- do.call(rbind, matches)
matches_df <- matches_df %>% filter(!is.na(TreeNum))

core_trees <- matches_df %>% dplyr::filter(!is.na(mean_bai_5))
