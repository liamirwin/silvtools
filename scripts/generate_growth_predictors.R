# Match Postex and Lidar Stems - Assess Crown Relationships
# May 2023
# Liam Irwin

# ---- Load Packages ----
library(sf) # For vector data processing
library(tidyverse) # For data manipulation
library(terra) # For raster processing
library(tmap) # For thematic maps
library(plotly) # For interactive plots
library(silvtools) # For custom functions
library(lidR) # For lidar processing
library(siplab) # For heygi and afree calculations
library(nngeo) # For nearest neighbour analysis
# ---- Processing Setup ----
plots_dir <- 'G:/Quesnel_2022/plot_summary/plots'
process_dir <- 'G:/Quesnel_2022/process'
blocks_dir <- list.dirs(plots_dir, recursive = FALSE)
# Area around plot stems to crop CHM to
chm_radius <- 10
coords <- 26910
# Output lists for matching and plotting
matches <- list()
matches_sf <- list()
scores <- list()
tmap_plots <- list()
# Should results be written to disk?
write_results <- TRUE
# Generate tmap validation plots? - Takes some time
plot_result <- FALSE
# Should irradiance rasters be generated for each plot? - Takes some time
generate_insol = FALSE
# Extract and include solar irradiance?
extract_insol = TRUE
# Extract and include topographic wetness?
extract_twi = TRUE
# Generate Pairwise competition index?
generate_heygi <- TRUE
# Generate Area Potentially Avaliable Index
generate_apa <- TRUE
# Calculate distance to nearest neighbouring tree top
generate_nn <- TRUE
# Decimate points before alphashapes are computed?
decimate <- FALSE
v_res <- 1
v_n <- 3
# Set crown height threshold to limit watershed segmentation extent
crown_height_threshold <- 0.70

# Set up competition index inputs
# Tree size metric(s) used to compute cindex
comp_input <- c('vol_concave','vol_convex')
maxR <- c(2,3,4,5,6,7,8,9,10,11,12,13,14,15)
# Match Reference Stems with Closest Detected Treetops

# Loop over all directories
for(i in 1:length(blocks_dir)){

  # Set up directories and output paths
  proj_dir <- blocks_dir[i]
  plot_name <- basename(proj_dir)
  block_name <- str_extract(plot_name, "CT\\d+")
  block_dir <- glue::glue('{process_dir}/{block_name}')
  raster_output <- glue::glue('{block_dir}/output/raster')
  vector_output <- glue::glue('{block_dir}/output/vector')

  print(glue::glue('Loaded plot directory for {plot_name}'))

  # Set up ALS directories with silvtools package
  silvtools::setup_als_dirs(proj_dir)

  # Load and filter reference trees
  area_stems <- st_read(glue::glue('{dirname(file.path(plots_dir))}/postex_w_cores_2022.gpkg'), quiet = TRUE)
  plot_stems <- area_stems %>% filter(PlotID == plot_name)
  print(glue::glue('Filtered {nrow(plot_stems)} reference trees for {plot_name} from overall {nrow(area_stems)} in dataset...'))

  # Filter tree cores
  core_stems <- plot_stems %>% filter(!is.na(mean_bai_5))
  print(glue::glue('Of the {nrow(plot_stems)} reference trees for {plot_name}, {nrow(core_stems)} had tree cores taken...'))

  # Load CHM, treetops and crowns files
  chm_file <- stringr::str_subset(list.files(glue::glue('{block_dir}/output/raster/chm'), pattern = '.tif$',full.names = T), pattern = 'fill')
  #ttops_file <- stringr::str_subset(list.files(glue::glue('{block_dir}/output/vector/treetops'), pattern = '.gpkg', full.names = T), pattern = 'ashape')
  ttops_file <- 'G:/Quesnel_2022/plot_summary/ct_ttops_ws1.gpkg'
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


  # Set up bounding box and adjust to include radius around plot
  bbox <- plot_stems %>% st_bbox() %>% st_as_sfc() %>% terra::vect() %>%
    terra::ext()
  bbox[1,] <- bbox[1,] - chm_radius
  bbox[2,] <- bbox[2,] + chm_radius
  bbox[3,] <- bbox[3,] - chm_radius
  bbox[4,] <- bbox[4,] + chm_radius
  bbox_sf <- st_as_sf(vect(bbox))
  st_crs(bbox_sf) <- st_crs(ttops)

  # Crop CHM, treetops and crowns to plot extent
  plot_chm <- terra::crop(chm, bbox)
  plot_ttops <- ttops %>% st_intersection(bbox_sf)
  plot_crowns <- crowns %>% st_intersection(bbox_sf)
  # Generate watershed segmentations limited to crown height threshold
  plot_crowns <- silvtools::crown_mask(chunk = plot_chm,
                                       ttops = plot_ttops,
                                       crown_height_threshold = crown_height_threshold,
                                       vis = FALSE) %>%
    terra::as.polygons(.) %>% sf::st_as_sf(.) %>%
  silvtools::convert_multi_to_single_polygons(polygons = ., fill_holes = TRUE)
  ctg_norm <- catalog(glue::glue('{block_dir}/input/las/norm'))
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


  if(decimate == TRUE){
  tree_las <- decimate_points(tree_las, random_per_voxel(res = v_res, n = v_n))
  print(glue::glue('Decimated Tree LAS with 5cm voxel'))

  }

  plot_ashapes <- get_alphashape_metrics(tree_las)

  # Get the column names from plot_ashapes that are not in plot_ttops
  new_cols <- setdiff(names(plot_ashapes), names(plot_ttops))

  # Join plot_ttops with only the new columns from plot_ashapes
  plot_ttops <- plot_ttops %>%
    left_join(plot_ashapes %>% select(all_of(c("treeID", new_cols))), by = "treeID")


  if(generate_insol == TRUE){

    # Setup SOlar Position Dataframe
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

  # Extract Solar Irradiance
  insol_file <- stringr::str_subset(list.files(glue::glue('{proj_dir}/output/raster/irradiance'), pattern = '.tif', full.names = T), pattern = 'mean')
  # If mean irradiance raster found and extraction requested, perform extraction
  if(length(insol_file) == 0){
    print(glue::glue('No mean irradiance raster was found in {proj_dir}/output/raster/irradiance, consider rerunning with generate_insol = TRUE'))
    extract_insol <- FALSE
    print(glue::glue('Skipping insolation extraction since no raster was found'))

  }


  if(extract_insol == TRUE & length(insol_file) == 1){

    print(glue::glue('Found a mean irradiance raster for {basename(proj_dir)} clipping to crown extents'))

    insol <- rast(insol_file) %>% terra::crop(., bbox)

    insol_crowns <- exactextractr::exact_extract(insol, plot_crowns, fun = 'mean', append_cols = TRUE)

    names(insol_crowns) <- c('treeID','irr_mean')

    plot_ttops <- merge(plot_ttops, insol_crowns, by = 'treeID')

  }

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

  plot_ttops <- plot_ttops %>% filter(!is.na(X) | !is.na(Y))

  if (generate_heygi == TRUE) {
      # if more than one comp_input value generate unique metric for each
      for(k in 1:length(comp_input)){
        # Generate Heygi cindex for each comp_input value
        c_input <- comp_input[k]
        # Generate Heygi cindex for each maxR value
        for (r in 1:length(maxR)) {
          heygi_ttops <-
            heygi_cindex(plot_ttops, comp_input = comp_input, maxR = maxR[r]) %>%
            st_drop_geometry() %>%
            select(cindex, treeID) %>% filter(!is.na(treeID))

          colnames(heygi_ttops)[1] <- glue::glue('cindex_{c_input}_{maxR[r]}m')

          plot_ttops <- merge(plot_ttops, heygi_ttops, by = 'treeID')
          print(glue::glue('Calculated Heygi cindex based on {c_input} for maxR = {maxR[r]}m radius sphere of influence'))
        }
        print(glue::glue('Calculated Heygi cindex based on {c_input} for all maxR values'))
      }
  }

  if(generate_apa == TRUE){

    for(n in 1:length(comp_input)){
      # Generate APA cindex for each comp_input value
      afree_ttops <- apa_cindex(plot_ttops, comp_input = comp_input[n]) %>%
        st_drop_geometry() %>%
        select(afree)
      colnames(afree_ttops)[1] <- glue::glue('afree_{comp_input[n]}')
      plot_ttops <- cbind(plot_ttops, afree_ttops)
      print(glue::glue('Calculated APA cindex based on {comp_input[n]}'))
    }
  }

  if(generate_nn){
    # Generate nearest neighbour distance
    nn_ttops <- nngeo::st_nn(plot_ttops, plot_ttops, k = 2, returnDist = TRUE)
    nn_dist <- nn_ttops[[2]] %>% do.call(rbind, .) %>% as.data.frame() %>% select(V2) %>%
      rename(nn_dist = V2)
    plot_ttops$nn_dist <- nn_dist
    print('Calculated nearest neighbour distance')
  }

  # Perform tree matching and score results
  matches[[i]] <- tree_matching(plot_stems, plot_ttops, plot_name)

  matches_sf[[i]] <- matches[[i]] %>% filter(!is.na(X) &!is.na(Y)) %>%
    st_as_sf(coords = c('X', 'Y'), crs = st_crs(plot_ttops))

  scores[[i]] <- tree_matching_scores(matches[[i]]) %>% mutate(PlotID = plot_name)

  # Return the score
  score[[i]]

  # Write matches df locally
  if(write_results){
  dir.create(glue::glue('{proj_dir}/output/matches'), showWarnings = F)
  dir.create(glue::glue('{proj_dir}/output/matches/scores'), showWarnings = F)
  dir.create(glue::glue('{proj_dir}/output/matches/sf'), showWarnings = F)
  dir.create(glue::glue('{proj_dir}/output/matches/csv'), showWarnings = F)

  write.csv(scores[[i]], glue::glue('{proj_dir}/output/matches/scores/{plot_name}_scores_cht{crown_height_threshold}.csv'), append = F)
  st_write(matches_sf[[i]], glue::glue('{proj_dir}/output/matches/sf/{plot_name}_matches_cht{crown_height_threshold}.gpkg'), append = F)
  }

}



matches_df <- do.call(rbind, matches)
matches_df <- matches_df %>% filter(!is.na(TreeNum))

core_trees <- matches_df %>% dplyr::filter(!is.na(mean_bai_5) & !is.na(vol_concave))

write.csv(st_drop_geometry(core_trees),'G:/Quesnel_2022/Modelling/core_trees_fix_review.csv')
