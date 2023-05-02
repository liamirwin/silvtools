# Pixel Metric Mass Generation


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

# ---- Processing switches ----
# ULS or DAP?
is_dap <- FALSE
# Run in parallel?
run_parallel <- T
num_cores <- 3L

chunk_buf <- 5
# Calculate Metrics?
make_mets <- T
met_res <- 10
# Size of voxel for voxel metrics
vox_size <- 1
# Is ALS?
is_als <- F
is_mls <- F
# List directories (each is one acquisiton of ULS/DAP)
blocks_dir <- 'I:/NZ_2023/Cass/ULS_Subset'
################################################################################
# START BUTTON
################################################################################

for(i in 1:length(blocks_dir)){

  timings_df <- data.frame(Metric = character(), Time = numeric(), stringsAsFactors = FALSE)

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

  # ---- Initialize Parallel Processing ----

  if(run_parallel == TRUE){
    future::plan(future::multisession, workers = num_cores)
    print(glue::glue('Parallel processing initiated using {num_cores} cores'))
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

    tictoc::tic()

    print(glue::glue('Generating basic metrics for {acq} at {met_res}m'))
    basic <- pixel_metrics(ctg_norm, ~lidRmetrics::metrics_basic(Z), res = met_res)
    terra::writeRaster(x = basic, filename = glue::glue('{raster_output}/metrics/{acq}_basic_{met_res}m.tif'), overwrite = TRUE)

    time_basic <- tictoc::toc(log = TRUE)
    timings_df <- rbind(timings_df, data.frame(Metric = "Basic", Time = time_basic$toc - time_basic$tic, stringsAsFactors = FALSE))

    # Percentiles metrics suite
    # -------------------------

    tictoc::tic()

    print(glue::glue('Generating percentiles metrics for {acq} at {met_res}m'))
    percentiles_metrics <- pixel_metrics(ctg_norm, ~lidRmetrics::metrics_percentiles(Z), res = met_res)
    terra::writeRaster(x = percentiles_metrics, filename = glue::glue('{raster_output}/metrics/{acq}_percentiles_{met_res}m.tif'), overwrite = TRUE)

    time_percentiles <- tictoc::toc(log = TRUE)
    timings_df <- rbind(timings_df, data.frame(Metric = "Percentiles", Time = time_percentiles$toc - time_percentiles$tic, stringsAsFactors = FALSE))


    # Percentage of returns above a threshold metrics suite
    # -----------------------------------------------------
    tictoc::tic()

    print(glue::glue('Generating percabove metrics for {acq} at {met_res}m'))
    percabove_metrics <- pixel_metrics(ctg_norm, ~lidRmetrics::metrics_percabove(Z), res = met_res)
    terra::writeRaster(x = percabove_metrics, filename = glue::glue('{raster_output}/metrics/{acq}_percabove_{met_res}m.tif'), overwrite = TRUE)

    time_percabove <- tictoc::toc(log = TRUE)
    timings_df <- rbind(timings_df, data.frame(Metric = "Percent Above", Time = time_percabove$toc - time_percabove$tic, stringsAsFactors = FALSE))


    # Dispersion metrics suite
    # ------------------------
    tictoc::tic()

    print(glue::glue('Generating dispersion metrics for {acq} at {met_res}m'))
    dispersion_metrics <- pixel_metrics(ctg_norm, ~lidRmetrics::metrics_dispersion(Z), res = met_res)
    terra::writeRaster(x = dispersion_metrics, filename = glue::glue('{raster_output}/metrics/{acq}_dispersion_{met_res}m.tif'), overwrite = TRUE)

    time_dispersion <- tictoc::toc(log = TRUE)
    timings_df <- rbind(timings_df, data.frame(Metric = "Dispersion", Time = time_dispersion$toc - time_dispersion$tic, stringsAsFactors = FALSE))


    # Canopy density metrics suite
    # ----------------------------

    tictoc::tic()

    print(glue::glue('Generating canopy density metrics for {acq} at {met_res}m'))
    canopydensity_metrics <- pixel_metrics(ctg_norm, ~lidRmetrics::metrics_canopydensity(Z), res = met_res)
    terra::writeRaster(x = canopydensity_metrics, filename = glue::glue('{raster_output}/metrics/{acq}_canopydensity_{met_res}m.tif'), overwrite = TRUE)

    time_canopydensity <- tictoc::toc(log = TRUE)
    timings_df <- rbind(timings_df, data.frame(Metric = "Canopy Density", Time = time_canopydensity$toc - time_canopydensity$tic, stringsAsFactors = FALSE))


    # L-moments metrics suite
    # ------------------------

    tictoc::tic()

    # Check if Lmoments package is installed and loaded
    if (!"Lmoments" %in% rownames(installed.packages())) {
      cat("Lmoments package is not installed.\n")
      install_option <- readline(prompt = "Do you want to install Lmoments package now? (y/n): ")

      if (tolower(install_option) == "y") {
        install.packages("Lmoments")
        library(Lmoments)
      } else {
        cat("Skipping L-moments metrics generation as the Lmoments package is not installed.\n")
      }
    }

    # If Lmoments package is installed and loaded, generate L-moments metrics
    if ("Lmoments" %in% .packages()) {
      library(Lmoments)
      print(glue::glue('Generating L-moments metrics for {acq} at {met_res}m'))
      Lmoments_metrics <- pixel_metrics(ctg_norm, ~lidRmetrics::metrics_Lmoments(Z), res = met_res)
      terra::writeRaster(x = Lmoments_metrics, filename = glue::glue('{raster_output}/metrics/{acq}_Lmoments_{met_res}m.tif'), overwrite = TRUE)
    } else {
      cat("Skipping L-moments metrics generation as the Lmoments package is not installed.\n")
    }

    time_Lmoments <- tictoc::toc(log = TRUE)
    timings_df <- rbind(timings_df, data.frame(Metric = "Lmoments", Time = time_Lmoments$toc - time_Lmoments$tic, stringsAsFactors = FALSE))



    # Leaf Area Density metrics suite
    # -------------------------------

    tictoc::tic()

    print(glue::glue('Generating leaf area density metrics for {acq} at {met_res}m'))
    lad_metrics <- pixel_metrics(ctg_norm, ~lidRmetrics::metrics_lad(Z), res = met_res)
    terra::writeRaster(x = lad_metrics, filename = glue::glue('{raster_output}/metrics/{acq}_lad_{met_res}m.tif'), overwrite = TRUE)

    time_lad <- tictoc::toc(log = TRUE)
    timings_df <- rbind(timings_df, data.frame(Metric = "Leaf Area Density", Time = time_lad$toc - time_lad$tic, stringsAsFactors = FALSE))


    # Interval metrics suite
    # ----------------------
    tictoc::tic()

    print(glue::glue('Generating interval metrics for {acq} at {met_res}m'))
    interval_metrics <- pixel_metrics(ctg_norm, ~lidRmetrics::metrics_interval(Z), res = met_res)
    terra::writeRaster(x = interval_metrics, filename = glue::glue('{raster_output}/metrics/{acq}_interval_{met_res}m.tif'), overwrite = TRUE)

    time_interval <- tictoc::toc(log = TRUE)
    timings_df <- rbind(timings_df, data.frame(Metric = "Intervals", Time = time_interval$toc - time_interval$tic, stringsAsFactors = FALSE))


    # Rumple metrics suite
    # --------------------

    tictoc::tic()

    # Check if geometry package is installed and loaded
    if (!"geometry" %in% rownames(installed.packages())) {
      cat("geometry package is not installed.\n")
      install_option <- readline(prompt = "Do you want to install geometry package now? (yes/no): ")

      if (tolower(install_option) == "yes") {
        install.packages("geometry")
        library(geometry)
      } else {
        cat("Skipping rumple metrics generation as the geometry package is not installed.\n")
      }
    }

    # If geometry package is installed and loaded, generate rumple metrics
    if ("geometry" %in% .packages()) {
      library(geometry)
      print(glue::glue('Generating rumple metrics for {acq} at {met_res}m'))
      rumple_metrics <- pixel_metrics(ctg_norm, ~lidRmetrics::metrics_rumple(x = X, y = Y, z = Z, pixel_size = met_res/10), res = met_res)
      terra::writeRaster(x = rumple_metrics, filename = glue::glue('{raster_output}/metrics/{acq}_rumple_{met_res}m.tif'), overwrite = TRUE)
    } else {
      cat("Skipping rumple metrics generation as the geometry package is not installed.\n")
    }

    time_rumple <- tictoc::toc(log = TRUE)
    timings_df <- rbind(timings_df, data.frame(Metric = "Rumple", Time = time_rumple$toc - time_rumple$tic, stringsAsFactors = FALSE))


    # Voxel metrics suite
    # -------------------
    tictoc::tic()

    print(glue::glue('Generating voxel metrics for {acq} at {met_res}m'))
    if(!exists('vox_size')){
      vox_size <- 1
    }
    voxel_metrics <- pixel_metrics(ctg_norm, ~lidRmetrics::metrics_voxels(x= X, y = Y, z = Z, vox_size = vox_size), res = met_res)
    terra::writeRaster(x = voxel_metrics, filename = glue::glue('{raster_output}/metrics/{acq}_voxels_{met_res}m.tif'), overwrite = TRUE)

    time_voxel <- tictoc::toc(log = TRUE)
    timings_df <- rbind(timings_df, data.frame(Metric = "Voxel", Time = time_voxel$toc - time_voxel$tic, stringsAsFactors = FALSE))


    # KDE metrics suite
    # -----------------

    tictoc::tic()

    print(glue::glue('Generating kernel density estimation (KDE) metrics for {acq} at {met_res}m'))
    kde_metrics <- pixel_metrics(ctg_norm, ~lidRmetrics::metrics_kde(Z), res = met_res)
    terra::writeRaster(x = kde_metrics, filename = glue::glue('{raster_output}/metrics/{acq}_kde_{met_res}m.tif'), overwrite = TRUE)

    time_KDE <- tictoc::toc(log = TRUE)
    timings_df <- rbind(timings_df, data.frame(Metric = "KDE", Time = time_KDE$toc - time_KDE$tic, stringsAsFactors = FALSE))


    ggplot(timings_df, aes(x = Metric, y = Time)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      labs(title = "Time Taken for Each Metric Suite", x = "Metric Suite", y = "Time (seconds)") +
      theme_minimal()


  }

  print(glue::glue('Finished point cloud processing for {acq}'))
  tictoc::toc()

}
