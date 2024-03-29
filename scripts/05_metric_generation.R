# Pixel Metric Mass Generation

# ----- Load Packages ----- #
library(lidR)
library(future)
library(glue)
library(lidRmetrics)
library(terra)
library(stringr)
library(geometry) # Required for rumple metrics
library(Lmoments) # Required for Lmoments metrics
library(ggplot2)

# ---- Processing switches ----
# ULS or DAP?
is_dap <- FALSE
# Run in parallel?
run_parallel <- T
num_cores <- 6L
chunk_buf <- 5
# Pixel metric spatial resolution (m)
met_res <- 0.25
# Generate time consuming voxel metrics?
make_vox_mets <- F
# Size of voxels for voxel metrics
vox_size <- 1
# Is ALS?
is_als <- F
is_mls <- F
# List directories (each is one acquisiton of ULS/DAP)
blocks_dir <- 'I:/NZ_2023/Cass/ULS'
# Acqusition Name
acq <- 'Cass_ULS'

################################################################################
# START BUTTON
################################################################################

for(i in 1:length(blocks_dir)){

  timings_df <- data.frame(Metric = character(), Time = numeric(),
                           stringsAsFactors = FALSE)

  if(length(blocks_dir) == 1){
    i <- 1
  }

  proj_dir <- blocks_dir[i]

  # ----- Output directories -----

  raster_output <- glue::glue('{proj_dir}/output/raster')
  vector_output <- glue::glue('{proj_dir}/output/vector')


  # ---- Initialize Parallel Processing ----
  if(run_parallel == TRUE){
    future::plan(future::multisession, workers = num_cores)
    print(glue::glue('Parallel processing initiated using {num_cores} cores'))
  }

  # ---- Metrics -----

    if(!dir.exists(glue::glue('{raster_output}/metrics'))){
      dir.create(glue::glue('{raster_output}/metrics'), recursive = T)
      print(glue::glue('Created a directory for metrics at {raster_output}/metrics'))
    }

    met_res <- met_res
    ctg_norm <- catalog(glue::glue('{proj_dir}/input/las/norm'))
    opt_progress(ctg_norm) <- T
    opt_chunk_buffer(ctg_norm) <- chunk_buf
    opt_select(ctg_norm) <- "xyz"
    opt_filter(ctg_norm) <- "-drop_withheld -drop_z_below 0 -drop_z_above 25"
    opt_progress(ctg_norm) <- T


    # Basic metrics
    tictoc::tic()
    print(glue::glue('Generating basic metrics for {acq} at {met_res}m'))
    basic <- pixel_metrics(ctg_norm, ~lidRmetrics::metrics_basic(Z), res = met_res)
    terra::writeRaster(x = basic, filename = glue::glue('{raster_output}/metrics/{acq}_basic_{met_res}m.tif'), overwrite = TRUE)
    time_basic <- tictoc::toc(log = TRUE)
    timings_df <- rbind(timings_df, data.frame(Metric = "Basic", Time = time_basic$toc - time_basic$tic, stringsAsFactors = FALSE))

    # Percentiles metrics
    tictoc::tic()
    print(glue::glue('Generating percentiles metrics for {acq} at {met_res}m'))
    percentiles_metrics <- pixel_metrics(ctg_norm, ~lidRmetrics::metrics_percentiles(Z), res = met_res)
    terra::writeRaster(x = percentiles_metrics, filename = glue::glue('{raster_output}/metrics/{acq}_percentiles_{met_res}m.tif'), overwrite = TRUE)
    time_percentiles <- tictoc::toc(log = TRUE)
    timings_df <- rbind(timings_df, data.frame(Metric = "Percentiles", Time = time_percentiles$toc - time_percentiles$tic, stringsAsFactors = FALSE))


    # Percentage of returns above a threshold metrics
    tictoc::tic()
    print(glue::glue('Generating percabove metrics for {acq} at {met_res}m'))
    percabove_metrics <- pixel_metrics(ctg_norm, ~lidRmetrics::metrics_percabove(Z), res = met_res)
    terra::writeRaster(x = percabove_metrics, filename = glue::glue('{raster_output}/metrics/{acq}_percabove_{met_res}m.tif'), overwrite = TRUE)
    time_percabove <- tictoc::toc(log = TRUE)
    timings_df <- rbind(timings_df, data.frame(Metric = "Percent Above", Time = time_percabove$toc - time_percabove$tic, stringsAsFactors = FALSE))


    # Dispersion metrics
    tictoc::tic()
    print(glue::glue('Generating dispersion metrics for {acq} at {met_res}m'))
    dispersion_metrics <- pixel_metrics(ctg_norm, ~lidRmetrics::metrics_dispersion(Z), res = met_res)
    terra::writeRaster(x = dispersion_metrics, filename = glue::glue('{raster_output}/metrics/{acq}_dispersion_{met_res}m.tif'), overwrite = TRUE)
    time_dispersion <- tictoc::toc(log = TRUE)
    timings_df <- rbind(timings_df, data.frame(Metric = "Dispersion", Time = time_dispersion$toc - time_dispersion$tic, stringsAsFactors = FALSE))


    # Canopy density metrics
    tictoc::tic()
    print(glue::glue('Generating canopy density metrics for {acq} at {met_res}m'))
    canopydensity_metrics <- pixel_metrics(ctg_norm, ~lidRmetrics::metrics_canopydensity(Z), res = met_res)
    terra::writeRaster(x = canopydensity_metrics, filename = glue::glue('{raster_output}/metrics/{acq}_canopydensity_{met_res}m.tif'), overwrite = TRUE)
    time_canopydensity <- tictoc::toc(log = TRUE)
    timings_df <- rbind(timings_df, data.frame(Metric = "Canopy Density", Time = time_canopydensity$toc - time_canopydensity$tic, stringsAsFactors = FALSE))


    # L-moments metrics
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



    # Leaf Area Density metrics
    tictoc::tic()
    print(glue::glue('Generating leaf area density metrics for {acq} at {met_res}m'))
    lad_metrics <- pixel_metrics(ctg_norm, ~lidRmetrics::metrics_lad(Z), res = met_res)
    terra::writeRaster(x = lad_metrics, filename = glue::glue('{raster_output}/metrics/{acq}_lad_{met_res}m.tif'), overwrite = TRUE)
    time_lad <- tictoc::toc(log = TRUE)
    timings_df <- rbind(timings_df, data.frame(Metric = "Leaf Area Density", Time = time_lad$toc - time_lad$tic, stringsAsFactors = FALSE))


    # Interval metrics
    tictoc::tic()
    print(glue::glue('Generating interval metrics for {acq} at {met_res}m'))
    interval_metrics <- pixel_metrics(ctg_norm, ~lidRmetrics::metrics_interval(Z), res = met_res)
    terra::writeRaster(x = interval_metrics, filename = glue::glue('{raster_output}/metrics/{acq}_interval_{met_res}m.tif'), overwrite = TRUE)
    time_interval <- tictoc::toc(log = TRUE)
    timings_df <- rbind(timings_df, data.frame(Metric = "Intervals", Time = time_interval$toc - time_interval$tic, stringsAsFactors = FALSE))


    # Rumple metrics
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


    # Voxel metrics
    if(make_vox_mets == T){
    tictoc::tic()
    if(!exists('vox_size')){
      vox_size <- 1
    }
    print(glue::glue('Generating voxel metrics for {acq} at {met_res}m with {vox_size}m voxels'))
    voxel_metrics <- pixel_metrics(ctg_norm, ~lidRmetrics::metrics_voxels(x= X, y = Y, z = Z, vox_size = vox_size), res = met_res)
    terra::writeRaster(x = voxel_metrics, filename = glue::glue('{raster_output}/metrics/{acq}_voxels_{met_res}m.tif'), overwrite = TRUE)
    time_voxel <- tictoc::toc(log = TRUE)
    timings_df <- rbind(timings_df, data.frame(Metric = "Voxel", Time = time_voxel$toc - time_voxel$tic, stringsAsFactors = FALSE))
    } else{
      print('make_vox_mets FALSE; skipping time consuming voxel metric generation')
    }

    # KDE metrics
    tictoc::tic()
    print(glue::glue('Generating kernel density estimation (KDE) metrics for {acq} at {met_res}m'))
    kde_metrics <- pixel_metrics(ctg_norm, ~lidRmetrics::metrics_kde(Z), res = met_res)
    terra::writeRaster(x = kde_metrics, filename = glue::glue('{raster_output}/metrics/{acq}_kde_{met_res}m.tif'), overwrite = TRUE)
    time_KDE <- tictoc::toc(log = TRUE)
    timings_df <- rbind(timings_df, data.frame(Metric = "KDE", Time = time_KDE$toc - time_KDE$tic, stringsAsFactors = FALSE))

    # Visualize time taken

    # Order the data frame by time taken
    timings_df <- timings_df[order(timings_df$Time),]

    # Reorder the levels of the Metric factor based on the ordering of the data frame
    timings_df$Metric <- factor(timings_df$Metric, levels = timings_df$Metric)

    # Create a color scale based on time taken
    color_scale <- scale_fill_gradientn(colors = c("steelblue", "darkorange", "firebrick"))

    # Plot the timings using ggplot2 with the requested changes
    p <- ggplot(timings_df, aes(x = Metric, y = Time, fill = Time)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = round(Time, 1)), vjust = -0.5) +
      labs(title = "Time Taken for Each Metric", x = "Metric", y = "Time (seconds)") +
      theme_light() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      color_scale

    print(glue::glue('Metric generation for {acq} complete...'))
    print(p)
}
