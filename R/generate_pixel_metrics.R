#' Generate lidar pixel metrics
#'
#' This function generates a pixel metrics from tiled LAS files.
#'
#' @param proj_dir Project directory containing the normalized LAS files.
#' @param res CHM resolution in meters. Default is 1.
#' @param metrics Sets of metrics to generate. Default is 'basic'. 'percentiles' and 'scan angle' are also available.
#' @param zmin Minimum Z value of points to consider in metrics. Default NA but EFI's often use 2m
#' @param num_cores Number of cores to use for parallel processing. Default is 1.
#' @param chunk_buf Chunk buffer size in meters. Default is 5% of tile size.
#' @param acq Acquisition name. Default is NULL.
#' @param ctg_norm Catalog of normalized LAS files. Default is NULL which will use the normalized LAS files in the proj_dir.
#' @param mets_output_dir Directory to save metrics. Default is NULL which will save to proj_dir.
#' @param abs Boolean to determine if absolute values should be used for Scan Angle metrics. Default is TRUE.
#' @param SA_attribute Scan angle attribute name. Default is 'ScanAngle'.
#' @param is_rads Boolean to determine if Scan Angle values are in radians. Default is FALSE.
#'
#' @return Saves lidar pixel metrics locally within proj_dir
#'
#' @import lidR
#' @import terra
#' @import dplyr
#' @import glue
#' @importFrom future plan
#'
#' @export
generate_pixel_metrics <- function(proj_dir, res = 1, metrics = c('basic'), zmin = NULL,
                         num_cores = 1L, chunk_buf = NULL,
                         acq = NULL, ctg_norm = NULL, mets_output_dir = NULL, abs = TRUE,
                         SA_attribute = 'ScanAngle', is_rads = FALSE) {
  # Handle parallelization
  if (num_cores == 1L) {
    future::plan("sequential")
  } else {
    future::plan("multisession", workers = num_cores)
  }

  # Define acquisition name dynamically if not supplied
  if (is.null(acq)) {
    acq <- basename(proj_dir)
  }

  # Set options for the catalog
  if(is.null(ctg_norm)) {
    ctg_norm <- lidR::catalog(glue::glue("{proj_dir}/input/las/norm"))
  } else {
    ctg_norm <- lidR::catalog(ctg_norm)
  }

  lidR::opt_progress(ctg_norm) <- TRUE

  # Auto-set the chunk buffer size (5% of tile_size) if not provided
  if (is.null(chunk_buf)) {
    tile_size <- get_tile_size(ctg_norm)
    chunk_buf <- tile_size * 0.05
  }
  # Set processing chunk buffer size
  lidR::opt_chunk_buffer(ctg_norm) <- chunk_buf


  if(is.null(mets_output_dir)) {
    # Create the metrics output directory if it doesn't exist
    mets_output_dir <- glue::glue("{proj_dir}/output/raster/metrics")
    if (!dir.exists(mets_output_dir)) {
      dir.create(mets_output_dir, showWarnings = FALSE, recursive = TRUE)
      print(glue::glue("Created a directory for metrics outputs {mets_output_dir}"))
    }
  } else{
    if (!dir.exists(mets_output_dir)) {
      dir.create(mets_output_dir, showWarnings = FALSE, recursive = TRUE)
      print(glue::glue("Created a directory for metrics outputs {mets_output_dir}"))
    }
  }

  # Set overwrite options
  ctg_norm@output_options$drivers$SpatRaster$param$overwrite <- TRUE

  # If zmin is not NA print/save in output names
  if (!is.null(zmin)) {
    zmin_msg <- glue::glue("zmin{zmin}_")
  } else {
    zmin_msg <- ''
  }

  # Generate Standard Metrics
  if('basic' %in% metrics) {
    tictoc::tic()
    print(glue::glue("Generating basic pixel metrics for {acq} at {res}m resolution {zmin_msg}"))
    # Set output file options for saving CHM tiles
    lidR::opt_output_files(ctg_norm) <- glue::glue("{mets_output_dir}/tiles/{acq}_basic_{zmin_msg}{res}m_{{XLEFT}}_{{YBOTTOM}}")
    metrics <-  pixel_metrics(ctg_norm, func = ~lidRmetrics::metrics_basic(z = Z, zmin = zmin), res = res)
    # Load Metrics Tiles as a virtual raster dataset
    mets_tiles_dir <- glue::glue("{mets_output_dir}/tiles")
    mets_tiles <- list.files(mets_tiles_dir, pattern = '.tif$', full.names = TRUE)
    mets <- terra::vrt(mets_tiles, glue::glue("{mets_tiles_dir}/{acq}_basic.vrt"), overwrite = TRUE)
    names(mets) <- names(terra::rast(mets_tiles[1]))
    # Save Metrics to disk
    terra::writeRaster(mets, glue::glue("{mets_output_dir}/{acq}_basic_{zmin_msg}_{res}m.tif"), overwrite = TRUE)
    # Delete intermediate tiles
    file.remove(mets_tiles)
    print('Generated standard pixel metrics at {res}m resolution for {acq} {zmin_msg}')
    tictoc::toc()
  }

  if('percentiles' %in% metrics) {
    tictoc::tic()
    print(glue::glue("Generating percentile pixel metrics for {acq} at {res}m resolution {zmin_msg}"))
    # Set output file options for saving CHM tiles
    lidR::opt_output_files(ctg_norm) <- glue::glue("{mets_output_dir}/tiles/{acq}_percentiles_{zmin_msg}_{res}m_{{XLEFT}}_{{YBOTTOM}}")
    metrics <-  pixel_metrics(ctg_norm, func = ~lidRmetrics::metrics_percentiles(z = Z, zmin = zmin), res = res)
    # Load Metrics Tiles as a virtual raster dataset
    mets_tiles_dir <- glue::glue("{mets_output_dir}/tiles")
    mets_tiles <- list.files(mets_tiles_dir, pattern = '.tif$', full.names = TRUE)
    mets <- terra::vrt(mets_tiles, glue::glue("{mets_tiles_dir}/{acq}_percentiles.vrt"), overwrite = TRUE)
    # Save Metrics to disk
    terra::writeRaster(mets, glue::glue("{mets_output_dir}/{acq}_percentiles_{zmin_msg}_{res}m.tif"), overwrite = TRUE)
    # Delete intermediate tiles
    file.remove(mets_tiles)
    print('Generated percentile pixel metrics at {res}m resolution for {acq} {zmin_msg}')
    tictoc::toc()
  }

  if('scan angle' %in% metrics){

    if(is_rads){
      print('is_rads = TRUE... Scan angle metrics will be converted from radians to degrees')
    }

    scan_angle_metrics <- function(ScanAngle, abs = TRUE, rads = FALSE){

      # If not integer, convert to integer
      if (is.integer(ScanAngle)) {
        ScanAngle <- as.numeric(ScanAngle)
      }

      # If radians, convert to degrees
      if (rads) {
        # Convert from radians to degrees
        ScanAngle <- ScanAngle * (180 / pi)
      }

      if(abs == TRUE){
        # Take absolute value of scan angle
        ScanAngle <- abs(ScanAngle)
      }

      list(
        mean <- mean(ScanAngle),
        median <- median(ScanAngle),
        sd <- sd(ScanAngle),
        cv <- sd(ScanAngle) / mean(ScanAngle),
        iqr <- IQR(ScanAngle),
        min <- min(ScanAngle),
        max <- max(ScanAngle),
        n <- length(ScanAngle)
      )

    }

    tictoc::tic()
    print(glue::glue("Generating scan angle metrics for {acq} at {res}m resolution {zmin_msg}"))
    # Set output file options for saving CHM tiles
    lidR::opt_output_files(ctg_norm) <- glue::glue("{mets_output_dir}/SA/tiles/{acq}_scan_angle_{zmin_msg}_{res}m_{{XLEFT}}_{{YBOTTOM}}")
    lidR::opt_select(ctg_norm) <- 'sa'
    # Create SA metrics and SA tile directories
    dir.create(glue::glue("{mets_output_dir}/SA"), showWarnings = FALSE, recursive = TRUE)
    dir.create(glue::glue("{mets_output_dir}/SA/tiles"), showWarnings = FALSE, recursive = TRUE)
    # Generate SA metrics
    metrics <-  pixel_metrics(ctg_norm, func = ~scan_angle_metrics(ScanAngle = get(SA_attribute), abs = abs, rads = is_rads), res = res)
    # Load Metrics Tiles as a virtual raster dataset
    sa_met_tiles <- list.files(glue::glue("{mets_output_dir}/SA/tiles"), pattern = '.tif$', full.names = TRUE)
    sa_met <- terra::vrt(sa_met_tiles, glue::glue("{mets_output_dir}/SA/tiles/{acq}_scan_angle.vrt"), overwrite = TRUE)
    names(sa_met) <- c('mean','median','sd','cv','iqr','min','max','n')
    # Save Metrics to disk
    terra::writeRaster(sa_met, glue::glue("{mets_output_dir}/{acq}_SA_{zmin_msg}_{res}m.tif"), overwrite = TRUE)
    # Delete intermediate tiles
    file.remove(sa_met_tiles)
    print(glue::glue('Generated scan angle metrics at {res}m resolution for {acq} {zmin_msg}'))
    tictoc::toc()
  }

}
