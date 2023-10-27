#' Generate lidar pixel metrics
#'
#' This function generates a pixel metrics from tiled LAS files.
#'
#' @param proj_dir Project directory containing the normalized LAS files.
#' @param res CHM resolution in meters. Default is 1.
#' @param metrics Sets of metrics to generate. Default is 'stdmetrics'.
#' @param num_cores Number of cores to use for parallel processing. Default is 1.
#' @param chunk_buf Chunk buffer size in meters. Default is 5% of tile size.
#' @param acq Acquisition name. Default is NULL.
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
generate_pixel_metrics <- function(proj_dir, res = 1, metrics = c('stdmetrics'),
                         num_cores = 1L, chunk_buf = NULL,
                         acq = NULL) {
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
  ctg_norm <- lidR::catalog(glue::glue("{proj_dir}/input/las/norm"))
  lidR::opt_progress(ctg_norm) <- TRUE

  # Auto-set the chunk buffer size (5% of tile_size) if not provided
  if (is.null(chunk_buf)) {
    tile_size <- get_tile_size(ctg_norm)
    chunk_buf <- tile_size * 0.05
  }
  # Set processing chunk buffer size
  lidR::opt_chunk_buffer(ctg_norm) <- chunk_buf

  # Create the metrics output directory if it doesn't exist
  mets_output_dir <- glue::glue("{proj_dir}/output/raster/metrics")
  if (!dir.exists(mets_output_dir)) {
    dir.create(mets_output_dir, showWarnings = FALSE, recursive = TRUE)
    print(glue::glue("Created a directory for metrics outputs {mets_output_dir}"))
  }

  # Set overwrite options
  ctg_norm@output_options$drivers$SpatRaster$param$overwrite <- TRUE

  # Generate Standard Metrics
  if('stdmetrics' %in% metrics) {
    tictoc::tic()
    print(glue::glue("Generating standard pixel metrics for {acq} at {res}m resolution"))
    # Set output file options for saving CHM tiles
    lidR::opt_output_files(ctg_norm) <- glue::glue("{mets_output_dir}/tiles/{acq}_stdmetrics_{res}m_{{XLEFT}}_{{YBOTTOM}}")
    metrics <-  pixel_metrics(ctg_norm, ~stdmetrics(X,Y,Z,Intensity,ReturnNumber,Classification,dz=1), res = met_res)
    # Load Metrics Tiles as a virtual raster dataset
    mets_tiles_dir <- glue::glue("{mets_output_dir}/tiles")
    mets_tiles <- list.files(mets_tiles_dir, pattern = '.tif$', full.names = TRUE)
    mets <- terra::vrt(mets_tiles, glue::glue("{mets_tiles_dir}/{acq}_stdmetrics.vrt"), overwrite = TRUE)
    # Save Metrics to disk
    terra::writeRaster(mets, glue::glue("{mets_output_dir}/{acq}_stdmetrics_{res}m.tif"), overwrite = TRUE)
    # Delete intermediate tiles
    file.remove(mets_tiles)
    print('Generated standard pixel metrics at {res}m resolution for {acq}')
    tictoc::toc()
  }

}
