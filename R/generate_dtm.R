#' Generate Digital Terrain Model (DTM)
#'
#' This function generates a Digital Terrain Model (DTM) from classified LAS files.
#'
#' @param proj_dir Project directory containing the classified LAS files.
#' @param dtm_res DTM resolution in meters.
#' @param num_cores Number of cores for parallel processing. Default is 1.
#' @param chunk_buf Chunk buffer size in meters. Default is NULL (auto-set to 5% of tile_size).
#' @param dtm_algorithm Custom DTM algorithm. Default is NULL (TIN).
#' @param acq Acquisition name. Default is NULL.
#'
#' @return A SpatRaster representing the DTM.
#'
#' @import lidR
#' @import terra
#' @import dplyr
#' @import glue
#' @importFrom future plan
#'
#' @export
generate_dtm <- function(proj_dir, dtm_res, num_cores = 1L, chunk_buf = NULL, dtm_algorithm = NULL, acq = NULL) {
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

  # Set the catalog for classified LAS files
  ctg_class <- lidR::catalog(glue::glue("{proj_dir}/input/las/class"))

  # Set options for the catalog
  lidR::opt_progress(ctg_class) <- TRUE

  # Auto-set the chunk buffer size (5% of tile_size) if not provided
  if (is.null(chunk_buf)) {
    tile_size <- get_tile_size(ctg_class)
    chunk_buf <- tile_size * 0.05
  }
  lidR::opt_chunk_buffer(ctg_class) <- chunk_buf

  opt_stop_early(ctg_class) <- FALSE

  # Create the DTM output directory if it doesn't exist
  dtm_output_dir <- glue::glue("{proj_dir}/output/dtm")
  if (!dir.exists(dtm_output_dir)) {
    dir.create(dtm_output_dir, showWarnings = FALSE, recursive = TRUE)
    print(glue::glue("Created a directory for DTM {dtm_output_dir}"))
  }

  # Set output files and overwrite options
  opt_output_files(ctg_class) <- glue::glue("{dtm_output_dir}/tiles/{acq}_dtm_{dtm_res}_{{XLEFT}}_{{YBOTTOM}}")
  ctg_class@output_options$drivers$SpatRaster$param$overwrite <- TRUE

  # Choose the DTM algorithm based on user input or default to TIN
  if (is.null(dtm_algorithm)) {
    dtm_algorithm <- lidR::tin()
  }

  # Generate the terrain model
  lidR::rasterize_terrain(ctg_class, res = dtm_res, algorithm = dtm_algorithm)

  # Load DTM Tiles as a virtual raster dataset
  dtm_tiles_dir <- glue::glue("{dtm_output_dir}/tiles")
  dtm_tiles <- list.files(dtm_tiles_dir, pattern = '.tif$', full.names = TRUE)
  dtm <- terra::vrt(dtm_tiles, glue::glue("{dtm_tiles_dir}/{acq}_dtm.vrt"), overwrite = TRUE)

  # Write the DTM to a raster file
  terra::writeRaster(dtm, filename = glue::glue("{dtm_output_dir}/{acq}_dtm_tin_{dtm_res}m.tif"), overwrite = TRUE)

  # Smooth the DTM
  dtm_smooth <- dtm %>%
    terra::focal(w = matrix(1, 25, 25),
                 fun = mean,
                 na.rm = TRUE,
                 pad = TRUE)

  # Write the smoothed DTM to a raster file
  terra::writeRaster(dtm_smooth, filename = glue::glue("{dtm_output_dir}/{acq}_dtm_tin_smooth_{dtm_res}m.tif"), overwrite = TRUE)

  print(glue::glue("DTM generation process complete for {acq}"))

  return(dtm) # Return the DTM as a SpatRaster
}
