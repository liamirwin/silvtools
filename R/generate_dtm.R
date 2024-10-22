#' Generate Digital Terrain Model (DTM)
#'
#' This function generates a Digital Terrain Model (DTM) from classified LAS files.
#'
#' @param proj_dir Project directory containing the classified LAS files.
#' @param res DTM resolution in meters.
#' @param dtm_algorithm Custom DTM algorithm. Default is NULL (TIN). c('tin', 'knnidw', 'kriging')
#' @param num_cores Number of cores for parallel processing. Default is 1.
#' @param chunk_buf Chunk buffer size in meters. Default is NULL (auto-set to 5% of tile_size).
#' @param k Number of neighbours of IDW/KNN. Default is 8
#' @param acq Acquisition name. Default is NULL.
#'
#' @return Saves DTMs locally to proj_dir/output/raster/dtm
#'
#' @import lidR
#' @import terra
#' @import dplyr
#' @import glue
#' @importFrom future plan
#'
#' @export
generate_dtm <- function(proj_dir, res = 1, dtm_algorithm = 'tin',
                         num_cores = 1L, chunk_buf = NULL, k = 8,
                         acq = basename(proj_dir)) {
  tictoc::tic()
  # Handle parallelization
  if (num_cores == 1L) {
    future::plan("sequential")
  } else {
    future::plan("multisession", workers = num_cores)
  }

  # Set the catalog for classified LAS files
  ctg_class <- lidR::catalog(glue::glue("{proj_dir}/input/las/class"))

  # Set options for the catalog
  lidR::opt_progress(ctg_class) <- TRUE

  # Filter out non-ground points
  lidR::opt_filter(ctg_class) <- 'keep_class 2'

  print('filtered out non-ground points')

  # Auto-set the chunk buffer size (5% of tile_size) if not provided
  if (is.null(chunk_buf)) {
    tile_size <- get_tile_size(ctg_class)
    chunk_buf <- tile_size * 0.05
  }
  lidR::opt_chunk_buffer(ctg_class) <- chunk_buf

  opt_stop_early(ctg_class) <- FALSE

  # Create the DTM output directory if it doesn't exist
  dtm_output_dir <- glue::glue("{proj_dir}/output/raster/dtm")
  if (!dir.exists(dtm_output_dir)) {
    dir.create(dtm_output_dir, showWarnings = FALSE, recursive = TRUE)
    print(glue::glue("Created a directory for DTM {dtm_output_dir}"))
  }

  # Set output files and overwrite options
  opt_output_files(ctg_class) <- glue::glue("{dtm_output_dir}/tiles/{acq}_dtm_{res}_{{XLEFT}}_{{YBOTTOM}}")
  ctg_class@output_options$drivers$SpatRaster$param$overwrite <- TRUE


  print(glue::glue("Generating DTM for {acq} at {res}m using {dtm_algorithm} algorithm(s)"))

  # Choose the DTM algorithm based on user input or default to TIN

  for (algo in dtm_algorithm) {
    algo_name <- algo

    # Select the DTM algorithm based on user input
    if (algo_name == 'tin') {
      algo_fun <- lidR::tin()
    } else if (algo_name == 'knnidw') {
      algo_fun <- lidR::knnidw(k = k)
    } else if (algo_name == 'kriging') {
      algo_fun <- lidR::kriging(k = k)
    }

    dtm_tiles_dir <- glue::glue("{dtm_output_dir}/tiles")

    # Delete intermediate raster tiles
    if (dir.exists(dtm_tiles_dir)) {
      unlink(dtm_tiles_dir, recursive = TRUE)
      print("Deleted intermediate DTM tiles")
    }


  # Generate the terrain model
  lidR::rasterize_terrain(ctg_class, res = res, algorithm = algo_fun)


  # Load DTM Tiles as a virtual raster dataset
  dtm_tiles_dir <- glue::glue("{dtm_output_dir}/tiles")
  dtm_tiles <- list.files(dtm_tiles_dir, pattern = '.tif$', full.names = TRUE)
  dtm <- terra::vrt(dtm_tiles, glue::glue("{dtm_tiles_dir}/{acq}_dtm.vrt"), overwrite = TRUE)

  # Adjust filenames to reflect the algorithm used
  base_filename <- glue::glue("{dtm_output_dir}/{acq}_dtm_{algo_name}_{res}m")
  terra::writeRaster(dtm, filename = paste0(base_filename, ".tif"), overwrite = TRUE)

  # Smooth the DTM
  dtm_smooth <- dtm %>%
    terra::focal(w = matrix(1, 25, 25),
                 fun = mean,
                 na.rm = TRUE,
                 pad = TRUE)

  # Write the smoothed DTM to a raster file
  terra::writeRaster(dtm_smooth, filename = paste0(base_filename, "_smooth.tif"), overwrite = TRUE)

  # Delete intermediate raster tiles
  if (dir.exists(dtm_tiles_dir)) {
    unlink(dtm_tiles_dir, recursive = TRUE)
    print("Deleted intermediate DTM tiles")
  }

  print(glue::glue("DTM generation process complete for {acq} using {algo_name} algorithm"))
  print(glue::glue('Wrote DTM rasters to {dtm_output_dir}'))

  }

  tictoc::toc()
}
