#' Generate Canopy Height Model (CHM)
#'
#' This function generates a Canopy Height Model (CHM) from normalized LAS files.
#'
#' @param proj_dir Project directory containing the normalized LAS files.
#' @param res CHM resolution in meters. Default is 1.
#' @param algorithm Algorithm to use for CHM generation. Default is 'p2r'. Options are 'p2r' and 'pitfree'
#' @param subcircle_p2r Subcircle argument for p2r CHM generation, default is 0.
#' @param subcircle_pitfree Subcircle argument for pitfree CHM generation, default is 0.
#' @param num_cores Number of cores for parallel processing. Default is 1.
#' @param chunk_buf Chunk buffer size in meters. Default is NULL (auto-set to 5% of tile_size).
#' @param acq Acquisition name. Default is NULL.
#' @param clamp Logical switch to clamp CHM values. Default is FALSE.
#' @param min_clamp Minimum value to clamp CHM at. Default is 0.
#' @param max_clamp Maximum value to clamp CHM at. Default is 50.
#'
#' @return Saves a CHM locally as well as smoothed and filled versions of the CHM.
#'
#' @import lidR
#' @import terra
#' @import dplyr
#' @import glue
#' @importFrom future plan
#'
#' @export
generate_chm <- function(proj_dir, res = 1, algorithm = c('p2r','pitfree'),
                         subcircle_p2r = 0, subcircle_pitfree = 0,
                         num_cores = 1L, chunk_buf = NULL,
                         acq = NULL, clamp = FALSE,
                         min_clamp = 0, max_clamp = 50) {
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
  lidR::opt_chunk_buffer(ctg_norm) <- chunk_buf

  # Create the CHM output directory if it doesn't exist
  chm_output_dir <- glue::glue("{proj_dir}/output/raster/chm")
  if (!dir.exists(chm_output_dir)) {
    dir.create(chm_output_dir, showWarnings = FALSE, recursive = TRUE)
    print(glue::glue("Created a directory for CHM outputs {chm_output_dir}"))
  }

  # Delete intermediate raster tiles
  if (dir.exists(glue::glue("{chm_output_dir}/tiles"))) {
    unlink(glue::glue("{chm_output_dir}/tiles"), recursive = TRUE)
    print("Deleted intermediate CHM tiles")
  }

  # Set overwrite options
  ctg_norm@output_options$drivers$SpatRaster$param$overwrite <- TRUE

  # Generate the canopy height model
  if('p2r' %in% algorithm) {
    tictoc::tic()
    # Set output file options for saving CHM tiles
    lidR::opt_output_files(ctg_norm) <- glue::glue("{chm_output_dir}/tiles/{acq}_chm_{res}m_p2r_sc{subcircle_p2r}_{{XLEFT}}_{{YBOTTOM}}")
    # Generate CHM
    print(glue::glue('Generating {res}m CHM using p2r algorithm with subcircle {subcircle_p2r} for {acq}'))
    lidR::rasterize_canopy(ctg_norm, res = res,
                           algorithm = lidR::p2r(subcircle = subcircle_p2r,
                           na.fill = lidR::knnidw()))
    # Load CHM Tiles as a virtual raster dataset
    chm_tiles_dir <- glue::glue("{chm_output_dir}/tiles")
    chm_tiles <- list.files(chm_tiles_dir, pattern = '.tif$', full.names = TRUE)
    chm <- terra::vrt(chm_tiles, glue::glue("{chm_tiles_dir}/{acq}_chm.vrt"), overwrite = TRUE)
    if (clamp) {
      # Clamp CHM values within specified range
      chm <- terra::clamp(chm, min_clamp, max_clamp, values = TRUE)
    }
    # Set layer name to Z
    names(chm) <- 'Z'
    # Fill CHM - fill NA values with mean of 3x3 neighborhood
    chm_filled <- terra::focal(chm, w = 3, fun = "mean", na.policy = "only", na.rm = TRUE)
    names(chm_filled) <- 'Z'
    # Smooth CHM
    fgauss <- function(sigma, n = ws) {
      m <- matrix(ncol = n, nrow = n)
      col <- rep(1:n, n)
      row <- rep(1:n, each = n)
      x <- col - ceiling(n/2)
      y <- row - ceiling(n/2)
      m[cbind(row, col)] <- 1/(2 * pi * sigma^2) * exp(-(x^2 + y^2)/(2 * sigma^2))
      m/sum(m)
    }
    chm_smooth <- terra::focal(chm_filled, w = fgauss(1, n = 5))
    names(chm_smooth) <- 'Z'
    # Write CHM files
    terra::writeRaster(chm, filename = glue::glue("{chm_output_dir}/{acq}_chm_{res}m_p2r_sc{subcircle_p2r}.tif"), overwrite = TRUE)
    terra::writeRaster(chm_filled, filename = glue::glue("{chm_output_dir}/{acq}_chm_fill_p2r_sc{subcircle_p2r}_{res}m.tif"), overwrite = TRUE)
    terra::writeRaster(chm_smooth, filename = glue::glue("{chm_output_dir}/{acq}_chm_smooth_p2r_sc{subcircle_p2r}_{res}m.tif"), overwrite = TRUE)
    # Delete intermediate raster tiles
    if (dir.exists(glue::glue("{chm_output_dir}/tiles"))) {
      unlink(glue::glue("{chm_output_dir}/tiles"), recursive = TRUE)
      print("Deleted intermediate CHM tiles")
    }
  print(glue::glue('Completed generation of {res}m CHM using p2r algorithm with subcircle {subcircle_p2r} for {acq}'))
  tictoc::toc()
  }

  if('pitfree' %in% algorithm){
    tictoc::tic()
    # Set output file options for saving CHM tiles
    lidR::opt_output_files(ctg_norm) <- glue::glue("{chm_output_dir}/tiles/{acq}_chm_{res}m_pitfree_sc{subcircle_pitfree}_{{XLEFT}}_{{YBOTTOM}}")
    # Generate CHM
    print(glue::glue('Generating {res}m CHM using pitfree algorithm with subcircle {subcircle_pitfree} for {acq}'))
    lidR::rasterize_canopy(ctg_norm, res = res,
                           algorithm = lidR::pitfree(),
                           subcircle = subcircle_pitfree)
    # Load CHM Tiles as a virtual raster dataset
    chm_tiles_dir <- glue::glue("{chm_output_dir}/tiles")
    chm_tiles <- list.files(chm_tiles_dir, pattern = '.tif$', full.names = TRUE)
    chm <- terra::vrt(chm_tiles, glue::glue("{chm_tiles_dir}/{acq}_chm.vrt"), overwrite = TRUE)
    if (clamp) {
      # Clamp CHM values within specified range
      chm <- terra::clamp(chm, min_clamp, max_clamp, values = TRUE)
    }
    # Set layer name to Z
    names(chm) <- 'Z'
    # Fill CHM - fill NA values with mean of 3x3 neighborhood
    chm_filled <- terra::focal(chm, w = 3, fun = "mean", na.policy = "only", na.rm = TRUE)
    names(chm_filled) <- 'Z'
    # Smooth CHM
    fgauss <- function(sigma, n = ws) {
      m <- matrix(ncol = n, nrow = n)
      col <- rep(1:n, n)
      row <- rep(1:n, each = n)
      x <- col - ceiling(n/2)
      y <- row - ceiling(n/2)
      m[cbind(row, col)] <- 1/(2 * pi * sigma^2) * exp(-(x^2 + y^2)/(2 * sigma^2))
      m/sum(m)
    }
    chm_smooth <- terra::focal(chm_filled, w = fgauss(1, n = 5))
    names(chm_smooth) <- 'Z'
    # Write CHM files
    terra::writeRaster(chm, filename = glue::glue("{chm_output_dir}/{acq}_chm_{res}m_pitfree_sc{subcircle_pitfree}.tif"), overwrite = TRUE)
    terra::writeRaster(chm_filled, filename = glue::glue("{chm_output_dir}/{acq}_chm_{res}m_fill_pitfree_sc{subcircle_pitfree}.tif"), overwrite = TRUE)
    terra::writeRaster(chm_smooth, filename = glue::glue("{chm_output_dir}/{acq}_chm_{res}m_smooth_pitfree_sc{subcircle_pitfree}.tif"), overwrite = TRUE)
    # Delete intermediate raster tiles
    if (dir.exists(glue::glue("{chm_output_dir}/tiles"))) {
      unlink(glue::glue("{chm_output_dir}/tiles"), recursive = TRUE)
      print("Deleted intermediate CHM tiles")
    }
    print(glue::glue('Completed generation of {res}m CHM using pitfree algorithm with subcircle {subcircle_pitfree} for {acq}'))
    tictoc::toc()
  }
}
