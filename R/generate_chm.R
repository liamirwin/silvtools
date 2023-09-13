#' Generate Canopy Height Model (CHM)
#'
#' This function generates a Canopy Height Model (CHM) from normalized LAS files.
#'
#' @param proj_dir Project directory containing the normalized LAS files.
#' @param chm_res CHM resolution in meters. Default is 1.
#' @param num_cores Number of cores for parallel processing. Default is 1.
#' @param chunk_buf Chunk buffer size in meters. Default is NULL (auto-set to 5% of tile_size).
#' @param acq Acquisition name. Default is NULL.
#' @param clamp Logical switch to clamp CHM values. Default is FALSE.
#' @param min_clamp Minimum value to clamp CHM at. Default is 0.
#' @param max_clamp Maximum value to clamp CHM at. Default is 50.
#'
#' @return A SpatRaster representing the CHM.
#'
#' @import lidR
#' @import terra
#' @import dplyr
#' @import glue
#' @importFrom future plan
#'
#' @export
generate_chm <- function(proj_dir, chm_res = 1, num_cores = 1L, chunk_buf = NULL, acq = NULL, clamp = FALSE, min_clamp = 0, max_clamp = 50) {
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
  chm_output_dir <- glue::glue("{proj_dir}/output/chm")
  if (!dir.exists(chm_output_dir)) {
    dir.create(chm_output_dir, showWarnings = FALSE, recursive = TRUE)
    print(glue::glue("Created a directory for CHM {chm_output_dir}"))
  }

  # Set output files and overwrite options
  lidR::opt_output_files(ctg_norm) <- glue::glue("{chm_output_dir}/tiles/{acq}_chm_{chm_res}_{{XLEFT}}_{{YBOTTOM}}")
  ctg_norm@output_options$drivers$SpatRaster$param$overwrite <- TRUE

  # Generate the canopy height model
  lidR::rasterize_canopy(ctg_norm, res = chm_res, algorithm = lidR::p2r(na.fill = lidR::knnidw()))

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

  # Fill CHM
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
  terra::writeRaster(chm, filename = glue::glue("{chm_output_dir}/{acq}_chm_{chm_res}m_sub0_p2r.tif"), overwrite = TRUE)
  terra::writeRaster(chm_filled, filename = glue::glue("{chm_output_dir}/{acq}_chm_fill_p2r_{chm_res}m.tif"), overwrite = TRUE)
  terra::writeRaster(chm_smooth, filename = glue::glue("{chm_output_dir}/{acq}_chm_smooth_p2r_{chm_res}m.tif"), overwrite = TRUE)

  # Delete intermediate raster tiles
  if (dir.exists(glue::glue("{chm_output_dir}/tiles"))) {
    unlink(glue::glue("{chm_output_dir}/tiles"), recursive = TRUE)
    print("Deleted intermediate CHM tiles")
  }

  print(glue::glue("CHM generation process complete for {acq}"))

  return(chm) # Return the CHM as a SpatRaster
}
