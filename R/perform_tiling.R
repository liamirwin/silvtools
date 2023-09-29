#' Perform Tiling on LAS files
#'
#' @param proj_dir Project directory containing the raw las files.
#' @param tile_size Size of each tile in meters.
#' @param chunk_buffer Buffer area around each tile.
#' @param index_tiles Logical, should tiles be indexed?
#' @param output_laz Logical, should the output be in LAZ format?
#' @param num_cores Number of cores for parallel processing.
#' @param acq Acquisition name.
#'
#' @return A catalog of tiled LAS files.
#' @export
#'
perform_tiling <- function(proj_dir,
                           tile_size = 250,
                           chunk_buffer = 0,
                           index_tiles = TRUE,
                           output_laz = TRUE,
                           num_cores = 1L,
                           acq = NULL) {

  if (num_cores == 1L) {
    future::plan("sequential")
  } else {
    future::plan("multisession", workers = num_cores)
  }

  if (is.null(acq)) {
    acq <- basename(proj_dir)
  }

  if (is.null(proj_dir)) {
    stop("proj_dir must be specified.")
  }

  tile_dir <- glue::glue("{proj_dir}/input/las/tile")

  ctg_tile <- lidR::catalog(glue::glue('{proj_dir}/input/las/raw'))

  lidR::opt_chunk_size(ctg_tile) <- tile_size
  lidR::opt_chunk_buffer(ctg_tile) <- chunk_buffer
  lidR::opt_progress(ctg_tile) <- TRUE
  lidR::opt_output_files(ctg_tile) <- glue::glue("{tile_dir}/{acq}_{{XLEFT}}_{{YBOTTOM}}")



  if (output_laz) {
    ctg_tile@output_options$drivers$LAS$param = list(index = index_tiles, compression = "LAZ")
  } else {
    ctg_tile@output_options$drivers$LAS$param = list(index = index_tiles, compression = "NONE")
  }

  if (!dir.exists(tile_dir)) {
    dir.create(tile_dir, recursive = TRUE)
  }

  ctg_tile <- lidR::catalog_retile(ctg_tile)

  if (index_tiles) {
    lidR:::catalog_laxindex(ctg_tile)
  }

  return(ctg_tile)
}
