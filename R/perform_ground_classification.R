#' Perform Ground Classification on a LAScatalog
#'
#' Classify the ground points in a LAScatalog using a specified algorithm.
#' Automatically sets tile size and buffer if not supplied.
#'
#' @param proj_dir The project directory where the LAS files are stored.
#' @param tile_size The size of each tile for processing. Default is NULL, which will use the size of the middle tile in the catalog.
#' @param chunk_buf The buffer around each chunk. Default is NULL, set to 5% of tile_size.
#' @param acq The name of the acquisition, default is NULL which uses the basename of the project directory.
#' @param ground_algorithm The ground classification algorithm to use. Default is NULL which uses the CSF algorithm.
#' @param num_cores The number of cores for parallel processing. Default is 1.
#' @param index_tiles Boolean, whether to index the classified tiles. Default is TRUE.
#' @param output_laz Logical, should the output be in LAZ format? Default is TRUE.
#'
#' @return A LAScatalog containing the ground classified points.
#'
#' @export
#' @importFrom lidR classify_ground catalog csf readLASheader
#' @importFrom glue glue
#' @importFrom future plan
perform_ground_classification <- function(proj_dir, tile_size = NULL, chunk_buf = NULL, acq = NULL, ground_algorithm = NULL, num_cores = 1L, index_tiles = TRUE, output_laz = TRUE) {

  # Handle parallelization
  if (num_cores == 1L) {
    future::plan("sequential")
  } else {
    future::plan("multisession", workers = num_cores)
  }

  # Handle acquisition name
  if (is.null(acq)) {
    acq <- basename(proj_dir)
  }

  # Initialize catalog
  tile_dir <- glue::glue("{proj_dir}/input/las/tile")
  ctg_tile <- lidR::catalog(tile_dir)

  # If tile_size is NULL, use the size of the middle tile in the catalog
  if (is.null(tile_size)) {
    tile_size <- get_tile_size(ctg_tile)
  }

  # If chunk_buf is NULL, set it to 5% of tile_size
  if (is.null(chunk_buf)) {
    chunk_buf <- tile_size * 0.05
  }

  # Set options
  lidR::opt_progress(ctg_tile) <- TRUE
  lidR::opt_chunk_size(ctg_tile) <- tile_size
  lidR::opt_chunk_buffer(ctg_tile) <- chunk_buf

  class_dir <- glue::glue("{proj_dir}/input/las/class")
  lidR::opt_output_files(ctg_tile) <- glue::glue("{class_dir}/{acq}_{{XLEFT}}_{{YBOTTOM}}_class")

  # Set default ground classification algorithm if not provided
  if (is.null(ground_algorithm)) {
    ground_algorithm <- lidR::csf(class_threshold = 0.25, cloth_resolution = 0.25, rigidness = 2)
  }

  # Classify ground using the specified algorithm
  ctg_class <- lidR::classify_ground(ctg_tile, algorithm = ground_algorithm)

  # Optionally index the classified tiles
  if (index_tiles) {
    lidR:::catalog_laxindex(ctg_class)
  }

  # Set the output format (LAS or LAZ)
  if (output_laz) {
    ctg_class@output_options$drivers$LAS$param = list(index = index_tiles, compression = "LAZ")
  } else {
    ctg_class@output_options$drivers$LAS$param = list(index = index_tiles, compression = "NONE")
  }

  return(ctg_class)
}
