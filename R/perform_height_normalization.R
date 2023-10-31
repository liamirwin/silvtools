#' Perform Height Normalization on Classified LAS Files
#'
#' @param proj_dir Project directory containing the classified las files.
#' @param chunk_buf Buffer area around each tile in meters.
#' @param tile_size The size of each tile for processing. Default is NULL, which will use the size of the middle tile in the catalog.
#' @param acq Acquisition name. Default is NULL.
#' @param norm_algorithm Normalization algorithm to use. Default is tin().
#' @param num_cores Number of cores for parallel processing. Default is 1.
#' @param index_tiles Logical, should tiles be indexed? Default is TRUE.
#' @param output_laz Logical, should the output be in LAZ format? Default is TRUE.
#' @param k Number of neighbours for knnidw algorithm
#' @param norm_dir Directory to write normalized tiles to; default is '{proj_dir}/input/las/norm'
#'
#' @return A catalog of height-normalized LAS files.
#' @export
#'
perform_height_normalization <- function(proj_dir,
                                         chunk_buf = NULL,
                                         tile_size = NULL,
                                         acq = NULL,
                                         norm_algorithm = 'tin',
                                         num_cores = 1L,
                                         index_tiles = TRUE,
                                         output_laz = TRUE,
                                         k = 8,
                                         norm_dir = NULL) {
  tictoc::tic()

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
  class_dir <- glue::glue("{proj_dir}/input/las/class")
  ctg_class <- lidR::catalog(class_dir)


  # If tile_size is NULL, use the size of the middle tile in the catalog
  if (is.null(tile_size)) {
    tile_size <- get_tile_size(ctg_class)
  }

  # If chunk_buf is NULL, set it to 5% of tile_size
  if (is.null(chunk_buf)) {
    chunk_buf <- tile_size * 0.05
  }

  # Set options
  lidR::opt_progress(ctg_class) <- TRUE
  lidR::opt_chunk_buffer(ctg_class) <- chunk_buf
  lidR::opt_laz_compression(ctg_class) <- TRUE

  if(is.null(norm_dir)){
  norm_dir <- glue::glue("{proj_dir}/input/las/norm")
  }

  # Check if norm_dir exists, if not create it
  if (!dir.exists(norm_dir)) {
    dir.create(norm_dir, recursive = TRUE)
  }

  # Set Normalization Algorithm
  if(norm_algorithm == 'tin') {
    algo <- lidR::tin()
  } else if (norm_algorithm == 'knnidw'){
    algo <- lidR::knnidw(k = k)
  }


  if (output_laz) {
    lidR::opt_output_files(ctg_class) <- glue::glue("{norm_dir}/{acq}_{{XLEFT}}_{{YBOTTOM}}_norm_{norm_algorithm}")
  } else {
    lidR::opt_output_files(ctg_class) <- glue::glue("{norm_dir}/{acq}_{{XLEFT}}_{{YBOTTOM}}_norm_{norm_algorithm}")
  }

  print(glue::glue('Beginning normalization of {acq} lidar tiles'))

  # Normalize point cloud using specified algorithm
  ctg_norm <- lidR::normalize_height(ctg_class, algorithm = algo)

  # Optionally index the normalized tiles
  if (index_tiles) {
  print(glue::glue('Indexing normalized lidar tiles for {acq}'))
    lidR:::catalog_laxindex(ctg_norm)
  }

  tictoc::toc()

}
