#' Add Crown IDs to LAS Files
#'
#' This function processes the LAS files in a catalog to add tree crown IDs based on the supplied crown vectors or rasters
#'
#' @param proj_dir Project directory containing the input LAS files and the output directories.
#' @param acq Acquisition name. If NULL (default), the basename of the `proj_dir` will be used as the acquisition name.
#' @param num_cores Number of cores to use for parallel processing. Default is 1.
#' @param chunk_buf Chunk buffer size. If NULL (default), it is auto-set to 5% of tile_size.
#' @param crown_dir Directory containing the crown shapefiles. Default is NULL.
#' @param crown_pattern Pattern to match the crown files. Default is 'lmfws2'.
#' @param crown_ext File extension for the crown files. Default is '.gpkg$'.
#' @param filter Logical to decide whether to apply any filtering while adding tree crowns. Default is FALSE.
#'
#' @return Creates new LAS files with tree crown IDs attached as attribute.
#'
#' @import lidR
#' @importFrom glue glue
#' @importFrom future plan
#'
#' @export
add_crowns_to_las <- function(proj_dir, acq = NULL,
                              num_cores = 1L, chunk_buf = NULL,
                              crown_dir = NULL, crown_pattern = 'lmfws2',
                              crown_ext = '.gpkg$', filter = FALSE) {
  tictoc::tic()
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
  ctg <- lidR::catalog(glue::glue("{proj_dir}/input/las/norm"))
  lidR::opt_progress(ctg) <- TRUE
  lidR::opt_output_files(ctg) <- "{proj_dir}/output/tree_las/{*}_treeid"
  lidR::opt_laz_compression(ctg) <- T

  # Create tree_las directory if it doesn't exist
  if (!dir.exists(glue::glue("{proj_dir}/output/tree_las"))) {
    dir.create(glue::glue("{proj_dir}/output/tree_las"))
  }

  # Auto-set the chunk buffer size (5% of tile_size) if not provided
  if (is.null(chunk_buf)) {
    tile_size <- get_tile_size(ctg)
    chunk_buf <- tile_size * 0.05
  }
  lidR::opt_chunk_buffer(ctg) <- chunk_buf

  # Load Crowns of Choice

  # If Crown Directory is not provided, use the default crown directory
  if (is.null(crown_dir)) {
    crown_dir <- glue::glue("{proj_dir}/output/vector/crowns")
  }

  crowns_file <- select_file_path(crown_dir, pattern = crown_pattern, ext = crown_ext)
  crowns <- sf::st_read(crowns_file, quiet = TRUE)


  apply_treeid_to_ctg <- function(chunk, crowns, filter = TRUE, buf = 5){

    # Read in chunk and check if empty
    las <- lidR::readLAS(chunk)
    box <- sf::st_as_sfc(st_bbox(chunk))
    if (lidR::is.empty(las)) return(NULL)

    box_buf <- sf::st_buffer(box, dist = buf)

    # Select crowns relevant to chunk of interest
    chunk_crowns <- crowns[sf::st_intersects(crowns, box_buf, sparse = FALSE),]
    # If no crowns intersect with chunk, return NULL
    if(nrow(chunk_crowns) == 0) return(NULL)
    # Rename columns if necessary
    if("Z" %in% names(chunk_crowns)){
      names(chunk_crowns) <- c('treeID','geometry')
    }
    # Merge tree IDs with las
    tree_las <- lidR::merge_spatial(las, chunk_crowns, attribute = 'treeID')
    # Add treeID attribute to las
    tree_las = add_lasattribute(tree_las, name="treeID", desc="ID of a tree")
    # Filter out points that do not belong to a tree
    if(filter){
      tree_las <- filter_poi(tree_las, !is.na(treeID))
    }
    glue::glue('Merged {nrow(chunk_crowns)} crowns with {length(las@data$Z)} points in chunk')
    return(tree_las)
  }

  # Apply treeIDs from polygons
  print(glue::glue('Starting to apply tree ids of {nrow(crowns)} to {nrow(ctg)} las tiles for {acq}'))
  tree_ctg <- lidR::catalog_apply(ctg, apply_treeid_to_ctg, crowns = crowns, filter = filter)
  print(glue::glue('Finished applying tree ids for {acq}'))

  # Index newly created LAS files
  lidR:::catalog_laxindex(tree_ctg)
  print(glue::glue('Finished indexing tree id las for {acq}'))
  tictoc::toc()

}
