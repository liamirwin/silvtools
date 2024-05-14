#' Applies tree IDs to LAS file based on crown intersections
#'
#' This function reads in a LAS file chunk and intersects it with a set of crowns
#' to assign tree IDs to each point in the LAS file that belongs to a tree. The
#' resulting LAS file will have an additional attribute called "treeID" that
#' contains the ID of the tree to which each point belongs.
#'
#' @param chunk A LAS file chunk to be processed.
#' @param crowns A set of crowns represented as an sf object.
#' @param filter Logical indicating whether to filter out points that do not correspond to a tree
#' @param buf Buffer distance to use when selecting crowns that intersect with the chunk
#'
#' @return A LAS file with an additional attribute called "treeID" that contains
#' the ID of the tree to which each point belongs.
#'
#' @import lidR
#' @import sf
#' @importFrom glue glue
#' @importFrom tictoc tic toc
#' @export
#'
#' @examples
#'
#' lasfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
#' las <- readLAS(lasfile)
#' norm_las <- normalize_height(las)
#'
#'
#'
#'
apply_treeid <- function(chunk, crowns, filter = TRUE, buf = 5){

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
