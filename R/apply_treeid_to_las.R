#' Applies tree IDs to LAS file based on crown intersections
#'
#' This function reads in a LAS file chunk and intersects it with a set of crowns
#' to assign tree IDs to each point in the LAS file that belongs to a tree. The
#' resulting LAS file will have an additional attribute called "treeID" that
#' contains the ID of the tree to which each point belongs.
#'
#' @param chunk A LAS file chunk to be processed.
#' @param crowns A set of crowns represented as an sf object.
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
apply_treeid_to_las <- function(chunk, crowns){
  tictoc::tic()
  las <- readLAS(chunk)
  box <- st_as_sfc(st_bbox(chunk))
  if (is.empty(las)) return(NULL)
  # if (!'treeID' %in% names(las@data)) return(NULL)
  box_buf <- st_buffer(box, dist = 5)
  # Select crowns relevant to chunk of interest
  chunk_crowns <- crowns[sf::st_intersects(crowns, box_buf, sparse = FALSE),]
  if(nrow(chunk_crowns) == 0) return(NULL)
  if("Z" %in% names(chunk_crowns)){
  names(chunk_crowns) <- c('treeID','geometry')
  }
  # Merge tree IDs with las
  tree_las <- merge_spatial(las, chunk_crowns, attribute = 'treeID')
  tree_las = add_lasattribute(tree_las, name="treeID", desc="ID of a tree")
  glue::glue('Merged {nrow(chunk_crowns)} crowns with {length(las@data$Z)} points in chunk')
  tictoc::toc()
  return(tree_las)
}
library(future)

# List directories (each is one acquisiton of ULS/DAP)
blocks_dir <- list.dirs('H:/Quesnel_2022/process', recursive = FALSE)
# Omit these already processed blocks from processing stream
processed <- c('CT2','CT3','CT4','CT5','CT1-T-DAP')
blocks_dir <- blocks_dir[!basename(blocks_dir) %in% processed]
i = 1
is_dap <- FALSE

proj_dir <- blocks_dir[i]

if(stringr::str_detect(basename(proj_dir), pattern = 'DAP')){
  is_dap == TRUE
  # Set acquisition name (DAPYY_blockname)
  acq <- glue::glue('DAP22_{stringr::str_replace(basename(proj_dir), pattern = "-DAP", replacement = "")}')
  print(glue::glue('Set acqusition type as DAP named {acq}'))
} else{
  is_dap == FALSE
  # Set acquisition name (ULSYY_blockname)
  acq <- glue::glue('ULS22_{basename(proj_dir)}')
  print(glue::glue('Set acqusition type as lidar (ULS) named {acq}'))
}


ctg <- catalog(glue::glue('{proj_dir}/input/las/norm'))
crowns <- st_read(glue::glue('{proj_dir}/output/vector/crowns/{acq}_lmf_ws2_watershed_crowns.shp'))
if(!dir.exists(glue::glue('{proj_dir}/output/tree_las'))){
  dir.create(glue::glue('{proj_dir}/output/tree_las'))
  glue::glue('Created tree las directory at {proj_dir}/output/tree_las')
}

opt_output_files(ctg) <- "{proj_dir}/output/tree_las/{*}_treeid"
opt_stop_early(ctg) <- F
opt_progress(ctg) <- T
opt_laz_compression(ctg) <- T
opt_chunk_buffer(ctg) <- 25
plan(multisession, workers = 3L)
tree_ctg <- catalog_apply(ctg, apply_treeid_to_las, crowns = crowns)

