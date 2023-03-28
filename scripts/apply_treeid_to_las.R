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

library(glue)
library(lidR)
library(future.apply)
library(future)
library(sf)
library(silvtools)
library(dplyr)
# List directories (each is one acquisiton of ULS/DAP)
blocks_dir <- list.dirs('H:/Quesnel_2022/process', recursive = FALSE)
# Omit these already processed blocks from processing stream
processed <- c('CT1','CT5')
blocks_dir <- blocks_dir[!basename(blocks_dir) %in% processed]
blocks_dir <- c("H:/Quesnel_2022/process/CT1-T-DAP")
is_dap <- TRUE

for(i in 1:length(blocks_dir)){

proj_dir <- blocks_dir[i]

if(stringr::str_detect(basename(proj_dir), pattern = 'DAP')){
  is_dap = TRUE
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
opt_chunk_buffer(ctg) <- 5
opt_filter(ctg) <- "-thin_with_voxel 0.05"
plan(multisession, workers = 3L)
plan(sequential)
tree_ctg <- catalog_apply(ctg, apply_treeid_to_las, crowns = crowns)
print(glue::glue('Finished applying tree ids for {acq}, beginning alphashape computation at {Sys.time()}'))
generate_alphashape(proj_dir, acq, num_cores = 3L)
print(glue::glue('Finished applying tree ids for {acq}, beginning alphashape computation at {Sys.time()}'))
}

# Clean duplicate alphashapes by filtering same tree ID and taking tree with greater n_points

# List directories (each is one acquisiton of ULS/DAP)
blocks_dir <- list.dirs('H:/Quesnel_2022/process', recursive = FALSE)
# Omit these already processed blocks from processing stream
processed <- c('CT1','CT2','CT3','CT4','CT5')
blocks_dir <- blocks_dir[!basename(blocks_dir) %in% processed]
i = 1
is_dap <- TRUE

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



# Functionalize

generate_alphashape <- function(proj_dir, acq, num_cores = 1) {
  # Generate Alphashapes

  # Define output directories
  raster_output <- glue::glue('{proj_dir}/output/raster')
  vector_output <- glue::glue('{proj_dir}/output/vector')
  image_output <- glue::glue('{proj_dir}/output/png/ashape_snapshots')

  # Get a list of all .laz files in the input directory
  tree_las_list <- list.files(glue::glue('{proj_dir}/output/tree_las'), pattern = '.laz', full.names = T)

  # Initialize a list to store alphashape metrics
  ashape_mets <- list()

  # Check if the output directory for alphashapes exists, and create it if it doesn't
  if(!dir.exists(glue::glue('{proj_dir}/output/crowns/ashapes'))){
    dir.create(glue::glue('{proj_dir}/output/crowns/ashapes'), recursive = T)
    print(glue::glue('Created alphashape save directory for {acq}'))
  }

  # Define a function to get alphashape metrics for a single LAS file
  get_alphashape_metrics_fun <- function(las_filename) {
    las <- readLAS(las_filename, filter = "-thin_with_voxel 0.05")
    ashape_metrics <- get_alphashape_metrics(las)
    readr::write_csv(
      ashape_metrics,
      file = glue::glue('{proj_dir}/output/crowns/ashapes/{acq}_chunk_5cmvoxel_ashapes.csv'),
      append = T,
      col_names = T)
    return(ashape_metrics)
  }
  print(glue::glue('Beginning parallel processing for {acq} Alphashapes at {Sys.time()}'))
  # Use future_lapply to process the LAS files in parallel if more than one core is specified
  if(num_cores > 1){
    library(future.apply)
    plan(multisession, workers = num_cores)

    results <- future_lapply(tree_las_list, get_alphashape_metrics_fun)

  } else {
    # Otherwise, use lapply to process the LAS files sequentially
    results <- lapply(tree_las_list, get_alphashape_metrics_fun)
  }

  # Define a function to print progress updates
  print_progress <- function(i, num_files) {
    message(glue::glue('Processed {i} out of {num_files} files.'))
  }

  # Loop through the results and print progress updates
  for(i in seq_along(results)) {
    print_progress(i, length(results))
  }

  # Return NULL invisibly
  return(invisible(NULL))
}



generate_alphashape(proj_dir, acq, num_cores = 4L)


# NON PARRALEL FUNCTION

generate_alphashapes <- function(proj_dir, acq) {
  # Define paths to output directories
  raster_output <- glue::glue('{proj_dir}/output/raster')
  vector_output <- glue::glue('{proj_dir}/output/vector')
  image_output <- glue::glue('{proj_dir}/output/png/ashape_snapshots')

  # Get a list of input files
  tree_las_list <- list.files(glue::glue('{proj_dir}/input/tree_las'), pattern = '.laz', full.names = T)

  # Check if output directory for alphashapes exists, create it if it doesn't
  if(!dir.exists(glue::glue('{proj_dir}/output/crowns/ashapes'))){
    dir.create(glue::glue('{proj_dir}/output/crowns/ashapes'), recursive = T)
    print(glue::glue('Created alphashape save directory for {acq}'))
  }

  # Define function to generate alphashape metrics for a single file
  generate_alphashape <- function(las_file) {
    las <- readLAS(las_file, filter = "-thin_with_voxel 0.05")
    print(glue::glue('Loaded las {las_file}'))
    ashape_metrics <- get_alphashape_metrics(las)
    # Write alphashape
    readr::write_csv(ashape_metrics,
                     file = glue::glue('{proj_dir}/output/crowns/ashapes/{acq}_chunk_5cmvoxel_ashapes.csv'),
                     append = T)
    print(glue::glue('Done processing {las_file}'))
    return(ashape_metrics)
  }



  # Apply generate_alphashape function to each input file
  ashape_mets <- lapply(tree_las_list, generate_alphashape)



  return(ashape_mets)
}
