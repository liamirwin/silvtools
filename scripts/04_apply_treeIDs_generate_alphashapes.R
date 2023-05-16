library(future)
library(glue)
library(lidR)
library(future.apply)
library(future)
library(sf)
library(silvtools)
library(dplyr)

# Apply segmented crowns to las point clouds

# List directories (each is one acquisiton of ULS/DAP)
blocks_dir <- list.dirs('G:/Quesnel_2022/process', recursive = FALSE)

# Omit these already processed blocks from processing stream
processed <- c('CT3','CT4')
processed <- c('CT2','CT4')
blocks_dir <- blocks_dir[basename(blocks_dir) %in% processed]
is_dap <- F
# Run in parallel?
run_parallel <- T
num_cores <- 2L
# Chunk buffer
chunk_buf <- 5

# ---- Apply TreeID from Crowns to LAS Files


for(i in 1:length(blocks_dir)){
  tictoc::tic()
  proj_dir <- blocks_dir[i]

  if(stringr::str_detect(basename(proj_dir), pattern = 'DAP') | is_dap == T){
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

  # Load Normalized LAS tiles
  ctg <- catalog(glue::glue('{proj_dir}/input/las/norm'))
  # Load pre-generated cleaned tree crown polygons
  crowns <- st_read(glue::glue('{proj_dir}/output/vector/crowns/{acq}_lmf_ws2_watershed_crowns.shp'))

  # Create directory for output treeID las files
  if(!dir.exists(glue::glue('{proj_dir}/output/tree_las'))){
    dir.create(glue::glue('{proj_dir}/output/tree_las'))
    glue::glue('Created tree las directory at {proj_dir}/output/tree_las')
  }

  # Catalog setup
  opt_output_files(ctg) <- "{proj_dir}/output/tree_las/{*}_treeid"
  opt_stop_early(ctg) <- F
  opt_progress(ctg) <- T
  opt_laz_compression(ctg) <- T
  opt_chunk_buffer(ctg) <- 0
  opt_filter(ctg) <- "-thin_with_voxel 0.05"

  if(run_parallel){
    # Enable parallel processing
    plan(multisession, workers = 3L)
    print(glue::glue('Parallel processing initiated for {acq} with {num_cores} cores'))
  } else{
    plan(sequential)
    print(glue::glue('Beginning applying treeIDs for {acq} with one core'))
  }
  # Apply treeIDs from polygons
  tree_ctg <- catalog_apply(ctg, apply_treeid_to_las, crowns = crowns)
  print(glue::glue('Finished applying tree ids for {acq}'))
  # Index newly created LAS files
  lidR:::catalog_laxindex(tree_ctg)
  print(glue::glue('Finished indexing tree id las for {acq}'))
  tictoc::toc()

}

# Generate alphashapes for treeID las files

for(i in 1:length(blocks_dir)){

  tictoc::tic()
  proj_dir <- blocks_dir[i]

  # Setup acquisiton name

  if(stringr::str_detect(basename(proj_dir), pattern = 'DAP') | is_dap == T){
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

  # Load treeID LAS files as catalog
  tree_ctg <- catalog(glue::glue('{proj_dir}/output/tree_las'))
  opt_progress(tree_ctg) <- T
  opt_chunk_buffer(tree_ctg) <- chunk_buf
  # Save alphashape dataframe for each tile with XLEFT and YBOTTOM position
  opt_output_files(tree_ctg) <- '{proj_dir}/output/crowns/{acq}_ashapes_{XLEFT}_{YBOTTOM}'
  # Output dataframes as CSV
  tree_ctg@output_options$drivers$data.frame$extension <- '.csv'

  if(run_parallel){
    # Enable parallel processing
    plan(multisession, workers = 3L)
    print(glue::glue('Parallel processing initiated for {acq} with {num_cores} cores'))
  } else{
    plan(sequential)
    print(glue::glue('Beginning generating alphashapes for {acq} with one core'))
  }
  # Generate alphashape dataframes
  ashapes <- catalog_apply(tree_ctg, get_alphashape_metrics)

  print(glue::glue('Finished alphashape computation for {acq}'))
  tictoc::toc()

}

