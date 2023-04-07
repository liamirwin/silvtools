
library(lidR)
library(glue)
library(future)
# Define a filter function for catalog_apply
filter_tree <- function(chunk, output_dir) {
  # Read the LAS data
  las <- readLAS(chunk)

  # Check if the data is not empty
  if (is.empty(las)) return(NULL)

  # Get unique treeIDs
  tree_ids <- unique(las$treeID)

  # Loop through each treeID and save them as separate LAZ files
  for (tree_id in tree_ids) {
    tree_las <- filter_poi(las, treeID == tree_id)
    output_file <- file.path(output_dir, glue("treeID_{tree_id}_{round(mean(las@data$X))}X_{round(mean(las@data$Y))}Y.laz"))
    if(is.empty(tree_las)){
      print('Empty las file')
    } else{
    writeLAS(tree_las, output_file)
    }
  }
}

blocks_dir <- list.dirs('H:/Quesnel_2022/process', recursive = FALSE)
processed <- c('CT1-T-DAP','CT1','CT5')
# target <- c('CT1')
blocks_dir <- blocks_dir[!basename(blocks_dir) %in% processed]
is_dap <- FALSE


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

tree_dir <- glue::glue('{proj_dir}/output/individual_trees')

if (!dir.exists(tree_dir)) {
  dir.create(tree_dir)
}

ctg <- catalog(glue::glue('{proj_dir}/output/tree_las'))
#opt_filter(ctg) = "-thin_with_voxel 0.1"
opt_laz_compression(ctg) <- T
opt_progress(ctg) <- T
plan(multisession, workers = 3L)

# Process the LAScatalog using catalog_apply
catalog_apply(ctg, filter_tree, output_dir = tree_dir)

}
