# Load required packages

library(sf)
library(tidyverse)
library(gridExtra)
library(silvtools)
library(patchwork)

# Internal Functions

# Define function to detect acquisition type of a project
acquisition_type <- function(proj_dir){
  # Initialize acquisition as NULL
  acq <- NULL
  # If the project name contains 'DAP', it's a DAP type acquisition
  if(stringr::str_detect(basename(proj_dir), pattern = 'DAP')){
    is_dap <- TRUE
    # Construct the acquisition name based on the project name
    acq <- glue::glue('DAP22_{stringr::str_replace(basename(proj_dir), pattern = "-DAP", replacement = "")}')
  } else{
    # If the project name does not contain 'DAP', it's a ULS type acquisition
    is_dap <- FALSE
    # Construct the acquisition name based on the project name
    acq <- glue::glue('ULS22_{basename(proj_dir)}')
  }
  return(acq)
}

# Define a function to read alphashape tree tops from a given directory
read_alphashape_tree_tops <- function(proj_dir, acq){
  # List all CSV files in the directory
  csvs <- list.files(glue::glue('{proj_dir}/output/crowns/ashapes'), pattern = '.csv', full.names = T)

  # If there are multiple CSV files, read and combine them all
  ashapes <- if(length(csvs) > 1){
    print(glue::glue('More than one alphashape CSV is present; reading and binding all {length(csvs)} tables'))
    ashapes <- list()
    for(i in 1:length(csvs)){
      ashapes[[i]] <- read.table(csvs[i], header = TRUE, sep = ",")
    }
    do.call(rbind, ashapes)
  } else{
    # If there is only one CSV file, read it
    read.table(glue::glue('{proj_dir}/output/crowns/ashapes/{acq}_chunk_5cmvoxel_ashapes.csv'), header = TRUE, sep = ",")
  }
  return(ashapes)
}

# Function for cleaning up the treetops dataset
clean_up_treetops <- function(ashapes){
  ttops <- ashapes %>%
    # Filter out tree ID's that are not in the names of the shapes
    filter(!ashapes$treeID %in% names(ashapes)) %>%
    # Convert data to simple feature object using X and Y as coordinates
    sf::st_as_sf(coords = c('X','Y'), crs = coords, remove = FALSE) %>%
    # Convert character variables to numeric
    mutate_if(is.character, as.numeric) %>%
    # Convert treeID to integer and create a new column 'point_id' as row number
    mutate(treeID = as.integer(treeID), point_id = row_number()) %>%
    # Group by treeID
    group_by(treeID) %>%
    # Filter to keep only rows where n_points is maximum within each group
    filter(n_points == max(n_points)) %>%
    # Keep only the first row in case of multiple rows with the same max n_points
    filter(row_number() == 1) %>%
    # Ungroup the data
    ungroup()
  # Return the cleaned up treetops
  ttops
}

# Function to cap outliers based on percentiles
cap_outliers <- function(data, metric, upper_cutoff_percentile, lower_cutoff_percentile) {
  # Calculate the upper and lower cutoff values based on the percentiles
  upper_cutoff_value <- quantile(data[[metric]], upper_cutoff_percentile / 100)
  lower_cutoff_value <- quantile(data[[metric]], lower_cutoff_percentile / 100)

  # Cap the values using the cutoff values
  capped_data <- data %>%
    # Cap values above the upper cutoff value
    mutate(across({{metric}}, ~ifelse(. > upper_cutoff_value, upper_cutoff_value, .), .names = "{col}")) %>%
    # Cap values below the lower cutoff value
    mutate(across({{metric}}, ~ifelse(. < lower_cutoff_value, lower_cutoff_value, .), .names = "{col}"))

  # Count the number of outliers
  num_up_outliers <- sum(data[[metric]] > upper_cutoff_value)
  num_low_outliers <- sum(data[[metric]] < lower_cutoff_value)

  # Print the number of capped values
  print(glue::glue("For {metric}, there were {num_up_outliers} rows above the upper cap value, and their values have been altered,
       there were {num_low_outliers} rows below the lower cap value, and their values have been altered,
       in total {num_up_outliers + num_low_outliers} trees had their values capped."))

  # Return the capped dataset
  return(capped_data)
}

filter_unlikely_trees <- function(ttops){
  # Filter out trees that have vol_concave less than 10, n_points less than 50, and Zmax greater than 40m
  ttops %>% filter(vol_concave > 10) %>% filter(n_points > 50) %>% filter(Zmax < 40)
}

# Function to chunk data into smaller pieces
chunk_data <- function(df, size = 10000){
  # Calculate the number of rows and the number of chunks needed
  n <- nrow(df)
  chunks <- ceiling(n / size)
  # Divide the row indices into chunks
  indices <- cut(1:n, breaks = chunks, labels = FALSE)
  # Split the data into chunks using the indices
  split(df, indices)
}

# Function to create a buffer around each chunk of data
create_buffer <- function(chunk, df, distance){
  # Create a buffer around each point in the chunk
  buffer <- st_buffer(chunk, distance)
  # Dissolve borders between points in the buffer
  buffer <- st_union(buffer)
  # Add a new column 'buffer_point' to the dataframe and initialize it as FALSE
  df$buffer_point <- FALSE
  # Identify which points from the original data fall within the buffer
  intersects_buffer <- st_intersection(df, buffer, sparse = TRUE)
  # Update 'buffer_point' to TRUE for points that are in the buffer but not in the original chunk
  intersects_buffer$buffer_point[!intersects_buffer$treeID %in% chunk$treeID] <- TRUE
  # Return the updated dataframe
  return(intersects_buffer)
}


# Apply Buffer Function to Each Chunk with a Buffer
calculate_heygi_for_chunks <- function(chunks, df, buffer_distance = 10, comp_input = 'vol_concave', maxR = 6) {
  # Apply the create_buffer function to each chunk, then apply the heygi_cindex function
  results <- chunks %>%
    purrr::map(~ create_buffer(.x, df, buffer_distance)) %>%
    purrr::map(~ heygi_cindex(.x, comp_input = comp_input, maxR = maxR))

  # Combine the results into a single dataframe
  df_processed <- do.call(rbind, results)
  # Filter out the rows where buffer_point is FALSE
  df_filtered <- filter(df_processed, buffer_point == FALSE)
  # Return the processed dataframe
  return(df_processed)
}


# Calculate Linear BAI Models
calculate_linear_bai_models <- function(treetops) {
  # Compute several new variables related to Basal Area Increment (BAI)
  treetops <- treetops %>% mutate(log_sum_bai_5 = (4.75962 + 0.75351 * log(vol_concave)), sum_bai_5 = exp(log_sum_bai_5),
                                  norm_cindex = ((cindex - min(cindex))/(max(cindex) - min(cindex))),
                                  norm_sum_bai_5 = ((sum_bai_5 - min(sum_bai_5))/(max(sum_bai_5) - min(sum_bai_5))),
                                  GCR_sbai5 = norm_sum_bai_5/norm_cindex)
  # Return the updated treetops data
  return(treetops)
}

# Function to normalize a numeric vector
normalize <- function(x) {
  # Subtract the minimum and divide by the range of the vector
  (x - min(x)) / (max(x) - min(x))
}

# Calculate GC Ratio
calculate_gc_ratio <- function(treetops, weight_growth = 0.5) {
  # Weight for competition is one minus the weight for growth
  weight_competition = 1 - weight_growth
  # Compute the G_C_score as a weighted sum of normalized sum_bai_5 and normalized cindex
  treetops$G_C_score <- weight_growth * normalize(treetops$sum_bai_5) + weight_competition * normalize(treetops$cindex)
  # Return the updated treetops data
  return(treetops)
}

# Function to calculate Normalized Difference Growth Competition Index (NDGCI)
calculate_ndgci <- function(norm_growth, norm_competition) {
  ndgci = (norm_growth - norm_competition) / (norm_growth + norm_competition)
  return(ndgci)
}

# New Index
calculate_gci <- function(norm_growth, norm_competition){
  # New index should range from 0 to 1
  gci = ((1 - norm_growth) + norm_competition)/2
  # Tree with value of 1 - slow growing, high competition
  # Tree with value of 0 - fast growing, low competition
  # Middle range of 0.5 - both low growth/comp and high growth/comp
}

# Function to read the boundary file
read_boundary_file <- function(proj_dir, clip_to_bdy = T){
  # If clip_to_bdy is FALSE, return NULL
  if(!clip_to_bdy){
    return(NULL)
  }

  # Find the boundary file in the input/vector/bdy directory
  bdy_fn <- list.files(glue::glue('{proj_dir}/input/vector/bdy'),
                       pattern = 'bdy.gpkg',
                       full.names = T)
  # If more than one boundary file is found, stop the function
  if(length(bdy_fn) > 1){
    stop("More than one boundary gpkg found in input/vector/bdy")
  }

  # If no boundary file is found and clip_to_bdy is TRUE, stop the function
  if(length(bdy_fn) == 0){
    stop("No boundary gpkg found in input/vector/bdy, but clip_to_bdy = T")
  }

  # Read the boundary file
  bdy <- st_read(bdy_fn, quiet = T)

  # Return the boundary data
  return(bdy)
}

# Function to clip the grids to the boundary
clip_to_boundary <- function(grids, bdy){
  # Intersect the grids with the boundary
  grids <- st_intersection(grids, st_geometry(bdy))

  # Print a message
  print('Grid cropped to block boundary')

  # Return the cropped grids
  return(grids)
}


# Function to compute percentile breaks
percentile_breaks <- function(data, percentiles) {
  # Apply the quantile function to each percentile
  sapply(percentiles, function(p) quantile(data, probs = p/100, na.rm = TRUE))
}

# Function to generate custom labels for breaks
custom_labels <- function(breaks) {
  # For each break, generate a label showing the range of values covered
  sapply(seq_along(breaks)[-1], function(i) {
    lower_break <- ifelse(breaks[i-1] < 1, round(breaks[i-1], 3), round(breaks[i-1], 0))
    upper_break <- ifelse(breaks[i] < 1, round(breaks[i], 3), round(breaks[i], 0))
    paste0(format(lower_break, big.mark = ","), " - ", format(upper_break, big.mark = ","))
  })
}


create_plot <- function(data, fill_var, breaks, color_scale, scale_name, crs, decimal_places = 2) {
  # Define a custom labeling function
  custom_labels <- function(breaks) {
    labels <- scales::label_number(big.mark = ",", accuracy = 10^(-decimal_places))
    labels(breaks)
  }


  plot <- ggplot() +
    geom_sf(data = data, aes_string(fill = paste0("cut(", fill_var, ", breaks = breaks)")), color = "black", size = 0.1) +
    scale_fill_manual(name = scale_name,
                      values = color_scale,
                      labels = custom_labels(breaks), drop = TRUE) +
    theme_minimal() +
    theme(legend.position = "bottom", legend.direction = "horizontal") +
    coord_sf(crs = crs)

  return(plot)
}


# Set up Project

# List directories (each is one acquisiton of ULS/DAP)
blocks_dir <- list.dirs('G:/Quesnel_2022/process', recursive = FALSE)
# Omit these already processed blocks from processing stream
processed <- c('CT1','CT2','CT3','CT4','CT5')
blocks_dir <- blocks_dir[basename(blocks_dir) %in% processed]
# Global variables
coords <- 26910
# Clip resulting grid to stand shape?
clip_to_bdy <- TRUE
# Metrics which need to be capped
metrics <- c('vol_concave','vol_convex')
# Percentiles to cap values by
upper_cutoff_percentile <- 99
lower_cutoff_percentile <- 1


# Area of summary grid tiles




generate_grid <- function(proj_dir, clip_trees_to_bdy = TRUE ,clip_grid_to_bdy = TRUE, grid_area = 2500){

tictoc::tic()

# Define Acquistion Name

acq <- acquisition_type(proj_dir = proj_dir)

print(glue::glue('Summarizing tree tops across {grid_area} m2 hexagonal grid for {acq}'))

# Load tree tops with computed alphashapes

ashapes <- read_alphashape_tree_tops(proj_dir, acq)

# Filter duplicate alphashapes

ttops <- ashapes %>% clean_up_treetops()

print(glue::glue('Filtered out {nrow(ashapes) - nrow(ttops)} duplicate alphashapes ({round(((nrow(ashapes) - nrow(ttops))/nrow(ashapes)) * 100)}%)'))

# Filter erroneous segmentations

ttops_filt <- ttops %>% filter_unlikely_trees()

print(glue::glue('Filtered out {nrow(ttops) - nrow(ttops_filt)} erroneous treetops ({round(((nrow(ttops) - nrow(ttops_filt))/nrow(ttops_filt)) * 100)}%)'))

# Cap extreme values

ttops_capped <- ttops_filt

for (metric in metrics) {
  ttops_capped <- cap_outliers(ttops_capped, metric, upper_cutoff_percentile, lower_cutoff_percentile)
}

ttops <- ttops_capped

# Generate heygi competition index

chunks <- ttops %>% chunk_data()

ttops_c <- calculate_heygi_for_chunks(chunks = chunks,
                                      df = ttops)

# If clip_trees_to_bdy is TRUE, clip the tree tops to the boundary

if(clip_trees_to_bdy){
  bdy <- read_boundary_file(proj_dir)
  ttops_c <- clip_to_boundary(ttops_c, bdy)
  clip <- 'tree_clipped'
  print('Tree tops clipped to boundary')
} else{
  clip <- 'tree_unclipped'
}

# Cap competition index values

ttops <- cap_outliers(ttops_c, metric = 'cindex', upper_cutoff_percentile = 99, lower_cutoff_percentile = 1)

# Now treetops all have crown volume and a cindex value we can use to calculate BAI

ttops <- ttops %>%
  calculate_linear_bai_models()

# Calculate Growth:Competition Indicies

# Original Index
ttops$NDGCI <- calculate_ndgci(ttops$norm_sum_bai_5, ttops$norm_cindex)
# New Index
ttops$GCI <- calculate_gci(ttops$norm_sum_bai_5, ttops$norm_cindex)

# Write final treetops to disk

# Drop extra X and Y columns
ttops <- ttops %>% select(-x, -y)

st_write(ttops, glue::glue('{proj_dir}/output/vector/treetops/{acq}_ttops_lmf2_cut{upper_cutoff_percentile}-{lower_cutoff_percentile}_{clip}_bai_cindex_gci.gpkg'),  append = FALSE)

print(glue::glue('Wrote {nrow(ttops)} treetops to disk for {acq}'))

# Summarize tree tops across regular hexagonal grids
print(glue::glue('Summarizing tree tops across {grid_area} m2 hexagonal grids'))
# Basal area increment grid
bai_grid <- silvtools::summarize_grid(points_sf = ttops, grid_area = grid_area, grid_shape = 'hexagon', summary_var = 'sum_bai_5', summary_fun = 'sum')
# Competition Index grid
cindex_grid <- silvtools::summarize_grid(points_sf = ttops, grid_area = grid_area, grid_shape = 'hexagon', summary_var = 'cindex', summary_fun = 'sum')
# Old Growth:Competition Index grid
NDGCI_grid_sum <- silvtools::summarize_grid(points_sf = ttops, grid_area = grid_area, grid_shape = 'hexagon', summary_var = 'NDGCI', summary_fun = 'sum')
NDGCI_grid_mean <- silvtools::summarize_grid(points_sf = ttops, grid_area = grid_area, grid_shape = 'hexagon', summary_var = 'NDGCI', summary_fun = 'mean')
# New Growth:Competition Index grid
GCI_grid_sum <- silvtools::summarize_grid(points_sf = ttops, grid_area = grid_area, grid_shape = 'hexagon', summary_var = 'GCI', summary_fun = 'sum')
GCI_grid_mean <- silvtools::summarize_grid(points_sf = ttops, grid_area = grid_area, grid_shape = 'hexagon', summary_var = 'GCI', summary_fun = 'mean')

# Join grids together
grids <- st_join(bai_grid, cindex_grid, left = TRUE, largest = TRUE) %>% mutate(sum_sum_bai_5_m2 = sum_sum_bai_5/1000000)
grids <- st_join(grids, NDGCI_grid_sum, left = TRUE, largest = TRUE)
grids <- st_join(grids, NDGCI_grid_mean, left = TRUE, largest = TRUE)
grids <- st_join(grids, GCI_grid_sum, left = TRUE, largest = TRUE)
grids <- st_join(grids, GCI_grid_mean, left = TRUE, largest = TRUE)
grids <- grids %>% dplyr::select(-contains("id"))

# Clip to block boundary if required
if(clip_to_bdy){

bdy <- read_boundary_file(proj_dir)


dir.create(glue::glue('{proj_dir}/output/vector/grids'), showWarnings = FALSE)
st_write(grids, glue::glue('{proj_dir}/output/vector/grids/{acq}_{round(grid_area)}m_cut{upper_cutoff_percentile}-{lower_cutoff_percentile}_{clip}_grids.gpkg'), driver = 'GPKG', append = FALSE)

grids <- clip_to_boundary(grids, bdy)
st_write(grids, glue::glue('{proj_dir}/output/vector/grids/{acq}_{round(grid_area)}m_cut{upper_cutoff_percentile}-{lower_cutoff_percentile}_{clip}_grids_clip.gpkg'), driver = 'GPKG', append = FALSE)

print(glue::glue('Clipped grids to block boundary and wrote to disk for {acq}'))

}

grids <- grids %>%
  mutate(norm_bai = sum_sum_bai_5 / max(sum_sum_bai_5),
         norm_cindex = sum_cindex / max(sum_cindex),
         norm_diff = abs(norm_bai - norm_cindex))

grids <- grids %>%
  mutate(geom_mean_index = sqrt(norm_bai * norm_cindex),
         diff_geom_mean_index = sqrt(norm_bai * norm_cindex * norm_diff))

tictoc::toc()

return(grids)

}


ct1 <- generate_grid("G:/Quesnel_2022/process/CT1", clip_trees_to_bdy = T, clip_grid_to_bdy = T, grid_area = 2500)
ct2 <- generate_grid("G:/Quesnel_2022/process/CT2", clip_trees_to_bdy = T, clip_grid_to_bdy = T, grid_area = 2500)
ct3 <- generate_grid("G:/Quesnel_2022/process/CT3", clip_trees_to_bdy = T, clip_grid_to_bdy = T, grid_area = 2500)
ct4 <- generate_grid("G:/Quesnel_2022/process/CT4", clip_trees_to_bdy = T, clip_grid_to_bdy = T, grid_area = 2500)
ct5 <- generate_grid("G:/Quesnel_2022/process/CT5", clip_trees_to_bdy = T, clip_grid_to_bdy = T, grid_area = 2500)

utm_crs <- st_crs("+init=EPSG:26910")

grids <- ct1

grids <- grids %>% mutate(area_ha = st_area(geometry)/10000, bai_m2_per_ha_5yr = sum_sum_bai_5_m2/area_ha, bai_m2_per_ha = bai_m2_per_ha_5yr/5)
#grids <- grids %>% filter(as.numeric(area_ha) >= 0.20)


# break_percentiles <- c(0, 10, 50, 75, 85, 95, 100)
break_percentiles <- c(0, 10, 25, 50, 75, 90, 100)

# Create a plot for sum_sum_bai_5
bai_breaks <- percentile_breaks(grids$sum_sum_bai_5_m2, break_percentiles)
plot_bai <- create_plot(data = grids,
                        fill_var = "sum_sum_bai_5_m2",
                        breaks = bai_breaks,
                        color_scale = colorRampPalette(c(alpha("white", 0.9), alpha("darkgreen", 0.9)), alpha = TRUE)(length(bai_breaks) - 1),
                        scale_name = "Cumulative Basal\nArea of Last 5 Years (m2/5 yr)",
                        crs = utm_crs,
                        decimal_places = 2)

# Create a plot for sum_cindex
cindex_breaks <- percentile_breaks(grids$sum_cindex, break_percentiles)
plot_cindex <- create_plot(data = grids,
                           fill_var = "sum_cindex",
                           breaks = cindex_breaks,
                           color_scale = colorRampPalette(c("white", "darkred"))(length(cindex_breaks) - 1),
                           scale_name = "Cumulative Competition Index",
                           crs = utm_crs,
                           decimal_places = 2)


dgindex_breaks <- percentile_breaks(grids$mean_NDGCI, break_percentiles)
color_scale_viridis <- viridisLite::viridis(length(dgindex_breaks) - 1)

plot_ndgci <- create_plot(data = grids,
                          fill_var = "mean_NDGCI",
                          breaks = dgindex_breaks,
                          color_scale = color_scale_viridis,
                          scale_name = "Mean NDGCI",
                          crs = utm_crs,
                          decimal_places = 2)

gciindex_breaks <- percentile_breaks(grids$mean_GCI, break_percentiles)
color_scale_viridis <- viridisLite::viridis(length(gciindex_breaks) - 1)
plot_gci <- create_plot(data = grids,
                        fill_var = "mean_GCI",
                        breaks = gciindex_breaks,
                        color_scale = color_scale_viridis,
                        scale_name = "Mean GCI",
                        crs = utm_crs,
                        decimal_places = 2)

# Plot GCI Sum

gci_sum_breaks <- percentile_breaks(grids$sum_GCI, break_percentiles)
color_scale_viridis <- viridisLite::viridis(length(gci_sum_breaks) - 1)

plot_gci_sum <- create_plot(data = grids,
                            fill_var = "sum_GCI",
                            breaks = gci_sum_breaks,
                            color_scale = color_scale_viridis,
                            scale_name = "Sum GCI",
                            crs = utm_crs,
                            decimal_places = 2)

library(patchwork)

plot_bai + plot_cindex + plot_gci


# Summary Metrics for Grids

grid <- ct2

calculate_summary_stats <- function(grid) {
  # Calculate summary statistics for each column
  summary_stats <- data.frame(
    Variable = colnames(grid),
    Mean = colMeans(grid, na.rm = TRUE),
    StdDev = sapply(grid, sd, na.rm = TRUE),
    Min = sapply(grid, min, na.rm = TRUE),
    Max = sapply(grid, max, na.rm = TRUE)
  )

  return(summary_stats)
}

# Call the function with your grid
grid_summary <- calculate_summary_stats(grid)

# Call the function with your list of grids

block_grids <- list(ct1, ct2, ct3, ct4, ct5)

grids_all <- do.call(rbind, block_grids) %>% st_drop_geometry()

grid_summary <- calculate_summary_stats(grids_all)

write.csv(grid_summary, 'D:/Proposal_2022/Thinning Paper/Figures/grid_summary_stats_oct.csv')

