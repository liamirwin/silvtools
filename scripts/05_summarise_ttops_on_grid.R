# Summarize point data on grid

library(sf)
library(tidyverse)
library(gridExtra)
library(silvtools)
library(patchwork)

# List directories (each is one acquisiton of ULS/DAP)
blocks_dir <- list.dirs('G:/Quesnel_2022/process', recursive = FALSE)
# Omit these already processed blocks from processing stream
processed <- c('CT1','CT2','CT3','CT4','CT5')
blocks_dir <- blocks_dir[basename(blocks_dir) %in% processed]
# ULS or DAP?
is_dap <- FALSE
# ---- Project setup ----
grid_area <- 2500
grid_shape <- 'hexagon' # or 'square'
plot_list <- list()
grid_list <- list()

for(i in 1:length(blocks_dir)){

  tictoc::tic()

  if(length(blocks_dir) == 1){
    i <- 1
  }

acq = NULL
proj_dir <- blocks_dir[i]
coords <- 26910
clip_to_bdy <- TRUE


if(is.null(acq)){
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
}


# Read Alphashape tree tops

if(length(list.files(glue::glue('{proj_dir}/output/crowns/ashapes'), pattern = '.csv')) > 1){
  csvs <- list.files(glue::glue('{proj_dir}/output/crowns/ashapes'), pattern = '.csv', full.names = T)
  print(glue::glue('More than one alphashape CSV is present; reading and binding all {length(csvs)} tables'))
  ashapes <- list()
  for(i in 1:length(csvs)){
    ashapes[[i]] <- read.table(csvs[i], header = TRUE, sep = ",")
  }
  ashapes <- do.call(rbind, ashapes)
} else{
ashapes <- read.table(glue::glue('{proj_dir}/output/crowns/ashapes/{acq}_chunk_5cmvoxel_ashapes.csv'), header = TRUE, sep = ",")

}

# Clean up treetops
ttops <- ashapes %>%
  # Remove headers
  filter(!ashapes$treeID %in% names(ashapes)) %>%
  # Convert to SF points
  sf::st_as_sf(coords = c('X','Y'), crs = coords, remove = FALSE) %>%
  # Convert value data types from character to numeric
  mutate_if(is.character, as.numeric) %>%
  mutate(treeID = as.integer(treeID), point_id = row_number()) %>%
  # Group alphashapes with common treeIDs
  group_by(treeID) %>%
  # Extract tree with highest number of returns
  filter(n_points == max(n_points)) %>%
  filter(row_number() == 1) %>%
  # Ungroup alphashapes
  ungroup()

print(glue::glue('Filtered out {nrow(ashapes) - nrow(ttops)} duplicate alphashapes ({round(((nrow(ashapes) - nrow(ttops))/nrow(ashapes)) * 100)}%)'))

# Filter out unlikely trees
ttops_filt <- ttops %>% filter(vol_concave > 10) %>% filter(n_points > 50) %>% filter(Zmax < 40)

print(glue::glue('Filtered out {nrow(ttops) - nrow(ttops_filt)} erroneous treetops ({round(((nrow(ttops) - nrow(ttops_filt))/nrow(ttops_filt)) * 100)}%)'))

ttops <- ttops_filt

# Clamp outliers and Convert Alphashape Attributed Treetops to Summary Grid for BAI and Competition
# Adjust Outliers based on quantiles
# Set the desired metrics to cap
metrics <- c('vol_concave','vol_convex')

# Replace XX with your chosen percentile (e.g., 95, 99, etc.)
upper_cutoff_percentile <- 99
lower_cutoff_percentile <- 1

# Function to cap outliers based on chosen quantile
cap_outliers <- function(data, metric,  upper_cutoff_percentile, lower_cutoff_percentile) {
  # Calculate the upper and lower cutoff values based on the percentiles
  upper_cutoff_value <- quantile(data[[metric]], upper_cutoff_percentile / 100)
  lower_cutoff_value <- quantile(data[[metric]], lower_cutoff_percentile / 100)

  # Cap the values using the cutoff values
  capped_data <- data %>%
    mutate(across({{metric}}, ~ifelse(. > upper_cutoff_value, upper_cutoff_value, .), .names = "{col}")) %>%
    mutate(across({{metric}}, ~ifelse(. < lower_cutoff_value, lower_cutoff_value, .), .names = "{col}"))

  # Count the number of outliers
  num_up_outliers <- sum(data[[metric]] > upper_cutoff_value)
  num_low_outliers <- sum(data[[metric]] < lower_cutoff_value)

  # Print the number of capped values
  print(glue::glue("For {metric}, there were {num_up_outliers} rows above the upper cap value, and their values have been altered,
       there were {num_low_outliers} rows below the lower cap value, and their values have been altered,
       in total {num_up_outliers + num_low_outliers} trees had their values capped."))

  # Create pre-outlier and post-outlier capping histograms
  p1 <- plot_outlier_distribution(data, metric, percentiles = c(upper_cutoff_percentile, lower_cutoff_percentile)) + labs(title = glue::glue('Pre-outlier capping histogram of {metric}'))
  p2 <- plot_outlier_distribution(capped_data, metric, percentiles = c(upper_cutoff_percentile, lower_cutoff_percentile)) + labs(title = glue::glue('Post-outlier capping histogram of {metric}'))

  # Print the histograms
  grid.arrange(p1, p2, ncol = 2)
  # Return the capped dataset
  return(capped_data)
}

# Cap outliers for each metric
ttops_capped <- ttops

for (metric in metrics) {
  ttops_capped <- cap_outliers(ttops_capped, metric, upper_cutoff_percentile, lower_cutoff_percentile)
}

ttops <- ttops_capped

st_write(ttops, glue::glue('{proj_dir}/output/vector/treetops/{acq}_ashape_ttops.gpkg'), append = FALSE)

# Calculate Heygi index for each tree
print(glue::glue('Calculating Heygi index for {nrow(ttops)} tree tops'))

# Large Tree top datasets need to be split up

# Define the maximum chunk size (number of trees)
chunk_size <- 10000

# Create a function to chunk the data
chunk_data <- function(df, size){
  n <- nrow(df)
  chunks <- ceiling(n / size)
  indices <- cut(1:n, breaks = chunks, labels = FALSE)
  split(df, indices)
}

# Chunk the data
chunks <- chunk_data(ttops, chunk_size)

# Create a buffer function
create_buffer <- function(chunk, df, distance){
  # Create buffer around each point in the chunk
  buffer <- st_buffer(chunk, distance)
  # Dissolve borders between points
  buffer <- st_union(buffer)
  # Add buffer_point column to all points in the dataframe and set it to FALSE
  df$buffer_point <- FALSE
  # Get a logical matrix indicating which points from the original data fall within the buffer
  intersects_buffer <- st_intersection(df, buffer, sparse = TRUE)
  # Update buffer_point to TRUE for points that are in the buffer but not in the original chunk
  intersects_buffer$buffer_point[!intersects_buffer$treeID %in% chunk$treeID] <- TRUE
  return(intersects_buffer)
}

# Define the buffer distance (in the units of the CRS)
buffer_distance <- 10  # Adjust as needed

# Apply the function to each chunk with a buffer
results <- chunks %>%
  purrr::map(~ create_buffer(.x, ttops, buffer_distance)) %>%
  purrr::map(~ heygi_cindex(.x, comp_input = 'vol_concave', maxR = 6))

# Combine the chunks back together
ttops_processed <- do.call(rbind, results)

# Filter out the buffer points
ttops_final <- filter(ttops_processed, buffer_point == FALSE)

ttops <- ttops_final

plot_outlier_distribution(ttops, metric = 'cindex', percentiles = c(1, 5, 99.5))
plot_outlier_distribution(ttops, metric = 'vol_concave', percentiles = c(1, 5, 99.5))

ttops_capped <- cap_outliers(ttops, metric = 'cindex', upper_cutoff_percentile = 99, lower_cutoff_percentile = 1)

ttops <- ttops_capped

# Apply linear bai models and calculate growth competition ratios/normalize

ttops <- ttops %>% mutate(log_sum_bai_5 = (4.76 + 0.75 * log(vol_concave)), sum_bai_5 = exp(log_sum_bai_5),
                          norm_cindex = ((cindex - min(cindex))/(max(cindex) - min(cindex))),
                          norm_sum_bai_5 = ((sum_bai_5 - min(sum_bai_5))/(max(sum_bai_5) - min(sum_bai_5))),
                          GCR_sbai5 = norm_sum_bai_5/norm_cindex)

normalize <- function(x) { (x - min(x)) / (max(x) - min(x)) }
ttops$G_C_difference <- normalize(ttops$sum_bai_5) - normalize(ttops$cindex)
ttops$G_C_ratio <- normalize(ttops$sum_bai_5) / normalize(ttops$cindex)
ttops$NDGCI <- (normalize(ttops$sum_bai_5) - normalize(ttops$cindex))/normalize(ttops$sum_bai_5) + normalize(ttops$cindex)
ttops$G_C_ratio[is.infinite(ttops$G_C_ratio)] <- 1


weight_growth = 0.5  # adjust as needed
weight_competition = 1 - weight_growth
ttops$G_C_score <- weight_growth * normalize(ttops$sum_bai_5) + weight_competition * normalize(ttops$cindex)

ggplot(ttops, aes(x = cindex, y = sum_bai_5, color = G_C_difference)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Competition", y = "Growth", color = "Difference Index", title = "Growth vs Competition with Difference Index") +
  scale_color_gradient(low = "blue", high = "red")

ggplot(ttops, aes(x = cindex, y = sum_bai_5, color = G_C_ratio)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Competition", y = "Growth", color = "Ratio Index", title = "Growth vs Competition with Ratio Index") +
  scale_color_gradient(low = "blue", high = "red")

ggplot(ttops, aes(x = cindex, y = sum_bai_5, color = G_C_score)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Competition", y = "Growth", color = "Composite Score", title = "Growth vs Competition with Composite Score") +
  scale_color_gradient(low = "blue", high = "red")

ggplot(ttops, aes(x = cindex, y = sum_bai_5, color = NDGCI)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Competition", y = "Growth", color = "Composite Score", title = "Growth vs Competition with Composite Score") +
  scale_color_gradient(low = "blue", high = "red")



bai_grid <- silvtools::summarize_grid(points_sf = ttops, grid_area = grid_area, grid_shape = 'hexagon', summary_var = 'sum_bai_5', summary_fun = 'sum')
cindex_grid <- silvtools::summarize_grid(points_sf = ttops, grid_area = grid_area, grid_shape = 'hexagon', summary_var = 'cindex', summary_fun = 'sum')
sum_GCR_grid <- silvtools::summarize_grid(points_sf = ttops, grid_area = grid_area, grid_shape = 'hexagon', summary_var = 'NDGCI', summary_fun = 'sum')
mean_GCR_grid <- silvtools::summarize_grid(points_sf = ttops, grid_area = grid_area, grid_shape = 'hexagon', summary_var = 'NDGCI', summary_fun = 'mean')

grids <- st_join(bai_grid, cindex_grid, left = TRUE, largest = TRUE) %>% mutate(sum_sum_bai_5_m2 = sum_sum_bai_5/1000000)
grids <- st_join(grids, GCR_grid, left = TRUE, largest = TRUE)
grids <- grids %>% select(-contains("id"))

if(clip_to_bdy == T){

  bdy_fn <- list.files(glue::glue('{proj_dir}/input/vector/bdy'),
                       pattern = 'bdy.gpkg',
                       full.names = T)

  if(length(bdy_fn) > 1){
    stop("More than one boundary gpkg found in input/vector/bdy")
  }

  if(length(bdy_fn) == 0){
    stop("No boundary gpkg found in input/vector/bdy, but clip_to_bdy = T")
  }

  bdy <- st_read(bdy_fn, quiet = T)

  # Clip the grids with the boundary
  grids <- st_intersection(grids, st_geometry(bdy))

  print('Grid cropped to block boundary')

}


# Normalize Grids with Index
# First Normalize BAI and CINDEX values based on their maximum
grids <- grids %>%
  mutate(norm_bai = sum_sum_bai_5 / max(sum_sum_bai_5),
         norm_cindex = sum_cindex / max(sum_cindex),
         norm_diff = abs(norm_bai - norm_cindex))
# Then take the geometric mean - useful for dealing with variables that are different orders of magnitude
# High values = High Growth/Competition (Don't Thin)
# Intermediate Value =
# Low Values = Low Growth/Competition (Don't Thin)
grids <- grids %>%
  mutate(geom_mean_index = sqrt(norm_bai * norm_cindex),
         diff_geom_mean_index = sqrt(norm_bai * norm_cindex * norm_diff))


# Plot Grids

percentile_breaks <- function(data, percentiles) {
  sapply(percentiles, function(p) quantile(data, probs = p/100, na.rm = TRUE))
}

percentiles <- c(0, 5, 15, 25, 50, 75, 90, 95, 100)
percentiles <- c(0, 5, 25, 50, 75, 90, 95, 100)
# Calculate breaks for sum_mean_bai_5 and sum_cindex
bai_breaks <- percentile_breaks(grids$sum_sum_bai_5, percentiles)
bai_m2_breaks <- percentile_breaks(grids$sum_sum_bai_5_m2, percentiles)
cindex_breaks <- percentile_breaks(grids$sum_cindex, percentiles)
gindex_breaks <- percentile_breaks(grids$geom_mean_index, percentiles)
dgindex_breaks <- percentile_breaks(grids$diff_geom_mean_index, percentiles)

# Outlier Distributions of Grid Cells

# plot_outlier_distribution(grids, metric = 'sum_sum_bai_5', percentiles = percentiles)
# plot_outlier_distribution(grids, metric = 'sum_cindex', percentiles = percentiles)
# plot_outlier_distribution(grids, metric = 'sum_sum_bai_5_m2', percentiles = percentiles)
# plot_outlier_distribution(grids, metric = 'norm_bai', percentiles = percentiles)
# plot_outlier_distribution(grids, metric = 'norm_cindex', percentiles = percentiles)
# plot_outlier_distribution(grids, metric = 'norm_diff', percentiles = percentiles)
# plot_outlier_distribution(grids, metric = 'geom_mean_index', percentiles = percentiles)
# plot_outlier_distribution(grids, metric = 'diff_geom_mean_index', percentiles = percentiles)

# Function to generate custom labels

custom_labels <- function(breaks) {
  sapply(seq_along(breaks)[-1], function(i) {
    lower_break <- ifelse(breaks[i-1] < 1, round(breaks[i-1], 3), round(breaks[i-1], 0))
    upper_break <- ifelse(breaks[i] < 1, round(breaks[i], 3), round(breaks[i], 0))
    paste0(format(lower_break, big.mark = ","), " - ", format(upper_break, big.mark = ","))
  })
}

# # Prepare CHM for Plotting
#
# chm <- terra::rast(glue::glue('{proj_dir}/output/raster/chm/{acq}_chm_fill_p2r_0.1m.tif')) %>%
#   terra::mask(terra::vect(bdy)) %>% terra::crop(terra::vect(bdy))
#
# chm_df <- as.data.frame(chm, xy = T)

# Define the desired CRS
utm_crs <- st_crs("+init=EPSG:26910")


# Create a plot for the CHM raster layer
# plot_chm <- ggplot() +
#   geom_raster(data = chm_df, aes(x = x, y = y, fill = Z)) +
#   scale_fill_gradientn(name = "Canopy Height Model",
#                        colours = c("black", "white")) +
#   theme_void()
#
# plot_chm <- ggplot() +
#   geom_raster(data = chm_df, aes(x = x, y = y, fill = Z)) +
#   scale_fill_gradientn(name = "Canopy Height Model",
#                        colours = c("black", "white")) +
#   theme_minimal() +
#   theme(legend.position = "bottom",
#         legend.direction = "horizontal") +
#   coord_sf(crs = utm_crs)


# Create a plot for sum_mean_bai_5 with UTM 10N coordinates
plot_bai <- ggplot() +
  geom_sf(data = grids, aes(fill = cut(sum_sum_bai_5, breaks = bai_breaks)), color = "black", size = 0.1) +
  scale_fill_manual(name = "Cumulative Basal\nArea of Last 5 Years (mm/5 yr)",
                    values = colorRampPalette(c(alpha("white", 0.9), alpha("darkgreen", 0.9)), alpha = TRUE)(length(bai_breaks) - 1),
                    labels = custom_labels(bai_breaks)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal") +
  coord_sf(crs = utm_crs)

# Create a plot for sum_cindex with UTM 10N coordinates
plot_cindex <- ggplot() +
  geom_sf(data = grids, aes(fill = cut(sum_cindex, breaks = cindex_breaks)), color = "black", size = 0.1) +
  scale_fill_manual(name = "Cumulative Competition Index",
                    values = colorRampPalette(c("white", "darkred"))(length(cindex_breaks) - 1),
                    labels = custom_labels(cindex_breaks)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal") +
  coord_sf(crs = utm_crs)

# Try plotting index values

threshold_bai <- 0.4
threshold_cindex <- 0.2

grids <- grids %>%
  mutate(scenario = case_when(
    norm_bai >= threshold_bai & norm_cindex >= threshold_cindex ~ "High/High",
    norm_bai < threshold_bai & norm_cindex >= threshold_cindex ~ "Low/High",
    norm_bai >= threshold_bai & norm_cindex < threshold_cindex ~ "High/Low",
    norm_bai < threshold_bai & norm_cindex < threshold_cindex ~ "Low/Low"
  ))

summary_table <- grids %>% st_drop_geometry() %>%
  group_by(scenario) %>%
  summarize(mean_diff_geom_mean_index = mean(diff_geom_mean_index),
            min_diff_geom_mean_index = min(diff_geom_mean_index),
            max_diff_geom_mean_index = max(diff_geom_mean_index),
            .groups = "drop")


plot_index <- ggplot() +
  geom_sf(data = grids, aes(fill = cut(diff_geom_mean_index, breaks = dgindex_breaks)), color = "black", size = 0.1) +
  scale_fill_manual(name = "Differential Geometric Growth/Competition Index",
                    values = colorRampPalette(c("white", "darkblue"))(length(dgindex_breaks) - 1),
                    labels = custom_labels(dgindex_breaks)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal") +
  coord_sf(crs = utm_crs)


plot_index <- ggplot() +
  geom_sf(data = grids, aes(fill = cut(diff_geom_mean_index, breaks = dgindex_breaks)), color = "black", size = 0.1) +
  scale_fill_viridis_d(name = "Differential Geometric Growth/Competition Index",
                     labels = custom_labels(dgindex_breaks),
                     option = "viridis") + # You can change the option to "magma", "plasma", or "inferno" for other Viridis color ramps
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal") +
  coord_sf(crs = utm_crs)

# Arrange the two plots side-by-side


#final_plot <- grid.arrange(plot_bai, plot_cindex, plot_index, ncol = 3)
plots <- list()
plots[[1]] <- plot_bai
plots[[2]] <- plot_cindex
plots[[3]] <- plot_index

plot_list[[i]] <- plots
grid_list[[i]] <- grids

}

grids <- grid_list[sapply(grid_list, function(x) !is.null(x))]
plots <- plot_list[sapply(grid_list, function(x) !is.null(x))]
plots <- discard(plot_list, is.null)

for(i in 1:length(plots)){

  p1 <- plots[[i]][[1]]
  p2 <- plots[[i]][[2]]
  p3 <- plots[[i]][[3]]

  p1 + p2 + p3

  # Prompt the user to press "Enter" to continue
  cat("Press Enter to continue...")
  readline()

}

p1 <- plots[[2]][[1]]
p2 <- plots[[2]][[2]]
p3 <- plots[[2]][[3]]

p1 + p2 + p3



for(i in 1:length(grid_list)){
  print(grid_list[[i]])
}

# Summarize values of grid cells across all five stands





# Compare Sum BAI 5 and sum_cindex at grid level

ggplot(x, aes(x = bai_m2_per_ha, y = sum_cindex)) +
  geom_point() +
  labs(x = "sum_sum_bai_5_m2", y = "sum_cindex") +
  ggtitle("Scatterplot of sum_sum_bai_5 vs sum_cindex")


x <- grid_list[[1]] %>% mutate(area = st_area(geometry), bai_m2_per_ha = (sum_sum_bai_5_m2/5)/(as.numeric(area)/10000))

plot_bai <- ggplot() +
  geom_sf(data = x, aes(fill = cut(sum_sum_bai_5, breaks = bai_breaks)), color = "black", size = 0.1) +
  scale_fill_manual(name = "Cumulative Basal\nArea of Last 5 Years (mm/5 yr)",
                    values = colorRampPalette(c(alpha("white", 0.9), alpha("darkgreen", 0.9)), alpha = TRUE)(length(bai_breaks) - 1),
                    labels = custom_labels(bai_breaks)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal") +
  coord_sf(crs = utm_crs)

plot_bai <- ggplot() +
  geom_sf(data = x, aes(fill = bai_m2_per_ha), color = "black", size = 0.1) +
  scale_fill_manual(name = "Cumulative Basal\nArea of Last 5 Years (mm/5 yr)",
                    values = colorRampPalette(c(alpha("white", 0.9), alpha("darkgreen", 0.9)), alpha = TRUE)(length(bai_breaks) - 1),
                    labels = custom_labels(bai_breaks)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal") +
  coord_sf(crs = utm_crs)



st_write(grid_cv_hex, glue::glue('{proj_dir}/output/vector/{acq}_heygi_bai_ashape_hexgrid_{cellsize}m.gpkg'), append = FALSE)

st_write(grid_cv, glue::glue('{proj_dir}/output/vector/{acq}_heygi_bai_ashape_sqgrid_{cellsize}m.gpkg'), append = FALSE)
x <- ttops_to_grid(proj_dir = blocks_dir[1])
st_write(x, 'H:/Quesnel_2022/scratch/CT1_heygi_ashape_grid.gpkg')
ttops <- st_read('H:/Quesnel_2022/process/CT1-T-DAP/output/crowns/ashapes/DAP22_CT1-T_heygi_ashapes.gpkg')
ttops <- st_read('H:/Quesnel_2022/process/CT1/output/crowns/ashapes/ULS22_CT1_heygi_ashapes.gpkg')
ttops <- ttops %>% filter(vol_concave < 800) %>% filter(n_points > 100) %>% filter(n_points < 100000)


