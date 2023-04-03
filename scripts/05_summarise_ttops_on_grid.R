# Summarize point data on grid

library(sf)
library(tidyverse)
library(gridExtra)


# List directories (each is one acquisiton of ULS/DAP)
blocks_dir <- list.dirs('H:/Quesnel_2022/process', recursive = FALSE)
# Omit these already processed blocks from processing stream
processed <- c('CT3','CT4','CT5')
blocks_dir <- blocks_dir[!basename(blocks_dir) %in% processed]
# ULS or DAP?
is_dap <- FALSE
i = 1
# ---- Project setup ----
acq = NULL
proj_dir <- blocks_dir[i]

# Clamp outliers and Convert Alphashape Attributed Treetops to Summary Grid for BAI and Competition

ttops_to_grid <- function(proj_dir, cellsize = 10, square = TRUE, acq = NULL, clip_to_bdy = FALSE){

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
ashapes <- read.table(glue::glue('{proj_dir}/output/crowns/ashapes/{acq}_chunk_5cmvoxel_ashapes.csv'), header = TRUE, sep = ",")

# Clean up treetops
ttops <- ashapes %>% filter(!ashapes$treeID %in% names(ashapes)) %>% st_as_sf(coords = c('X','Y'), crs = 26910, remove = FALSE) %>% mutate_if(is.character, as.numeric) %>%
  mutate(treeID = as.integer(treeID), point_id = row_number()) %>% group_by(treeID) %>%
  filter(n_points == max(n_points)) %>% filter(row_number() == 1) %>% ungroup()

# Filter out unlikely trees

ttops <- ttops %>% filter(vol_concave > 10) %>% filter(n_points > 50)


# Adjust Outliers based on quantiles
# Set the desired metrics to cap
metrics <- c('vol_concave','vol_convex')

# Replace XX with your chosen percentile (e.g., 95, 99, etc.)
upper_cutoff_percentile <- 99.5
lower_cutoff_percentile <- 0.5

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

# Calculate Heygi index for each tree
print(glue::glue('Calculating Heygi index for {nrow(ttops)} tree tops'))

ttops <- silvtools::heygi_cindex(ttops, comp_input = 'vol_concave', maxR = 6)

plot_outlier_distribution(ttops, metric = 'cindex', percentiles = c(1, 5, 99.5))
plot_outlier_distribution(ttops, metric = 'vol_concave', percentiles = c(1, 5, 99.5))

ttops_capped <- cap_outliers(ttops, metric = 'cindex', upper_cutoff_percentile = 99, lower_cutoff_percentile = 1)

ttops <- ttops_capped

# Apply linear bai models and calculate growth competition ratios/normalize

ttops <- ttops %>% mutate(mean_bai_5 = (369.924 + 3.676 * vol_concave), mean_bai_10 = (450.285 + 3.803 * vol_concave),
                          sum_bai_5 = (4289.062 + 43.452 * vol_concave), sum_bai_10 = (9727.30 + 81.75 * vol_concave),
                          norm_cindex = ((cindex - min(cindex))/(max(cindex) - min(cindex))),
                          norm_mean_bai_5 = ((mean_bai_5 - min(mean_bai_5))/(max(mean_bai_5) - min(mean_bai_5))),
                          norm_mean_bai_10 = ((mean_bai_10 - min(mean_bai_10))/(max(mean_bai_10) - min(mean_bai_10))),
                          norm_sum_bai_5 = ((sum_bai_5 - min(sum_bai_5))/(max(sum_bai_5) - min(sum_bai_5))),
                          norm_sum_bai_10 = ((sum_bai_10 - min(sum_bai_10))/(max(sum_bai_10) - min(sum_bai_10))),

                          GCR_mbai5 = norm_mean_bai_5/norm_cindex,
                          GCR_mbai10 = norm_mean_bai_10/norm_cindex,
                          GCR_sbai5 = norm_sum_bai_5/norm_cindex,
                          GCR_sbai10 = norm_sum_bai_10/norm_cindex)


bai_grid <- silvtools::summarize_grid(points_sf = ttops, grid_area = 1000, grid_shape = 'hexagon', summary_var = 'mean_bai_5', summary_fun = 'sum')
cindex_grid <- silvtools::summarize_grid(points_sf = ttops, grid_area = 1000, grid_shape = 'hexagon', summary_var = 'cindex', summary_fun = 'sum')

grids <- st_join(bai_grid, cindex_grid)

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

return(grids)

}

# Whats going on with our treetops...

# Select attributes for clustering
attributes <- c("Zmax", "Zq999", "Zq99", "Z_mean", "n_points", "vol_concave", "vol_convex", "CV_Z", "CRR", "cindex")
ttops_clustering <- ttops[, attributes] %>% st_drop_geometry()
ttops_clustering <- ttops_clustering[complete.cases(ttops_clustering) & !apply(ttops_clustering, 1, function(x) any(is.infinite(x))), ]


# Perform k-means clustering
set.seed(42)  # For reproducibility
num_clusters <- 3
kmeans_model <- kmeans(ttops_clustering, centers = num_clusters)

# Add cluster labels to the dataframe
ttops_clustering$cluster <- as.factor(kmeans_model$cluster)

# Visualize clusters using a scatter plot (using two attributes, e.g., Zmax and vol_concave)
ggplot(ttops_clustering, aes(x = Zmax, y = vol_concave, color = cluster)) +
  geom_point() +
  labs(title = "K-means Clustering of Top Cindex Trees",
       x = "Zmax",
       y = "Vol Concave") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))




# Normalize Grids with Index
# First Normalize BAI and CINDEX values based on their maximum
grids <- grids %>%
  mutate(norm_bai = sum_mean_bai_5 / max(sum_mean_bai_5),
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
bai_breaks <- percentile_breaks(grids$sum_mean_bai_5, percentiles)
cindex_breaks <- percentile_breaks(grids$sum_cindex, percentiles)
gindex_breaks <- percentile_breaks(grids$geom_mean_index, percentiles)
dgindex_breaks <- percentile_breaks(grids$diff_geom_mean_index, percentiles)

# Outlier Distributions of Grid Cells

plot_outlier_distribution(grids, metric = 'sum_mean_bai_5', percentiles = percentiles)
plot_outlier_distribution(grids, metric = 'sum_cindex', percentiles = percentiles)
plot_outlier_distribution(grids, metric = 'norm_bai', percentiles = percentiles)
plot_outlier_distribution(grids, metric = 'norm_cindex', percentiles = percentiles)
plot_outlier_distribution(grids, metric = 'norm_diff', percentiles = percentiles)
plot_outlier_distribution(grids, metric = 'geom_mean_index', percentiles = percentiles)
plot_outlier_distribution(grids, metric = 'diff_geom_mean_index', percentiles = percentiles)

# Function to generate custom labels

custom_labels <- function(breaks) {
  sapply(seq_along(breaks)[-1], function(i) {
    lower_break <- ifelse(breaks[i-1] < 1, round(breaks[i-1], 3), round(breaks[i-1], 0))
    upper_break <- ifelse(breaks[i] < 1, round(breaks[i], 3), round(breaks[i], 0))
    paste0(format(lower_break, big.mark = ","), " - ", format(upper_break, big.mark = ","))
  })
}

# Prepare CHM for Plotting

chm <- terra::rast(glue::glue('{proj_dir}/output/raster/chm/{acq}_chm_fill_p2r_0.05m.tif')) %>%
  terra::mask(terra::vect(bdy)) %>% terra::crop(terra::vect(bdy))

chm_df <- as.data.frame(chm, xy = T)

# Define the desired CRS
utm_crs <- st_crs("+init=EPSG:26910")


# Create a plot for the CHM raster layer
plot_chm <- ggplot() +
  geom_raster(data = chm_df, aes(x = x, y = y, fill = Z)) +
  scale_fill_gradientn(name = "Canopy Height Model",
                       colours = c("black", "white")) +
  theme_void()

plot_chm <- ggplot() +
  geom_raster(data = chm_df, aes(x = x, y = y, fill = Z)) +
  scale_fill_gradientn(name = "Canopy Height Model",
                       colours = c("black", "white")) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal") +
  coord_sf(crs = utm_crs)


# Create a plot for sum_mean_bai_5 with UTM 10N coordinates
plot_bai <- ggplot() +
  geom_sf(data = grids, aes(fill = cut(sum_mean_bai_5, breaks = bai_breaks)), color = "black", size = 0.1) +
  scale_fill_manual(name = "Cumulative Mean Basal\nArea of Last 5 Years (mm/yr)",
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

# Arrange the two plots side-by-side
grid.arrange(plot_bai, plot_cindex, ncol = 2)

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
grid.arrange(plot_bai, plot_cindex, plot_index, ncol = 3)















st_write(grid_cv_hex, glue::glue('{proj_dir}/output/vector/{acq}_heygi_bai_ashape_hexgrid_{cellsize}m.gpkg'), append = FALSE)

st_write(grid_cv, glue::glue('{proj_dir}/output/vector/{acq}_heygi_bai_ashape_sqgrid_{cellsize}m.gpkg'), append = FALSE)
x <- ttops_to_grid(proj_dir = blocks_dir[1])
st_write(x, 'H:/Quesnel_2022/scratch/CT1_heygi_ashape_grid.gpkg')
ttops <- st_read('H:/Quesnel_2022/process/CT1-T-DAP/output/crowns/ashapes/DAP22_CT1-T_heygi_ashapes.gpkg')
ttops <- st_read('H:/Quesnel_2022/process/CT1/output/crowns/ashapes/ULS22_CT1_heygi_ashapes.gpkg')
ttops <- ttops %>% filter(vol_concave < 800) %>% filter(n_points > 100) %>% filter(n_points < 100000)

ttops_attr <- ttops %>% mutate(mean_bai_5 = (369.924 + 3.676 * vol_concave), mean_bai_10 = (450.285 + 3.803 * vol_concave),
                          sum_bai_5 = (4289.062 + 43.452 * vol_concave), sum_bai_10 = (9727.30 + 81.75 * vol_concave))


hist(ttops_attr$n_points)

x <- ttops_attr %>% select(mean_bai_5)
ttops_attr2 <- st_read('H:/Quesnel_2022/process/CT1-T-DAP/output/crowns/CT1-T-DAP_ttops_attr_incomplete.gpkg')
st_write(ttops_attr, 'H:/Quesnel_2022/process/CT1-T-DAP/output/crowns/ashapes/CT1-T-DAP_ttops_attr_incomplete_2.gpkg', append = FALSE)

ttops_attrx <- rbind(ttops_attr, ttops_attr2)
