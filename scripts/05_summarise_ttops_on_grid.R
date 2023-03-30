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

ttops_to_grid <- function(proj_dir, cellsize = 10, square = TRUE, acq = NULL){

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
  print(p1)
  print(p2)

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

ttops_capped <- ttops %>% mutate(cindex = ifelse(cindex > 25, 25, cindex))

ttops <- ttops_capped

plot_outlier_distribution(ttops_capped, metric = 'cindex')

# Apply linear bai models

ttops <- ttops %>% mutate(mean_bai_5 = (369.924 + 3.676 * vol_concave), mean_bai_10 = (450.285 + 3.803 * vol_concave),
                          sum_bai_5 = (4289.062 + 43.452 * vol_concave), sum_bai_10 = (9727.30 + 81.75 * vol_concave))

bai_grid <- silvtools::summarize_grid(points_sf = ttops, grid_area = 1000, grid_shape = 'hexagon', summary_var = 'mean_bai_5', summary_fun = 'sum')
cindex_grid <- silvtools::summarize_grid(points_sf = ttops, grid_area = 1000, grid_shape = 'hexagon', summary_var = 'cindex', summary_fun = 'sum')

grids <- st_join(bai_grid, cindex_grid)

return(grids)

}

# Plot Grids

percentile_breaks <- function(data, percentiles) {
  sapply(percentiles, function(p) quantile(data, probs = p/100, na.rm = TRUE))
}

percentiles <- c(0, 5, 15, 25, 50, 75, 90, 95, 100)

# Calculate breaks for sum_mean_bai_5 and sum_cindex
bai_breaks <- percentile_breaks(grids$sum_mean_bai_5, percentiles)
cindex_breaks <- percentile_breaks(grids$sum_cindex, percentiles)

plot_outlier_distribution(grids, metric = 'sum_mean_bai_5', percentiles = percentiles)
plot_outlier_distribution(grids, metric = 'sum_cindex', percentiles = percentiles)

# Function to generate custom labels
custom_labels <- function(breaks) {
  sapply(seq_along(breaks)[-1], function(i) {
    paste0(format(round(breaks[i-1], 0), big.mark = ","), " - ", format(round(breaks[i], 0), big.mark = ","))
  })
}

library(sf)
library(ggplot2)
library(gridExtra)

# Define the desired CRS
utm_crs <- st_crs("+init=EPSG:26910")

# Create a plot for sum_mean_bai_5 with UTM 10N coordinates
plot_bai <- ggplot() +
  geom_sf(data = grids, aes(fill = cut(sum_mean_bai_5, breaks = bai_breaks)), color = "black", size = 0.1) +
  scale_fill_manual(name = "Cumulative Mean Basal\nArea of Last 5 Years (mm/yr)",
                    values = colorRampPalette(c("white", "darkgreen"))(length(bai_breaks) - 1),
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
