# Summarize point data on grid

library(sf)
library(tidyverse)
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

# Adjust Volume Outliers based on quantiles
# Load the necessary libraries
library(dplyr)

# Assuming your dataset is in a data frame called `tree_data`
# Set the desired percentiles
lower_percentile <- 0.02  # 1st percentile
upper_percentile <- 0.98  # 99th percentile

# Calculate the percentile-based cutoffs
vol_concave_cutoffs <- quantile(ttops$vol_concave, c(lower_percentile, upper_percentile))
vol_convex_cutoffs <- quantile(ttops$vol_convex, c(lower_percentile, upper_percentile))
n_points_cutoffs <- quantile(ttops$n_points, c(lower_percentile, upper_percentile))


# Visualize potential outliers

# Replace values outside the cutoffs with the cutoff values
tree_data_adjusted <- ttops %>%
  mutate(vol_concave = ifelse(vol_concave < vol_concave_cutoffs[1], vol_concave_cutoffs[1],
                              ifelse(vol_concave > vol_concave_cutoffs[2], vol_concave_cutoffs[2], vol_concave)),
         vol_convex = ifelse(vol_convex < vol_convex_cutoffs[1], vol_convex_cutoffs[1],
                             ifelse(vol_convex > vol_convex_cutoffs[2], vol_convex_cutoffs[2], vol_convex)),
         n_points = ifelse(n_points < n_points_cutoffs[1], n_points_cutoffs[1],
                           ifelse(n_points > n_points_cutoffs[2], n_points_cutoff)))



# Calculate Heygi index for each tree
print(glue::glue('Calculating Heygi index for {nrow(ttops)} tree tops'))
ttops <- silvtools::heygi_cindex(ttops, comp_input = 'vol_concave', maxR = 6)

ttops <- ttops %>% mutate(cindex = ifelse(cindex > 40, 40, cindex)) %>% mutate(vol_concave = ifelse(vol_concave > 1000, 1000, vol_concave))

# Apply linear bai models

ttops <- ttops %>% mutate(mean_bai_5 = (369.924 + 3.676 * vol_concave), mean_bai_10 = (450.285 + 3.803 * vol_concave),
                          sum_bai_5 = (4289.062 + 43.452 * vol_concave), sum_bai_10 = (9727.30 + 81.75 * vol_concave))

# Create a square grid with 100m2 cells
square_grid <- st_make_grid(ttops, cellsize = cellsize, square = TRUE) %>% st_as_sf()
# Create a hexagonal grid with 100m2 cells
hex_grid <- st_make_grid(ttops, cellsize = cellsize, square = FALSE) %>% st_as_sf()

print(glue::glue('Intersecting {cellsize}x{cellsize}m grid with tree tops'))

plot(st_geometry(ttops))
plot(st_geometry(square_grid), border = 'red', add = T)

# Add grid ID of cell intersecting with each tree top
contained  <- st_contains(square_grid, ttops) %>%
  as.data.frame() %>%
  rename(cell_id = row.id, point_id = col.id) %>%
  merge(ttops, . , by = 'point_id')

# Generate summary metrics for each grid cell
print(glue::glue('Generating summary metrics for ashapes on {cellsize}x{cellsize}m'))
stats <- contained %>% group_by(cell_id) %>% summarise(mean_cv = mean(vol_concave), sd_cv = sd(vol_concave), n = n(),
                                                       mean_ci = mean(cindex), sd_ci = sd(cindex),
                                                       mean_bai_5 = mean(mean_bai_5), sd_bai_5 = sd(mean_bai_5),
                                                       sum_mean_bai_5 = sum(mean_bai_5), sum_mean_bai_10 = sum(mean_bai_10),
                                                       mean_bai_10 = mean(mean_bai_10), sd_bai_10 = sd(mean_bai_10), n_trees = n())
# Join summary metrics to original grid
grid_cv <- st_join(square_grid, stats) #%>% filter(!is.na(cell_id))

# Hexagons

# Add grid ID of cell intersecting with each tree top
contained  <- st_contains(hex_grid, ttops) %>%
  as.data.frame() %>%
  rename(cell_id = row.id, point_id = col.id) %>%
  merge(ttops, . , by = 'point_id')

# Generate summary metrics for each grid cell
print(glue::glue('Generating summary metrics for ashapes on {cellsize}x{cellsize}m'))

stats <- contained %>% group_by(cell_id) %>% summarise(mean_cv = mean(vol_concave), sd_cv = sd(vol_concave), n = n(),
                                                       mean_ci = mean(cindex), sd_ci = sd(cindex),
                                                       mean_bai_5 = mean(mean_bai_5), sd_bai_5 = sd(mean_bai_5),
                                                       sum_mean_bai_5 = sum(mean_bai_5), sum_mean_bai_10 = sum(mean_bai_10),
                                                       mean_bai_10 = mean(mean_bai_10), sd_bai_10 = sd(mean_bai_10)) %>% ungroup()
# Join summary metrics to original grid
grid_cv_hex <- st_join(hex_grid, stats) %>% filter(!is.na(cell_id))


grids <- list(grid_cv, grid_cv_hex)
return(grids)

}

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
