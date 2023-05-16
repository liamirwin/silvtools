# Format Postex Stem Maps and match with detected tree tops

# Include; postex formatting, stem matching, accuracy metrics


library(sf)
library(tidyverse)
library(silvtools)
library(gridExtra)

map_dir <- 'F:/Quesnel_2022/Quesnel_2022_PosTex/georeferenced_adjusted'

stem_maps <- list.files(map_dir, pattern = 'shp$', full.names = T) %>%
  map_df(~sf::st_read(., quiet = TRUE)) %>% mutate(CC = replace(CC, CC == 'C/D', 'C')) %>%
  mutate(CC = factor(CC, levels = c('S','I','C','D'), labels = c('Suppressed','Intermediate', 'Co-dominant', 'Dominant'),
                     ordered = TRUE)) %>% mutate(TreeNum = str_pad(TreeNum, 3, pad = '0'),
                                                 tree_id = paste0(PlotID, '-', Species, TreeNum)) %>%
  relocate(tree_id, .before = PlotID) %>%
  rename(DBH = Diametr) %>% filter(DBH > 0)

# Clip Tree Tops to Plot Extent

# List directories (each is one block in the experiment)
blocks_dir <- list.dirs('H:/Quesnel_2022/process', recursive = FALSE)
# Omit these blocks from processing stream
target <- c('CT1','CT2','CT3','CT4','CT5')
# Final block directories to assess
blocks_dir <- blocks_dir[basename(blocks_dir) %in% target]
# Maxima type ('ws2', 'lmfv', 'lmfauto')
lmf <- 'ws2'
# Plot buffer - distance from plot centre to extract tree tops
plot_radius <- 11.28

accuracy_scores <- list()
match_tables <- list()

for(i in 1:length(blocks_dir)){

  if(length(blocks_dir) == 1){
    i <- 1
  }

  proj_dir <- blocks_dir[i]

  # ----- Output directories -----

  raster_output <- glue::glue('{proj_dir}/output/raster')
  vector_output <- glue::glue('{proj_dir}/output/vector')

  # Acquisition name (block)
  acq <- basename(proj_dir)

  # Load tree tops for acquisition

  ttops_files <- list.files(glue::glue('{vector_output}/treetops'),
                            pattern = 'gpkg$',
                            full.names = T)

  ttops <- sf::st_read(ttops_files[str_detect(ttops_files, lmf)], quiet = TRUE)

  # Plot Centre

  plot_centre <- sf::st_read('F:/Quesnel_2022/shp/CT2022_plots.shp', quiet = TRUE) %>%
    st_transform('EPSG:26910') %>% filter(str_detect(PlotID, acq))

  print(glue::glue('Found {nrow(plot_centre)} stem map plots within {acq}'))

  # Buffer Plot Centre

  plot_buf <- plot_centre %>% sf::st_buffer(dist = plot_radius)

  # Extract Postex Stems for accuracy assessment

  plot_stem_maps <- stem_maps %>% filter(str_detect(tree_id, acq))

  print(glue::glue('{nrow(plot_stem_maps)} reference trees mapped in {nrow(plot_buf)} plots of {acq}'))

  # Summarize data by Species and CC
  stem_maps_summary <- plot_stem_maps %>%
    group_by(PlotID, CC, Species) %>%
    summarize(count = n()) %>%
    ungroup()

  # Create the bar graph
  ggplot(stem_maps_summary, aes(x = CC, y = count, fill = Species)) +
    geom_col(position = "stack") +
    labs(title = "Distribution of Crown Classes by Species",
         x = "Crown Class",
         y = "Number of Trees",
         fill = "Species") +
    theme_minimal() +
    facet_wrap(~PlotID)

  print(glue::glue('Discarded {nrow(filter(plot_stem_maps, CC == "Suppressed"))} suppressed trees from accuracy assessment'))

  plot_stem_maps <- plot_stem_maps %>% filter(CC != 'Suppressed')

  accuracy_stats <- list()

  matches <- list()

  for(n in 1:nrow(plot_buf)){

    plot <- plot_buf[n,]

    reference_trees <- plot_stem_maps %>% filter(PlotID == plot$PlotID) %>% st_intersection(plot_buf)

    # Extract treetops within plot boundaries for assessment

    detected_trees <- st_intersection(plot, ttops)

    print(glue::glue('{nrow(detected_trees)} detected maxima fell within the {plot_radius}m plot {plot$PlotID}'))

    matches[[n]] <- tree_matching(reference_trees, detected_trees, plot$PlotID)

    accuracy_stats[[n]] <- tree_matching_scores(matches[[n]]) %>% mutate(PlotID = plot$PlotID)

  }

  match_tables[[i]] <- do.call(rbind, matches)

  accuracy_scores[[i]] <- do.call(rbind, accuracy_stats)

}

match_tables <- do.call(rbind, match_tables)

accuracy_scores <- do.call(rbind, accuracy_scores)


filtered_trees <- match_tables %>% filter(Height > 0)

# Fit a linear model
lm_model <- lm(Z.detected ~ Height, data = filtered_trees)

# Calculate the R-squared value
r_squared <- summary(lm_model)$r.squared

# Create a label for the R-squared value
r_squared_label <- paste("R² =", round(r_squared, 2))


ggplot(filtered_trees, aes(x = Height, y = Z.detected)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue",  linetype = "dashed") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(title = "Relationship between Field and Lidar Tree Top Height",
       x = "Field Measured Height",
       y = "Lidar Measured Height") +
  annotate("text", x = 0.95 * max(filtered_trees$Height), y = 0.95 * max(filtered_trees$Z.detected), label = r_squared_label, hjust = 1)

# Reshape the data frame to a long format
accuracy_scores_long <- accuracy_scores %>%
  select(PlotID, Precision, Recall, Fscore) %>%
  gather(key = "Statistic", value = "Value", -PlotID)

# Create bar plots for Precision, Recall, and F-score, faceted by PlotID
accuracy_plot <- ggplot(accuracy_scores_long, aes(x = Statistic, y = Value, fill = Statistic)) +
  geom_col(position = "dodge") +
  facet_wrap(~PlotID, ncol = 3) +
  labs(title = "Accuracy Statistics by Plot ID",
       x = "Statistic",
       y = "Value",
       fill = "Statistic") +
  ylim(0, 1) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Display the plot
accuracy_plot

tree_matching = function(reference, detected, PlotID)
{
  stopifnot(is(detected, "sf"))
  stopifnot(is(reference, "sf"))

  reference <- reference %>%
    dplyr::rename(X_postex = X, Y_postex = Y) %>%
    dplyr::mutate(X = unlist(purrr::map(.$geometry,1)), Y = unlist(purrr::map(.$geometry,2)),
                  PLOTID = PlotID) %>% sf::st_drop_geometry()

  detected <- detected %>%
    dplyr::mutate(X = unlist(purrr::map(.$geometry,1)), Y = unlist(purrr::map(.$geometry,2))) %>% sf::st_drop_geometry()

  xy_truth    = reference %>% dplyr::select(X, Y)
  xy_detected = detected %>% dplyr::select(X, Y)
  x_truth     = xy_truth[,1]
  y_truth     = xy_truth[,2]
  x_detected  = xy_detected[,1]
  y_detected  = xy_detected[,2]
  z_detected  = detected$Z

  # Attribution of nearest and 2nd neareast referenced tree index for each detected tree
  tree <- SearchTrees::createTree(xy_truth)
  # Finds 2 nearest neighbours using xy values of truth/detected trees
  knn  <- SearchTrees::knnLookup(tree, newdat = xy_detected, k = 2L)
  # 1st nearest neighbour ID
  inds1 <- knn[,1]
  # 2nd nearest neighbour ID
  inds2 <- knn[,2]

  detected$PLOTID <- reference$PLOTID[inds1]

  match_table <- data.table::data.table(index_detected = 1:length(inds1),
                                        index_truth1   = inds1,
                                        index_truth2   = inds2)

  # Horizontal distance between detected tree and two truth neighbours (m)
  match_table$distance1 <- sqrt((x_truth[inds1] - x_detected)^2 + (y_truth[inds1] - y_detected)^2)
  match_table$distance2 <- sqrt((x_truth[inds2] - x_detected)^2 + (y_truth[inds2] - y_detected)^2)

  # Takes 10% of tree height as maximum distance, if this is less than 2m make 2m
  dist_max = z_detected*0.10
  dist_max[dist_max < 2] = 2
  # assign the column index_truth a value of NA if distance more than dist max
  match_table <- match_table %>% mutate(index_truth1 = ifelse(distance1 > dist_max, NA, index_truth1),
                                        index_truth2 = ifelse(distance2 > dist_max, NA, index_truth2))
  # Get IDs of closest trees
  id = match_table[, .I[which.min(distance1)], by = index_truth1]
  # Set match table truth index to NA
  match_table$index_truth1 = NA_integer_
  # Set index truth to tree value
  match_table[id$V1, index_truth1 := id$index_truth1]
  # Get IDs of second closest trees
  id = match_table[, .I[which.min(distance2)], by = index_truth2]
  match_table$index_truth2 = NA_integer_
  match_table[id$V1, index_truth2 := id$index_truth2]
  # If the same index value appears in both columns take the one where it is a shorter distance (index 1)
  match_table[index_truth2 %in% index_truth1, index_truth2 := NA_integer_]
  # Collate the two index truth ID columns into one
  match_table$index_truth = ifelse(is.na(match_table$index_truth1), match_table$index_truth2, match_table$index_truth1)

  ###

  # Create a rowID column for the ground truth trees
  reference$num_tree = 1:nrow(reference)
  # Apply that ID to the row of the detected tree matched to it
  detected$num_tree = match_table$index_truth
  detected$distance1 = match_table$distance1
  detected$distance2 = match_table$distance2
  # Convert to data table format
  dt_reference = data.table::as.data.table(reference)
  dt_detected     = data.table::as.data.table(detected)

  X = dplyr::full_join(dt_reference, dt_detected, by = "num_tree")
  X$PLOTID <- ifelse(is.na(X$PLOTID.x), X$PLOTID.y, X$PLOTID.x)
  X = dplyr::select(X, -PLOTID.x, -PLOTID.y)

  data.table::setDT(X)

  if ("Z" %in% names(dt_reference)){
    data.table::setnames(X, c("X.x", "Y.x", "Z.x", "X.y", "Y.y", "Z.y"), c("X", "Y", "Z", "X.detected", "Y.detected", "Z.detected"))
  }else {
    data.table::setnames(X, c("X.x", "Y.x", "X.y", "Y.y", "Z"), c("X", "Y", "X.detected", "Y.detected", "Z.detected"))
  }

  X$status = "FN"
  X$status[is.na(X$num_tree)] = "FP"
  X$status[!is.na(X$Z.detected) & !is.na(X$num_tree)] = "TP"
  class(X) <- append("TreeMatching", class(X))
  return(X)
}
