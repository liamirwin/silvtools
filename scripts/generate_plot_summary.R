# Plot Summary App
# May 2023
# Liam Irwin

library(sf)
library(dplyr)
library(tidyverse)
library(terra)
library(tmap)
library(plotly)

plots_dir <- 'G:/Quesnel_2022/plot_summary/plots'
process_dir <- 'G:/Quesnel_2022/process'
blocks_dir <- list.dirs(plots_dir, recursive = FALSE)

# Area around plot stems to crop CHM to

chm_radius <- 2.5
coords <- 26910
matches <- list()
plot_result <- TRUE
tmap_plots <- list()

# Match Reference Stems with Closest Detected Treetops

for(i in 1:length(blocks_dir)){

  proj_dir <- blocks_dir[i]
  plot_name <- basename(proj_dir)

  # Directory where CHM etc is stored
  block_name <- str_extract(plot_name, "CT\\d+")
  block_dir <- glue::glue('{process_dir}/{block_name}')
  raster_output <- glue::glue('{block_dir}/output/raster')
  vector_output <- glue::glue('{block_dir}/output/vector')

  print(glue::glue('Loaded plot directory for {plot_name}'))

  silvtools::setup_als_dirs(proj_dir)

  area_stems <- st_read(glue::glue('{dirname(file.path(plots_dir))}/postex_w_cores_2022.gpkg'), quiet = TRUE)

  plot_stems <- area_stems %>% filter(PlotID == plot_name)

  print(glue::glue('Filtered {nrow(plot_stems)} reference trees for {plot_name} from overall {nrow(area_stems)} in dataset...'))

  core_stems <- plot_stems %>% filter(!is.na(mean_bai_5))

  print(glue::glue('Of the {nrow(plot_stems)} reference trees for {plot_name}, {nrow(core_stems)} had tree cores taken...'))

  chm_file <- stringr::str_subset(list.files(glue::glue('{block_dir}/output/raster/chm'), pattern = '.tif$',full.names = T), pattern = 'fill')

  ttops_file <- stringr::str_subset(list.files(glue::glue('{block_dir}/output/vector/treetops'), pattern = '.gpkg', full.names = T), pattern = 'ashape')

  crown_file <- stringr::str_subset(list.files(glue::glue('{block_dir}/output/vector/crowns'), pattern = '.shp$', full.names = T), pattern = 'ws2')

  if(length(chm_file) != 1){

    print(glue::glue("{length(chm_file)} CHM files matched expected input in {raster_output}/chm"))

    for(k in 1:length(chm_file)){
      print(glue::glue('{k}.  {chm_file[k]}'))
    }

    repeat {
      user_choice <- as.integer(readline("Enter the number corresponding to the correct CHM file: "))
      if (user_choice >= 1 && user_choice <= length(chm_file)) {
        break
      } else {
        cat("Invalid choice. Please enter a number between 1 and", length(chm_file), ".\n")
      }
    }

    chm_file <- chm_file[user_choice]
    cat(glue::glue("\nSelected CHM file: {chm_file}\n"))
  }

  chm <- rast(chm_file)

  ttops <- st_read(ttops_file, quiet = TRUE)

  crowns <- st_read(crown_file, quiet = TRUE)

  bbox <- plot_stems %>% st_bbox() %>% st_as_sfc() %>% terra::vect() %>%
    terra::ext()

  bbox[1,] <- bbox[1,] - chm_radius
  bbox[2,] <- bbox[2,] + chm_radius
  bbox[3,] <- bbox[3,] - chm_radius
  bbox[4,] <- bbox[4,] + chm_radius

  bbox_sf <- st_as_sf(vect(bbox))
  st_crs(bbox_sf) <- st_crs(ttops)

  plot_chm <- terra::crop(chm, bbox)

  plot_ttops <- ttops %>% st_intersection(bbox_sf)

  plot_crowns <- crowns %>% st_intersection(bbox_sf)

  matches[[i]] <- tree_matching(plot_stems, plot_ttops, plot_name)

  score <- tree_matching_scores(matches[[i]]) %>% mutate(PlotID = plot_name)

  if(plot_result == TRUE){

  tmap_mode("plot")

  matched_trees <- matches[[i]] %>% filter(!is.na(X.detected) & !is.na(Y.detected) & !is.na(X) & !is.na(Y))

  create_lines <- function(matched_trees) {
    # Create an sf object from the dataframe
    matched_trees_sf <- st_as_sf(matched_trees, coords = c("X", "Y"), crs = coords)

    # Create a second sf object for detected points
    detected_points_sf <- st_as_sf(matched_trees, coords = c("X.detected", "Y.detected"), crs = coords)

    # Combine the two sf objects into one, creating lines between the reference and detected points
    lines_sf <- sapply(1:nrow(matched_trees), function(i) {

      ref <- st_as_sf(matched_trees[i,], coords = c("X", "Y"), crs = coords)
      det <- st_as_sf(matched_trees[i,], coords = c("X.detected", "Y.detected"), crs = coords)
      line <- nngeo::st_connect(ref, det, progress = FALSE)

      return(line)

    }, simplify = FALSE)

    # Convert the list to an sf object
    lines_sf <- do.call("rbind", lapply(lines_sf, st_sf))

    return(lines_sf)
  }

  # Use the function with the matched_trees dataframe
  lines_sf <- create_lines(matched_trees)

  tmap_plots[[i]] <- tm_shape(plot_chm) +
    tm_raster(palette = viridis::viridis(150), style = "cont") +
    tm_shape(plot_ttops) +
    tm_dots(col = "red", size = 'vol_concave', shape = 21, border.col = "black",
            labels = "Detected Tree Tops") + # Add legend label for the red dots
    tm_shape(plot_stems) +
    tm_dots(col = "Species", size = "Diametr", shape = 21, border.col = "black") + # Add legend label for the dots colored by Species and sized by Diametr
    tm_shape(lines_sf) +
    tm_lines(col = 'red', lwd = 1, lty = 'solid') +
    tm_layout(main.title = glue::glue("Canopy Height Model with Detected Maxima - {plot_name}"),
              legend.outside = TRUE)

  }

  score

}

matches_df <- do.call(rbind, matches)

matches_df <- matches_df %>% filter(!is.na(TreeNum))

create_scatter_plot <- function(data, x_var, y_var, title) {
  # Calculate the linear regression model
  model <- lm(data[[y_var]] ~ data[[x_var]])

  # Get the coefficients (intercept and slope) of the linear regression model
  intercept <- coef(model)[1]
  slope <- coef(model)[2]

  # Calculate the R-squared value
  r_squared <- summary(model)$r.squared

  # Create the scatter plot with the linear regression line, equation, and R-squared value
  scatter_plot <- ggplot(data, aes_string(x = x_var, y = y_var, color = title)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = "blue") +
    annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1,
             label = paste("y =", round(intercept, 2), "+", round(slope, 2), "x\n",
                           "R^2 =", round(r_squared, 2)),
             size = 4, color = "blue", family = "mono") +
    labs(title = title, x = x_var, y = y_var) +
    theme_minimal()

  return(scatter_plot)
}

core_trees <- matches_df %>% dplyr::filter(!is.na(mean_bai_5) & !is.na(vol_concave))

# Plots
create_scatter_plot(core_trees, "mean_bai_5", "vol_concave", "Species") %>% ggplotly() %>% print()
create_scatter_plot(core_trees, "mean_bai_5", "Diametr", "Species") %>% print()
create_scatter_plot(core_trees, "mean_bai_5", "Zmax", "Species") %>% print()
create_scatter_plot(core_trees, "Diametr", "Zmax", "Species") %>% print()
matches_df %>% filter(Height > 6 & Height < 22) %>% create_scatter_plot("Height", "Zmax", "Species") %>% print()
matches_df %>% filter(Height > 6 & Height < 22) %>% create_scatter_plot("Height", "Diametr", "Species") %>% print()

scores_df <- matches_df %>% filter(DC == 1 & Species %in% c('Pl','Sx','Fd','Bf') & !is.na(CC))

score_summary <- list()

for(i in 1:length(unique(scores_df$CC))){
  crown_class <- unique(scores_df$CC)[i]
  score_summary[[i]] <- scores_df %>%
    filter(CC == crown_class) %>%
    tree_matching_scores() %>%
    mutate(CC = crown_class)
}


library(gt)

score_summary_df <- do.call(rbind, score_summary)

score_summary_df <- score_summary_df %>%
  dplyr::relocate(CC, .before = 1) %>%
  dplyr::relocate(Fscore, .before = max_dist)

# Define the desired order
levels <- c("Dominant", "Co-dominant", "Intermediate", "Suppressed")

# Change the 'Crown Class' column to a factor and specify the order
score_summary_df$CC <- factor(score_summary_df$CC, levels = levels)

# Sort the dataframe by the 'Crown Class' column
score_summary_df <- score_summary_df[order(score_summary_df$CC), ]

colnames(score_summary_df) <- c("Crown Class","N reference", "N detected", "FP", "TP", "FN", "Precision", "Recall","F-Score",
                             "Max Distance (m)", "Min Distance (m)", "Standard Deviation (m)", "Mean Distance (m)" )

tbl <- gt(score_summary_df)

# Set the table title
tbl <- tbl %>%
  tab_header(
    title = "Accuracy Assessment of Tree Detection"
  )

# Format the table
tbl <- tbl %>%
  fmt_number(
    columns = c("N reference", "N detected", "FP", "TP", "FN"),
    decimals = 0
  ) %>%
  fmt_number(
    columns = c("Precision", "Recall", "Max Distance (m)", "Standard Deviation (m)", "Mean Distance (m)", "F-Score"),
    decimals = 2
  ) %>%
  fmt_scientific(
    columns = "Min Distance (m)",
    decimals = 2
  )

# Display the table
tbl
