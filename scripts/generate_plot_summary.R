# Match Postex and Lidar Stems - Assess Crown Relationships
# May 2023
# Liam Irwin

# ---- Load Packages ----

library(sf)
library(dplyr)
library(tidyverse)
library(terra)
library(tmap)
library(plotly)
library(silvtools)
library(lidR)
library(siplab)
# ---- Processing Setup ----
plots_dir <- 'G:/Quesnel_2022/plot_summary/plots'
process_dir <- 'G:/Quesnel_2022/process'
blocks_dir <- list.dirs(plots_dir, recursive = FALSE)

# Area around plot stems to crop CHM to
chm_radius <- 10
coords <- 26910
# Output lists for matching and plotting
matches <- list()
tmap_plots <- list()
# Generate tmap validation plots? - Takes some time
plot_result <- FALSE
# Should irradiance rasters be generated for each plot? - Takes some time
generate_insol = FALSE
# Extract and include solar irradiance?
extract_insol = TRUE
# Extract and include topographic wetness?
extract_twi = TRUE
# Generate Pairwise competition index?
generate_heygi <- TRUE
# Generate Area Potentially Avaliable Index
generate_apa <- TRUE

comp_input <- 'vol_concave'
# Compute with multiple maxR values
maxR <- c(2,3,4,5,6,7,8,9,10,11,12,13,14,15)
# Match Reference Stems with Closest Detected Treetops

# Loop over all directories
for(i in 1:length(blocks_dir)){

  # Set up directories and output paths
  proj_dir <- blocks_dir[i]
  plot_name <- basename(proj_dir)
  block_name <- str_extract(plot_name, "CT\\d+")
  block_dir <- glue::glue('{process_dir}/{block_name}')
  raster_output <- glue::glue('{block_dir}/output/raster')
  vector_output <- glue::glue('{block_dir}/output/vector')

  print(glue::glue('Loaded plot directory for {plot_name}'))

  # Set up ALS directories with silvtools package
  silvtools::setup_als_dirs(proj_dir)

  # Load and filter reference trees
  area_stems <- st_read(glue::glue('{dirname(file.path(plots_dir))}/postex_w_cores_2022.gpkg'), quiet = TRUE)
  plot_stems <- area_stems %>% filter(PlotID == plot_name)
  print(glue::glue('Filtered {nrow(plot_stems)} reference trees for {plot_name} from overall {nrow(area_stems)} in dataset...'))

  # Filter tree cores
  core_stems <- plot_stems %>% filter(!is.na(mean_bai_5))
  print(glue::glue('of the {nrow(plot_stems)} reference trees for {plot_name}, {nrow(core_stems)} had tree cores taken...'))

  # Load CHM, treetops and crowns files
  chm_file <- stringr::str_subset(list.files(glue::glue('{block_dir}/output/raster/chm'), pattern = '.tif$',full.names = T), pattern = 'fill')
  ttops_file <- stringr::str_subset(list.files(glue::glue('{block_dir}/output/vector/treetops'), pattern = '.gpkg', full.names = T), pattern = 'ashape')
  crown_file <- stringr::str_subset(list.files(glue::glue('{block_dir}/output/vector/crowns'), pattern = '.shp$', full.names = T), pattern = 'ws2')

  # If multiple CHM files found, allow user to select correct one
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


  # Set up bounding box and adjust to include radius around plot
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

  if(generate_insol == TRUE){

    # Setup SOlar Position Dataframe
    start_date <- as.POSIXct("2022-05-15 00:00:00", tz = "America/Los_Angeles")
    end_date <- as.POSIXct("2022-09-15 00:00:00", tz = "America/Los_Angeles")
    interval <- '1 hour'
    lat <- 53.371759251761766
    lon <- -122.76808019195623
    time_df <- get_solar_pos(start_date, end_date, interval, lat, lon)
    # Filter positions of interest
    solar_pos <- time_df %>% filter(wday == 2 & as.numeric(hour) %in% c(9,10,11,12,13,14,15,16))
    # Find the required file
    dsm_file <- list.files(glue::glue('{block_dir}/output/raster/dsm'), pattern = 'fill_p2r_0.1m.tif$', full.names = T)
    dsm <- rast(dsm_file) %>% terra::crop(., bbox)


    # Call solar simulator
    solar_simulator(dsm_file = dsm,
                    proj_dir = proj_dir,
                    site_name = plot_name,
                    sun_pos = solar_pos,
                    num_cores = 1
    )

    create_animation(glue::glue('{proj_dir}/output/raster/irradiance/{plot_name}'), sitename = plot_name)

    # Define function to calculate mean irradiance
    mean_irradiance <- function(proj_dir, sitename){
      irradiance_dir <- paste0(proj_dir, '/output/raster/irradiance/', sitename)
      irradiance_files <- list.files(irradiance_dir, pattern = '.tif$', full.names = T)
      print(glue::glue('Found {length(irradiance_files)} irradiance rasters in folder for {sitename}'))
      r <- terra::rast(irradiance_files)
      print(glue::glue('Calculating mean value across {length(irradiance_files)} rasters'))
      r_mean <- terra::app(r, mean)
      terra::writeRaster(r_mean, paste0(proj_dir, '/output/raster/irradiance/', sitename, '_mean_irradiance.tif'), overwrite = T)
    }

    # Call the mean_irradiance function
    mean_irradiance(proj_dir = proj_dir, sitename = plot_name)

  }


  # Extract Solar Irradiance
  insol_file <- stringr::str_subset(list.files(glue::glue('{proj_dir}/output/raster/irradiance'), pattern = '.tif', full.names = T), pattern = 'mean')
  # If mean irradiance raster found and extraction requested, perform extraction
  if(length(insol_file) == 0){

    print(glue::glue('No mean irradiance raster was found in {proj_dir}/output/raster/irradiance, consider rerunning with generate_insol = TRUE'))

  }


  if(extract_insol == TRUE & length(insol_file) == 1){

    print(glue::glue('Found a mean irradiance raster for {basename(proj_dir)} clipping to crown extents'))

    insol <- rast(insol_file) %>% terra::crop(., bbox)

    insol_crowns <- exactextractr::exact_extract(insol, plot_crowns, fun = 'mean', append_cols = TRUE)

    names(insol_crowns) <- c('treeID','irr_mean')

    plot_ttops <- merge(plot_ttops, insol_crowns, by = 'treeID')

  }

  # Extract Topographic Wetness
  # If TWI raster found and extraction requested, perform extraction
  twi_file <- stringr::str_subset(list.files(glue::glue('{block_dir}/output/raster/topography'), pattern = '.tif$', full.names = T), pattern = 'twi')

  if(length(twi_file) == 0){

    print(glue::glue('No topographic wetness index raster was found in {block_dir}/output/raster/topography'))

  }

  if(extract_twi == TRUE & length(twi_file) == 1){

    print(glue::glue('Found a topographic wetness index raster for {basename(block_dir)} clipping to crown extents'))

    twi <- rast(twi_file) %>% terra::crop(., bbox)

    twi_crowns <- exactextractr::exact_extract(twi, plot_crowns, fun = 'mean', append_cols = TRUE)

    names(twi_crowns) <- c('treeID','twi_mean')

    plot_ttops <- merge(plot_ttops, twi_crowns, by = 'treeID')

  }


  if (generate_heygi == TRUE) {
    # if more than one maxR value generate unique metric for each
    if (length(maxR) > 1) {
      for (i in 1:length(maxR)) {
        heygi_ttops <-
          heygi_cindex(plot_ttops, comp_input = comp_input, maxR = maxR[i]) %>%
          st_drop_geometry() %>%
          select(cindex, treeID) %>% filter(!is.na(treeID))

        plot_ttops <- merge(plot_ttops, heygi_ttops, by = 'treeID')

      }

    } else {
      heygi_ttops <-
        heygi_cindex(plot_ttops, comp_input = comp_input, maxR = maxR) %>%
        st_drop_geometry() %>%
        select(cindex, treeID) %>% filter(!is.na(treeID))

      plot_ttops <- merge(plot_ttops, heygi_ttops, by = 'treeID')
    }

  if(generate_apa == TRUE){
    apa_ttops <- apa_cindex(plot_ttops, comp_input = comp_input) %>%
      st_drop_geometry() %>%
      select(aindex, afree, treeID) %>% filter(!is.na(treeID))


  plot_ttops <- merge(plot_ttops, apa_ttops, by = 'treeID')

  }


  # Perform tree matching and score results
  matches[[i]] <- tree_matching(plot_stems, plot_ttops, plot_name)
  score <- tree_matching_scores(matches[[i]]) %>% mutate(PlotID = plot_name)

  # If plotting requested, generate plot
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
    tm_shape(plot_crowns) +
    tm_polygons(alpha = 0, border.col = "red", border.lwd = 1, border.lty = 'dashed') +
    tm_layout(main.title = glue::glue("Canopy Height Model with Detected Maxima - {plot_name}"),
              legend.outside = TRUE)

  }

  # Return the score
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
    geom_smooth(method = "lm", se = FALSE, linetype = "solid", aes(color = title)) +
    annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1,
             label = paste("y =", round(intercept, 2), "+", round(slope, 2), "x\n",
                           "R^2 =", round(r_squared, 2)),
             size = 4, color = "blue", family = "mono") +
    labs(title = title, x = x_var, y = y_var) +
    theme_minimal()

  return(scatter_plot)
}

core_trees <- matches_df %>% dplyr::filter(!is.na(mean_bai_5) & !is.na(vol_concave))

# Plotting the relationship between mean_bai_5 and vol_concave for different species
create_scatter_plot(core_trees, "sum_bai_5", "vol_concave", "Species") %>% print()

# Plotting the relationship between mean_bai_5 and irr_mean for different species
create_scatter_plot(core_trees, "mean_bai_5", "irr_mean", "Species") %>% print()

# Plotting the relationship between mean_bai_5 and twi_mean for different species
create_scatter_plot(core_trees, "mean_bai_5", "twi_mean", "Species") %>% print()

# Plotting the relationship between mean_bai_5 and vol_concave for different species with specific filtering
core_trees %>% filter(vol_concave < 500 & mean_bai_5 < 2500) %>% create_scatter_plot(., "mean_bai_5", "vol_concave", "Species") %>% print()

# Plotting the relationship between mean_bai_5 and Diametr for different species
create_scatter_plot(core_trees, "mean_bai_5", "Diametr", "Species") %>% print()

# Plotting the relationship between mean_bai_5 and Zmax for different species
create_scatter_plot(core_trees, "mean_bai_5", "Zmax", "Species") %>% print()

# Plotting the relationship between Diametr and Zmax for different species
create_scatter_plot(core_trees, "Diametr", "Zmax", "Species") %>% print()

# Plotting the relationship between Diametr and Zmax for different species
create_scatter_plot(core_trees, "vol_concave","sum_bai_5", "Species") %>% print()

# Plotting the relationship between Height and Zmax for different species with specific filtering
matches_df %>% filter(Height > 6 & Height < 22) %>% create_scatter_plot(., "Height", "Zmax", "Species") %>% print()

# Plotting the relationship between Height and Diametr for different species with specific filtering
matches_df %>% filter(Height > 6 & Height < 22) %>% create_scatter_plot(., "Height", "Diametr", "Species") %>% print()


core_trees <- core_trees %>% select(-c(Date, PlotRds, Dist1, Dist2, Dist3, R, Theta, X_postex, Y_postex, core_tree_id, tree_id.y, num_tree, treeID, point_id, Z.detected))

write.csv(core_trees, 'G:/Quesnel_2022/Modelling/core_trees.csv')
# ---- Linear Mixed Effects Model ----
# Nested structure of trees within plots
core_trees <- read.csv('G:/Quesnel_2022/Modelling/core_trees.csv')

library(lme4)
library(lmerTest)
library(performance)
# Normal to log transform BAI measurements

# Mixed effects model which incorporates variability between plots

# Model 1: Crown Volume and Transformations

core_trees %>% lm(log(sum_bai_5) ~ log(vol_concave), data = .) %>% summary()

# Plot relationships for each species...

m1 <- lmerTest::lmer(log(sum_bai_5) ~ log(vol_concave) + (1|PlotID), data = core_trees)
r2_nakagawa(m1)

m1 <- lmerTest::lmer(log(sum_bai_5) ~ log(vol_concave) + (1|PlotID), data = core_trees)
m1 <- lmerTest::lmer(log(sum_bai_5) ~ log(vol_concave) + (1|Species) + (1|PlotID), data = core_trees)
ranef(m1)
summary(m1)



# Volume + volume squared
# Uncertain which method makes more sense
m1 <- lmer(sum_bai_5 ~ log(vol_concave) + (log(vol_concave)|PlotID), data = core_trees)
summary(m1)
m2 <- lmerTest::lmer(sum_bai_5 ~ cindex + (1|PlotID), data = core_trees)
summary(m2)
m3 <- lmer(sum_bai_5 ~ vol_concave + cindex + (1|PlotID), data = core_trees)
summary(m3)
m4 <- lmer(sum_bai_5 ~ vol_concave + cindex + irr_mean + twi_mean + (1|PlotID), data = core_trees)

library(lmerTest)

m1 <- m3

# Generate predicted values
core_trees$predicted_mean_bai_5 <- predict(m1)

# Plot actual vs predicted values
ggplot(core_trees, aes(x = cindex, y = log(sum_bai_5))) +
  geom_point() +  # Plot the actual data
  geom_line(aes(y = predicted_mean_bai_5), color = "blue") +  # Plot the model predictions
  labs(x = "vol_concave", y = "mean_bai_5",
       title = "Actual vs predicted mean_bai_5 based on vol_concave") +
  theme_minimal()

# Log Transformed
ggplot(core_trees, aes(x = log(vol_concave), y = mean_bai_5)) +
  geom_point() +  # Plot the actual data
  geom_line(aes(y = predicted_mean_bai_5), color = "blue") +  # Plot the model predictions
  labs(x = "vol_concave", y = "mean_bai_5",
       title = "Actual vs predicted mean_bai_5 based on ln of vol_concave") +
  theme_minimal()

# Plot actual vs predicted values
ggplot(core_trees, aes(x = cindex, y = mean_bai_5)) +
  geom_point() +  # Plot the actual data
  geom_line(aes(y = predicted_mean_bai_5), color = "blue") +  # Plot the model predictions
  labs(x = "Competition Index", y = "Mean Basal Area Increment (5 Years)",
       title = "Actual vs predicted Mean BAI 5 based on Competition Index") +
  theme_minimal()

# Residuals
# R2
# Cross Validation - Nested
# Predicted vs Observed
# Basal Area Competition Index - Belowground - One based on Crown Volume
# Rumple of Individual Tree Crown

# Model 1: Volume of Crown, Log Volume, Volume Squared, Log Volume Squared
# Different Competition Indexes

# Marginal and conditional predictions;
# Crown Area based index
# Apex angle?
# Species Stuff- Significant predictor of BAI but we need to map it
# Species Classificaiton?
# Potential Growth Quantile Regresssion
# Modifiers based on Competition - Pretzch
# GCC (Green/Green + blue + red) BCC RCC - Mean SD for points in Crown (1st returns)


# ---- Model Assumptions ----

# Check Residuals
df <- data.frame(resid = residuals(m1), fitted = fitted(m1))
ggplot(df, aes(x = fitted, y = resid)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(x = "Fitted values", y = "Residuals", title = "Residuals vs Fitted values")

# QQ Plot
qqnorm(residuals(m1))
qqline(residuals(m1))

# Check Variance of Random Effects for each Plot
ranef_df <- as.data.frame(ranef(m1)$PlotID)
ranef_df$PlotID <- rownames(ranef_df)  # Add a column for PlotID
ggplot(ranef_df, aes(x = PlotID, y = `(Intercept)`)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(x = "PlotID", y = "Random effect", title = "Random effects for each Plot") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")  # Add a horizontal line at y = 0











# ---- Samuel Cross Validation Code ----

# variable-wise cross validation
crossv_strat =
  function (data, var, id = '.id') {
    n <- nrow(data)
    folds = data[var][[1]]
    idx <- seq_len(n)
    fold_idx <- split(idx, folds)
    fold <- function(test) {
      list(train = modelr::resample(data, setdiff(idx, test)),
           test = modelr::resample(data, test))
    }
    cols <- purrr::transpose(purrr::map(fold_idx, fold))
    tibble::as_tibble(cols)
  }

# cross validation by POPULATION
val_mods_prov = blups_mod %>%
  crossv_strat(var = 'Prov') %>%
  # fit each model with the training data
  # mutate(model_complex = purrr::map(train,
  #                                   ~lm(noquote(paste0('DBH16 ~ ', candidate_list[[1]])), data=.))) %>%
  mutate(model_simple = purrr::map(train,
                                   ~lm(noquote(paste0('DBH16 ~ ', candidate_list[[1]])), data=.))) %>%
  mutate(model_vol_zq999 = purrr::map(train,
                                      ~lm(noquote(paste0('DBH16 ~ ', candidate_list[[2]])), data=.))) %>%
  mutate(model_struct = purrr::map(train,
                                   ~lm(noquote(paste0('DBH16 ~ ', candidate_list[[3]])), data=.))) %>%
  mutate(model_spec = purrr::map(train,
                                 ~lm(noquote(paste0('DBH16 ~ ', candidate_list[[4]])), data=.))) %>%
  mutate(model_ndre1 = purrr::map(train,
                                  ~lm(noquote(paste0('DBH16 ~ ', candidate_list[[5]])), data=.))) %>%
  mutate(model_ndvi = purrr::map(train,
                                 ~lm(noquote(paste0('DBH16 ~ ', candidate_list[[6]])), data=.))) %>%
  mutate(model_cci = purrr::map(train,
                                ~lm(noquote(paste0('DBH16 ~ ', candidate_list[[7]])), data=.))) %>%
  mutate(model_nirvndre = purrr::map(train,
                                     ~lm(noquote(paste0('DBH16 ~ ', candidate_list[[8]])), data=.))) %>%
  mutate(model_sipi = purrr::map(train,
                                 ~lm(noquote(paste0('DBH16 ~ ', candidate_list[[9]])), data=.))) %>%
  mutate(model_3 = purrr::map(train,
                              ~lm(noquote(paste0('DBH16 ~ ', candidate_list[[10]])), data=.))) %>%
  mutate(model_4 = purrr::map(train,
                              ~lm(noquote(paste0('DBH16 ~ ', candidate_list[[11]])), data=.))) %>%
  mutate(model_5 = purrr::map(train,
                              ~lm(noquote(paste0('DBH16 ~ ', candidate_list[[12]])), data=.))) %>%
  mutate(model_6 = purrr::map(train,
                              ~lm(noquote(paste0('DBH16 ~ ', candidate_list[[13]])), data=.))) %>%
  mutate(model_8 = purrr::map(train,
                              ~lm(noquote(paste0('DBH16 ~ ', candidate_list[[14]])), data=.))) %>%
  mutate(model_9 = purrr::map(train,
                              ~lm(noquote(paste0('DBH16 ~ ', candidate_list[[15]])), data=.))) %>%
  # mutate(model_10 = purrr::map(train,
  #                             ~lm(noquote(paste0('DBH16 ~ ', candidate_list[[17]])), data=.))) %>%
  mutate(mod_z = purrr::map(train,
                            ~lm(noquote(paste0('HT16 ~ ', candidate_list[[16]])), data=.))) %>%
  pivot_longer(cols = c(#model_complex,
    model_simple,
    model_vol_zq999,
    model_struct,
    model_spec,
    model_ndre1,
    model_ndvi,
    model_cci,
    model_nirvndre,
    model_sipi,
    model_3,
    model_4,
    model_5,
    model_6,
    model_8,
    model_9,
    #model_10,
    mod_z),
    names_to = 'model_name', values_to = 'model') %>%
  mutate(predicted = map2(model, test, ~ augment(.x, newdata = .y))) %>%
  unnest(predicted) %>%
  dplyr::select(-train, -test, -model) %>%
  mutate(phen = if_else(model_name == 'mod_z', 'Height', 'DBH'),
         cross_v = 'Population')





# ---- Generate Treetops for each plot area ----

# Use finer window size to enable matching with all stems (eg avoid false negative especially with core trees)

chm_radius <- 10
ttops_list <- list()

# Loop over all directories
for(i in 1:length(blocks_dir)){

  # Set up directories and output paths
  proj_dir <- blocks_dir[i]
  plot_name <- basename(proj_dir)
  block_name <- str_extract(plot_name, "CT\\d+")
  block_dir <- glue::glue('{process_dir}/{block_name}')
  raster_output <- glue::glue('{block_dir}/output/raster')
  vector_output <- glue::glue('{block_dir}/output/vector')

  print(glue::glue('Loaded plot directory for {plot_name}'))

  # Set up ALS directories with silvtools package
  silvtools::setup_als_dirs(proj_dir)

  # Load CHM, treetops and crowns files
  chm_file <- stringr::str_subset(list.files(glue::glue('{block_dir}/output/raster/chm'), pattern = '.tif$',full.names = T), pattern = 'smooth')

  # If multiple CHM files found, allow user to select correct one
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

  # Load and filter reference trees
  area_stems <- st_read(glue::glue('{dirname(file.path(plots_dir))}/postex_w_cores_2022.gpkg'), quiet = TRUE)
  plot_stems <- area_stems %>% filter(PlotID == plot_name)
  print(glue::glue('Filtered {nrow(plot_stems)} reference trees for {plot_name} from overall {nrow(area_stems)} in dataset...'))

  # Filter tree cores
  core_stems <- plot_stems %>% filter(!is.na(mean_bai_5))

  # Set up bounding box and adjust to include radius around plot
  bbox <- plot_stems %>% st_bbox() %>% st_as_sfc() %>% terra::vect() %>%
    terra::ext()
  bbox[1,] <- bbox[1,] - chm_radius
  bbox[2,] <- bbox[2,] + chm_radius
  bbox[3,] <- bbox[3,] - chm_radius
  bbox[4,] <- bbox[4,] + chm_radius
  bbox_sf <- st_as_sf(vect(bbox))
  st_crs(bbox_sf) <- st_crs(ttops)

  plot_chm <- terra::crop(chm, bbox)

  ttops <- locate_trees(plot_chm, lmf(ws = 1, hmin = 5))
  ttops$PlotID <- plot_name
  plot(plot_chm)
  plot(st_geometry(ttops),add = T)

  ttops_list[[i]] <- ttops

}

ttops_ws1 <- do.call(rbind, ttops_list)
#st_write(ttops_ws1, 'G:/Quesnel_2022/plot_summary/ct_ttops_ws1.gpkg')


# Rematch with adjusted ttops - less false negatives?

# Output lists for matching and plotting
matches <- list()
tmap_plots <- list()
chm_radius <- 10

# Loop over all directories
for(i in 1:length(blocks_dir)){

  # Set up directories and output paths
  proj_dir <- blocks_dir[i]
  plot_name <- basename(proj_dir)
  block_name <- str_extract(plot_name, "CT\\d+")
  block_dir <- glue::glue('{process_dir}/{block_name}')
  raster_output <- glue::glue('{block_dir}/output/raster')
  vector_output <- glue::glue('{block_dir}/output/vector')

  print(glue::glue('Loaded plot directory for {plot_name}'))

  # Set up ALS directories with silvtools package
  silvtools::setup_als_dirs(proj_dir)

  # Load and filter reference trees
  area_stems <- st_read(glue::glue('{dirname(file.path(plots_dir))}/postex_w_cores_2022.gpkg'), quiet = TRUE)
  plot_stems <- area_stems %>% filter(PlotID == plot_name)
  print(glue::glue('Filtered {nrow(plot_stems)} reference trees for {plot_name} from overall {nrow(area_stems)} in dataset...'))

  # Filter tree cores
  core_stems <- plot_stems %>% filter(!is.na(mean_bai_5))
  print(glue::glue('Of the {nrow(plot_stems)} reference trees for {plot_name}, {nrow(core_stems)} had tree cores taken...'))

  # Load CHM, treetops and crowns files
  chm_file <- stringr::str_subset(list.files(glue::glue('{block_dir}/output/raster/chm'), pattern = '.tif$',full.names = T), pattern = 'fill')
  #ttops_file <- stringr::str_subset(list.files(glue::glue('{block_dir}/output/vector/treetops'), pattern = '.gpkg', full.names = T), pattern = 'ashape')
  ttops_file <- 'G:/Quesnel_2022/plot_summary/ct_ttops_ws1.gpkg'
  crown_file <- stringr::str_subset(list.files(glue::glue('{block_dir}/output/vector/crowns'), pattern = '.shp$', full.names = T), pattern = 'ws2')

  # If multiple CHM files found, allow user to select correct one
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


  # Set up bounding box and adjust to include radius around plot
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
  plot_crowns <- silvtools::crown_mask(chunk = plot_chm, ttops = ttops, crown_height_threshold = 0.25, vis = FALSE) %>% terra::as.polygons(.) %>% sf::st_as_sf(.) %>%
    silvtools::convert_multi_to_single_polygons(polygons = ., fill_holes = TRUE)
  ctg_norm <- catalog(glue::glue('{block_dir}/input/las/norm'))
  plot_las <- clip_roi(ctg_norm, bbox_sf)
  print('Clipped normalized las to plot extent')

  apply_treeid_to_las <- function(chunk, crowns){
    tictoc::tic()
    if('LAS' %in% class(chunk)){
      las <- chunk
    } else{
    las <- readLAS(chunk)
    }
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
    tree_las <- filter_poi(tree_las, !is.na(treeID))
    glue::glue('Merged {nrow(chunk_crowns)} crowns with {length(las@data$Z)} points in chunk')
    tictoc::toc()
    return(tree_las)
  }

  tree_las <- apply_treeid_to_las(plot_las, plot_crowns)
  print('Applied crown tree IDs to plot las')

  plot_ashapes <- get_alphashape_metrics(tree_las)
  plot_ashape_ttops <- merge(plot_ttops, plot_ashapes, by = 'treeID')
  # Perform tree matching and score results
  matches[[i]] <- tree_matching(plot_stems, plot_ashape_ttops, plot_name)
  score <- tree_matching_scores(matches[[i]]) %>% mutate(PlotID = plot_name)

  # Return the score
  score

}


matches_df <- do.call(rbind, matches)
matches_df <- matches_df %>% filter(!is.na(TreeNum))

core_trees <- matches_df %>% dplyr::filter(!is.na(mean_bai_5) & !is.na(vol_concave))

write.csv(st_drop_geometry(core_trees),'G:/Quesnel_2022/Modelling/core_trees_fix3.csv')

ct_sf <- st_as_sf(core_trees, coords = c('X','Y'),crs = 26910)
st_crs(ct_sf)
# False Negative Matches with Core Trees?

matches_df %>% filter(status == 'FN' & !is.na(mean_bai_5))

# ---- Compute Tree Detection Accuracy Scores ----

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
