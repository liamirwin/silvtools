# Solar Simulator for CT1P1, CT1P2, CT2P2
library(tidyverse)
library(silvtools)
library(sf)
library(lidR)
blocks <- list.dirs('H:/Quesnel_2022/blocks', recursive = FALSE)
blocks <- blocks[c(1,4,7,8,9)]
plot_bufs <- st_read('H:/Quesnel_2022/shp/ct_plot_20m_buffer.shp') #%>% filter(PlotID %in% plots) %>% st_buffer(dist = 10)


start_date <- as.POSIXct("2022-05-15 00:00:00", tz = "America/Los_Angeles")
end_date <- as.POSIXct("2022-09-15 00:00:00", tz = "America/Los_Angeles")
interval <- '10 min'
lat <- 53.371759251761766
lon <- -122.76808019195623

time_df <- get_solar_pos(start_date, end_date, interval, lat, lon)

solar_pos <- time_df %>% filter(wday == 2 & as.numeric(hour) %in% 7:19)


mean_irradiance <- function(proj_dir, sitename){
  irradiance_dir <- paste0(proj_dir, '/output/raster/irradiance/', sitename)
  irradiance_files <- list.files(irradiance_dir, pattern = '.tif$', full.names = T)
  print(glue::glue('Found {length(irradiance_files)} irradiance rasters in folder for {sitename}'))
  r <- terra::rast(irradiance_files)
  print(glue::glue('Calculating mean value across {length(irradiance_files)} rasters'))
  r_mean <- terra::app(r, mean)
  terra::writeRaster(r_mean, paste0(proj_dir, '/output/raster/irradiance/', sitename, '_mean_irradiance.tif'), overwrite = T)
}

mean_irradiance(proj_dir = 'H:/Quesnel_2022/blocks/CT1', sitename = 'CT1P1')
mean_irradiance(proj_dir = 'H:/Quesnel_2022/blocks/CT1', sitename = 'CT1P2')
mean_irradiance(proj_dir = 'H:/Quesnel_2022/blocks/CT2', sitename = 'CT2P1')
mean_irradiance(proj_dir = 'H:/Quesnel_2022/blocks/CT3', sitename = 'CT2P1')


map_dir <- 'F:/Quesnel_2022/Quesnel_2022_PosTex/georeferenced_adjusted'
stem_maps <- list.files(map_dir, pattern = 'shp$', full.names = T) %>%
  map_df(~sf::st_read(., quiet = TRUE)) %>% mutate(CC = replace(CC, CC == 'C/D', 'C')) %>%
  mutate(CC = factor(CC, levels = c('I','C','D'), labels = c('Intermediate', 'Co-dominant', 'Dominant'),
                     ordered = TRUE)) %>% mutate(TreeNum = str_pad(TreeNum, 3, pad = '0'),
                                                 tree_id = paste0(PlotID, '-', Species, TreeNum)) %>% relocate(tree_id, .before = PlotID)

cored_trees <- stem_maps %>% filter(Diametr > 0)
plots <- unique(cored_trees$PlotID)
plot_sf <- st_read('F:/Quesnel_2022/trimble_corrected/05_SHP/quesnel_2022_plots_utm10n.shp')
plots <- plots[6:9]
for(i in 1:length(plots)){
  plot <- plots[i]
  print(plot)
  blk <- str_extract(plot, pattern = 'CT[0-9]')
  proj_dir <- glue::glue('H:/Quesnel_2022/blocks/{blk}')
  site_name <- plot
  plot_location <- plot_sf %>% filter(PlotID == plot)
  aoi <- plot_bufs %>% filter(PlotID == plot)
  dsm <- list.files(glue::glue('{proj_dir}/output/raster/dsm'), pattern = 'fill_p2r_0.1m.tif$', full.names = T)

  silvtools::solar_simulator(dsm_file = dsm,
                  proj_dir = proj_dir,
                  site_name = site_name,
                  aoi_file = aoi,
                  sun_pos = solar_pos
  )

  print(glue::glue('Successfully created solar simulation for {site_name}, averaging results'))

  mean_irradiance(proj_dir = proj_dir, sitename = site_name)

  print(glue::glue('Generating TWI for {site_name}'))

  dem_file <- list.files(glue::glue('{proj_dir}/output/raster/dtm'), pattern = '0.25m.tif$', full.names = T)
  calc_twi(dem_file = dem_file, proj_dir = proj_dir )

}







plot_las <- list.files('H:/Quesnel_2022/plot_las/ULS22/las', pattern = '_norm.laz', full.names = T)

plot_ashape_list <- list()

for(i in 1:length(plots)){
plot <- plots[i]
print(plot)
las_dir <- str_subset(plot_las, pattern = plot)
las <- readLAS(las_dir)
ttops <- locate_trees(las, lmf(ws = 2, hmin = 10))
reference <- cored_trees %>% filter(PlotID == plot)

chm <- rasterize_canopy(las, res = 0.1, p2r(0.05))
crowns <- silvtools::mcwatershed(chm = chm, treetops = ttops, th_tree = 10)()
tree_las <- merge_spatial(las, source = crowns, attribute = 'treeID')
tree_las = add_lasattribute(tree_las, name="treeID", desc="ID of a tree")

crowns_p <- crowns %>% terra::as.polygons()
blk <- str_extract(plot, pattern = 'CT[0-9]')
irr_rasts <- glue::glue('H:/Quesnel_2022/blocks/{blk}/output/raster/irradiance/{plot}_mean_irradiance.tif')
#irr_rasts <- str_subset(list.files(glue::glue('H:/Quesnel_2022/blocks/{blk}/output/raster/irradiance/'), full.names = T), pattern = '.tif$')
#irr_rasts <- terra::rast(str_subset(irr_rasts, pattern = glue::glue('{plot}_mean')))
irr_rast <- terra::rast(irr_rasts)
irr_crowns <- terra::extract(irr_rast, crowns_p, fun = 'mean')
names(irr_crowns) <- c('treeID','irr_mean')
ttops_irr <- merge(ttops, irr_crowns, by = 'treeID')

ttops <- ttops_irr

twi_rasts <- str_subset(list.files(glue::glue('H:/Quesnel_2022/blocks/{blk}/output/raster/topography/'), full.names = T), pattern = 'twi.tif$')
twi_rast <- terra::rast(twi_rasts)
twi_crowns <- terra::extract(twi_rast, crowns_p, fun = 'mean')
names(twi_crowns) <- c('treeID','twi_mean')
ttops_twi <- merge(ttops, twi_crowns, by = 'treeID')

ttops <- ttops_twi

match_table <- t_match(reference, ttops)

ashape_list <- list()

for(n in 1:length(unique(na.omit(tree_las@data$treeID)))){
  tictoc::tic()
  tree_id <- unique(na.omit(tree_las@data$treeID))[n]

  tree <- filter_poi(tree_las, treeID == tree_id)

  ashape_list[[n]] <- alphashape_metrics(tree)
  ashape_list[[n]]$treeID <- tree_id
  print(glue::glue('Generated alphashape metrics for tree {tree_id}'))
  tictoc::toc()
}

ashape_df <- do.call('rbind', ashape_list)

ashape_ttops <- left_join(match_table, ashape_df, by = 'treeID')

plot_ashape_list[[i]] <- ashape_ttops

}

plot_ashape_list

master_df <- do.call(rbind, plot_ashape_list)
write_csv(master_df, 'D:/Proposal_2022/Thinning Paper/Analyses/master_tree_df_lmf2_irr_twi.csv')

library(tidyverse)

master_df <- read.csv('D:/Proposal_2022/Thinning Paper/Analyses/master_tree_df_lmf2_irr_twi.csv')

master_df %>% filter(Diametr > 0 & distance1 < 1) %>% filter(CC %in% c('Dominant','Co-dominant')) %>% mutate(log_vol = log(vol_convex)) %>%
  lm(Diametr ~ vol_convex, data = .) %>% ggplotRegression(.) + theme_classic() + labs(y = 'Diameter at Breast Height (cm)', x = 'Convex Hull Volume (m3)')

master_df %>% filter(Diametr > 0) %>% filter(Height > 0 ) %>%  filter(status == 'TP' & distance1 < 0.5) %>% filter(CC %in% c('Dominant','Co-dominant')) %>% mutate(log_vol = log(vol_convex)) %>%
  lm(Height ~ vol_convex, data = .) %>% ggplotRegression(.)

master_df %>% ggplot(aes(x = distance1, fill = Diametr)) + geom_histogram()

master_df %>% filter(Height > 0 ) %>% ggplot(aes(x = Z.detected, y = Height)) + geom_point()

count(master_df, status)
library(sf)
cored_trees <- sf::st_read('F:/Quesnel_2022/Dendrochronology/cores_to_crowns/postex_w_cores_2022.shp')
sum(is.na(master_df$Y.detected))
master_sf <- st_as_sf(filter(master_df, !is.na(X.detected)), coords = c('X.detected','Y.detected'), remove = F) %>% mutate(PLOTID = PlotID)

master_cindex <- heygi_cindex(master_sf, comp_input = 'vol_convex', maxR = 6)

reference <- cored_trees
detected <- master_cindex %>% filter(PLOTID %in% na.omit(unique(reference$PlotID)))

cored_match <- match(cored_trees, master_cindex)

write.csv(cored_match, 'D:/Proposal_2022/Thinning Paper/Analyses/matched_cores_jan30.csv')


cored_match <- cored_match %>% filter(sum_b_5 > 0)

cored_match %>% mutate(log_vol = log(vol_convex)) %>% filter(men_b_5 > 220) %>% filter(distance1 < 1 | distance2 < 1) %>%
  lm(men_b_5 ~ vol_convex, data = .) %>% ggplotRegression(.) + theme_classic() +
  labs(y = 'Mean 5 Year Basal Area Increment (cm2)', x = 'Convex Hull Volume (m3)')

cored_match %>% mutate(log_vol = log(vol_convex)) %>% filter(men_b_5 > 220 )  %>%
  lm(men_b_5 ~ cindex, data = .) %>% ggplotRegression(.) + theme_classic() +
  labs(y = 'Mean 5 Year Basal Area Increment (cm2)', x = 'Convex Hull Volume Competition Index')

cored_match %>% mutate(log_vol = log(vol_convex)) %>% filter(men_b_5 > 220) %>%
  lm(men_b_5 ~ irr_mean, data = .) %>% ggplotRegression(.) + theme_classic() +
  labs(y = 'Mean 5 Year Basal Area Increment (cm2)', x = 'Mean Solar Irradiance (Percent of average growing day)')

cored_match %>% mutate(log_vol = log(vol_convex)) %>% filter(men_b_5 > 220 )  %>%
  lm(men_b_5 ~ twi_mean, data = .) %>% ggplotRegression(.) + theme_classic() +
  labs(y = 'Mean 5 Year Basal Area Increment (cm2)', x = 'Mean Topographic Wetness Index')

cored_match %>% mutate(log_vol = log(vol_convex)) %>% filter(men_b_5 > 220) %>%
  lm(men_b_5 ~ Z.detected, data = .) %>% ggplotRegression(.) + theme_classic() +
  labs(y = 'Mean 5 Year Basal Area Increment (cm2)', x = 'Tree top Height (m)')

cored_match %>% mutate(log_vol = log(vol_convex)) %>% filter(men_b_5 > 220) %>%
  lm(men_b_5 ~ irr_mean, data = .) %>% ggplotRegression(.) + theme_classic() +
  labs(y = 'Mean 5 Year Basal Area Increment (cm2)', x = 'Mean Solar Irradiance (Percent of average growing day)')

cored_match %>% mutate(log_vol = log(vol_convex)) %>% filter(men_b_5 > 220) %>%
  lm(men_b_5 ~ vol_convex + cindex + irr_mean, data = .) %>% summary(.)


# ---- Random Forest Attempt ----

library(randomForest)
cored_match %>% filter(sum_b_5 > 0)
# Define the predictors (features) and response variables
predictors <- cored_match %>% filter(sum_b_5 > 0) %>% select(c("CC.x",
                "cindex", "Species.x", "vol_convex", "vol_concave", "vol_a05",
                "CV_Z", "rumple", "CRR", "irr_mean", "twi_mean", "Z.detected"))

names(predictors) <- c("Crown class",
                       "Competition index (vol_convex)", "Species", "Convex hull volume", "Concave hull volume", "Alpha 0.5 hull volume",
                       "Coefficient of Variation for Z", "Rumple index", "Crown Representation Rate",
                       "Mean Solar Irradiance", "Mean Topographic Wetness Index", 'Tree Height')


response <- cored_match$men_b_5


# Fit the random forest regression model
rf_model <- randomForest(x = predictors, y = response, data = cored_match, ntree = 500)

# Print the variable importance scores
importance(rf_model)

varImpPlot(rf_model, main = 'Random Forest Model Variable Importance for Mean Basal Area Increment Growth (5 year)')
oob_error <- rf_model$err.rate[rf_model$ntree]

match <- function(reference, detected)
{
  stopifnot(is(detected, "sf"))
  stopifnot(is(reference, "sf"))

  plot_id <- na.omit(unique(reference$PlotID))

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
  dist_max[dist_max < 5] = 5
  # := operator assigns the column index_truth a value of NA if distance is dist max
  match_table[distance1 > dist_max, index_truth1 := NA]
  match_table[distance2 > dist_max, index_truth2 := NA]
  # Get IDs of closest trees
  id = match_table[, .I[which.min(distance1)], by = index_truth1]
  match_table$index_truth1 = NA_integer_
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

  # if ("Z" %in% names(dt_reference)){
  #   data.table::setnames(X, c("X.x", "Y.x", "Z.x", "X.y", "Y.y", "Z.y"), c("X", "Y", "Z", "X.detected", "Y.detected", "Z.detected"))
  # }else {
  #   data.table::setnames(X, c("X.x", "Y.x", "X.y", "Y.y", "Z"), c("X", "Y", "X.detected", "Y.detected", "Z.detected"))
  # }

  X$status = "FN"
  X$status[is.na(X$num_tree)] = "FP"
  X$status[!is.na(X$Z.detected) & !is.na(X$num_tree)] = "TP"
  class(X) <- append("TreeMatching", class(X))
  return(X)
}
