# Connect cored trees to lidar trees and calculate competition metrics/solar simulation
# Correlate with ring data

library(sf)
library(lidR)
library(tidyverse)

cored_trees <- st_read('F:/Quesnel_2022/Dendrochronology/cores_to_crowns/postex_w_cores_2022.shp')


plot_las <- list.files('H:/Quesnel_2022/plot_las/ULS22/las', pattern = '_norm.laz', full.names = T)

plots <- unique(cored_trees$PlotID)


# Stem map directory (adjusted postex shapefiles)

map_dir <- 'F:/Quesnel_2022/Quesnel_2022_PosTex/georeferenced_adjusted'

# Read in stem map as one big dataframe and clean columns

stem_maps <- list.files(map_dir, pattern = 'shp$', full.names = T) %>%
  map_df(~sf::st_read(., quiet = TRUE)) %>% mutate(CC = replace(CC, CC == 'C/D', 'C')) %>%
  mutate(CC = factor(CC, levels = c('I','C','D'), labels = c('Intermediate', 'Co-dominant', 'Dominant'),
                     ordered = TRUE)) %>% mutate(TreeNum = str_pad(TreeNum, 3, pad = '0'),
                                                 tree_id = paste0(PlotID, '-', Species, TreeNum)) %>% relocate(tree_id, .before = PlotID)

cored_trees <- stem_maps %>% filter(Diametr > 0)
plots <- unique(cored_trees$PlotID)
plot_ashape_list <- list()

for(i in 1:length(plots))
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
  irr_rasts <- str_subset(list.files(glue::glue('H:/Quesnel_2022/blocks/{blk}/output/raster/irradiance/'), full.names = T), pattern = '.tif$')
  irr_rasts <- terra::rast(str_subset(irr_rasts, pattern = glue::glue('{plot}_mean')))
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



for(i in 1:length(plots)){
  plot <- plots[i]
  las_dir <- str_subset(plot_las, pattern = plot)
  las <- readLAS(las_dir)
  ttops <- locate_trees(las, lmf(ws = 2, hmin = 10))
  chm <- rasterize_canopy(las, res = 0.1, p2r(0.05))
  crowns <- silvtools::mcwatershed(chm = chm, treetops = ttops, th_tree = 10)()
  crowns_p <- crowns %>% terra::as.polygons()
  blk <- str_extract(plot, pattern = 'CT[0-9]')
  irr_rasts <- list.files(glue::glue('H:/Quesnel_2022/blocks/{blk}/output/raster/irradiance/'), full.names = T)
  irr_rasts <- terra::rast(str_subset(irr_rasts, pattern = glue::glue('{plot}_mean')))
  irr_crowns <- terra::extract(irr_rast, crowns_p, fun = 'mean') %>% rename(treeID = ID)
  ttops_irr <- merge(ttops, irr_crowns, by = 'treeID')
  ttops <- ttops_irr

  t_match = function(reference, detected)
  {
    stopifnot(is(detected, "sf"))
    stopifnot(is(reference, "sf"))

    plot_id <- na.omit(unique(reference$PlotID))

    reference <- reference %>%
      dplyr::rename(X_postex = X, Y_postex = Y) %>%
      dplyr::mutate(X = unlist(purrr::map(.$geometry,1)), Y = unlist(purrr::map(.$geometry,2)),
                    PLOTID = plot_id) %>% sf::st_drop_geometry()

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
    match_table <- match_table %>% mutate(index_truth1 = ifelse(distance1 > dist_max, NA, index_truth1))
    match_table <- match_table %>% mutate(index_truth2 = ifelse(distance2 > dist_max, NA, index_truth2))
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



ggplotRegression <- function (fit) {

  require(ggplot2)

  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}


plot_ashape_list

ashape_merge <- do.call('rbind', plot_ashape_list)


ashape_merge <- read.csv('D:/Proposal_2022/Thinning Paper/Analyses/master_tree_df_lmf2_irr_twi.csv')
ashape_merge <- left_join(ashape_merge, cored_trees, by = 'tree_id', suffix = c("", ".y")) %>% # merge the two data frames by the "id" column and avoid duplicate column names
  select(-ends_with(".y"))

ashape_merge %>% ggplot(aes(x = distance1 )) + geom_histogram()

ashape_merge %>% lm(men_b_5 ~ vol_convex) %>% summary()


ashape_merge %>% filter(men_b_5 > 250) %>% lm(men_b_5 ~ irr_mean, data = .) %>% ggplotRegression(.)
ashape_merge %>% filter(men_b_5 > 500 & men_b_5 < 3000) %>% lm(men_b_5 ~ irr_mean + vol_convex + Z.detected, data = .) %>% ggplotRegression(.)
ashape_merge %>% filter(men_b_5 > 500 & men_b_5 < 3000) %>% lm(men_b_5 ~ vol_convex, data = .) %>% summary()
ashape_merge %>% filter(men_b_5 > 500 & men_b_5 < 3000) %>% lm(men_b_5 ~ vol_convex, data = .) %>% ggplotRegression(.)
ashape_merge %>% filter(men_b_5 > 500 & men_b_5 < 3000) %>% lm(mn_b_10 ~ vol_convex, data = .) %>% ggplotRegression(.)
ashape_merge %>% filter(men_b_5 > 500 & men_b_5 < 3000) %>% lm(men_b_5 ~ vol_convex, data = .) %>% ggplotRegression(.)

ashape_merge %>% lm(Diametr ~ irr_mean + vol_convex + Z.detected, data = .) %>% ggplotRegression(.)
ashape_merge %>% filter(Diametr > 0) %>% nrow()
cor(ashape_merge$)

ashape_merge %>% ggplot(aes(x = vol_convex, y = Z.detected)) + geom_point()
ashape_merge %>% filter(men_b_5 > 250) %>% ggplot(aes(x = vol_convex, y = Z.detected)) + geom_point()

ashape_merge %>% filter(men_b_5 > 250) %>% lm(men_b_5 ~ twi_mean, data = .) %>% ggplotRegression(.)
ashape_merge %>% filter(Diametr > 0) %>% nrow()
ashape_merge %>% filter(men_b_5 > 250) %>% pairs(.)
### Calculate Heygi style competition index for each tree in dataset

library(spatstat)

ashape_merge

comp_input <- 'vol_convex'

ttops <- ashape_merge %>% filter(!is.na(X.detected) & !is.na(Y.detected)) %>% st_as_sf(coords = c('X.detected','Y.detected'), remove = FALSE) %>% filter(!st_is_empty(geometry))


heygi_cindex <- function(ttops, comp_input = 'vol_convex', maxR = 6){

  tictoc::tic()

  # Convert ttops to ppp object
  trees_ppp <- ttops %>%
    mutate(X = X.detected, Y = Y.detected) %>%
    as.data.frame() %>%
    select(c(X, Y, comp_input, PLOTID))

  names(trees_ppp) <- c('X','Y','comp_value')

  trees_ppp <- trees_ppp %>%
    spatstat.geom::ppp(
      x = .$X,
      y = .$Y,
      window = spatstat.geom::owin(range(.$X),
                                   range(.$Y)),
      marks = .$comp_value)

  # Calculate heygi competition index for each tree (Z instead of dbh)
  heygi <- trees_ppp %>%
    siplab::pairwise(., maxR=maxR, kernel=siplab::powers_ker,
                     kerpar=list(pi=1, pj=1, pr=1, smark=1))

  # Join new cindex with original ttops

  trees_cindex <- heygi %>%
    as.data.frame() %>% select(x, y, cindex) %>% st_as_sf(coords = c('x', 'y'))

  ttops <- ttops %>% st_as_sf(coords = c('X.detected', 'Y.detected'))

  trees_cindex <- trees_cindex %>%
    st_join(., ttops)


  print(glue::glue('Calculated Heygi style competition for {nrow(trees_cindex)} trees assesing their {comp_input} within a {maxR}m radius'))

  tictoc::toc()

  return(trees_cindex)

}


att_trees <- ttops %>% mutate(ln_vol_convex = log(vol_convex))

c <- heygi_cindex(att_trees, comp_input = 'ln_vol_convex', maxR = 6)
c <- heygi_cindex(att_trees, comp_input = 'vol_convex', maxR = 6)
c %>% filter(men_b_5 > 500 & men_b_5 < 3000) %>%  lm(men_b_5 ~ cindex + vol_concave, data = .) %>% ggplotRegression(.)
c %>% filter(men_b_5 > 500 & men_b_5 < 3000) %>%  lm(men_b_5 ~ cindex + vol_convex + Z.detected, data = .) %>% ggplotRegression(.)
c %>% filter(men_b_5 > 500 & men_b_5 < 3000) %>%  lm(men_b_5 ~ Z.detected + Diametr, data = .) %>% summary()
c %>% filter(cindex < 20) %>% lm(Diametr ~ cindex, data = .) %>% ggplotRegression(.)

c %>% filter(cindex < 10) %>% ggplot(aes(x = cindex, y = men_b_5)) + geom_point()
