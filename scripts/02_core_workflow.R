'''

Tree Core RWL files to Attributed Tree Crowns

November 15 2022

'''


#---- Packages ----

library(tidyverse)
library(dplR)
library(sf)
library(terra)

theme_set(theme_light())

#---- Load in Cores ----

proj_dir <- 'tree_rings'



#---- Deprecated code ----

csv <- list.files(paste0(proj_dir, '/rwl/csv'), pattern = '.csv', full.names = T)
rwl <- list.files(paste0(proj_dir, '/rwl'), pattern = '.rwl', full.names = T)
# Read in RWL files, merge into list

rwl_list <- list()

for(i in 1:length(rwl)){
  # Get core id without extension
  core_id <- tools::file_path_sans_ext(basename(rwl[i]))
  # Read in RWL file
  rwl_list[[i]] <- dplR::read.rwl(rwl[i], "tucson")
  # Force proper core id as rwl column
  colnames(rwl_list[[i]]) <- core_id
}

rwl_df <- dplR::combine.rwl(rwl_list)

#----

# Read RWL from CooRecorder

rwl_df <- dplR::read.rwl('F:/Quesnel_2022/Dendrochronology/tree_rings/quesnel_cores_oct.rwl', format = "tucson")


# Calculate basal area increment for each year

rwl_bai <- bai.in(rwl_df) %>%
  rownames_to_column(var = "year") %>%
  pivot_longer(cols = !year, names_to = 'core_id', values_to = 'bai') %>% filter(!is.na(bai))

# Format ring widths as long

rwl_width <- rwl_df %>% rownames_to_column(var = "year") %>%
  pivot_longer(cols = !year, names_to = 'core_id', values_to = 'ring.width') %>% filter(!is.na(ring.width))

core_id_df <- read.csv('tree_rings/core_label_guide_2023.csv') %>% select(old_code, code) %>% rename(core_id = code, tree_id = old_code)
core_id_df <- read.csv('F:/Quesnel_2022/Dendrochronology/tree_rings/quesnel_core_list_2023.csv')
rwl_width <- merge(rwl_width, core_id_df, by = 'core_id')

# Merge BAI and widths for each core, each year (long format)

core_df <- left_join(rwl_width, rwl_bai, by = c('core_id', 'year')) %>% mutate(tree_id = str_replace(tree_id, pattern = '-[A-Z]-', replacement = '')) %>%
  mutate(year = lubridate::year(as.Date(year, format = '%Y'))) %>% relocate(tree_id, .before = ring.width)

# BAI Quality Check

core_df %>% ggplot(aes(x = bai, fill = core_id)) +
  geom_histogram() + facet_wrap(~year) + theme(legend.position = 'none')

core_df %>% ggplot(aes(x = bai, y = ring.width,color = core_id)) +
  geom_point() + facet_wrap(~year) + theme(legend.position = 'none')

core_df %>% ggplot(aes(x = bai, y = ring.width,color = core_id)) +
  geom_point() + facet_wrap(~year) + theme(legend.position = 'none') + labs(x = 'Basal Area Increment (mm2)', y = 'Ring Width (mm)' )

# Summarise across last 5 and 10 years (and average both E and S cores)

last_5 <- core_df %>% filter(year >= 2017) %>% group_by(tree_id) %>%
  summarise(mean_bai_5 = mean(bai), sum_bai_5 = sum(bai),
             mean_rw_5 = mean(ring.width), sum_rw_5 = sum(ring.width))

last_10 <- core_df %>% filter(year >= 2012) %>% group_by(tree_id) %>%
  summarise(mean_bai_10 = mean(bai), sum_bai_10 = sum(bai),
            mean_rw_10 = mean(ring.width), sum_rw_10 = sum(ring.width))

core_summary <- left_join(last_5, last_10, by = 'tree_id') %>%
  mutate(tree_id_nosp = str_remove(tree_id, pattern = '[A-Z][a-z]'))



#---- Load in Stem Maps ----

# Stem map directory (adjusted postex shapefiles)

map_dir <- 'F:/Quesnel_2022/Quesnel_2022_PosTex/georeferenced_adjusted'

# Read in stem map as one big dataframe and clean columns

stem_maps <- list.files(map_dir, pattern = 'shp$', full.names = T) %>%
  map_df(~sf::st_read(., quiet = TRUE)) %>% mutate(CC = replace(CC, CC == 'C/D', 'C')) %>%
  mutate(CC = factor(CC, levels = c('I','C','D'), labels = c('Intermediate', 'Co-dominant', 'Dominant'),
                     ordered = TRUE)) %>% mutate(TreeNum = str_pad(TreeNum, 3, pad = '0'),
                                                 tree_id = paste0(PlotID, '-', Species, TreeNum)) %>% relocate(tree_id, .before = PlotID)
# Link stem map data to cores

core_stems <- stem_maps %>% filter(tree_id %in% core_df$tree_id)

core_stems <- left_join(core_stems, core_summary, by = 'tree_id')

st_write(core_stems, 'F:/Quesnel_2022/Dendrochronology/cores_to_crowns/postex_w_cores_2022_coorecorder.shp', append = FALSE )

# Generate Summary of Cores

dbh <- core_stems %>% summarise(mean_dbh = mean(Diametr), sd_dbh = sd(Diametr), min_dbh = min(Diametr), max_dbh = max(Diametr))
bai <- core_df %>% summarise(mean = mean(bai), sd = sd(bai), min = min(bai), max = max(bai))
rwi <- core_df %>% summarise(mean = mean(ring.width), sd = sd(ring.width), min = min(ring.width), max = max(ring.width))
dbh <- core_stems %>% summarise(mean_dbh = mean(Diametr), sd_dbh = sd(Diametr), min_dbh = min(Diametr), max_dbh = max(Diametr))



#---- Load tree crowns ----

crown_dir <- 'F:/Quesnel_2022/Dendrochronology/cores_to_crowns/'
crown_dir <- 'H:/Quesnel_2022/crowns/shp/'
plot_id <- 'CT1P1'

crowns_1 <- st_read(paste0(crown_dir, plot_id, '_crowns.shp')) %>%
  mutate(plot_id = plot_id, tree_id = paste0(plot_id, '-', str_pad(treeID, 3, pad = '0'))) %>% select(!treeID)

plot_id <- 'CT1P2'

crowns_2 <- st_read(paste0(crown_dir, plot_id, '_crowns.shp'))  %>%
  mutate(plot_id = plot_id, tree_id = paste0(plot_id, '-', str_pad(treeID, 3, pad = '0'))) %>% select(!treeID)

plot_id <- 'CT2P1'

crowns_3 <- st_read(paste0(crown_dir, plot_id, '_crowns.shp'))  %>%
  mutate(plot_id = plot_id, tree_id = paste0(plot_id, '-', str_pad(treeID, 3, pad = '0'))) %>% select(!treeID)


crowns <- rbind(crowns_1, crowns_2)
crowns <- rbind(crowns, crowns_3)

core_crowns <- crowns %>% filter(tree_id %in% core_summary$tree_id_nosp)

core_stems <- core_stems %>% mutate(tree_id = tree_id_nosp) %>% filter(tree_id %in% core_crowns$tree_id) %>%
  mutate(tree_id_nosp = str_remove(tree_id, pattern = '[A-Z][a-z]')) %>% st_drop_geometry()

core_trees <- right_join(core_crowns, core_stems, by = 'tree_id') %>% mutate(crown_area = as.numeric(st_area(.))) %>%
  mutate(CW = (CWmaj + CWmin)/2, CA = pi*(CW/2)^2)

st_write(core_trees, 'F:/Quesnel_2022/Dendrochronology/cores_to_crowns/measured_crowns.shp', append = FALSE)

# ---- Cored Tree Crown Quality Assessment ----

block_dir <- 'H:/Quesnel_2022/blocks'

blocks <- basename(list.dirs(block_dir, recursive = FALSE))
blocks <- blocks[2]

# Visualizing Cored Trees by Plot

for(i in 1:length(blocks)){
  # Visualze by CT Block
  block <- blocks[i]
  # Filter out core trees in block
  block_trees <- core_trees[str_detect(core_trees$PlotID, pattern = block),]
  # Rank trees by BAI and Crown Area
  block_trees <- block_trees %>% mutate(bai_rank = dense_rank(desc(mean_bai_5)),
                                        ca_rank = dense_rank(desc(crown_area)),)
  # Load CHM for block
  chm <- rast(str_subset(list.files(glue::glue('{block_dir}/{block}/output/raster/chm'),
                                     pattern = '.tif$', full.names = T), pattern = 'fill'))
  for(i in 1:nrow(block_trees)){
    # Select i tree
    tree <- block_trees[i,]
    # Crop block CHM to buffer around crown
    tree_chm <- terra::crop(chm, ext(vect(st_geometry(tree) %>% st_buffer( dist = 3))))
    # Plotting
    plot(tree_chm, col = viridis::viridis(50))
    plot(st_geometry(tree), add = T, border = 'red', col = NA)
    # Print BAI and ranking, crown area and ranking
    print(glue::glue('{tree$tree_id}-{tree$Species} Mean BAI 5 = {round(tree$mean_bai_5)}mm2 Crown Area = {round(tree$crown_area)}m2,
                     BAI Rank = {tree$bai_rank}/{nrow(block_trees)}, Crown Area = {tree$ca_rank}/{nrow(block_trees)}'))
    # User input to continue to next tree
    readline(prompt="Press [enter] to continue")
  }
  print(glue::glue('Visualized {nrow(block_trees)} crowns for {block}'))

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

core_trees%>% filter(plot_id == 'CT1P2')%>% lm(mean_bai_5 ~ crown_area + Z + Diametr, data = .) %>% ggplotRegression(.)

mess <- core_trees %>% filter(plot_id == 'CT2P1') %>%  slice_max(mean_bai_5, n = 1)
mess2 <- core_trees %>% filter(plot_id == 'CT1P2') %>%  slice_min(mean_rw_5, n = 1)
mess <- rbind(mess, mess2)
st_write(mess, 'F:/Quesnel_2022/Dendrochronology/cores_to_crowns/messy_crowns.shp', append = FALSE)

ggplotRegression(m1)

#---- Random Forest -----

library(randomForest)
library(rfviz)
set.seed(71)

ct <- core_trees  %>% st_drop_geometry() %>% select(Species, CC, crown_area, mean_bai_5, Diametr) %>% filter(crown_area > 1)

rf <-randomForest(mean_rw_5~.,data=ct, ntree=500)

rf <- rf_prep(ct, core_trees$mean_bai_5)

bcrf <- rf_viz(rf, input=TRUE, imp=TRUE, cmd=TRUE)

importance(rf)
varImpPlot(rf)

mean(predict(rf) - ct$mean_bai_5)

core_trees %>% filter(crown_area > 1) %>% ggplot(aes(x = crown_area, y = mean_bai_5)) + geom_point() + geom_smooth()
