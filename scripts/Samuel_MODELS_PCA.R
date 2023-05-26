library(ggfortify)
library(ggbiplot)
library(numform)
library(tidyverse)
library(broom)
library(modelr)
library(sf)
library(correlation)
library(ggcorrplot)
library(lmerTest)
library(MuMIn)
library(caret)
library(tidytext)
library(yardstick)
library(lsmeans)
library(ggpmisc)
library(ggExtra)
library(flextable)
library(qpcR)
library(dplyr)


windowsFonts(Times = windowsFont("Times New Roman"))
windowsFonts(Cambria = windowsFont("Cambria"))

pols_dat = read_rds("D://Sync//_Sites//Sx_Genecology//all_metrics_no_outliers.rds") %>%
  filter(Edge != 1 & Dead != 1) %>%
  mutate(Tree = as.numeric(Tree),
         HT16 = HT16 /100) %>%
  mutate(Class = if_else(region == "BC" & class != "B", "Orchard", "Wildstand")) %>%
  as.data.frame() %>%
  dplyr::select(-geometry) %>%
  drop_na(mean.NDVI:CRR) %>%
  group_by(Site) %>%
  dplyr::mutate(total_spec = n())

pols_dat_dbh = pols_dat %>%
  filter(no_dbh != 1) %>%
  mutate(Tree = as.numeric(Tree)) %>%
  mutate(Class = if_else(region == "BC" & class != "B", "Orchard", "Wildstand")) %>%
  group_by(Site) %>%
  dplyr::mutate(total_dbh = n(),
         total_weev = sum(Weev16)) %>%
  ungroup() %>%
  mutate(Prov_f = factor(Prov),
         species_f = factor(species),
         Weev16_f = factor(Weev16))

pols_dat_dbh %>%
  dplyr::select(Site, total_spec, total_dbh, total_weev) %>%
  distinct()


pols_dat_dbh = mutate(pols_dat_dbh, Site = recode_factor(Site, "JordanRiver" = "Jordan River",
                      "Harrison" = "Harrison",
                      "Skimikin" = "Skimikin",
                      "Kalamalka" = "Kalamalka",
                      "TJC" = "Tête Jaune Cache",
                      "Whitecourt" = "Whitecourt"))

class_cols =  c("Wildstand" = "#1752A2",
                "Orchard" = "#0A8C35"
)

################################################################################
# TESTING CLIMATE CORRELATIONS
cutoff = 0.8

# sx_pops_clim = read.csv("D:\\Sync\\_Sites\\Sx_Genecology\\_Climate\\output\\Sx_CC_Seedlots_Normal_1961_1990_Y_S.csv") %>%
#   dplyr::select(Latitude:DD1040, -Longitude) %>%
#    relocate("Latitude", "AHM", "MCMT", "TD", "MAP") %>%
#   cor(use = "complete.obs", method = "pearson") %>%
#   data.frame() %>%
#   filter(abs(Latitude) < cutoff &
#            abs(AHM) < cutoff &
#            abs(MCMT) < cutoff &
#            abs(TD) < cutoff &
#            abs(MAP) < cutoff)
#
#
#
# sx_pops_clim %>%
#   cor(use = "complete.obs", method = "pearson",
#   y = dplyr::select(., "MAP", "TD", "MCMT", "AHM", "Latitude"),
#   x = .) %>%
#   ggcorrplot(hc.order = FALSE, type = "full",
#              outline.col = "white", lab = TRUE)
# %>%
#   dplyr::select(Seedlot, Latitude, TD, EMT, MAP, AHM)

################################################################################
# check that trends in reflectance are similar across sites.
# if not, perhaps the calibration was poor.
#
pols_dat_dbh %>%
  pivot_longer(cols = c(mean.R444:mean.R842), names_to = "band", values_to = "reflectance") %>%
  ggplot(aes(x = reflectance, color = Site)) +
  geom_density(size = 1.5) +
  theme_bw(base_size = 20) +
  facet_wrap(. ~ band, ncol = 5,
             scales = "free")


pols_dat_dbh %>%
  pivot_longer(cols = c(mean.R444:mean.R842), names_to = "band", values_to = "reflectance") %>%
  ggplot(aes(x = band, y = reflectance, fill = Site)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  facet_wrap(. ~ Site, ncol = 3)

pols_dat_dbh %>%
  ggplot(aes(x = Site, y = mean.NDVI, fill = Site)) +
  geom_boxplot() +
  theme_bw(base_size = 20)


################################################################################
# a function for fitting a mixed effects model

#df = filter(blups, variable == "HT16")
blup_mod = function(df) {
  # fit a random effects model
  mm = lmer(var_value ~ (1|Prov_f:Site_f) + (1|Blk_f:(Rep_f:Site_f)) + (1|Rep_f:Site_f)
            #+ (1|Rep_f:Site_f)
            + Site_f, data = df, REML = TRUE)

  # ran = ranef(mm)[2]

  mm_weev = lmer(var_value ~ (1|Prov_f:Site_f) + (1|Blk_f:(Rep_f:Site_f)) + (1|Rep_f:Site_f)
                 #+ (1|Rep_f:Site_f)
                 + Site_f + (1|Weev16_f:(Prov_f:Site_f)), data = df, REML = TRUE)

  # ran_weev = ranef(mm_weev)[3]
  #
  # ran_2 = rownames_to_column(ran$`Prov_f:Site_f`) %>%
  #   left_join(rownames_to_column(ran_weev$`Prov_f:Site_f`),
  #             by = "rowname")
  #
  # ggplot(ran_2, aes(x = `(Intercept).x`, y = `(Intercept).y`)) +
  #   geom_point() +
  #   theme_bw()

  site_means =
    fixef(mm) %>%
    as.data.frame() %>%
    dplyr::rename('start' = '.') %>%
    rownames_to_column() %>%
    mutate(rel_val = if_else(rowname == '(Intercept)', 0, start),
           rowname = if_else(rowname == '(Intercept)', 'Site_fJordan River', rowname)) %>%
    group_by(rowname) %>%
    mutate(Site = str_split(rowname, pattern = '_f')[[1]][2]) %>%
    ungroup() %>%
    dplyr::select(-rowname) %>%
    mutate(start_val = start[1],
           site_mean = start_val + rel_val)

  site_means_weev =
    fixef(mm_weev) %>%
    as.data.frame() %>%
    dplyr::rename('start' = '.') %>%
    rownames_to_column() %>%
    mutate(rel_val = if_else(rowname == '(Intercept)', 0, start),
           rowname = if_else(rowname == '(Intercept)', 'Site_fJordan River', rowname)) %>%
    group_by(rowname) %>%
    mutate(Site = str_split(rowname, pattern = '_f')[[1]][2]) %>%
    ungroup() %>%
    dplyr::select(-rowname) %>%
    mutate(start_val = start[1],
           site_mean_weev = start_val + rel_val)

  test = ranef(mm)$'Prov_f:Site_f' %>%
    rownames_to_column() %>%
    as_tibble() %>%
    separate(rowname, sep = ':', into = c('Prov', 'Site')) %>%
    dplyr::rename('BLUP' = '(Intercept)') %>%
    left_join(site_means) %>%
    mutate(TRAIT = BLUP + site_mean)

  test_weev = ranef(mm_weev)$'Prov_f:Site_f' %>%
    rownames_to_column() %>%
    as_tibble() %>%
    separate(rowname, sep = ':', into = c('Prov', 'Site')) %>%
    dplyr::rename('BLUP_weev' = '(Intercept)') %>%
    left_join(site_means_weev) %>%
    mutate(TRAIT_weev = BLUP_weev + site_mean_weev)

  test2 = left_join(test, test_weev, by = c("Prov", "Site")) %>%
    dplyr::select(!contains("."))
}

# produce BLUPs across variables and sites
blups = pols_dat_dbh %>%
  mutate(Prov_f = factor(Prov),
         species_f = factor(species),
         Rep_f = factor(Rep),
         Blk_f = factor(Blk),
         Site_f = factor(Site, ordered = FALSE),
         Weev16_f = factor(Weev16)) %>%
  mutate(vol_convex_raw = vol_convex,
         vol_concave_raw = vol_concave,
         vol_convex = log1p(vol_convex),
         vol_concave = log1p(vol_concave)) %>%
  pivot_longer(cols = c(HT16, DBH16, stem_vol, mean.NDVI:CRR),
               names_to = 'variable',
               values_to = 'var_value') %>%
  group_by(variable) %>%
  nest() %>%
  mutate(mm = map(data, blup_mod)) %>%
  # mutate(mm_weev = map(data, blup_mod_weev)) %>%
  # mutate(mm_weev = map(mm_weev, dplyr::select, BLUP_weev)) %>%
  unnest(c(mm)) %>%
  dplyr::select(-data)



# produce BLUPs across variables and sites
# # EXCLUDING WEEVIL
# blups_no_weevil = pols_dat_dbh %>%
#   filter(Weev16 != 1) %>%
#   mutate(Prov_f = as_factor(Prov),
#          species_f = as_factor(species),
#          Rep_f = as_factor(Rep),
#          Blk_f = as_factor(Blk),
#          Site_f = factor(Site, ordered = FALSE)) %>%
#   mutate(vol_convex_raw = vol_convex,
#          vol_concave_raw = vol_concave,
#          vol_convex = log1p(vol_convex),
#          vol_concave = log1p(vol_concave)) %>%
#   pivot_longer(cols = c(HT16, DBH16, stem_vol, mean.NDVI:CRR),
#                names_to = 'variable',
#                values_to = 'var_value') %>%
#   group_by(variable) %>%
#   nest() %>%
#   mutate(mm = map(data, blup_mod)) %>%
#   # mutate(mm_weev = map(data, blup_mod_weev)) %>%
#   # mutate(mm_weev = map(mm_weev, dplyr::select, BLUP_weev)) %>%
#   unnest(c(mm)) %>%
#   rename(BLUP_no_weevil = BLUP,
#          TRAIT_no_weevil = TRAIT,
#          site_mean_no_weevil = site_mean) %>%
#   dplyr::select(Site, variable, Prov, BLUP_no_weevil,
#                 TRAIT_no_weevil,
#                 site_mean_no_weevil)


by_pop = pols_dat_dbh %>%
  pivot_longer(cols = c(HT16, DBH16, mean.NDVI:CRR),
               names_to = 'variable',
               values_to = 'var_mean') %>%
  dplyr::select(Site, variable, var_mean, Euc, zone, Prov, Rep, Blk, Class, species, Weev16) %>%
  mutate(Group = case_when(variable %in% c('DBH16', 'HT16') ~ 'Field Traits',
                           str_detect(variable, 'mean.') ~ 'Spectral Indices'),
         Group = if_else(is.na(Group), 'Structural Metrics', Group)) %>%
  mutate(Group_lab = case_when(Group == 'Field Traits' ~ 'Field',
                               Group == 'Structural Metrics' ~ 'Structural',
                               Group == 'Spectral Indices' ~ 'Spectral')) %>%
  #filter(class == 'A+') %>%
  group_by(Prov, Site, variable) %>%
  summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.))) %>%
  mutate(Prov = as.factor(Prov)) %>%
  left_join(blups, by = c('Site', 'variable', 'Prov')) %>%
  #left_join(blups_no_weevil, by = c('Site', 'variable', 'Prov')) %>%
  pivot_longer(cols = c(TRAIT, BLUP, TRAIT_weev, BLUP_weev),
               names_to = 'mean_type',
               values_to = 'var_value')
# %>%
#   pivot_longer(cols = c(TRAIT_no_weevil, BLUP_no_weevil),
#                names_to = 'mean_type_no_weevil',
#                values_to = 'var_value_no_weevil')

################################################################################
# climate PCA

clim = read.csv("D:\\Sync\\_Sites\\Sx_Genecology\\_Climate\\output\\Sx_CC_Sites_ReferencePeriod_2005_2020_Y_S.csv") %>%
  #dplyr::select(Site, Latitude, TD, EMT, MAP, AHM) %>%
  # mutate(logMAP = log(MAP)) %>%
  # dplyr::select(-MAP)  %>%
  mutate(type = "Site",
         name = as.character(Site)) %>%
  dplyr::select(-Site)


sx_pops_clim = read.csv("D:\\Sync\\_Sites\\Sx_Genecology\\_Climate\\output\\Sx_CC_Seedlots_Normal_1961_1990_Y_S.csv") %>%
  #dplyr::select(Seedlot, Latitude, TD, EMT, MAP, AHM) %>%
  # mutate(logMAP = log(MAP)) %>%
  # dplyr::select(-MAP) %>%
  mutate(type = "Seedlot",
         name = as.character(Seedlot)) %>%
  dplyr::select(-Seedlot)

# sx_pops = sx_pops_clim %>%
#   filter(!is.na(name)) %>%
#   left_join(read.csv("D:\\Sync\\_Sites\\Sx_Genecology\\_Climate\\Sx_CC_Seedlots_info.csv") %>%
#               dplyr::select(-number, -donor), by = c("Seedlot" = "population")) %>%
# mutate(zone = case_when(region == "NM" | region == "AZ" ~ "US (south)",
#                         region == "ID" | region == "MT" | region == "WA" ~ "US (north)",
#                         region == "BC" ~ "British Columbia",
#                         region == "AB" ~ "Alberta",
#                         region == "YK" | region == "NWT" ~ "Canada (north)",
#                         region == "ON" ~ "Ontario"),
#        class = if_else(region == "AB", "B", class)) %>%
#   left_join(read.csv("D:\\Sync\\_Sites\\Sx_Genecology\\SxGenecology_predicted_Euc.csv") %>%
#               dplyr::select(sxProv, species, propGla, propEng, propSit),
#             by = c("Seedlot" = "sxProv"))

clim_dat = sx_pops_clim %>%
  bind_rows(clim) %>%
  dplyr::select(-X, -Longitude, -MAR)

# clim_dat_pca = clim_dat %>%
#   dplyr::select(-Site, -Seedlot)

pca_clim = prcomp(~ ., clim_dat[1:26], scale = TRUE)
pca_clim$x[,1:3]

pca_plot = bind_cols(clim_dat, pca_clim$x[,1:3])

#autoplot(pca_clim, colour = 'name', label = TRUE)
(bp1 = ggbiplot(pca_clim, choices = c(1,2), groups = pca_plot$type,
         varname.size = 3, color = "black", labels = clim_dat$name) +
  theme_bw(base_size = 20))

(bp2 = ggbiplot(pca_clim, choices = c(1,3), groups = pca_plot$type,
         varname.size = 3, color = "black", labels = clim_dat$name) +
  theme_bw(base_size = 20))

#s = summary(pca_clim)
#barplot(pca_clim$rotation[,1], main="PC 1 Loadings Plot", las=2)

# PC1: temperature
(pc1_plot = pca_clim$rotation %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  ggplot(aes(x = reorder(rowname, desc(abs(PC1))), y = PC1)) +
  geom_bar(stat = "identity") +
  labs(x = "Climate variables") +
  theme_bw(base_size = 24) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5)))

ggsave(plot = pc1_plot,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\PC1.tiff',
       device = tiff,
       width = 15,
       height = 7,
       units = 'in',
       dpi = 600,
       bg = 'white')

ggsave(plot = pc1_plot,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\PC1.jpeg',
       device = jpeg,
       width = 15,
       height = 7,
       units = 'in',
       dpi = 600,
       bg = 'white')

# PC2: growing season aridity
(pc2_plot = pca_clim$rotation %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  ggplot(aes(x = reorder(rowname, desc(abs(PC2))), y = PC2)) +
  geom_bar(stat = "identity") +
  labs(x = "Climate variables") +
  theme_bw(base_size = 24) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5)))

ggsave(plot = pc2_plot,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\PC2.tiff',
       device = tiff,
       width = 15,
       height = 10,
       units = 'in',
       dpi = 600,
       bg = 'white')


# PC3 growing season length
(pc3_plot = pca_clim$rotation %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  ggplot(aes(x = reorder(rowname, desc(abs(PC3))), y = PC3)) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 24) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5)))

ggsave(plot = pc3_plot,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\PC3.tiff',
       device = tiff,
       width = 15,
       height = 10,
       units = 'in',
       dpi = 600,
       bg = 'white')


################################################################################

sites_pc = pca_plot %>%
  filter(type == "Site") %>%
  mutate(Site = name) %>%
  dplyr::select(-type) %>%
  mutate(Site = recode_factor(Site, "JordanRiver" = "Jordan River",
                              "Harrison" = "Harrison",
                              "Skimikin" = "Skimikin",
                              "Kalamalka" = "Kalamalka",
                              "TJC" = "Tête Jaune Cache",
                              "Whitecourt" = "Whitecourt")) %>%
  dplyr::select(Site, name, PC1, PC2, PC3) %>%
  dplyr::rename(c(PC1_site = PC1,
           PC2_site = PC2,
           PC3_site = PC3))

pops_pc = pca_plot %>%
  filter(type == "Seedlot") %>%
  mutate(Prov = name) %>%
  dplyr::select(-name, -type) %>%
  dplyr::select(Prov, PC1, PC2, PC3)


(by_pop_plot = blups %>%
    left_join(sites_pc, by = c("Site" = "Site")) %>%
    left_join(pops_pc, by = c("Prov" = "Prov")) %>%
    mutate(`PC1____<warm_cold>` = PC1_site - PC1,
           `PC2____<arid_moist>` = PC2_site - PC2,
           `PC3____<long_short>` = PC3_site - PC3) %>%
    pivot_longer(cols = `PC1____<warm_cold>`:`PC3____<long_short>`,
                 names_to = "PC", values_to = "distance") )

saveRDS(by_pop_plot, "D:\\Sync\\_Sites\\Sx_Genecology\\by_pop_pc.rds")

pops_plot =
  sx_pops_clim %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  right_join(pops_pc, by = c("name" = "Prov"))

sites_plot = clim %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  right_join(sites_pc, by = "name") %>%
  dplyr::rename(c(PC1 = PC1_site,
                  PC2 = PC2_site,
                  PC3 = PC3_site))


all_plot = bind_rows(pops_plot, sites_plot)

saveRDS(all_plot, "D:\\Sync\\_Sites\\Sx_Genecology\\all_plot_pc.rds")

states = st_read("D:\\Sync\\_Sites\\Sx_Genecology\\_rnaturalearth\\ne_50m_admin_1_states_provinces.shp") %>%
  filter(admin %in% c("United States of America", "Canada"))

ocean = st_read("D:\\Sync\\_Sites\\Sx_Genecology\\_rnaturalearth\\ne_50m_ocean.shp")

rast_df1 = read_rds("D:\\Sync\\_Sites\\Sx_Genecology\\rast_df1")

# sx range
(pc1_map = ggplot() +
    geom_raster(data = rast_df1,  aes(x = x, y = y, fill = elev)) +
    scale_fill_gradient(high = "white", low = "grey20", guide = "none") +
  geom_sf(data = countries, fill = NA, color = "grey20") +
  geom_sf(data = states, fill = NA) +
  geom_sf(data = ocean, fill = "white", color = "grey20") +
  geom_sf(data = all_plot, aes(color = PC1), shape = 16, size = 3.5) +
    scale_color_viridis_c(option = "C", labels = scales::number_format(accuracy = 0.5)) +
  geom_sf(data = sites, color = "black", size = 3.5, shape = 21) +
  coord_sf(xlim = c(225000, 2800000), ylim = c(-768333, 2150000),
           crs = st_crs(3005)) +
  theme(panel.border = element_rect(fill = NA, color = "grey25", size = .5),
        legend.position = c(0.94, 0.5),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(family = "Cambria", size = 16),
        legend.title = element_text(family = "Cambria", size = 24),
        legend.margin = margin(8, 8, 14, 8),
        panel.background = element_rect(fill = "white", color = "black", size = .5),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()))

ggsave(plot = pc1_map,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\map_PC1.jpeg',
       device = tiff,
       width = 12.5,
       height = 14,
       units = 'in',
       dpi = 300,
       bg = 'white')

(pc2_map = ggplot() +
    geom_raster(data = rast_df1,  aes(x = x, y = y, fill = elev)) +
    scale_fill_gradient(high = "white", low = "grey20", guide = "none") +
    geom_sf(data = countries, fill = NA, color = "grey20") +
    geom_sf(data = states, fill = NA) +
    geom_sf(data = ocean, fill = "white", color = "grey20") +
    geom_sf(data = all_plot, aes(color = PC2), shape = 16, size = 3.5) +
    scale_color_viridis_c(option = "C", labels = scales::number_format(accuracy = 0.5)) +
    geom_sf(data = sites, color = "black", size = 3.5, shape = 21) +
    coord_sf(xlim = c(225000, 2800000), ylim = c(-768333, 2150000),
             crs = st_crs(3005)) +
    theme(panel.border = element_rect(fill = NA, color = "grey25", size = .5),
          legend.position = c(0.94, 0.5),
          legend.background = element_rect(fill = "white", color = "black"),
          legend.text = element_text(family = "Cambria", size = 16),
          legend.title = element_text(family = "Cambria", size = 24),
          legend.margin = margin(8, 8, 14, 8),
          panel.background = element_rect(fill = "white", color = "black", size = .5),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank()))

ggsave(plot = pc2_map,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\map_PC2.jpeg',
       device = tiff,
       width = 12.5,
       height = 14,
       units = 'in',
       dpi = 300,
       bg = 'white')

(pc3_map = ggplot() +
    geom_raster(data = rast_df1,  aes(x = x, y = y, fill = elev)) +
    scale_fill_gradient(high = "white", low = "grey20", guide = "none") +
    geom_sf(data = countries, fill = NA, color = "grey20") +
    geom_sf(data = states, fill = NA) +
    geom_sf(data = ocean, fill = "white", color = "grey20") +
    geom_sf(data = all_plot, aes(color = PC3), shape = 16, size = 3.5) +
    scale_color_viridis_c(option = "C", labels = scales::number_format(accuracy = 0.5)) +
    geom_sf(data = sites, color = "black", size = 3.5, shape = 21) +
    coord_sf(xlim = c(225000, 2800000), ylim = c(-768333, 2150000),
             crs = st_crs(3005)) +
    theme(panel.border = element_rect(fill = NA, color = "grey25", size = .5),
          legend.position = c(0.94, 0.5),
          legend.background = element_rect(fill = "white", color = "black"),
          legend.text = element_text(family = "Cambria", size = 16),
          legend.title = element_text(family = "Cambria", size = 24),
          legend.margin = margin(8, 8, 14, 8),
          panel.background = element_rect(fill = "white", color = "black", size = .5),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank()))

ggsave(plot = pc3_map,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\map_PC3.jpeg',
       device = tiff,
       width = 12.5,
       height = 14,
       units = 'in',
       dpi = 300,
       bg = 'white')


(by_pop = pols_dat_dbh %>%
  pivot_longer(cols = c(HT16, DBH16, mean.NDVI:CRR),
               names_to = 'variable',
               values_to = 'var_mean') %>%
  dplyr::select(Site, variable, var_mean, Euc, zone, Prov, Rep, Blk, Class, species, Weev16) %>%
  mutate(Group = case_when(variable %in% c('DBH16', 'HT16') ~ 'Field Traits',
                           str_detect(variable, 'mean.') ~ 'Spectral Indices'),
         Group = if_else(is.na(Group), 'Structural Metrics', Group)) %>%
  mutate(Group_lab = case_when(Group == 'Field Traits' ~ 'Field',
                               Group == 'Spectral Indices' ~ 'Spectral',
                               Group == 'Structural Metrics' ~ 'Structural')) %>%
  #filter(class == 'A+') %>%
  group_by(Prov, Site, variable) %>%
  summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.))) %>%
  mutate(Prov = as.factor(Prov)) %>%
  left_join(blups, by = c('Site', 'variable', 'Prov')) %>%
  #left_join(blups_no_weevil, by = c('Site', 'variable', 'Prov')) %>%
  pivot_longer(cols = c(TRAIT, BLUP, TRAIT_weev, BLUP_weev),
               names_to = 'mean_type',
               values_to = 'var_value')  %>%
  left_join(sites_pc, by = c("Site" = "Site")) %>%
  left_join(pops_pc, by = c("Prov" = "Prov")) %>%
  mutate(`PC1 (temperature)` = PC1_site - PC1,
         `PC2 (moisture)` = PC2_site - PC2,
         `PC3 (growing season)` = PC3_site - PC3) %>%
  pivot_longer(cols = `PC1 (temperature)`:`PC3 (growing season)`,
               names_to = "PC", values_to = "distance"))



by_site_corr = by_pop %>%
  mutate(sign = if_else(distance > 0, "pos", "neg")) %>%
  group_by(Site, PC, variable, sign, mean_type) %>%
  mutate(r = cor(distance, var_value, method = 'pearson',),
         r_sign = if_else(r > 0, "pos", "neg")) %>%
  mutate(cline = case_when(PC == "PC1 (temperature)" & sign == "neg" ~ "warm",
                           PC == "PC1 (temperature)" & sign == "pos" ~ "cold",
                           PC == "PC2 (moisture)" & sign == "neg" ~ "arid",
                           PC == "PC2 (moisture)" & sign == "pos" ~ "moist",
                           PC == "PC3 (growing season)" & sign == "neg" ~ "long",
                           PC == "PC3 (growing season)" & sign == "pos" ~ "short"),
         r_p = cor.test(distance, var_value, method = 'pearson')$p.value,
         cline = ordered(cline,
                                    levels = c("warm", "cold",
                                               "arid", "moist",
                                               "long", "short"))) %>%
  mutate(p_plot = case_when(r_p < 0.001 ~ '***',
                            r_p >=  0.001 & r_p < 0.01 ~ '**',
                            r_p >=  0.01 & r_p < 0.05 ~ '*',
                            r_p >=  0.05 & r_p < 0.1 ~ '•',
                            r_p >  0.1 ~ 'NS')) %>%
  mutate(r_plot = as.character(f_num(r, retain.leading.zero = FALSE, digits = 2)))

by_site_corr_sum = summarise_all(by_site_corr, funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.))) %>%
  mutate(r_abs = abs(r)) %>%
  mutate(r_abs_pos = if_else(sign == "pos", r_abs, NA_real_),
         r_abs_neg = if_else(sign == "neg", r_abs * -1, NA_real_),
         r_abs_sig = if_else(p_plot > 0.1, 0, r_abs),
         sig = if_else(p_plot > 0.1, 0, 1)) %>%
  group_by(variable, mean_type, PC) %>%
  mutate(r_mean = mean(r_abs),
         r_mean_sig = mean(r_abs_sig),
         n_sig = sum(sig)) %>%
  ungroup()


by_site_best_plot = by_site_corr %>%
  filter(variable == 'mean.NIRvNDRE' & mean_type == 'BLUP')

summary_stats = by_site_best_plot %>%
  group_by(Site, PC) %>%
  mutate(mean_d = mean(distance),
            std_d = sd(distance)) %>%
  pivot_wider(names_from = PC, values_from = c(mean_d, std_d))

by_site_best_sig = filter(by_site_best_plot, p_plot != "NS") %>%
  group_by(PC, Site, sign) %>%
  mutate(lab_dist = case_when(sign == "neg" ~ min(distance) - 1.4,
                              sign == "pos" & max(distance) < 5 ~ max(distance) - 3,
                              sign == "pos" & max(distance) >= 5 ~ max(distance) * .1))


(g5 = by_site_best_plot %>%
    ggplot(aes(x = distance, y = var_value, group = cline)) +
    geom_vline(xintercept = 0, color = "grey20") +
    geom_point(size = 4, alpha = .5,
               color = "steelblue4"
               #,aes(color = Class)
    ) +
    geom_smooth(data = by_site_best_sig, method = 'lm', color = '#B50010', se = FALSE, size = 1.2, span = 3) +
    #geom_smooth(data = filter(by_site_best_plot, p_plot == "NS"), method = 'lm', color = 'black', se = FALSE, size = 1.2, span = 3, linetype = 3) +
    #geom_line(aes(x = Euc, y = pred), color = 'red3', size = 1.2) +
    theme_bw(base_size = 15) +
    labs(x = 'PC transfer distance',
         #y = 'CCI') +
         y = expression('Near-infrared' ~ reflectance ~ of ~ vegetation ~
                          (NIR[V ^ NDRE]) ~ BLUPs)) +
    #      # expression(bold(Slope~of~the~upper~red~edge~(RE[UPPER]))))+
    scale_y_continuous(expand = expansion(add = c(.008,.012)), breaks = seq(-.03, .03, .01)) +
    scale_color_manual(values = class_cols) +
      geom_text(data = by_site_best_sig,
                aes(x = lab_dist),
                y = Inf,
                label = paste0(expression(italic('r'))),
                size = 6.5,
                parse = TRUE,
                vjust = 2.5,
                hjust = 0,
                family = "Cambria") +
      geom_text(data = by_site_best_sig,
                aes(x = lab_dist + 1),
                y = Inf,
                label = paste0('= ', by_site_best_sig$r_plot, " ", by_site_best_sig$p_plot),
                size = 6.5,
                parse = FALSE,
                vjust = 1.9,
                hjust = 0,
                family = "Cambria") +
      # geom_text(data = by_site_best_sig,
      #           aes(x = lab_dist + 4.3),
      #           y = Inf,
      #           label = by_site_best_sig$p_plot,
      #           size = 5.5,
      #           parse = FALSE,
      #           vjust = 2,
      #           hjust = 0,
      #           family = "Cambria") +
    theme(legend.position = 'none',
          strip.background = element_rect(fill = 'white'),
          axis.text = element_text(family='Cambria', size=18, color = 'black'),
          plot.background = element_rect(color = 'white', fill = NA, size = 1),
          #axis.ticks.y = element_blank(),
          axis.text.y = element_text(margin = margin(r = 10, l = -10)),
          plot.margin = margin(.4, .4, .4, .4, 'cm'),
          panel.border = element_rect(color = 'black', fill = NA, linewidth = .5),
          axis.title.y = element_text(size = 25, margin = margin(t = 0, r = 30, b = 0, l = 0), family = 'Cambria'),
          axis.title.x = element_text(size = 25, margin = margin(t = 12, r = 0, b = 0, l = 0), family = 'Cambria'),
          strip.text = element_text(size = 19, face = "bold", family = "Cambria")) +
    facet_grid(Site ~ PC,
               #ncol = 3,
               scales = "free"))


ggsave(plot = g5,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\figure_5.tiff',
       device = tiff,
       width = 14,
       height = 15,
       units = 'in',
       dpi = 600,
       bg = 'white')

ggsave(plot = g5,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\figure_5.jpeg',
       device = jpeg,
       width = 14,
       height = 14,
       units = 'in',
       dpi = 300,
       bg = 'white')

# keeps_list = readRDS('D:\\Sync\\_Sites\\Sx_Genecology\\_R_scripts\\output\\Models\\keeps.rds') %>%
#   filter(!(Parameter1 %in% c('mean.PRI', 'mean.Gcc', 'CRR', 'CV_Z'))) %>%
#   .$Parameter1 %>%
#   append('DBH16')


var_labs = c(DBH16 = 'DBH',
             #stem_vol = expression(Vol['stem']),
             HT16 = 'Height',
             mean.CCI = 'CCI',
             mean.SIPI = 'SIPI',
             mean.RE_upper = expression(RE['upper']),
             mean.NDVI = 'NDVI',
             mean.NIRvNDVI = expression(NIR['V'^'NDVI']),
             mean.NIRvNDRE = expression(NIR['V'^'NDRE']),
             mean.NIRvCCI = expression(NIR['V'^'CCI']),
             mean.RE_total = expression(RE['total']),
             mean.NDRE1 = expression(NDRE['705']),
             mean.EVI = 'EVI',
             mean.NDRE2 = expression(NDRE['717']),
             mean.NDRE3 = expression(NDRE['740']),
             mean.RE_lower = expression(RE['lower']),
             mean.PRI = 'PRI',
             mean.Gcc = 'GCC',
             Zq99 = expression(Z['Q99']),
             Z_mean = expression(mu['Z']),
             Zq999 = expression(Z['Q999']),
             vol_convex = expression(italic('')~V['convex']),
             n_points = expression(n[points]),
             vol_concave = expression(italic('')~V['concave']),
             area = 'Area', #expression(Area['crown']),
             rumple = 'Rumple',
             apex_angle = expression(theta['apex']),
             apex_sd = expression(sigma['apex']),
             CRR = 'CRR',
             CV_Z = expression(CV['Z']))


(g6 = by_site_corr_sum %>%
    filter(PC == "PC1 (temperature)"
           #PC == "PC2 (moisture)"
           #PC == "PC3 (growing season)"
          &
             mean_type == "BLUP" & p_plot != "NS") %>%
    mutate(Group_lab = ordered(Group_lab,
                               levels = c("Field", "Structural", "Spectral"))) %>%
    ggplot(aes(y = r_abs_neg, x = reorder(variable, r_mean_sig))) +
    geom_bar(aes(y = r_abs_neg, fill = Group_lab), stat = 'identity', width = .8) +
    geom_bar(aes(y = r_abs_pos, fill = Group_lab), stat = 'identity', width = .8,
             alpha = 1) +
    scale_fill_manual(values = c('#415080', '#774A70',
                                 '#A24E5C')) +
    geom_text(aes(y = r_abs_pos, label = p_plot), family = "Cambria", size = 5, hjust = -.18) +
    geom_text(aes(y = r_abs_neg, label = p_plot), family = "Cambria", size = 5, hjust = 1.18) +
    geom_hline(yintercept = 0, size = .8) +
    expand_limits(y = c(-.9, .87)) +
    coord_flip() +
    labs(y = expression('Correlation with transfer distance along PC1 (|'~italic(r)~'|)'),
         #y = expression('Correlation with transfer distance along PC2 <arid (|'~italic(r)~'|) moist>'),
         #y = expression('Correlation with transfer distance along PC3 <long (|'~italic(r)~'|) short>'),
         x = 'Variable') +
    theme_bw(base_size = 18) +
    theme(legend.position = 'none',
          axis.text.y = element_text(size = 18, color = 'black',
                                     family = "Cambria",
                                     margin = margin(r = 10, l = 10), vjust = .5),
          axis.text.x = element_text(size = 15, color = 'black', family = "Cambria"),
          panel.background = element_rect(color = 'white', fill = NA, size = 1),
          axis.ticks.y = element_blank(),
          plot.margin = margin(.4, .4, .4, .4, 'cm'),
          panel.border = element_rect(color = 'black', fill = NA, size = .5),
          axis.title.y = element_blank(),
          # element_text(size = 22, margin = margin(t = 0, r = 4, b = 0, l = 0),
          #                           family = 'Cambria'),
          axis.title.x = element_text(size = 22, margin = margin(t = 12, r = 0, b = 0, l = 0),
                                      family = 'Cambria'),
          strip.text = element_text(size = 18, face = "bold", family = "Cambria"),
          strip.background = element_blank(),
          strip.placement = "outside") +
    scale_x_discrete(labels = var_labs) +
    scale_y_continuous(breaks = seq(-.8, .8, .4), labels = c(".8",".4","0",".4",".8")) +#labels = numform::ff_num(zero = 0)) +
    facet_grid(Group_lab ~ Site, scales = 'free_y',
               space = 'free'))

ggsave(plot = g6,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\figure_6_P1.tiff',
       device = tiff,
       width = 15,
       height = 10,
       units = 'in',
       dpi = 600,
       bg = "white")

ggsave(plot = g6,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\figure_6_P1.jpeg',
       device = jpeg,
       width = 15,
       height = 10,
       units = 'in',
       dpi = 300,
       bg = "white")

#------------------------------------------------------------------------------#

(g6 = by_site_corr_sum %>%
    filter(#PC == "PC1 (temperature)"
           PC == "PC2 (moisture)"
           #PC == "PC3 (growing season)"
           &
             mean_type == "BLUP" & p_plot != "NS") %>%
    mutate(Group_lab = ordered(Group_lab,
                               levels = c("Field", "Structural", "Spectral"))) %>%
    ggplot(aes(y = r_abs_neg, x = reorder(variable, r_mean_sig))) +
    geom_bar(aes(y = r_abs_neg, fill = Group_lab), stat = 'identity', width = .8) +
    geom_bar(aes(y = r_abs_pos, fill = Group_lab), stat = 'identity', width = .8,
             alpha = 1) +
    scale_fill_manual(values = c('#415080', '#774A70',
                                 '#A24E5C')) +
    geom_text(aes(y = r_abs_pos, label = p_plot), family = "Cambria", size = 5, hjust = -.18) +
    geom_text(aes(y = r_abs_neg, label = p_plot), family = "Cambria", size = 5, hjust = 1.18) +
    geom_hline(yintercept = 0, size = .8) +
    expand_limits(y = c(-.9, .87)) +
    coord_flip() +
    labs(y = expression('Correlation with transfer distance along PC2 (|'~italic(r)~'|)'),
         x = 'Variable') +
    theme_bw(base_size = 18) +
    theme(legend.position = 'none',
          axis.text.y = element_text(size = 18, color = 'black',
                                     family = "Cambria",
                                     margin = margin(r = 10, l = 10), vjust = .5),
          axis.text.x = element_text(size = 15, color = 'black', family = "Cambria"),
          panel.background = element_rect(color = 'white', fill = NA, size = 1),
          axis.ticks.y = element_blank(),
          plot.margin = margin(.4, .4, .4, .4, 'cm'),
          panel.border = element_rect(color = 'black', fill = NA, size = .5),
          axis.title.y = element_blank(),
          # element_text(size = 22, margin = margin(t = 0, r = 4, b = 0, l = 0),
          #                           family = 'Cambria'),
          axis.title.x = element_text(size = 22, margin = margin(t = 12, r = 0, b = 0, l = 0),
                                      family = 'Cambria'),
          strip.text = element_text(size = 18, face = "bold", family = "Cambria"),
          strip.background = element_blank(),
          strip.placement = "outside") +
    scale_x_discrete(labels = var_labs) +
    scale_y_continuous(breaks = seq(-.8, .8, .4), labels = c(".8",".4","0",".4",".8")) +#labels = numform::ff_num(zero = 0)) +
    facet_grid(Group_lab ~ Site, scales = 'free_y',
               space = 'free'))

ggsave(plot = g6,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\figure_6_P2.tiff',
       device = tiff,
       width = 15,
       height = 10,
       units = 'in',
       dpi = 600,
       bg = "white")

ggsave(plot = g6,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\figure_6_P2.jpeg',
       device = jpeg,
       width = 15,
       height = 10,
       units = 'in',
       dpi = 300,
       bg = "white")

#------------------------------------------------------------------------------#

(g6 = by_site_corr_sum %>%
    filter(#PC == "PC1 (temperature)"
           #PC == "PC2 (moisture)"
           PC == "PC3 (growing season)"
           &
             mean_type == "BLUP" & p_plot != "NS") %>%
    mutate(Group_lab = ordered(Group_lab,
                               levels = c("Field", "Structural", "Spectral"))) %>%
    ggplot(aes(y = r_abs_neg, x = reorder(variable, r_mean_sig))) +
    geom_bar(aes(y = r_abs_neg, fill = Group_lab), stat = 'identity', width = .8) +
    geom_bar(aes(y = r_abs_pos, fill = Group_lab), stat = 'identity', width = .8,
             alpha = 1) +
    scale_fill_manual(values = c('#415080', '#774A70',
                                 '#A24E5C')) +
    geom_text(aes(y = r_abs_pos, label = p_plot), family = "Cambria", size = 5, hjust = -.18) +
    geom_text(aes(y = r_abs_neg, label = p_plot), family = "Cambria", size = 5, hjust = 1.18) +
    geom_hline(yintercept = 0, size = .8) +
    expand_limits(y = c(-.9, .87)) +
    coord_flip() +
    labs(y = expression('Correlation with transfer distance along PC3 (|'~italic(r)~'|)'),
         #y = expression('Correlation with transfer distance along PC2 <arid (|'~italic(r)~'|) moist>'),
         #y = expression('Correlation with transfer distance along PC3 <long (|'~italic(r)~'|) short>'),
         x = 'Variable') +
    theme_bw(base_size = 18) +
    theme(legend.position = 'none',
          axis.text.y = element_text(size = 18, color = 'black',
                                     family = "Cambria",
                                     margin = margin(r = 10, l = 10), vjust = .5),
          axis.text.x = element_text(size = 15, color = 'black', family = "Cambria"),
          panel.background = element_rect(color = 'white', fill = NA, size = 1),
          axis.ticks.y = element_blank(),
          plot.margin = margin(.4, .4, .4, .4, 'cm'),
          panel.border = element_rect(color = 'black', fill = NA, size = .5),
          axis.title.y = element_blank(),
          # element_text(size = 22, margin = margin(t = 0, r = 4, b = 0, l = 0),
          #                           family = 'Cambria'),
          axis.title.x = element_text(size = 22, margin = margin(t = 12, r = 0, b = 0, l = 0),
                                      family = 'Cambria'),
          strip.text = element_text(size = 18, face = "bold", family = "Cambria"),
          strip.background = element_blank(),
          strip.placement = "outside") +
    scale_x_discrete(labels = var_labs) +
    scale_y_continuous(breaks = seq(-.8, .8, .4), labels = c(".8",".4","0",".4",".8")) +#labels = numform::ff_num(zero = 0)) +
    facet_grid(Group_lab ~ Site, scales = 'free_y',
               space = 'free'))

ggsave(plot = g6,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\figure_6_P3.tiff',
       device = tiff,
       width = 15,
       height = 10,
       units = 'in',
       dpi = 600,
       bg = "white")

ggsave(plot = g6,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\figure_6_P3.jpeg',
       device = jpeg,
       width = 15,
       height = 10,
       units = 'in',
       dpi = 300,
       bg = "white")



################################################################################
# pivot the transfer fits wider and look at the best R^2 across sites

# keeps_list = readRDS("D:\\Sync\\_Sites\\Sx_Genecology\\_R_scripts\\output\\Models\\keeps.rds") %>%
#   filter(!(Parameter1 %in% c("mean.PRI", "mean.Gcc", "CRR", "CV_Z"))) %>%
#   .$Parameter1 %>%
#   append("DBH16")


pols_dat_dbh %>%
  dplyr::select(#HT16, DBH16, vol_per_ha, growth_ratio,
                mean.NDVI:CRR,
                Site) %>%
  ungroup() %>%
  filter(Site == "Whitecourt") %>%
  dplyr::select(-Site) %>%
  cor(use = "complete.obs", method = "pearson") %>%
  ggcorrplot(hc.order = FALSE, type = "lower",
             outline.col = "white", lab = TRUE)


# correlation among variables
# select variables
prop_var = summary(pca_clim)$importance[2,1:3] %>%
  unname()

# rankings across PC1 and PC2 for variable selection
by_site_corr_rank = by_site_corr_sum %>%
  filter(mean_type == 'BLUP') %>%
  mutate(pc_prop = case_when(PC == "PC1 (temperature)" ~ prop_var[1],
                             PC == "PC2 (moisture)" ~ prop_var[2],
                             PC == "PC3 (growing season)" ~ prop_var[3])) %>%
  group_by(variable, PC) %>%
  mutate(r_sum_sig = sum(r_abs_sig),
         r_sum_sig_weighted = r_sum_sig * pc_prop) %>%
  group_by(variable) %>%
  mutate(score = sum(r_sum_sig_weighted)) %>%
  dplyr::select(variable, score) %>%
  distinct()
#
# keep_sig = filter(by_site_corr_rank, n_sig_mean > 3 &
#                     variable != "HT16" &
#                     variable != "DBH16") %>%
#   dplyr::select(variable) %>%
#   distinct %>%
#   .$variable
#

# drop variables with very low importance
scores_plot = ggplot(by_site_corr_rank,
       aes(x = reorder(variable, score, decreasing = TRUE),
           y = score)) +
  geom_point(size = 3)  +
  theme_bw(base_size = 25) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5),
        axis.title.x = element_blank()) +
  scale_x_discrete(labels = var_labs)

ggsave(plot = scores_plot,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\variable_scores.jpeg',
       device = jpeg,
       width = 15,
       height = 10,
       units = 'in',
       dpi = 300,
       bg = "white")

################################################################################
# table for paper

(cor_table = by_site_corr_sum %>%
   filter(mean_type == "BLUP") %>%
   dplyr::select(Site:sign, r, p_plot, Group_lab, r_mean_sig) %>%
   mutate(Group_lab = ordered(Group_lab,
                              levels = c("Field", "Structural", "Spectral")),
          p_plot = na_if(p_plot, "NS"),
          r = round(r, 2)) %>%
   unite(r_tab, c("r", "p_plot"), sep = "", na.rm = TRUE) %>%
   pivot_wider(names_from = c(Site, sign), values_from = r_tab) %>%
   left_join(by_site_corr_rank, by = "variable") %>%
   group_by(variable, PC, Group_lab) %>%
   mutate(r_mean_sig_plot = mean(r_mean_sig)) %>%
   group_by(PC, Group_lab) %>%
   arrange(desc(r_mean_sig_plot), .by_group = TRUE) %>%
   ungroup() %>%
   dplyr::select(-r_mean_sig_plot, -r_mean_sig, -score) %>%
   add_row(.before = 59) %>%
   add_row(.before = 30) %>%
   add_row(.before = 1) %>%
   flextable::flextable())

################################################################################

save_as_docx(
  'Correlations by Site, PC, and sign' = cor_table,
  path = 'D:\\Sync\\Figures\\_TABLES_paper3\\paper3_table_corr_GCB.docx')


keep_score = by_site_corr_rank %>%
  filter(score > 10.1 &
           !(variable %in% c("HT16", "DBH16"))) %>%
  .$variable

cors = pols_dat_dbh %>%
  dplyr::select(#HT16, DBH16, vol_per_ha, growth_ratio,
    all_of(keep_score), Site) %>%
  group_by(Site) %>%
  nest() %>%
  mutate(cor_df = map(data, correlation, method = 'pearson', redundant = TRUE)) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(cor_df)) %>%
  left_join(by_site_corr_rank, by = c('Parameter1' = 'variable')) %>%
  dplyr::rename("score_1" = "score") %>%
  left_join(by_site_corr_rank, by = c('Parameter2' = 'variable')) %>%
  dplyr::rename('score_2' = "score") %>%
  filter(Parameter1 != Parameter2) %>%
  group_by(Parameter1) %>%
  dplyr::mutate(r_max = max(abs(r))) %>%
  ungroup()

drops = cors %>%
  dplyr::filter(abs(r) > .9 &
                  (Parameter2 == "vol_convex" |
                  Parameter1 == "vol_convex")) %>%
  distinct(Parameter1) %>%
  filter(Parameter1 != "vol_convex") %>%
  .$Parameter1

# # which is most closely correlated to DBH?
# test_dbh = filter(cors, Parameter2 == "DBH16" & Parameter1 %in% drops$Parameter1) %>%
#   dplyr::select(Site, Parameter1, r, score_1) %>%
#   group_by(Parameter1) %>%
#   summarise(mean_r = mean(r))
#   distinct()

cors_use = cors %>%
  dplyr::filter(!(Parameter1 %in% drops | Parameter2 %in% drops)) %>%
  group_by(Parameter1) %>%
  dplyr::mutate(r_max = max(abs(r))) %>%
  ungroup()

# vars which aren't highly correlated with any others
keeps1 = cors_use %>%
  filter(r_max <= .9) %>%
  ungroup() %>%
  dplyr::select(Parameter1, score_1, r_max) %>%
  distinct()


# for groups are vars that are highly correlated,
# choose the one with the best transfer fit
keeps2 = cors_use %>%
  filter(abs(r) > .9) %>%
  dplyr::select(Parameter1, Parameter2, score_1, score_2, r_max, r) %>%
  dplyr::filter(score_1 > score_2) %>%
  dplyr::mutate(test = if_else(Parameter1 %in% if_any(Parameter2), 1, 0)) %>%
  dplyr::filter(test != 1) %>%
  dplyr::distinct()

keeps = bind_rows(keeps1, keeps2) %>%
  group_by(Parameter1) %>%
  #tally() %>%
  #filter(n >= 6) %>%
  left_join(dplyr::select(cors, Parameter1, score_1)) %>%
  #dplyr::select(-Site) %>%
  distinct(Parameter1, score_1)

saveRDS(keeps, 'D:\\Sync\\_Sites\\Sx_Genecology\\_R_scripts\\output\\Models\\keeps_pca.rds')
keeps = readRDS('D:\\Sync\\_Sites\\Sx_Genecology\\_R_scripts\\output\\Models\\keeps_pca.rds')

keeps_list = keeps %>%
  # filter(!(Parameter1 %in% c('mean.PRI', 'mean.Gcc',
  #                            'CRR', 'CV_Z'))) %>%
  .$Parameter1

pols_dat_dbh %>%
  ungroup() %>%
  dplyr::select(#HT16, DBH16, vol_per_ha, growth_ratio,
    matches(keeps_list), Site) %>%
  filter(Site == 'Skimikin') %>% #'Tête Jaune Cache') %>%
  dplyr::select(-Site, -(vol_con_Zq999:growth_ratio_Zq999)) %>%
  cor(use = 'complete.obs', method = 'pearson') %>%
  ggcorrplot(hc.order = FALSE, type = 'lower',
             outline.col = 'white', lab = TRUE)



################################################################################
blups_mod = by_pop %>%
  filter(mean_type == 'BLUP') %>%
  dplyr::select(Prov:Euc, Class, species, var_value, -var_mean) %>%
  distinct() %>%
  pivot_wider(names_from = 'variable', values_from = 'var_value') %>%
  dplyr::select(matches(keeps_list),
                Site,
                DBH16, HT16,
                Euc, Prov, Class) %>%
  mutate(lag = as.factor(case_when(
    #Site == 'Jordan River' ~ 'lag_19_20',
    Site %in% c('Harrison', 'Kalamalka', 'Skimikin') ~ 'lag_0',
    Site %in% c('Tête Jaune Cache', 'Whitecourt', 'Jordan River') ~ 'lag_1'))) %>%
  distinct()


# Give lm with all variables; after filtering to colinearity

mod_all = lm(blups_mod,
             formula = DBH16 ~
               vol_convex +
               lag:vol_convex +
               apex_angle +
               CRR +
               mean.NDVI +
               mean.SIPI +
               mean.NDRE1 +
               mean.NIRvNDRE +
               Zq999 +
               mean.CCI,
             na.action = 'na.fail')

# Conv to dataframe

dd = dredge(mod_all) %>%
  tibble::rownames_to_column()

# Gives AIC orders models

### best performing model
# dd_complex = dd %>%
#   subset(delta == 0)
# dd_best_complex = MuMIn::get.models(dd, subset = 1)[[1]]

### simple
dd_simple = dd %>%
  subset((MuMIn::has('vol_convex', 'lag:vol_convex'))) %>%
  subset(df == 4)

# Best model with 3, 4, 5; compare models with similar levels of complexity(DF can be used to subset models with less varibles)

dd_best_simple = MuMIn::get.models(dd_simple, subset = 1)[[1]]

dd_vol_zq999 = dd %>%
  subset((has('vol_convex', 'lag:vol_convex', 'Zq999'))) %>%
  subset(df == 5)

dd_best_vol_zq999 = MuMIn::get.models(dd_vol_zq999, subset = 1)[[1]]

### structural
dd_struct = dd %>%
  subset(has(!'mean.NDRE1', !'mean.NDVI', !'mean.NIRvNDRE', !'mean.CCI', !'mean.SIPI'))

dd_best_struct = MuMIn::get.models(dd_struct, subset = 1)[[1]]

### vol + spectral
dd_spec = dd %>%
  subset(has(!'apex_angle', !'CRR', !'Zq999', 'vol_convex', 'lag:vol_convex'))

dd_best_spec = MuMIn::get.models(dd_spec, subset = 1)[[1]]

### vol + ndre1
dd_vol_ndre1 = dd %>%
  subset((has('vol_convex', 'lag:vol_convex'))) %>%
  subset(has('mean.NDRE1')) %>%
  subset(df == 5)

dd_best_vol_ndre1 = MuMIn::get.models(dd_vol_ndre1, subset = 1)[[1]]

### vol + ndvi
dd_vol_ndvi = dd %>%
  subset((has('vol_convex', 'lag:vol_convex'))) %>%
  subset(has('mean.NDVI')) %>%
  subset(df == 5)

dd_best_vol_ndvi = MuMIn::get.models(dd_vol_ndvi, subset = 1)[[1]]

### vol + NIRvNDRE
dd_vol_nirvndre = dd %>%
  subset((has('vol_convex', 'lag:vol_convex'))) %>%
  subset(has('mean.NIRvNDRE')) %>%
  subset(df == 5)

dd_best_vol_nirvndre = MuMIn::get.models(dd_vol_nirvndre, subset = 1)[[1]]

# ### vol + RE_lower
# dd_vol_re_lower = dd %>%
#   subset((has('vol_convex', 'lag:vol_convex'))) %>%
#   subset(has('mean.RE_lower')) %>%
#   subset(df == 5)
#
# dd_best_vol_re_lower = MuMIn::get.models(dd_vol_re_lower, subset = 1)[[1]]

### vol + NIRvCCI
dd_vol_cci = dd %>%
  subset((has('vol_convex', 'lag:vol_convex'))) %>%
  subset(has('mean.CCI')) %>%
  subset(df == 5)

dd_best_vol_cci = MuMIn::get.models(dd_vol_cci, subset = 1)[[1]]

### vol + SIPI
dd_vol_sipi = dd %>%
  subset((has('vol_convex', 'lag:vol_convex'))) %>%
  subset(has('mean.SIPI')) %>%
  subset(df == 5)

dd_best_vol_sipi = MuMIn::get.models(dd_vol_sipi, subset = 1)[[1]]


### simple
dd_3 = dd %>%
  subset((has('vol_convex', 'lag:vol_convex'))) %>%
  subset(df == 6)

dd_best_3 = MuMIn::get.models(dd_3, subset = 1)[[1]]

### simple
dd_4 = dd %>%
  subset((has('vol_convex', 'lag:vol_convex'))) %>%
  subset(df == 7)

dd_best_4 = MuMIn::get.models(dd_4, subset = 1)[[1]]

dd_5 = dd %>%
  subset((has('vol_convex', 'lag:vol_convex'))) %>%
  subset(df == 8)

dd_best_5 = MuMIn::get.models(dd_5, subset = 1)[[1]]


dd_6 = dd %>%
  subset((has('vol_convex', 'lag:vol_convex'))) %>%
  subset(df == 9)

dd_best_6 = MuMIn::get.models(dd_6, subset = 1)[[1]]


dd_8 = dd %>%
  subset((has('vol_convex', 'lag:vol_convex'))) %>%
  subset(df == 11)

dd_best_8 = MuMIn::get.models(dd_8, subset = 1)[[1]]


dd_9 = dd %>%
  subset((has('vol_convex', 'lag:vol_convex'))) %>%
  subset(df == 12)

dd_best_9 = MuMIn::get.models(dd_9, subset = 1)[[1]]

# dd_10 = dd %>%
#   subset((has('vol_convex', 'lag:vol_convex'))) %>%
#   subset(df == 13)
#
# dd_best_10 = MuMIn::get.models(dd_10, subset = 1)[[1]]


# model for height wih height lag
mod_z = lm(blups_mod,
           formula = HT16 ~
             Zq999 +
             lag:Zq999 + 1,
           na.action = 'na.fail')


candidate_list = c(#dd_best_complex$call$formula,
                   dd_best_simple$call$formula,
                   dd_best_vol_zq999$call$formula,
                   dd_best_struct$call$formula,
                   dd_best_spec$call$formula,
                   dd_best_vol_ndre1$call$formula,
                   dd_best_vol_ndvi$call$formula,
                   dd_best_vol_cci$call$formula,
                   dd_best_vol_nirvndre$call$formula,
                   dd_best_vol_sipi$call$formula,
                   dd_best_3$call$formula,
                   dd_best_4$call$formula,
                   dd_best_5$call$formula,
                   dd_best_6$call$formula,
                   dd_best_8$call$formula,
                   dd_best_9$call$formula,
                   #dd_best_10$call$formula,
                   mod_z$call$formula) %>%
  lapply(function(x) substr(unlist(x), 1, nchar(x) - 4)[3]) %>%
  unlist()

################################################################################
# table of model coefficients for manuscript

cand_mod_list = c(#dd_complex$rowname[1],
                  dd_simple$rowname[1],
                  dd_vol_zq999$rowname[1],
                  dd_struct$rowname[1],
                  dd_spec$rowname[1],
                  dd_vol_ndre1$rowname[1],
                  dd_vol_ndvi$rowname[1],
                  dd_vol_cci$rowname[1],
                  dd_vol_nirvndre$rowname[1],
                  dd_vol_sipi$rowname[1],
                  dd_3$rowname[1],
                  dd_4$rowname[1],
                  dd_5$rowname[1],
                  dd_6$rowname[1],
                  dd_8$rowname[1],
                  dd_9$rowname[1])

lag_coeff_list = data.frame(lag = c(#tail(dd_best_complex$coefficients, n = 1),
                                    tail(dd_best_simple$coefficients, n = 1),
                                    tail(dd_best_vol_zq999$coefficients, n = 1),
                                    tail(dd_best_struct$coefficients, n = 1),
                                    tail(dd_best_spec$coefficients, n = 1),
                                    tail(dd_best_vol_ndre1$coefficients, n = 1),
                                    tail(dd_best_vol_ndvi$coefficients, n = 1),
                                    tail(dd_best_vol_cci$coefficients, n = 1),
                                    tail(dd_best_vol_nirvndre$coefficients, n = 1),
                                    tail(dd_best_vol_sipi$coefficients, n = 1),
                                    tail(dd_best_3$coefficients, n = 1),
                                    tail(dd_best_4$coefficients, n = 1),
                                    tail(dd_best_5$coefficients, n = 1),
                                    tail(dd_best_6$coefficients, n = 1),
                                    tail(dd_best_8$coefficients, n = 1),
                                    tail(dd_best_9$coefficients, n = 1)),
                                    #tail(dd_best_10$coefficients, n = 1)),
                            rowname = cand_mod_list)


(dd_table = dd %>%
  subset(rowname %in% cand_mod_list) %>%
  data.frame() %>%
  dplyr::select(rowname:Zq999, logLik, AICc) %>%
  left_join(lag_coeff_list, by = "rowname") %>%
  dplyr::select(vol_convex, lag,
                any_of(keeps_list),
                logLik, AICc) %>%
  mutate(across(vol_convex:AICc, round, 2)) %>%
  add_row(.before = 1) %>%
  flextable::flextable())

# save_as_docx(
#   'Table of model parameters' = dd_table,
#   path = 'D:\\Sync\\Figures\\_TABLES_paper3\\paper3_model_params.docx')

################################################################################
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

# cross validation by SITE
val_mods_site = blups_mod %>%
  crossv_strat(var = 'Site') %>%
  # fit each model with the training data
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
         cross_v = 'Site')


val_mods = val_mods_prov %>%
  bind_rows(val_mods_site) %>%
  distinct()

val_mods_sum = val_mods %>%
  group_by(model_name, Site, cross_v) %>%
  mutate(R2 = if_else(model_name == 'mod_z',
                         yardstick::rsq_trad_vec(.fitted, HT16),
                         yardstick::rsq_trad_vec(.fitted, DBH16)),
            rmse = if_else(model_name == 'mod_z',
                           caret::RMSE(.fitted, HT16),
                           caret::RMSE(.fitted, DBH16))) %>%
  dplyr::select(model_name, Site, cross_v, R2, rmse) %>%
  distinct()




################################################################################
# HALF NORMAL functions

by_pop_site_mean = by_pop %>%
  ungroup() %>%
  filter(variable %in% c('DBH16', 'HT16')) %>%
  dplyr::select(Site, site_mean, variable) %>%
  mutate(phen = if_else(variable == 'DBH16', 'DBH', 'Height')) %>%
  distinct()

val_mods_blup = val_mods %>%
  # add rows with the field measurements
  bind_rows(blups_mod) %>%
  replace_na(list(model_name = "Field height",
                  phen = "Height")) %>%
  mutate(.fitted = if_else(is.na(.fitted), HT16, .fitted)) %>%
  bind_rows(blups_mod) %>%
  replace_na(list(model_name = "Field DBH",
                  phen = "DBH",
                  cross_v = "Population")) %>%
  mutate(.fitted = if_else(is.na(.fitted), DBH16, .fitted)) %>%
  bind_rows(blups_mod) %>%
  replace_na(list(model_name = "Field DBH",
                  phen = "DBH",
                  cross_v = "Site")) %>%
  mutate(.fitted = if_else(is.na(.fitted), DBH16, .fitted)) %>%
  left_join(by_pop_site_mean) %>%
  distinct() %>%
  mutate(TRAIT = .fitted + site_mean)

halfnorm = val_mods_blup %>%
  group_by(Site, Class, model_name, cross_v) %>%
  mutate(max_y = max(TRAIT),
         min_x = max(Euc)) %>%
  nest() %>%
  mutate(model = purrr::map(data,
                            ~nls(TRAIT ~ a * exp(-0.5 * Euc^2/b^2),
                                 start = list(a = first(.$max_y), b = first(.$min_x)),
                                 control = nls.control(maxiter = 5000),
                                 data = ., trace = TRUE))) %>%
  mutate(preds = map2(data, model, add_predictions)) %>%
  mutate(preds = map(preds, dplyr::select, pred)) %>%
  mutate(resids = map2(data, model, add_residuals)) %>%
  mutate(resids = map(resids, dplyr::select, resid)) %>%
  mutate(glance = map(model, broom::glance))  %>%
  # get the a and sigma parameters of the function
  mutate(par_a = map(model, function(x){x$m$getPars()[[1]]}),
         par_b = map(model, function(x){x$m$getPars()[[2]]})) %>%
  #mutate(rmse_halfnorm = map(model, qpcR::RMSE)) %>%
  unnest(c(data, preds, resids, glance, par_a, par_b)) %>%
  dplyr::select(model_name, Site, Class, Euc, TRAIT, model, pred, resid, sigma, phen, logLik:deviance, par_a, par_b) %>%
  distinct() %>%
  mutate(legend_id = case_when(
    Class == "Orchard" & model_name %in% c("Field DBH", "Field height") ~ "orch_meas",
    Class == "Orchard" & !(model_name %in% c("Field DBH", "Field height")) ~ "orch_pred",
    Class == "Wildstand" & model_name %in% c("Field DBH", "Field height") ~ "wild_meas",
    Class == "Wildstand" & !(model_name %in% c("Field DBH", "Field height")) ~ "wild_pred",
  ))


################################################################################
# Model selection


dbh_table_rank = halfnorm %>%
  dplyr::select(model_name, Site, cross_v, TRAIT, AIC, BIC, phen, legend_id) %>%
  left_join(val_mods_sum,
            by = c("model_name", "Site", "cross_v")) %>%
  distinct() %>%
  #select(order(colnames(.))) %>%
  group_by(Site, Class, cross_v, phen) %>%
  mutate(rmse_rank = percent_rank(rmse),
         rmse_mean = mean(rmse, na.rm = TRUE),
         rmse_min = min(rmsena.rm = TRUE),
         BIC_rank = percent_rank(BIC),
         BIC_mean = mean(BIC),
         BIC_min = min(BIC),
         AIC_rank = percent_rank(AIC),
         AIC_mean = mean(AIC),
         AIC_min = min(AIC),
         R2_rank = percent_rank(desc(R2))) %>%
  group_by(model_name) %>%
  mutate(rmse_rel = rmse - rmse_min,
         BIC_rel = BIC - BIC_min,
         AIC_rel = AIC - AIC_min,
         mean_AIC = mean(AIC),
         mean_rank_AIC = mean(AIC_rank),
         mean_rank_rmse = mean(rmse_rank),
         mean_rank_combined = (mean_rank_AIC + mean_rank_rmse) / 2,
         mean_rmse_rel = mean(rmse_rel),
         mean_AIC_mean = mean(AIC_mean)) %>%
  # what are you selecting by?
  ungroup() %>%
  dplyr::select(-TRAIT) %>%
  distinct()


# compare models' fit to choose the best one


mod_labs = c(
    model_6 = expression(~V['convex']~'+'#~italic('lag')~V['convex']~'+'
                       ~Z['Q999']~'+'
                       ~CCI~'+'~NDRE['705']~'+'~NDVI~'+'~NIR['V'^'NDRE']),
    model_4 = expression(~V['convex']~'+'#~italic('lag')~V['convex']~'+'
                         ~Z['Q999']~'+'
                         ~NDRE['705']~'+'~NIR['V'^'NDRE']),
    model_8 = expression(~V['convex']~'+'#~italic('lag')~V['convex']~'+'
                       ~Z['Q999']~'+'~CRR~'+'
                       ~CCI~'+'~NDRE['705']~'+'~NDVI~'+'~NIR['V'^'NDRE']~'+'
                       ~SIPI),
    model_5 = expression(~V['convex']~'+'#~italic('lag')~V['convex']~'+'
                         ~Z['Q999']~'+'
                         ~NDRE['705']~'+'~NDVI~'+'~NIR['V'^'NDRE']),
    model_spec = expression(~V['convex']~'+'#~italic('lag')~V['convex']~'+'
                            ~CCI~'+'~NDRE['705']~'+'~NDVI~'+'~NIR['V'^'NDRE']~'+'
                            ~SIPI),
    model_9 = expression(~V['convex']~'+'#~italic('lag')~V['convex']~'+'
                       ~Z['Q999']~'+'~CRR~'+'~theta['apex']~'+'
                       ~CCI~'+'~NDRE['705']~'+'~NDVI~'+'~NIR['V'^'NDRE']~'+'
                       ~SIPI),
  model_3 = expression(~V['convex']~'+'#~italic('lag')~V['convex']~'+'
                       ~Z['Q999']~'+'
                       ~NIR['V'^'NDRE']),
  model_nirvndre = expression(~V['convex']~'+'#~italic('lag')~V['convex']~'+'
                              ~NIR['V'^'NDRE']),
  model_ndre1 = expression(~V['convex']~'+'#~italic('lag')~V['convex']~'+'
                              ~NDRE['705']),
  model_cci = expression(~V['convex']~'+'#~italic('lag')~V['convex']~'+'
                           ~CCI),
  model_sipi = expression(~V['convex']~'+'#~italic('lag')~V['convex']~'+'
                              ~SIPI),
  model_vol_zq999 = expression(~V['convex']~'+'#~italic('lag')~V['convex']~'+'
                              ~Z['Q999']),
  model_ndvi = expression(~V['convex']~'+'#~italic('lag')~V['convex']~'+'
                          ~NDVI),
  model_struct = expression(~V['convex']~'+'#~italic('lag')~V['convex']~'+'
                            ~CRR~'+'~theta['apex']),
  model_simple = expression(~V['convex']#~italic('lag')~V['convex']
                            )
)

(s4 = dbh_table_rank %>%
   filter(model_name != 'mod_z' &
            model_name != 'Field height' &
            model_name != 'Field DBH' &
            Class != "Orchard") %>%
    pivot_longer(cols = c(rmse_rel, AIC_rel), values_to = 'value', names_to = 'metric') %>%
    mutate(metric = recode(metric, AIC_rel = "ΔAIC of half-normal\ntransfer function",
                                     rmse_rel = "ΔRMSE of cross-validation\nregression (mm)"),
           cross_v = recode(cross_v, Population = "Population-wise cross validation",
                            Site = "Site-wise cross validation")) %>%
    group_by(metric, cross_v, model_name) %>%
    mutate(mean_val = mean(value)) %>%
    ggplot(aes(x = reorder(model_name, desc(mean_rank_combined)), y = value, color = Site, group = Site)) +
    coord_flip() +
   geom_point(size = 4, alpha = .6) +
   geom_line(size = 1.7, alpha = .5) +
    #geom_point(aes(y = mean_val), color = "grey20", size = 3) +
    geom_line(aes(y = mean_val), color = "grey20", size = 1.1) +
   theme_bw(base_size = 18) +
    scale_x_discrete(labels = mod_labs) +
    theme(axis.title = element_blank(),
          axis.text.y = element_text(hjust = 0, family = "Cambria", size = 20),
          axis.text.x = element_text(family = "Cambria", size = 16),
          legend.position = c(.87, .935),
          legend.margin = margin(0, 3, 3, 3),
          legend.box.background = element_rect(size = .5),
          legend.text = element_text(size = 15, family = "Cambria"),
          strip.background = element_rect(fill = "white", color = "white"),
          strip.text = element_text(family = "Cambria", size = 23),
          legend.title = element_blank(),
          panel.border = element_rect(color = "black", size = 1)) +
    scale_color_brewer(palette = "Dark2") +
   facet_grid(cross_v ~ metric,
              scales = "free"))


ggsave(plot = s4,
       filename = "D:\\Sync\\Figures\\_FIGURES_GCB\\figure_S4.tiff",
       device = tiff,
       width = 19,
       height = 15,
       units = "in",
       dpi = 600)

ggsave(plot = s4,
       filename = "D:\\Sync\\Figures\\_FIGURES_GCB\\figure_S4.jpeg",
       device = jpeg,
       width = 19,
       height = 15,
       units = "in",
       dpi = 600)



  (table_pub = val_mods_sum %>%
      ungroup() %>%
     filter(model_name %in% c('mod_z', 'model_6')) %>%
     mutate(across(R2:rmse, round, 3)) %>%
     dplyr::select(model_name:rmse) %>%
     distinct() %>%
     pivot_longer(cols = c('R2', 'rmse'), names_to = 'fit_metric', values_to = 'fit_value') %>%
     pivot_wider(names_from = c('cross_v', 'fit_metric'), values_from = 'fit_value', names_glue = '{cross_v}_{fit_metric}') %>%
     add_row(.before = 1) %>%
     mutate(model_name = replace_na(model_name, 'Height')) %>%
     add_row(.before = 8) %>%
     mutate(model_name = replace_na(model_name, 'DBH')) %>%
     # mutate(across(2:13, .fns = ~ ifelse(is.na(.x) == TRUE, cur_column(), .x))) %>%
     # separate_rows(2:13, sep = '_') %>%
     flextable::flextable())

save_as_docx(
  'Table of cross validation model fits by site' = table_pub,
  path = 'D:\\Sync\\Figures\\_TABLES_paper3\\paper3_table_pub_GCB.docx')
################################################################################

val_mods_plot = filter(val_mods,
                       model_name == 'model_6' &
                         cross_v == "Population") %>%
  distinct()

val_mods_sum_plot = filter(val_mods_sum, model_name == 'model_6' &
                             cross_v == "Population")

# predictions vs DBH
# predictions vs DBH
# predictions vs DBH
(g7 = ggplot(val_mods_sum_plot) +
  geom_abline(slope = 1, intercept = 0, color = 'black') +
    geom_point(data = arrange(val_mods_plot, desc(Prov)), aes(y = .fitted, x = DBH16, color = Class), size = 3.4, alpha = .5
               #color = '#607080'
    ) +
    geom_smooth(data = val_mods_plot, aes(y = .fitted, x = DBH16),
                method = 'lm',
                color = '#B50010',
                se = FALSE,
                size = 1.2) +
  theme_bw(base_size = 18) +
    scale_color_manual(values = class_cols) +
  labs(x = 'Field-assessed DBH BLUPs (mm)',
       y = 'Model-predicted DBH values (mm)') +
  # expression(bold(Slope~of~the~upper~red~edge~(RE[UPPER])))
    geom_text(x = min(val_mods_plot$.fitted) +
                ((max(val_mods_plot$.fitted) - min(val_mods_plot$.fitted)) * .08),
              y = max(val_mods_plot$DBH16) * .8,
              label = paste0('R^2'),
              size = 6,
              parse = TRUE,
              vjust = .5,
              hjust = 0,
              family = "Cambria") +
    geom_text(x = min(val_mods_plot$.fitted) +
                ((max(val_mods_plot$.fitted) - min(val_mods_plot$.fitted)) * .145),
              y = max(val_mods_plot$DBH16) * .785,
              label = paste0("= ", numform::f_num(val_mods_sum_plot$R2, digits = 2, retain.leading.zero = FALSE)),
              size = 6,
              parse = FALSE,
              vjust = .5,
              hjust = 0,
              family = "Cambria") +
  geom_text(x = min(val_mods_plot$.fitted) +
              ((max(val_mods_plot$.fitted) - min(val_mods_plot$.fitted)) * .08),
            y = max(val_mods_plot$DBH16) * .62,
            label = paste0('RMSE = ', numform::f_num(val_mods_sum_plot$rmse, digits = 2, retain.leading.zero = TRUE)),
            size = 6,
            parse = FALSE,
            vjust = .5,
            hjust = 0,
            family = "Cambria") +
  # geom_text(x = min(val_mods_plot$.fitted) +
  #             ((max(val_mods_plot$.fitted) - min(val_mods_plot$.fitted)) * .08),
  #           y = max(val_mods_plot$DBH16) * .6,
  #           label = paste0('MAE == ', round(val_mods_sum_plot$mae, digits = 2)),
  #           size = 5.5,
  #           parse = TRUE,
  #           vjust = .5,
  #           hjust = 0) +
  theme(legend.position = 'none',
        strip.background = element_rect(fill = 'white'),
        axis.text = element_text(family='Cambria', size=16, color = 'black'),
        plot.background = element_rect(color = 'white', fill = 'white', size = 1),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(margin = margin(r = 10, l = -10), vjust = 1),
        plot.margin = margin(.4, .4, .4, .4, 'cm'),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        axis.title.y = element_text(size = 22, margin = margin(t = 0, r = 30, b = 0, l = 0), family = 'Cambria'),
        axis.title.x = element_text(size = 22, margin = margin(t = 12, r = 0, b = 0, l = 0), family = 'Cambria'),
        strip.text = element_text(size = 24, face = "bold", family = "Cambria")) +
  facet_wrap(. ~ Site,
             nrow = 2))

ggsave(plot = g7,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\figure_7.tiff',
       device = tiff,
       width = 15,
       height = 10,
       units = 'in',
       dpi = 600)

ggsave(plot = g7,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\figure_7.jpeg',
       device = jpeg,
       width = 15,
       height = 10,
       units = 'in',
       dpi = 300)

################################################################################
# same thing but for height, SI

val_mods_plot = filter(val_mods,
                       model_name == 'mod_z' &
                         cross_v == "Population") %>%
  distinct()

val_mods_sum_plot = filter(val_mods_sum, model_name == 'mod_z' &
                             cross_v == "Population")
# predictions vs DBH
# predictions vs DBH
# predictions vs DBH
(s1 = ggplot(val_mods_sum_plot) +
    geom_abline(slope = 1, intercept = 0, color = 'black') +
    geom_point(data = arrange(val_mods_plot, desc(Prov)), aes(y = .fitted, x = HT16, color = Class), size = 3.4, alpha = .5
               #color = '#607080'
                 ) +
    geom_smooth(data = val_mods_plot, aes(y = .fitted, x = HT16),
                method = 'lm',
                color = '#B50010',
                se = FALSE,
                size = 1.2) +
    theme_bw(base_size = 18) +
    scale_color_manual(values = class_cols) +
    labs(x = 'Field-assessed height BLUPs (m)',
         y = 'Model-predicted height values (m)') +
    # expression(bold(Slope~of~the~upper~red~edge~(RE[UPPER])))
    geom_text(x = min(val_mods_plot$.fitted) +
                ((max(val_mods_plot$.fitted) - min(val_mods_plot$.fitted)) * .08),
              y = max(val_mods_plot$HT16) * .8,
              label = paste0('R^2'),
              size = 6,
              parse = TRUE,
              vjust = .5,
              hjust = 0,
              family = "Cambria") +
    geom_text(x = min(val_mods_plot$.fitted) +
                ((max(val_mods_plot$.fitted) - min(val_mods_plot$.fitted)) * .145),
              y = max(val_mods_plot$HT16) * .785,
              label = paste0("= ", numform::f_num(val_mods_sum_plot$R2, digits = 3, retain.leading.zero = FALSE)),
              size = 6,
              parse = FALSE,
              vjust = .5,
              hjust = 0,
              family = "Cambria") +
    geom_text(x = min(val_mods_plot$.fitted) +
                ((max(val_mods_plot$.fitted) - min(val_mods_plot$.fitted)) * .08),
              y = max(val_mods_plot$HT16) * .62,
              label = paste0('RMSE = ', f_num(val_mods_sum_plot$rmse, digits = 2, retain.leading.zero = TRUE)),
              size = 6,
              parse = FALSE,
              vjust = .5,
              hjust = 0,
              family = "Cambria") +
    # geom_text(x = min(val_mods_plot$.fitted) +
    #             ((max(val_mods_plot$.fitted) - min(val_mods_plot$.fitted)) * .08),
    #           y = max(val_mods_plot$DBH16) * .6,
    #           label = paste0('MAE == ', round(val_mods_sum_plot$mae, digits = 2)),
    #           size = 5.5,
    #           parse = TRUE,
    #           vjust = .5,
    #           hjust = 0) +
    theme(legend.position = 'none',
          strip.background = element_rect(fill = 'white'),
          axis.text = element_text(family='Cambria', size=16, color = 'black'),
          plot.background = element_rect(color = 'white', fill = "white", size = 1),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(margin = margin(r = 10, l = -10), vjust = 1),
          plot.margin = margin(.4, .4, .4, .4, 'cm'),
          panel.border = element_rect(color = 'black', fill = NA, size = 1),
          axis.title.y = element_text(size = 22, margin = margin(t = 0, r = 30, b = 0, l = 0), family = 'Cambria'),
          axis.title.x = element_text(size = 22, margin = margin(t = 12, r = 0, b = 0, l = 0), family = 'Cambria'),
          strip.text = element_text(size = 24, face = "bold", family = "Cambria")) +
    facet_wrap(. ~ Site,
               nrow = 2))

ggsave(plot = s1,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\figure_S1.tiff',
       device = tiff,
       width = 15,
       height = 10,
       units = 'in',
       dpi = 600)

ggsave(plot = s1,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\figure_S1.jpeg',
       device = jpeg,
       width = 15,
       height = 10,
       units = 'in',
       dpi = 300)

################################################################################

#class_cols = c("Wildstand" = "#366B9A", "Orchard" = "#30A275")

halfnorm_plot = halfnorm %>%
  filter(model_name %in% c("model_6", "mod_z", "Field DBH", "Field height") &
           cross_v != "Site")

halfnorm_table = halfnorm_plot %>%
  # mutate(x = max(Euc),
  #        y = max(BLUP)) %>%
  dplyr::select(Site, Class, model_name, AIC, par_a, par_b) %>%
  distinct() %>%
  mutate(across(AIC:par_b, round, 2)) %>%
  arrange(Class) %>%
  group_by(Site) %>%
  nest()

data.tb = tibble(tb = list(halfnorm_table))

plot_labs = c("orch_meas" = "",
              "orch_pred" = "",
              "wild_meas" = "Field assessment",
              "wild_pred" = "Model prediction")

legend_title = "  Orchard     Wildstand"


(g8 = halfnorm_plot %>% filter(Site == 'Tête Jaune Cache'
                                 ) %>%
    mutate(phen = factor(phen, levels = c('Height', 'DBH'))) %>%
    mutate(phen = recode_factor(phen, Height = "Height (m)", DBH = "DBH (mm)")) %>%
  ggplot(aes(x = Euc, y = TRAIT, color = legend_id, alpha = legend_id, shape = legend_id, linetype = legend_id)) +
  geom_point(aes(size = legend_id)) +
  labs(x = 'Climate transfer distance (Euclidean units)',
       y = 'DBH (mm)') +
  #geom_smooth(size = 1.3, alpha = 1, se = FALSE, method = 'lm') +
  geom_line(aes(y = pred), size = 1.1) +
  #geom_hline(aes(yintercept = par_a * .9)) +
  theme_bw(base_size = 20) +
  scale_color_manual(name = legend_title,
                     values = c('orch_meas' = class_cols[[2]], #'#D67101',
                                'orch_pred' = class_cols[[2]], #'#D67101',
                                'wild_meas' = class_cols[[1]], #'#1B7588',
                                'wild_pred' = class_cols[[1]]), #'#1B7588'),
                     labels = plot_labs) +
  scale_shape_manual(name = legend_title,
                     values = c('orch_meas' = 16,
                                'orch_pred' = 16,
                                'wild_meas' = 16,
                                'wild_pred' = 16),
                     labels = plot_labs) +
  scale_alpha_manual(name = legend_title,
                     values = c('orch_meas' = 1,
                                'orch_pred' = .45,
                                'wild_meas' = 1,
                                'wild_pred' = .45),
                     labels = plot_labs) +
  scale_linetype_manual(name = legend_title,
                        values = c('orch_meas' = 1,
                                   'orch_pred' = 1,
                                   'wild_meas' = 1,
                                   'wild_pred' = 1),
                        labels = plot_labs) +
    scale_size_manual(name = legend_title,
                          values = c('orch_meas' = 2.5,
                                     'orch_pred' = 4.4,
                                     'wild_meas' = 2.5,
                                     'wild_pred' = 4.4),
                          labels = plot_labs) +
  theme(legend.position = c(0.84, 0.925),
        legend.background = element_rect(fill = "white", color = "black", size = .5),
        panel.border = element_rect(fill = NA, color = "black", size = 1),
        legend.key.width = unit(2.75,"cm"),legend.key.height = unit(.7,"cm"),
        legend.margin = margin(4, 4, 4, 4),
        legend.spacing.x = unit(.1, 'cm'),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 24, family = "Cambria"),
        strip.placement = "left",
        strip.background.y = element_blank(),
        strip.switch.pad.wrap = unit(.3, 'cm'),
        legend.title = element_text(size = 17, family = "Cambria"),
        legend.text = element_text(size = 17, family = "Cambria"),
        axis.text = element_text(family='Cambria', size = 18, color = 'black'),
        axis.title.y = element_blank(),
        axis.title.x = element_text(vjust = -.5, size = 24, family = "Cambria")) +
  guides(color = guide_legend(ncol = 2)) +
    #scale_y_continuous(expand = expansion(mult = c(.06, .16))) + #Kalamalka
    #scale_y_continuous(expand = expansion(mult = c(.06, .13))) + #TJC
    #scale_y_continuous(expand = expansion(mult = c(.06, .08))) + #Skim, Whit
  facet_wrap(. ~ phen, scales = 'free', strip.position = "left"))

ggsave(plot = g8,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\figure_8_TJC.tiff',
       device = tiff,
       width = 15,
       height = 10,
       units = 'in',
       dpi = 600)

ggsave(plot = g8,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\figure_8_TJC.jpeg',
       device = jpeg,
       width = 15,
       height = 10,
       units = 'in',
       dpi = 300)


# COLUMN - improves transfer fit yes/no


(table_halfnorm_pub = halfnorm_plot %>%
    ungroup() %>%
    dplyr::select(Site:Class, phen, legend_id, AIC, par_a, par_b) %>%
    distinct() %>%
  arrange(desc(phen)) %>%
  mutate(across(AIC:par_b, round, 2)) %>%
    distinct() %>%
    pivot_longer(cols = c('AIC', 'par_a', 'par_b'), names_to = 'fit_metric', values_to = 'fit_value') %>%
    pivot_wider(names_from = c('Class', 'legend_id', 'fit_metric'), values_from = 'fit_value', names_glue = '{Class}_{legend_id}_{fit_metric}') %>%
  dplyr::select(order(colnames(.))) %>%
  dplyr::select(phen:Site, everything()) %>%
    add_row(.before = 1) %>%
    #mutate(model_name = replace_na(model_name, 'Height')) %>%
    add_row(.before = 8) %>%
    #mutate(model_name = replace_na(model_name, 'DBH')) %>%
    # mutate(across(2:13, .fns = ~ ifelse(is.na(.x) == TRUE, cur_column(), .x))) %>%
    # separate_rows(2:13, sep = '_') %>%
    flextable::flextable())

save_as_docx(
  'Table of half-normal fits by site' = table_halfnorm_pub,
  path = 'D:\\Sync\\Figures\\_TABLES_paper3\\paper3_table_halfnorm_pub_GCB.docx')


library(ggh4x)

(s9 = halfnorm_plot %>%
  pivot_longer(cols = c(par_a, par_b), values_to = 'value', names_to = 'par') %>%
  mutate(par = factor(par, labels = c("italic(a)", "italic(sigma)^2"))) %>%
  #filter(par == "par_a") %>%
  # mutate(metric = recode(metric, AIC_rel = "ΔAIC of half-normal\ntransfer function",
  #                        rmse_rel = "ΔRMSE of cross-validation\nregression (mm)"),
  #        cross_v = recode(cross_v, Population = "Population-wise cross validation",
  #                         Site = "Site-wise cross validation")) %>%
  # group_by(metric, cross_v, model_name) %>%
  # mutate(mean_val = mean(value)) %>%
  ggplot(aes(x = Site, y = value, color = legend_id, alpha = legend_id,
             shape = legend_id, linetype = legend_id, group = legend_id)) +
 # geom_point(aes(size = legend_id)) +
  coord_flip() +
  geom_line(size = 1) +
  #geom_point(aes(y = mean_val), color = "grey20", size = 3) +
  #geom_line(aes(y = mean_val), color = "grey20", size = 1.1) +
  theme_bw(base_size = 18) +
  #scale_x_discrete(labels = mod_labs) +
  scale_color_manual(name = legend_title,
                     values = c('orch_meas' = class_cols[[2]], #'#D67101',
                                'orch_pred' = class_cols[[2]], #'#D67101',
                                'wild_meas' = class_cols[[1]], #'#1B7588',
                                'wild_pred' = class_cols[[1]]), #'#1B7588'),
                     labels = plot_labs) +
  scale_shape_manual(name = legend_title,
                     values = c('orch_meas' = 16,
                                'orch_pred' = 16,
                                'wild_meas' = 16,
                                'wild_pred' = 16),
                     labels = plot_labs) +
  scale_alpha_manual(name = legend_title,
                     values = c('orch_meas' = 1,
                                'orch_pred' = .45,
                                'wild_meas' = 1,
                                'wild_pred' = .45),
                     labels = plot_labs) +
  scale_linetype_manual(name = legend_title,
                        values = c('orch_meas' = 1,
                                   'orch_pred' = 1,
                                   'wild_meas' = 1,
                                   'wild_pred' = 1),
                        labels = plot_labs) +
  scale_size_manual(name = legend_title,
                    values = c('orch_meas' = 2.5,
                               'orch_pred' = 4.4,
                               'wild_meas' = 2.5,
                               'wild_pred' = 4.4),
                    labels = plot_labs) +
  theme(legend.position = c(0.83, 0.434),
        legend.background = element_rect(fill = "white", color = "black", size = .5),
        panel.border = element_rect(fill = NA, color = "black", size = 1),
        legend.key.width = unit(2.55,"cm"),legend.key.height = unit(.7,"cm"),
        legend.margin = margin(4, 4, 4, 4),
        legend.spacing.x = unit(.1, 'cm'),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 24, family = "Cambria"),
        strip.placement = "right",
        strip.background.y = element_blank(),
        strip.switch.pad.wrap = unit(.3, 'cm'),
        legend.title = element_text(size = 16, family = "Cambria"),
        legend.text = element_text(size = 16, family = "Cambria"),
        axis.text = element_text(family='Cambria', size = 18, color = 'black'),
        axis.title = element_blank()) +
  guides(color = guide_legend(ncol = 2)) +
  facet_grid2(par ~ phen,
             scales = "free",
             independent = "x",
             labeller = label_parsed))



ggsave(plot = s9,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\figure_s9.tiff',
       device = tiff,
       width = 15,
       height = 14,
       units = 'in',
       dpi = 600)

ggsave(plot = s9,
       filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\figure_s9.jpeg',
       device = jpeg,
       width = 15,
       height = 14,
       units = 'in',
       dpi = 300)



























#
#
# ################################################################################
#
# # predictions vs DBH
# g6 = ggplot(val_mods_sum_plot) +
#   geom_point(data = val_mods_plot, aes(y = .fitted, x = DBH16), color = "steelblue4", size = 2, alpha = .3) +
#   geom_smooth(data = val_mods_plot, aes(y = .fitted, x = DBH16),
#               method = "lm", size = 1.2, color = "red3", se = FALSE) +
#   geom_text(x = min(val_mods_plot$.fitted) +
#                    ((max(val_mods_plot$.fitted) - min(val_mods_plot$.fitted)) * .08),
#             y = max(val_mods_plot$DBH16) * .9,
#             label = paste0("R^2 == ", round(val_mods_sum_plot$R2, digits = 2)),
#             size = 5,
#             parse = TRUE,
#             vjust = 1,
#             hjust = 0) +
#   geom_text(x = min(val_mods_plot$.fitted) +
#               ((max(val_mods_plot$.fitted) - min(val_mods_plot$.fitted)) * .08),
#             y = max(val_mods_plot$DBH16) * .82,
#             label = paste0("RMSE == ", round(val_mods_sum_plot$rmse, digits = 2)),
#             size = 5,
#             parse = TRUE,
#             vjust = 1,
#             hjust = 0) +
#   # geom_text(x = min(val_mods_plot$.fitted) +
#   #             ((max(val_mods_plot$.fitted) - min(val_mods_plot$.fitted)) * .08),
#   #           y = max(val_mods_plot$DBH16) * .76,
#   #           label = paste0("MAE == ", round(val_mods_sum_plot$mae, digits = 2)),
#   #           size = 5,
#   #           parse = TRUE,
#   #           vjust = 1,
#   #           hjust = 0) +
#   # labs(x = "Euclidean climatic transfer distance",
#   #      y = "BLUP") +
#   theme_bw(base_size = 18) +
#   theme(strip.background = element_rect(fill = "white"),
#         strip.text = element_text(size = 22),
#         legend.title = element_text(size = 12),
#         legend.text = element_text(size = 12),
#         axis.title.y = element_text(vjust = 2, size = 20),
#         axis.title.x = element_text(vjust = -.5, size = 20)) +
#   scale_x_continuous(breaks = seq(0, 150, 50)) +
#   facet_wrap(. ~ Site,
#              nrow = 2)
#
#
#
# ggsave(plot = g6,
#        filename = "D:\\Sync\\Figures\\_FIGURES_paper3\\figure_6.tiff",
#        device = "tiff",
#        width = 15,
#        height = 10,
#        units = "in",
#        dpi = 600)
#
# ggsave(plot = g6,
#        filename = "D:\\Sync\\Figures\\_FIGURES_paper3\\figure_6.jpeg",
#        device = "jpeg",
#        width = 15,
#        height = 10,
#        units = "in",
#        dpi = 300)
#
#
#
# # predictions vs DBH
# ggplot(val_mods_sum_plot) +
#   geom_point(data = val_mods_plot, aes(y = .resid, x = DBH16), color = "steelblue4", size = 3, alpha = .3) +
#   geom_smooth(data = val_mods_plot, aes(y = .resid, x = DBH16),
#               method = "loess", size = 1.5, color = "indianred", se = FALSE) +
# theme_bw(base_size = 20) +
#   facet_wrap(Site ~ ., nrow = 2)
#
#
# ###############################################################################
#
#
# val_mods_blup = val_mods %>%
#   group_by(model_name, Site) %>%
#   rename(var_value = .fitted) %>%
#   nest() %>%
#   mutate(mm = map(data, blup_mod)) %>%
#   mutate(mm_weev = map(data, blup_mod_weev)) %>%
#   mutate(mm_weev = map(mm_weev, dplyr::select, BLUP_weev)) %>%
#   unnest(c(mm, mm_weev)) %>%
#   dplyr::select(-data) %>%
#   left_join(
#     # get the Euc column for each pop at each site
#     dplyr::select(by_pop, Prov, Site, Euc, Class), by = c("Prov", "Site")) %>%
#   distinct() %>%
#     # then append the BLUPS for the census data
#   bind_rows(filter(by_pop, variable == "DBH16", mean_type != "var_mean") %>%
#       pivot_wider(names_from = mean_type, values_from = var_value) %>%
#       dplyr::select(-variable, -Rep, -Blk, -Weev16) %>%
#       rename("model_name" = Group))
#
#
#
# # plot the BLUPS for the model predictions
#
# val_mods_blup_plot = val_mods_blup %>%
#   filter(model_name %in% c(#"model_ndre1",
#                            "model_spec_site",
#                            "Census Traits"))
#
# # plot the model at each site
# ggplot(val_mods_blup_plot,
#        aes(x = Euc, y = BLUP_weev, color = model_name)) +
#   geom_point(size = 4, alpha = .7) +
#   geom_smooth(span = 2, size = 1) +
#   theme_bw(base_size = 20) +
#   # geom_text(x = max(by_site$Euc) * .7,
#   #           y = Inf,
#   #           label = paste0("R^2 == ", round(by_site_best_plot$adj.r.squared, digits = 2)),
#   #           size = 5,
#   #           parse = TRUE,
#   #           vjust = 1.5,
#   #           hjust = 0) +
#   labs(x = "Euclidean climatic transfer distance",
#        y = "BLUP") +
#   # geom_text(x = max(by_site$Euc) * .88,
#   #           y = Inf,
#   #           label = paste0("n == ", by_site_best$nobs),
#   #           size = 5,
#   #           parse = TRUE,
#   #           vjust = 3.6,
#   #           hjust = 0) +
#   # geom_text(x = max(by_site$Euc) * .88,
#   #           y = Inf,
#   #           label = paste0("AIC == ", round(by_site_best$AIC, digits = 1)),
#   #           size = 5,
# #           parse = TRUE,
# #           vjust = 5,
# #           hjust = 0) +
# # geom_text(x = max(by_site$Euc) * .88,
# #           y = Inf,
# #           label = paste0("type == ", by_site_best$type),
# #           size = 5,
# #           parse = TRUE,
# #           vjust = 5.5,
# #           hjust = 0) +
# facet_wrap(. ~ Site,
#            nrow = 2,
#            scales = "free_y")
#
#
# ################################################################################
# # fit a half-normal transfer function
#
# # nls_mod = nls(BLUP ~ a * exp(-0.5 * Euc^2/b^2), start = list(a = 120, b = 2),
# #               data = val_mods_test, trace = TRUE)
#
# halfnorm = val_mods_blup %>%
#   # bind_rows((val_mods_blup %>% mutate(Class = "All"))) %>%
#   group_by(Site, Class, model_name) %>%
#   nest() %>%
#   mutate(model = purrr::map(data,
#                                    ~nls(BLUP ~ a * exp(-0.5 * Euc^2/b^2),
#                                        start = list(a = 90, b = 3), data = .))) %>%
#   mutate(preds = map2(data, model, add_predictions)) %>%
#   mutate(preds = map(preds, dplyr::select, pred)) %>%
#   mutate(resids = map2(data, model, add_residuals)) %>%
#   mutate(resids = map(resids, dplyr::select, resid)) %>%
#   mutate(glance = map(model, broom::glance)) %>%
#   unnest(c(data, preds, resids, glance)) %>%
#   mutate(legend_id = case_when(
#     Class == "Orchard" & model_name == "Census Traits" ~ "orch_meas",
#     Class == "Orchard" & model_name != "Census Traits" ~ "orch_pred",
#     Class == "Wildstand" & model_name == "Census Traits" ~ "wild_meas",
#     Class == "Wildstand" & model_name != "Census Traits" ~ "wild_pred"))
#
#
# halfnorm2 = halfnorm %>%
#   dplyr::select(model_name, AIC, Site) %>%
#   filter(model_name == "Census Traits") %>%
#   group_by(Site) %>%
#   mutate(AIC_census = AIC) %>%
#   dplyr::select(-AIC, -model_name)
#
# halfnorm3 = left_join(halfnorm, halfnorm2) %>%
#   mutate(delta_AIC = AIC - AIC_census)
#
# ################################################################################
#
# #------------------------------------------------------------------------------#
# # tables for the paper
#
#
# dbh_compare_site = val_mods_sum %>%
#   filter(grepl("site", model_name)) %>%
#   select(-rmse) %>%
#   rename("MAE" = "mae") %>%
#   #select(order(colnames(.))) %>%
#   group_by(Site) %>%
#   mutate(mae_rank = percent_rank(desc(MAE)),
#          R2_rank = percent_rank(R2)) %>%
#   group_by(model_name) %>%
#   mutate(rank_mae_median = median(mae_rank),
#          rank_mae_mean = mean(mae_rank),
#          rank_R2_median = median(R2_rank),
#          rank_R2_mean = mean(R2_rank)) %>%
#   ungroup() %>%
#   pivot_longer(contains("rank"),
#                names_to = "fit",
#                values_to = "value")
#
# s4 = dbh_compare_site %>%
#   filter(fit == "R2_rank") %>%
#   ggplot(aes(x = reorder(model_name, value, mean), y = value, color = Site, group = Site)) +
#   geom_point(size = 7) +
#   geom_line(size = 1.5) +
#   labs(title = "Rank plot comparison of predictive DBH models, by site") +
#   scale_y_reverse() +
#   theme_bw(base_size = 25) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   facet_grid(fit ~ .)
#
#
# ggsave(plot = s4,
#        filename = "D:\\Sync\\Figures\\_FIGURES_paper3\\figure_S4.tiff",
#        device = "tiff",
#        width = 15,
#        height = 10,
#        units = "in",
#        dpi = 600)
#
# # DBH models
#
# dbh_table_rank = val_mods_sum %>%
#   filter(grepl("site", model_name)) %>%
#   select(-rmse) %>%
#   rename("MAE" = "mae") %>%
#   #select(order(colnames(.))) %>%
#   group_by(Site) %>%
#   mutate(mae_rank = percent_rank(MAE),
#          R2_rank = percent_rank(desc(R2))) %>%
#   group_by(model_name) %>%
#   mutate(mae_rank_median = median(mae_rank),
#          mae_rank_mean = mean(mae_rank),
#          R2_rank_median = median(R2_rank),
#          R2_rank_mean = mean(R2_rank)) %>%
#   arrange(R2_rank_mean) %>%
#   ungroup()
#
#
#
# (dbh_table_pub = dbh_table_rank %>%
#     mutate(across(MAE:R2, round, 2)) %>%
#     dplyr::select(model_name:R2) %>%
#     pivot_longer(cols = c("MAE", "R2"), names_to = "fit_metric", values_to = "fit_value") %>%
#     pivot_wider(names_from = c("Site", "fit_metric"), values_from = "fit_value", names_glue = "{Site}_{fit_metric}") %>%
#     add_row(.before = 1) %>%
#     mutate(across(2:13, .fns = ~ ifelse(is.na(.x) == TRUE, cur_column(), .x))) %>%
#     separate_rows(2:13, sep = "_") %>%
#     flextable::flextable())
#
#
# ################################################################################
#
# halfnorm_compare_site = halfnorm3 %>%
#   dplyr::select(Site, Class, model_name, sigma:delta_AIC) %>%
#   distinct() %>%
#   filter(grepl("site", model_name) | model_name == "Census Traits") %>%
#   distinct() %>%
#   group_by(Site, Class) %>%
#   mutate(sign = if_else(delta_AIC > 0, 1, -1),
#          delta_AIC_rank = (percent_rank(abs(delta_AIC))) * sign)
#   #        #deviance_rank = dense_rank(deviance),
#   #        AIC_rank = dense_rank(AIC)) %>%
#   # ungroup() %>%
#   # mutate(model_name_new = reorder_within(model_name, AIC, Class))
# # %>%
# #   pivot_longer(contains("rank"),
# #                names_to = "fit",
# #                values_to = "value")
#
# s2 = halfnorm_compare_site %>%
#   #filter(fit == "AIC_rank") %>%
#   ggplot(aes(x = Site, y = delta_AIC, color = model_name, group = model_name)) +
#   geom_point(size = 7) +
#   geom_line(size = 1.5) +
#   labs(title = "Rank plot comparison of half-normal models, validated across dataset") +
#   scale_y_reverse() +
#   theme_bw(base_size = 20) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   facet_grid(. ~ Class, scales = "free_x")
#
#
#
# ggsave(plot = s2,
#        filename = "D:\\Sync\\Figures\\_FIGURES_paper3\\figure_S4.tiff",
#        device = "tiff",
#        width = 15,
#        height = 10,
#        units = "in",
#        dpi = 600)
#
# ggsave(plot = s1,
#        filename = "D:\\Sync\\Figures\\_FIGURES_paper3\\figure_S4.jpeg",
#        device = "jpeg",
#        width = 15,
#        height = 10,
#        units = "in",
#        dpi = 300)
#
#
#
# # halfnorm
# (halfnorm_table_pub = halfnorm %>%
#     select(-logLik) %>%
#     group_by(Site) %>%
#     dplyr::select(Site, Class, model_name, AIC, deviance) %>%
#     distinct() %>%
#     filter(grepl("site", model_name) | model_name == "Census Traits") %>%
#     group_by(Class, Site) %>%
#     mutate(dev_rank = dense_rank(deviance),
#            aic_rank = dense_rank(AIC)) %>%
#     group_by(model_name) %>%
#     mutate(rank_median = median(aic_rank),
#            rank_mean = mean(aic_rank)) %>%
#     arrange(rank_mean) %>%
#     ungroup() %>%
#     mutate(across(AIC:deviance, round, 2)) %>%
#     dplyr::select(-(dev_rank:rank_mean)) %>%
#     select(order(colnames(.))) %>%
#     pivot_longer(cols = c("AIC", "deviance"), names_to = "fit_metric", values_to = "fit_value") %>%
#     pivot_wider(names_from = c("Class", "Site", "fit_metric"), values_from = "fit_value", names_glue = "{Site}_{Class}_{fit_metric}") %>%
#     add_row(.before = 1) %>%
#     mutate(across(2:25, .fns = ~ ifelse(is.na(.x) == TRUE, cur_column(), .x))) %>%
#     separate_rows(2:25, sep = "_") %>%
#     flextable::flextable())
#
#
# save_as_docx(
#   "dbh_table" = dbh_table_pub,
#   "halfnorm table" = halfnorm_table_pub,
#   path = "D:\\Sync\\Figures\\_TABLES_paper3\\paper3_tables.docx")
#
# #------------------------------------------------------------------------------#
#
#
# halfnorm_plot = filter(halfnorm, model_name %in% c("Census Traits", "model_spec_site"))
#
# halfnorm_table = halfnorm_plot %>%
#   group_by(Site) %>%
#   # mutate(x = max(Euc),
#   #        y = max(BLUP)) %>%
#   dplyr::select(Site, Class, model_name, logLik, AIC, deviance) %>%
#   distinct() %>%
#   mutate(across(logLik:deviance, round, 2)) %>%
#   arrange(Class) %>%
#   group_by(Site) %>%
#   nest()
#
# data.tb = tibble(tb = list(halfnorm_table))
#
# plot_labs = c("orch_meas" = "",
#               "orch_pred" = "",
#               "wild_meas" = "Measured",
#               "wild_pred" = "Predicted")
#
# legend_title = " Orchard    Wildstand"
#
# # predictions vs DBH
# g7 = ggplot(val_mods_sum_plot) +
#   geom_abline(slope = 1, intercept = 0, color = 'black') +
#   geom_point(data = val_mods_plot, aes(y = .fitted, x = DBH16), size = 3, alpha = .7,
#              color = 'steelblue4',) +
#   geom_smooth(data = val_mods_plot, aes(y = .fitted, x = DBH16),
#               method = 'lm',
#               color = 'red4',
#               se = FALSE) +
#   theme_bw(base_size = 18) +
#   labs(x = 'Field-measured DBH BLUPs (cm)',
#        y = 'Model-predicted DBH values (cm)') +
#   # expression(bold(Slope~of~the~upper~red~edge~(RE[UPPER])))
#   geom_text(x = min(val_mods_plot$.fitted) +
#               ((max(val_mods_plot$.fitted) - min(val_mods_plot$.fitted)) * .08),
#             y = max(val_mods_plot$DBH16) * .9,
#             label = paste0('R^2 == ', round(val_mods_sum_plot$R2, digits = 2)),
#             size = 5.5,
#             parse = TRUE,
#             vjust = .5,
#             hjust = 0) +
#   geom_text(x = min(val_mods_plot$.fitted) +
#               ((max(val_mods_plot$.fitted) - min(val_mods_plot$.fitted)) * .08),
#             y = max(val_mods_plot$DBH16) * .75,
#             label = paste0('RMSE == ', round(val_mods_sum_plot$rmse, digits = 2)),
#             size = 5.5,
#             parse = TRUE,
#             vjust = .5,
#             hjust = 0) +
#   # geom_text(x = min(val_mods_plot$.fitted) +
#   #             ((max(val_mods_plot$.fitted) - min(val_mods_plot$.fitted)) * .08),
#   #           y = max(val_mods_plot$DBH16) * .6,
#   #           label = paste0('MAE == ', round(val_mods_sum_plot$mae, digits = 2)),
#   #           size = 5.5,
#   #           parse = TRUE,
#   #           vjust = .5,
#   #           hjust = 0) +
#   theme(text=element_text(family='Cambria', size=14, face = 'bold'),
#         legend.position = 'none',
#         strip.background = element_rect(fill = 'white'),
#         axis.text = element_text(size = 18, face = 'bold', color = 'black'),
#         plot.background = element_rect(color = 'black', fill = NA, size = 1),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_text(margin = margin(r = 10, l = -10), vjust = 1),
#         plot.margin = margin(.4, .4, .4, .4, 'cm'),
#         panel.border = element_rect(color = 'black', fill = NA, size = .5),
#         axis.title.y = element_text(size = 22, margin = margin(t = 0, r = 30, b = 0, l = 0), family = 'Cambria'),
#         axis.title.x = element_text(size = 22, margin = margin(t = 12, r = 0, b = 0, l = 0), family = 'Cambria'),
#         strip.text = element_text(size = 24)) +
#   facet_wrap(. ~ Site,
#              nrow = 2)
#
# g7
# ggsave(plot = g7,
#        filename = 'D:\\Sync\\Figures\\_FIGURES_paper3\\figure_7.tiff',
#        device = 'tiff',
#        width = 15,
#        height = 10,
#        units = 'in',
#        bg = 'white',
#        dpi = 600)
# ggsave(plot = g7,
#        filename = 'D:\\Sync\\Figures\\_FIGURES_paper3\\figure_7.tiff',
#        device = tiff,
#        width = 15,
#        height = 10,
#        units = 'in',
#        bg = 'white',
#        dpi = 600)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
