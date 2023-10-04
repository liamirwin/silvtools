# Produce a filled CHM for 2020 and 2021, snap treetops to local maxima, produce
# supercells and updated segment-based polygons for further cleaning. 

library(tidyverse)
library(sf)
library(sp)
library(spatial)
library(rgdal)
library(rgeos)
library(raster)
library(terra)
library(lidR)
library(sp)
library(nngeo)
library(future)
library(rmapshaper)
library(concaveman)
library(parallel)
library(foreach)
library(smoothr)
library(ForestTools)
library(rmapshaper)
library(gdalUtilities)
library(exifr)
library(lubridate)
library(devtools)
library(exactextractr)
library(lwgeom)
library(terra)


################################################################################
devtools::source_gist('https://gist.github.com/JoshOBrien/7cf19b8b686e6d6230a78a1a9799883b')

dir = "D:\\Sync\\_Sites\\Skimikin_spectral\\"

# # snap all rasters to 2020_08_25 MS
ms_snap = rast("D:\\Sync\\_Sites\\Skimikin_spectral\\input\\MS_ORTHO\\Skimikin_2020_08_25_MS.tif")

polys_manual = st_read(paste0(dir, "GIS\\Skimikin_z50_spectral.shp")) %>% 
  st_set_crs(26910) %>% 
  st_snap_to_grid(terra::res(ms_snap)) 

# treetops to be delineated, within 40 cm of actual treetop
grid_sf = st_read(paste0(dir, "GIS\\Skimikin_grid_spectral.shp")) %>% 
  st_set_crs(26910) %>% 
  mutate(Absent = replace_na(Absent, 0)) %>% 
  dplyr::filter(Absent != 1 & Row) %>% 
  rowid_to_column()

# a buffer to clip out pixels which are not of interest,
# which are far from treetops to be delineated 
grid_buf = st_buffer(grid_sf, 3) %>% 
  st_union() %>% # unite to a geometry object
  st_sf() %>% # make the geometry a data frame object
  mutate(centrum = T) # return back the data value 

# the last one here is the one that will be used for thresholding
seasons = c("2021_CHM_max", #"2021_CHM_mean", "2021_CHM_min",
            "2020_CHM_max"#, "2020_CHM_mean", "2020_CHM_min"
            )

bound = st_read("D:\\Sync\\_Sites\\Skimikin_spectral\\GIS\\bound.shp") %>% 
  st_transform(crs = 26910)


for (i in seq_along(seasons)) {
  
  season = seasons[i]
  
  # load the smoothed CHM
  CHM_max = rast(paste0(dir, "\\output\\RASTER\\Crowns\\season_", season, ".tif")) %>% 
    terra::focal(w = 3, fun = "median", na.policy = "only", na.rm = TRUE)
  
  CHM_max_fill = CHM_max %>% 
    terra::focal(w = 3, fun = "median", na.policy = "only", na.rm = TRUE) %>% 
    resample(ms_snap, method = "bilinear",
                 filename = paste0(dir, "output\\RASTER\\Crowns\\fill_", season, ".tif"),
                 overwrite = TRUE)
    
  ###
  buffers = st_buffer(grid_sf, dist = .4) %>%
    vect()
  
  polys_buff = st_buffer(polys_manual, dist = -.10) %>% 
    vect()
  
  CHM_mask = CHM_max_fill %>% 
    mask(buffers) %>% 
    mask(polys_buff) 
  # %>% 
  #   terra::focal(w = focalWeight(., .02, "Gauss"), na.rm = FALSE, NAonly = FALSE)
  
  CHM_delin = CHM_max_fill %>%
    ifel(!is.na(CHM_mask), CHM_mask, .) %>%
    terra::focal(w = focalWeight(., .01, "Gauss"), na.policy = "omit", na.rm = TRUE, NAonly = FALSE) %>%
    ifel(!is.na(.), ., CHM_max_fill,
         filename = paste0(dir, "output\\RASTER\\Crowns\\delin_", season, ".tif"),
         overwrite = TRUE)
  
  CHM_mask2 = CHM_mask %>% 
  ifel(!is.na(.), CHM_delin, .,
       filename = paste0(dir, "output\\RASTER\\Crowns\\mask_", season, ".tif"),
              overwrite = TRUE)
  
  
  # we will search for local maxima within a buffer of these points,
  # and "snap" to the correct nearby maximum for each point
  tops = find_trees(CHM_mask2, algorithm = lmf(ws = .4,
                                              hmin = 0,
                                              shape = "circular")) %>% 
    remove.duplicates() %>% 
    st_as_sf() %>% 
    st_set_crs(26910)
  
  # here we snap the census data grid to the true local maxima
  joined1 = st_join(tops, grid_sf, join = st_nearest_feature) %>%
    group_by(rowid) %>% 
    filter(Z == max(Z))
  
  
  # save the new treetops if you want
  st_write(obj = joined1,
           dsn = paste0(dir, "output\\HULLS\\Skimikin_treetops_", season, ".shp"),
           append = FALSE)
  
  # read in the treetops as an sp object
  ttops_sp = st_read(paste0(dir, "output\\HULLS\\Skimikin_treetops_", season, ".shp")) %>% 
    as_Spatial()
}



################################################################################


library(supercells)
library(RStoolbox)

dir = "D:\\Sync\\_Sites\\Skimikin_spectral\\"

bound = grid_buf

CHM_max_2020 = rast(paste0(dir, "output\\RASTER\\\\Crowns\\delin_2020_CHM_max.tif")) %>% 
  # terra::focal(w = 3, fun = "median", na.policy = "only", na.rm = TRUE) %>% 
  # terra::focal(w = 3, fun = "median", na.policy = "only", na.rm = TRUE) %>% 
  crop(bound) %>% 
  trim()

CHM_max_2021 = rast(paste0(dir, "output\\RASTER\\\\Crowns\\delin_2021_CHM_max.tif")) %>% 
  # terra::focal(w = 3, fun = "median", na.policy = "only", na.rm = TRUE) %>% 
  # terra::focal(w = 3, fun = "median", na.policy = "only", na.rm = TRUE) %>% 
  #resample(CHM_max_2020) %>% 
  crop(bound) %>% 
  trim()

ortho = c(#rast("D:\\Sync\\_Sites\\Skimikin_spectral\\input\\MS_ORTHO\\Skimikin_2020_03_22_MS.tif")[[1:10]],
               rast("D:\\Sync\\_Sites\\Skimikin_spectral\\input\\MS_ORTHO\\Skimikin_2020_07_02_MS.tif")[[c(1,6,9)]],
               # rast("D:\\Sync\\_Sites\\Skimikin_spectral\\input\\MS_ORTHO\\Skimikin_2020_08_07_MS.tif")[[1:10]],
               # rast("D:\\Sync\\_Sites\\Skimikin_spectral\\input\\MS_ORTHO\\Skimikin_2021_07_29_MS.tif")[[1:10]],
               rast("D:\\Sync\\_Sites\\Skimikin_spectral\\input\\MS_ORTHO\\Skimikin_2021_08_14_MS.tif")[[c(1,6,9)]]) %>% 
  #resample(CHM_max_2020, method = "near") %>% 
  crop(bound) %>% 
  trim() 
# %>% 
#   terra::focal(w = focalWeight(., .01, "Gauss"), na.policy = "omit", na.rm = TRUE, NAonly = FALSE) 


# ws = rast(paste0(dir, "output\\RASTER\\Crowns\\ws_raw_", season, ".tif")) %>% 
#   #resample(CHM_max_2020) %>% 
#   mask(bbox) %>% 
#   trim()

rgb = c(rast("D:\\Sync\\_Sites\\Skimikin_spectral\\input\\ORTHO\\Skimikin_2020_07_02_RGB.tif")[[1:3]] %>% 
          resample(CHM_max_2020),
        rast("D:\\Sync\\_Sites\\Skimikin_spectral\\input\\ORTHO\\Skimikin_2021_08_14_RGB.tif")[[1:3]] %>% 
          resample(CHM_max_2020)) %>% 
  crop(bound) %>% 
  trim()
# %>% 
  # terra::focal(w = focalWeight(., .01, "Gauss"), na.policy = "omit", na.rm = TRUE, NAonly = FALSE) 

use = c(CHM_max_2020, CHM_max_2021,
        ortho, rgb)

use[use > 60000] = NA

###

# tree_rgb = brick("D:\\Sync\\_Sites\\Skimikin_spectral\\input\\ORTHO\\Skimikin_2020_07_02_RGB.tif") %>% 
#   crop(bbox)

# ttops = st_read(paste0(dir, "output\\HULLS\\Skimikin_treetops_2020_CHM_max.shp")) %>% 
#   st_crop(bbox)

use2 = terra::scale(use)
  
writeRaster(use2, paste0(dir, "output\\RASTER\\Crowns\\scaled_for_segments.tif"),
            overwrite = TRUE)

use3 = use2 %>% 
  ifel(CHM_max_2021 < .25 & CHM_max_2020 < .25, NA, .)

#plan(multisession, workers = 8L)

seg = supercells(use3, step = 6, compactness = 5, iter = 50)

#SCALE

#seg2 = dplyr::select(seg, supercells, geometry)

st_write(obj = dplyr::select(seg, supercells, geometry), 
         dsn = paste0(dir, "\\output\\HULLS\\SEGS_step6_c5.shp"),
         driver = "ESRI Shapefile",
         append = FALSE)


################################################################################
# filter out polygons that don't touch the 2020 Skimikin manual edits 

seg = st_read(paste0(dir, "\\output\\HULLS\\SEGS_step6_c5.shp"))

CHM_max_2020 = rast(paste0(dir, "output\\RASTER\\\\Crowns\\fill_2020_CHM_max.tif")) %>% 
  crop(bound) %>% 
  trim()

CHM_max_2021 = rast(paste0(dir, "output\\RASTER\\\\Crowns\\fill_2021_CHM_max.tif")) %>% 
  crop(bound) %>% 
  trim()

chm_dat = c(CHM_max_2020, CHM_max_2021)
names(chm_dat) = c("CHM_2020", "CHM_2021")

# read in the treetops as an sp object
ttops_2020 = st_read(paste0(dir, "output\\HULLS\\Edits\\Skimikin_treetops_2020_CHM_max_edits.shp")) %>% 
  vect() %>% 
  terra::extract(chm_dat[[1]], ., bind = TRUE) %>% 
  st_as_sf() %>% 
  dplyr::select(CHM_2020, geometry)

ttops_2021 = st_read(paste0(dir, "output\\HULLS\\Edits\\Skimikin_treetops_2021_CHM_max_edits.shp")) %>% 
  vect() %>% 
  terra::extract(chm_dat[[2]], ., bind = TRUE) %>% 
  st_as_sf() %>% 
  dplyr::select(CHM_2021, geometry)

polys_manual = st_read(paste0(dir, "GIS\\Skimikin_z50_spectral.shp")) %>% 
  st_set_crs(26910) %>% 
  st_join(ttops_2020) %>% 
  st_join(ttops_2021) %>% 
  dplyr::select(Obs, CHM_2020, CHM_2021, geometry) %>% 
  arrange(Obs)

seg_extract = seg %>% 
  exact_extract(x = chm_dat,
                y = .,
                fun = "max",
                force_df = TRUE,
                append_cols = c("supercells"))


seg1 = st_join(seg, polys_manual) %>% 
  drop_na() %>% 
  group_by(supercells) %>% 
  mutate(count = n()) %>% 
  rename("max_2020_poly" = "CHM_2020",
         "max_2021_poly" = "CHM_2021") %>% 
  left_join(seg_extract, by = c("supercells"))

saveRDS(seg1, paste0(dir, "\\output\\HULLS\\SEG_intersection.rds"))

seg2 = seg1 %>% 
  mutate(z50_2020 = max_2020_poly * .5,
         z50_2021 = max_2021_poly * .5) %>% 
  mutate(Obs = if_else(max.CHM_2020 > z50_2020 &
                      max.CHM_2021 > z50_2021,
                    Obs, 0)) %>% 
  group_by(supercells) %>% 
  mutate(count = n()) %>% 
  mutate(Obs = if_else(count > 1, 0, Obs))
  

st_write(obj = seg2, 
         dsn = paste0(dir, "\\output\\HULLS\\SEGS_initial.shp"),
         driver = "ESRI Shapefile",
         append = FALSE)

# seg3 = seg1 %>% 
#   mutate(z50_2020 = max_2020_poly * .5,
#          z50_2021 = max_2021_poly * .5) %>% 
#   filter(max.CHM_2020 > z50_2020 & max.CHM_2021 > z50_2021) %>% 
#   group_by(supercells) %>% 
#   mutate(count = n()) %>% 
#   filter(count < 2)
# 
# st_write(obj = seg3, 
#          dsn = paste0(dir, "\\output\\HULLS\\SEGS_filter.shp"),
#          driver = "ESRI Shapefile",
#          append = FALSE)



seg3 = seg2 %>% 
  filter(Obs != 0) %>% 
  st_centroid() %>% 
  st_join(y = polys_manual, join = st_within) %>% 
  drop_na() %>% 
  dplyr::select(supercells, geometry)

st_write(obj = seg3, 
         dsn = paste0(dir, "\\output\\HULLS\\SEGS_centroids.shp"),
         driver = "ESRI Shapefile",
         append = FALSE)

seg4 = seg2 %>% 
  st_join(seg3) %>% 
  drop_na() 

seg5 = seg4 %>% 
  group_by(Obs) %>% 
  summarise()

# seg6 = seg4 %>% 
#   filter(is.na(Obs.y))
  

st_write(obj = seg4, 
         dsn = paste0(dir, "\\output\\HULLS\\SEGS_keep.shp"),
         driver = "ESRI Shapefile",
         append = FALSE)

st_write(obj = seg5, 
         dsn = paste0(dir, "\\output\\HULLS\\SEGS_merge.shp"),
         driver = "ESRI Shapefile",
         append = FALSE)

################################################################################
# read in edited polys, merge and clean up 

ttops_2020_snap = st_read(paste0(dir, "\\output\\HULLS\\Edits\\Skimikin_treetops_2020_CHM_max_edits.shp")) %>% 
  dplyr::select(geometry)

segs_edited = st_read(paste0(dir, "\\output\\HULLS\\Edits\\SEGS_keep_edited.shp")) %>% 
  group_by(Obs) %>% 
  summarize() %>% 
  st_cast("MULTIPOLYGON") %>% 
  st_cast("POLYGON") %>% 
  st_join(st_buffer(ttops_2020_snap, .02), left = FALSE) %>% 
  drop_na() %>% 
  vect() %>% 
  fillHoles() %>% 
  st_as_sf() %>% 
  #dplyr::filter(Obs != 0 & Row != 128) %>% 
  left_join(read.csv("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\Sx_MASTER_ver_21_ADDITIONAL_TRAITS.csv") %>% 
              filter(Site == "Skim") %>% 
              dplyr::select(Meas_seq, Prov, Rep, Tree, HT16, DBH16, DBH18:COMM18),
            by = c("Obs" = "Meas_seq"))


#test = st_join(segs_edited, ttops_2020_snap, left = FALSE)

st_write(obj = segs_edited, 
         dsn = paste0(dir, "\\output\\HULLS\\Edits\\Skimikin_z50_updated.shp"),
         driver = "ESRI Shapefile",
         append = FALSE)


################################################################################


# st_write(obj = seg6, 
#          dsn = paste0(dir, "\\output\\HULLS\\SEGS_merge.shp"),
#          driver = "ESRI Shapefile",
#          append = FALSE)

# seg_inner = seg4 %>% 
#   ms_innerlines() %>% 
#   vect() %>% 
#   buffer(res(use3))
# 
# bounds = mask(CHM_max_2020, buffer(seg_inner, res(CHM_max_2020)))
# 
# use_new = use3 %>%
#   mask(seg4) %>% 
#   ifel(is.na(bounds), ., NA)
# 
# 
# #CHM_delin[!is.na(bounds)] = NA
#   
# seg_new = supercells(use_new, k = ttops_2020, step = 80, compactness = 5, iter = 5, minarea = 10)
# 
# st_write(obj = dplyr::select(seg_new, supercells, geometry), 
#          dsn = paste0(dir, "\\output\\HULLS\\SEGS_round2.shp"),
#          driver = "ESRI Shapefile",
#          append = FALSE)

################################################################################




ggRGB(tree_rgb, r = 1, g = 2, b = 3,
            scale = 255, 
            stretch = "lin",
            quantiles = c(0, 1)) +
  geom_sf(data = seg, fill = NA, color = "red3", linewidth = 1) +
  theme_void() +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = .9),
        legend.position = "none",
        legend.justification = c(0, 1),
        legend.title.align = 0,
        legend.direction = "vertical",
        legend.key.width = unit(.28, "in"),
        legend.key.height = unit(.26, "in"),
        legend.text = element_text(family = "Cambria", size = 17),
        legend.title = element_blank())


#future:::ClusterRegistry("stop")

seg_mean = seg %>% 
  dplyr::select(supercells, geometry) %>% 
  exact_extract(x = CHM_max_2020, 
                y = ., 
                fun = "mean",
                force_df = TRUE,
                append_cols = c("supercells")) %>% 
  left_join(x = seg, y = .) %>% 
  filter(max > .25) %>% 
  dplyr::select(supercells, max, geometry)

plot(seg_mean)

seg_trees = st_join(seg_mean, st_buffer(ttops, .3)) %>% 
  drop_na()

st_write(obj = seg_mean, 
         dsn = paste0(dir, "\\output\\HULLS\\SEGMENTS.shp"),
         delete_dsn = TRUE,
         driver = "ESRI Shapefile",
         append = FALSE)



################################################################################
ws_all = rast(c(paste0(dir, "output\\RASTER\\Crowns\\ws_2020_CHM_max.tif"),
                # paste0(dir, "output\\RASTER\\Crowns\\ws_2020_CHM_mean.tif"),
                # paste0(dir, "output\\RASTER\\Crowns\\ws_2020_CHM_min.tif"),
                paste0(dir, "output\\RASTER\\Crowns\\ws_2021_CHM_max.tif")
                # paste0(dir, "output\\RASTER\\Crowns\\ws_2021_CHM_mean.tif"),
                # paste0(dir, "output\\RASTER\\Crowns\\ws_2021_CHM_min.tif")
))

ws_mean = app(ws_all, "mean",
              filename = paste0(dir, "output\\RASTER\\Crowns\\ws_mean.tif"),
              overwrite = TRUE)

# ws_mode = app(ws_all, "modal",
#                 filename = paste0(dir, "output\\RASTER\\Crowns\\ws_mode.tif"),
#                 overwrite = TRUE)

CHM_vote = CHM_threshold %>% 
  ifel(ws_all[[1]] == ws_mean, ., NA,
       filename = paste0(dir, "output\\RASTER\\Crowns\\vote_CHM.tif"),
       overwrite = TRUE)


ws2 = ForestTools::mcws(CHM = raster(CHM_vote),
                        treetops = ttops_sp,
                        OSGeoPath = "C:\\OSGeo4W64",
                        format = "polygons",
                        minHeight = 0.2)

ws2_smooth = st_as_sf(ws2) %>% 
  fill_holes(threshold = units::set_units(2000, cm^2)) 

ws3 = ws2_smooth %>%
  dplyr::filter(Obs != 0 & Row != 128) %>% 
  left_join(read.csv("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\Sx_MASTER_ver_21_ADDITIONAL_TRAITS.csv") %>% 
              filter(Site == "Skim") %>% 
              dplyr::select(Meas_seq, DBH18:COMM18),
            by = c("Obs" = "Meas_seq"))

st_write(obj = ws3, 
         dsn = paste0(dir, "\\output\\HULLS\\Skimikin_z50.shp"),
         # layer = name,
         driver = "ESRI Shapefile",
         append = FALSE)



# inters1 <- st_difference(bfr2, bfr1)    
# final_inters1 <- inters1[which(inters1$id == inters1$id.1),]
