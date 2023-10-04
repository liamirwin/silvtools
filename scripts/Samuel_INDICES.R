
library(tidyverse)
library(sf)
library(sp)
library(spatial)
library(rgdal)
library(rgeos)
library(raster)
library(terra)
library(rayshader)
library(lidR)
library(sp)
library(nngeo)
library(future)
library(rmapshaper)
library(concaveman)
library(parallel)
library(foreach)
library(doParallel)
library(smoothr)
library(ForestTools)
library(rmapshaper)
library(gdalUtilities)
library(exactextractr)
library(alphashape3d)
library(lubridate)
library(modelr)
library(lmerTest)

##---------------------------------------------------------------------------##

dir = "D:\\Sync\\_Sites\\Skimikin_spectral"

date_list = c(
  "Skimikin_2020_03_22",
  "Skimikin_2020_07_02",
  "Skimikin_2020_08_07",
  "Skimikin_2020_08_25",
  "Skimikin_2021_03_31",
  "Skimikin_2021_06_29",
  "Skimikin_2021_07_29",
  "Skimikin_2021_08_14",
  "Skimikin_2023_02_24")

#------------------------------------------------------------------------------#
# produce a hillshade mask 

library(rayshader)
library(suncalc)
library(exifr)
library(lubridate)
library(rgdal)
library(LaplacesDemon)

# make a mean values CHM
# fill NA with a minimum filter,
# smooth aggressively
# read in norm cloud
#plan(multisession, workers = 3L)

################################################################################

## Function to find index values where derivative switches from neg to post (aka local min)----------:
find_closest_to_zero_and_index <- function(values) {
  neg_slope <- numeric()
  pos_slope <- numeric()
  index_of_closest <- numeric()
  definition = numeric()
  neg_sum = numeric()
  pos_sum = numeric()
  
  k <- 1
  
  for (i in 2:length(values)) {
    #print(i)
    if (values[i - 1] < 0 & values[i] >= 0 & i > 20) { #i>15 because local min will not be in first 15 and this stops an error occuring where a local min is found in teh first 15 values and the def_positive indexing doe snot work
      # Check if the absolute value of the previous and current values is closer to zero
      # if (abs(values[i - 1]) < abs(values[i])) {
      #   closest_to_zero[k] <- values[i - 1]
      #   index_of_closest[k] <- i - 1
      # } else {
      #   closest_to_zero[k] <- values[i]
      #   index_of_closest[k] <- i
      # }
      def_positive <- c(values[i:i+15]) # vector of next 7 gradient values
      def_pos <- def_positive[def_positive > 0] #only taking positives
      pos_sum[k] <- sum(def_pos)
      
      def_negative <- c(values[(i - 10):i-15]) # vector of next 7 gradient values # 7 was used and found first min
      def_neg <- def_negative[def_negative < 0] #only taking negatives
      neg_sum[k] <- sum(def_neg)
      
      neg_slope[k] <- values[i - 1]
      pos_slope[k] <- values[i]
      index_of_closest[k] <- i - 1
      definition[k] <- sum(abs(def_neg), abs(def_pos)) # higher the value, more pronounced the local min
      k <- k + 1
    }
  }
  
  return(data.frame(Neg_Value = neg_slope, Pos_Value = pos_slope, Index = index_of_closest, definition = definition, neg_def = neg_sum, pos_def = pos_sum))
}


pols = st_read(paste0(dir, "\\output\\HULLS\\Edits\\Skimikin_z50_updated.shp")) %>%
  filter(!st_is_empty(.))

pols_spat = pols %>% 
  terra::vect()

ratios = read_rds(paste0(dir, "\\CSV\\Corrected_values//ratios.rds"))

for (x in 6:9#seq_along(date_list)
) {
  
  print(date_list[x])
  
  ms_temp = rast(paste0(dir, "\\input\\MS_ORTHO\\", date_list[x], "_MS.tif"))
  
  #make an NDVI mask
  NDVI = (ms_temp[[10]] - ms_temp[[6]]) / (ms_temp[[10]] + ms_temp[[6]])
  
  terra::writeRaster(NDVI, paste0(dir,"\\output\\RASTER\\", date_list[x], "\\",  date_list[x], "_NDVI.tif"),
                     overwrite = TRUE)
  
  NDVI_mask = NDVI
  threshold = minmax(NDVI_mask)[2] * .75
  NDVI_mask[NDVI_mask >= threshold] = NA
  NDVI_mask[NDVI_mask < threshold] = 1
  
  
  terra::writeRaster(NDVI_mask, paste0(dir,"\\output\\RASTER\\", date_list[x], "\\",  date_list[x], "_NDVI_mask.tif"),
                     overwrite = TRUE)
  
  #ratio = ratios[x,]
  
  # make a shadow mask from NIR reflectance values
  NIR = ms_temp[[10]] %>% 
    ifel(NDVI_mask > threshold, ., NA) %>% 
    crop(pols_spat) %>% 
    clamp(upper = 50000, values = FALSE) %>% 
    as.vector()
  
  #saveRDS(NIR, "D:\\Sync\\_Sites\\Skimikin_spectral\\Skimikin_2020_07_02_NIR.rds")
  
  NIR_na <- na.omit(NIR)
  density_values <- density(NIR_na)
  
  if (is.multimodal(NIR_na)){
    
    dy_dt <- pracma::gradient(density_values$y) #list of first derivatives
    
    zeros <- find_closest_to_zero_and_index(dy_dt) #list of values and their indicies where slope switches from neg to post (local min)
    zeros_filtered <- zeros[(zeros$pos_def > 0),]#filters rows with pos_def > 0, filtering out small minimums on negative slopes
    zeros_local_min <- zeros_filtered[which.max(zeros_filtered$definition), ]
    
    x_zeros <- density_values$x[zeros_local_min$Index] #selecting NIR (aka x_mid) values that correspond to index values where the slope switches from neg to positive (indicating a local min)
    #threshold <- max(x_zeros[x_zeros<0.6])
    threshold = x_zeros
    mode <- "multimodal"
    
  }else{
    mean_illum = mean(NIR, na.rm = TRUE)
    threshold = mean_illum * (2/3)
    mode <- "Unimodal"
  }
  #plot(density_values$x, density_values$y)
  
  
  
  width = 8
  height = 8
  
  (hist = NIR %>% 
      as_tibble() %>% 
      ggplot() +
      geom_histogram(aes(x = NIR), bins = 150) +
      geom_vline(xintercept = threshold, color = "red3") +
      labs(title = paste0(mode, " , threshold: ",round(threshold, digits = 2) ))+
      theme_bw()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 13,hjust = 0.75, vjust = -28)))
  
  
  ggsave(hist,
         filename = paste0(dir, "\\", date_list[x], "_NIR_shadow_hist.jpeg"),
         device = jpeg,
         width = 8,
         height = 8)
  
  saveRDS(as_tibble(NIR), 
          paste0(dir, "\\CSV\\NIR_hist\\", date_list[x], "_NIR.rds"))
  
  #hist(NIR, breaks = 100) + abline(v = threshold, col='red', lwd = 3)
  
  shadow_mask = ms_temp[[10]]
  shadow_mask[shadow_mask > threshold] = NA
  shadow_mask[shadow_mask <= threshold] = 1
  
  terra::writeRaster(shadow_mask, paste0(dir, "\\output\\RASTER\\", date_list[x], "\\",  date_list[x], "_NIR_shadow_mask.tif"),
                     overwrite = TRUE)
  
  
  shadow_patches = raster::clump(raster::raster(shadow_mask), directions = 4) %>% 
    rast()
  
  clumps = data.frame(freq(shadow_patches))
  # threshold 200 cm^2
  num_pix = 0.02 / (res(shadow_patches)[1]^2)
  flecks = clumps[clumps$count > num_pix,] #remove clump observations with frequency smaller than 9
  flecks = as.vector(flecks$value) # record IDs from clumps which met the criteria in previous step
  
  new_mask = shadow_patches %in% flecks
  new_mask[new_mask == 0] = NA
  
  terra::writeRaster(new_mask, paste0(dir, "\\output\\RASTER\\", date_list[x], "\\",  date_list[x], "_NIR_shadow_mask2.tif"),
                     overwrite = TRUE)
  
}


##---------------------------------------------------------------------------##
# extract reflectance values


pols = st_read(paste0(dir, "\\output\\HULLS\\Edits\\Skimikin_z50_updated.shp")) %>%
  filter(!st_is_empty(.))

pols_spat = pols %>% 
  terra::vect()

foreach(x = 6:9, 
        .combine = 'c',
        .packages = c("dplyr", "raster", "terra", "sf", "exactextractr", "lidR")
) %do% {
  
  ms = rast(paste0(dir, "\\input\\MS_ORTHO\\", date_list[x], "_MS.tif"))
  
  ms_mask = terra::mask(ms, pols_spat)
  
  terra::writeRaster(ms_mask, paste0(dir, "\\output\\RASTER\\", date_list[x], "\\",  date_list[x],  "_MS_mask.tif"),
                     overwrite = TRUE)
  
  mask_shadow = rast(paste0(dir,"\\output\\RASTER\\", date_list[x], "\\",  date_list[x], "_NIR_shadow_mask2.tif"))
  mask_ndvi = rast(paste0(dir,"\\output\\RASTER\\", date_list[x], "\\",  date_list[x], "_NDVI_mask.tif"))
  mask_combined = cover(mask_shadow, mask_ndvi)
  
  ms_mask_trim = trim(ms_mask)
  mask_trim = terra::crop(mask_combined, ms_mask_trim)
  
  # may not work with larger datasets. was crashing before being trimmed 
  ms_mask2 = terra::mask(ms_mask_trim, mask_trim, maskvalues = 1, updatevalue = NA)
  
  terra::writeRaster(ms_mask2, paste0(dir, "\\output\\RASTER\\", date_list[x], "\\",  date_list[x], "_MS_use.tif"),
                     overwrite = TRUE)
  
  ms_mask = rast(paste0(dir, "\\output\\RASTER\\", date_list[x], "\\",  date_list[x], "_MS_use.tif"))
  
  # ms_sum = (ms_mask[[1]] +
  #             ms_mask[[2]] +
  #             ms_mask[[3]] +
  #             ms_mask[[4]] +
  #             ms_mask[[5]] +
  #             ms_mask[[6]] +
  #             ms_mask[[7]] +
  #             ms_mask[[8]] +
  #             ms_mask[[9]] +
  #             ms_mask[[10]])
  
  rast_list = list(# reflectance values
    R444 = ms_mask[[1]] / 32752,
    R475 = ms_mask[[2]] / 32752,
    R531 = ms_mask[[3]] / 32752,
    R560 = ms_mask[[4]] / 32752,               
    R650 = ms_mask[[5]] / 32752,               
    R668 = ms_mask[[6]] / 32752,               
    R705 = ms_mask[[7]] / 32752,               
    R717 = ms_mask[[8]] / 32752,               
    R740 = ms_mask[[9]] / 32752,               
    R842 = ms_mask[[10]] / 32752,                              
    
    # chlorophyll
    NDVI = (ms_mask[[10]] - ms_mask[[6]]) / (ms_mask[[10]] + ms_mask[[6]]), 
    NIRv = ((ms_mask[[10]] - ms_mask[[6]]) / (ms_mask[[10]] + ms_mask[[6]])) * ms_mask[[10]],                          
    NDRE1 = (ms_mask[[10]] - ms_mask[[7]]) / (ms_mask[[10]] + ms_mask[[7]]),     
    #NDRE2 = (ms_mask[[10]] - ms_mask[[8]]) / (ms_mask[[10]] + ms_mask[[8]]), 
    #NDRE3 = (ms_mask[[10]] - ms_mask[[9]]) / (ms_mask[[10]] + ms_mask[[9]]), 
    EVI = 2.5 * ((ms_mask[[10]] - ms_mask[[6]]) / 
                   (ms_mask[[10]] + (6 * ms_mask[[6]]) - (7.5 * ms_mask[[1]]) + 1)),             
    Gcc = ms_mask[[4]] / (ms_mask[[2]] + ms_mask[[4]] + ms_mask[[5]]),             
    
    # carotenoids             
    SIPI = (ms_mask[[10]] - ms_mask[[1]]) / (ms_mask[[10]] - ms_mask[[6]]),             
    PRI = (ms_mask[[3]] - ms_mask[[4]]) / (ms_mask[[3]] + ms_mask[[4]]),             
    CCI = (ms_mask[[3]] - ms_mask[[5]]) / (ms_mask[[3]] + ms_mask[[5]]),             
    
    # Red edge  
    RE_lower = (ms_mask[[8]] - ms_mask[[7]]) / 12,             
    RE_total = (ms_mask[[9]] - ms_mask[[7]]) / 35, 
    RE_upper = (ms_mask[[9]] - ms_mask[[8]]) / 23)             
  
  rast_all = rast(rast_list)
  
  df_spectral = exact_extract(rast_all, pols, fun = "mean", append_cols = "Obs")
  
  saveRDS(df_spectral, paste0(dir, "\\CSV\\TRAITS\\", date_list[x], "_spectral.rds"))
}

################################################################################

# all ITC spectral values
df_spectral_all = bind_rows(
  read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2020_03_22_spectral.rds") %>% 
    mutate(Date = ymd("2020-03-22"),
           timepoint = "LW"),
  read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2020_07_02_spectral.rds") %>% 
    mutate(Date = ymd("2020-07-02"),
           timepoint = "ES"),
  read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2020_08_07_spectral.rds") %>% 
    mutate(Date = ymd("2020-08-07"),
           timepoint = "MS"),
  read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2020_08_25_spectral.rds") %>% 
    mutate(Date = ymd("2020-08-25"),
           timepoint = "LS"),
  read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2021_03_31_spectral.rds") %>% 
    mutate(Date = ymd("2021-03-31"),
           timepoint = "LW"),
  read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2021_06_29_spectral.rds") %>% 
    mutate(Date = ymd("2021-06-29"),
           timepoint = "ES"),
  read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2021_07_29_spectral.rds") %>% 
    mutate(Date = ymd("2021-07-29"),
           timepoint = "MS"),
  read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2021_08_14_spectral.rds") %>% 
    mutate(Date = ymd("2021-08-14"),
           timepoint = "LS"),
  read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2023_02_24_spectral.rds") %>% 
    mutate(Date = ymd("2023-02-24"),
           timepoint = "LW")) %>% 
  rename_with(~ gsub("mean.", "", .x, fixed = TRUE)) %>% 
  pivot_longer(R444:RE_upper, names_to = "index", values_to = "value") %>% 
  dplyr::select(-timepoint) %>% 
  mutate(value = if_else(is.nan(value), NA_real_, value)) %>% 
  pivot_wider(names_from = "Date", values_from = "value") %>% 
  drop_na() %>% 
  mutate(mean_summer = (`2020-07-02`+`2020-08-07`+
                          `2021-06-29`+`2021-07-29`)/4,
         mean_winter = (`2020-03-22`+`2021-03-31`+`2023-02-24`)/3,
         mean_LS = (`2020-08-25`+`2021-08-14`)/2,
         greenup = mean_summer - mean_winter,
         decline = mean_LS - mean_summer)

saveRDS(df_spectral_all, paste0(dir, "\\CSV\\TRAITS\\Skimikin_all_spectral.rds"))

################################################################################
# 
# # all ITC spectral values
# df_spectral_all = 
#   read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2020_03_22_spectral.rds") %>% 
#     rename_with(~ gsub("mean.", "", .x, fixed = TRUE)) %>% 
#     rename_with(~paste0(.,"_2020_LW"), -Obs) %>% 
#   left_join(read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2020_07_02_spectral.rds") %>% 
#     rename_with(~ gsub("mean.", "", .x, fixed = TRUE)) %>% 
#     rename_with(~paste0(.,"_2020_ES"), -Obs)) %>% 
#   left_join(read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2020_08_07_spectral.rds") %>% 
#     rename_with(~ gsub("mean.", "", .x, fixed = TRUE)) %>% 
#     rename_with(~paste0(.,"_2020_MS"), -Obs)) %>% 
#   left_join(read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2020_08_25_spectral.rds") %>% 
#     rename_with(~ gsub("mean.", "", .x, fixed = TRUE)) %>% 
#     rename_with(~paste0(.,"_2020_LS"), -Obs)) %>% 
#   left_join(read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2021_03_31_spectral.rds") %>% 
#     rename_with(~ gsub("mean.", "", .x, fixed = TRUE)) %>% 
#     rename_with(~paste0(.,"_2021_LW"), -Obs)) %>% 
#   left_join(read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2021_06_29_spectral.rds") %>% 
#     rename_with(~ gsub("mean.", "", .x, fixed = TRUE)) %>% 
#     rename_with(~paste0(.,"_2021_ES"), -Obs)) %>% 
#   left_join(read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2021_07_29_spectral.rds") %>% 
#     rename_with(~ gsub("mean.", "", .x, fixed = TRUE)) %>% 
#     rename_with(~paste0(.,"_2021_MS"), -Obs)) %>% 
#   left_join(read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2021_08_14_spectral.rds") %>% 
#     rename_with(~ gsub("mean.", "", .x, fixed = TRUE)) %>% 
#     rename_with(~paste0(.,"_2021_LS"), -Obs)) %>% 
#   left_join(read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2023_02_24_spectral.rds") %>% 
#     rename_with(~ gsub("mean.", "", .x, fixed = TRUE)) %>% 
#     rename_with(~paste0(.,"_2023_LW"), -Obs)) %>% 
#   
#   
# 
# saveRDS(df_spectral_all, paste0(dir, "\\CSV\\TRAITS\\Skimikin_all_spectral.rds"))

