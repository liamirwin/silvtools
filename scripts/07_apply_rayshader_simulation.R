# Generate Solar Irradiance and TWI Rasters for CT Blocks
# May 2023

# ---- Project Setup ----
library(silvtools)
# List directories (each is one acquisiton of ULS/DAP)
blocks_dir <- list.dirs('G:/Quesnel_2022/process', recursive = FALSE)
# Omit these blocks from processing stream
processed <- c('CT1','CT2','CT3','CT4','CT5')
# target <- c('CT1')
blocks_dir <- blocks_dir[basename(blocks_dir) %in% processed]

# Setup SOlar Position Dataframe
start_date <- as.POSIXct("2022-05-15 00:00:00", tz = "America/Los_Angeles")
end_date <- as.POSIXct("2022-09-15 00:00:00", tz = "America/Los_Angeles")
interval <- '10 min'
interval <- '1 hour'
# 400 Road Lat Lon
lat <- 53.371759251761766
lon <- -122.76808019195623
# Block 18 Lat Lon
lat <- 48.2304193
lon <- -81.3747929

time_df <- get_solar_pos(start_date, end_date, interval, lat, lon)
# Block 18 filtering
solar_pos <- time_df %>% dplyr::filter(wday == 2 & as.numeric(hour) %in% c(6,7,8,9,10,11,12,13,14,15,16,17,18))
# Filter positions of interest
solar_pos <- time_df %>% dplyr::filter(wday == 2 & as.numeric(hour) %in% c(11,13,15))
solar_pos <- time_df %>% dplyr::filter(date == "2022-06-21" & as.numeric(hour) %in% c(9, 14))
# Processing Switches
# Generate solar irradiance rasters? - Very slow (~1 hr/timepoint/site)
generate_irr <- TRUE
# Generate TWI and topography rasters?
generate_twi <- FALSE
# Is site DAP? (not ULS)
is_dap <- FALSE
num_cores <- 1
################################################################################
# START BUTTON
################################################################################

# Loop through each directory in blocks_dir
for(i in 1:length(blocks_dir)){

  tictoc::tic() # Begin timer

  # If there's only one directory, use the first index
  if(length(blocks_dir) == 1){
    i <- 1
  }

  # Assign current directory to proj_dir
  proj_dir <- blocks_dir[i]

  # Define output directories
  raster_output <- glue::glue('{proj_dir}/output/raster')
  vector_output <- glue::glue('{proj_dir}/output/vector')

  # Check the base name of the directory for certain patterns and set flags/variables accordingly
  if(stringr::str_detect(basename(proj_dir), pattern = 'DAP') | is_dap){
    is_dap = TRUE
    acq <- paste0('DAP22_', stringr::str_replace(basename(proj_dir), pattern = "-DAP", replacement = ""))
    print(paste0('Set acqusition type as DAP named ', acq))
  } else if(stringr::str_detect(basename(proj_dir), pattern = 'MLS')){
    is_mls = TRUE
    acq <- paste0('MLS22_', stringr::str_replace(basename(proj_dir), pattern = "-MLS", replacement = ""))
    print(paste0('Set acqusition type as MLS named ', acq))
  } else{
    is_dap = FALSE
    acq <- paste0('ULS22_',basename(proj_dir))
    print(paste0('Set acqusition type as lidar (ULS) named ', acq))
  }

  # Set site_name as acq
  site_name <- acq

  # Generate irradiance if requested
  if(generate_irr == TRUE){

    # Find the required file
    dsm <- list.files(glue::glue('{proj_dir}/output/raster/dsm'), pattern = 'fill_p2r', full.names = T)

    # Call solar simulator
    solar_simulator(dsm_file = dsm,
                               proj_dir = proj_dir,
                               site_name = acq,
                               sun_pos = solar_pos,
                      num_cores = num_cores
    )

    print(glue::glue('Successfully created solar simulation for {site_name}, averaging results'))

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
    mean_irradiance(proj_dir = proj_dir, sitename = site_name)

  }

  # Generate TWI if requested
  if(generate_twi){

    print(glue::glue('Generating TWI for {site_name}'))

    # Find the required file
    dem_file <- list.files(glue::glue('{proj_dir}/output/raster/dtm'), pattern = '0.25m.tif$', full.names = T)

    if(length(dem_file) != 1){

      print(glue::glue("{length(dem_file)} DEM files matched expected input in {raster_output}/dtm"))

      for(k in 1:length(dem_file)){
        print(glue::glue('{k}.  {dem_file[k]}'))
      }

      repeat {
        user_choice <- as.integer(readline("Enter the number corresponding to the correct DEM file: "))
        if (user_choice >= 1 && user_choice <= length(dem_file)) {
          break
        } else {
          cat("Invalid choice. Please enter a number between 1 and", length(dem_file), ".\n")
        }
      }

      dem_file <- dem_file[user_choice]
      cat(glue::glue("\nSelected DEM file: {dem_file}\n"))
    }


    # Call TWI calculator function
    calc_twi(dem_file = dem_file, proj_dir = proj_dir )

  }

  tictoc::toc() # End timer after completing each directory
}


create_animation(proj_dir, sitename = acq)

proj_dir <- 'G:/Ontario_2023/Block_18/scratch/Solar_Irradiance_test/output/raster/irradiance/ULS22_Solar_Irradiance_test'
sitename <- acq
filetype <- '.tif$'
fps <- 10
label <- TRUE

animate <- function (proj_dir, fps = 10, label = TRUE, filetype = ".tif$",
          sitename)
{
  tictoc::tic()
  gif_dir <- glue::glue("{proj_dir}/gifs")
  print(glue::glue("Creating animation for {sitename}"))
  image_files <- list.files(path = proj_dir, pattern = filetype,
                            full.names = TRUE)
  image_times <- file.info(image_files)$ctime
  images_sorted <- image_files[order(image_times)]
  image_list <- lapply(images_sorted, magick::image_read)
  if (label == TRUE) {
    for (i in 1:length(image_list)) {
      file <- basename(images_sorted[i])
      date <- stringr::str_extract(file, "([0-9]{4}-[0-9]{2}-[0-9]{2})")
      hour <- stringr::str_replace(stringr::str_extract(file,
                                                        "([0-9]{2})hr"), pattern = "hr", replacement = "")
      hour <- ifelse(nchar(hour) == 1, paste0("0", hour),
                     hour)
      min <- stringr::str_replace(stringr::str_extract(file,
                                                       "([0-9]{1,2})min"), pattern = "min", replacement = "")
      min <- ifelse(nchar(min) == 1, paste0("0", min),
                    min)
      date_time_formatted <- paste(date, hour, sep = " ")
      date_time_formatted <- paste(date_time_formatted,
                                   min, sep = ":")
      date_time_formatted <- as.POSIXct(date_time_formatted,
                                        format = "%Y-%m-%d %H:%M")
      image_list[[i]] <- magick::image_annotate(image_list[[i]],
                                                date_time_formatted, color = "red", boxcolor = "white",
                                                location = "+10+10", size = 25, font = "Arial",
                                                weight = 500, gravity = "northwest")
      print(paste0(i, "/", length(image_list), " images annotated"))
    }
    print(glue::glue("Image annotation complete, creating animation with {length(image_list)} frames; this could take some time..."))
  }
  else {
    print("Label == FALSE, skipping annotation process")
  }
  animated_gif <- magick::image_join(image_list)
  img_animated = magick::image_animate(animated_gif, fps = fps)
  gif_dir <- glue::glue("{proj_dir}/gifs")
  if (!dir.exists(gif_dir)) {
    dir.create(gif_dir)
  }
  magick::image_write(image = img_animated, path = glue::glue("{gif_dir}/{sitename}_{fps}FPS.gif"))
  print(glue::glue("Successfully created and wrote irradiance gif with\n                 length(image_list) frames at {fps}FPS\n                 for {sitename} to {proj_dir}/gifs"))
  tictoc::toc()
}
