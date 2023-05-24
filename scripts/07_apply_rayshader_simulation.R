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
interval <- '1 hour'
lat <- 53.371759251761766
lon <- -122.76808019195623
time_df <- get_solar_pos(start_date, end_date, interval, lat, lon)
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
    dsm <- list.files(glue::glue('{proj_dir}/output/raster/dsm'), pattern = 'fill_p2r_0.1m.tif$', full.names = T)

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

