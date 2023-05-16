# List directories (each is one acquisiton of ULS/DAP)
blocks_dir <- list.dirs('G:/Quesnel_2022/process', recursive = FALSE)
# Omit these blocks from processing stream
processed <- c('CT4','CT5')
# target <- c('CT1')
blocks_dir <- blocks_dir[basename(blocks_dir) %in% processed]

################################################################################
# START BUTTON
################################################################################

is_dap <- FALSE

start_date <- as.POSIXct("2022-05-15 00:00:00", tz = "America/Los_Angeles")
end_date <- as.POSIXct("2022-09-15 00:00:00", tz = "America/Los_Angeles")
interval <- '1 hour'
lat <- 53.371759251761766
lon <- -122.76808019195623

time_df <- get_solar_pos(start_date, end_date, interval, lat, lon)

solar_pos <- time_df %>% filter(wday == 2 & as.numeric(hour) %in% c(11,13,15))

for(i in 1:length(blocks_dir)){

  tictoc::tic()

  if(length(blocks_dir) == 1){
    i <- 1
  }

  proj_dir <- blocks_dir[i]

  # ----- Output directories -----

  raster_output <- glue::glue('{proj_dir}/output/raster')
  vector_output <- glue::glue('{proj_dir}/output/vector')


  if(stringr::str_detect(basename(proj_dir), pattern = 'DAP') | is_dap){
    is_dap = TRUE
    # Set acquisition name (DAPYY_blockname)
    acq <- paste0('DAP22_', stringr::str_replace(basename(proj_dir), pattern = "-DAP", replacement = ""))
    print(paste0('Set acqusition type as DAP named ', acq))
  } else if(stringr::str_detect(basename(proj_dir), pattern = 'MLS')){
    is_mls = TRUE
    # Set acquisition name (DAPYY_blockname)
    acq <- paste0('MLS22_', stringr::str_replace(basename(proj_dir), pattern = "-MLS", replacement = ""))
    print(paste0('Set acqusition type as MLS named ', acq))
  } else{
    is_dap = FALSE
    # Set acquisition name (ULSYY_blockname)
    acq <- paste0('ULS22_',basename(proj_dir))
    print(paste0('Set acqusition type as lidar (ULS) named ', acq))
  }

    site_name <- acq

    dsm <- list.files(glue::glue('{proj_dir}/output/raster/dsm'), pattern = 'fill_p2r_0.1m.tif$', full.names = T)

    silvtools::solar_simulator(dsm_file = dsm,
                               proj_dir = proj_dir,
                               site_name = acq,
                               sun_pos = solar_pos
    )

    print(glue::glue('Successfully created solar simulation for {site_name}, averaging results'))


    mean_irradiance <- function(proj_dir, sitename){
      irradiance_dir <- paste0(proj_dir, '/output/raster/irradiance/', sitename)
      irradiance_files <- list.files(irradiance_dir, pattern = '.tif$', full.names = T)
      print(glue::glue('Found {length(irradiance_files)} irradiance rasters in folder for {sitename}'))
      r <- terra::rast(irradiance_files)
      print(glue::glue('Calculating mean value across {length(irradiance_files)} rasters'))
      r_mean <- terra::app(r, mean)
      terra::writeRaster(r_mean, paste0(proj_dir, '/output/raster/irradiance/', sitename, '_mean_irradiance.tif'), overwrite = T)
    }

    mean_irradiance(proj_dir = proj_dir, sitename = site_name)

    print(glue::glue('Generating TWI for {site_name}'))

    dem_file <- list.files(glue::glue('{proj_dir}/output/raster/dtm'), pattern = '0.25m.tif$', full.names = T)

    calc_twi(dem_file = dem_file, proj_dir = proj_dir )

}
