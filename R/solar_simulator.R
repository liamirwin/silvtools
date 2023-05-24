#' Rayshading solar simulator function
#'
#' @param dsm_file File path to DSM file of interest typically 'path/dsm.tif'
#' @param proj_dir Directory for block/acquisition creates an /output/irradiance folder within proj_dir
#' @param site_name Site name for file naming purposes ex 'CT1'
#' @param aoi_file Optional area of interest for internal cropping of DSM either file pointing to shp/gpkg etc or loaded sf object
#' @param sun_pos Sun position data frame from get_sun_pos, zenith, azimuth for each timepoint of interest
#'
#' @return Writes a raster for each timepoint recording the irradiance for each given pixel
#' @export
#'
#' @examples
#' \dontrun{
#'
#' dsm_file <- 'H:/Quesnel_2022/blocks/CT1-DAP/output/raster/dsm/CT1-DAP_dsm_fill_p2r_0.05m.tif'
#' proj_dir <- 'H:/Quesnel_2022/blocks/CT1-DAP'
#' site_name <- 'CT1-DAP'
#' aoi_file <- 'H:/Quesnel_2022/blocks/CT1-DAP/input/vector/CT1P1_buffer.shp'
#'
#' sun_pos <- get_solar_pos(start_date, end_date, interval, lat, lon)
#' solar_simulator(dsm_file, proj_dir, site_name, aoi_file, sun_pos)
#' }
#'
solar_simulator <- function(dsm_file,
                            proj_dir,
                            site_name,
                            aoi_file = NULL,
                            sun_pos,
                            num_cores = NULL){
  # Determine Simulation Start Time
  sim_start <- Sys.time()
  # Load DSM raster and crop if aoi provided (useful for testing/benchmarking)
  if (!is.null(aoi_file)) {
    if ('sf' %in% class(aoi_file)) {
      aoi <- aoi_file
    } else{
      aoi <- sf::st_read(aoi_file)
    }
    dsm <- terra::rast(dsm_file) %>% terra::crop(terra::vect(aoi)) %>% raster::raster(.)
    print('DSM cropped to provided area of interest')
  } else{
    dsm <- raster::raster(dsm_file)
    print('DSM raster successfully loaded')
  }
  # Convert to matrix for rayshader
  rasM <- rayshader::raster_to_matrix(dsm)
  # Create output directory for irradiance models if it doesnt exist
  if(!dir.exists(glue::glue('{proj_dir}/output/raster/irradiance/{site_name}'))) {
    dir.create(glue::glue('{proj_dir}/output/raster/irradiance/{site_name}'),
               recursive = T)
    print(glue::glue('Created directory for saving irradiance files /output/irradiance/{site_name}'))
  }

  print(glue::glue('Beginning rayshading process of {site_name} DSM for {nrow(sun_pos)} timepoints at {Sys.time()}'))

  for(k in 1:nrow(sun_pos)){
    tictoc::tic()
    timepoint <- sun_pos[k,]
    zenith <- timepoint$alt_deg
    azimuth <- timepoint$az_deg

    print(glue::glue('Beginning rayshading for timepoint {timepoint$date_posixct} ({k}/{nrow(sun_pos)}) for {site_name}'))

    if(is.null(num_cores)){
      num_cores = 1
    }

    if(num_cores > 1){
      run_parallel = TRUE
      options(cores = num_cores)
      print(glue::glue('{num_cores} enabled for parallel processing'))
    } else{
      run_parallel = FALSE
    }


    if(num_cores > 1){
      run_parallel = TRUE
    } else{
      run_parallel = FALSE
    }

    zscale = terra::res(dsm)[1]

    rshade= rayshader::ray_shade(
      rasM,
      sunaltitude = zenith,
      sunangle = azimuth,
      zscale = zscale,
      progbar = FALSE,
      multicore = run_parallel)

    rshade <- t(rshade)
    rshade<-rshade[,c(ncol(rshade):1),drop = FALSE]

    ras2 = raster::raster(rshade)
    ras2@file <- raster::raster(dsm)@file # Change metadata
    ras2@legend <- raster::raster(dsm)@legend
    ras2@extent <- raster::raster(dsm)@extent
    ras2@rotation <- raster::raster(dsm)@rotation
    ras2@crs <- raster::raster(dsm)@crs
    raster::writeRaster(
      ras2,
      glue::glue(
        '{proj_dir}/output/raster/irradiance/{site_name}/{site_name}_{timepoint$date}-{timepoint$hour}hr-{timepoint$min}min_rayshade.tif'
      ),
      overwrite = T
    )
    tictoc::toc()
  }
}
