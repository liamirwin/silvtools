# Rayshading solar simulator function
#' Title
#'
#' @param dsm_file
#' @param weekday
#' @param time_df
#' @param proj_dir
#' @param site_name
#' @param aoi_file
#' @param sun_pos
#' @param interval
#' @param start_time
#' @param end_time
#' @param filt
#'
#' @return
#' @export
#'
#' @examples
solar_simulator <- function(dsm_file,
                            weekday = 'Mon',
                            time_df,
                            proj_dir,
                            site_name,
                            aoi_file = NULL,
                            sun_pos,
                            interval,
                            start_time,
                            end_time,
                            filt = FALSE){

  sim_start <- Sys.time()


  # Filter time dataframe for target day of the week/intervals

  filter_by_interval <- function(data,
                                 interval, # Options: "10min", "30min", "hour"
                                 weekday, # Options: "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday","Sunday"
                                 start_time, # Starting hour (1:24)
                                 end_time # Ending hour (1:24)
  ) {
    tictoc::tic()

    if (interval %in% c("10min", "30min", "60min")) {
      if (weekday %in% c("Mon",
                         "Tue",
                         "Wed",
                         "Thu",
                         "Fri",
                         "Sat",
                         "Sun")) {
        if (interval == "10min") {
          filtered_data <-
            data %>% filter(
              as.numeric(hour) %in% c(start_time:end_time),
              as.numeric(min) %% 10 == 0,
              weekday == lubridate::wday(date, label = TRUE)
            )
        } else if (interval == "30min") {
          filtered_data <- data %>% filter(
            as.numeric(hour) %in% c(start_time:end_time),
            as.numeric(min) %in% c(0, 30),
            weekday == lubridate::wday(date, label = TRUE)
          )
        } else if (interval == "60min") {
          filtered_data <- data %>% filter(
            as.numeric(min) == 0,
            as.numeric(hour) %in% c(start_time:end_time),
            weekday == lubridate::wday(date, label = TRUE)
          )
        }

      } else {
        stop("Invalid weekday, please choose from: Mon, Tue, Wed, Thu, Fri, Sat, Sun")
      }
    }
    print(
      glue::glue(
        'Filtered out {nrow(data) - nrow(filtered_data)} ({round((nrow(data) - nrow(filtered_data))/nrow(data)*100)}%)
                   timepoints, keeping every solar positions for {weekday} every {interval} interval'
      )
    )
    tictoc::toc()
    return(filtered_data)
  }

  if (filt == TRUE) {
    filtered_times <- filter_by_interval(
      data = sun_pos,
      interval = interval,
      weekday = weekday,
      start_time = start_time,
      end_time = end_time
    )
    print('hoho')
  } else {
    filtered_times <- sun_pos
    print('Filter == FALSE so timepoints not filtered')
  }
  # Load DSM raster and crop if aoi provided (useful for testing/benchmarking)
  if(!is.null(aoi_file)) {
    aoi <- st_read(aoi_file)
    dsm <- terra::rast(dsm_file) %>% terra::crop(aoi)
  } else{
    dsm <- terra::rast(dsm_file)
  }
  # Convert to matrix for rayshader
  rasM <- rayshader::raster_to_matrix(dsm)

  # Create output directory for irradiance models if it doesnt exist
  if(!dir.exists(glue::glue('{proj_dir}/output/raster/irradiance'))) {
    dir.create(glue::glue('{proj_dir}/output/raster/irradiance'),
               recursive = T)
  }

  print(glue::glue('Beginning rayshading process of {site_name} DSM for {nrow(filtered_times)} timepoitns at {Sys.time()}'))

  for(k in 1:nrow(filtered_times)){
    tictoc::tic()
    timepoint <- filtered_times[k,]
    zenith <- timepoint$alt_deg
    azimuth <- timepoint$az_deg

    print(glue::glue('Beginning rayshading for timepoint {timepoint$date_posixct} ({k}/{nrow(filtered_times)}) for {site_name}'))

    rshade= rayshader::ray_shade(
      rasM,
      sunaltitude = zenith,
      sunangle = azimuth,
      zscale = 0.05,
      progbar = interactive())

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
        '{proj_dir}/output/raster/irradiance/{site_name}_{timepoint$date}-{timepoint$hour}hr-{timepoint$min}min_rayshade.tif'
      ),
      overwrite = T
    )
    tictoc::toc()

  }
}
