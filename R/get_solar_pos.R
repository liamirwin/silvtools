#Get solar positions ' Title ' ' @param start_date ' @param end_date ' @param
#interval ' @param lat ' @param lon ' ' @return ' @export '
#' @examples{
#'
#' # Start and end date for observations
#' start_date <- as.POSIXct("2022-05-15 00:00:00", tz = "America/Los_Angeles")
#' end_date <- as.POSIXct("2022-09-15 00:00:00", tz = "America/Los_Angeles")
#' Observation interval ('10 min', '30 min', 'hour')
#' interval <- '10 min'
#' Latitude and longitude of area of interest
#' lat <- 53.371759251761766
#' lon <- -122.76808019195623
#'
#' solar_pos <- get_solar_pos(start_date, end_date, interval, lat, lon)
#'
#'
# #
get_solar_pos <- function(start_date,
                          end_date,
                          interval,
                          lat,
                          lon){
  tictoc::tic()
  # Create sequence of time objects with 10 min intervals
  time_vec <- seq(from = start_date, to = end_date, by = interval)

  time_df <- data.frame(
    date_posixct = time_vec,
    year = format(time_vec, "%Y"),
    month = format(time_vec, "%m"),
    date = format(time_vec, "%Y-%m-%d"),
    day = format(time_vec, "%d"),
    week = format(time_vec, "%U"),
    hour = format(time_vec, "%H"),
    wday = lubridate::wday(time_vec),
    min = lubridate::minute(time_vec)
  )

  print(glue::glue('Getting solar positions for {nrow(time_df)} timepoints between:
                 {start_date} and {end_date} every {interval}'))

  # Get sun zenith and azimuth for each timepoint at given latitude and longitude

  solarPos = suncalc::getSunlightPosition(date = time_df$date_posixct,
                                          lat = lat,
                                          lon = lon) %>%
    mutate(alt_deg = altitude * 180 / pi,
           az_deg = 180 + azimuth * 180 / pi) %>% rename(date_posixct = date)

  time_df <- merge(time_df, solarPos, by = 'date_posixct')

  tictoc::toc()
  return(time_df)

}
