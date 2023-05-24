#' Create Animation
#' Animates solar irradiance tifs (or any images) into gif and saves locally
#' @param proj_dir Directory where rasters from solar_simulator are stored (.tifs), creates a /gifs output folder for saving
#' @param fps Chosen frames per second for final gif; change depending on number of images default is 10 fps
#' @param label Logical; TRUE/FALSE. Determines if final animation is annotated with YYYY-MM-DD HH:MM:SS time stamp for each frame
#' @param filetype File type of animation images; default assumes raster in .tif format
#' @param sitename string name of site ex 'CT1P1'
#'
#' @return Saves a .gif file locally
#' @export
#'
#' @examples
#' \dontrun{
#' proj_dir <- 'H:/Quesnel_2022/blocks/CT1-DAP'
#' create_animation(proj_dir, sitename = 'CT1P1')
#' }
create_animation <- function(proj_dir,
                             fps = 10,
                             label = TRUE,
                             filetype = '.tif$',
                             sitename){
  tictoc::tic()
  # Directory should be where irradiance tifs are stored
  gif_dir <- glue::glue('{proj_dir}/gifs')
  print(glue::glue('Creating animation for {sitename}'))
  # Get list of images
  image_files <- list.files(path = proj_dir, pattern = filetype, full.names = TRUE)
  # Get creation times of image files
  image_times <- file.info(image_files)$ctime
  # Sort image files by creation time
  images_sorted <- image_files[order(image_times)]
  # Read in images with magick in order of creation
  image_list <- lapply(images_sorted, magick::image_read)
  if(label == TRUE){
    # Apply date and time label to each image
    for (i in 1:length(image_list)){
      file <- basename(images_sorted[i])
      # Extract date and time from filename
      date <- stringr::str_extract(file, "([0-9]{4}-[0-9]{2}-[0-9]{2})")
      hour <- stringr::str_replace(stringr::str_extract(file, "([0-9]{2})hr"), pattern = 'hr', replacement = "")
      hour <- ifelse(nchar(hour)==1, paste0("0",hour), hour)
      min <- stringr::str_replace(stringr::str_extract(file, "([0-9]{1,2})min"), pattern = 'min', replacement = "")
      min <- ifelse(nchar(min)==1, paste0("0",min), min)
      # Format date, hour, minute into YYYY-MM-DD HH:MM:SS TZ format
      date_time_formatted <- paste(date,hour,sep = " ")
      date_time_formatted <- paste(date_time_formatted,min,sep = ":")
      date_time_formatted <- as.POSIXct(date_time_formatted, format = "%Y-%m-%d %H:%M")
      # Annotate image with date and time label
      image_list[[i]] <- magick::image_annotate(image_list[[i]], date_time_formatted, color = 'red', boxcolor = 'white',
                                                location = "+10+10", size = 25, font = "Arial", weight = 500,
                                                gravity = "northwest")
      print(paste0(i, '/', length(image_list), ' images annotated'))
    }
    print(glue::glue('Image annotation complete, creating animation with {length(image_list)} frames; this could take some time...'))
  } else{
    print('Label == FALSE, skipping annotation process')
  }


  # Create an animated GIF with target frames per second
  animated_gif <- magick::image_join(image_list)
  img_animated = magick::image_animate(animated_gif, fps = fps)
  # Output directory for gif
  gif_dir <- glue::glue('{proj_dir}/gifs')
  if(!dir.exists(gif_dir)){dir.create(gif_dir)}
  # Write gif, lengthly process...
  magick::image_write(image = img_animated,
                      path = glue::glue('{gif_dir}/{sitename}_{fps}FPS.gif'))
  print(glue::glue('Successfully created and wrote irradiance gif with
                 length(image_list) frames at {fps}FPS
                 for {sitename} to {proj_dir}/gifs'))
  tictoc::toc()
}
