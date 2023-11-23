#' Segment tree crowns from a CHM
#'
#' This function segments tree crowns using a canopy height model (CHM) based on the tree tops located earlier.
#'
#' @param proj_dir A string specifying the project directory.
#' @param chm_dir A string specifying the raster output directory. Default is NULL.
#' @return A list containing the results of each segmentation method.
#' @examples
#' \dontrun{
#' approximate_chm_crowns("path/to/project")
#' }
#' @export
approximate_chm_crowns <- function(proj_dir,
                               chm_dir = NULL,
                               vector_output = NULL,
                               raster_output = NULL,
                               crown_methods = c('fixed', 'auto', 'variable'),
                               hmin = 1,
                               crown_height_threshold = 0.25,
                               vis = FALSE,
                               chm_ext = 'smooth',
                               window_size = NULL) {

  tictoc::tic()

  assertthat::assert_that(
    is.character(proj_dir), length(proj_dir) == 1,
    is.character(chm_dir) | is.null(chm_dir),
    is.character(vector_output) | is.null(vector_output),
    is.character(raster_output) | is.null(raster_output),
    is.numeric(crown_height_threshold),
    is.logical(vis)
  )

  # ----- Output directories -----
  if(is.null(chm_dir)) {
    chm_dir <- glue::glue('{proj_dir}/output/raster/chm')
  }
  if(is.null(vector_output)) {
    vector_output <- glue::glue('{proj_dir}/output/vector')
  }
  if(is.null(raster_output)) {
    raster_output <- glue::glue('{proj_dir}/output/raster')
  }

  # Create output directories
  if(!dir.exists(vector_output)) {
    dir.create(vector_output, showWarnings = FALSE, recursive = TRUE)
  }
  if(!dir.exists(raster_output)) {
    dir.create(raster_output, showWarnings = FALSE, recursive = TRUE)
  }

    acq <- basename(proj_dir)

    print(glue::glue('Beginning tree segmentation for {acq} - {Sys.time()}'))

    # ----- Output directories -----
    if(is.null(chm_dir)) {
      chm_dir <- glue::glue('{proj_dir}/output/raster/chm')
    }


    # ----- Output directories -----

    raster_output <- glue::glue('{proj_dir}/output/raster')
    vector_output <- glue::glue('{proj_dir}/output/vector')

    # ---- Raster Watershed Segmentation ----

    if (!dir.exists(chm_dir)) {
      stop(glue::glue("The specified CHM directory {chm_dir} does not exist."))
    }

    chm <- terra::rast(silvtools::select_file_path(chm_dir, pattern = chm_ext, ext = '.tif$'))

    ttops_files <- list.files(glue::glue('{vector_output}/treetops'), pattern = '.gpkg', full.names = T)

    ct <- round(crown_height_threshold * 100)

    # ---- Segment Crowns ----

    # Create Crown Save Directory

    if(!dir.exists(glue::glue('{raster_output}/crowns'))) {
      dir.create(glue::glue('{raster_output}/crowns'), showWarnings = FALSE, recursive = TRUE)
    }

    if(!dir.exists(glue::glue('{vector_output}/crowns'))) {
      dir.create(glue::glue('{vector_output}/crowns'), showWarnings = FALSE, recursive = TRUE)
    }

    # Start with Fixed Window Size Crowns

    if('fixed' %in% crown_methods) {
      tictoc::tic()
      # If not supplied; default window size to 2m, else take supplied argument
      if (is.null(window_size)) {
        window_size <- 2
      }
      print(glue::glue('Beginning {window_size}m fixed window size watershed lmf crown segmentation for {acq} with {ct}% crown height threshold'))
      # Load Tree Tops
      ttops_ws <- sf::st_read(stringr::str_subset(ttops_files, pattern = glue::glue('lmfws{window_size}')), quiet = TRUE)
      print(glue::glue('Successfully loaded {nrow(ttops_ws)} - lmf {window_size}m window treetops for {acq}'))
      # Segment Crowns
      crowns_ws <- silvtools::crown_mask(chunk = chm, ttops = ttops_ws, hmin = hmin, crown_height_threshold = crown_height_threshold, vis = FALSE)
      # Clean up crowns
      print('Cleaning fixed window size lmf crowns...')
      crowns_ws_p <- sf::st_as_sf(terra::as.polygons(crowns_ws))
      print('Converted fixed window crowns to polygons... cleaning them now')
      crowns_ws_p <- silvtools::convert_multi_to_single_polygons(polygons = crowns_ws_p, fill_holes = TRUE)
      # Save crowns
      terra::writeRaster(crowns_ws, glue::glue('{raster_output}/crowns/{acq}_lmf_ws{window_size}_watershed_crowns_{chm_ext}_ct{ct}.tif'), overwrite = T)
      sf::st_write(crowns_ws_p, glue::glue('{vector_output}/crowns/{acq}_lmf_ws{window_size}_watershed_crowns_{chm_ext}_ct{ct}.gpkg'), append = FALSE)
      print('fixed {window_size} m ws lmf crown segmentation complete...')
      tictoc::toc()
    } else {
      print(glue::glue('Skipping fixed window size lmf crown segmentation for {acq}'))
    }

    if('auto' %in% crown_methods) {
      tictoc::tic()
      print(glue::glue('Beginning lmfauto crown segmentation for {acq} with {ct}% crown height threshold'))
      # Load Tree Tops
      ttops_auto <- sf::st_read(stringr::str_subset(ttops_files, pattern = 'auto'), quiet = TRUE)
      print(glue::glue('Successfully loaded {nrow(ttops_auto)} - lmfauto treetops for {acq}'))
      # Segment Crowns
      crowns_auto <- silvtools::crown_mask(chunk = chm, ttops = ttops_auto, hmin = hmin, crown_height_threshold = crown_height_threshold, vis = FALSE)
      # Clean up crowns
      print('Cleaning lmfauto crowns...')
      crowns_auto_p <- sf::st_as_sf(terra::as.polygons(crowns_auto))
      print('Converted lmfauto crowns to polygons... cleaning them now')
      crowns_auto_p <- silvtools::convert_multi_to_single_polygons(polygons = crowns_auto_p, fill_holes = TRUE)
      # Save crowns
      terra::writeRaster(crowns_auto, glue::glue('{raster_output}/crowns/{acq}_lmfauto_watershed_crowns_{chm_ext}_ct{ct}.tif'), overwrite = T)
      sf::st_write(crowns_auto_p, glue::glue('{vector_output}/crowns/{acq}_lmfauto_watershed_crowns_{chm_ext}_ct{ct}.gpkg'), append = FALSE)
      print('lmfauto crown segmentation complete...')
      tictoc::toc()
      } else{
        print(glue::glue('Skipping lmfauto crown segmentation for {acq}'))
      }

    if('variable' %in% crown_methods) {
      tictoc::tic()
      print(glue::glue('Beginning lmfv crown segmentation for {acq} with {ct}% crown height threshold'))
      # Load Tree Tops
      ttops_v <- sf::st_read(stringr::str_subset(ttops_files, pattern = 'lmfv'), quiet = TRUE)
      print(glue::glue('Successfully loaded {nrow(ttops_v)} - lmfv treetops for {acq}'))
      # Segment Crowns
      crowns_v <- silvtools::crown_mask(chunk = chm, ttops = ttops_v, hmin = hmin, crown_height_threshold)
      # Clean up crowns
      print('Cleaning lmfv crowns...')
      crowns_v_p <- sf::st_as_sf(terra::as.polygons(crowns_v))
      print('Converted lmfv crowns to polygons... cleaning them now')
      crowns_v_p <- silvtools::convert_multi_to_single_polygons(polygons = crowns_v_p, fill_holes = TRUE)
      # Save crowns
      terra::writeRaster(crowns_v, glue::glue('{raster_output}/crowns/{acq}_lmfv_watershed_crowns_{chm_ext}_ct{ct}.tif'), overwrite = T)
      sf::st_write(crowns_v_p, glue::glue('{vector_output}/crowns/{acq}_lmfv_watershed_crowns_{chm_ext}_ct{ct}.gpkg'), append = FALSE)
      print('lmfv crown segmentation complete...')
      tictoc::toc()
    } else {
      print(glue::glue('Skipping lmfv crown segmentation for {acq}'))
    }

  print(glue::glue('Completed Crown Segmentation for {acq}'))
  tictoc::toc()
}







