#' Segment tree crowns from a CHM
#'
#' This function segments tree crowns using a canopy height model (CHM) based on the tree tops located earlier.
#'
#' @param proj_dir A string specifying the project directory.
#' @param chm_dir A string specifying the raster output directory. Default is NULL.
#' @return A list containing the results of each segmentation method.
#' @examples
#' \dontrun{
#' segment_chm_crowns("path/to/project")
#' }
#' @export
segment_chm_crowns <- function(proj_dir,
                               chm_dir = NULL,
                               vector_output = NULL,
                               raster_output = NULL,
                               crown_methods = c('fixed', 'auto', 'variable'),
                               hmin = 1,
                               crown_height_threshold = 0.25,
                               vis = FALSE,
                               chm_ext = 'smooth') {

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
  dir.create(glue::glue('{vector_output}/crowns'), showWarnings = FALSE, recursive = TRUE)
  dir.create(glue::glue('{raster_output}/crowns'), showWarnings = FALSE, recursive = TRUE)


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

    # ---- Load Tree tops ----

    # Three different lmf methods lmf(ws = 2) lmfauto( ), variable ws lmf
    ttops_ws <- st_read(str_subset(ttops_files, pattern = 'ws'), quiet = TRUE)
    ttops_auto <- st_read(str_subset(ttops_files, pattern = 'auto'), quiet = TRUE)
    ttops_v <- st_read(str_subset(ttops_files, pattern = 'lmfv'), quiet = TRUE)

    print(glue::glue('Successfully loaded all sets of treetops for {acq}
                 lmfws2 = {nrow(ttops_ws)} trees
                 lmfauto = {nrow(ttops_auto)} trees
                 lmfv = {nrow(ttops_v)}'))

    # ---- Segment Crowns ----

    print(glue::glue('Begnning Crown Segmentation for {acq}'))

    # lmfws
    crowns <- silvtools::crown_mask(chunk = chm, ttops = ttops_ws, hmin = hmin, crown_height_threshold = crown_height_threshold, vis = FALSE)
    print('fixed ws lmf crown segmentation complete...')
    # lmfauto( )
    crowns_auto <- silvtools::crown_mask(chunk = chm, ttops = ttops_auto, hmin = hmin, crown_height_threshold = crown_height_threshold, vis = FALSE)
    print('lmfauto( ) crown segmentation complete...')
    # lmfv()
    crowns_v <- silvtools::crown_mask(chunk = chm, ttops = ttops_v, hmin = hmin, crown_height_threshold = crown_height_threshold, vis = FALSE)
    print('lmfv( ) crown segmentation complete...')

    # ---- Clean Crowns ----

    # Clean crowns, take largest polygon for each treeID; fill holes within crowns

    print('Cleaning fixed window size lmf crowns...')
    # lmf ws
    crowns_p <- sf::st_as_sf(terra::as.polygons(crowns)) %>%
      convert_multi_to_single_polygons(polygons = ., fill_holes = TRUE)

    print('Cleaning lmfauto( ) crowns...')
    # lmfauto( )
    crowns_auto_p <- sf::st_as_sf(terra::as.polygons(crowns_auto)) %>%
      convert_multi_to_single_polygons(polygons = ., fill_holes = TRUE)

    print('Cleaning lmfv( ) crowns...')
    # lmfv()
    crowns_v_p <- sf::st_as_sf(terra::as.polygons(crowns_v)) %>%
      convert_multi_to_single_polygons(polygons = ., fill_holes = TRUE)

    dir.create(glue::glue('{vector_output}/crowns'), showWarnings = FALSE, recursive = T)
    dir.create(glue::glue('{raster_output}/crowns'), showWarnings = FALSE, recursive = T)

    # ---- Write Crowns ----

    print(glue::glue('Segmentation process complete; attempting to write raster and vector crowns for {acq}'))

    ct <- round(crown_height_threshold * 100)

    # Write lmf(ws = 2) crowns
    terra::writeRaster(crowns, glue::glue('{raster_output}/crowns/{acq}_lmf_ws2_watershed_crowns_ct{ct}.tif'), overwrite = T)
    sf::st_write(crowns_p, glue::glue('{vector_output}/crowns/{acq}_lmf_ws2_watershed_crowns_ct{ct}.shp'), append = FALSE)
    # Write lmfauto( ) crowns
    terra::writeRaster(crowns_auto, glue::glue('{raster_output}/crowns/{acq}_lmf_auto_watershed_crowns_ct{ct}.tif'), overwrite = T)
    sf::st_write(crowns_auto_p, glue::glue('{vector_output}/crowns/{acq}_lmf_auto_watershed_crowns_ct{ct}.shp'), append = FALSE)
    # Write lmfv( ) crowns
    terra::writeRaster(crowns_v, glue::glue('{raster_output}/crowns/{acq}_lmf_v_watershed_crowns_ct{ct}.tif'), overwrite = T)
    sf::st_write(crowns_v_p, glue::glue('{vector_output}/crowns/{acq}_lmf_v_watershed_crowns_ct{ct}.shp'), append = FALSE)

    print(glue::glue('Wrote crown rasters for {acq}'))

    tictoc::toc()


}







