#' Approximate tree tops from a CHM
#'
#' This function approximates tree locations using a canopy height model (CHM).
#' It performs tree detection with fixed, automatic, and variable window sizes.
#'
#' @param proj_dir A string specifying the project directory.
#' @param chm_dir A string specifying the raster output directory. Default is NULL.
#' @param vector_output A string specifying the vector output directory. Default is NULL.
#' @param fixed_window A logical indicating whether to perform fixed window tree detection. Default is TRUE.
#' @param auto_window A logical indicating whether to perform automatic window tree detection. Default is TRUE.
#' @param variable_window A logical indicating whether to perform variable window tree detection. Default is TRUE.
#' @param fix_ws A numeric specifying the fixed window size for tree detection. Default is 2.
#' @param hmin A numeric specifying the minimum tree height. Default is 2m.
#' @param mod A numeric modifier for the variable window size function. Default is 0.07.
#' @param chm_ext A string related to the desired chm file selection pattern (default smooth) if multiple options; user chooses
#' @param save_output A logical indicating whether to save the output to disk. Default is FALSE.
#' @param vector_ext A string specifying the vector file extension. Default is '.gpkg', options '.shp', '.gpkg', '.geojson'.
#' @param is_vrt A logical indicating whether the CHM is a virtual raster. Default is FALSE.
#' @return A list containing the results of each tree detection method.
#'
#' @export
approximate_chm_treetops <- function(proj_dir,
                              chm_dir = NULL,
                              vector_output = NULL,
                              fixed_window = TRUE,
                              auto_window = TRUE,
                              variable_window = TRUE,
                              fix_ws = 2,
                              hmin = 2,
                              mod = 0.07,
                              chm_ext = 'smooth_p2r',
                              save_output = TRUE,
                              vector_ext = '.gpkg',
                              is_vrt = FALSE) {

  tictoc::tic()

  assertthat::assert_that(
    is.character(proj_dir), length(proj_dir) == 1,
    is.character(chm_dir) | is.null(chm_dir),
    is.character(vector_output) | is.null(vector_output),
    is.logical(fixed_window), is.logical(auto_window), is.logical(variable_window),
    is.numeric(fix_ws), is.numeric(hmin), is.numeric(mod)
  )

  acq <- basename(proj_dir)

  print(glue::glue('Beginning tree approximation for {acq} - {Sys.time()}'))

  # ----- Output directories -----
  if(is.null(chm_dir)) {
    chm_dir <- glue::glue('{proj_dir}/output/raster/chm')
  }

  # ---- Find Treetops ----

  if (!dir.exists(chm_dir)) {
    stop(glue::glue("The specified CHM directory {chm_dir} does not exist."))
  }



  if(is_vrt){
    vrt_file <- list.files(chm_dir, pattern = '.vrt$', full.names = T)
    chm <- vrt(vrt_file)
  } else{
  chm <- terra::rast(select_file_path(chm_dir, pattern = chm_ext, ext = '.tif$'))
  }

  results <- list()

  # Fixed window size
  if(fixed_window) {
    tictoc::tic()
    print(glue::glue('Locating Fixed Window Treetops ({fix_ws}m) - {Sys.time()}'))
    results$lmf_ws <- locate_trees(chm, lmf(ws = fix_ws, hmin = hmin), uniqueness = 'incremental')
    tictoc::toc()
  }

  # Lmf auto
  if(auto_window) {
    tictoc::tic()
    print(glue::glue('Locating Automatic Window Treetops - {Sys.time()}'))
    results$lmf_auto <- locate_trees(chm, lidRplugins::lmfauto(), uniqueness = 'incremental')
    tictoc::toc()
  }

  # Variable window size
  if(variable_window) {
    tictoc::tic()
    print(glue::glue('Locating Variable Window Treetops - {Sys.time()}'))
    f <- function(x) {x * mod + 1}
    results$lmf_v <- locate_trees(chm, lmf(f, hmin = hmin), uniqueness = 'incremental')
    tictoc::toc()
  }


  if (save_output == TRUE) {
    if (is.null(vector_output)) {
      vector_output <- glue::glue('{proj_dir}/output/vector')
    }

    dir.create(glue::glue('{proj_dir}/output/vector/treetops'),
               showWarnings = FALSE)

    if (fixed_window == TRUE) {
      sf::st_write(
        results$lmf_ws,
        dsn =  glue::glue(
          '{proj_dir}/output/vector/treetops/{acq}_lmf_ws_{fix_ws}m_hmin{hmin}{vector_ext}'
        ) ,
        append = FALSE
      )
      print(glue::glue(
        'Wrote {nrow(results$lmf_ws)} fixed window ({fix_ws}m) tree tops'
      ))
    }
    if (auto_window == TRUE) {
      sf::st_write(
        results$lmf_auto,
        dsn =  glue::glue('{proj_dir}/output/vector/treetops/{acq}_lmfauto_hmin{hmin}{vector_ext}'),
        append = FALSE
      )
      print(glue::glue(
        'Wrote {nrow(results$lmf_auto)} fixed window ({fix_ws}m) tree tops'
      ))
    }
    if (variable_window == TRUE) {
      sf::st_write(
        results$lmf_v,
        dsn =  glue::glue('{proj_dir}/output/vector/treetops/{acq}_lmfv_hmin{hmin}{vector_ext}'),
        append = FALSE
      )
      print(glue::glue('Wrote {nrow(results$lmf_v)} variable window tree tops'))
    }
  }

  print('Tree top detection complete...')

  return(results)
}
