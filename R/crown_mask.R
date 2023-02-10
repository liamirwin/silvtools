#' Watershed segmentation of crowns with threshold
#'
#' Height limited watershed segmentation; takes a las/chunk/chm and performs watershed segmentation
#' if treetops are provided it will use those; function masks crowns by crown_height_threshold defined as
#' percentage of treetop height; this will determine if pixels are included in final crown (eg if CHM pixel less than 50% of tree height; exclude)
#'
#' @param chunk lasCatalog chunk (catalog_apply) or las file loaded with readLAS or CHM loaded as SpatRaster
#' @param ttops optional treetops argument, if provided lmf process is skipped
#' @param chm_res desired resolution of chm for segmentation
#' @param crown_height_threshold desired threshold of treetop height to mask crown
#' @param ws local maxima finding window size
#' @param hmin minimum height of local maxima to identify as treetops
#' @param uniqueness see ?lidR::locate_trees argument defines how tree IDs are determined
#' @param vis LOGICAL; choose whether or not to plot result (default = FALSE)
#'
#'
#' @export
crown_mask <- function(chunk, ttops = NULL, chm_res = 0.25, crown_height_threshold = 0.5, ws = 2, hmin = 2, uniqueness = 'incremental', vis = FALSE){


  # Check that 'chunk' is not missing or NULL
  if (is.null(chunk) || missing(chunk)) {
    stop("Error: 'chunk' argument must be specified.")
  }

  # Check if input is a SpatRaster object
  if("SpatRaster" == class(chunk)){
    chm <- chunk
    chm_res <- terra::res(chm)[1]
    print('Input is SpatRaster (Presuming CHM) proceeding with crown segmentation')
  } else{
    # Check if input is a las file
    if("LAS" %in% class(chunk)){
      las <- chunk
      print('Input is las file proceeding with crown segmentation')
    }
    else{
      # Check if input is a las catalog chunk
      las <- lidR::readLAS(chunk)
      if (lidR::is.empty(las)) {
        return(NULL)
      }
      print('Input is catalog tile proceeding with crown segmentation')
    }

    # Generate canopy height
    chm <- lidR::rasterize_canopy(las, res = chm_res , lidR::p2r(na.fill = lidR::knnidw()))
  }

  # Check that 'ttops' is not NULL
  if (is.null(ttops)){
    # Locate treetops based on defined ws
    ttops <- lidR::locate_trees(las, algorithm = lidR::lmf(ws = ws, hmin = hmin), uniqueness = uniqueness)
  }


  # Initial crown segmentation

  crowns <- silvtools::mcwatershed(chm, ttops, th_tree = hmin, tree_id = 'treeID')()
  # Replace crown IDs with maximum canopy height within crown
  crowns_max <- terra::classify(x = crowns,
                                rcl = as.data.frame(terra::zonal(chm,
                                                                 crowns,
                                                                 fun = 'max')))
  # Make CHM mask raster where height is less than threshold percent
  chm_mask <- chm
  chm_mask[chm_mask < (crown_height_threshold * crowns_max)] = NA
  # Re-run segmentation on masked CHM
  crowns_masked <- silvtools::mcwatershed(chm_mask, ttops, th_tree = hmin, tree_id = 'treeID')()


  if(vis == TRUE){
    # Crown polygons
    crowns_p <- sf::st_as_sf(terra::as.polygons(crowns))
    crowns_mask_p <- sf::st_as_sf(terra::as.polygons(crowns_masked))
    # Plotting (Optional)
    terra::plot(chm)
    plot(sf::st_geometry(crowns_p), add = T, border = 'blue',lty = 2)
    plot(sf::st_geometry(crowns_mask_p), add = T, border = 'red',lty = 2)
    plot(sf::st_geometry(ttops), add = T)
  }

  return(crowns_masked)
}
