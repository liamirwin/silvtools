#' Watershed segmentation of crowns with threshold
#'
#' @param chunk lasCatalog chunk (catalog_apply) or las file loaded with readLAS or CHM loaded as SpatRaster
#' @param ttops optional treetops argument, if provided lmf process is skipped
#' @param chm_res desired resolution of chm for segmentation
#' @param chm_cutoff_percent desired threshold of treetop height to mask crown
#' @param ws local maxima finding window size
#' @param uniqueness see ?lidR::locate_trees argument defines how tree IDs are determined
#' @param vis LOGICAL; choose whether or not to plot result (default = FALSE)
#'
#'
#' @export
crown_mask <- function(chunk, ttops = NULL, chm_res = 0.25, chm_cutoff_percent = 0.5, ws = 2, uniqueness = 'incremental', vis = FALSE){


  if("SpatRaster" == class(chunk)){
    chm <- chunk
    chm_res <- terra::res(chm)[1]
    print('Input is SpatRaster (Presuming CHM) proceeding with crown segmentation')
  } else{
  if("LAS" %in% class(chunk)){
    # Check if input is las instead of ctg chunk
    las <- chunk
    print('Input is las file proceeding with crown segmentation')
  }
  else{
    # Standard chunk check for catalog apply
    las <- lidR::readLAS(chunk)
    print('Input is catalog tile proceeding with crown segmentation')
  }
  if (lidR::is.empty(las)) return(NULL)
    # Generate canopy height
    chm <- lidR::rasterize_canopy(las, res = chm_res , lidR::p2r(na.fill = lidR::knnidw()))
  }


  if (is.null(ttops) == TRUE){
  # Locate treetops based on defined ws
  ttops <- lidR::locate_trees(las, algorithm = lidR::lmf(ws = ws, hmin = 2), uniqueness = uniqueness)
  }
  crowns <- silvtools::mcwatershed(chm, ttops, th_tree = 2, tree_id = 'treeID')()
  # Replace crown IDs with maximum canopy height within crown
  crowns_max <- terra::classify(x = crowns,
                                rcl = as.data.frame(terra::zonal(chm,
                                                                 crowns,
                                                                 fun = 'max')))
  # Make CHM mask raster where height is less than threshold percent
  chm_mask <- chm
  chm_mask[chm_mask < (chm_cutoff_percent * crowns_max)] = NA
  # Re-run segmentation on masked CHM
  crowns_masked <- silvtools::mcwatershed(chm_mask, ttops, th_tree = 2, tree_id = 'treeID')()


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
