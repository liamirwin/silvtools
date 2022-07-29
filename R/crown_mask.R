#' Watershed segmentation of crowns with threshold
#'
#' @param chunk lasCatalog chunk (catalog_apply) or las file loaded with readLAS
#' @param chm_res desired resolution of chm for segmentation
#' @param chm_cutoff_percent desired threshold of treetop height to mask crown
#' @param ws local maxima finding window size
#' @param uniqueness see ?lidR::locate_trees argument defines how tree IDs are determined
#' @param vis LOGICAL; choose whether or not to plot result (default = FALSE)
#'
#'
#' @export
crown_mask <- function(chunk, chm_res = 0.25, chm_cutoff_percent = 0.5, ws = 2, uniqueness = 'incremental', vis = FALSE){

  if("LAS" %in% class(las)){
    # Check if input is las instead of ctg chunk
    las <- chunk
  }
  else{
    # Standard chunk check for catalog apply
    las <- lidR::readLAS(chunk)
  }

  if (lidR::is.empty(las)) return(NULL)
  # Generate canopy height
  chm <- lidR::rasterize_canopy(las, res = chm_res , lidR::p2r(subcircle = 0.25, na.fill = lidR::knnidw()))
  # Locate treetops based on defined ws
  ttops <- lidR::locate_trees(las, algorithm = lidR::lmf(ws = ws, hmin = 2), uniqueness = uniqueness)
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
    plot(chm)
    plot(sf::st_geometry(crowns_p), add = T, border = 'blue',lty = 2)
    plot(sf::st_geometry(crowns_mask_p), add = T, border = 'red',lty = 2)
    plot(sf::st_geometry(ttops), add = T)
  }

  return(crowns_masked)
}
