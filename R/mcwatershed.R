#' Marker Controlled Watershed Segmentation
#'
#' @param chm - Canopy height model SpatRast
#' @param treetops - Local Maxima (Markers)
#' @param th_tree - Minimum height threshold for segmented trees (Default is 2m)
#' @param tree_id - String of field refering to tree ID attribute associated with treetops
#'
#' @return Returns SpatRast of segmented trees by treeID
#' @export
#'
#' @examples
#' \dontrun{
#' mcwatershed(chm, treetops, th_tree = 2, tree_id = 'treeID')()
#' }
mcwatershed <-function(chm, treetops, th_tree = 2, tree_id = "treeID")
{
  f = function(bbox)
  {
    if (!missing(bbox)) chm <- terra::crop(chm, bbox)

    # Convert the CHM to a matrix
    Canopy <- terra::as.matrix(chm, wide = TRUE)
    mask   <- Canopy < th_tree | is.na(Canopy)
    Canopy[mask] <- 0

    cells <- terra::cellFromXY(chm, sf::st_coordinates(treetops)[, c(1, 2)])

    # check cells for NA values this indicates that ttops and chm to not overlap completely
    if (any(is.na(cells))) {
      stop("treetops and chm do not overlap completely")
    }

    # ids <- seq_along(cells)
    ids <- dplyr::pull(treetops, tree_id)

    seeds = chm
    seeds[] = 0L
    seeds[cells] = ids
    treetops = terra::as.matrix(seeds, wide = TRUE)

    Canopy <- Canopy/max(Canopy)
    Canopy <- imager::as.cimg(Canopy)
    treetops  <- imager::as.cimg(treetops)
    Crowns <- imager::watershed(treetops, Canopy)
    Crowns <- Crowns[,,1,1]
    #storage.mode(Crowns) = "integer"

    Crowns[mask] <- NA_integer_
    out <- terra::setValues(chm, Crowns)

    return(out)
  }

  f <- lidR::plugin_its(f, raster_based = TRUE)
  return(f)
}
