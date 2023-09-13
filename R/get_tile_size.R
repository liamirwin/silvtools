#' Get Tile Size from a LAScatalog
#'
#' This internal function calculates the size of the middle tile in a LAScatalog.
#' The size is rounded to the nearest whole number.
#'
#' @param ctg A LAScatalog object.
#'
#' @return The size of the middle tile, rounded to the nearest whole number.
#'
#' @keywords internal
#' @importFrom lidR readLASheader
get_tile_size <- function(ctg) {
  tile_count <- length(ctg$filename)
  middle_tile_index <- ceiling(tile_count / 2)
  middle_tile <- ctg$filename[middle_tile_index]

  hdr <- lidR::readLASheader(middle_tile)
  phb <- hdr@PHB
  x_size <- phb$`Max X` - phb$`Min X`
  y_size <- phb$`Max Y` - phb$`Min Y`

  return(round(max(x_size, y_size)))
}
