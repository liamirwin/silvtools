#' Convert MultiPolygons to Single Polygons
#'
#' This function takes in a set of 'sf' polygon objects and separates them into
#' single polygon objects and multi-polygon objects. The function then takes
#' the largest polygon from each multi-polygon and fills any holes in the
#' resulting polygon (optional). The final output is a set of single polygon
#' objects.
#'
#' @param polygons a sf object of polygon geometries.
#' @param fill_holes a logical indicating whether to fill any holes in the
#' resulting polygon (default is TRUE).
#'
#' @return A sf object of single polygon geometries.
#'
#' @examples
#'
#' # Load example data
#' library(sf)
#' library(nngeo)
#'
#' # Load example MULTIPOLYGON tree crown dataset
#' multi_crowns <- st_read(system.file("extdata", "multi_crowns.gpkg", package = "silvtools"))
#'
#' # Apply function to convert example data to only largest single polygons
#' crowns <- convert_multi_to_single_polygons(polygons = multi_crowns, fill_holes = TRUE)
#'
#' @export
#'
convert_multi_to_single_polygons <- function(polygons, fill_holes = TRUE){
tictoc::tic()
if(!"MULTIPOLYGON" %in% unique(sf::st_geometry_type(polygons))){
  print('ERROR: Input sf polygon df contained zero MUTLIPOLYGONS; conversion not neccessary')
  stop()
}
# Seperate out multipolygons
mp <- polygons %>% dplyr::filter(sf::st_geometry_type(polygons) == 'MULTIPOLYGON')
# Seperate out polygons
sp <- polygons %>% dplyr::filter(sf::st_geometry_type(polygons) == 'POLYGON')

largest_polygon <- function(x, Z) {
  # Turn multipolygon into vector of single polygons
  x <- sf::st_combine(x)
  x <- sf::st_cast(x, "POLYGON")
  # Calculate area of each single polygon
  areas <- sf::st_area(x)
  # Take largest polygon
  max_area_index <- which.max(areas)
  # Re-add the attribute column
  largest <- x[max_area_index] %>% sf::st_as_sf(crs = sf::st_crs(x))
  # Re-name geometry column
  sf::st_geometry(largest) <- 'geometry'
  largest <- cbind(largest, Z) %>% dplyr::relocate(geometry, .after = tidyselect::last_col())
  if(fill_holes){
  # Fill any holes in the resulting polygon
  largest <- nngeo::st_remove_holes(largest)
  }
  return(largest)
}

# Go through multipolygons and apply largest_polygon function
# Empty list for fixed polygons
p_polygons <- list()
# Progress bar
pb <- progress::progress_bar$new(total = nrow(mp))
print('Taking largest polygon of each MULTIPOLYGON crown')
for (i in 1:nrow(mp)) {
  p_polygons[[i]] <- largest_polygon(sf::st_geometry(mp[i,]), Z = sf::st_drop_geometry(mp[i,]))
  pb$tick()
}
# Bind list of corrected polygons
print(glue::glue('Binding together {length(p_polygons)} cleaned polygons'))
p_polygons <- do.call(rbind, p_polygons)
sf::st_crs(p_polygons) <- sf::st_crs(polygons)
if(fill_holes == TRUE & nrow(sp) > 0){
# Fill holes in single polygon polygons
sp <- nngeo::st_remove_holes(sp)
}
# Re-join all polygons together
polygons <- rbind(sp, p_polygons)
print(glue::glue('Finished cleaning {nrow(polygons)} polygons'))
tictoc::toc()
return(polygons)
}
