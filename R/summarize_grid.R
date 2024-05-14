#' Summarize an SF object with tree top locations across a regular grid
#'
#' This function takes an SF object with tree top locations and summarizes a given variable
#' across a regular grid (either square or hexagonal) of a specified area.
#'
#' @param points_sf An SF object containing the tree top locations with at least one numeric variable for summarization.
#' @param grid_area The area of the grid cells, in square meters (default is 100 m²).
#' @param grid_shape The shape of the grid, either 'square' or 'hexagon' (default is 'square').
#' @param summary_var The name of the variable to summarize (default is 'cindex').
#' @param summary_fun The summary function to apply (default is 'mean').
#' @param count_trees Logical argument to count tree tops per cell or not (default is false)
#' @param by_species Logical argument; if TRUE, summarization will be grouped by a cateogrical column: 'species' result will be summary variable for each species
#' @return An SF object with the same CRS as the input data, containing the summarized variable for each grid cell.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Create an example point dataset
#' toy_data <- data.frame(
#'   x = c(1, 2, 3, 4, 5, 6),
#'   y = c(1, 2, 3, 4, 5, 6),
#'   cindex = c(10, 20, 30, 40, 50, 60)
#' )
#'
#' # Convert to an SF object
#' points_sf <- sf::st_as_sf(toy_data, coords = c("x", "y"), crs = 4326)
#'
#' # Call the summarize_grid function
#' result <- summarize_grid(points_sf, grid_area = 4, grid_shape = 'square', summary_var = 'cindex', summary_fun = 'mean')
#'
#' # Print the result
#' print(result)
#' }
summarize_grid <- function(points_sf, grid_area = 400, grid_shape = 'square',
                           summary_var = 'cindex', summary_fun = 'mean',
                           count_trees = FALSE, by_species = FALSE) {

  # Set cell size based on grid_area and grid_shape
  if (grid_shape == 'square') {
    cell_size <- sqrt(grid_area)
  } else if (grid_shape == 'hexagon') {
    # Calculate the edge length based on the given area
    edge_length <- sqrt((2 * grid_area) / (3 * sqrt(3)))

    # Calculate the cell size (distance between opposite edges)
    cell_size <- edge_length * sqrt(3)
  } else {
    stop("Invalid grid_shape: must be 'square' or 'hexagon'")
  }

  # Create the grid
  if (grid_shape == 'square') {
    # Create square grid
    grid <- sf::st_make_grid(points_sf, cellsize = c(cell_size, cell_size), square = TRUE)
  } else {
    # Calculate the number of hexagons per row
    xbnds <- range(st_coordinates(points_sf)[, 1])
    ybnds <- range(st_coordinates(points_sf)[, 2])
    hex_cols <- ceiling(diff(xbnds) / cell_size)

    # Create hexagonal grid
    grid <- sf::st_make_grid(points_sf, cellsize = c(cell_size, cell_size), square = FALSE, n = c(hex_cols, NA))
  }

  # Calculate summaries
  # Assign grid cell id to each point
  points_sf$grid_id <- as.character(st_within(points_sf, grid))

  
  
  if(by_species){
    # Calculate summary metric for each species group
    summaries <- points_sf %>%
    group_by(grid_id, species) %>%
    summarize(across(all_of(summary_var), match.fun(summary_fun), .names = paste0(summary_fun, "_", summary_var)),
              n_trees = if(count_trees) n() else NULL,
              .groups = 'drop') %>% sf::st_drop_geometry()
  } else{
  # Group points by grid_id and calculate summary statistics
  summaries <- points_sf %>%
    group_by(grid_id) %>%
    summarize(across(all_of(summary_var), match.fun(summary_fun), .names = paste0(summary_fun, "_", summary_var)),
              n_trees = if(count_trees) n() else NULL,
              .groups = 'drop') %>% sf::st_drop_geometry()
  }

  # Convert grid to sf data frame and add id column
  grid_sf <- sf::st_sf(geometry = grid)
  grid_sf$id <- as.character(seq_len(nrow(grid_sf)))

  # Merge summaries back into the grid
  grid_summary <- left_join(grid_sf, summaries, by = c("id" = "grid_id"))

  # Remove grid cells with no points within them
  grid_summary <- grid_summary[!is.na(grid_summary[[paste0(summary_fun, "_", summary_var)]]), ]

  return(grid_summary)
}
