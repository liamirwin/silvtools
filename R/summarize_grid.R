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
#' @param count_trees Logical argument to count tree tops per cell or not (default is FALSE).
#' @param group_var The name of the variable to group by (e.g., "species", "crown_class", "status"). If NULL, no grouping is performed (default is NULL).
#' @param return_points_and_grid Logical argument. If TRUE, the function returns a list containing the original points with grid IDs and the grid itself. If FALSE (default), it returns the summarized grid.
#' @return An SF object with the same CRS as the input data, containing the summarized variable for each grid cell.
#' @export
#'
#' @examples
#' \dontrun{
#' # Create an example point dataset
#' toy_data <- data.frame(
#'   x = c(1, 2, 3, 4, 5, 6),
#'   y = c(1, 2, 3, 4, 5, 6),
#'   cindex = c(10, 20, 30, 40, 50, 60),
#'   species = c("A", "B", "A", "B", "A", "B")
#' )
#'
#' # Convert to an SF object
#' points_sf <- sf::st_as_sf(toy_data, coords = c("x", "y"), crs = 4326)
#'
#' # Call the summarize_grid function, grouped by species
#' result_species <- summarize_grid(points_sf, grid_area = 4, grid_shape = 'square', summary_var = 'cindex', summary_fun = 'mean', group_var = "species")
#'
#' #Call the function without grouping
#' result_no_group <- summarize_grid(points_sf, grid_area = 4, grid_shape = 'square', summary_var = 'cindex', summary_fun = 'mean')
#'
#' # Print the result
#' print(result_species)
#' print(result_no_group)
#' }
summarize_grid <- function(points_sf, grid_area = 400, grid_shape = 'square',
                           summary_var = 'cindex', summary_fun = 'mean',
                           count_trees = FALSE, group_var = NULL,
                           return_points_and_grid = FALSE) {

  tictoc::tic()

  cli::cli_alert_info("Summarizing grid data...")

  # Input Validation
  if (!inherits(points_sf, "sf")) {
    cli::cli_abort("points_sf must be an sf object.")
  }


  if(!return_points_and_grid){

  if (!summary_var %in% names(points_sf)) {
    cli::cli_abort("summary_var '{summary_var}' not found in points_sf.")
  }

  if (!(summary_fun %in% c("mean", "median", "sum", "min", "max"))) {
    cli::cli_abort("Invalid summary_fun. Must be one of: mean, median, sum, min, max")
  }

  if (!is.null(group_var) && !group_var %in% names(points_sf)) {
    cli::cli_abort("group_var '{group_var}' not found in points_sf.")
  }

  }

  # Calculate cell size
  cell_size <- switch(grid_shape,
                      square = sqrt(grid_area),
                      hexagon = sqrt((2 * grid_area) / (3 * sqrt(3))) * sqrt(3),
                      cli::cli_abort("Invalid grid_shape: must be 'square' or 'hexagon'"))


  cli::cli_alert("Creating {grid_shape} grid with cell size: {round(cell_size, 1)}m")

  # Create the grid
  if (grid_shape == 'square') {
    grid <- sf::st_make_grid(points_sf, cellsize = c(cell_size, cell_size), square = TRUE)
  } else {
    xbnds <- range(st_coordinates(points_sf)[, 1])
    ybnds <- range(st_coordinates(points_sf)[, 2])
    hex_cols <- ceiling(diff(xbnds) / cell_size)
    grid <- sf::st_make_grid(points_sf, cellsize = c(cell_size, cell_size), square = FALSE, n = c(hex_cols, NA))
  }

  grid <- grid %>% st_as_sf() %>% rowid_to_column(var = 'grid_id')


  # Calculate summaries

  # Add grid ID to points
  points_sf <- st_join(points_sf, st_as_sf(grid), join = st_within)


  if(return_points_and_grid){

    out <- list(points_sf, grid)
    names(out) <- c("points", "grid")

    return(out)
  }


  if (!is.null(group_var)) {
    cli::cli_alert(glue("Summarizing {summary_var} by {group_var} across {grid_shape} grid ({grid_area})m²..."))

    summaries <- points_sf %>%
      group_by(grid_id, !!sym(group_var)) %>% # Use !!sym to unquote the group_var
      summarize(across(all_of(summary_var), match.fun(summary_fun), .names = paste0(summary_fun, "_", summary_var)),
                n_trees = if (count_trees) n() else NULL,
                .groups = 'drop') %>% sf::st_drop_geometry()
  } else {
    cli::cli_alert(glue("Summarizing {summary_var} across {grid_shape} grid ({grid_area})m²..."))
    summaries <- points_sf %>%
      group_by(grid_id) %>%
      summarize(across(all_of(summary_var), match.fun(summary_fun), .names = paste0(summary_fun, "_", summary_var)),
                n_trees = if (count_trees) n() else NULL,
                .groups = 'drop') %>% sf::st_drop_geometry()
  }

  grid_sf <- sf::st_sf(geometry = grid)

  grid_summary <- left_join(grid_sf, summaries, by = c("grid_id" = "grid_id"))

  # Remove rows where the summary variable is NA (meaning no points in that grid cell)
  grid_summary <- grid_summary %>%
    filter(!is.na(.data[[paste0(summary_fun, "_", summary_var)]]))

  cli::cli_alert_success("Grid summarization complete.")

  tictoc::toc()

  return(grid_summary)
}
