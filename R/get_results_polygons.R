#' Retrieve RESULTS Polygons Intersecting a Given Extent
#'
#' This function queries the BC Data Catalogue for RESULTS polygons that
#' intersect a specified spatial extent. It performs CRS validation and
#' transformation, retrieves data using `bcdata`, intersects the data with
#' the extent, and optionally plots the resulting polygons using `ggplot2`.
#'
#' @param extent An `sf` object representing the spatial extent for the query.
#'   It must be a simple feature or a simple feature collection, and will be
#'   converted to a simple bounding box polygon. If the CRS is not BC Albers
#'   (EPSG:3005), it will be transformed.
#' @param plot Logical, whether to plot the resulting polygons. Default is
#'   `TRUE`.
#' @param symbology_field Character, the field to use for symbology in the plot.
#'   If `NULL` (default), the function will attempt to use "FOREST_COVER_ID" if
#'   it exists in the data. If not provided and "FOREST_COVER_ID" is not found,
#'   an error is thrown. This field is converted to a factor for plotting.
#'
#' @return An `sf` object containing the RESULTS polygons that intersect the
#'   given extent. If no polygons are found, an empty `sf` object is returned.
#'   If `plot` is `TRUE`, a `ggplot` object is also displayed, showing the
#'   polygons filled according to the `symbology_field`.
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # Example extent (replace with your actual extent)
#' extent <- st_read("path/to/your/extent.gpkg")
#'
#' # Retrieve and plot RESULTS polygons
#' results <- get_results_polygons(extent)
#'
#' # Retrieve without plotting, using a custom symbology field
#' results_no_plot <- get_results_polygons(extent, plot = FALSE, symbology_field = "MY_FIELD")
#' }
#'
#' @details
#' The function performs several steps:
#' 1. Converts the input `extent` to its bounding box.
#' 2. Validates and transforms the CRS of the extent to BC Albers (EPSG:3005).
#' 3. Queries the BC Data Catalogue for RESULTS polygons intersecting the extent.
#' 4. Performs a spatial intersection between the retrieved polygons and the extent.
#' 5. Optionally plots the resulting polygons using `ggplot2`, with symbology
#'   based on the specified field.
#'
#' The function uses `cli` for console output, providing informative messages,
#' warnings, and errors during execution.
#' @export
get_results_polygons <- function(extent, plot = TRUE, symbology_field = NULL) {

  # Check and install bcdata if necessary
  if (!requireNamespace("bcdata", quietly = TRUE)) {
    cli::cli_alert_warning("Package 'bcdata' is not installed. Installing...")
    tryCatch({
      install.packages("bcdata")
      if (!requireNamespace("bcdata", quietly = TRUE)) { # Double-check after attempt
        cli::cli_abort("Failed to install 'bcdata'. Please install it manually.")
      }
      cli::cli_alert_success("Package 'bcdata' installed successfully.")
    }, error = function(e) {
      cli::cli_abort("Error installing 'bcdata': {e$message}")
    })
  }


  cli::cli_h1("Retrieving RESULTS Polygons")

  # Ensure extent is a simple bounding box polygon
  extent <- st_as_sfc(st_bbox(extent))

  # Input Validation and CRS Transformation
  cli::cli_alert_info("Checking and transforming CRS...")
  tryCatch({
    if (st_crs(extent) != st_crs(3005)) {
      extent <- st_transform(extent, 3005)
      cli::cli_alert_warning("CRS was not BC Albers (3005), transformed to match.")
    }
    cli::cli_alert_success("CRS check and transformation complete.")
  }, error = function(e) {
    cli::cli_abort("Error transforming CRS: {e$message}")
  })


  # Data Retrieval with bcdata and Error Handling
  cli_alert_info("Querying bcdata...")
  tryCatch({
    results_polys <-
      bcdata::bcdc_query_geodata("258bb088-4113-47b1-b568-ce20bd64e3e3") %>%
      bcdata::filter(INTERSECTS(extent)) %>%
      bcdata::collect()
    cli_alert_success("Data retrieved from bcdata.") # Tick after successful operation
  }, error = function(e) {
    cli_abort("Error retrieving data from bcdata: {e$message}")
  })


  # Intersection and Casting
  cli::cli_alert_info("Performing spatial intersection...")
  tryCatch({
    results_polys <- results_polys %>%
      sf::st_intersection(extent) %>%
      sf::st_cast("POLYGON")
    cli::cli_alert_success("Spatial intersection complete.")
  }, error = function(e) {
    cli::cli_abort("Error during spatial intersection: {e$message}")
  })

  n_results <- nrow(results_polys)
  cli::cli_alert_success("Found {n_results} RESULTS polygons in the extent.")

  if (n_results == 0) {
    cli::cli_alert_warning("No results found. Returning an empty sf object.")
    return(sf::st_sf(geometry = sf::st_sfc()))
  }

  # Plotting
  if (plot) {
    cli::cli_alert_info("Preparing plot...")
    # Symbology Field Determination and Validation
    if (is.null(symbology_field)) {
      if ("FOREST_COVER_ID" %in% colnames(results_polys)) {
        symbology_field <- "FOREST_COVER_ID"
        cli::cli_alert_info("Using 'FOREST_COVER_ID' for symbology.")
      } else {
        cli::cli_abort("No symbology field provided and 'FOREST_COVER_ID' not found in data.")
      }
    }

    cli::cli_alert_info("Plotting with symbology field: {symbology_field}")
    tryCatch({
      # Convert symbology field to factor
      results_polys[[symbology_field]] <- as.factor(results_polys[[symbology_field]])

      p <- results_polys %>%
        terra::vect() %>%  # Use terra for plotting sf objects with ggplot2
        ggplot2::ggplot(data = .) +
        ggplot2::geom_spatvector(ggplot2::aes(fill = .data[[symbology_field]], color = "black", size = 0.1)) +
        ggplot2::coord_sf(expand = FALSE) +
        ggplot2::theme_bw() +
        ggplot2::labs(title = 'RESULTS Polygons for extent') +
        ggplot2::theme(legend.position = 'none')

      print(p)
      cli::cli_alert_success("Plotting complete.")
    }, error = function(e) {
      cli::cli_abort("Error during plotting: {e$message}")
    })

  }

  return(results_polys)
}
