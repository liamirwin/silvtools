#' Interactive Transect Clipping with lidR
#'
#' This function interactively clips a transect from a LAS dataset. It first
#' generates a surface height model and displays it. The user is then prompted
#' to click two points (start and end of the transect) on the plot. A transect polygon
#' is constructed from these points and used to clip the LAS dataset.
#'
#' @param las A character string specifying the path to a LAS/LAZ file or a LAS object.
#' @param width Numeric. The width of the transect. Default is 1.
#' @param res Numeric. The resolution for rasterizing the canopy height model. Default is 0.5.
#' @param subcircle Numeric. Parameter for the \code{p2r} function used in canopy rasterization.
#'        Default is 0.
#'
#' @return A LAS object corresponding to the clipped transect.
#'
#' @import lidR
#' @import sf
#' @import cli
#' @importFrom glue glue
#'
#' @examples
#' \dontrun{
#'   # Using the Megaplot example data from lidR:
#'   lasFile <- system.file("extdata", "Megaplot.laz", package = "lidR")
#'
#'   las_clipped <- clip_transect_interactive(
#'       las       = lasFile,
#'       width     = 5,
#'       res       = 1,
#'       subcircle = 0.25
#'     )
#'
#' }
#'
#' @export
clip_transect_interactive <- function(las, width = 5, res = 1, subcircle = 0) {

  # Check if 'las' is a file path or a LAS object, and load the file if needed
  if (is.character(las)) {
    las <- tryCatch(
      lidR::readLAS(las),
      error = function(e) {
        cli::cli_progress_message(NULL)
        cli::cli_alert_danger(glue::glue("Error reading LAS file: {e$message}"))
        stop(e)
      }
    )

    if (is.null(las)) {
      cli::cli_progress_message(NULL)
      cli::cli_alert_danger("LAS file could not be read or is empty.")
      stop("LAS file could not be read or is empty.")
    }
    cli::cli_progress_step("Load LAS file")
  } else if (!inherits(las, "LAS")) {
    cli::cli_alert_danger("Input is neither a LAS object nor a valid file path.")
    stop("Input is neither a LAS object nor a valid file path.")
  }

  # Generate the canopy height model (CHM) at the specified resolution
  cli::cli_progress_message(glue::glue("Rasterizing canopy height model at {res}m spatial resolution..."))
  chm <- tryCatch(
    lidR::rasterize_canopy(las, res = res, lidR::p2r(subcircle = subcircle)),
    error = function(e) {
      cli::cli_progress_message(NULL)
      cli::cli_alert_danger(glue::glue("Error during CHM generation: {e$message}"))
      stop(e)
    }
  )
  cli::cli_progress_step("Generate surface model")

  # Close extra graphics devices if any are open
  while (length(dev.list()) > 1) {
    dev.off()
  }

  # Plot the CHM and prompt user for two clicks (the start and end points of the transect)
  plot(chm, main = "Click on two points to define the transect")
  cli::cli_progress_message("Please click on two points in the plot (start and end of transect)...")
  clicks <- locator(2)
  cli::cli_progress_step("User selected points")

  # Validate that two points were selected
  if (length(clicks$x) < 2 || length(clicks$y) < 2) {
    cli::cli_alert_danger("You must click on two points to define the transect.")
    stop("You must click on two points to define the transect.")
  }

  # Compute the transect endpoints and check that they are distinct
  p1 <- c(clicks$x[1], clicks$y[1])
  p2 <- c(clicks$x[2], clicks$y[2])
  dx <- p2[1] - p1[1]
  dy <- p2[2] - p1[2]
  line_length <- sqrt(dx^2 + dy^2)

  if (line_length == 0) {
    cli::cli_alert_danger("The two clicked points are identical. Please select two distinct points.")
    stop("The two clicked points are identical. Please select two distinct points.")
  }

  # Compute perpendicular vectors to define the transect's width
  half_width <- width / 2
  ux <- -dy / line_length
  uy <- dx / line_length

  p1_left  <- c(p1[1] + half_width * ux, p1[2] + half_width * uy)
  p2_left  <- c(p2[1] + half_width * ux, p2[2] + half_width * uy)
  p2_right <- c(p2[1] - half_width * ux, p2[2] - half_width * uy)
  p1_right <- c(p1[1] - half_width * ux, p1[2] - half_width * uy)

  polygon_coords <- rbind(p1_left, p2_left, p2_right, p1_right, p1_left)

  # Create an sf polygon from the clicked points
  transect_poly <- tryCatch(
    {
      poly <- sf::st_polygon(list(polygon_coords))
      # Retrieve CRS from the CHM or fallback to LAS
      crs_info <- tryCatch(
        sf::st_crs(chm),
        error = function(e) sf::st_crs(las)
      )
      sf::st_sfc(poly, crs = crs_info)
    },
    error = function(e) {
      cli::cli_progress_message(NULL)
      cli::cli_alert_danger(glue::glue("Error creating transect polygon: {e$message}"))
      stop(e)
    }
  )
  cli::cli_progress_step("Polygon created")

  # Overlay the transect polygon on the CHM plot
  plot(transect_poly, add = TRUE, border = "red", lwd = 2)

  # Clip the LAS data using the transect polygon
  las_clipped <- tryCatch(
    lidR::clip_roi(las, transect_poly),
    error = function(e) {
      cli::cli_progress_message(NULL)
      cli::cli_alert_danger(glue::glue("Error clipping the LAS object: {e$message}"))
      stop(e)
    }
  )
  cli::cli_progress_step("LAS object clipped")

  # Final message showing the number of points in the clipped LAS
  message(glue::glue(
    "Finished clipping {width}m transect from the LAS object (",
    "{nrow(las_clipped@data)} points)."
  ))

  return(las_clipped)
}
