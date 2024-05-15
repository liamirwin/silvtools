#' Create Buffer Around Points in an sf DataFrame
#'
#' This function creates a buffer around each chunk of points in a given sf dataframe,
#' dissolves borders between the buffered points, and updates the dataframe to indicate
#' which points fall within the buffer but are not in the original chunk.
#'
#' @param chunk An sf dataframe containing the chunk of points to be buffered.
#' @param df An sf dataframe containing the original set of points.
#' @param distance Numeric. The buffer distance to be applied around each point in the chunk.
#'
#' @return An sf dataframe with an additional column `buffer_point` indicating which points
#' fall within the buffer and are not in the original chunk.
#'
#' @examples
#' \dontrun{
#'   library(sf)
#'   # Create example data
#'   df <- st_as_sf(data.frame(treeID = 1:10, geometry = st_sfc(st_point(c(0, 0)), st_point(c(1, 1)),
#'   st_point(c(2, 2)), st_point(c(3, 3)), st_point(c(4, 4)), st_point(c(5, 5)), st_point(c(6, 6)),
#'   st_point(c(7, 7)), st_point(c(8, 8)), st_point(c(9, 9)))))
#'   chunk <- df[1:5, ]
#'   # Apply buffer
#'   result <- create_buffer(chunk, df, 2)
#'   print(result)
#' }
#'
#' @importFrom sf st_buffer st_union st_intersection
#' @export
create_buffer <- function(chunk, df, distance){
  # Create a buffer around each point in the chunk
  buffer <- st_buffer(chunk, distance)
  # Dissolve borders between points in the buffer
  buffer <- st_union(buffer)
  # Add a new column 'buffer_point' to the dataframe and initialize it as FALSE
  df$buffer_point <- FALSE
  # Identify which points from the original data fall within the buffer
  intersects_buffer <- st_intersection(df, buffer)
  # Update 'buffer_point' to TRUE for points that are in the buffer but not in the original chunk
  intersects_buffer$buffer_point <- !intersects_buffer$treeID %in% chunk$treeID
  # Return the updated dataframe
  return(intersects_buffer)
}

#' Chunk Dataframe into Smaller Chunks
#'
#' This function splits a dataframe into smaller chunks of a specified size.
#'
#' @param df A dataframe to be chunked.
#' @param size Numeric. The size of each chunk (number of rows). Default is 10,000.
#'
#' @return A list of dataframes, each representing a chunk of the original dataframe.
#'
#' @examples
#' \dontrun{
#'   # Create example data
#'   df <- data.frame(x = 1:100000, y = rnorm(100000))
#'   # Split data into chunks of 10,000 rows each
#'   chunks <- chunk_data(df, size = 10000)
#'   print(length(chunks))  # Should print 10
#' }
#'
#' @export
chunk_data <- function(df, size = 10000){
  # Calculate the number of rows and the number of chunks needed
  n <- nrow(df)
  chunks <- ceiling(n / size)
  # Divide the row indices into chunks
  indices <- cut(1:n, breaks = chunks, labels = FALSE)
  # Split the data into chunks using the indices
  split(df, indices)
}
