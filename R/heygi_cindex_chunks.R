#' Apply Buffer Function and Heygi Index Calculation to Each Chunk of Data
#'
#' This function chunks the input dataframe, applies the `create_buffer` function to each chunk,
#' then calculates the Heygi competition index for the buffered chunks.
#'
#' @param ttops tree tops sf object
#' @param comp_input string value indicating competition index input; traditionally DBH/Height, we test crown structural metrics
#' @param maxR radius of sphere of influence; treetops within this radius will be included in the calculation of cindex
#' @param quiet logical; if TRUE, suppresses progress messages
#' @param size Numeric. The size of each chunk (number of rows). Default is 10,000.
#' @param buffer_distance Numeric. The buffer distance to be applied around each point in the chunk. Default is 10.
#' @param comp_input Character. The name of the competition input variable. Default is 'vol_concave'.
#' @param maxR Numeric. The maximum radius for the Heygi competition index calculation. Default is 6.
#' @param num_cores Numeric. The number of cores to use for parallel processing. Default is 1.
#'
#' @return A dataframe with the Heygi competition index calculated for each chunk, and buffer points identified.
#'
#' @examples
#' \dontrun{
#' lasfile <- system.file("extdata", "Megaplot.laz", package = "lidR")
#' las <- readLAS(lasfile)
#' norm_las <- normalize_height(las, tin())
#' chm <- rasterize_canopy(norm_las, p2r(subcircle = 0.25), res = 1)
#' ttops <- locate_trees(chm, algorithm = lmf(ws = 1, hmin = 5))
#'
#' c <- calculate_heygi_for_chunks(ttops, size = 5000, buffer_distance = 5, comp_input = 'Z', maxR = 6, num_cores = 3)
#'
#' }
#'
#' @importFrom sf st_buffer st_union st_intersection
#' @importFrom purrr map
#' @importFrom furrr future_map
#' @importFrom future plan multisession
#' @importFrom dplyr filter
#' @export
calculate_heygi_for_chunks <- function(ttops, size = 10000, buffer_distance = 10, comp_input = 'vol_concave', maxR = 6, num_cores = 1) {

  tictoc::tic()

  # Chunk the data
  chunks <- chunk_data(ttops, size)

  if (num_cores == 1) {
    # Apply the create_buffer function to each chunk, then apply the heygi_cindex function
    results <- chunks %>%
      purrr::map( ~ create_buffer(., ttops, buffer_distance)) %>%
      purrr::map( ~ heygi_cindex(., comp_input = comp_input, maxR = maxR), .progress = T)

  } else {
    print(glue::glue("Processing {length(chunks)} chunks with {num_cores} cores"))
    future::plan(future::multisession, workers = num_cores)

    results <- chunks %>%
      furrr::future_map( ~ create_buffer(., ttops, buffer_distance)) %>%
      furrr::future_map( ~ heygi_cindex(., comp_input = comp_input, maxR = maxR), .progress = T)
  }

  # Combine the results into a single dataframe
  ttops_processed <- do.call(rbind, results)
  # Filter out the rows where buffer_point is FALSE
  ttops_filtered <- filter(ttops_processed, buffer_point == FALSE)

  tictoc::toc()
  # Return the processed dataframe
  return(ttops_filtered)
}
