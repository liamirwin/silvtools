#' Calculate KKL competition Indicies
#'
#' @param ttops tree tops object; expects an X.detected, Y.detected column, typically from silvtools::tree_matching() function
#' @param comp_input string value indicating competition index input specifically related to size; traditionally DBH/Height, crown volume etc
#'
#' @return returns original ttops dataset (sf/dataframe) with new cindex column
#' @export
#'
#' @examples
#' \dontrun{
#' c <- kkl_cindex(ttops, comp_input = 'Zmax')
#' }
kkl_cindex <- function(ttops, comp_input = 'Zmax'){

  tictoc::tic()

  # Convert ttops to ppp object
  trees_ppp <- ttops %>%
    dplyr::mutate(X = X, Y = Y) %>%
    as.data.frame() %>%
    dplyr::select(c(X, Y, comp_input))

  names(trees_ppp) <- c('X','Y','comp_value')

  trees_ppp <- trees_ppp %>%
    spatstat.geom::ppp(
      x = .$X,
      y = .$Y,
      window = spatstat.geom::owin(range(.$X),
                                   range(.$Y)),
      marks = .$comp_value)

  # Calculate KKL competition index for each tree (Z instead of dbh)
  kkl <- trees_ppp %>%
    di_steneker(dbh = .$comp_value, area = 1) # assuming that comp_value is equivalent to dbh and area is 1

  # Join new cindex with original ttops
  trees_cindex <- kkl %>%
    as.data.frame() %>% dplyr::select(x, y, aindex, afree) %>% sf::st_as_sf(coords = c('x', 'y'), crs = st_crs(ttops))

  ttops <- ttops %>% sf::st_as_sf(coords = c('X', 'Y'))

  trees_cindex <- trees_cindex %>%
    sf::st_join(., ttops)

  print(glue::glue('Calculated KKL style competition for {nrow(trees_cindex)} trees using size-dependent variable: {comp_input}'))

  tictoc::toc()

  return(trees_cindex)

}
