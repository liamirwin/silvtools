#' Calculate Heygi Style Competition Indicies
#'
#' @param ttops tree tops object; expects an X.detected, Y.detected column, typically from silvtools::tree_matching() function
#' @param comp_input string value indicating competition index input; traditionally DBH/Height, we test crown structural metrics
#' @param maxR radius of sphere of influence; treetops within this radius will be included in the calculation of cindex
#' @param quiet logical; if TRUE, suppresses progress messages
#'
#' @return returns original ttops dataset (sf/dataframe) with new cindex column
#' @export
#'
#' @examples
#' \dontrun{
#' c <- heygi_cindex(ttops, comp_input = 'vol_convex', maxR = 6)
#' }
heygi_cindex <- function(ttops, comp_input = 'vol_convex', maxR = 6, quiet = FALSE){
  # Start timer if not quiet
  if(!quiet){
  tictoc::tic()
  }
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

  # Calculate heygi competition index for each tree (Z instead of dbh)
  heygi <- trees_ppp %>%
    siplab::pairwise(., maxR=maxR, kernel=siplab::powers_ker,
                     kerpar=list(pi=1, pj=1, pr=1, smark=1))


  # Join new cindex with original ttops

  trees_cindex <- heygi %>%
    as.data.frame() %>% dplyr::select(x, y, cindex)

  ttops <- ttops %>% sf::st_as_sf(coords = c('X', 'Y'))

  trees_cindex <- dplyr::bind_cols(ttops, trees_cindex)

  # Print progress if not quiet
  if(!quiet){
  print(glue::glue('Calculated Heygi style competition for {nrow(trees_cindex)} trees assesing their {comp_input} within a {maxR}m radius'))
  tictoc::toc()
  }

  return(trees_cindex)
}
