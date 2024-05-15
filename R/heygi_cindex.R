#' Calculate Heygi Style Competition Indicies
#'
#' @param ttops tree tops object
#' @param comp_input string value indicating competition index input; traditionally DBH/Height, we test crown structural metrics
#' @param maxR radius of sphere of influence; treetops within this radius will be included in the calculation of cindex
#' @param quiet logical; if TRUE, suppresses progress messages
#'
#' @return returns original ttops dataset (sf/dataframe) with new cindex column
#' @export
#'
#' @examples
#' \dontrun{
#'
#' lasfile <- system.file("extdata", "MixedConifer.laz", package = "lidR")
#' las <- readLAS(lasfile)
#' norm_las <- normalize_height(las, tin())
#' chm <- rasterize_canopy(norm_las, p2r(subcircle = 0.25), res = 1)
#' ttops <- locate_trees(chm, algorithm = lmf(ws = 2, hmin = 5))
#'
#' c <- heygi_cindex(ttops, comp_input = 'Z', maxR = 6)
#' }
heygi_cindex <- function(ttops, comp_input = 'vol_convex', maxR = 6, quiet = FALSE){
  # Start timer if not quiet
  if(!quiet){
  tictoc::tic()
  }

  # If X column does not exist; extract from geometry
  if(!('X' %in% names(ttops))){
    ttops <- ttops %>%
      dplyr::mutate(X = unlist(purrr::map(.$geometry,1)),
                    Y = unlist(purrr::map(.$geometry,2)))
  }

  # Print warning if there are more than 30000 rows
  if(nrow(ttops) > 30000){
    warning('This function may take a long time to run with more than 30,000 tree tops - use heygi_cindex_chunks instead')
  }

  # Convert ttops to ppp object
  trees_ppp <- ttops %>%
    dplyr::mutate(X = X, Y = Y) %>%
    as.data.frame() %>%
    dplyr::select(X, Y, comp_input)

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
