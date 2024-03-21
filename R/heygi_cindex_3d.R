#' Calculate Heygi Style Competition Indicies
#'
#' @param ttops tree tops object; expects an X.detected, Y.detected column, typically from silvtools::tree_matching() function
#' @param comp_input string value indicating competition index input; traditionally DBH/Height, we test crown structural metrics
#' @param maxR radius of sphere of influence; treetops within this radius will be included in the calculation of cindex
#' @param inverted_search_cone logical; if TRUE, the function will calculate the competition index using the inverted search cone method
#' @param cone_height numeric; height of the search cone defined as percentage of tree top (apex) height
#' @param cone_angle numeric; angle of the search cone defined as degrees from the apex (default 60 degrees)
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

  if(!inverted_search_cone){
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

  if(inverted_search_cone){
    # Convert ttops to ppp object but keep tree height
    trees_ppp <- ttops %>%
      dplyr::mutate(X = X, Y = Y) %>%
      as.data.frame() %>%
      dplyr::select(X, Y, comp_input, Z)

    # Rename columns
    names(trees_ppp) <- c('X','Y','comp_value', 'Z')

    trees_ppp <- trees_ppp %>%
      spatstat.geom::ppp(
        x = .$X,
        y = .$Y,
        window = spatstat.geom::owin(range(.$X),
                                     range(.$Y)),
        marks = .$comp_value)




}
}



search_cone <- function(ttops, start_height_threshold = 0.6, base_height = 50, angle_degrees = 60){

  angle_radians <- angle_degrees * (pi/180)

  # Calculate the radius of the cone's base
  ttops <- ttops %>% mutate(
    cone_start_height = Z * start_height_threshold,
    cone_base_radius = cone_start_height / tan(angle_radians),
    cone_angle = angle_degrees,
    cone_base_height = base_height,
    # Cone base radius at height above tree
    cone_radius_at_height = cone_base_radius + (base_height / tan(angle_radians))
  )
  return(ttops)
}

#ttops_cone <- search_cone(trees_ppp, start_height_threshold = 0.6, angle_degrees = 60, base_height = 30)

# Figure out neighbouring trees within search radius

find_neighbours <- function(ttops, k = nrow(ttops)){

  # If not sf object, convert to sf
  if(!inherits(ttops, 'sf')){
    ttops <- ttops %>% sf::st_as_sf(coords = c('X', 'Y'))
  }

  nn <- nngeo::st_nn(ttops, ttops, k = nrow(ttops), returnDist = TRUE)

  # Collapse to list of dataframes with tree ID and distance
  df_list <- map2(nn[[1]], nn[[2]], ~data.frame(nn_ID = .x, nn_dist = .y))
  # Add tree ID index to each dataframe
  df_list <- imap(df_list, ~ mutate(.x, treeID = .y))
  # Extract the X, Y, Z information for joining
  coords_data <- ttops %>%
    mutate(ID = row_number()) %>%
    select(ID, X, Y, Z)
  # Attach x y and z of neighbours to dataframe
  df_list <- map(df_list, ~left_join(.x, coords_data, by = c('nn_ID' = 'ID')))

  return(df_n)
}

#ttops_n <- find_neighbours(ttops_cone)

# Function to annotate neighbors within search cone
annotate_within_cone <- function(ttops_cone, ttops_n) {
  # For each tree's neighbors, check if they are within the search cone
  ttops_cone_with_radius <- ttops_cone %>% select(treeID, cone_radius_at_height)

  ttops_n_annotated <- map2(ttops_n, ttops_cone_with_radius$treeID, function(neighbors, id) {
    cone_radius <- filter(ttops_cone_with_radius, treeID == id)$cone_radius_at_height
    mutate(neighbors, within_cone = nn_dist <= cone_radius)
  })

  return(ttops_n_annotated)
}

# Annotate neighbors within search cone
# ttops_n_annotated <- annotate_within_cone(ttops_cone, ttops_n)
#
# # Example output for the first tree
# ttops_n_annotated[[1]]



plot_neighbours <- function(df) {
  # Find the target tree (the one with nn_dist = 0.0)
  target_tree <- df[1,]
  # Calculate relative X and Y coordinates for neighbouring trees
  df <- df %>% mutate(X = X - target_tree$X, Y = Y - target_tree$Y)
  # Set Target tree X and Y to 0
  target_tree$X <- 0
  target_tree$Y <- 0
  # Create the plot
  p <- ggplot() +
    geom_point(data = df, aes(x = X, y = Y, color = nn_dist), size = 3) +
    geom_point(data = target_tree, aes(x = X, y = Y), color = 'red', size = 3) +
    scale_color_gradient(low = "blue", high = "red") +
    coord_fixed() +
    theme_minimal() +
    labs(title = "Neighbouring Trees", x = "X coordinate", y = "Y coordinate", color = "Distance") +
    geom_text(data = target_tree, aes(x = X, y = Y, label = nn_ID), vjust = -1)

  # Print the plot
  print(p)
}

# Call the function with the first dataframe in your list as an example
#plot_neighbours(df_list[[1]])
