#' Match ground truth with detected tree tops
#' Calculates distance between ground truth and nearest 2 detected tree tops, returns a matched dataset containing the
#' tree ID of the ground truth and detected matches
#' @param reference sf object - points of recorded tree locations gathered in field
#' @param detected sf object - tree top locations typically from lidR::locate_trees
#' @param PlotID character string; ID of plot where accuracy is being assessed
#'
#' @return returns a tree matching data table
#' @export
#'
#' @examples
#' \dontrun{
#' detected <- lidR::locate_trees(las, algorithm = lmf(ws = 2, hmin = 10))
#'
#' reference <- sf::st_read('F:/Quesnel_2022/Quesnel_2022_PosTex/
#' georeferenced_adjusted/CT1P1_postex_utm10n.shp')
#'
#' matched_trees <- tree_matching(reference, detected)
#' }
tree_matching = function(reference, detected, PlotID)
{
  stopifnot(is(detected, "sf"))
  stopifnot(is(reference, "sf"))

  reference <- reference %>%
    dplyr::rename(X_postex = X, Y_postex = Y) %>%
    dplyr::mutate(X = unlist(purrr::map(.$geometry,1)), Y = unlist(purrr::map(.$geometry,2)),
           PLOTID = PlotID) %>% sf::st_drop_geometry()

  detected <- detected %>%
    dplyr::mutate(X = unlist(purrr::map(.$geometry,1)), Y = unlist(purrr::map(.$geometry,2))) %>% sf::st_drop_geometry()

  xy_truth    = reference %>% dplyr::select(X, Y)
  xy_detected = detected %>% dplyr::select(X, Y)
  x_truth     = xy_truth[,1]
  y_truth     = xy_truth[,2]
  x_detected  = xy_detected[,1]
  y_detected  = xy_detected[,2]
  z_detected  = detected$Z

  # Attribution of nearest and 2nd neareast referenced tree index for each detected tree
  tree <- SearchTrees::createTree(xy_truth)
  # Finds 2 nearest neighbours using xy values of truth/detected trees
  knn  <- SearchTrees::knnLookup(tree, newdat = xy_detected, k = 2L)
  # 1st nearest neighbour ID
  inds1 <- knn[,1]
  # 2nd nearest neighbour ID
  inds2 <- knn[,2]

  detected$PLOTID <- reference$PLOTID[inds1]

  match_table <- data.table::data.table(index_detected = 1:length(inds1),
                                        index_truth1   = inds1,
                                        index_truth2   = inds2)

  # Horizontal distance between detected tree and two truth neighbours (m)
  match_table$distance1 <- sqrt((x_truth[inds1] - x_detected)^2 + (y_truth[inds1] - y_detected)^2)
  match_table$distance2 <- sqrt((x_truth[inds2] - x_detected)^2 + (y_truth[inds2] - y_detected)^2)

  # Takes 10% of tree height as maximum distance, if this is less than 2m make 2m
  dist_max = z_detected*0.10
  dist_max[dist_max < 5] = 5
  # := operator assigns the column index_truth a value of NA if distance is dist max
  match_table[distance1 > dist_max, index_truth1 := NA]
  match_table[distance2 > dist_max, index_truth2 := NA]
  # Get IDs of closest trees
  id = match_table[, .I[which.min(distance1)], by = index_truth1]
  match_table$index_truth1 = NA_integer_
  match_table[id$V1, index_truth1 := id$index_truth1]
  # Get IDs of second closest trees
  id = match_table[, .I[which.min(distance2)], by = index_truth2]
  match_table$index_truth2 = NA_integer_
  match_table[id$V1, index_truth2 := id$index_truth2]
  # If the same index value appears in both columns take the one where it is a shorter distance (index 1)
  match_table[index_truth2 %in% index_truth1, index_truth2 := NA_integer_]
  # Collate the two index truth ID columns into one
  match_table$index_truth = ifelse(is.na(match_table$index_truth1), match_table$index_truth2, match_table$index_truth1)

  ###

  # Create a rowID column for the ground truth trees
  reference$num_tree = 1:nrow(reference)
  # Apply that ID to the row of the detected tree matched to it
  detected$num_tree = match_table$index_truth
  detected$distance1 = match_table$distance1
  detected$distance2 = match_table$distance2
  # Convert to data table format
  dt_reference = data.table::as.data.table(reference)
  dt_detected     = data.table::as.data.table(detected)

  X = dplyr::full_join(dt_reference, dt_detected, by = "num_tree")
  X$PLOTID <- ifelse(is.na(X$PLOTID.x), X$PLOTID.y, X$PLOTID.x)
  X = dplyr::select(X, -PLOTID.x, -PLOTID.y)

  data.table::setDT(X)

  if ("Z" %in% names(dt_reference)){
    data.table::setnames(X, c("X.x", "Y.x", "Z.x", "X.y", "Y.y", "Z.y"), c("X", "Y", "Z", "X.detected", "Y.detected", "Z.detected"))
  }else {
    data.table::setnames(X, c("X.x", "Y.x", "X.y", "Y.y", "Z"), c("X", "Y", "X.detected", "Y.detected", "Z.detected"))
  }

  X$status = "FN"
  X$status[is.na(X$num_tree)] = "FP"
  X$status[!is.na(X$Z.detected) & !is.na(X$num_tree)] = "TP"
  class(X) <- append("TreeMatching", class(X))
  return(X)
}
