#' @title Convert LAS file to Alphashape and compute crown metrics
#'
#' @description The get_alphashape_metrics function takes a tree-attributed LAS object and converts it to a 3D alphashape.
#' Crown metrics such as X, Y, Zmax, Zq999, Zq99, Z_mean, n_points, vol_convex, vol_concave, vol_a05, CV_Z, and CRR are then computed for each tree.
#'
#' @param chunk an object of class LAS or a catalog chunk which has points attributed with treeID values, typically computed from lidR::segment_trees.
#' @param prog_bar logical. Indicates whether to display a progress bar for the computation (defaults to TRUE).
#'
#' @return A data frame with columns treeID, X, Y, Zmax, Zq999, Zq99, Z_mean, n_points, vol_convex, vol_concave, vol_a05, CV_Z, CRR.
#'
#' @examples
#' \dontrun{
#' library(alphashape3d)
#' library(lidR)
#' las <- readLAS(system.file("extdata", "uls.laz", package = "silvtools"))
#' tree_las <- segment_trees(las, mcwatershed(chm = rasterize_canopy(las, res = 1, p2r()),
#' treetops = locate_trees(las, lmf(ws = 5, hmin = 5))))
#' ashape_df <- get_alphashape_metrics(tree_las, prog_bar = TRUE)
#' }
#' @export
get_alphashape_metrics <- function(chunk, prog_bar = TRUE){

if ("LAS" %in% class(chunk)) {
    tree_las <- fchunkile
    print('Individual LAS object input into function')
} else{
tree_las <- lidR::readLAS(chunk)
}

if (is.empty(tree_las)) return(NULL)

if(prog_bar){
pb <- progress::progress_bar$new(
    format = "Generating crown metrics (:percent )| :elapsed elapsed | :eta eta",
    total = length(na.omit(unique((tree_las@data$treeID)))),
    width = 60)
} else{
  pb = NULL
}

print(glue::glue('Beginning crown metric generation for chunk'))

tree_las <- filter_duplicates(tree_las)

obs <- tree_las@data %>%
  dplyr::filter(!is.na(treeID)) %>%
  dplyr::select(X,Y,Z,treeID) %>%
  dplyr::group_by(treeID) %>%
  summarise(n = n()) %>% ungroup() %>% filter(n <= 4)

if(nrow(obs) > 0){
print(glue::glue('{nrow(obs)} treeIDs had 4 or fewer points and were discarded'))
}

mets <- tree_las@data %>%
  dplyr::filter(!is.na(treeID)) %>%
  dplyr::filter(!treeID %in% obs$treeID) %>%
  dplyr::select(X,Y,Z,treeID) %>%
  dplyr::group_by(treeID) %>%
  dplyr::summarise(ashape_metrics = get_crown_attributes(X, Y, Z, pb))

mets <- cbind(mets$treeID, mets$ashape_metrics)

colnames(mets)[1] <- "treeID"

return(mets)

}
