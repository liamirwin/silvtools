#' Calculate Tree Matching Scores
#'
#' @param match_table output from silvtools::tree_matching
#'
#' @return returns dataframe with metrics including number of FP/TP/FN as well as
#' precision, recall, Fscore, and statistics on distances between matches
#' @export
#'
#' @examples
#' \dontrun{
#' reference <- st_read('ground_reference_stems.shp')
#' detected <- locate_trees(chm, lmf(ws = 2))
#' match_table <- tree_matching(reference, detected)
#' scores <- tree_matching_scores(match_table)
#' }
tree_matching_scores <- function(match_table){

  # Quantify matching effectiveness

  status = match_table$status

  FP = sum(status == "FP")
  TP = sum(status == "TP")
  FN = sum(status == "FN")

  scores = data.frame(
    n_trees = TP + FN,
    n_detected = TP + FP,
    FP = FP,
    TP = TP,
    FN = FN,
    Precision = TP / (TP + FP),
    Recall = TP / (TP + FN),
    max_dist = max(na.omit(c(match_table$distance1, match_table$distance2))),
    min_dist = min(na.omit(c(match_table$distance1, match_table$distance2))),
    sd_dist = sd(na.omit(c(match_table$distance1, match_table$distance2))),
    mean_dist = mean(na.omit(c(match_table$distance1, match_table$distance2))))

  scores$Fscore = with(scores, 2 * Recall * Precision / (Recall + Precision))


  print(glue::glue('Out of {scores$n_trees} measured trees {scores$n_detected} were detected,
                   Matching resulted in {scores$TP} True Positives, {scores$FP} False Positives, and {scores$FN} False Negatives,
                   matching achieved a precision of {round(scores$Precision, 2)}, a recall of {round(scores$Recall, 2)} and an F-Score of {round(scores$Fscore, 2)}'))


  return(scores)

}
