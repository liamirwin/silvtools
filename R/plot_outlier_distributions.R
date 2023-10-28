#' Create a histogram with percentile vertical lines
#'
#' This function creates a histogram of a specified metric from a given dataframe and adds vertical lines at the specified percentiles. It supports optional faceting by a categorical variable.
#'
#' @param data A dataframe containing the metric of interest.
#' @param metric A character string specifying the column name of the metric to plot.
#' @param percentiles A numeric vector of percentiles to display as vertical lines (default: c(1, 2, 90, 95, 99)).
#'
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' # Create example data
#' example_data <- data.frame(id = 1:100,
#'                            height = rnorm(100, mean = 10, sd = 2),
#'                            species = factor(rep(c("A", "B", "C", "D"),
#'                             each = 25)))
#' # Plot the histogram
#' plot <- plot_outlier_distribution(example_data, metric = "height",
#' percentiles = c(1, 2, 90, 95, 99))
#' print(plot)
#' }
#' @export
plot_outlier_distribution <- function(data, metric, percentiles = c(1, 2, 90, 95, 99)) {
  # Determine cutoffs for given percentiles
  metric_data <- data %>% pull({{metric}})
  cutoff_values <- quantile(metric_data, probs = percentiles / 100)

  # Create a dataframe to store the x positions and labels for the vertical lines
  annot_df <- data.frame(x = cutoff_values, quantile = names(cutoff_values))

  # Create the histogram plot
  plot <- ggplot(data, aes(x = !!sym(metric))) +
    geom_histogram(bins = 30, color = "black", fill = "dodgerblue", alpha = 0.7) + # Histogram
    labs(title = paste0("Histogram of ", metric), x = metric, y = "Frequency") + # Axis labels and title

    # Add vertical lines at specified percentiles
    geom_vline(data = annot_df,
               mapping = aes(xintercept=x, color=quantile), linetype = "dashed", size = 0.8) +

    # Set the color scale for the vertical lines and order the legend items
    scale_color_manual(values = viridis::viridis(length(percentiles)),
                       breaks = names(sort(cutoff_values))) +

    # Apply a minimal theme to the plot
    theme_minimal()

  names(annot_df) <- c(metric, 'quantile')
  print(annot_df)

  return(plot)
}

