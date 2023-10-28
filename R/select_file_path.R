#' Select file path from directory
#'
#' This function checks a specified directory for files matching a specified pattern,
#' and either loads a file if there's only one match, or gives the user an option
#' to choose which file to load if there are multiple matches.
#' If no files match the specified pattern, an error is returned.
#'
#' @param file_dir A string specifying the directory to check.
#' @param pattern A string specifying the pattern to match file names against.
#' @param ext A string specifying the file extension to filter on. Default is '.tif$'.
#' @return A string containing the full path of the selected file.
#' @examples
#' \dontrun{
#' select_file_path("path/to/directory", "pattern")
#' }
#' @export
select_file_path  <- function(file_dir, pattern, ext = '.tif$') {

  assertthat::assert_that(!is.null(pattern))
  assertthat::assert_that(dir.exists(file_dir))
  assertthat::assert_that(is.character(file_dir), length(file_dir) == 1)
  assertthat::assert_that(is.character(pattern), length(pattern) == 1)
  assertthat::assert_that(is.character(ext), length(ext) == 1)

  files <- list.files(file_dir, pattern = pattern, full.names = TRUE)

  # Filter files that meet the raster subset criterion
  files <- stringr::str_subset(files, pattern = ext)

  if (length(files) > 1) {
    print(glue::glue("Multiple files matched the expected input"))

    for (k in 1:length(files)) {
      print(glue::glue("{k}. {basename(files[k])}"))
    }

    repeat {
      user_choice <- as.integer(readline("Enter the number corresponding to the correct file: "))
      if (user_choice >= 1 && user_choice <= length(files)) {
        cat(glue::glue("\nSelected file: {basename(files[user_choice])}\n"))
        break
      } else {
        cat("Invalid choice. Please enter a number between 1 and", length(files), ".\n")
      }
    }

    file <- files[user_choice]
  } else if (length(files) == 1) {
    file <- files[1]
  } else {
    stop(glue::glue("No files matched the expected input in directory {file_dir} with pattern {pattern}. Please check the directory and try again."))

  }

  return(file)
}
