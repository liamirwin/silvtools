#' Run PyInstaller-packaged Python Script to Calculate Crown Volumes (convex)
#'
#' This function runs a PyInstaller-packaged Python script to calculate the convex volume of 3D points from a Parquet file.
#'
#' @param input_parquet Path to the input Parquet file containing 3D points.
#' @param output_parquet Path to the output Parquet file where the results will be saved.
#' @param py_executable Path to the PyInstaller-packaged Python executable.
#' @return A character vector containing the output from the Python script execution.
#' @examples
#' \dontrun{
#' input_path <- "G:/calc_concave_vol/tree_las.parquet"
#' output_path <- "G:/calc_concave_vol/tree_las_output_2.parquet"
#' executable_path <- "G:/calc_concave_vol/dist/calc_concave_volume/calc_concave_volume.exe"
#' result <- calc_convex_vol_py(input_path, output_path, executable_path)
#' print(result)
#' }
#' @export
calc_convex_vol_py <- function(input_parquet, output_parquet, py_executable){

  command <- paste(py_executable, input_parquet, output_parquet)
  result <- system(command, intern = TRUE)

  return(result)

}

