#' DEM to Topographic Wetness Index
#'
#' @param dem_file location of elevation model raster
#' @param proj_dir directory where /output/raster/topography dir will be created and outputs stored
#'
#' @return saves a median filter, slope, filled dem, flow accumulation raster, and TWI raster in out_dir
#' @export
#'
#' @examples
#' \dontrun{
#' proj_dir <- 'H:/Quesnel_2022/blocks/CT1'
#' dem_file <- file.path(proj_dir, 'output/raster/dtm/ULS22_CT1_dtm_tin_0.25m.tif')
#' }
#'
calc_twi <- function(dem_file, proj_dir){

  out_dir <- paste0(proj_dir, '/output/raster/topography')

  if(!dir.exists(out_dir)){
    dir.create(out_dir, recursive = T)
  }

  # Median filter DEM

  whitebox::wbt_median_filter(dem_file,
                    file.path(out_dir, paste0(tools::file_path_sans_ext(basename(dem_file)), '_median_filt.tif')))

  # Generate slope raster from DEM

  whitebox::wbt_slope(
    file.path(out_dir, paste0(tools::file_path_sans_ext(basename(dem_file)), '_median_filt.tif')),
    file.path(out_dir, paste0(tools::file_path_sans_ext(basename(dem_file)), '_slope.tif')),
    zfactor=NULL,
    units="degrees")


  # Remove depressions (holes) in dem

  whitebox::wbt_breach_depressions_least_cost(
    file.path(out_dir, paste0(tools::file_path_sans_ext(basename(dem_file)), '_median_filt.tif')),
    file.path(out_dir, paste0(tools::file_path_sans_ext(basename(dem_file)), '_fill.tif')),
    dist=11,
    max_cost=NULL,
    min_dist=T,
    flat_increment=NULL,
    fill=T
  )

  # Generate flow accumulation raster

  whitebox::wbt_d8_flow_accumulation(
    file.path(out_dir, paste0(tools::file_path_sans_ext(basename(dem_file)), '_fill.tif')),
    file.path(out_dir, paste0(tools::file_path_sans_ext(basename(dem_file)), '_flow_acc.tif')),
    out_type="cells",
    log=F,
    clip=F,
    pntr=F,
    esri_pntr=F
  )

  # Calculate topographic wetness index

  whitebox::wbt_wetness_index(
    file.path(out_dir, paste0(tools::file_path_sans_ext(basename(dem_file)), '_flow_acc.tif')),
    file.path(out_dir, paste0(tools::file_path_sans_ext(basename(dem_file)), '_slope.tif')),
    file.path(out_dir, paste0(tools::file_path_sans_ext(basename(dem_file)), '_twi.tif'))
  )

}
