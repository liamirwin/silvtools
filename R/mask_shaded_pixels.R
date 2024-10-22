library(terra)
library(ggplot2)

# CHM should be pre-loaded by user as a SpatRaster to speed up processing

chm_rast <- rast('G:/Ontario_2023/Block_18/LAS/DLS23/output/raster/chm/B18_DLS23_chm_smooth_p2r_sc0.025_0.1m.tif')

mask_shaded_pixels <- function(input_dir, output_dir, nir_band_num = 10,
                               output_ext, chm_rast = NULL, chm_threshold = 0.5,
                               plot_histogram = FALSE) {

  # Load MX Raster
  mx_raster <- rast(input_dir)

  # Subset only the NIR band (Band 10 aka R842)
  nir_band <- mx_raster[[nir_band_num]]

  # Initial Mask out pixels based on coincident CHM < certain height (optional)
  if (!is.null(chm_rast)) {
    chm_crop <- terra::crop(chm_rast, nir_band)
    mask_rast <- chm_crop >= chm_threshold
    aligned_mask <- terra::extend(mask_rast, nir_band)
    aligned_mask <- terra::resample(mask_rast, nir_band)
    # If CRS doesnt match, warp aligned mask to match the NIR band
    if (crs(aligned_mask) != crs(nir_band)) {
      aligned_mask <- project(aligned_mask, nir_band)
    }
    nir_band <- mask(nir_band, mask = aligned_mask, maskvalues = 0)


  }

  mx_raster <- mx_raster/32752

  max_vals <- global(mx_raster, 'max')

  # Create an empty SpatRaster object to store the normalized layers
  normalized_raster <- mx_raster

  # Normalize each layer by dividing by its max value
  for (i in 1:nlyr(mx_raster)) {
    normalized_raster[[i]] <- mx_raster[[i]] / max_vals[i,]
  }

  mx_raster <- normalized_raster

  # Sum all band values for each pixel
  pan <- mx_raster[1] + mx_raster[2] +
    mx_raster[3] + mx_raster[4] + mx_raster[5] +
    mx_raster[6] + mx_raster[7] + mx_raster[8] + mx_raster[9] + mx_raster[10]

  pan <- app(mx_raster, fun = mean)
  pan <- pan/10


  writeRaster(pan, 'G:/Ontario_2023/Block_18/scratch/pan_3.tif')
  pan

  # Normalize each band by its maximum value


  nir_band <- nir_band/32752
  nir_mask <- nir_band > 0.33
  writeRaster(nir_band, 'G:/Ontario_2023/Block_18/scratch/mask33.tif')
  # Calculate the local minima between modes of the distribution
  nir_values <- values(nir_band, na.rm = TRUE)
  nir_values <- nir_values/32752
  density_nir <- density(nir_values)
  local_minima <- density_nir$x[which.min(density_nir$y)]


  # Generate histogram plot (optional by parameter)
  if (plot_histogram) {
    df <- data.frame(NIR = nir_values)
    names(df) <- "NIR"

   ggplot(df, aes(x = NIR)) +
      geom_density(fill = 'blue', color = 'black') +
      geom_vline(xintercept = local_minima, linetype = "dotted", color = "red") +
      annotate("text", x = local_minima, y = Inf, label = round(local_minima, 2), vjust = 2, color = "red") +
      ggtitle("NIR Band Reflectance Distribution") +
      xlab("Reflectance") +
      ylab("Frequency")

    # Save histogram image to subdirectory of the output directory (PLOT) as .png with ggsave
    plot_dir <- file.path(output_dir, "PLOT")
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE)
    }
    ggsave(filename = file.path(plot_dir, "nir_histogram.png"), plot = hist_plot, width = 8, height = 6)
  }

  # Mask pixels below the local minima threshold
  mask_threshold <- nir_band > local_minima
  masked_nir_band <- mask(nir_band, mask_threshold)

  # Save the masked NIR band raster
  output_file <- file.path(output_dir, paste0("masked_nir_band", output_ext, ".tif"))
  writeRaster(masked_nir_band, output_file, format = "GTiff", overwrite = TRUE)
}

# Example usage:
# mask_shaded_pixels(input_dir = "path/to/mx_raster.tif", output_dir = "path/to/output", output_ext = "_masked", chm_path = "path/to/chm_raster.tif", chm_threshold = 0.5, plot_histogram = TRUE)
