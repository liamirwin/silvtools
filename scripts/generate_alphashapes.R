

blocks_dir <- list.dirs('H:/Quesnel_2022/process', recursive = FALSE)

is_dap <- FALSE

for(i in 1:length(blocks_dir))

tictoc::tic()

proj_dir <- blocks_dir[i]

if(stringr::str_detect(basename(proj_dir), pattern = 'DAP')){
  is_dap == TRUE
  # Set acquisition name (DAPYY_blockname)
  acq <- paste0('DAP22_', stringr::str_replace(basename(proj_dir), pattern = "-DAP", replacement = ""))
  print(paste0('Set acqusition type as DAP named ', acq))
} else{
  is_dap == FALSE
  # Set acquisition name (ULSYY_blockname)
  acq <- paste0('ULS22_',basename(proj_dir))
  print(paste0('Set acqusition type as lidar (ULS) named ', acq))
}

raster_output <- glue::glue('{proj_dir}/output/raster')
vector_output <- glue::glue('{proj_dir}/output/vector')
image_output <- glue::glue('{proj_dir}/output/png/ashape_snapshots')

tree_las_list <- list.files(glue::glue('{proj_dir}/input/tree_las'), pattern = '.laz', full.names = T)
ashape_mets <- list()

if(!dir.exists(glue::glue('{proj_dir}/output/crowns/ashapes'))){
  dir.create(glue::glue('{proj_dir}/output/crowns/ashapes'), recursive = T)
  print(glue::glue('Created alphashape save directory for {acq}'))
}

for(i in 1:length(tree_las_list)){
  tictoc::tic()
  las <- readLAS(tree_las_list[i], filter = "-thin_with_voxel 0.05")
  print(glue::glue('Loaded las {i} of {length(tree_las_list)}'))
  ashape_mets[[i]] <- get_alphashape_metrics(las)
  write.csv(ashape_mets[[i]], file = glue::glue('{proj_dir}/output/crowns/ashapes/{acq}_chunk_5cmvoxel_{i}_ashapes.csv'))
  print(glue::glue('Done processing for chunk {i}'))
  tictoc::toc()
}

ashape_df <- do.call(rbind, ashape_mets)
ashape_sf <- ashape_df %>% st_as_sf(coords = c('X','Y'), crs = 26910, remove = FALSE)

write.csv(ashape_df, glue::glue('{proj_dir}/output/crowns/ashapes/{acq}_chunk_5cmvoxel_all_ashapes.csv'))
st_write(ashape_sf, glue::glue('{proj_dir}/output/crowns/ashapes/{acq}_chunk_5cmvoxel_ashape_ttops.gpkg'))


# Create a square grid with 100m2 cells
square_grid <- st_make_grid(ashape_sf, cellsize = 10, square = TRUE)

# Create a hexagonal grid with 100m2 cells
hex_grid <- st_make_grid(st_as_sfc(st_bbox(ashape_sf)), cellsize = 10, square = FALSE)

bdy <- st_read(list.files(glue::glue('{proj_dir}/input/vector/'), pattern = '.gpkg', full.names = T))
a <- sf::st_intersection(square_grid, bdy)
b <- sf::st_intersection(a, ashape_sf)
b <- sf::st_intersects(a, ashape_sf)
b <- st_within(ashape_sf, a, sparse = FALSE)
c <- ashape_sf[,]
pointsID <- st_join(ashape_sf, a)
# Assign each point to a cell in the hexagonal grid
hex_cells <- tile.which(hex_grid, ashape_sf)
hex_means <- aggregate(ashape_sf$vol_concave, list(hex_cells), mean, na.rm = TRUE)

# Convert the square grid to a raster using terra
square_raster <- rast(square_grid, values = square_means$x)


c <- st_contains(a, ashape_sf) %>% as.data.frame() %>% rename(cell_id = row.id, point_id = col.id)

sf <- ashape_sf %>% mutate(row_num = row_number())

sf <- sf %>%
  left_join(c, by = c("row_num" = "point_id")) %>%
  mutate(cell_id = ifelse(is.na(cell_id), "NA", as.character(cell_id)))

x <- sf %>% filter(!is.na(cell_id)) %>%
  group_by(cell_id) %>%
  summarise(mean_cv = mean(vol_concave),
            sum_cv = sum(vol_concave))

