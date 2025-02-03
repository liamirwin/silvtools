#' Run the Alphacrown Shiny App
#'
#' This function starts the Alphacrown Shiny App, which allows users to analyze tree crowns using the specified folder path and alpha value.
#' Simply input a folder containing individual tree las files and the app will allow you to visualize different trees with different alpha values and provide a set of summary metrics.
#'
#' @return No return value. This function runs a Shiny app.
#' @examples
#' \dontrun{
#'
#' library(silvtools)
#' # Get the example dataset path
#' example_data_path <- system.file("extdata", "individual_trees", package = "silvtools")
#'
#' # Launch the Shiny app with the example dataset
#' run_alphacrown(example_data_path)
#' }
#'
#' @import shiny rgl glue lidR alphashape3d
#' @export
#'
run_alphacrown <- function(){


ui <- fluidPage(
  titlePanel("Alphacrown"),
  tags$style(type = 'text/css', '
  #metrics_table {width: 100%; overflow-x: auto;}
  #tree_plot {border: 1px solid #ccc; padding: 5px;}
  .btn {margin-right: 4px;}
  .btn-primary {background-color: #3399ff; border-color: #3399ff;}
  .btn-secondary {background-color: #C4C7CC; border-color: #cccccc;} # Lighten up the gray color
  .btn-info {background-color: #17a2b8; border-color: #17a2b8;}
'),
  fluidRow(
    column(3,
           textInput("folder", "Enter Tree Folder Path"),
           numericInput("alpha", "Alpha Value:", 1),
           actionButton("run_analysis", "Run Analysis", class = "btn btn-primary"),
           uiOutput("tree_slider_ui"),
           actionButton("prev_tree", "Previous Tree", class = "btn btn-secondary"),
           actionButton("next_tree", "Next Tree", class = "btn btn-secondary"),
           checkboxInput("decimate", "Decimate LAS Data", value = TRUE),
           checkboxInput("alpha_inf", "Infinite Alpha Value", value = FALSE),
           textOutput("timer_log"), # Add the timer log output here
    ),
    column(9,
           rglwidgetOutput("tree_plot"),
           div(style = "overflow-x: auto;", tableOutput("metrics_table"))
    )
  )
)


# Server
server <- function(input, output, session) {
  # Reactive value for the folder path
  folder_path <- reactive({
    req(input$folder)
    # Replace backslashes with forward slashes for R compatibility
    gsub("\\\\", "/", input$folder)
  })

  # Reactive value for the tree index
  tree_index <- reactive(input$tree_slider)

  # Function to analyze the ashape
  analyze_ashape <- function(las_file, alpha_value, decimate = FALSE) {

    snapshot_dir <- glue::glue('{folder_path()}/snapshots')
    if (!dir.exists(snapshot_dir)) {
      dir.create(snapshot_dir)
    }
    las <- lidR::readLAS(las_file)

    if(decimate == TRUE){
      n1 <- length(las@data$Z)
      las <- lidR::decimate_points(las, random_per_voxel(res = 0.25, n = 1))
      n2 <- length(las@data$Z)
    }

    X <- las@data$X
    Y <- las@data$Y
    Z <- las@data$Z
    ID <- unique(las@data$treeID)

    if (length(X) <= 3 || length(Y) <= 3 || length(Z) <= 3) {
      print('Cannot compute a 3D hull from 3 or fewer points')
      return(NULL)
    }

    # alphashadep3d
    a3d = cbind(X, Y, Z)
    a3d <- unique(a3d)
    # get treetop location
    top_x <- a3d[which.max(a3d[,3]),1]
    top_y <- a3d[which.max(a3d[,3]),2]
    # normalize X and Y
    a3d[,1] = a3d[,1] - mean(a3d[,1]) #center points around 0,0,0
    a3d[,2] = a3d[,2] - mean(a3d[,2]) #center points around 0,0,0

    # Replace this line:
    # alpha <- c(Inf, 1)
    # with:
    alpha <- c(Inf, alpha_value)

    ashape <- alphashape3d::ashape3d(x = a3d, alpha = alpha)

    # calculate crown metrics
    df <- data.frame(
      # Crown height
      Zmax = max(Z),
      Zq999 = as.numeric(quantile(Z, 0.999)),
      Zq99 = as.numeric(quantile(Z, 0.990)),
      Z_mean = mean(Z),
      # Crown size
      n_points = length(Z),
      vol_concave = alphashape3d::volume_ashape3d(ashape, indexAlpha = 1),
      vol_convex = alphashape3d::volume_ashape3d(ashape, indexAlpha = 2),
      # Crown complexity
      CV_Z = sd(Z) / mean(Z),
      CRR = (mean(Z) - min(Z)) / (max(Z) - min(Z)),
      # Append georeferenced tree top coordinates
      X = top_x,
      Y = top_y)

    open3d() # Open a new RGL device

    # Plot Concave
    plot(ashape, indexAlpha = 2, transparency = 0.4, axes = FALSE)
    # Add Axes
    axes3d(
      edges=c('x--', 'y+-', 'z--'),
      labels=T
    )
    # Add Points
    points3d(a3d, color = 'black')
    in3d <- inashape3d(ashape, points = a3d)
    # Rotate Plot
    delta_angle <- 2
    angle <- rep(delta_angle * pi / 8, 8 / delta_angle)[1]
    i = 1
    while(i < 8) {
      view3d(userMatrix = rotate3d(par3d("userMatrix"), angle, 0, 0, 1))
      i = i + 1
    }
    # Add Title
    title3d(main = glue::glue('ID: {ID} - v_conc: {round(df$vol_concave)} v_conv: {round(df$vol_convex)} n: {df$n_points} CRR: {round(df$CRR, 3)}'))
    par3d(windowRect = c(0, 0, 600, 600))


    rglwidget <- suppressWarnings(rglwidget(controllers = "none"))
    close3d()

    return(list(plot = rglwidget, metrics = df, snapshot_dir = snapshot_dir))

  }

  # Reactive value for the list of tree files
  tree_list <- reactive({
    # Check if folder is valid
    validate(
      need(file.exists(folder_path()), "The folder path does not exist or is invalid. Please provide a valid path.")
    )
    # List laz/las  files
    list.files(folder_path(), pattern = '\\.(laz|las)$', full.names = T, ignore.case = TRUE)
  })

  output$tree_slider_ui <- renderUI({
    req(length(tree_list()))
    sliderInput("tree_slider", "Select Tree ID:", min = 1, max = length(tree_list()), value = 1)
  })

  # Render rgl plot in Shiny

  output$tree_plot <- renderRglwidget({
    req(plot_data())

    try(close3d(), silent = TRUE)

    plot_data()$plot
  })

  output$metrics_table <- renderTable({
    req(plot_data())
    plot_data()$metrics
  })

  observeEvent(input$prev_tree, {
    if (tree_index() > 1) {
      updateSliderInput(session, "tree_slider", value = tree_index() - 1)
    }
  })

  observeEvent(input$next_tree, {
    if (tree_index() < length(tree_list())) {
      updateSliderInput(session, "tree_slider", value = tree_index() + 1)
    }
  })

  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("plot-", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      rgl.snapshot(file, fmt = "png")
    },
    contentType = "image/png"
  )

  start_time <- reactiveVal()
  end_time <- reactiveVal()

  # Update the plot and metrics when the tree_slider value changes

  observe({
    req(input$tree_slider)
    tree_las <- tree_list()[tree_index()]

    alpha_value <- if (input$alpha_inf) Inf else input$alpha # Conditionally set the alpha value

    validate(
      need(!is.na(alpha_value), "Please provide a valid alpha value.")
    )

    if (is.na(alpha_value)) {
      showNotification("Warning: Alpha value is missing.", type = "warning", duration = 5, closeButton = TRUE)
      return()
    }

    start_time(Sys.time())

    withProgress(message = 'Processing...', value = 0, {
      result <- analyze_ashape(tree_las, alpha_value, input$decimate)
      incProgress(0.5) # Increment progress to 50%

      Sys.sleep(1) # Sleep for 1 second to demonstrate progress bar
      incProgress(0.5) # Increment progress to 100%
    })

    end_time(Sys.time())

    if (!is.null(result)) {
      plot_data(result)
    } else {
      showNotification("Error: Could not compute the 3D hull.", type = "error", duration = 5, closeButton = TRUE)
    }
  })



  plot_data <- reactiveVal(NULL)

  output$timer_log <- renderText({
    req(end_time(), start_time())
    elapsed_time <- round(difftime(end_time(), start_time(), units = "secs"), 2)
    paste("Processing time:", elapsed_time, "seconds")
  })
}





# Run the app
shinyApp(ui = ui, server = server)


}
