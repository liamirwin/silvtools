#' Configure Lidar Processing Switches
#'
#' This Shiny application allows users to set processing switches for lidar data processing.
#' It provides an interactive user interface to select various options and displays the selected options in a table.
#'
#' @import shiny
#' @import shinythemes
#' @import crayon
#'
#' @export
#'
#' @title Lidar Processing Switches
#'
#' @description A Shiny application for setting lidar processing switches.
#'
#' @seealso \code{\link[shiny]{shinyApp}}
#'
#' @examples
#' \dontrun{
#'   # Run the Shiny app
#'   configure_las_process()
#' }
configure_las_process <- function() {

  # Define UI for application
  ui <- fluidPage(
    theme = shinytheme("darkly"), # Add a theme
    titlePanel("Lidar Processing Switches"),
    fluidRow(
      splitLayout(
        cellWidths = c("50%", "50%"),
        wellPanel(
          radioButtons("dap", "Is Acqusition lidar or DAP?", choices = list("lidar" = FALSE, "DAP" = TRUE), selected = FALSE),
          checkboxInput("parallel", "Run in parallel?", FALSE),
          numericInput("cores", "Number of cores:", 1, min = 1, max = 16),
          checkboxInput("tile", "Tile area?", FALSE),
          numericInput("tile_size", "Tile size (m):", 10, min = 1, max = 1000),
          numericInput("chunk_buf", "Chunk buffer:", 10, min = 1, max = 100),
          checkboxInput("ground_classify", "Classify ground points?", FALSE),
          checkboxInput("normalize", "Normalize points?", FALSE),
          checkboxInput("filter_normalize", "Filter outlier normalized returns?", FALSE)
        ),
        wellPanel(
          checkboxInput("dsm", "Create DSM?", FALSE),
          numericInput("dsm_res", "DSM resolution (m):", 1, min = 0.01, max = 100),
          checkboxInput("chm", "Create CHM?", FALSE),
          numericInput("chm_res", "CHM resolution (m):", 1, min = 0.01, max = 100),
          checkboxInput("dtm", "Create DTM?", FALSE),
          numericInput("dtm_res", "DTM resolution (m):", 1, min = 0.01, max = 100),
          checkboxInput("mets", "Calculate Metrics?", TRUE),
          numericInput("met_res", "Metrics resolution (m):", 1, min = 1, max = 100),
          checkboxInput("is_als", "Is Acqusition ALS?", FALSE),
          checkboxInput("is_mls", "Is Acqusition MLS?", FALSE)
        )
      ),
      actionButton("confirm", "Confirm Settings"),
      tableOutput("result")
    )
  )

  # Define server logic
  server <- function(input, output, session) {
    options <- reactive({
      list(
        is_dap = as.logical(input$dap), run_parallel = as.logical(input$parallel),
        num_cores = as.numeric(input$cores), make_tile = as.logical(input$tile),
        tile_size = as.numeric(input$tile_size), chunk_buf = as.numeric(input$chunk_buf),
        ground_classify = as.logical(input$ground_classify), normalize = as.logical(input$normalize),
        filter_normalize = as.logical(input$filter_normalize), make_dsm = as.logical(input$dsm),
        dsm_res = as.numeric(input$dsm_res), make_chm = as.logical(input$chm),
        chm_res = as.numeric(input$chm_res), make_dtm = as.logical(input$dtm),
        dtm_res = as.numeric(input$dtm_res), make_mets = as.logical(input$mets),
        met_res = as.numeric(input$met_res), is_als = as.logical(input$is_als),
        is_mls = as.logical(input$is_mls)
      )
    })

    observeEvent(input$confirm, {
      # Assign options to global environment
      lapply(names(options()), function(opt) assign(opt, options()[[opt]], envir = .GlobalEnv))

      # Pretty print settings in the console with colors
      sapply(names(options()), function(opt) {
        if(is.logical(options()[[opt]])) {
          if(options()[[opt]]) {
            cat(crayon::green(opt), "Enabled\n")
          } else {
            cat(crayon::red(opt), "Disabled\n")
          }
        } else {
          cat(opt, options()[[opt]], "\n")
        }
      })

      # Close the app
      stopApp()
    })

    output$result <- renderTable({
      opt = names(options())
      vals = ifelse(unlist(options()), "Enabled", "Disabled")
      len = length(opt)
      half = ceiling(len / 2)
      data.frame(
        Option1 = opt[1:half], Value1 = vals[1:half],
        Option2 = c(opt[(half + 1):len], rep(NA, 2 * half - len)),
        Value2 = c(vals[(half + 1):len], rep(NA, 2 * half - len)),
        stringsAsFactors = FALSE
      )
    }, rownames = FALSE)
  }

  # Run the application
  shinyApp(ui = ui, server = server)
}
