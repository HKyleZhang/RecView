#' Run RecView ShinyApp
#'
#' @description Initiate the RecView ShinyApp. The RecView ShinyApp is intended to visualize and provide options automatically locate recombinations.
#' @usage run_RecView_App()
#' @note Before execute run_RecView_App(), it is recommmended to navigate the working directory to the directory of input files. This makes it convenient for saving result output files from the app.
#'
#' @export
run_RecView_App <- function() {
  options(shiny.maxRequestSize = 1000*1024^2)
  #shiny::shinyApp(RecViewR:::recview_ui, RecViewR:::recview_server)
  shiny::runGadget(RecViewR:::recview_ui, RecViewR:::recview_server, viewer = shiny::paneViewer(minHeight = 800))
}

