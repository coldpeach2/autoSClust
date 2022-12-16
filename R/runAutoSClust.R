#' @importFrom shiny runApp

runAutoSClust <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "autoSClust")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}
