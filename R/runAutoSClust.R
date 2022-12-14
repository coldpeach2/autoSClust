#' @importFrom shiny runApp

runClustREval <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "autoSClust")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}
