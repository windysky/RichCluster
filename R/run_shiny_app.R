#' Launch RichStudio shiny app
#'
#' @import Rcpp
#' @export
launch_shiny <- function() {
  appDir <- system.file("application", package = "RichCluster")
  if (appDir == "") {
    stop("Could not find application. Try re-installing `RichCluster`.", call. = FALSE)
  }

  # Source all R scripts from subdirectories
  sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "\\.[Rr]$", full.names = TRUE)) {
      if(trace) cat(nm,":")
      source(nm, ...)
      if(trace) cat("\n")
    }
  }

  # Source scripts in subdirectories
  sourceDir(system.file("R/modules", package = "RichCluster"))

  shiny::runApp(appDir, display.mode = "normal")
}
