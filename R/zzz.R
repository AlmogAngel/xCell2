.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("pracma", quietly = TRUE)) {
    message("Installing 'pracma' from R-Forge...")
    install.packages("pracma", repos = "http://R-Forge.R-project.org")
  }
  
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("The 'quadprog' package is required but not installed. It should install automatically from CRAN. Please install it manually if this message persists.")
  }
}