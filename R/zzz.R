.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("pracma", quietly = TRUE)) {
    packageStartupMessage("Installing 'pracma' from R-Forge...")
    utils::install.packages("pracma", repos = "http://R-Forge.R-project.org")
  }
  
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("The 'quadprog' package is required but not installed. Please install it manually from CRAN.")
  }
}

.onAttach <- function(libname, pkgname) {
  if (!requireNamespace("pracma", quietly = TRUE)) {
    packageStartupMessage("The 'pracma' package is required and will be installed from R-Forge.")
    utils::install.packages("pracma", repos = "http://R-Forge.R-project.org")
  }
}