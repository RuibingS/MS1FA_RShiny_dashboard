cran_packages <- c(
  "shiny", "DT", "shinyWidgets", "dplyr", "igraph", "stringr",
  "readxl", "purrr", "readr", "plyr", "data.table",
  "tidyverse", "hrbrthemes", "viridis", "viridisLite", "ggplot2", "roxygen2",
  "rlang", "RcppArmadillo", "webshot", "htmlwidgets", "profvis", "shinythemes",
  "shinyjs", "visNetwork", "bs4Dash","magick","chromote"
)
bioconductor_packages <- c(
  "MetaboCoreUtils", "enviPat","MSnbase"
)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

install_cran_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg)
      print(paste0("Please install the required package: ", pkg))
      # After installation, load the package
      library(pkg, character.only = TRUE)
      
    } else {
      # If already installed, just make sure it's loaded
      library(pkg, character.only = TRUE)
    
    }
  }
}
install_cran_packages(packages=cran_packages)

# Function to install Bioconductor packages
install_bioconductor_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      print(paste0("Please install the required package: ", pkg)) 
      BiocManager::install(pkg)
      library(pkg, character.only = TRUE)
      
    }
    else {
      # If already installed, just make sure it's loaded
      library(pkg, character.only = TRUE)
     
    }
  }
}
install_bioconductor_packages(packages=bioconductor_packages)


