install.packages("devtools")

# Install bioconductor packages used in chipVis:
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite('bamsignals')
biocLite('rtracklayer')
biocLite('GenomicRanges')
biocLite('Gviz')
biocLite('annotate')

library(devtools)

# Library used in chipVis from github:
devtools::install_github("AliciaSchep/iheatmapr")

# Libraries used in chipVis from CRAN:
install.packages("shiny")
install.packages("plotly")
install.packages("fastcluster")
install.packages("stringr")

# Libraries used in development:
install.packages("roxygen2")
install.packages("testthat")
install.packages("lintr")
install.packages("formatR")
install.packages("covr")
install.packages("DT")
install.packages("remotes")