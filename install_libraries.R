install.packages("devtools")

# Install bioconductor packages used in chipVis:
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite('bamsignals')
biocLite('rtracklayer')
biocLite('GenomicRanges')
biocLite('Gviz')
biocLite('annotate')
biocLite("S4Vectors")
biocLite("biovizBase")
biocLite("OrganismDbi")

library(devtools)

# Library used in chipVis from github:
devtools::install_github("AliciaSchep/iheatmapr@e4c4bdd54bde09643fc350f0d4df951f040b5267")

# Libraries used in chipVis from CRAN:
install.packages("shiny")
install.packages("plotly")
install.packages("fastcluster")
install.packages("stringr")
install.packages("RColorBrewer")
install.packages("dplyr")

# Libraries used in development:
install.packages("roxygen2")
install.packages("testthat")
install.packages("lintr")
install.packages("formatR")
install.packages("covr")
install.packages("DT")
install.packages("remotes")