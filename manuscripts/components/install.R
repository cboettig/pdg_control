#!/usr/bin/Rscript
library("methods")

# Install packrat if not available
if(!require("devtools")) install.packages("devtools")
if(!require("packrat")) devtools::install_github("rstudio/packrat")

# Let packrat set up the manuscript dependencies
packrat::packify()
source(".Rprofile"); readRenviron(".Renviron")
packrat::restore()
source(".Rprofile"); readRenviron(".Renviron")


