#!/usr/bin/Rscript
library("methods")
# Let packrat set up the manuscript dependencies
if(file.exists("packrat.sources") && file.exists(".Rprofile") && file.exists(".Renviron")){
  if(exists("initPackrat"))
    initPackrat()
  else 
    NULL
} else {
# Install packrat if not available
  if(!require("devtools")) install.packages("devtools")
  if(!require("packrat")) devtools::install_github("rstudio/packrat")
  packrat::packify()
  packrat::restore()
}
