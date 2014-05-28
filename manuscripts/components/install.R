#!/usr/bin/Rscript

if(!require("pdgControl")){
  install.packages("devtools")
  library("devtools")
  install_github("reshape")
  install_github("rmarkdown", "rstudio")
  install_github("cboettig/cboettigR")


  install.packages("..", repos = NULL, dependencies = c("Depends", "Imports", "Suggests"))
}
