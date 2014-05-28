#!/usr/bin/Rscript

if(!require("pdgControl")){
  install.packages("devtools")
  library("devtools")
  install_github("reshape")
  install_github("rmarkdown", "rstudio")
  install_github("cboettig/cboettigR")
  install_github("cboettig/pdg_control", dependencies = c("Depends", "Imports", "Suggests"))
}
