Pretty Darn Good Control
========================

This repository provides an R package with general purpose stochastic dynamic programming code developed
and used in this (draft) manuscript.


Getting started editing manuscript
----------------------------------

- Download and unzip the source file locally
- Change into the manuscripts directory (That is, set `manuscripts/` as your working directory)
- Run the install script
- Restore the cache (optional, avoids waiting for R code to re-run locally)
- Compile the `manuscript.Rmd` document into a pdf (excutes the R code and applies a LaTeX template)[^1].

These steps can be done from within R by pasting in the following code[^2]:

```coffee
download.file("https://github.com/cboettig/pdg_control/archive/master.zip", "pdg_control.zip", "wget")
unzip("pdg_control.zip")
setwd("pdg_control-master/manuscripts")

source("components/install.R")
source("components/restore-cache.R")
rmarkdown::render("manuscript.Rmd")

```

[^1]: Requires [RStudio >= 0.98b](http://www.rstudio.com/ide/download/preview)

[^2]: Windows users may have to download the zip file manually.


Editing
-------

Use RStudio to edit the file `manuscript.Rmd` in the `pdg_control/manuscripts` directory.
