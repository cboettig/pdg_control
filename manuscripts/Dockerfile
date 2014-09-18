## Launch an R terminal: docker run -it cboettig/pdg_control /usr/bin/R
## Launch via RStudio: docker run -d -p 8787:8787 cboettig/pdg_control
##   and then navigate to http://localhost:8787 

FROM cboettig/ropensci:latest
MAINTAINER Carl Boettiger cboettig@ropensci.org

ENV PACKAGE pdg_control

## Remain current & install binaries for dependencies that don't compile
RUN apt-get update -qq && apt-get dist-upgrade -qy && apt-get install -qy r-cran-snowfall r-cran-expm 

## Install package, with dependencies and suggests
RUN Rscript -e 'devtools::install_github(c("cboettig/cboettigR", "jrnold/ggthemes", "hadley/reshape", "cboettig/pdg_control"), dependencies=TRUE)'

# Add manuscripts directory onto image 
ADD . /host/pdg_control/manuscripts

# Check image
WORKDIR /host/pdg_control/manuscripts
RUN rm -rf components/cache/*
RUN Rscript -e 'rmarkdown::render("manuscript.Rmd")'

## Note that cboettig/rstudio's CMD uses supervisord to first run a shell script which will make sure 
## that the user's whole home directory, including the subdirectories we just added (as root) 
## have ownership permissions assigned to the user specified in --e USER=<username> (rstudio by default). 
