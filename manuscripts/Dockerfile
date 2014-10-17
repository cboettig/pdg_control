## Launch an R terminal: docker run -it cboettig/pdg_control /usr/bin/R
## Launch via RStudio: docker run -d -p 8787:8787 cboettig/pdg_control
##   and then navigate to http://localhost:8787 

FROM rocker/hadleyverse
MAINTAINER Carl Boettiger cboettig@ropensci.org

ENV PACKAGE pdg_control

## Remain current & install binaries for dependencies that don't compile
RUN apt-get update \ 
  && apt-get install -y r-cran-snowfall r-cran-expm

## Install package, with dependencies and suggests
RUN installGithub.r --deps TRUE \ 
  cboettig/cboettigR \
  jrnold/ggthemes \ 
  hadley/reshape \ 
  cboettig/pdg_control

# Add manuscripts directory onto image 
COPY . /host/pdg_control/manuscripts

# Check image
WORKDIR /host/pdg_control/manuscripts


## Test build but leave image clean
#RUN nice -n 19 Rscript -e 'rmarkdown::render("manuscript.Rmd")' \
#  && rm -rf components/cache/*

## Note that cboettig/rstudio's CMD uses supervisord to first run a shell script which will make sure 
## that the user's whole home directory, including the subdirectories we just added (as root) 
## have ownership permissions assigned to the user specified in --e USER=<username> (rstudio by default). 
