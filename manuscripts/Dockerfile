## Launch an R terminal: docker run -it cboettig/pdg_control /usr/bin/R
## Launch via RStudio: docker run -d -p 8787:8787 cboettig/pdg_control
##   and then navigate to http://localhost:8787 

FROM cboettig/ropensci:latest
MAINTAINER Carl Boettiger cboettig@ropensci.org

ENV PACKAGE pdg_control

## Remain current & install dependencies
RUN apt-get update -qq && apt-get dist-upgrade -qy && apt-get install -qy texlive-humanities lmodern texlive-fonts-recommended texlive-latex-extra r-cran-snowfall r-cran-expm

## Install package, with dependencies and suggests
RUN Rscript -e 'devtools::install_github(c("cboettig/cboettigR", "jrnold/ggthemes", "cboettig/pdg_control"))'

RUN apt-get install -y imagemagick

ADD . /home/rstudio/pdg_control/manuscripts
RUN chmod -R 777 /home/rstudio/pdg_control/manuscripts
WORKDIR /home/rstudio/pdg_control/manuscripts

