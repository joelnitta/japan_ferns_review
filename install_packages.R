# This script writes renv.lock for installing R packages to a docker image.
# It should be run from within the rocker/verse:3.6.0 container.
# e.g., after launching the container with
# docker run -it --rm -e DISABLE_AUTH=true -v /Users/joelnitta/repos/japan_ferns_review:/home/rstudio/project rocker/verse:3.6.0 bash
#
# Then build the image with
# docker build . -t joelnitta/japan_ferns_review

### Initialize renv ###

# Install renv
install.packages("remotes", repos = "https://cran.rstudio.com/")
remotes::install_github("rstudio/renv")

# Initialize renv, but don't let it try to find packages to install itself.
renv::init(
  bare = TRUE,
  force = TRUE,
  restart = FALSE)

renv::activate()

### Setup repositories ###

# Install packages that install packages.
install.packages("remotes", repos = "https://cran.rstudio.com/")

# Set repos.
my_repos <- vector()
my_repos["CRAN"] <- "https://cran.rstudio.com/"
options(repos = my_repos)

### Install CRAN packages ###
cran_packages <- c(
  "conflicted",
  "ape",
  "assertr",
  "assertthat",
  "broom",
  "drake",
  "fs",
  "here",
  "janitor",
  "magrittr",
  "multcomp",
  "maps",
  "picante",
  "readxl",
  "scales",
  "scico",
  "tidyverse")

install.packages(cran_packages)

### Install github packages ###
github_packages <- c(
  "joelnitta/jntools",
  "thomasp85/patchwork",
  "danlwarren/RWTY"
)

remotes::install_github(github_packages)

### Take snapshot ###

renv::snapshot(type = "simple")
