# This script writes packrat.lock for installing R packages to a docker image.
# It should be run from within the rocker/verse:3.6.0 container.
# e.g., after launching the container with
# docker run -d --rm -e DISABLE_AUTH=true -v /Users/joelnitta/repos/japan_ferns_review:/home/rstudio/project rocker/verse:3.6.0
#
# Then build the image with
# docker build . -t joelnitta/japan_ferns_review
#
# For more info on installing R packages to docker images with
# packrat, see https://www.joelnitta.com/post/docker-and-packrat/

### Initialize packrat ###

# Don't let packrat try to find
# packages to install itself.

install.packages("packrat", repos = "https://cran.rstudio.com/")
packrat::init(
  infer.dependencies = FALSE,
  enter = TRUE,
  restart = FALSE)

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
  "thomasp85/patchwork"
)

remotes::install_github(github_packages)

### Take snapshot ###

packrat::snapshot(
  snapshot.sources = FALSE,
  ignore.stale = TRUE,
  infer.dependencies = FALSE)
