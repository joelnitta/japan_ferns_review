# make.R

# Master script for running analyses.
# This project uses drake to manage workflows.
# For more information about drake, see
# https://ropensci.github.io/drake/

# Setup ----

# For reproducibility
set.seed(9130)

# Set working directory
setwd(here::here())

# Load packages
source("code/packages.R")

# Update drake settings
pkgconfig::set_config("drake::strings_in_dots" = "literals")

# Load functions and plans  ----
source("code/functions.R")
source("code/plan.R")

# Run analyses ----
make(plan) 
