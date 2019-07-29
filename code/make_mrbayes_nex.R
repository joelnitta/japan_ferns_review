# Read in rbcL alignment, set up for analysis with MrBayes, and write
# out nexus file formatted for MrBayes. The nexus file can then be
# uploaded to CIPRES for phylogenetic analysis.
#
# This is performed separately from the rest of the drake plan because
# requires manual uploading and downloading.

library(ape)
library(glue)
library(tidyverse)

# Read in alignment ----

# The original alignment file doesn't seem to be in proper phylip format.
# Read in as raw text after skipping first two lines (matrix dimensions and
# a blank line), then write out as fasta and read back in with ape.

rbcl_df <- read_delim(
  "data/JpFern_rbcL.txt", 
  skip = 2, delim = " ", 
  col_names = c("name", "seq", "extra")) %>%
  select(name, seq) %>%
  mutate(
    name = paste0(">", name),
    line = map2(name, seq, c))

temp_file <- tempfile()

write_lines(unlist(rbcl_df$line), temp_file)

rbcl <- ape::read.FASTA(temp_file)

# Set up constriaints ----

# Make tibble of constraint commands for MrBayes.
# Each constraint is for a PPGI fern family with three or more taxa.
# (the first part of the sequence name separated by
# underscores is the family).
constraints <-
  tibble(
    tip = names(rbcl)
  ) %>%
  mutate(
    family = map_chr(tip, ~ str_split(., "_") %>% map_chr(1))
  ) %>%
  group_by(family) %>%
  summarize(
    tips = paste(tip, collapse = " "),
    n = n()
  ) %>%
  filter(n > 2) %>%
  mutate(
    family = paste0("f", family),
    constraint = glue("constraint {family} = {tips};")
  )

# Define MrBayes block ----
mrbayes_block <- c(
  "BEGIN mrbayes;",
  "set autoclose=no nowarn=yes;",
  constraints$constraint,
  "outgroup 601_1;",
  "prset statefreqpr=dirichlet(1.0,1.0,1.0,1.0);",
  paste0(
    "prset topologypr = constraints (", 
    paste(constraints$family, collapse = ", "), 
    ");"),
  "lset nst=6 rates=invgamma; [GTR+I+G]",
  "mcmcp ngen=1000000 nruns=2 nchains=4 samplefreq = 1000 printfreq = 1000 temp = 0.005 relburnin = yes burninfrac = 0.25;",
  "mcmc;",
  "sumt conformat = simple;",
  "END;"
)

# Write out nexus file for MrBayes ----

# Write out the alignment in nexus format
write.nexus.data(rbcl, "data/rbcl_mrbayes.nex", interleaved = FALSE)

# Append mrbayes block to nexus
write_lines(mrbayes_block, "data/rbcl_mrbayes.nex", append = TRUE)
