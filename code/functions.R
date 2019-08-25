# Data processing ----

#' Clean up reproductive mode data for pteridophytes of Japan
#'
#' @param data Tibble. Raw data read-in from Electronic Supp. Mat. 1
#'
#' @return Tibble
process_repro_data <- function (data) {
  
  data %>%
    clean_names %>%
    mutate(
      reproductive_mode = case_when(
        reproductive_mode == 0 ~ "unknown",
        reproductive_mode == 1 ~ "sexual", 
        reproductive_mode == 2 ~ "apomictic",
        reproductive_mode == 3 ~ "sex_apo",
        TRUE ~ "unknown"
      ) %>% as.factor,
      sexual_diploid = case_when(
        sexual_diploid == 1 ~ TRUE,
        TRUE ~ FALSE
      ),
      sexual_polyploid = case_when(
        sexual_polyploid == 1 ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    select(taxon_id, reproductive_mode, sexual_diploid, sexual_polyploid) 
  
}

#' Read in a nexus file contained in a zipped archive
#'
#' @param zip_folder Path to zip file
#' @param nexus_file Name of nexus file within zip file
#'
#' @return List
#' 
read_nexus_in_zip <- function (zip_folder, nexus_file) {
  
  temp_dir <- tempdir()
  
  unzip(zip_folder, exdir = temp_dir)
  
  ape::read.nexus(fs::path(temp_dir, nexus_file))

}

# Basic stats ----

#' Count species per grid cell excluding hybrids
#'
#' @param occ_data Occurrence data, with one row per
#' grid cell per taxon, including hybrids.
#' @param repro_data Reproductive mode data, with
#' one row per taxon, excluding hybrids.
#' @return tibble
count_species_per_cell <- function (occ_data, repro_data) {
  occ_data %>%
    filter(taxon_id %in% repro_data$taxon_id) %>%
    group_by(secondary_grid_code) %>%
    count(sort = TRUE) %>%
    ungroup
}

#' Count grid cells per species
#'
#' @param occ_data Occurrence data, with one row per
#' grid cell per taxon, including hybrids.
#' @return tibble
count_cells_per_species <- function (occ_data) {
  occ_data %>%
    group_by(taxon_name) %>%
    count(sort = TRUE) %>%
    ungroup
}

#' Count number of grid cells per species by reproductive mode
#'
#' @param occ_data Occurrence data, with one row per
#' grid cell per taxon, including hybrids.
#' @param repro_data Reproductive mode mata, with
#' one row per taxon, excluding hybrids.
#' @return tibble
count_cells_per_species_by_repro <- function(occ_data, repro_data) {
  occ_data %>%
    group_by(taxon_id) %>%
    summarize(
      n_grids = n()
    ) %>%
    inner_join(select(repro_data, taxon_id, reproductive_mode)) %>%
    ungroup() %>%
    filter(reproductive_mode != "unknown")
}

#' Count number of grid cells per species by ploidy level
#' 
#'
#' @param occ_data Occurrence data, with one row per
#' grid cell per taxon, including hybrids.
#' @param repro_data Reproductive mode mata, with
#' one row per taxon, excluding hybrids.
#' @return tibble
count_cells_per_species_by_ploidy <- function(occ_data, repro_data) {
  occ_data %>%
    group_by(taxon_id) %>%
    summarize(
      n_grids = n()
    ) %>%
    inner_join(select(repro_data, taxon_id, reproductive_mode, sexual_diploid, sexual_polyploid)) %>%
    ungroup()
}

#' Count number of grid cells per species by growth type
#' 
#' @param occ_data Occurrence data, with one row per
#' grid cell per taxon, including hybrids.
#' @param growth_data Growth mode mata, with
#' one row per taxon, excluding hybrids.
#' @param cells_per_species Number of grid cells per species
#' 
#' @return tibble
count_cells_per_species_by_growth <- function(occ_data, growth_data, cells_per_species) {
  
  # Make table of taxon IDs and names for occurrence data
  occ_taxa <-
    occ_data %>%
    select(taxon_id, taxon_name) %>%
    unique
  
  # Optional: check for missing species.
  # These are only hybrid taxa
  in_occ_missing_from_growth <-
    anti_join(occ_taxa, growth_data)
  
  # Join growth data with cells per species.
  # Note that some species are both, so these 
  # will be repeated.
  growth_data %>%
    # Make sure all taxon IDs are in the occurrence data
    verify(all(taxon_id %in% occ_taxa$taxon_id)) %>%
    left_join(growth_data) %>%
    left_join(occ_taxa) %>%
    left_join(cells_per_species) %>%
    rename(n_grids = n)
  
}

#' Get mean number of grid cells per species by reproductive mode
#'
#' @param cells_per_species_by_repro Tibble
#'
#' @return Tibble
avg_cells_per_species_by_repro <- function (cells_per_species_by_repro) {
  cells_per_species_by_repro %>%
    group_by(reproductive_mode) %>%
    summarize(
      mean = mean(n_grids, na.rm = TRUE),
      n = n(),
      sd = sd(n_grids, na.rm = TRUE)
    )
}

#' Get mean number of grid cells per species by ploidy level
#'
#' @param cells_per_species_by_ploidy Tibble 
#'
#' @return Tibble
avg_cells_per_species_by_ploidy <- function (cells_per_species_by_repro) {
  bind_rows(
    cells_per_species_by_repro %>%
    filter(reproductive_mode == "sexual") %>%
    group_by(sexual_diploid) %>%
    summarize(
      mean = mean(n_grids, na.rm = TRUE),
      n = n(),
      sd = sd(n_grids, na.rm = TRUE)
    ),
    cells_per_species_by_repro %>%
      filter(reproductive_mode == "sexual") %>%
      group_by(sexual_polyploid) %>%
      summarize(
        mean = mean(n_grids, na.rm = TRUE),
        n = n(),
        sd = sd(n_grids, na.rm = TRUE)
      ),
  )
}

#' Determine latitudinal breadth for each species,
#' including reproductive mode
#'
#' @param occ_data Occurrence data, with one row per
#' grid cell per taxon, including hybrids.
#' @param repro_data Reproductive mode mata, with
#' one row per taxon, excluding hybrids.
#'
#' @return Tibble. One species per row with latitudinal
#' breadth (max - min) and reproductive mode, excluding
#' hybrids and repro mode unknown.
count_lat_by_repro <- function(occ_data, repro_data) {
  occ_data %>%
    group_by(taxon_id) %>%
    summarize(
      lat_breadth = max(latitude) - min(latitude)
    ) %>%
    inner_join(select(repro_data, taxon_id, reproductive_mode)) %>%
    ungroup() %>%
    filter(reproductive_mode != "unknown")
}

#' Get mean latitudinal breadth per species by reproductive mode
#'
#' @param cells_per_species_by_repro 
#'
#' @return Tibble
avg_lat_by_repro <- function(lat_by_repro) {
  lat_by_repro %>%
    group_by(reproductive_mode) %>%
    summarize(
      mean = mean(lat_breadth, na.rm = TRUE),
      n = n(),
      sd = sd(lat_breadth, na.rm = TRUE)
    )
}

# Taxonomy ----

#' Modify Pteridophyte Phylogeny Group I taxonomy
#' to match the version used for Pteridophytes of Japan
#'
#' @param ppgi Original PPGI taxonomy
#'
#' @return tibble
modify_ppgi <- function (ppgi) {
  
  # Use normal "e" for Isoetes
  ppgi_mod <- 
    ppgi %>%
    mutate_all(~str_replace_all(., "ë", "e"))
  
  # Add genera missing in PPGI that are included in Japan pteridophyte checklist
  # Use the sister (or encompassing) genus for each, so other higher-order
  # taxonomy will be correct
  anisocampium_dat <-
    ppgi_mod %>% filter(genus == "Athyrium") %>%
    mutate(genus = "Anisocampium")
  
  humata_dat <-
    ppgi_mod %>% filter(genus == "Davallia") %>%
    mutate(genus = "Humata")
  
  drynaria_dat <-
    ppgi_mod %>% filter(genus == "Aglaomorpha") %>%
    mutate(genus = "Drynaria")
  
  bind_rows(ppgi_mod, anisocampium_dat, humata_dat, drynaria_dat) %>%
    select(genus, family, order, class)
  
}

#' Add taxonomy data to occurrence data
#'
#' @param occ_data Occurence data, including
#' at least one column called "taxon_name" where
#' the first word separated by spaces is the genus name.
#' @param taxonomy_data Taxonomy data at the
#' genus level and higher
#'
#' @return tibble
add_taxonomy <- function(occ_data, taxonomy_data) {
  occ_data %>%
    mutate(genus = str_split(taxon_name, " ") %>% map_chr(1)) %>%
    filter(!is.na(genus)) %>%
    left_join(taxonomy_data)
}

#' Summarize taxonomic data
#'
#' @param repro_data Reproductive mode data, with
#' one row per taxon, excluding hybrids.
#' @param occ_data_with_taxonomy Occurrence data, with one row per
#' grid cell per taxon, including hybrids. Should already have
#' higher-level taxonomy (genus, family, etc) included.
#' 
#' @return tibble: one row per taxon with higher level taxonomy.
summarize_taxonomy <- function(repro_data, occ_data_with_taxonomy) {
  
  # Just get relevant parts of taxonomic data for joining to repro data
  occ_data_with_taxonomy <-
  select(occ_data_with_taxonomy, 
         taxon_id, taxon_name, 
         genus, family, order, class) %>% unique
  
  taxonomic_data <- left_join(repro_data, occ_data_with_taxonomy) %>%
    # Add column for lowest taxonomic rank
    mutate(lowest_taxonomic_rank = case_when(
      str_detect(taxon_name, "var\\.") ~ "infraspecies",
      str_detect(taxon_name, "subsp\\.") ~ "infraspecies",
      TRUE ~ "species"
    )) %>%
    # Add column for species (without infra sp. taxon)
    mutate(species_name = str_split(taxon_name, " ") %>% 
             map_chr(., ~magrittr::extract(., 1:2) %>% jntools::paste3(collapse = " ")))
  
  # Double check that our assessment of species vs infraspecies was correct:
  # all "species" should have exactly one space in their name
  assert_that(
    taxonomic_data %>% 
      filter(lowest_taxonomic_rank == "species") %>%
      pull(taxon_name) %>%
      map_dbl(~str_count(., " ")) %>%
      unique == 1
  )
  
  return(taxonomic_data)
}

# Breeding system ----

#' Calculate the percentage of apomictic taxa in a given grid cell
#'
#' @param occ_data Occurence data, including
#' a column called "taxon_name" where
#' the first word separated by spaces is the genus name,
#' and "secondary_grid_cell" where that species occurs.
#' @param repro_data Reproductive mode data, with
#' one row per taxon, excluding hybrids.
#' @param all_cells List of all 1km2 grid cells with
#' grid cell code and long/lat.
#'
#' @return tibble
calc_percent_apomictic <- function (occ_data, repro_data, all_cells) {
  inner_join(
    occ_data,
    repro_data
  ) %>%
    filter(!is.na(secondary_grid_code)) %>%
    mutate(is_apo = case_when(
      reproductive_mode == "apomictic" ~ TRUE,
      TRUE ~ FALSE
    )) %>%
    group_by(secondary_grid_code, latitude, longitude) %>%
    summarize(
      num_apomictic = sum(is_apo),
      num_total = n()
    ) %>%
    ungroup %>%
    mutate(percent_apomictic = num_apomictic / num_total) %>%
    select(secondary_grid_code, percent_apomictic) %>%
    # Add 0s for cells with no occurrence data
    right_join(all_cells) %>%
    mutate(percent_apomictic = replace_na(percent_apomictic, 0))
}

# Community diversity ----

#' Calculate species richness for all grid cells
#' 
#' The long / lats of the occurrence data
#' are already cell centroids, so it is trivial
#' to compute species richness.
#'
#' @param occ_data Occurrence data, where each row is the 
#' occurrence of a species in a grid cell
#' @param all_cells Dataframe with long/lat of all 10km2 grid
#' cells across Japan
#' @return tibble
make_richness_matrix <- function (occ_data, all_cells) {
  occ_data %>%
    group_by(secondary_grid_code) %>%
    summarize(
      richness = n()
    ) %>%
    # Add all cells, including those missing from occ_data
    right_join(all_cells) %>%
    mutate(richness = replace_na(richness, 0))
}


#' Format tip labels in Japanese pteridophytes rbcL tree
#'
#' The tips in original phylogeny file are coded with two numbers 
#' separated by an underscore,
#' e.g., 601_14.
#' The second part is the taxon_id in occurrence and reproductive 
#' mode data. Keep only this as the tip label.
#'
#' @param phy Phylogeny of Japanese pteridophytes with tip labels
#' formatted as two numbers separated by an underscore
#'
#' @return List of class "phylo"
format_tip_labels <- function (phy) {
  
  phy$tip.label <- str_split(phy$tip.label, "_") %>% map_chr(2)
  
  phy
  
}

#' Make a community matrix
#'
#' @param occ_data Occurrence data, with one row per
#' grid cell per taxon, including hybrids.
#' 
#' @return tibble. One column for species then the rest
#' of the columns are presence/absence of that species in 
#' each site, where a "site" is a 1km2 grid cell. Names
#' of sites are grid-cell codes.
#' Species are stored as taxon_id values.
make_comm_matrix <- function (occ_data) {
  occ_data %>%
    select(species = taxon_id, site = secondary_grid_code) %>%
    mutate(
      abundance = 1,
      site = as.character(site)) %>%
    spread(site, abundance) %>%
    mutate_at(vars(-species), ~replace_na(., 0)) %>%
    mutate(species = as.character(species))
}

#' Calculate Faith's PD for a single community
#' 
#' Does not include the root when summing branch lengths.
#' 
#' At least two taxa in the community must be present in the tree; otherwise,
#' communities with 0 species will return 0,
#' communities with 1 species will return NA.
#' 
#' @param single_comm Single community, formatted as tibble with
#' two columns: 'species' (character) and 'abundance' (numeric)
#' @param phy Phylogenetic tree
#' @param shuffle_tips Logical; should tips of the phylogeny be
#' randomized? Used for generating null distributions.
#' 
#' @return a number: phylogenetic diversity (sum of branch lengths)
#' for that community
calc_pd <- function(single_comm, phy, shuffle_tips = FALSE) {
  
  # Filter to only species present in the focal community
  single_comm <- filter(single_comm, abundance > 0)
  
  # Return 0 if there are zero species present
  if(nrow(single_comm) == 0) return (0)
  
  # The phylogenetic distance for a single species without using the root is
  # undefined.
  if(nrow(single_comm) == 1) return (NA)
  
  # Optionally shuffle tips when generating null distributions
  if(isTRUE(shuffle_tips)) phy <- picante::tipShuffle(phy)
  
  # Prune tree to only species present in the community
  phy <- ape::keep.tip(phy, single_comm$species)
  
  # Get sum of branches remaining
  sum(phy$edge.length)
}

#' Calculate Faith's PD across multiple communities.
#' 
#' Does not include the root when summing branch lengths. At least two taxa in each
#' community must be present in the tree.
#' 
#' @param comm community matrix. One column must be named
#' 'species', and the rest should correspond to presence or absence of species
#' in communities (sites).
#' @param phy phylogenetic tree
#' 
#' @return a data frame with observed Faith's PD (pd_obs)
calc_pd_comm <- function (comm, phy) {
  
  assert_that(isTRUE(all.equal(
    sort(comm$species), sort(phy$tip.label) )),
    msg = "Species don't match exactly between 'comm' and 'phy'"
  )
  
  # Nest by community, then calculate PD for each
  comm %>%
    gather(secondary_grid_code, abundance, -species) %>%
    nest(-secondary_grid_code) %>%
    mutate(
      pd_obs = map_dbl(data, ~ calc_pd(., phy = phy))
    ) %>%
    select(-data)
  
}


#' Match community data and tree
#' 
#' Order of species in comm will be rearranged to match the
#' phylogeny.
#'
#' @param comm Community data frame, with one column for sites and
#' the rest for species.
#' @param phy Phylogeny (list of class "phylo")
#' @param return Type of object to return
#'
#' @return Either a dataframe or a list of class "phylo"; the tree or
#' the community, pruned so that only species occurring in both datasets
#' are included.
#' @export
#'
#' @examples
match_comm_and_tree <- function (comm, phy, return = c("comm", "tree")) {
  
  assert_that("species" %in% colnames(comm))
  
  # Keep only species in phylogeny
  comm <- comm %>%
    filter(species %in% phy$tip.label) 
  
  # Trim to only species with trait data
  phy <- drop.tip(phy, setdiff(phy$tip.label, comm$species))
  
  # Get comm in same order as tips
  comm <- left_join(
    tibble(species = phy$tip.label),
    comm
  )
  
  # Make sure that worked
  assert_that(isTRUE(all.equal(comm$species, phy$tip.label)))
  
  # Return comm or tree
  assert_that(return %in% c("tree", "comm"))
  
  if(return == "tree") { 
    return (phy) 
  } else {
    return (comm)
  }
  
}

# Plotting ----

#' Define ggplot theme
#' 
#' BW theme with no gridlines, black axis text, main font size 11,
#' axis ticks size 9.
#'
standard_theme2 <- function () {
  ggplot2::theme_bw() + 
    theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(colour="black"),
      axis.text.y = ggplot2::element_text(colour="black")
    )
}

#' Get the lower, upper, or absolute maximum value
#' of a variable in a dataframe
#' 
#' For setting plotting limits manually
#'
#' @param data Dataframe
#' @param var Name of variable (column) in dataframe
#' @param digits Number of digits desired in output
#' @param type Type of limit to calculate: "min" is lower,
#' "max" is upper, and "abs" is the absolute greatest value.
#'
#' @return Number
#'
#' @examples
#' get_limit(mtcars, disp, "max")
get_limit <- function (data, var, type = c("min", "max", "abs"), digits = 2) {
  
  var <- enquo(var)
  
  switch(type,
         
         max = data %>% 
           pull(!!var) %>% 
           max(na.rm = TRUE) %>% 
           multiply_by(10^digits) %>% ceiling %>% divide_by(10^digits),
         
         min = data %>% 
           pull(!!var) %>% 
           min(na.rm = TRUE) %>% 
           multiply_by(10^digits) %>% floor %>% divide_by(10^digits),
         
         abs = c(
           data %>% pull(!!var) %>% max(na.rm = TRUE),
           data %>% pull(!!var) %>% min(na.rm = TRUE)) %>% 
           abs %>% max %>%
           multiply_by(10^digits) %>% ceiling %>% divide_by(10^digits)
  )
  
}

#' Make a plot showing selected alpha diversity metric on a map of Japan
#'
#' @param div_data Alpha diversity matrix; rows are communities
#' (1km2 grid cells), and columns are various alpha diversity metrics.
#' @param world_map Background world mapp
#' @param occ_data Occurrence data, with one row per
#' grid cell per taxon, including hybrids.
#' @param div_metric Selected diversity metric to plot.
#' Must one of the column names of div_dat.
#' @param metric_title Character: title to use for legend for
#' diversity metric.
#'
#' @return ggplot object
make_diversity_map <- function (div_data, world_map, occ_data, div_metric, metric_title, label) {
  
  div_metric <- sym(div_metric)
  
  ggplot(world_map, aes(x = longitude, y = latitude)) +
    geom_polygon(aes(group = group), fill = "light grey") +
    geom_tile(data = div_data,
              aes(fill = !!div_metric)) + 
    coord_quickmap(
      xlim = c(pull(occ_data, longitude) %>% min %>% floor, 
               pull(occ_data, longitude) %>% max %>% ceiling),
      ylim = c(pull(occ_data, latitude) %>% min %>% floor, 
               pull(occ_data, latitude) %>% max %>% ceiling)
    )  +
    labs(
      fill = metric_title
    ) +
    annotate(
      "text", label = label, size = 12/.pt, fontface = "bold",
      x = -Inf, 
      y = Inf,
      vjust = 1.2, hjust = -0.5) +
    jntools::blank_x_theme() +
    jntools::blank_y_theme() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.title = element_text(size = 20/.pt),
      legend.text = element_text(size = 16/.pt),
      legend.justification=c(1,0), 
      legend.position=c(1.1,0))
}

#' Make a plot showing selected SES of PD on map of Japan with
#' only significant communities colored.
#'
#' @param div_data Alpha diversity matrix; rows are communities
#' (1km2 grid cells), and columns are various alpha diversity metrics.
#' @param world_map Background world mapp
#' @param occ_data Occurrence data, with one row per
#' grid cell per taxon, including hybrids.
#'
#' @return ggplot object
make_pd_highlight_map <- function (div_data, world_map, occ_data) {
  
  # gghighlight only "knows about" data in the most recent layer
  # So make plot with world map on top of highlighted SES of PD,
  # then rearrange layers.
  plot <-
    ggplot(div_data, aes(x = longitude, y = latitude)) +
    geom_tile(aes(fill = ses_pd), color = "black") +
    gghighlight( (pd_obs_p > 0.975 | pd_obs_p < 0.025) & !is.na(ses_pd) ) +
    coord_quickmap(
      xlim = c(pull(occ_data, longitude) %>% min %>% floor, 
               pull(occ_data, longitude) %>% max %>% ceiling),
      ylim = c(pull(occ_data, latitude) %>% min %>% floor, 
               pull(occ_data, latitude) %>% max %>% ceiling)
    ) +
    geom_polygon(data = world_map, aes(group = group), fill = "light grey")
  
  # See
  # https://stackoverflow.com/questions/20249653/insert-layer-underneath-existing-layers-in-ggplot2-object
  plot$layers <- plot$layers[c(1,3,2)]
  
  plot +
    scale_fill_scico(palette = "vik", na.value="transparent") +
    jntools::blank_x_theme() +
    jntools::blank_y_theme() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA)
    ) +
    labs(
      fill = "SES of PD"
    )
  
}

#' Assemble a set of jitter plots showing differences in
#' number of grid cells or latitudinal breadth by different
#' categories
#'
#' @param cps_by_repro Dataframe of number of cells per species
#' by reproductive mode
#' @param lat_by_repro Dataframe of latitudinal breadth per species
#' by reproductive mode
#' @param cps_by_ploidy Dataframe of number of cells per species
#' by ploidy level
#'
#' @return GGplot object
#'
#' @examples
assemble_jitter_plots <- function(cps_by_repro, lat_by_repro, cps_by_ploidy) {
  
  b <- cps_by_repro %>%
    mutate(
      reproductive_mode = forcats::fct_recode(
        reproductive_mode,
        Sexual = "sexual",
        `Sex. apo.` = "sex_apo",
        Apomictic = "apomictic"
      )
    ) %>%
    ggplot(aes(x = reproductive_mode, y = n_grids, color = reproductive_mode)) +
    geom_jitter(alpha = 0.7) +
    geom_boxplot(fill = "transparent", color = "dark grey", outlier.shape = NA) +
    labs(
      y = "No. Grid cells",
      x = "Reproductive mode",
      subtitle = "b"
    ) +
    standard_theme2() +
    theme(
      legend.position = "none",
      plot.subtitle = element_text(face = "bold")
    )
  
  c <- lat_by_repro %>%
    mutate(
      reproductive_mode = forcats::fct_recode(
        reproductive_mode,
        Sexual = "sexual",
        `Sex. apo.` = "sex_apo",
        Apomictic = "apomictic"
      )
    ) %>%
    ggplot(aes(x = reproductive_mode, y = lat_breadth, color = reproductive_mode)) +
    geom_jitter(alpha = 0.7) +
    geom_boxplot(fill = "transparent", color = "grey", outlier.shape = NA) +
    labs(
      y = "Lat. breadth (°)",
      x = "Reproductive mode",
      subtitle = "c"
    ) +
    standard_theme2() +
    theme(
      legend.position = "none",
      plot.subtitle = element_text(face = "bold")
    )
  
  d <- cps_by_ploidy %>%
    mutate(ploidy = case_when(
      sexual_diploid == TRUE ~ "Diploid",
      sexual_polyploid == TRUE ~ "Polyploid",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(ploidy)) %>%
    ggplot(aes(x = ploidy, y = n_grids, color = ploidy)) +
    geom_jitter(alpha = 0.7) +
    geom_boxplot(fill = "transparent", color = "grey", outlier.shape = NA) +
    labs(
      y = "No. Grid cells",
      x = "Ploidy",
      subtitle = "d"
    ) +
    standard_theme2() +
    theme(
      legend.position = "none",
      plot.subtitle = element_text(face = "bold")
    )
  
  b + c + d + plot_layout(ncol = 2, nrow = 2)
  
}
