# Define analysis plan
plan <- drake_plan (
  
  # Load and process raw data ----
  
  # Load Pteridophyte Phylogeny Group I (PPGI) taxonomy.
  # - original version
  ppgi_raw = read_csv("data/ppgi_taxonomy.csv"),
  
  # - modify slightly for Pteridophytes of Japan
  ppgi = modify_ppgi(ppgi_raw),
  
  # Load Fern Green List, with conservation status for each species.
  green_list = read_excel("data/FernGreenListV1.01.xls") %>%
    select(taxon_id = ID20160331, scientific_name = `GreenList学名`,
           endemic = `固有`, conservation_status = `RL2012`) %>%
    mutate(taxon_id = as.character(taxon_id)),
  
  # Load reproductive mode data, with one row per species.
  repro_data_raw = read_csv(
    "data/ESM1.csv",
    col_types = "cccccnnn"),
  
  repro_data = process_repro_data(repro_data_raw),
  
  # Load occurrence data, with multiple rows per species.
  # Occurrences are presences in a set of 10km2 grid 
  # cells across Japan, not actual occurrence points of specimens.
  occ_data_raw = read_csv(
    "data/ESM2.csv",
    col_types = "cccnnccc"),
  
  # - occurrence data including ferns and lycophytes
  occ_data_pteridos = clean_names(occ_data_raw) %>%
    add_taxonomy(ppgi),
  
  # - occurrence data including ferns only
  occ_data_ferns = 
    occ_data_pteridos %>%
    assert(not_na, class) %>%
    filter(class == "Polypodiopsida"),
  
  # - occurrence data including apomictic ferns only
  occ_data_apos = 
    occ_data_ferns %>%
    left_join(repro_data) %>%
    filter(reproductive_mode == "apomictic"),
  
  # - occurrence data including endangered pteridophytes only
  occ_data_endangered = 
    occ_data_pteridos %>%
    left_join(green_list) %>%
    filter(!is.na(conservation_status)),
  
  # Load phylogenetic tree of all non-hybrid pteridophyte
  # taxa based on rbcL gene from phylogenetic analysis with
  # mrBayes on CIPRES.
  japan_pterido_tree_raw = read_nexus_in_zip(
    "data/japan_pterido_rbcl_cipres.zip", "infile.nex.con.tre")[[2]],
  
  # Process trees.
  # - tree including ferns and lycophtyes
  japan_pterido_tree = format_tip_labels(japan_pterido_tree_raw),
  # - tree including ferns only
  japan_fern_tree = drop.tip(
    japan_pterido_tree, 
    setdiff(japan_pterido_tree$tip.lab, occ_data_ferns$taxon_id)
  ),
  
  # Load basic world map.
  world_map = ggplot2::map_data("world") %>%
    rename(longitude = long, latitude = lat),
  
  # Load list of all 10km2 grid cells across Japan.
  all_cells = read_csv(
    "data/2_grid_cells_all.csv",
    col_types = "ccnn") %>%
    select(secondary_grid_code = id, longitude = x, latitude = y),
  
  # Analyze basic statistics ----
  
  # Make taxonomic data table for all pteridophytes.
  taxonomic_data = summarize_taxonomy(repro_data, occ_data_pteridos),
  
  # Count species per grid cell excluding hybrids.
  species_per_cell = count_species_per_cell(occ_data_pteridos, repro_data),
  
  # Count grid cells per species (CPS).
  cells_per_species = count_cells_per_species(occ_data_pteridos),
  
  # Count CPS by reproductive mode.
  cps_by_repro = count_cells_per_species_by_repro(
    occ_data_pteridos, repro_data
  ),
  
  # Calculate mean CPS by reproductive mode.
  cps_by_repro_means = avg_cells_per_species_by_repro(
    cps_by_repro
  ),
  
  # Count CPS by ploidy level.
  cps_by_ploidy = count_cells_per_species_by_ploidy(
    occ_data_pteridos, repro_data
  ),
  
  # Calculate mean CPS by ploidy level.
  cps_by_ploidy_means = avg_cells_per_species_by_ploidy(
    cps_by_ploidy
  ),
  
  # Run analysis of variance (AOV) on CPS by reproductive mode.
  cps_by_repro_model_summary = aov(
    n_grids ~ reproductive_mode, 
    data = cps_by_repro) %>% tidy,
  
  # Make tibble with latitudinal breadth and
  # reproductive mode for all non-hybrids.
  lat_by_repro = count_lat_by_repro(occ_data_pteridos, repro_data),
  
  # Calculate mean latitudinal breadth per species by reproductive mode.
  lat_by_repro_means = avg_lat_by_repro(lat_by_repro),
  
  # Run analysis of variance (AOV) on 
  # latitudinal breadth by reproductive mode.
  lat_by_repro_model = aov(
    lat_breadth ~ reproductive_mode, 
    data = lat_by_repro),
  
  lat_by_repro_model_summary = tidy(lat_by_repro_model),
  
  # AOV of latidudinal breadth by reproductive mode 
  # showed a significant difference, so
  # run Tukey HSD test on results.
  lat_by_repro_tukey = TukeyHSD(lat_by_repro_model) %>% tidy,
  
  # Analyze species richness ----
  
  # - Ferns and lycophytes
  richness_pteridos = make_richness_matrix(occ_data_pteridos, all_cells),
  
  # - Ferns only
  richness_ferns = make_richness_matrix(occ_data_ferns, all_cells),
  
  # - Apomictic ferns only
  richness_apos = make_richness_matrix(occ_data_apos, all_cells),
  
  # - Endangered ferns only
  richness_endangered = make_richness_matrix(occ_data_endangered, all_cells),
  
  ### Analyze phylogenetic diversity ----
  
  # Make community matrix (presence/absence of each species in
  # 1km2 grid cells), only including species in tree.
  # - Ferns and lycophytes
  comm_pteridos = make_comm_matrix(occ_data_pteridos) %>% 
    match_comm_and_tree(japan_pterido_tree, "comm"),
  
  # - Ferns only
  comm_ferns = make_comm_matrix(occ_data_ferns) %>% 
    match_comm_and_tree(japan_fern_tree, "comm"),
  
  # Calculate PD
  # - Ferns and lycophytes
  pd_pteridos = calc_pd_comm(comm_pteridos, japan_pterido_tree),
  
  # - Ferns only
  pd_ferns = calc_pd_comm(comm_ferns, japan_fern_tree),
  
  # Calculate percentage of apomictic ferns ----
  # percentage is out of all species that are in both
  # occurrence data and reproductive data, so excludes hybrids.
  # Set 0s to NA for plotting on log-scale
  percent_apomictic_ferns = calc_percent_apomictic(
    occ_data_ferns,
    repro_data,
    all_cells
  ),
  
  # Combine richness and PD metrics for plotting
  alpha_diversity_pteridos = left_join(
    richness_pteridos, pd_pteridos
  ),
  
  alpha_diversity_ferns = left_join(
    richness_ferns, pd_ferns
  ),
  
  # Make figures ----
  
  # - Richness of ferns and lycophytes
  pterido_richness_map = make_diversity_map(
    div_data = alpha_diversity_pteridos, 
    world_map = world_map, 
    occ_data = occ_data_pteridos, 
    div_metric = "richness", 
    metric_title = "No. species ",
    label = "a") +
    scale_fill_scico(palette = "bamako", na.value="grey"),
  
  # - PD of ferns and lycophytes
  pterido_pd_map = make_diversity_map(
    div_data = alpha_diversity_pteridos, 
    world_map = world_map, 
    occ_data = occ_data_pteridos, 
    div_metric = "pd_obs", 
    metric_title = "Phylo. \n diversity",
    label = "b") +
    scale_fill_scico(palette = "bamako", na.value="grey"),

  # - PD of ferns only
  fern_pd_map = make_diversity_map(
    div_data = alpha_diversity_ferns, 
    world_map = world_map, 
    occ_data = occ_data_pteridos, 
    div_metric = "pd_obs", 
    metric_title = "Phylo. \n diversity",
    label = "c") +
    scale_fill_scico(palette = "bamako", na.value="grey"),
  
  # - Richness of apomictic ferns
  fern_apo_map = make_diversity_map(
    div_data = richness_apos, 
    world_map = world_map, 
    occ_data = occ_data_pteridos, 
    div_metric = "richness", 
    metric_title = "No. species",
    label = "d") +
    scale_fill_scico(palette = "bamako", na.value="grey"),
  
  # - Percent apomictic taxa (on log-scale)
  fern_apo_frac_map = make_diversity_map(
    # Replace 0 with NA for log-transform
    div_data = percent_apomictic_ferns %>%
      mutate(percent_apomictic = case_when(
        percent_apomictic == 0 ~ NaN,
        TRUE ~ percent_apomictic * 100
      )), 
    world_map = world_map, 
    occ_data = occ_data_pteridos, 
    div_metric = "percent_apomictic", 
    metric_title = "% apomictic",
    label = "e") +
    scale_fill_scico(
      name = "% apomictic",
      palette = "bamako", na.value="#003F4C", # set NA value to be same color as 0
      trans = "log", breaks = c(6,12,25,50,100)),

  # - Richness of endangered ferns and lycophytes (on log-scale)
  fern_endangered_map = make_diversity_map(
    # Replace 0 with NA for log-transform
    div_data = richness_endangered %>%
      mutate(richness = case_when(
        richness == 0 ~ NaN,
        TRUE ~ richness
      )), 
    world_map = world_map, 
    occ_data = occ_data_pteridos, 
    div_metric = "richness", 
    metric_title = "No. species",
    label = "f") +
    scale_fill_scico(
      palette = "bamako", na.value="#003F4C", # set NA value to be same color as 0
      trans = "log", breaks = c(3,6,12,24,48)),
  
  # Combine subplots.
  fig_1 = patchwork::wrap_plots(
    list(
      pterido_richness_map, pterido_pd_map,
      fern_pd_map, fern_apo_map,
      fern_apo_frac_map, fern_endangered_map),
    ncol = 2, nrow = 3
  ),
  
  # Write out final plot at full-page size.
  fig_1_out = ggsave(
    plot = fig_1, 
    filename = file_out(here("manuscript/fig_1.pdf")),
    height = 234, width = 174, units = "mm"),

  # Write out manuscript ----
  ms = rmarkdown::render(
    knitr_in(here::here("manuscript/japan_ferns_diversity_ms.Rmd")),
    output_file = file_out(here::here("manuscript/japan_ferns_diversity_ms.docx")),
    quiet = TRUE)
)
