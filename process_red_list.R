library(tidyverse)
library(drake)
library(readxl)

download.file(
  "http://www.rdplants.org/gl/FernGreenListV1.01.xls",
  "data/FernGreenListV1.01.xls")

green_list <- read_excel("data/FernGreenListV1.01.xls") %>%
  select(id = ID20160331, hybrid_variety = `雑種・品種`, 
         scientific_name = `GreenList学名`, family = `PPG科名`,
         endemic = `固有`, cons_status = `RL2012`)

green_list %>%
  filter(!is.na(cons_status))




red_list <- read_csv("data/redlist_sy21.csv", locale=locale(encoding="CP932")) %>%
  select(rank = `ランク`, species = `学名`)

world_ferns <- read_csv("data/world_ferns.csv")

loadd(occ_data_pteridos)

loadd(ppgi)


# Make vector of pteridophyte genera. Note this includes
# many NON pteridophyte synonyms, but they could
# potential be names of pteridophytes once synonymy
# is accounted for.
wf_pterido_genera <- unique(world_ferns$genericName)

# check genera not in ppgi
red_list %>%
  mutate(genus = str_split(species, " ") %>% map_chr(1)) %>%
  anti_join(select(ppgi, genus, family, order)) %>% 
  count(genus, sort = TRUE) %>%
  filter(genus %in% wf_pterido_genera)


# Filter only to species in official list
pterido_taxa_df <- occ_data_pteridos %>%
  mutate(genus = str_split(taxon_name, " ") %>% map_chr(1)) %>%
  mutate(species = str_split(taxon_name, " ") %>% map_chr(2)) %>%
  mutate(genus_species = paste0(genus, " ", species)) %>%
  select(taxon_name, genus_species) %>%
  unique()
  
pterido_taxa <- c(pull(pterido_taxa_df, taxon_name), pull(pterido_taxa_df, genus_species)) %>% unique

red_list %>%
  filter(species %in% pterido_taxa_df)

red_list %>%
  mutate(genus = str_split(species, " ") %>% map_chr(1)) %>%
  filter(genus %in% wf_pterido_genera)

# fuzzyjoin::stringdist_inner_join