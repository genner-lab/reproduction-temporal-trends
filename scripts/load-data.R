#!/usr/bin/env Rscript

# load libs and funs
source(here::here("scripts/funs-stats.R"))

### LOAD DATA

# load data
edna.all <- read_csv(file=here("temp/results/results-by-marker.csv"))

# filter events without data, recode the date, and remove the locality stuff (will add back later), add taxonRank
edna.filt <- edna.all %>% filter(!grepl("Blank",replicateFilter)) %>% mutate(eventDate=ymd(eventDate)) %>% select(-localityID,-localitySite,-minimumDepthInMeters) %>% mutate(taxonRank="species")

# load primer bias
primer.bias <- read_csv(file=here("temp/results/primer-efficiency-results.csv"))

# load uk.species to get taxonomy
uk.species <- read_csv("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/assets/species-table.csv",guess_max=99999)
tax.table <- uk.species %>% rename(sciNameValid="validName") %>% select(class,order,family,genus,sciNameValid) %>% distinct()

# load events
events <- read_csv(file=here("assets/events-master.csv"))

# load the trawl data
trad.master <- read_csv(file=here("assets/trad-data-master.csv")) %>% 
    rename(species=scientificName) %>% 
    mutate(partnerID=str_split_fixed(eventID,"-",3)[,1]) %>% # add partnerID
    mutate(fieldNumber=if_else(is.na(fieldNumber),eventID,fieldNumber)) # add a dummy field number for those missing
