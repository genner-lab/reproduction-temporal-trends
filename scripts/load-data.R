#!/usr/bin/env Rscript

# load libs and funs
source(here::here("scripts/funs-stats.R"))

### LOAD DATA

# load hashes
events.hashes <- bind_rows(lapply(list.files(here("meta-fish-pipe/temp/processing"),full.names=TRUE,recursive=TRUE,pattern="events-hashes.csv"),read_csv)) %>% 
    rename(sampleHash=hashLabel)

# load data
edna.all <- read_csv(file=here("meta-fish-pipe/results/fishes-by-sample.csv"))

# add replicate info from hashes
edna.all %<>% left_join(events.hashes,by=c("primerSet","library","eventID","sampleHash")) %>% 
    mutate(replicatePCR="pcr1")

# add sample details to read data
events.df <- read_csv(file=here("assets/events-master.csv")) %>% 
    mutate(eventDate=ymd(paste(year,month,day,sep="-"))) %>% 
    filter(sampleType=="Water") %>%
    distinct(eventID,partnerID,localityID,localitySite,minimumDepthInMeters,eventDate)

# join
edna.all %<>% left_join(events.df,by="eventID")

# filter events without data, recode the date, and remove the locality stuff (will add back later), add taxonRank
edna.filt <- edna.all %>% 
    select(-localityID,-localitySite,-minimumDepthInMeters) %>% 
    mutate(taxonRank="species") %>% 
    rename(species=assignedName,nreads=nReads)

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
