#!/usr/bin/env Rscript

# load libs and funs
source(here::here("scripts/funs-bioinformatics.R"))

# set the primer
primer.code <- "tele02"

# set the name of the folder that the file will be added to 
lib <- "lib1"
lib <- "lib2"
lib <- "lib3"
lib <- "lib4"

# make paths
proj.path <- here("temp/processing",paste0(primer.code,"-",lib))
# confirm proj path 
print(proj.path)


#################################################### 
#################################################### 

### STEPS 1-3: RUN `scripts/taxonomic-assignment.sh` (HOMOLOGY SEARCH AND SINTAX) ### 
### STEPS 1-3: RUN `scripts/taxonomic-assignment.sh` (HOMOLOGY SEARCH AND SINTAX) ### 



#################################################### 
#################################################### 

### STEP 4: FOR EACH LIB SUMMARISE THE CLEANED ASVS ### 
### STEP 4: FOR EACH LIB SUMMARISE THE CLEANED ASVS ### 

# read asv table back in written table to get format correct
otu.tab <- as_tibble(read.table(file=paste0(proj.path,"/results/asv-table.tsv"), sep="\t", header=TRUE, as.is=TRUE, row.names=1, check.names=FALSE), rownames="asv")
dnas.curated <- read.FASTA(file=paste0(proj.path,"/results/asvs.fna"))

# load up data cleaned by the HMM
dnas.curated.clean <- read.FASTA(file=paste0(proj.path,"/results/asvs-clean.fasta"))
dnas.curated.dirty <- dnas.curated[which(!(names(dnas.curated) %in% names(dnas.curated.clean)))]
#write.FASTA(dnas.curated.dirty, file=paste0(proj.path,"/results/otus-dirty.fasta"))
# subset the clean sequences from the ASV table
otu.tab.clean <- otu.tab %>% filter(asv %in% names(dnas.curated.clean))
# write out
write_tsv(otu.tab.clean,file=paste0(proj.path,"/results/asv-table-clean.tsv"))

# get numbers of seqs lost during hmm search
total <- otu.tab %>% summarise_if(is.numeric, sum, na.rm=TRUE) %>% rowSums()
lost <- otu.tab %>% filter(!asv %in% names(dnas.curated.clean)) %>% summarise_if(is.numeric, sum, na.rm=TRUE) %>% rowSums()
retained <- otu.tab.clean %>% summarise_if(is.numeric, sum, na.rm=TRUE) %>% rowSums()
lost+retained == total
print(retained)



#################################################### 
#################################################### 

### STEP 5: FOR EACH LIB SUMMARISE SINTAX RESULTS AND ISOLATE FISH READS ### 
### STEP 5: FOR EACH LIB SUMMARISE SINTAX RESULTS AND ISOLATE FISH READS ### 

# read in taxonomy assignment
tax.ass.df <- read_tsv(file=paste0(proj.path,"/results/sintax-output.tsv"),col_names=c("asv","idsProbs","strand","ids"),guess_max=999999)

tax.ass.df %<>% 
    mutate(sintaxBS=str_replace_all(map(str_split(pull(tax.ass.df,idsProbs),"\\("), last),"\\)","")) %>%
    mutate(idsProbsTemp=str_replace_all(idsProbs,"[a-z]:",""),idsProbsTemp=str_replace_all(idsProbsTemp,"\\([0-9].[0-9][0-9]\\)","")) %>%
    separate(idsProbsTemp,into=c("kingdom","phylum","class","order","family","genus","sintaxSpeciesID"),sep=",") %>%
    mutate(sintaxSpeciesID=str_replace_all(sintaxSpeciesID,"_", " ")) %>%
    mutate(isFish=if_else(class=="Cephalaspidomorphi" | class=="Elasmobranchii" | class=="Actinopterygii", TRUE, FALSE)) %>%
    select(asv,isFish,kingdom,phylum,class,order,family,genus,sintaxSpeciesID,sintaxBS)

# remake the otus frame with DNAs
otus.df <- tibble(asv=names(dnas.curated.clean), dnas=sapply(as.character(dnas.curated.clean), paste, collapse="")) %>% mutate(asvLength=str_length(dnas))

# to make the annotated ASV table
# add read numbers to OTU table
# get number of reads per sample
otu.tab.clean %>% select_if(is.numeric) %>% summarise_all(list(sum))
# number reads per OTU
otu.tab.clean.sums <- otu.tab.clean %>% select_if(is.numeric) %>% mutate(sums=rowSums(.))
otu.tab.clean.sums %>% select(sums)

# make new table 
otu.reads <- tibble(asv=pull(otu.tab.clean,asv), totalReads=pull(otu.tab.clean.sums,sums))
otu.reads %<>% mutate(propTotal=totalReads/sum(totalReads))

# merge with sintax results and reads tables
sintax.results <- purrr::reduce(list(otus.df,tax.ass.df,otu.reads),left_join, by="asv")

# rearrange
sintax.results %<>% 
    select(asv,kingdom,phylum,class,order,family,genus,sintaxSpeciesID,asvLength,sintaxBS,totalReads,propTotal,dnas) %>% 
    mutate(lib=paste0(primer.code,"-",lib)) %>% 
    arrange(desc(totalReads))

# write out original data
write_tsv(sintax.results, path=paste0(proj.path,"/results/sintax-results-sorted.tsv"))

# prep for blast search - generate fish queries file
# import reference library and exclusions file
custom.refs <- read_csv(here("temp/reference-library/custom-refs.csv"),guess_max=99999,col_types=cols())

# reload
sintax.results <- read_tsv(file=paste0(proj.path,"/results/sintax-results-sorted.tsv"))

# filter out fish reads
sintax.results %<>% filter(class=="Actinopterygii" | class=="Cephalaspidomorphi" | class=="Elasmobranchii")
print(sintax.results,n=Inf)

# write out
write.FASTA(tab2fas(df=sintax.results, seqcol="dnas", namecol="asv"), file=paste0(proj.path,"/results/fishqueries.fasta"))



#################################################### 
#################################################### 

### STEP 6: RUN `scripts/taxonomic-assignment.sh` (BLAST ASSIGNMENT) ### 
### STEP 6: RUN `scripts/taxonomic-assignment.sh` (BLAST ASSIGNMENT) ### 



#################################################### 
#################################################### 

### STEP 7: PROCESS BLAST RESULTS ### 
### STEP 7: PROCESS BLAST RESULTS ### 

# load blast result
local.db.blast <- read_tsv(paste0(proj.path,"/results/fish-blast-result.tsv"))

# choose "best" hit based on bitscore
# also add scinames
local.db.blast.sorted <- local.db.blast %>%
    group_by(asv) %>%
    arrange(desc(blastBitscore),.by_group=TRUE) %>%
    filter(blastBitscore==max(blastBitscore)) %>%
    mutate(blastSpeciesID=pull(custom.refs,sciNameValid)[match(blastDbid,pull(custom.refs,dbid))]) %>%
    arrange(blastSpeciesID,.by_group=TRUE) %>%
    mutate(blastSpeciesID=paste(unique(blastSpeciesID),collapse="; ")) %>%
    slice(1) %>% 
    ungroup()
local.db.blast.sorted %>% arrange(desc(blastBitscore)) %>% print(n=Inf)

# write out for later 
local.db.blast.sorted %>% write_csv(paste0(proj.path,"/results/fish-blast-result-sorted.csv"))



#################################################### 
#################################################### 

### STEP 8: PREP INPUT FOR EPA ASSIGNMENT ###
### STEP 8: PREP INPUT FOR EPA ASSIGNMENT ###

# prefixes for hmm
prefix <- "12s.taberlet.noprimers"

# subset the marker from the custom reference library reflib - using 50% length cutoff
reflib.sub <- subset_by_marker(prefix=prefix,df=custom.refs,thresh=0.5)

# collapse haps
reflib.red <- reflib.sub %>% group_by(sciNameValid) %>% group_modify(~ hap_collapse_df(df=.x,lengthcol=paste0("lengthFrag.",prefix),nuccol=paste0("nucleotidesFrag.",prefix),cores=1)) %>% ungroup()

# load up the queries
fish.queries <- read.FASTA(file=paste0(proj.path,"/results/fishqueries.fasta"))

# concatenate with the reference and align
combined.fas <- c(fish.queries,tab2fas(reflib.red, seqcol=paste0("nucleotidesFrag.",prefix),namecol="dbid"))

# remove Ns
combined.fas <- rm_ns(bin=combined.fas)

# align
combined.fas.aligned <- ips::mafft(combined.fas,exec="mafft",method="retree 2",maxiterate=2)

# make a subdir
dir.create(paste0(proj.path,"/results/epa"))

# make a reference and query alignment
references.aligned <- combined.fas.aligned[grep("asv[0-9]+",rownames(combined.fas.aligned),invert=TRUE),]
queries.aligned <- combined.fas.aligned[grep("asv[0-9]+",rownames(combined.fas.aligned)),]

write.FASTA(references.aligned, file=paste0(proj.path,"/results/epa/references.aligned.fasta"))
write.FASTA(queries.aligned, file=paste0(proj.path,"/results/epa/queries.aligned.fasta"))

# generate taxonomy file
reflib.red %>% mutate(sciNameValid=str_replace_all(sciNameValid," ", "_")) %>% 
    mutate(taxonomy=paste(subphylum,class,order,family,genus,sciNameValid,sep=";")) %>% 
    select(dbid,taxonomy) %>% 
    write_tsv(file=paste0(proj.path,"/results/epa/references.taxonomy.tsv"),col_names=FALSE)



#################################################### 
#################################################### 

### STEP 9: RUN `scripts/taxonomic-assignment.sh` (EPA ASSIGNMENT) ###
### STEP 9: RUN `scripts/taxonomic-assignment.sh` (EPA ASSIGNMENT) ### 



#################################################### 
#################################################### 

### STEP 10: PROCESS EPA RESULTS ###
### STEP 10: PROCESS EPA RESULTS ###

# load up epa results
epa.results <- read_tsv(file=paste0(proj.path,"/results/epa/assign_per_query.tsv"))
epa.results %<>% rename(asv=name)
# LWR: likelihood weight that was assigned to this exact taxonomic path
# fract: LWR divided by the global total likelihood weight
# aLWR: accumulated likelihood weights that were assigned either to this taxonomic path or any taxonomic path below this
# afract: aLWR divided by the global total likelihood weight
# taxopath: the taxonomic path

# best the best id
epa.results.best <- epa.results %>% 
    group_by(asv) %>%
    arrange(desc(LWR),.by_group=TRUE) %>% 
    slice(1) %>% 
    mutate(epaID=sapply(str_split(taxopath,";"),last), epaID=str_replace_all(epaID,"_"," ")) %>% 
    ungroup()

# get species IDs
epa.results.species <- epa.results %>% 
    group_by(asv) %>% 
    filter(grepl("_",taxopath)) %>% 
    arrange(taxopath,.by_group=TRUE) %>% 
    mutate(epaBestSppID=sapply(str_split(taxopath,";"),last), epaBestSppID=str_replace_all(epaBestSppID,"_"," ")) %>%  
    mutate(epaAllSpp=paste(unique(epaBestSppID),collapse="; ")) %>% 
    arrange(desc(LWR),.by_group=TRUE) %>% 
    slice(1) %>% 
    ungroup() %>%
    select(asv,epaBestSppID,epaAllSpp)

# join and tidy
epa.results.filtered <- left_join(epa.results.best,epa.results.species) %>% 
    select(asv,LWR,epaID,epaBestSppID,epaAllSpp) %>% 
    rename(epaIdLWR=LWR)

# write out 
epa.results.filtered %>% write_csv(file=paste0(proj.path,"/results/epa/epa-results-filtered.csv"))



#################################################### 
#################################################### 

### STEP 11: COMBINE ASSIGNMENT RESULTS ###
### STEP 11: COMBINE ASSIGNMENT RESULTS ###

# reload global blast, local blast, epa results
fish.otus <- read_tsv(file=paste0(proj.path,"/results/sintax-results-sorted.tsv")) %>% filter(class=="Actinopterygii" | class=="Cephalaspidomorphi" | class=="Elasmobranchii") %>% select(-lib)
local.db.blast.sorted <- read_csv(paste0(proj.path,"/results/fish-blast-result-sorted.csv"))
epa.results.filtered <- read_csv(file=paste0(proj.path,"/results/epa/epa-results-filtered.csv"))
# read in contamination exclusions
contam.table <- read_csv(here("assets/contaminants-exclude.csv"))

# combine
taxonomy.results <- purrr::reduce(list(fish.otus,local.db.blast.sorted,epa.results.filtered),left_join, by="asv")

# get mode
minLen <- ceiling(mode_avg(pull(taxonomy.results,asvLength))*0.9)
# plot(sort(taxonomy.results$bitscoreLocal))
# plot(sort(taxonomy.results$pidentLocal))

# process
# two step
# either (a) or (b) must be satisfied
# (a) epaAssign = highest likelihood EPA id and blast id the same, EPA likelihood > 0.9, and pident > 0.9
# (b) blastAssign = the best species-level EPA id (not necesarily best overall id) and blast id the same, and pident >0.97 and lengthLocal >0.9 of modal ASV length 
taxonomy.results %<>% 
    mutate(epaAssign=if_else(epaID==blastSpeciesID & epaIdLWR>=0.9 & blastPident>=90,TRUE,FALSE)) %>% 
    mutate(blastAssign=if_else(str_detect(blastSpeciesID,epaBestSppID) & blastLength>=minLen & blastPident>=97,TRUE,FALSE)) %>% 
    mutate(assigned=if_else(epaAssign==TRUE | blastAssign==TRUE,TRUE,FALSE)) %>% 
    mutate(assignedName=if_else(assigned==TRUE,if_else(epaAssign==TRUE,epaID,epaBestSppID),"NA")) %>% # assigns name based on either assignment 
    mutate(lib=paste0(primer.code,"-",lib))

# filter contams and make a label
contam.table %<>% filter(primer==primer.code) %>% mutate(label=paste(otuCode,bestId,sep="|"))

# add contam status
taxonomy.results %<>% mutate(contaminationID=purrr::map_chr(dnas,find_matches,referencedb=contam.table), isContaminated=if_else(is.na(contaminationID),FALSE,TRUE)) 

# check
taxonomy.results %>% select(asv,sintaxSpeciesID,assigned,assignedName,contaminationID) %>% arrange(desc(assigned),asv) %>% print(n=Inf)

# take a look
taxonomy.results %>% 
    select(-kingdom,-phylum,-class,-order,-family,-genus,-propTotal,-dnas,-blastDbid,-blastEvalue,-blastScore,-blastBitscore,-lib) %>%
    arrange(desc(assigned),desc(epaAssign),desc(blastAssign)) %>% 
    write_csv(file=paste0(proj.path,"/results/sanity-check.csv"))

# write out
taxonomy.results %>% write_csv(file=paste0(proj.path,"/results/taxonomic-assignment.csv"))



#################################################### 
#################################################### 

### STEP 12: GET SAMPLE LEVEL ASSIGNMENTS ###
### STEP 12: GET SAMPLE LEVEL ASSIGNMENTS ###

# first run `prep-barcodes.R` to get the "plates" and "events.master" objects
# first run `prep-barcodes.R` to get the "plates" and "events.master" objects

# load results
taxonomy.results <- read_csv(file=paste0(proj.path,"/results/taxonomic-assignment.csv"),guess_max=99999)
otu.tab <- read.table(file=paste0(proj.path,"/results/asv-table.tsv"), sep="\t", header=TRUE, as.is=TRUE, row.names=1, check.names=FALSE)
curated.table.fish <- as_tibble(otu.tab, rownames="asv")

# check for dups  (should be nothing)
(plates %>% filter(duplicated(hashLabel)) %>% dim())[1]
# md5 not matching (should be nothing)
(curated.table.fish %>% pivot_longer(cols=-asv,names_to="sampleHash",values_to="nreads") %>% filter(!sampleHash %in% pull(plates,hashLabel)) %>% dim())[1]
# events not matching (should be negs)
plates %>% filter(!eventID %in% pull(events.master,eventID)) %>% select(eventID,hashLabel,primerSet,library,replicateFilter,replicatePCR) %>% print(n=Inf)

# filter, and replace the names in the otu table
taxonomy.results.assigned <- taxonomy.results %>% filter(assigned==TRUE & isContaminated==FALSE)
# remove the unassigned reads
curated.table.fish %<>% filter(asv %in% pull(taxonomy.results.assigned,asv)) %>% mutate(asv=pull(taxonomy.results.assigned,assignedName)[match(asv,pull(taxonomy.results.assigned,asv))])

# collapse by species
curated.table.fish %<>% group_by(asv) %>% summarise_all(list(sum))

# process the table 
curated.fish.col <- curated.table.fish %>% 
    # reshape into column format
    pivot_longer(cols=-asv,names_to="sampleHash",values_to="nreads") %>%
    # remove zero reads and rename species column
    filter(nreads > 0) %>% 
    rename(species=asv) %>% 
    # add eventID
    mutate(eventID=pull(plates,eventID)[match(sampleHash,pull(plates,hashLabel))]) %>% 
    # add PCR rep and filtration rep 
    mutate(replicateFilter=pull(plates,replicateFilter)[match(sampleHash,pull(plates,hashLabel))], replicatePCR=pull(plates,replicatePCR)[match(sampleHash,pull(plates,hashLabel))]) %>% 
    # add marker and lib
    mutate(primerSet=pull(plates,primerSet)[match(sampleHash,pull(plates,hashLabel))], library=pull(plates,library)[match(sampleHash,pull(plates,hashLabel))]) %>%
    # add partner
    mutate(partnerID=pull(events.master,partnerID)[match(eventID,pull(events.master,eventID))]) %>% 
    # add locality
    mutate(localityID=pull(events.master,localityID)[match(eventID,pull(events.master,eventID))]) %>% 
    # add site
    mutate(localitySite=pull(events.master,localitySite)[match(eventID,pull(events.master,eventID))]) %>% 
    # add the date
    mutate(eventDate=paste(pull(events.master,year)[match(eventID,pull(events.master,eventID))],pull(events.master,month)[match(eventID,pull(events.master,eventID))],pull(events.master,day)[match(eventID,pull(events.master,eventID))],sep="-")) %>% 
    # depth
    mutate(minimumDepthInMeters=pull(events.master,minimumDepthInMeters)[match(eventID,pull(events.master,eventID))]) %>% 
    # select cols and sort
    select(primerSet,library,partnerID,eventID,sampleHash,replicateFilter,replicatePCR,localityID,localitySite,minimumDepthInMeters,eventDate,species,nreads) %>%
    arrange(partnerID,eventDate,localityID,localitySite,replicateFilter,replicatePCR,desc(nreads))

# write out the final table
curated.fish.col %>% write_csv(file=paste0(proj.path,"/results/fish-final-table.csv"))


# UNASSIGNED READS
# load sintax
sintax.results <- read_tsv(file=paste0(proj.path,"/results/sintax-results-sorted.tsv"))
sintax.results.unassigned <- sintax.results %>% filter(!asv %in% pull(taxonomy.results.assigned,asv))

# add the contamination data
sintax.results.unassigned %<>% mutate(contaminationID=purrr::map_chr(dnas,find_matches,referencedb=contam.table), isContaminated=if_else(is.na(contaminationID),FALSE,TRUE)) 

# reload asv table
curated.table.fish <- as_tibble(otu.tab, rownames="asv")
# subset 
unassigned.otu.table <- curated.table.fish %>% filter(asv %in% pull(sintax.results.unassigned,asv))

# add data
unassigned.otu.df <- unassigned.otu.table %>% 
    pivot_longer(cols=-asv,names_to="sampleHash",values_to="nreads") %>%
    filter(nreads > 0) %>% 
    mutate(eventID=pull(plates,eventID)[match(sampleHash,pull(plates,hashLabel))]) %>% 
    mutate(replicateFilter=pull(plates,replicateFilter)[match(sampleHash,pull(plates,hashLabel))], replicatePCR=pull(plates,replicatePCR)[match(sampleHash,pull(plates,hashLabel))]) %>% 
    mutate(primerSet=pull(plates,primerSet)[match(sampleHash,pull(plates,hashLabel))], library=pull(plates,library)[match(sampleHash,pull(plates,hashLabel))]) %>%
    mutate(partnerID=pull(events.master,partnerID)[match(eventID,pull(events.master,eventID))]) %>% 
    mutate(localityID=pull(events.master,localityID)[match(eventID,pull(events.master,eventID))]) %>% 
    mutate(localitySite=pull(events.master,localitySite)[match(eventID,pull(events.master,eventID))]) %>% 
    mutate(eventDate=paste(pull(events.master,year)[match(eventID,pull(events.master,eventID))],pull(events.master,month)[match(eventID,pull(events.master,eventID))],pull(events.master,day)[match(eventID,pull(events.master,eventID))],sep="-")) %>% 
    mutate(minimumDepthInMeters=pull(events.master,minimumDepthInMeters)[match(eventID,pull(events.master,eventID))])  %>% 
    mutate(phylum=pull(sintax.results.unassigned,phylum)[match(asv,pull(sintax.results.unassigned,asv))],
        class=pull(sintax.results.unassigned,class)[match(asv,pull(sintax.results.unassigned,asv))],
        order=pull(sintax.results.unassigned,order)[match(asv,pull(sintax.results.unassigned,asv))],
        family=pull(sintax.results.unassigned,family)[match(asv,pull(sintax.results.unassigned,asv))],
        genus=pull(sintax.results.unassigned,genus)[match(asv,pull(sintax.results.unassigned,asv))],
        species=pull(sintax.results.unassigned,sintaxSpeciesID)[match(asv,pull(sintax.results.unassigned,asv))],
        length=pull(sintax.results.unassigned,asvLength)[match(asv,pull(sintax.results.unassigned,asv))],
        sintaxBS=pull(sintax.results.unassigned,sintaxBS)[match(asv,pull(sintax.results.unassigned,asv))],
        contaminationID=pull(sintax.results.unassigned,contaminationID)[match(asv,pull(sintax.results.unassigned,asv))],
        isContaminated=pull(sintax.results.unassigned,isContaminated)[match(asv,pull(sintax.results.unassigned,asv))]) %>%
    select(primerSet,library,partnerID,eventID,sampleHash,replicateFilter,replicatePCR,localityID,localitySite,minimumDepthInMeters,eventDate,phylum,class,order,family,genus,species,asv,nreads,length,sintaxBS,isContaminated,contaminationID) %>%
    arrange(partnerID,eventDate,localityID,localitySite,replicateFilter,replicatePCR,desc(nreads))
# writeout
unassigned.otu.df %>% write_csv(file=paste0(proj.path,"/results/unassigned-final-table.csv"))



#################################################### 
#################################################### 

### STEP 13: COMBINE LIBRARIES (RUN PREVIOUS STEPS FOR ALL LIBS FIRST) ###
### STEP 13: COMBINE LIBRARIES (RUN PREVIOUS STEPS FOR ALL LIBS FIRST) ###

# taxon assignments
bind_rows(read_csv(file=here("temp/processing/tele02-lib1/results/taxonomic-assignment.csv"),guess_max=999999),
    read_csv(file=here("temp/processing/tele02-lib2/results/taxonomic-assignment.csv"),guess_max=999999),
    read_csv(file=here("temp/processing/tele02-lib3/results/taxonomic-assignment.csv"),guess_max=999999),
    read_csv(file=here("temp/processing/tele02-lib4/results/taxonomic-assignment.csv"),guess_max=999999)) %>%
    tidyr::separate(col=lib,into=c("primerSet","library"),sep="-lib",extra="merge") %>% mutate(library=paste0("lib",library)) %>% # next time add this at the assembly stage and the select below
    select(primerSet,library,asv,totalReads,propTotal,assignedName,assigned,contaminationID,kingdom:sintaxSpeciesID,sintaxBS,asvLength,blastDbid:isContaminated,dnas) %>% 
    write_csv(file=here("temp/results/fish-assignments.csv"))

# otu tables
bind_rows(read_csv(file=here("temp/processing/tele02-lib1/results/fish-final-table.csv"),guess_max=999999),
    read_csv(file=here("temp/processing/tele02-lib2/results/fish-final-table.csv"),guess_max=999999),
    read_csv(file=here("temp/processing/tele02-lib3/results/fish-final-table.csv"),guess_max=999999),
    read_csv(file=here("temp/processing/tele02-lib4/results/fish-final-table.csv"),guess_max=999999)) %>%
    arrange(primerSet,library,partnerID,eventID,replicateFilter,replicatePCR) %>%
    write_csv(file=here("temp/results/results-by-marker.csv"))

# unnassigned table
bind_rows(read_csv(file=here("temp/processing/tele02-lib1/results/unassigned-final-table.csv"),guess_max=999999),
    read_csv(file=here("temp/processing/tele02-lib2/results/unassigned-final-table.csv"),guess_max=999999),
    read_csv(file=here("temp/processing/tele02-lib3/results/unassigned-final-table.csv"),guess_max=999999),
    read_csv(file=here("temp/processing/tele02-lib4/results/unassigned-final-table.csv"),guess_max=999999)) %>%
    arrange(primerSet,library,partnerID,eventID,replicateFilter,replicatePCR,desc(nreads)) %>%
    write_csv(file=here("temp/results/unassigned-results-all.csv"))
