#!/usr/bin/env Rscript

# load libs funs
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

# set primer length
trlens <- c(18,20)# tele02

# load up the data
plates <- read_csv(file=here("assets/sequencing-master.csv"))

# check samples are all present in events-master - missing ones are blanks/ntc
events.master <- read_csv(file=here("assets/events-master.csv"))
plates %>% filter(!eventID %in% pull(events.master,eventID)) %>% print(n=Inf)

# filter by marker
plates %<>% filter(primerSet==primer.code & library==lib)

# create the barcode tags
plates %<>% mutate(barcodesFwd=str_replace_all(oligoFwd,"N",""), 
    barcodesFwd=str_trunc(barcodesFwd, width=10, side="right", ellipsis=""),
    barcodesRev=str_replace_all(oligoRev,"N",""),
    barcodesRev=str_trunc(barcodesRev, width=10, side="right", ellipsis=""),
    primerFwd=str_trunc(oligoFwd, width=trlens[1], side="left", ellipsis=""),
    primerRev=str_trunc(oligoRev, width=trlens[2], side="left", ellipsis=""),
    labelFwd=str_trunc(str_replace_all(oligoFwd,"N",""), width=unique(str_length(str_replace_all(oligoFwd,"N",""))-trlens[1]), side="right", ellipsis=""),
    labelRev=str_trunc(str_replace_all(oligoRev,"N",""), width=unique(str_length(str_replace_all(oligoRev,"N",""))-trlens[2]), side="right", ellipsis="")
    )

# check oligos match (should be true)
table(pull(plates,labelFwd) == pull(plates,labelRev))
# cat the labels
plates %<>% mutate(senseLabel=paste(eventID,primerSet,library,replicateFilter,replicatePCR,labelFwd,sep="."), antisenseLabel=paste(eventID,primerSet,library,replicateFilter,replicatePCR,labelRev,sep="."))
# check labels match (should be true)
table(pull(plates,senseLabel) == pull(plates,antisenseLabel))

# create the hashes
plates %<>% mutate(senseLabelMD5=str_trunc(md5(senseLabel),width=12,side="right",ellipsis=""), antisenseLabelMD5=str_trunc(md5(antisenseLabel),width=12,side="right",ellipsis="")) 
# check md5s match (should be true)
table(pull(plates,senseLabelMD5) == pull(plates,antisenseLabelMD5))
# check for collisions (should be true)
length(pull(plates,senseLabelMD5)) == length(unique(pull(plates,senseLabelMD5)))

# rename as hashLabel
plates %<>% mutate(hashLabel=if_else(senseLabelMD5 == antisenseLabelMD5, as.character(senseLabelMD5), "NA"))

# reverse comp the barcodes
plates %<>% mutate(barcodesFwdRevComp=revcomp_dna(dnacol=barcodesFwd),barcodesRevRevComp=revcomp_dna(dnacol=barcodesRev))


# RUN TO HERE IF JUST WANT THE 'PLATES' DATAFRAME (STEP REQUIRED FOR TAXON ASSIGNMENT)

# make fasta 
bcF <- tab2fas(df=plates, seqcol="barcodesFwd", namecol="hashLabel")
bcR <- tab2fas(df=plates, seqcol="barcodesRev", namecol="hashLabel")

# create dirs if needed
dir.create(path=proj.path,recursive=TRUE)

# write out in folder for running analyses
write.FASTA(bcF, file=here(proj.path,"barcodes-sense.fas"))
write.FASTA(bcR, file=here(proj.path,"barcodes-antisense.fas"))
