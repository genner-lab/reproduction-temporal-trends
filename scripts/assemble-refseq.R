#!/usr/bin/env Rscript

# load libs/funs
source(here::here("scripts/funs-bioinformatics.R"))

# read in the mito refseq catalogue
# first download with `taxonomic-assignment.sh`
mito.cat <- read_tsv(here("temp/reference-library/refseq.mitochondrion.cat.tsv"),guess_max=999999,col_names=c("taxid","scientificName","accession","dir","status","len"))

# read in the mito refseq fasta
mt <- read.FASTA(here("temp/reference-library/refseq.mitochondrion.genomic.fna"))

# remove the long and short sequences leaving animal mitogenomes 
mt <- mt[lapply(mt,length) < 21000 & lapply(mt,length) > 12000]

# get accesssions from names
accs <- names(mt)

# filter the catalogue down
mito.cat %<>% filter(accession %in% accs)

# add the genus from sciName
mito.cat %<>% mutate(genus=str_split_fixed(scientificName," ",2)[,1])

# create a db from GBIF
td_create(provider="gbif",dbdir=here("temp/reference-library"),overwrite=TRUE)

# connect to db
td_connect(dbdir=here("temp/reference-library"))

# pull out all the db then subset and remove dups
gbif.taxonomy <- taxa_tbl(provider="gbif") %>% 
    collect() %>% 
    filter(scientificName %in% pull(mito.cat,genus)) %>% 
    select(kingdom,phylum,class,order,family,genus) %>% 
    distinct() %>% 
    add_count(genus) %>% 
    filter(n==1) %>% 
    arrange(kingdom,phylum,class,order,family,genus) %>% 
    mutate(classified=TRUE)

# annotate the mito.cat
mito.cat.annotated <- mito.cat %>% 
    left_join(gbif.taxonomy,by="genus") %>% 
    filter(classified==TRUE) %>% 
    filter(kingdom=="Animalia") %>% 
    arrange(kingdom,phylum,class,order,family,genus,scientificName) %>% 
    mutate(label=paste0(accession,";tax=k:",kingdom,",p:",phylum,",c:",class,",o:",order,",f:",family,",g:",genus,",s:",scientificName)) %>% 
    mutate(label=str_replace_all(label," ","_"))

# drop the unwanted seqs
mt.sub <- mt[names(mt) %in% pull(mito.cat.annotated,accession)]

# rename with taxonomy
names(mt.sub) <- pull(mito.cat.annotated,label)[match(names(mt.sub),pull(mito.cat.annotated,accession))]

# write out
write.FASTA(mt.sub,file=here("temp/reference-library/refseq-annotated.fasta"))

# disconnect from db
td_disconnect()
