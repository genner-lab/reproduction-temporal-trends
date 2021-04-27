#!/usr/bin/env Rscript

# load libs/funs
source(here::here("scripts/funs-bioinformatics.R"))

# load latest reflib from the meta-fishlinb project
# https://doi.org/10.5281/zenodo.4443447
source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-load-remote.R")
source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-clean.R")
# change reflib
ref.lib <- reflib.orig

###

# add the SEADNA reference reads
sea.refs.df <- read_csv(file=here("assets/tissues-master.csv"))
sea.refs.fa <- read.FASTA(file=here("assets/reference-library.fasta"))

# remove COI and clean names
sea.refs.fa <- sea.refs.fa[grep("12S",names(sea.refs.fa))]
names(sea.refs.fa) <- str_replace(names(sea.refs.fa),"12S\\|","")

# combine
ref.comb <- c(sea.refs.fa,tab2fas(ref.lib,seqcol="nucleotides",namecol="dbid"))

# write out to temp file for hmm
ape::write.FASTA(ref.comb,file=here("temp/reference-library/custom-refs-temp.fasta"))

# subset the markers with hmmer
# need to temporily paste hmm file 'reference-libraries/hmms' into '../hmms'
setwd(here("scripts"))
# no primers
taberlet.fas <- run_hmmer3(dir=here("temp/reference-library"), infile="custom-refs-temp.fasta", prefix="12s.taberlet.noprimers", evalue="10", coords="env")
miya.fas <- run_hmmer3(dir=here("temp/reference-library"), infile="custom-refs-temp.fasta", prefix="12s.miya.noprimers", evalue="10", coords="env")
# with primers (taberlet is inside Miya)
miya.fas.primers <- run_hmmer3(here("temp/reference-library"), infile="custom-refs-temp.fasta", prefix="12s.miya.primers", evalue="10", coords="env")

# convert to char
sea.refs.char <- lapply(as.character(sea.refs.fa),paste,collapse="")
taberlet.char <- lapply(as.character(taberlet.fas),paste,collapse="")
miya.char <- lapply(as.character(miya.fas),paste,collapse="")
miya.char.primers <- lapply(as.character(miya.fas.primers),paste,collapse="")

# add new subseqs to the seadna df along with lengths
sea.refs.df %<>% mutate(nucleotides=as.character(sea.refs.char[match(otherCatalogNumbers,names(sea.refs.char))]), nucleotides=str_replace_all(nucleotides,"NULL","NA"), 
    nucleotidesFrag.12s.taberlet.noprimers=as.character(taberlet.char[match(otherCatalogNumbers,names(taberlet.char))]), nucleotidesFrag.12s.taberlet.noprimers=str_replace_all(nucleotidesFrag.12s.taberlet.noprimers,"NULL","NA"),
    nucleotidesFrag.12s.miya.noprimers=as.character(miya.char[match(otherCatalogNumbers,names(miya.char))]), nucleotidesFrag.12s.miya.noprimers=str_replace_all(nucleotidesFrag.12s.miya.noprimers,"NULL","NA"),
    nucleotidesFrag.12s.miya.primers=as.character(miya.char.primers[match(otherCatalogNumbers,names(miya.char.primers))]), nucleotidesFrag.12s.miya.primers=str_replace_all(nucleotidesFrag.12s.miya.primers,"NULL","NA"),
    length=str_length(nucleotides),
    lengthFrag.12s.taberlet.noprimers=str_length(nucleotidesFrag.12s.taberlet.noprimers),
    lengthFrag.12s.miya.noprimers=str_length(nucleotidesFrag.12s.miya.noprimers),
    lengthFrag.12s.miya.primers=str_length(nucleotidesFrag.12s.miya.primers))

# filter only identified species with DNAs and clean up
sea.refs.df %<>% filter(taxonRank=="species" & nucleotides!="NA") %>% mutate(dbid=otherCatalogNumbers,sciNameValid=paste(genus,specificEpithet),subphylum="Vertebrata",source="SeaDNA")

# add the miya.primers frag to the gb reflib
ref.lib %<>% mutate(nucleotidesFrag.12s.miya.primers=as.character(miya.char.primers[match(dbid,names(miya.char.primers))]), nucleotidesFrag.12s.miya.primers=str_replace_all(nucleotidesFrag.12s.miya.primers,"NULL","NA"),
    lengthFrag.12s.miya.primers=str_length(nucleotidesFrag.12s.miya.primers),kingdom="Animalia",phylum="Chordata")

# merge the DBs
ref.lib <- bind_rows(ref.lib,sea.refs.df)

## write out

# write out df
write_csv(ref.lib,file=here("temp/reference-library/custom-refs.csv"))

# write out fasta for blastn and epa
ape::write.FASTA(tab2fas(df=ref.lib,seqcol="nucleotides",namecol="dbid"),file=here("temp/reference-library/custom-refs.fasta"))

# load refseq
refseq.refs <- read.FASTA(here("temp/reference-library/refseq-annotated.fasta"))
in.refseq.refs <- str_split(names(refseq.refs),"\\.",2,simplify=TRUE)[,1]

# filter and add taxonomy to refs
ref.lib.filt <- ref.lib %>% filter(!is.na(nucleotidesFrag.12s.miya.noprimers)) %>%
    mutate(label=paste0(dbid,";tax=k:",kingdom,",p:",phylum,",c:",class,",o:",order,",f:",family,",g:",genus,",s:",str_replace_all(sciNameValid," ","_"))) %>%
    filter(!gbAccession %in% in.refseq.refs)

# write out annotated 12S fasta for sintax
ape::write.FASTA(tab2fas(ref.lib.filt,seqcol="nucleotidesFrag.12s.miya.noprimers",namecol="label"),file=here("temp/reference-library/custom-refs-annotated.fasta"))


## quickly QC the 12S reference library

# load prefix
prefix <- "12s.taberlet.noprimers"

# load reference library
ref.lib <- read_csv(here("temp/reference-library/custom-refs.csv"),guess_max=99999,col_types=cols())

# get number haplotypes 
ref.lib.haps <- subset_nucs(pref=prefix,df=ref.lib) %>% haps2fas()

# make a tree
ref.lib.haps.tr <- phylogenize(fas=tab2fas(df=ref.lib.haps,seqcol="nucleotidesFrag",namecol="noms"), prefix=prefix, binLoc="~/Software/standard-RAxML/raxmlHPC-AVX",outdir=here("temp/reference-library"))

# plot the tree
plot_trees(tr=ref.lib.haps.tr,df=ref.lib.haps,prefix=prefix,outdir=here("temp/reference-library"))
