#!/usr/bin/env Rscript
# script to generate primer efficiency values using the DECIPHER 

# load libs
source(here::here("scripts/funs-bioinformatics.R"))

# make results dir
if(!dir.exists(here("temp/results"))){
    dir.create(here("temp/results"),recursive=TRUE)
}

### add primers to custom reference library

# load refs
custom.refs <- read_csv(here("meta-fish-pipe/temp/taxonomic-assignment/custom-reference-library.csv"),guess_max=99999,show_col_types=FALSE,lazy=FALSE)

# convert
writeLines("\nWriting out fasta and running hmmer3 - may take a few minutes.")
ref.fas <- tab2fas(custom.refs,seqcol="nucleotides",namecol="dbid")

# write out to temp file for hmm
ape::write.FASTA(ref.fas,file=here("temp/custom-refs-temp.fasta"))

# subset the markers with hmmer
setwd(here("temp"))
# with primers (taberlet is inside Miya)
miya.fas.primers <- run_hmmer3(dir=here("temp"), infile="custom-refs-temp.fasta", prefix="12s.miya.primers", evalue="10", coords="env")
setwd(here())

# convert to char
miya.char.primers <- lapply(as.character(miya.fas.primers),paste,collapse="")
custom.refs %<>% mutate(
    nucleotidesFrag.12s.miya.primers=as.character(miya.char.primers[match(dbid,names(miya.char.primers))]), 
    nucleotidesFrag.12s.miya.primers=str_replace_all(nucleotidesFrag.12s.miya.primers,"NULL","NA"),
    nucleotidesFrag.12s.miya.primers=dplyr::na_if(nucleotidesFrag.12s.miya.primers,"NA"),
    lengthFrag.12s.miya.primers=str_length(nucleotidesFrag.12s.miya.primers)
    )

# prefixes for subset
prefix <- "12s.miya.primers"

# subset the marker from the custom reference library reflib - using 95% length cutoff
reflib.sub <- subset_by_marker(prefix=prefix,df=custom.refs,thresh=0.95)

# collapse haps
reflib.red <- reflib.sub %>% group_by(sciNameValid) %>% group_modify(~ hap_collapse_df(df=.x,lengthcol=paste0("lengthFrag.",prefix),nuccol=paste0("nucleotidesFrag.",prefix),cores=1)) %>% ungroup()

# convert to fasta
reflib.fas <- tab2fas(reflib.red,seqcol=paste0("nucleotidesFrag.",prefix),namecol="dbid")

# remove Ns
reflib.fas <- rm_ns(bin=reflib.fas)

# align
writeLines("\nNow aligning sequences and running DECIPHER - may take a few minutes.")
reflib.fas.ali <- ips::mafft(reflib.fas,exec="mafft",method="retree 2",maxiterate=2)

# convert to char and filter out those with gap chars at the ends
reflib.fas.ali.char <- lapply(as.character(as.list(reflib.fas.ali)),paste,collapse="")
starts <- reflib.fas.ali.char[!grepl("^-",reflib.fas.ali.char)]
starts.ends <- starts[!grepl("-$",starts)]
#length(starts.ends)

# now filter those out of reflib
reflib.red %<>% filter(dbid %in% names(starts.ends))

# prepare the reference library as named character vector
seqs.dealign <- toupper(str_replace_all(starts.ends,"-",""))
names(seqs.dealign) <- names(starts.ends)
#length(seqs.dealign)

# set primers
# tele02
tele02 <- c("AAACTCGTGCCAGCCACC","GGGTATCTAATCCCAGTTTG")
# mifish-u-mod
mifish.u.mod <- c("GCCGGTAAAACTCGTGCCAGC","CATAGTGGGGTATCTAATCCCAGTTTG")
# mifish-u
mifish.u <- c("GTCGGTAAAACTCGTGCCAGC","CATAGTGGGGTATCTAATCCCAGTTTG")


# run AmplifyDNA for each primer
suppressPackageStartupMessages(
    library("DECIPHER")
    )
# run on 8 cores
tele02.res <- mcmapply(function(x) extract_eff(x,ppair=tele02,temp=54,proc=1), seqs.dealign, mc.cores=8, SIMPLIFY=TRUE, USE.NAMES=TRUE)
mifish.u.mod.res <- mcmapply(function(x) extract_eff(x,ppair=mifish.u.mod,temp=60,proc=1), seqs.dealign, mc.cores=8, SIMPLIFY=TRUE, USE.NAMES=TRUE)
mifish.u.res <- mcmapply(function(x) extract_eff(x,ppair=mifish.u,temp=60,proc=1), seqs.dealign, mc.cores=8, SIMPLIFY=TRUE, USE.NAMES=TRUE)
detach("package:DECIPHER")

# combine the dataframes
effic.df <- dplyr::bind_rows(
    enframe(tele02.res,name="dbid",value="efficiency") %>% mutate(primerSet="tele02"),
    enframe(mifish.u.mod.res,name="dbid",value="efficiency") %>% mutate(primerSet="mifish-u-mod"),
    enframe(mifish.u.res,name="dbid",value="efficiency") %>% mutate(primerSet="mifish-u")
    )

# format and write out
effic.df %>% mutate(species=pull(custom.refs,sciNameValid)[match(dbid,pull(custom.refs,dbid))]) %>% 
    arrange(primerSet,species,desc(efficiency)) %>% 
    select(primerSet,dbid,species,efficiency) %>%
    write_csv(here("temp/results/primer-efficiency-results.csv"))

# report
writeLines("\nPrimer efficiencies written to 'temp/results/primer-efficiency-results.csv'.")
