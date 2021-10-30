#!/usr/bin/env Rscript


# load libs
# install a bioc package with renv
# renv::install(packages="bioc::DECIPHER")
suppressPackageStartupMessages({
    library("here")
    library("parallel")
    library("tidyverse")
    library("magrittr")
    library("ape")
    library("ips")
})


# function to reverse complement DNA
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/revcomp_dna.R")

# function for making fasta files from tables
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/tab2fas.R")


# function to subset the reference dataframe by marker and filter by sequence length
subset_by_marker <- function(prefix,df,thresh){
    df1 <- df %>% filter(!is.na(!!as.name(paste0("nucleotidesFrag.",prefix))))
    df1 %<>% filter(!!as.name(paste0("lengthFrag.",prefix)) >= (median(!!as.name(paste0("lengthFrag.",prefix)))*thresh))
    return(df1)
}


# collapses haplotypes (from a dataframe format to a dataframe format)
hap_collapse_df <- function(df,lengthcol,nuccol,cores){
    odf <- df[order(df[[lengthcol]],decreasing=TRUE),]
    reps <- mcmapply(FUN=function(x) which(str_detect(string=odf[[nuccol]], pattern=x) == TRUE)[1], odf[[nuccol]], SIMPLIFY=TRUE, USE.NAMES=FALSE, mc.cores=cores)
    ind <- unique(reps)
    dat <- odf[ind,]
    dat[["nHaps"]] <- as.numeric(table(reps))
    return(dat)
}


# function to get retrieve species names of sequences with an identical haplotype as your query 
get_sames <- function(df,ids,nucs,sppVec,query){
    per.ind <- df[[sppVec]][str_detect(df[[nucs]], query)]
    return(per.ind)
}


# removes Ns from a DNAbin list
rm_ns <-function(bin){
    bin.char <- as.character(bin)
    bin.spread <- sapply(bin.char, paste, collapse="")
    bin.rep <- sapply(bin.spread, str_replace_all, "[^actg]", "")
    bin.split <- strsplit(bin.rep, "")
    bin.bin <- as.DNAbin(bin.split)
    return(bin.bin)
}


# R script to run a hidden markov model on a sequence
run_hmmer3 <- function(dir, infile, hmm, prefix, evalue, coords){#
    string.hmmer <- paste0("nhmmer -E ", evalue, " --incE ", evalue, " --dfamtblout ", dir, "/", prefix, ".hmmer.tbl ", "../assets/", prefix, ".hmm ", dir, "/", infile)
    system(command=string.hmmer, ignore.stdout=TRUE)
    hmm.tbl <- readr::read_table(file=paste0(dir, "/", prefix, ".hmmer.tbl"), col_names=FALSE, progress=FALSE, comment="#", col_types=cols(), guess_max=100000)
    names(hmm.tbl) <- c("targetName","acc","queryName","bits","eValue","bias","hmmStart","hmmEnd","strand","aliStart","aliEnd","envStart","envEnd","sqLen","descriptionTarget")
    hmm.tbl %<>% filter(strand=="+") %>% distinct(targetName, .keep_all=TRUE) %>% mutate(coords=paste(envStart,envEnd,sep=":"))
    mtdna <- read.FASTA(file=paste0(dir,"/",infile))
    mtdna.sub <- as.character(mtdna[match(hmm.tbl$targetName,names(mtdna))])
    if(coords=="env"){
    mtdna.sub.coords <- as.DNAbin(mapply(function(x,y,z) x[y:z], mtdna.sub, hmm.tbl$envStart, hmm.tbl$envEnd, SIMPLIFY=TRUE, USE.NAMES=TRUE))
    } else if(coords=="ali"){
    mtdna.sub.coords <- as.DNAbin(mapply(function(x,y,z) x[y:z], mtdna.sub, hmm.tbl$aliStart, hmm.tbl$aliEnd, SIMPLIFY=TRUE, USE.NAMES=TRUE))
    } else {
    stop("Please provide 'env' or 'ali' as arguments to coords")
    }
    return(mtdna.sub.coords)
}


# function to run AmplifyDNA to get primer efficiencies
extract_eff <- function(qseq,ppair,temp,proc){
    res <- DECIPHER::AmplifyDNA(primers=ppair,myDNAStringSet=qseq,maxProductSize=300,annealingTemp=temp,P=5e-7,minEfficiency=0,taqEfficiency=TRUE,maxDistance=0.4,maxGaps=2,processors=proc)
    if(length(res)==0) {out <- NA}
    else {out <- res@ranges@NAMES[1]}
    out <- as.numeric(str_replace_all(out,"%.*",""))/100
    names(out) <- names(qseq)
    return(out)
}
