#!/usr/bin/env Rscript

# install a bioc package with renv
# renv::install(packages="bioc::dada2")

# load libs
library("here")
library("parallel")
library("tidyverse")
library("magrittr")
library("lubridate")
library("openssl")
library("ape")
library("ips")
library("spider")
library("phangorn")
library("dada2")
library("taxadb")
# sudo R devtools::install_github("ropensci/taxadb")


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
# need to specify columns that contain sequence lengths, and nucleotides
# hap_collapse_df(df=mydataframe,lengthcol="lengthFrag",nuccol="nucleotidesFrag")
# add a number of each haplotype
hap_collapse_df <- function(df,lengthcol,nuccol,cores){
    odf <- df[order(df[[lengthcol]],decreasing=TRUE),]
    reps <- mcmapply(FUN=function(x) which(str_detect(string=odf[[nuccol]], pattern=x) == TRUE)[1], odf[[nuccol]], SIMPLIFY=TRUE, USE.NAMES=FALSE, mc.cores=cores)
    ind <- unique(reps)
    dat <- odf[ind,]
    dat[["nHaps"]] <- as.numeric(table(reps))
    return(dat)
}


# function to get retrieve species names of sequences with an identical haplotype as your query 
# works on a dataframe
# get_sames(df=mydataframe,ids="dbid",nucs="nucleotidesFrag",sppVec="sciNameValid",query=mydataframe$nucleotidesFrag[[1]])
get_sames <- function(df,ids,nucs,sppVec,query){
    per.ind <- df[[sppVec]][str_detect(df[[nucs]], query)]
    return(per.ind)
}


# mode function
mode_avg <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
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


# function to create short file paths
cpath <- function(sense,step,r){
    path <- paste0(proj.path,"/processed-reads/",sense,"/",step,"-",r)
    return(path)
}


# R script to run a hidden markov model on a sequence
# need to have "hmmer" and "biosquid" installed 
# if not, run 'sudo apt install hmmer biosquid'
# also needs the ape package
# requires a tempfile directory (e.g. "../temp")
# requires an infile in fasta format to be in the same dir as the tempfiles (e.g. "myfile.fas")
# requires the name of the hmm you want to use (e.g. "12s.miya.noprimers.hmm")
# requires a prefix for the hmmer output (e.g. "analysis1")
# assumes the hidden markov model is located in ../assets
# returns a DNAbin object of the sequences matched by hmmer 

run_hmmer3 <- function(dir, infile, hmm, prefix, evalue, coords){#
    string.hmmer <- paste0("nhmmer -E ", evalue, " --incE ", evalue, " --dfamtblout ", dir, "/", prefix, ".hmmer.tbl ", "../assets/", prefix, ".hmm ", dir, "/", infile)
    system(command=string.hmmer, ignore.stdout=TRUE)
    hmm.tbl <- read_delim(file=paste0(dir, "/", prefix, ".hmmer.tbl"), delim=" ", col_names=FALSE, trim_ws=TRUE, progress=FALSE, comment="#", col_types=cols(), guess_max=100000)
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



# fun to subset a reference lib for each marker
subset_nucs <- function(pref,df){
    df %<>% rename(nucleotidesFrag=!!as.name(paste0("nucleotidesFrag.",pref)), lengthFrag=!!as.name(paste0("lengthFrag.",pref)))
    df %<>% filter(!is.na(nucleotidesFrag))
    return(df)
}


# fun to annotate a reference library table with number haplotypes per species
haps2fas <- function(df){
    df <- bind_rows(mcmapply(FUN=function(x) hap_collapse_df(df=x,lengthcol="lengthFrag",nuccol="nucleotidesFrag",cores=1), split(df,pull(df,sciNameValid)), SIMPLIFY=FALSE,mc.cores=1))
    sames <- mclapply(FUN=function(x) get_sames(df=df,ids="dbid",nucs="nucleotidesFrag",sppVec="sciNameValid",query=x), pull(df,nucleotidesFrag), mc.cores=1)
    df %<>% mutate(nMatches=sapply(sames, function(x) length(unique(x))), matchTax=sapply(sames, function(x) paste(unique(x),collapse=" | ")))
    df %<>% mutate(noms=paste(dbid,str_replace_all(sciNameValid," ","_"),nHaps,sep="|")) %>% arrange(class,order,family,genus,sciNameValid,lengthFrag,dbid)
    return(df)
}


# fun to align seqs and make a phylogentic tree
phylogenize <- function(fas,prefix,binLoc,outdir){
    fas <- ips::mafft(fas,exec="mafft",method="retree 2",maxiterate=2)
    tr <- ips::raxml(fas, file=paste0("fromR-",prefix), m="GTRCAT", f="d", p=42, exec=binLoc, N=1)
    tr <- tr$bestTree
    tmp.path <- paste0(outdir,"/qc_",paste(month(ymd(Sys.Date()),label=TRUE),year(ymd(Sys.Date())),sep="-"))
    dir.create(path=tmp.path)
    flist <- list.files(pattern=prefix)
    file.copy(flist, paste0(tmp.path,"/",flist))
    file.remove(flist)
    write.tree(tr,file=paste0(tmp.path,"/",prefix,".nwk"))
    return(tr)
}


# fun to plot and annotate phylogenetic trees
plot_trees <- function(tr,df,prefix,outdir){
    tr <- ape::ladderize(phangorn::midpoint(tr))
    sppv <- pull(df,sciNameValid)[match(str_split_fixed(tr$tip.label,"\\|",3)[,1],pull(df,dbid))]
    monov <- spider::monophyly(tr,sppVector=sppv)
    allmono <- monov[match(sppv, unique(sppv))]
    cols <- rep("gray20",length(tr$tip.label))
    cols[which(allmono==FALSE)] <- "hotpink"
    cols[match(df$noms[which(df$nMatches>1)], tr$tip.label)] <- "green3"
    tmp.path <- paste0(outdir,"/qc_",paste(month(ymd(Sys.Date()),label=TRUE),year(ymd(Sys.Date())),sep="-"))
    dfs <- df %>% summarise(nSeqs=sum(nHaps),nHaps=length(nHaps),nSpp=length(unique(sciNameValid)))
    tit <- paste0(str_replace_all(prefix,"\\.noprimers",""),"\n(n=",pull(dfs,nSeqs),", n haps=",pull(dfs,nHaps),", n spp.=",pull(dfs,nSpp),")\npink = non-monophyletic species\ngreen = shared haplotypes\nscroll down for tree ...")
    pdf(file=paste0(tmp.path,"/RAxML_bestTree.",prefix,".pdf"), width=15, height=length(tr$tip.label)/10)
    plot.phylo(tr, tip.col=cols, cex=0.5, font=1, label.offset=0.01, no.margin=TRUE)
    title(tit, line=-10)
    dev.off()
}


# function to query for contamination in an otu table
find_matches <- function(query,referencedb){
    res <- pull(referencedb,label)[which(grepl(query,pull(referencedb,nucleotides)))]
    if(length(res)==0){
    res <- NA} else {res}
    return(res)
}
