#!/usr/bin/env Rscript

# load libs and funs
source(here::here("scripts/funs-bioinformatics.R"))
# check dada version
packageVersion("dada2")

# set the primer
primer.code <- "tele02"

# set the name of the folder that the file will be added to 
lib <- "lib1"
lib <- "lib2"
lib <- "lib3"
lib <- "lib4"

# make paths
proj.path <- here("temp/processing",paste0(primer.code,"-",lib))
#
dir.sense <- paste0(proj.path,"/processed-reads/sense")
dir.antisense <- paste0(proj.path,"/processed-reads/antisense")
# confirm proj path 
print(proj.path)


# trucLens
# for Miya MiFish - truncLen 105 gives at least 29 bp overlap for the longest amplicons (e.g. Raja clavata @ 181 bp), and 40 bp for the regular 170 bp
trucVal <- c(105,105)

# quality trim Ns and truncate
filterAndTrim(fwd=cpath("sense","trimmed","R1"), filt=cpath("sense","filtered","R1"), rev=cpath("sense","trimmed","R2"), filt.rev=cpath("sense","filtered","R2"), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, truncLen=trucVal, multithread=TRUE, verbose=TRUE, matchIDs=TRUE)
filterAndTrim(fwd=cpath("antisense","trimmed","R1"), filt=cpath("antisense","filtered","R1"), rev=cpath("antisense","trimmed","R2"), filt.rev=cpath("antisense","filtered","R2"), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, truncLen=trucVal, multithread=TRUE, verbose=TRUE, matchIDs=TRUE)

# learn errors (can restart from here)
set.seed(42)
sense.filt.R1.errs <- learnErrors(cpath("sense","filtered","R1"), multithread=TRUE, randomize=TRUE, nbases=1e+08, verbose=TRUE)
set.seed(42)
sense.filt.R2.errs <- learnErrors(cpath("sense","filtered","R2"), multithread=TRUE, randomize=TRUE, nbases=1e+08, verbose=TRUE)
set.seed(42)
antisense.filt.R1.errs <- learnErrors(cpath("antisense","filtered","R1"), multithread=TRUE, randomize=TRUE, nbases=1e+08, verbose=TRUE)
set.seed(42)
antisense.filt.R2.errs <- learnErrors(cpath("antisense","filtered","R2"), multithread=TRUE, randomize=TRUE, nbases=1e+08, verbose=TRUE)

# plot the errors
#plotErrors(dq=sense.filt.R1.errs, nominalQ=TRUE)
#plotErrors(dq=sense.filt.R2.errs, nominalQ=TRUE)
#plotErrors(dq=antisense.filt.R1.errs, nominalQ=TRUE)
#plotErrors(dq=antisense.filt.R2.errs, nominalQ=TRUE)


# make some fish priors
prefix <- "12s.taberlet.noprimers"

# load reflib - first need to have run 'join-references.R' to create
reflib <- read_csv(here("temp/reference-library/custom-refs.csv"),guess_max=99999,col_types=cols())

# subset the marker from the reflib and rev comp and trim
reflib.sub <- subset_by_marker(prefix=prefix,df=reflib,thresh=0.625) %>% 
    select(sciNameValid,!!as.name(paste0("nucleotidesFrag.",prefix))) %>% 
    rename(seqs=!!as.name(paste0("nucleotidesFrag.",prefix))) %>%
    filter(!str_detect(seqs,"[^actg]")) %>%
    mutate(fwd=toupper(str_trunc(seqs,width=trucVal[1],side="right",ellipsis=""))) %>%
    mutate(revcomp=toupper(str_trunc(revcomp_dna(seqs),width=trucVal[1],side="right",ellipsis=""))) 

# pull out the priors
fish.priors <- unique(c(pull(reflib.sub,fwd),pull(reflib.sub,revcomp)))

# run dada denoising - takes time with pool=TRUE
sense.filt.R1.dada <- dada(cpath("sense","filtered","R1"), err=sense.filt.R1.errs, multithread=TRUE, pool=TRUE, priors=fish.priors)
sense.filt.R2.dada <- dada(cpath("sense","filtered","R2"), err=sense.filt.R2.errs, multithread=TRUE, pool=TRUE, priors=fish.priors)
antisense.filt.R1.dada <- dada(cpath("antisense","filtered","R1"), err=antisense.filt.R1.errs, multithread=TRUE, pool=TRUE, priors=fish.priors)
antisense.filt.R2.dada <- dada(cpath("antisense","filtered","R2"), err=antisense.filt.R2.errs, multithread=TRUE, pool=TRUE, priors=fish.priors)

# look at dada objects
sense.filt.R1.dada[[1]]
sense.filt.R2.dada[[1]]
antisense.filt.R1.dada[[1]]
antisense.filt.R2.dada[[1]]

# derep
sense.filt.R1.derep <- derepFastq(cpath("sense","filtered","R1"))
sense.filt.R2.derep <- derepFastq(cpath("sense","filtered","R2"))
antisense.filt.R1.derep <- derepFastq(cpath("antisense","filtered","R1"))
antisense.filt.R2.derep <- derepFastq(cpath("antisense","filtered","R2"))

# merge the R1 and R2
sense.merged <- mergePairs(dadaF=sense.filt.R1.dada, derepF=sense.filt.R1.derep, dadaR=sense.filt.R2.dada, derepR=sense.filt.R2.derep, verbose=TRUE, maxMismatch=0)
antisense.merged <- mergePairs(dadaF=antisense.filt.R1.dada, derepF=antisense.filt.R1.derep, dadaR=antisense.filt.R2.dada,  derepR=antisense.filt.R2.derep, verbose=TRUE, maxMismatch=0)

head(sense.merged[[1]])
head(antisense.merged[[1]])

# make an OTU table
sense.seqtab <- makeSequenceTable(sense.merged)
antisense.seqtab <- makeSequenceTable(antisense.merged)
dim(sense.seqtab)
dim(antisense.seqtab)

# reverse comp the antisense
colnames(antisense.seqtab) <- dada2::rc(colnames(antisense.seqtab))

# fix the names before merging
rownames(sense.seqtab) <- str_split_fixed(rownames(sense.seqtab),"\\.",4)[,1]
rownames(antisense.seqtab) <- str_split_fixed(rownames(antisense.seqtab),"\\.",4)[,1]

# merge the tables
merged.seqtab <- mergeSequenceTables(table1=sense.seqtab, table2=antisense.seqtab, repeats="sum")
dim(merged.seqtab)

# remove chimaeras
merged.seqtab.nochim <- removeBimeraDenovo(merged.seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(merged.seqtab.nochim)

# stats to copy into 'logs/stats.csv'
# remember to run 'after filter' line in demultiplex.sh
table(nchar(getSequences(merged.seqtab.nochim)))
sum(merged.seqtab)# merge
sum(merged.seqtab.nochim)# chim
sum(merged.seqtab.nochim)/sum(merged.seqtab)

# make df and fasta for IDs
otus.df <- tibble(names=paste0("asv",str_pad(seq_along(colnames(merged.seqtab.nochim)),width=4,side="left",pad="0")), dnas=colnames(merged.seqtab.nochim)) %>% mutate(len=str_length(dnas))

# write out
write.FASTA(tab2fas(df=otus.df, seqcol="dnas", namecol="names"), file=paste0(proj.path,"/results/asvs.fna"))

# save the OTU table as df
colnames(merged.seqtab.nochim) <- paste0("asv",str_pad(seq_along(colnames(merged.seqtab.nochim)),width=4,side="left",pad="0"))
write_tsv(as_tibble(t(merged.seqtab.nochim), rownames="asv"), file=paste0(proj.path,"/results/asv-table.tsv"))
