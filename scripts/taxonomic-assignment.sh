#!/usr/bin/env/ sh

# cd to directory
# set project
PROJ="temp/processing"



#################################################### 
#################################################### 

### STEP 1: PREP REFERENCE LIBRARIES - ONLY NEED TO RUN ONCE ### 
### STEP 1: PREP REFERENCE LIBRARIES - ONLY NEED TO RUN ONCE ### 

# combine the refseq and the custom libraries (only need to run once)
cat ../temp/reference-library/refseq-annotated.fasta ../temp/reference-library/custom-refs-annotated.fasta > ../temp/reference-library/reflib-combined-annotated.fasta

# make blast db (only need to do this step once)
makeblastdb -in ../temp/reference-library/custom-refs.fasta -dbtype nucl -blastdb_version 5



#################################################### 
#################################################### 

#### STEP 2: FOR EACH LIB DO SEQUENCE HOMOLOGY HMM SEARCH TO CLEAN DNA ####
#### STEP 2: FOR EACH LIB DO SEQUENCE HOMOLOGY HMM SEARCH TO CLEAN DNA ####

# set marker
# tele02
PRIMER="tele02"

# set libs
LIB="lib1"
LIB="lib2"
LIB="lib3"
LIB="lib4"

# create marker lib combo and check
DIR="$PRIMER""-""$LIB"
echo "$DIR"

# hmms
HMM="12s.taberlet"

# subset just the fragment using an HMM
hmmsearch -E 0.01 --incE 0.01 ../assets/"$HMM".noprimers.hmm ../"$PROJ"/"$DIR"/results/asvs.fna | grep ">> asv[0-9+]" | sed -e 's/>> //g' -e 's/[[:space:]]//g' -e 's/$/$/' | sort | uniq > ../"$PROJ"/"$DIR"/results/hmm-out.txt
#
grep -A 1 -f ../"$PROJ"/"$DIR"/results/hmm-out.txt ../"$PROJ"/"$DIR"/results/asvs.fna | grep -v "-" > ../"$PROJ"/"$DIR"/results/asvs-clean.fasta
rm ../"$PROJ"/"$DIR"/results/hmm-out.txt



#################################################### 
#################################################### 

### STEP 3: FOR EACH LIB DO SINTAX SEARCH WITH REFSEQ ### 
### STEP 3: FOR EACH LIB DO SINTAX SEARCH WITH REFSEQ ### 

# cd back
cd ../"$PROJ"

# run sintax
vsearch --threads 8 --sintax "$DIR"/results/asvs-clean.fasta --db ../reference-library/reflib-combined-annotated.fasta --sintax_cutoff 0.7 --tabbedout "$DIR"/results/sintax-output.tsv



#################################################### 
#################################################### 

### STEP 4-5: FOR EACH LIB RUN `scripts/taxonomic-assignment.R` (PROCESS AND SUMMARISE) ### 
### STEP 4-5: FOR EACH LIB RUN `scripts/taxonomic-assignment.R` (PROCESS AND SUMMARISE) ### 



#################################################### 
#################################################### 

### STEP 6: BLAST REFERENCE DATABASE LOCALLY FOR ONLY FISH READS ###
### STEP 6: BLAST REFERENCE DATABASE LOCALLY FOR ONLY FISH READS ###


# blast the blast db
# get better hits with smaller word size
blastn -task blastn -num_threads 8 -evalue 1000 -word_size 7 -max_target_seqs 500 -db ../reference-library/custom-refs.fasta -outfmt "6 qseqid sseqid evalue length pident nident score bitscore" -out "$DIR"/results/fish-blast.out -query "$DIR"/results/fishqueries.fasta

# join the header
echo -e "asv\tblastDbid\tblastEvalue\tblastLength\tblastPident\tblastNident\tblastScore\tblastBitscore" > "$DIR"/results/headers
cat "$DIR"/results/headers "$DIR"/results/fish-blast.out > "$DIR"/results/fish-blast-result.tsv
rm "$DIR"/results/fish-blast.out
rm "$DIR"/results/headers



#################################################### 
#################################################### 

### STEP 7: FOR EACH LIB RUN `scripts/taxonomic-assignment.R` (PROCESS BLAST) ### 
### STEP 7: FOR EACH LIB RUN `scripts/taxonomic-assignment.R` (PROCESS BLAST) ### 

### STEP 8: FOR EACH LIB RUN `scripts/taxonomic-assignment.R` (PREP EPA INPUT) ### 
### STEP 8: FOR EACH LIB RUN `scripts/taxonomic-assignment.R` (PREP EPA INPUT) ### 



#################################################### 
#################################################### 

### STEP 9: EPA ASSIGNMENT ###
### STEP 9: EPA ASSIGNMENT ###

# https://github.com/amkozlov/raxml-ng
# https://github.com/lczech/gappa
# https://github.com/Pbdas/epa-ng

# make RAxML trees for all markers
# run these in different terminals
# create a binary and set model
raxml-ng --msa "$DIR"/results/epa/references.aligned.fasta --model TN93+G --parse

# run a tree search
raxml-ng --msa "$DIR"/results/epa/references.aligned.fasta.raxml.rba --tree pars{1} --search --seed 42 --threads 1

# optimise params
raxml-ng --msa "$DIR"/results/epa/references.aligned.fasta.raxml.rba --evaluate --tree "$DIR"/results/epa/references.aligned.fasta.raxml.rba.raxml.bestTree --prefix "$DIR"/results/epa/opt

# run epa-ng
epa-ng --ref-msa "$DIR"/results/epa/references.aligned.fasta --tree "$DIR"/results/epa/references.aligned.fasta.raxml.rba.raxml.bestTree --query "$DIR"/results/epa/queries.aligned.fasta --outdir "$DIR"/results/epa --model "$DIR"/results/epa/opt.raxml.bestModel --redo --preserve-rooting off

# run gappa
gappa examine assign --per-query-results --jplace-path "$DIR"/results/epa/epa_result.jplace --taxon-file "$DIR"/results/epa/references.taxonomy.tsv --out-dir "$DIR"/results/epa

# rename best tree to keep
mv "$DIR"/results/epa/references.aligned.fasta.raxml.rba.raxml.bestTree "$DIR"/results/epa/references.tree.nwk

# clean up
rm "$DIR"/results/epa/*.log
rm "$DIR"/results/epa/*raxml*



#################################################### 
#################################################### 

### STEP 10-13: RUN `scripts/taxonomic-assignment.R` (PROCESS EPA RESULTS AND COMBINE ALL RESULTS) ### 
### STEP 10-13: RUN `scripts/taxonomic-assignment.R` (PROCESS EPA RESULTS AND COMBINE ALL RESULTS) ### 
