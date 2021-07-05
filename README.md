[![DOI](https://zenodo.org/badge/362065020.svg)](https://zenodo.org/badge/latestdoi/362065020)

```
# new code and instructions - july 2021

# clone repos
git clone https://github.com/genner-lab/reproduction-temporal-trends.git
cd reproduction-temporal-trends
git clone https://github.com/genner-lab/meta-fish-pipe.git
git clone https://github.com/genner-lab/refseq-reflib.git

# restore R libs
Rscript -e "renv::restore()"
# get data from ncbi
scripts/get-data.sh

# in R get the uk fish reference library (and add local seqs)
library("tidyverse")
library("ape")
source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-load-remote.R")
source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-clean.R")
locals <- read_csv(file="assets/local-12s.csv")
reflib.orig %>% bind_rows(locals) %>% write_csv(file="meta-fish-pipe/assets/meta-fish-lib-v243.csv")

# get refseq
cd refseq-reflib
mkdir temp references
Rscript -e "renv::restore()"
scripts/download.sh
scripts/extract.R -p tele02
scripts/annotate.R -s 42 -p tele02
rm temp/duckdb
cd ..
cp refseq-reflib/references/refseq206-annotated-tele02.csv meta-fish-pipe/assets/refseq206-annotated-tele02.csv


# copy across sample sheet and contam file to the pipeline lib
cp assets/sequencing-master.csv meta-fish-pipe/assets/sequencing-master.csv
cp assets/contaminants-exclude.csv meta-fish-pipe/assets/contaminants-exclude.csv


# set up pipeline
cd meta-fish-pipe
Rscript -e "renv::restore()"
scripts/session-info.sh  -r assets/refseq206-annotated-tele02.csv -c assets/meta-fish-lib-v243.csv

# set up libs
scripts/prepare-libraries.sh -p tele02 -l lib1
scripts/prepare-libraries.sh -p tele02 -l lib2
scripts/prepare-libraries.sh -p tele02 -l lib3
scripts/prepare-libraries.sh -p tele02 -l lib4

cd ..
# make symlinks LIB1
ln -s -r temp/data/SeaDNA_Teleo02_01_S1_L001_R1_001.fastq.gz meta-fish-pipe/temp/processing/tele02-lib1/fastq/R1.fastq.gz
ln -s -r temp/data/SeaDNA_Teleo02_01_S1_L001_R2_001.fastq.gz meta-fish-pipe/temp/processing/tele02-lib1/fastq/R2.fastq.gz
# make symlinks LIB2
ln -s -r temp/data/SeaDNA_Teleo02_02_S2_L001_R1_001.fastq.gz meta-fish-pipe/temp/processing/tele02-lib2/fastq/R1.fastq.gz
ln -s -r temp/data/SeaDNA_Teleo02_02_S2_L001_R2_001.fastq.gz meta-fish-pipe/temp/processing/tele02-lib2/fastq/R2.fastq.gz
# make symlinks LIB3
ln -s -r temp/data/SeaDNA_Tele02_Lib03v2_R1.fastq.gz.1 meta-fish-pipe/temp/processing/tele02-lib3/fastq/R1.fastq.gz
ln -s -r temp/data/SeaDNA_Tele02_Lib03v2_R2.fastq.gz.1 meta-fish-pipe/temp/processing/tele02-lib3/fastq/R2.fastq.gz
# make symlinks LIB4
ln -s -r temp/data/SeaDNA_Teleo02_Lib-04_S2_L001_R1_001.fastq.gz meta-fish-pipe/temp/processing/tele02-lib4/fastq/R1.fastq.gz
ln -s -r temp/data/SeaDNA_Teleo02_Lib-04_S2_L001_R2_001.fastq.gz meta-fish-pipe/temp/processing/tele02-lib4/fastq/R2.fastq.gz

# copy across sample sheet and contam file to the pipeline lib
cp assets/sequencing-master.csv meta-fish-pipe/assets/sequencing-master.csv
cp assets/contaminants-exclude.csv meta-fish-pipe/assets/contaminants-exclude.csv

# generate barcodes
cd meta-fish-pipe
scripts/generate-barcodes.R -p tele02 -l lib1 -f 18 -r 20 -m assets/sequencing-master.csv
scripts/generate-barcodes.R -p tele02 -l lib2 -f 18 -r 20 -m assets/sequencing-master.csv
scripts/generate-barcodes.R -p tele02 -l lib3 -f 18 -r 20 -m assets/sequencing-master.csv
scripts/generate-barcodes.R -p tele02 -l lib4 -f 18 -r 20 -m assets/sequencing-master.csv

# demultiplex
scripts/demultiplex.sh -p tele02 -l lib1 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18
scripts/demultiplex.sh -p tele02 -l lib2 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18
scripts/demultiplex.sh -p tele02 -l lib3 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18
scripts/demultiplex.sh -p tele02 -l lib4 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18

# denoise with dada2
scripts/dada2.R -p tele02 -l lib1
scripts/dada2.R -p tele02 -l lib2
scripts/dada2.R -p tele02 -l lib3
scripts/dada2.R -p tele02 -l lib4

# generate stats
scripts/generate-stats.sh -p tele02 -l lib1 -t 8
scripts/generate-stats.sh -p tele02 -l lib2 -t 8
scripts/generate-stats.sh -p tele02 -l lib3 -t 8
scripts/generate-stats.sh -p tele02 -l lib4 -t 8

# run taxonomic assignment
scripts/taxonomic-assignment.sh -t 8 -p tele02

# assemble results
scripts/assemble-results.R -c assets/contaminants-exclude.csv


```






# Reproduction explains marine eDNA variation

Code for: Collins, R.A., Baillie, C., Halliday, N.C., Rainbird, S., Sims, D.W., Mariani, S. & Genner, M.J. (2021) Reproduction explains marine eDNA variation. _Insert Journal_. [https://doi.org/xxx](https://doi.org/xxx).

To support this code, a full data description and summary can be found at [assets/data-description.pdf](assets/data-description.pdf).

## Important notes

* Unless otherwise stated, run R from the project base directory, but run bash terminal from the script location `scripts/`.

* It is generally best to run scripts line-by-line rather than excecuting them, in case of errors.

* Always restart R/bash before running each step.

## Organise and document

### 1: Clone repo, install software, report versions

* In bash terminal, clone repository with `git clone https://github.com/genner-lab/reproduction-temporal-trends.git`. 

* Change to base directory with `cd reproduction-temporal-trends`, and create a temporary output file with `mkdir temp`.

* Install the following software on your system, and make them available on your $PATH: [cutadapt](https://github.com/marcelm/cutadapt) v2.10, [vsearch](https://github.com/torognes/vsearch) v2.15, [fqtools](https://github.com/alastair-droop/fqtools) v2.3, [raxml-ng](https://github.com/amkozlov/raxml-ng) v1.0.2, [epa-ng](https://github.com/Pbdas/epa-ng) v0.3.7, [gappa](https://github.com/lczech/gappa) v0.6.1, [hmmer](http://hmmer.org/) v3.1b2, [ncbi-blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) v2.10.1, [mafft](http://mafft.cbrc.jp/alignment/software/) v7.271.

* In R, install R packages with renv:

```r
# install all from renv.lock
renv::restore()
```

* To report package versions to the `temp/sessionInfo.txt` file, in R run:

```r
libs <- sort(unique(unlist(lapply(list.files(path="scripts",pattern="\\.R"),function(x) grep("^library\\(.+\\)",readLines(here::here("scripts",x)),value=TRUE)))))
writeLines(libs,here::here("temp/libs.R"))
source(here::here("temp/libs.R"))
si <- sessionInfo()
date <- Sys.Date()
sink(here::here("temp/sessionInfo.txt"))
writeLines(paste("--------------------------------------------------\nR Package version list. Today's date is",date,"\n--------------------------------------------------\n\n"))
print(si)
sink()
quit(save="no")
```

* To report program versions to the `temp/sessionInfo.txt` file, in a bash terminal (in project base directory), run: 

```bash
{ printf "\n\n-----------------------\nLinux software versions\n-----------------------\n"; } >> temp/sessionInfo.txt
{ printf "\ncutadapt --version\n" & cutadapt --version; } >> temp/sessionInfo.txt
{ printf "\ngappa --version\n" & gappa --version; } >> temp/sessionInfo.txt
{ printf "\nepa-ng --version\n" & epa-ng --version; } >> temp/sessionInfo.txt
{ printf "\nraxml-ng --version\n" & raxml-ng --version | grep "RAxML-NG"; } >> temp/sessionInfo.txt
{ printf "\nhmmsearch -h\n" & hmmsearch -h | grep "HMMER"; } >> temp/sessionInfo.txt
{ printf "\nblastn -version\n" & blastn -version | grep -v "Package"; } >> temp/sessionInfo.txt
{ printf "\nvsearch --version\n" & vsearch --version | grep "linux"; } >> temp/sessionInfo.txt
{ printf "\nfqtools -v\n" & fqtools -v; } >> temp/sessionInfo.txt
```
(`vsearch --version` currently broken, so add manually)


### 2: Prepare reference library

Only needs to be done once at the beginning.

* Run `scripts/download-refseq.sh` to download the REFSEQ reference library.

* Run `scripts/assemble-refseq.R` to annotate the references. 

* Run `scripts/join-references.R` to combine the custom genbank ([doi.org/10.5281/zenodo.4443447](https://doi.org/10.5281/zenodo.4443447)) and the local SeaDNA reference libraries.

* Run `scripts/primer-fit.R` to generate primer efficiency scores.


### 3: Retrieve data

* Run `scripts/get-data.sh` to retrieve the fastq datasets from NCBI Sequence Read Archive (SRA) and verify the md5sums.

* Project URL on SRA is [https://www.ncbi.nlm.nih.gov/sra/PRJNA725897](https://www.ncbi.nlm.nih.gov/sra/PRJNA725897).


## Process libraries

These steps are carried out for each sequencing library individually. Each script  at the top, has a choice of libraries to run (lib1, lib2, lib3, lib4) , and they must be run one at a time, after restarting R or a new bash session.


### 4: Prepare barcodes for demultiplexing

* Run `scripts/prep-barcodes.R` to generate sample barcodes for demultiplexing.


### 5: Demultiplex reads with cutadapt

* Run `scripts/demultiplex.sh` to reorientate, demultiplex, trim and obtain stats.


### 6: Denoise reads

* Run `scripts/dada2.R` to denoise reads and generate OTU tables.


### 7: Assign taxonomy

The taxonomy assignment step is a bit of a pain because in requires jumping back and forth between a bash (`scripts/taxonomic-assignment.sh`) and an R script (`scripts/taxonomic-assignment.R`) to process the data at different steps, and for each library. Follow the step numbers in the scripts. Process each library one at a time through all steps. Eventually I will turn this into one executable script, and run assignment only once on the combined outputs for all libraries.

* [bash] Start in `scripts/taxonomic-assignment.sh`, and create Blast and Sintax databases (Step 1 in `scripts/taxonomic-assignment.sh`). This step only needs to be run once for all libraries.

* [bash] Run the sequence homology HMM filter (Step 2 in `scripts/taxonomic-assignment.sh`).

* [bash] Run the Sintax taxonomic assignment (Step 3 in `scripts/taxonomic-assignment.sh`).

* [R] Summarise the HMM cleaned ASVs (Step 4 in `scripts/taxonomic-assignment.R`).

* [R] Summarise Sintax results and isolate fish reads (Step 5 in `scripts/taxonomic-assignment.R`).

* [bash] Run Blast taxonomic assignment (Step 6 in `scripts/taxonomic-assignment.sh`).

* [R] Process blast results (Step 7 in `scripts/taxonomic-assignment.R`).

* [R] Prepare input for EPA assignment (Step 8 in `scripts/taxonomic-assignment.R`).

* [bash] Run EPA assignment (Step 9 in `scripts/taxonomic-assignment.sh`).

* [R] Process EPA results (Step 10 in `scripts/taxonomic-assignment.R`).

* [R] Combine the taxonomic assignment results (Step 11 in `scripts/taxonomic-assignment.R`). 

* [R] Get sample level assignments (Step 12 in `scripts/taxonomic-assignment.R`). This step additionally requires the 'plates' and 'events.master' objects to have been created by running the `scripts/prep-barcodes.R` script (for same libraries).


### 8: Combine assignment results for all libraries

* [R] Combine taxonomic assignment results for all libraries (Step 13 in `scripts/taxonomic-assignment.R`). Must run taxonomic assignment for all libraries first before combining results.


## Statistical analysis


### 9: Seasonal trends

* Run `scripts/seasonal-trends.R` to generate the seasonal PCoA plots for eDNA and demersal trawl.


### 10: Spawning 

* Run `scripts/spawning.R` to generate the logistic regression plots for fish spawning months and eDNA abundance.


### 11: Adundance correlations & models

* Run `scripts/adults.R` to generate the abundance correlation plots for eDNA and demersal trawl for adult fishes, as well the results of the  zero-inflated negative-binomial generalised-linear mixed-model.

* Run `scripts/larvae.R` to generate the abundance correlation plots for eDNA and larval ichthyoplankton.

* Run `scripts/eggs.R` to generate the abundance correlation plots for eDNA and egg ichthyoplankton.


### 12: Supporting information

* Run `scripts/supporting-info.R` to generate the supporting info plots and tables.

* Run `scripts/maps.R` to generate the supporting info maps.
