# Reproduction explains environmental DNA variation in a temperate marine fish community

Code for the marine eDNA fish reproduction paper published as:

Collins, R.A., Baillie, C., Halliday, N.C., Rainbird, S., Sims, D.W., Mariani, S. & Genner, M.J. (2021) Reproduction explains environmental DNA variation
in a temperate marine fish community. _Insert Journal_. [https://doi.org/xxx](https://doi.org/xxx).

## Organise and document

Important note: unless otherwise stated, run R from the project base directory, but run bash terminal from the script location `scripts`.

### 1: Install and record software installed

* In bash terminal (from base directory) create a temp project file with `mkdir temp`.

* In R, install R packages with renv.

```r
# after cloning the repository
renv::restore()
```

* In R run:

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

* In a bash terminal (in project base directory), run: 

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

Only needs to be done once at the beginning. Best to run scripts line by line rather than excecuting them, in case of errors. Restart R before running each step.

* Run `scripts/download-refseq.sh` to download the REFSEQ reference library.

* Run `scripts/assemble-refseq.R` to annotate the references. 

* Run `scripts/join-references.R` to combine the custom genbank ((https://doi.org/10.5281/zenodo.4443447)[https://doi.org/10.5281/zenodo.4443447]) and the local SeaDNA reference libraries.

* Run `scripts/primer-fit.R` to generate primer efficiency scores.


## Process libraries

These steps are carried out for each sequencing library individually. Each script  at the top, has a choice of libraries to run (lib1, lib2, lib3 etc) , and they must be run one at a time, after restarting R or a new bash session.

### 3: Prepare barcodes for demultiplexing

* Run `scripts/prep-barcodes.R` to generate sample barcodes for demultiplexing.


### 4: Demultiplex reads with cutadapt

* Run `scripts/demultiplex.sh` to reorientate, demultiplex, trim and obtain stats.


### 5: Denoise reads

* Run `scripts/dada2.R` to denoise reads and generate OTU tables.


### 6: Assign taxonomy

The taxonomy assignment step is a bit of a pain because in requires jumping back and forth between a bash (`scripts/taxonomic-assignment.sh`) and an R script (`scripts/taxonomic-assignment.R`) to process the data at different steps, and for each library. Follow the step numbers in the scripts. Process each library one at a time through all steps.

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

* [R] Get sample level assignments (Step 12 in `scripts/taxonomic-assignment.R`). This step additionally requires the `plates` and `events.master` objects to have been created by running the `scripts/prep-barcodes.R` script (for same libraries).


### 7: Combine assignment results for all libraries

* [R] Combine taxonomic assignment results for all libraries (Step 13 in `scripts/taxonomic-assignment.R`). Must run taxonomic assignment for all libraries first before combining results.


## Statistical analysis


### 8: Season trends

* Run `scripts/seasonal-trends.R` to generate the seasonal PCoA plots for eDNA and demersal trawl.

### 9: Spawning 

* Run `scripts/spawning.R` to generate the logistic regression plots for fish spawning months and eDNA abundance.

### 10: Adundance correlations & models

* Run `scripts/adults.R` to generate the abundance correlation plots for eDNA and demersal trawl for adult fishes, as well the results of the  zero-inflated negative-binomial generalised-linear mixed-model.

* Run `scripts/larvae.R` to generate the abundance correlation plots for eDNA and larval ichthyoplankton.

* Run `scripts/eggs.R` to generate the abundance correlation plots for eDNA and egg ichthyoplankton.

### 11: Supporting information

* Run `scripts/supporting-info.R` to generate the supporting info plots and tables.

* Run `scripts/maps.R` to generate the supporting info maps.
