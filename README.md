[![DOI](https://zenodo.org/badge/362065020.svg)](https://zenodo.org/badge/latestdoi/362065020)

# Reproduction explains marine eDNA variation

Code for: Collins, R.A., Baillie, C., Halliday, N.C., Rainbird, S., Sims, D.W., Mariani, S. & Genner, M.J. (202x) Reproduction explains marine eDNA variation. _Insert Journal_. [https://doi.org/xxx](https://doi.org/xxx).

To support this code, a full data description and summary can be found at [assets/data-description.pdf](assets/data-description.pdf).


### 1. Install software

* Ensure you have `git`, `mafft` and `hmmer` installed on your system (the [meta-fish-pipe](https://github.com/genner-lab/meta-fish-pipe) and [refseq-reflib](https://github.com/genner-lab/refseq-reflib) pipelines will require their own software also - please refer to their install instructions). For Ubuntu users, run:

```bash
# run in bash
sudo apt install git mafft hmmer
```


### 2. Clone pipeline modules

* In bash terminal clone the pipeline modules using git, and install required R packages in each module:

```bash
# run in bash
git clone https://github.com/genner-lab/reproduction-temporal-trends.git
cd reproduction-temporal-trends
Rscript -e "renv::restore()"
git clone https://github.com/genner-lab/meta-fish-pipe.git
cd meta-fish-pipe
git checkout v1.1
Rscript -e "renv::restore()"
cd ..
git clone https://github.com/genner-lab/refseq-reflib.git
cd refseq-reflib
git checkout v1.1
Rscript -e "renv::restore()"
cd ..
```


### 3. Prepare reference libraries

* Open an R session and download the UK custom fish reference library from [meta-fish-lib](https://github.com/genner-lab/meta-fish-lib):

```r
# run in R
library("tidyverse")
library("ape")
source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-load-remote.R")
source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-clean.R")
locals <- read_csv(file="assets/local-12s.csv")
reflib.orig %>% 
    bind_rows(locals) %>% 
    mutate(subphylum="Chordata") %>% 
    rename(phylum=subphylum) %>% 
    write_csv(file="meta-fish-pipe/assets/meta-fish-lib-v245.csv")
```

* Obtain the RefSeq general reference library from [refseq-reflib](https://github.com/genner-lab/refseq-reflib), running the following code in the bash terminal:

```bash
# run in bash
cd refseq-reflib
mkdir temp references
scripts/download.sh
scripts/extract.R -p tele02
scripts/annotate.R -s 42 -p tele02
cd ..
cp refseq-reflib/references/refseq208-annotated-tele02.csv meta-fish-pipe/assets/refseq208-annotated-tele02.csv
```


### 4. Retrieve fastq data and create required files

* Retrieve the fastq datasets from NCBI Sequence Read Archive (SRA project [PRJNA725897](https://www.ncbi.nlm.nih.gov/sra/PRJNA725897)) and verify the md5sums:

```bash
# run in bash
scripts/get-data.sh
```

* Copy files required for running [meta-fish-pipe](https://github.com/genner-lab/meta-fish-pipe) and set up session:

```bash
# run in bash
cd meta-fish-pipe
cp ../assets/sequencing-master.csv assets/sequencing-master.csv
cp ../assets/contaminants-exclude.csv assets/contaminants-exclude.csv
scripts/session-info.sh  -r assets/refseq208-annotated-tele02.csv -c assets/meta-fish-lib-v245.csv
```

* Prepare directory structures for each library:

```bash
# run in bash
scripts/prepare-libraries.sh -p tele02 -l lib1
scripts/prepare-libraries.sh -p tele02 -l lib2
scripts/prepare-libraries.sh -p tele02 -l lib3
scripts/prepare-libraries.sh -p tele02 -l lib4
```

* Create symlinks for the fastq datasets:

```bash
# run in bash
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
```

* Generate sample demultiplexing files:

```bash
# run in bash
cd meta-fish-pipe
scripts/generate-barcodes.R -p tele02 -l lib1 -f 18 -r 20 -m assets/sequencing-master.csv
scripts/generate-barcodes.R -p tele02 -l lib2 -f 18 -r 20 -m assets/sequencing-master.csv
scripts/generate-barcodes.R -p tele02 -l lib3 -f 18 -r 20 -m assets/sequencing-master.csv
scripts/generate-barcodes.R -p tele02 -l lib4 -f 18 -r 20 -m assets/sequencing-master.csv
```


### 5. Run processing pipeline:

* Demultiplex each library with cutadapt:

```bash
# run in bash
scripts/demultiplex.sh -p tele02 -l lib1 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18
scripts/demultiplex.sh -p tele02 -l lib2 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18
scripts/demultiplex.sh -p tele02 -l lib3 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18
scripts/demultiplex.sh -p tele02 -l lib4 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18
```

* Denoise each library with dada2:

```bash
# run in bash
scripts/dada2.R -p tele02 -l lib1
scripts/dada2.R -p tele02 -l lib2
scripts/dada2.R -p tele02 -l lib3
scripts/dada2.R -p tele02 -l lib4
```

* Generate pipeline statistics:

```bash
# run in bash
scripts/generate-stats.sh -p tele02 -l lib1 -t 8
scripts/generate-stats.sh -p tele02 -l lib2 -t 8
scripts/generate-stats.sh -p tele02 -l lib3 -t 8
scripts/generate-stats.sh -p tele02 -l lib4 -t 8
```

* Run taxonomic assignment:

```bash
# run in bash
scripts/taxonomic-assignment.sh -t 8 -p tele02
```

* Assemble results:

```bash
# run in bash
scripts/assemble-results.R -c assets/contaminants-exclude.csv
```


### 6. Run statistical analyses

* Calculate primer efficiencies:

```bash
# run in bash
cd ..
scripts/primer-fit.R
```

* Generate the seasonal PCoA plots for eDNA and demersal trawl:

```bash
# run in bash
scripts/seasonal-trends.R
```

* Generate the logistic regression plots for fish spawning months and eDNA abundance:

```bash
# run in bash
scripts/spawning.R
```

* Generate the abundance correlation plots for eDNA and demersal trawl for adult fishes, as well the results of the  zero-inflated negative-binomial generalised-linear mixed-model:

```bash
# run in bash
scripts/adults.R
```

* Generate the abundance correlation plots for eDNA and larval ichthyoplankton:

```bash
# run in bash
scripts/larvae.R
```

* Generate the abundance correlation plots for eDNA and egg ichthyoplankton:

```bash
# run in bash
scripts/eggs.R
```

* Generate the supporting information plots and tables:

```bash
# run in bash
scripts/supporting-info.R 
```
