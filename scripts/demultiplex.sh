#!/usr/bin/env sh
# add to .bashrc: `export PATH=~/.local/bin:$PATH`
# python3 -m pip install --upgrade pip
# python3 -m pip install cutadapt
# python3 -m pip install --upgrade cutadapt

# set project
PROJ="temp/processing"
# cd to directory
cd ../"$PROJ"

# set params #

# tele02
PRIMER="tele02"
FWD="AAACTCGTGCCAGCCACC"
REV="GGGTATCTAATCCCAGTTTG"
MINLEN="18"
TRUCVAL="105"

# set libs
LIB="lib1"
LIB="lib2"
LIB="lib3"
LIB="lib4"

# create marker lib combo
DIR="$PRIMER""-""$LIB"

# make a grepable degenerate primer
FWDGR="$(echo "$FWD" | sed -e 's/W/[A\|T]/g' -e 's/R/[A\|G]/g' -e 's/Y/[C\|T]/g' -e 's/M/[A\|C]/g' -e 's/N/[A\|C\|G\|T]/g' -e 's/$/\|\$/g')"
REVGR="$(echo "$REV" | sed -e 's/W/[A\|T]/g' -e 's/R/[A\|G]/g' -e 's/Y/[C\|T]/g' -e 's/M/[A\|C]/g' -e 's/N/[A\|C\|G\|T]/g' -e 's/$/\|\$/g')"


# make dirs for logs and sense/antisense, results etc
mkdir -p "$DIR"/fastq "$DIR"/logs "$DIR"/processed-reads/sense/dmplx "$DIR"/processed-reads/sense/trimmed-R1 "$DIR"/processed-reads/sense/trimmed-R2 "$DIR"/processed-reads/antisense/dmplx "$DIR"/processed-reads/antisense/trimmed-R1 "$DIR"/processed-reads/antisense/trimmed-R2 "$DIR"/results "$DIR"/trash

# make a stats file
# most of the stats can be found in cutadapt logs
echo -e "stat,reads\npf,\nprimer,\nbarcode,\ntrim,\nfilter,\nmerge,\nchim,\nhomol," > "$DIR"/logs/stats.csv

# make a symlink to the fastq files
# check md5sums
# check md5sums in data/sequencing-master.csv
md5sum "$DIR"/fastq/R1.fastq.gz
md5sum "$DIR"/fastq/R2.fastq.gz

# filter reads by orientation (need to wait for fwd to finish before doing rev)
# dumps all reads with fwd and rev primer on sense/antisense strands into "$DIR"/sense or "$DIR"/antisense

# sense
cutadapt --no-indels --error-rate 0 --overlap "$MINLEN" --pair-adapters --action=none -g senseF="$FWD" -G senseR="$REV" --untrimmed-output "$DIR"/trash/untrimmed.R1.fastq.gz --untrimmed-paired-output "$DIR"/trash/untrimmed.R2.fastq.gz -o "$DIR"/processed-reads/sense/R1.fastq.gz -p "$DIR"/processed-reads/sense/R2.fastq.gz "$DIR"/fastq/R1.fastq.gz "$DIR"/fastq/R2.fastq.gz > "$DIR"/logs/cutadapt.sense.log

# antisense
cutadapt --no-indels --error-rate 0 --overlap "$MINLEN" --pair-adapters --action=none -g antisenseF="$REV" -G antisenseR="$FWD" --untrimmed-output "$DIR"/trash/noprimer.R1.fastq.gz --untrimmed-paired-output "$DIR"/trash/noprimer.R2.fastq.gz -o "$DIR"/processed-reads/antisense/R1.fastq.gz -p "$DIR"/processed-reads/antisense/R2.fastq.gz "$DIR"/trash/untrimmed.R1.fastq.gz "$DIR"/trash/untrimmed.R2.fastq.gz > "$DIR"/logs/cutadapt.antisense.log

# check output and trash for primers (some of rev can be missing as rev primer longer than fwd)
gzip -cd "$DIR"/processed-reads/sense/R1.fastq.gz | sed -n '2~4p' | head -n 40 | GREP_COLOR="1;44" grep --color=always -E "$FWDGR" | GREP_COLOR="1;41" grep --color=always -E "$REVGR"
gzip -cd "$DIR"/processed-reads/sense/R2.fastq.gz | sed -n '2~4p' | head -n 40 | GREP_COLOR="1;44" grep --color=always -E "$FWDGR" | GREP_COLOR="1;41" grep --color=always -E "$REVGR"
gzip -cd "$DIR"/trash/noprimer.R1.fastq.gz | sed -n '2~4p' | head -n 40 | GREP_COLOR="1;44" grep --color=always -E "$FWDGR" | GREP_COLOR="1;41" grep --color=always -E "$REVGR"
gzip -cd "$DIR"/trash/noprimer.R2.fastq.gz | sed -n '2~4p' | head -n 40 | GREP_COLOR="1;44" grep --color=always -E "$FWDGR" | GREP_COLOR="1;41" grep --color=always -E "$REVGR"

# combinatorial paired demultiplex on sense
cutadapt --no-indels --error-rate 0 --overlap 10 --pair-adapters --action=none -g file:"$DIR"/barcodes-sense.fas -G file:"$DIR"/barcodes-antisense.fas --discard-untrimmed -o "$DIR"/processed-reads/sense/dmplx/{name}.R1.fastq.gz -p "$DIR"/processed-reads/sense/dmplx/{name}.R2.fastq.gz "$DIR"/processed-reads/sense/R1.fastq.gz "$DIR"/processed-reads/sense/R2.fastq.gz > "$DIR"/logs/cutadapt.dmplx.barcodes.sense.log

# combinatorial paired demultiplex on antisense
cutadapt --no-indels --error-rate 0 --overlap 10 --pair-adapters --action=none -g file:"$DIR"/barcodes-antisense.fas -G file:"$DIR"/barcodes-sense.fas --discard-untrimmed -o "$DIR"/processed-reads/antisense/dmplx/{name}.R1.fastq.gz -p "$DIR"/processed-reads/antisense/dmplx/{name}.R2.fastq.gz "$DIR"/processed-reads/antisense/R1.fastq.gz "$DIR"/processed-reads/antisense/R2.fastq.gz > "$DIR"/logs/cutadapt.dmplx.barcodes.antisense.log


# trim SENSE with cutadapt
sense="$(ls "$DIR"/processed-reads/sense/dmplx/*.fastq.gz | sed --expression='s/\.R.\.fastq\.gz//g' | sed --expression='s/unknown//g' | sed --expression="s/$DIR\/processed-reads\/sense\/dmplx\///g" | uniq)"
# check
for f in $sense; do echo "$f"; done
# now run
for i in $sense; do
cutadapt -n 5 --error-rate 0.15 --minimum-length "$TRUCVAL" -g "$FWD" -G "$REV" -o "$DIR"/processed-reads/sense/trimmed-R1/"$i".R1.fastq.gz -p "$DIR"/processed-reads/sense/trimmed-R2/"$i".R2.fastq.gz --discard-untrimmed "$DIR"/processed-reads/sense/dmplx/"$i".R1.fastq.gz "$DIR"/processed-reads/sense/dmplx/"$i".R2.fastq.gz >> "$DIR"/logs/cutadapt.sense.trimming.log &
done
echo "Trimming done!"

# trim ANTISENSE with cutadapt
antisense="$(ls "$DIR"/processed-reads/antisense/dmplx/*.fastq.gz | sed --expression='s/\.R.\.fastq\.gz//g' | sed --expression='s/unknown//g' | sed --expression="s/$DIR\/processed-reads\/antisense\/dmplx\///g" | uniq)"
# check
for f in $antisense; do echo "$f"; done
# now run
for i in $antisense; do
cutadapt -n 5 --error-rate 0.15 --minimum-length "$TRUCVAL" -g "$REV" -G "$FWD" -o "$DIR"/processed-reads/antisense/trimmed-R1/"$i".R1.fastq.gz -p "$DIR"/processed-reads/antisense/trimmed-R2/"$i".R2.fastq.gz --discard-untrimmed "$DIR"/processed-reads/antisense/dmplx/"$i".R1.fastq.gz "$DIR"/processed-reads/antisense/dmplx/"$i".R2.fastq.gz >> "$DIR"/logs/cutadapt.antisense.trimming.log &
done
echo "Trimming done!"


# check output for primers (should be none)
cat "$DIR"/processed-reads/sense/trimmed-R1/*.R1.fastq.gz | gzip -cd | sed -n '2~4p' | GREP_COLOR="1;44" grep --color=always -E "$FWDGR" | GREP_COLOR="1;41" grep --color=always -E "$REVGR" | head -n 400
cat "$DIR"/processed-reads/antisense/trimmed-R1/*.R1.fastq.gz | gzip -cd | sed -n '2~4p' | GREP_COLOR="1;44" grep --color=always -E "$FWDGR" | GREP_COLOR="1;41" grep --color=always -E "$REVGR" | head -n 400


#### get counts

# total pf
fqtools count "$DIR"/fastq/R1.fastq.gz

# with primer
expr `fqtools count "$DIR"/processed-reads/sense/R1.fastq.gz` + `fqtools count "$DIR"/processed-reads/antisense/R1.fastq.gz`

# with barcodes
expr `cat "$DIR"/processed-reads/sense/dmplx/*.R1.fastq.gz | fqtools count` + `cat "$DIR"/processed-reads/antisense/dmplx/*.R1.fastq.gz | fqtools count`

# after trim
expr `cat "$DIR"/processed-reads/sense/trimmed-R1/*.fastq.gz | fqtools count` + `cat "$DIR"/processed-reads/antisense/trimmed-R1/*.fastq.gz | fqtools count`


####################################
# NOW MOVE TO DADA2 IN R 
####################################
####################################

# after filter
expr `cat "$DIR"/processed-reads/sense/filtered-R1/*.fastq.gz | fqtools count` + `cat "$DIR"/processed-reads/antisense/filtered-R1/*.fastq.gz | fqtools count`
