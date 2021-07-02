#!/usr/bin/env sh

# make data data dir
mkdir -p temp/data

# cd
cd temp/data

# get the R1
tail -n +2 ../../assets/sequencing-master.csv | cut -d "," -f 9 | sort | uniq | xargs wget

# get the R2
tail -n +2 ../../assets/sequencing-master.csv | cut -d "," -f 12 | sort | uniq | xargs wget

# check the md5sums for R1
printf "checking md5sums"

# retrieve md5sums
tail -n +2 ../../assets/sequencing-master.csv | cut -d "," -f 7,8 | sort | uniq | sed 's/,/  /g' > checksums.md5
tail -n +2 ../../assets/sequencing-master.csv | cut -d "," -f 10,11 | sort | uniq | sed 's/,/  /g' >> checksums.md5

# check
md5sum -c checksums.md5
