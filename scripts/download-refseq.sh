#!/usr/bin/env/ sh

#################################################### 
#################################################### 


# DOWNLOAD REFSEQ
mkdir temp/reference-library
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER -P temp/reference-library
head temp/reference-library/RELEASE_NUMBER
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/README -P temp/reference-library
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz -P temp/reference-library
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.1.genomic.fna.gz -P temp/reference-library
# check the release number
file="temp/reference-library/RELEASE_NUMBER"
gbv=$(cat "$file")
echo "$gbv"
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/RefSeq-release"$gbv".catalog.gz -P temp/reference-library
gzip -cd temp/reference-library/RefSeq-release"$gbv".catalog.gz | grep "mitochondrion" > temp/reference-library/refseq.mitochondrion.cat.tsv
rm temp/reference-library/RefSeq-release"$gbv".catalog.gz

# unzip
gzip -d temp/reference-library/mitochondrion.1.1.genomic.fna.gz 
gzip -d temp/reference-library/mitochondrion.2.1.genomic.fna.gz 

# join
cat temp/reference-library/mitochondrion.1.1.genomic.fna temp/reference-library/mitochondrion.2.1.genomic.fna > temp/reference-library/refseq.mitochondrion.genomic.fna
rm temp/reference-library/mitochondrion.1.1.genomic.fna temp/reference-library/mitochondrion.2.1.genomic.fna

# clean
sed -i -e 's/ .*//g' temp/reference-library/refseq.mitochondrion.genomic.fna

# check 
head temp/reference-library/refseq.mitochondrion.genomic.fna
