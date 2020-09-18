#!/bin/bash

#This script is designed to filter intial methratio calls across all samples (from map.sh), retaining sites with at least 90% coverage (can easily be modified for different ranges of coverages depending on application)
#Run this script in parallel for each chromosome (CHROMNAME).
touch all_mratios_CHROMNAME.txt

#ids.txt is a list of sam file names from map.sh
for f in `cat ids.txt`; do awk 'BEGIN {OFS="\t"} {print $0,FILENAME}' $f | grep -P 'CHROMNAME\t' >> all_mratios_CHROMNAME.txt; done

# get sites covered in at least X% of the dataset
# Replace 250 with the desired number of samples with data per site for intitial filtering
awk '{print $1"_"$2}' all_mratios_CHROMNAME.txt | sort | uniq -c > temp_CHROMNAME.txt
awk '$1 > 250' temp_CHROMNAME.txt | awk '{print $2}' > all_mratios_90_CHROMNAME.txt
rm temp_CHROMNAME.txt

sed -e s/_/'\t'/g all_mratios_90_CHROMNAME.txt | awk '{OFS="\t"; print $1,$2,$2}' > all_mratios_90_CHROMNAME.bed

awk '{OFS="\t"; print $1,$2,$2,$5,$6,$7,$8,$13}' all_mratios_CHROMNAME.txt > all_mratios_CHROMNAME.bed
rm all_mratios_CHROMNAME.txt

# raw data for sites covered in at least 90% of the dataset
intersectBed -a all_mratios_CHROMNAME.bed -b all_mratios_90_CHROMNAME.bed -u > all_mratios_CHROMNAME_v2.txt

