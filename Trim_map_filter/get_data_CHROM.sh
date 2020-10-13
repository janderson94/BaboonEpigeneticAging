#!/bin/bash

##################################################################
#Filtering round 1
##################################################################
#Requires:
#Bedtools
#Awk

#This script is designed to filter initial methratio calls across all samples (from map.sh), retaining sites with at least 90% coverage. 
#It can easily be modified for different ranges of coverages depending on application.

#########################################
#Combining methratio calls across samples
#########################################
#Run this script in parallel for each chromosome (CHROMNAME).
touch all_mratios_CHROMNAME.txt

#ids.txt is the list of methratio file names output by methratio.py from BSMAP
#This iterates across methratio files and combines all CpG sites for a given chromosomes (CHROMNAME) into one file
for f in `cat ids.txt`; do awk 'BEGIN {OFS="\t"} {print $0,FILENAME}' $f | grep -P 'CHROMNAME\t' >> all_mratios_CHROMNAME.txt; done



#########################################
#Filtering based on coverage in dataset
#########################################
#From our list of CpG sites with methratio calls, we calculate the # of samples with call per site
awk '{print $1"_"$2}' all_mratios_CHROMNAME.txt | sort | uniq -c > temp_CHROMNAME.txt

# Then we retain only sites with calls in at least 90% of the dataset (250 samples /277 total samples)
# Replace 250 with the desired number of samples with data per site for initial filtering
awk '$1 > 250' temp_CHROMNAME.txt | awk '{print $2}' > all_mratios_90_CHROMNAME.txt
rm temp_CHROMNAME.txt

#Reformat the CpG sites of interest into bed format
sed -e s/_/'\t'/g all_mratios_90_CHROMNAME.txt | awk '{OFS="\t"; print $1,$2,$2}' > all_mratios_90_CHROMNAME.bed

#Reformat the list of full methratio calls at all sites
awk '{OFS="\t"; print $1,$2,$2,$5,$6,$7,$8,$13}' all_mratios_CHROMNAME.txt > all_mratios_CHROMNAME.bed
rm all_mratios_CHROMNAME.txt

#Subset the raw data for sites covered in at least 90% of the dataset
intersectBed -a all_mratios_CHROMNAME.bed -b all_mratios_90_CHROMNAME.bed -u > all_mratios_CHROMNAME_v2.txt

