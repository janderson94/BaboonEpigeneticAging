#!/bin/bash
#SBATCH --get-user-env

# get all data for CHROMNAME

touch all_mratios_CHROMNAME.txt

for f in `cat ids.txt`; do awk 'BEGIN {OFS="\t"} {print $0,FILENAME}' $f | grep -P 'CHROMNAME\t' >> all_mratios_CHROMNAME.txt; done

# get sites covered in at least 50% of the dataset
# replace 35 with 1/2 number of individuals in dataset (or whatever)

awk '{print $1"_"$2}' all_mratios_CHROMNAME.txt | sort | uniq -c > temp_CHROMNAME.txt
awk '$1 > 250' temp_CHROMNAME.txt | awk '{print $2}' > all_mratios_50_CHROMNAME.txt
rm temp_CHROMNAME.txt

sed -e s/_/'\t'/g all_mratios_50_CHROMNAME.txt | awk '{OFS="\t"; print $1,$2,$2}' > all_mratios_50_CHROMNAME.bed

awk '{OFS="\t"; print $1,$2,$2,$5,$6,$7,$8,$13}' all_mratios_CHROMNAME.txt > all_mratios_CHROMNAME.bed
rm all_mratios_CHROMNAME.txt

# raw data for sites covered in at least 50% of the dataset

intersectBed -a all_mratios_CHROMNAME.bed -b all_mratios_50_CHROMNAME.bed -u > all_mratios_CHROMNAME_v2.txt

