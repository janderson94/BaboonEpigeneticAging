#!/bin/bash
#SBATCH --get-user-env

module load TrimGalore

trim_galore -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --gzip --rrbs --length 15 --stringency 4 FILEINFO


# map with BSMAP
module load BSMAP
module load samtools

path_fastq=FILE_trimmed.fq.gz
path_genome=/data/tunglab/shared/genomes/pap_lam_genomes2.fasta
path_sam=FILE.sam

bsmap -a $path_fastq -d $path_genome -o $path_sam -v 0.1 -r 0

# get mratios

path_out=FILE.methratio.txt

python /data/tunglab/tpv/Programs/methratio.py -o $path_out -d $path_genome --combine-CpG --context=CG $path_sam
