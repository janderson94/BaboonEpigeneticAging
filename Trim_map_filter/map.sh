#!/bin/bash

##################################################
#Coding for trimming and mapping RRBS data and calling methratios
##################################################
#Requires:
#TrimGalore (A functional version of Cutadapt and FastQC are required for TrimGalore).
#BSMAP
#Python
#Samtools


###################
#Trimming
###################
#Path to TrimGalore software
trim_galore_path='path_to_trim_galore'

#-a can specify adaptors, if known, otherwise adaptors will be automatically identified from known lists of possible adaptor sequences (e.g. standard illumina adaptors shown here)
#--length specifies minimum length for retaining reads
#--stringency dictates the required overlap between the end of the read and adaptor sequence for trimming
#--rrbs is required for RRBS data to be correctly trimmed of artificially introduced bases from library construction following msp1 digestion.
#--gzip if files are in .fastq.gz format
#FILE.fastq.gz or path to fastq file of interest
trim_galore_path/trim_galore -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --gzip --rrbs --length 15 --stringency 4 FILE.fastq.gz


###################
#Mapping
###################
#Path to BSMAP software
bsmap_path='path_to_bsmap'

#Path to trimmed read file
path_fastq=FILE_trimmed.fq.gz

#Path to genome of interest. Here, we used a combined Panu2 + Lambda Phage genome so reads could be mapped and bisulfite conversion rate could be estimated in one step.
path_genome=/path_to_genome/genome.fasta

#Path to desired sam output file
path_sam=FILE.sam

#-r 0 reports unique hits only
#-v sets allowed # of mismatches. If between 0 and 1, it is interpreted as a % of the read length (i.e. here 10% of read length is allowed to be mismatched)
bsmap_path/bsmap -a $path_fastq -d $path_genome -o $path_sam -v 0.1 -r 0


###################
#Calling methratios
###################
#Path to methratio.py script from BSMAP
path_to_methratio.py='path_to_methratio.py'

#Path to desired methratio output file
path_out=FILE.methratio.txt

#--context=CG will only return mratio calls at CpG dinucleotides, the vast majority of the locations of DNA methylation in baboons (and many primates, as far as we know)
#--combine-CpG will combine mratio call information from both strands
python path_to_methratio.py/methratio.py -o $path_out -d $path_genome --combine-CpG --context=CG $path_sam
