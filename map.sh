#!/bin/bash
#SBATCH --get-user-env

#An updated version of TrimGalore is required.
#A functional version of Cutadapt and optionally FastQC are required for TrimGalore.
#module load TrimGalore
trim_galore_path='path_to_trimming_software/

#-a can specify adaptors, if known, otherwise adaptors will be automatically identified from known lists of possible adaptor sequences
#--length specifies minimum length for retaining reads
#--stringency dictates the required overlap between the end of the read and adaptor sequence for trimming
#--rrbs is required for RRBS data to be correctly trimming of artificially introduced bases from library construction following msp1 digestion.
#--gzip if files are in .fastq.gz format
#FILE.fastq.gz or path to file of interest
trim_galore_path/trim_galore -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --gzip --rrbs --length 15 --stringency 4 FILE.fastq.gz


# map with BSMAP
#module load BSMAP
#module load samtools
#BSMAP, pythoh, and samtools are required for mapping and calling mratios from trimmed reads

#Path to trimmed read file
path_fastq=FILE_trimmed.fq.gz

#Path to genome of interest. Here, we used a combined Panu2 + Lambda Phage genome so reads could be mapped and bisulfite conversion rate could be estimated in one step.
path_genome=/data/tunglab/shared/genomes/pap_lam_genomes2.fasta

#Path to desired sam output file
path_sam=FILE.sam

#-r 0 reports unique hits only
#-v sets allowed # of mismatches. If between 0 and 1, it is interpreted as a % of the read length (i.e. here 10% of read length is allowed to be mismatched)
bsmap -a $path_fastq -d $path_genome -o $path_sam -v 0.1 -r 0

# get mratios

path_out=FILE.methratio.txt
#Path to methratio.py script from BSMAP
path_to_methratio.py='path/to/methratio.py

#--context=CG will only return mratio calls at CpG dinucleotides, the vast majority of the locations of DNA methylation in baboons (and many primates, as far as we know)
#--combine-CpG will combine mratio call information from both strands

python path_to_methratio.py -o $path_out -d $path_genome --combine-CpG --context=CG $path_sam
