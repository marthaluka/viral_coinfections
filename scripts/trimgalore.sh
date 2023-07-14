#!/bin/bash

# Set the input directory where reads are stored
#read_dir="/home1/illumina/Externally_Sequenced/Fastq/220528_NS500215_0149_AHMT52BGXL"
read_dir="/home1/illumina/illumina-transferred/NextSeq500_backup/Fastq/220118_NB500959_0233_AHF2LMBGXK/fastq_files"

# Set the output directory for the trimmed reads
out_dir="/home3/2762681l/coinfections/data/trimmed"

# Make the output directory if it does not exist
if [ ! -d "$out_dir" ]; then
  mkdir "$out_dir"
fi


# Loop through all fastq files in the directory
for file in "$read_dir"/P*.fastq.gz
do
  # Run TrimGalore on the file
  trim_galore --length 40 --clip_R1 6 --three_prime_clip_R1 3 --output_dir $out_dir --quality 30 --illumina $file 
done
