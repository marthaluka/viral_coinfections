#!/bin/bash

# Set the input directory where reads are stored
#read_dir="/home1/illumina/Externally_Sequenced/Fastq/220528_NS500215_0149_AHMT52BGXL"
#read_dir="/home1/illumina/illumina-transferred/NextSeq500_backup/Fastq/220118_NB500959_0233_AHF2LMBGXK/fastq_files"
read_dir="/home3/2762681l/coinfections/data/trimmed"


# Set the output directory for the mapped reads
out_dir="/home3/2762681l/coinfections/data/mapped_trimmed"

# Set the reference genome fasta file
ref_genome="/home3/2762681l/coinfections/data/ref_seq/hg38_index"

# Make the output directory if it does not exist
if [ ! -d "$out_dir" ]; then
  mkdir "$out_dir"
fi

# Loop through all the fastq files in the input directory
#for read_file in "$read_dir"/P*.fastq.gz; do
for read_file in "$read_dir"/P*.fq.gz; do
  # Get the basename of the file without the path and extension
  #base=$(basename "$read_file".fastq.gz)
  base=$(basename "$read_file")
  # cleaned_basename=${base%_R1_001_trimmed_fq.gz}
  # echo $cleaned_basename
  # Map the reads to the reference genome
  hisat2 -x $ref_genome -U $read_file -S "$out_dir/$base.sam"
done
