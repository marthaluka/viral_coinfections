#!/bin/bash


#mapped_dir="/home3/2762681l/coinfections/data/mapped"
mapped_dir="/home3/2762681l/coinfections/data/mapped_trimmed"

cd $mapped_dir

# sam to bam
for file in "."/P*.sam; do

   #remove unnecessary extensions
   #base_name=${file%.fastq.gz.fastq.gz.sam}
   base_name=${file%._R1_001_trimmed.fq.gz.sam}

   #sam to bam
   samtools view -S -b $file > ${base_name}.bam

done


 # sort bam files 
for file2 in "."/*.bam; do
  
   samtools sort $file2 -o ${file2%.bam}_sorted.bam
  
done

# featureCounts
# #dont use this# featureCounts -O -t CDS -a ../ref_seq/hg38.refGene.gtf.gz -o ./${prefix}_count_table.txt *_sorted.bam
  # featureCounts -a ../ref_seq/hg38.refGene.gtf.gz -o count_table.txt *_sorted.bam

  
