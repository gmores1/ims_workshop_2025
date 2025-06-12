#!/bin/bash
#SBATCH --job-name="bacteria_assembly"
#SBATCH -A e32825
#SBATCH -p short
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH -o 
#SBATCH -e 

#open the proper conda environment
conda activate /projects/e32825/condaenvs/short_assembly_env

#loop through each file available in the desired folder and create a variable "sample" that will be named based on filename
for i in projects/e32825/data/bacteria_assembly/on_your_own_reads
  sample=$(basename -s _{1,2}.fastq.gz $i)

  #output sample name
  echo ${sample}
  
  #check read counts
  perl /projects/e32825/Scripts/fastq_stats.pl \
    -p /projects/e32825/data/bacteria_assembly/on_your_own/${sample} \
    > ${sample}.read_stats.untrimmed.txt
  
  #downsampling - reduce read coverage to 60-100x
  perl /projects/e32825/Scripts/downsample_paired.pl \
    -p ${sample}_ds \
    -g 1790000 \
    -f 40 \
    -o \
    /projects/e32825/data/bacteria_assembly/on_your_own/sample_1.fastq.gz \
    /projects/e32825/data/bacteria_assembly/on_your_own/sample_2.fastq.gz 
  
  export KRAKEN_DEFAULT_DB="/projects/e32825/Databases/kraken_db/minikraken_20171019_8GB"
  kraken \
      --preload --threads 1 --fastq-input --gzip-compressed \
      --paired ${sample}_ds_1.fastq.gz ${sample}_ds_2.fastq.gz \
      | kraken-report - > ${sample}.kraken_report.txt
  
  perl /projects/e32825/Scripts/kraken-topspecies.pl \
    ${sample}.kraken_report.txt ${sample} \
    > ${sample}.kraken_topspecies.txt
  
  fastp \
      --in1 ${sample}_ds_1.fastq.gz --in2 ${sample}_ds_2.fastq.gz \
      --out1 ${sample}_trimmed_paired_1.fastq.gz --unpaired1 ${sample}_trimmed_unpaired_1.fastq.gz \
      --out2 ${sample}_trimmed_paired_2.fastq.gz --unpaired2 ${sample}_trimmed_unpaired_2.fastq.gz \
      -h ${sample}_fastp.html -j ${sample}_fastp.json \
      -w 1
  
  perl /projects/e32825/Scripts/fastq_stats.pl -p ${sample}_trimmed_paired > ${sample}.read_stats.trimmed.txt
  
  done
