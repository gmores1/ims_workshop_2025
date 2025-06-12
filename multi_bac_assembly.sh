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

#loop through each file available in the desired folder and create a variable "sample" that will be named based on file name
for i in projects/e32825/data/bacteria_assembly/on_your_own_reads
  sample=$(basename -s _{1,2}.fastq.gz $i)

  #output sample name
  echo ${sample}
  
  #check read counts
  perl /projects/e32825/Scripts/fastq_stats.pl \
    #path and prefix for the reads
    -p /projects/e32825/data/bacteria_assembly/on_your_own/${sample} \
    > ${sample}.read_stats.untrimmed.txt
  
  #downsampling - reduce read coverage to 60-100x
  perl /projects/e32825/Scripts/downsample_paired.pl \
    #output file prefix
    -p ${sample}_ds \
    #estimated genome size, in bp
    -g 1790000 \
    #fold coverage to output
    -f 40 \
    #guarantees output, even if the input read set is less than the value given to f
    -o \
    /projects/e32825/data/bacteria_assembly/on_your_own/${sample}_1.fastq.gz \
    /projects/e32825/data/bacteria_assembly/on_your_own/${sample}_2.fastq.gz 

  #species estimation using Kraken
  export KRAKEN_DEFAULT_DB="/projects/e32825/Databases/kraken_db/minikraken_20171019_8GB"
  kraken \
      --preload --threads 1 --fastq-input --gzip-compressed \
      --paired ${sample}_ds_1.fastq.gz ${sample}_ds_2.fastq.gz \
      | kraken-report - > ${sample}.kraken_report.txt
  #distill kraken report down to the top species-level assignment and percentage of reads assigned to the species
  perl /projects/e32825/Scripts/kraken-topspecies.pl \
    ${sample}.kraken_report.txt ${sample} \
    > ${sample}.kraken_topspecies.txt
  #output the top species
  echo ${sample}.kraken_topspecies.txt

  #trim reads by removing Illumina adapter sequences that may be included in the reads due to short library fragments
  fastp \
      #input read files
      --in1 ${sample}_ds_1.fastq.gz --in2 ${sample}_ds_2.fastq.gz \
      #paired output read files and singleton read files if only one of the paired reads passed filters
      --out1 ${sample}_trimmed_paired_1.fastq.gz --unpaired1 ${sample}_trimmed_unpaired_1.fastq.gz \
      --out2 ${sample}_trimmed_paired_2.fastq.gz --unpaired2 ${sample}_trimmed_unpaired_2.fastq.gz \
      #filtering and trimming report in html format
      -h ${sample}_fastp.html -j ${sample}_fastp.json \
      #number of parallel threads to use
      -w 1

  #count paired reads
  perl /projects/e32825/Scripts/fastq_stats.pl -p ${sample}_trimmed_paired > ${sample}.read_stats.trimmed.txt

  #generate de novo whole assembly from trimmed reads using SPAdes
  spades.py \
    -o ${sample}.spades \
    -1 ${sample}_trimmed_paired_1.fastq.gz \
    -2 ${sample}_trimmed_paired_2.fastq.gz \
    --isolate \
    -t 12 -cov-cutoff auto \
      tee ${sample}_spades_output.txt

  sbatch spades.sh
  
  done
