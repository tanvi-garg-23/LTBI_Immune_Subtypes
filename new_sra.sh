#!/bin/bash
#SBATCH --time=25:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --array=1-200
#SBATCH --job-name=fastq_dump_GSE107993
#SBATCH --output=log/fastq_dump_test_%a.o
#SBATCH --error=log/fastq_dump_test_%a.e

export PATH=$PATH:/work/gr-fe/garg/sratoolkit.3.0.0-ubuntu64/bin 


sampleinfo="GSE79362.txt"
file_name=`sed -n "$SLURM_ARRAY_TASK_ID"p $sampleinfo |  awk '{print $1}'`
echo "file_name" $file_name

# for paired end reads only
#fastq-dump --split-3  $file_name
prefetch $file_name \\
fasterq-dump $file_name
# zip all resulting read files
#gzip *.fastq

