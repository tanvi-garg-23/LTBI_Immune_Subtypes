#!/bin/bash
#SBATCH --time=34:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
##SBATCH --array=1-200
#SBATCH --job-name=gzip_GSE79362
#SBATCH --output=log/gzip.o
#SBATCH --error=log/gzip.e

#module load sratoolkit/2.10.9


#!/bin/bash

# Loop through each file in the current directory
for file in *.fastq; do
    # Check if the file is a regular file
    if [ -f "$file" ]; then
        # Compress the file using gzip
        gzip "$file"
    fi
done

echo "Compression complete."

