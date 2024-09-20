#!/bin/bash
#SBATCH --chdir /scratch/garg
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --qos=serial
#SBATCH --mem 16G
#SBATCH --time 24:00:00
#SBATCH --mail-user=tanvig@kth.se
#SBATCH --mail-type=FAIL,END
#SBATCH -o /scratch/garg/temp/mat_c/ssGSEA_output.out # Standard output
#SBATCH -e /scratch/garg/temp/mat_c/ssGSEA_error.err # Standard error
#SBATCH --job-name=ssGSEApipeline

####Input files needed: counts.rds file , filtered_gene_set.rds, script.R
set -x
RLIBRARY="/work/gr-fe/R_4.3.1"
COUNTSM="/work/gr-fe/garg/counts_data/data/ssgsea/counts_protein_exprs93.rds"
GENES="/work/gr-fe/garg/counts_data/data/ssgsea/filtered_gene_sets93.rds"
OUTPUT="/work/gr-fe/garg/counts_data/data/ssgsea/ssgsea_result_protein_exprs93.rds"
SCRIPT="/work/gr-fe/garg/counts_data/scripts/ssgsea_function.R"

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"

Rscript ${SCRIPT} ${RLIBRARY} ${COUNTSM} ${GENES} ${OUTPUT}

# print end date and echo total runtime
end=`date +%s`
runtime=$((end-start))
echo Runtime: $runtime seconds
