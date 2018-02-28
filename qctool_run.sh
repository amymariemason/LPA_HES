#!/bin/bash
#SBATCH --time=8:0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mem_bind=verbose,local
#SBATCH --error=qctool_%a.err
#SBATCH --mem=64G
#SBATCH --array=1-22
#SBATCH --job-name=qct_%a
#SBATCH --output=qctool_%a.out


. /etc/profile.d/modules.sh      
module load default-cardio       
module load slurm
module load use.own

module load gcc/5.2.0 
module load qctool2

date

chr=${SLURM_ARRAY_TASK_ID}

qctool -g /scratch/curated_genetic_data/uk_biobank/imputed/full_release/HRC_subset/ukb_imp_chr${chr}_HRConly.bgen \
-og /scratch/am2609/Gwas/lpa_snps_chr${chr}.bgen -incl-positions /scratch/am2609/Gwas/test_snps.txt


date




