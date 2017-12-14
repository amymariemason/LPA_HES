#!/bin/bash
#SBATCH --time=8:0:0 
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH -n 1
#SBATCH -c 1 
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mem_bind=verbose,local
#SBATCH --output=qcrun.out
#SBATCH --error=qcrun.err
#SBATCH --job-name=13_12extract
#SBATCH --mem=64G

. /etc/profile.d/modules.sh
module load gcc/5.2.0
module load qctool2

qctool -g /scratch/curated_genetic_data/uk_biobank/imputed/full_release/HRC_subset/ukb_imp_chr6_HRConly.bgen \
-og /scratch/am2609/Gwas/lpa_snps.bgen -incl-positions /scratch/am2609/Gwas/lpa_vars.txt
