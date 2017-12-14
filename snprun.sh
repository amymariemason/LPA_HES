#!/bin/bash
#SBATCH --time=8:0:0 
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH -n 1
#SBATCH -c 1 
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mem_bind=verbose,local
#SBATCH --output=snptest_h.out
#SBATCH --error=snptest_h.err
#SBATCH --job-name=13_12extract
#SBATCH --mem=64G

. /etc/profile.d/modules.sh
module load snptest

snptest -data /scratch/am2609/Gwas/lpa_snps.bgen /scratch/am2609/Gwas/lpa_HESoutcomes.sample \
-o /scratch/am2609/Gwas/GWAS_HESoutcomes_haem.out \
-exclude_samples /scratch/curated_genetic_data/uk_biobank/imputed/full_release/HRC_EUR_subset/toExcl_Imp_nonEur_or_QCf.txt \
-frequentist 1 \
-method score \
-pheno haem_comb \
-cov_names PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10


