#!/bin/bash
#SBATCH --time=8:0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mem_bind=verbose,local
#SBATCH --output=GWAS1.out
#SBATCH --error=GWAS1.err
#SBATCH --job-name=13_12extract
#SBATCH --mem=64G

. /etc/profile.d/modules.sh
module load snptest

module load gcc/5.2.0
module load qctool2

qctool -g /scratch/curated_genetic_data/uk_biobank/imputed/full_release/HRC_subset/ukb_imp_chr6_HRConly.bgen \
-og /scratch/am2609/Gwas/lpa_snps.bgen -incl-positions /scratch/am2609/Gwas/lpa_vars.txt


for output in isch_comb pad_comb sah_comb ich_comb haem_comb
do
echo $output
snptest -data /scratch/am2609/Gwas/lpa_snps.bgen /scratch/am2609/Gwas/lpa_HESoutcomes.sample \
-o /scratch/am2609/Gwas/GWAS_HESoutcomes_$output-output.out \
-exclude_samples /scratch/am2609/Gwas/exclusion_list.txt \
-frequentist 1 \
-method score \
-pheno $output \
-cov_names PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
done




