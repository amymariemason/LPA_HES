#!/bin/bash
#SBATCH --time=8:0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mem_bind=verbose,local
#SBATCH --error=GWAS1.err
#SBATCH --mem=64G
#SBATCH --array=1-22
#SBATCH --job-name=qctool_%A_%a
#SBATCH --error=reports/qctool_%A_%a.err
#SBATCH --output=reports/qctool_%A_%a.out
#FILENAME: qctool_run
#AUTHOR : Amy Mason
#PURPOSE: create subset of the UKbiobank bgen files to enable running snptest
#INPUT: input_snps, a space seperated text file listing snps to extract. They must be in the form chr:pos with leading zeros if chr<10 ( eg 06:160493099)
#       output_name, what to call the output files
#       output_dir, where to put the output files (this will be set to the bgen folder in your project file)
#OUTPUT: one .bgen file for each chromosone, containing the extracted snps for that chr in .bgen format

. /etc/profile.d/modules.sh      
module load default-cardio       
module load slurm
module load use.own

module load gcc/5.2.0 
module load qctool2

date

chr=${SLURM_ARRAY_TASK_ID}


qctool -g "/scratch/curated_genetic_data/uk_biobank/imputed/full_release/HRC_subset/ukb_imp_chr${chr}_HRConly.bgen" \
-og "${output_dir}/${output_name}_${chr}.bgen" -incl-positions "${input_snps}"


date




