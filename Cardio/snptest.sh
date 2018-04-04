#!/bin/bash
#SBATCH --time=8:0:0 
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mem_bind=verbose,local
#SBATCH --output=snptest_%A.out
#SBATCH --error=snptest_%A.err
#SBATCH --job-name=snptest
#SBATCH --mem=64G
#FILENAME: snptest_run
#AUTHOR : Amy Mason
#PURPOSE: run snptest on inputed bgen file using inputed .sample file for a set of outcomes
#INPUT: bgen_file .bgen file containing genetic information for each snp/ participant 
#       sample_file .sample file containing outcome information for each outcome/participant
#       outcome_file:  space seperated file containing each outcome for which to run program
#       output_dir, where to put the output files (this will be set to the output folder in your project file)
#OUTPUT: csv file of associations with snps for each outcome

. /etc/profile.d/modules.sh
module load snptest

# load list of what outcomes to use
value=$(<${outcome_file})
echo "$value"


#for output in isch_comb pad_comb sah_comb ich_comb haem_comb
#do
#echo $output
#echo snptest -data /scratch/am2609/Gwas/test2.bgen /scratch/am2609/Gwas/lpa_HESoutcomes.sample \
#-o /scratch/am2609/Gwas/GWAS_HESoutcomes_$output-output.out \
#-exclude_samples /scratch/curated_genetic_data/uk_biobank/imputed/full_release/HRC_EUR_subset/toExcl_Imp_nonEur_or_QCf.txt \
#-include_samples /scratch/curated_genetic_data/uk_biobank/genotyped/full_release/QCd_data/QCed_Eur_unrelated.txt
#-frequentist 1 \
#-method score \
#-pheno $output \
#-cov_names PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
#done



