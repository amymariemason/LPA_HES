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
#SBATCH --mem=16G
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

# create an array of the covariates

IFS=' ' read -a covar_arr <<<"$covar"
echo "${covar_arr[@]}" 
# loop over outcomes to use

while read line
do
   printf "\n%s\n" "snptest running for outcome: ${line}"
snptest -data "${bgen_file}" "${sample_file}" \
-o "${output_dir}"/"${proj}_${line}.out" \
-exclude_samples /scratch/am2609/Gwas/input_all/exclusion_list.txt \
-frequentist 1 \
-method score \
-pheno "${line}" \
-cov_names "${covar_arr[@]}" 
done <"${outcome_file}"




#done



