#!/bin/bash
#SBATCH --time=0:5:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mem_bind=verbose,local
#SBATCH --mem=1G
#SBATCH --job-name=test
#SBATCH --output=test_%A.out
#SBATCH --error=test_%A.err
#Author: Amy Mason
#Filename: test
#Purpose: practise loops and variables in Unix

. /etc/profile.d/modules.sh      
module load default-cardio       
module load slurm
module load use.own

module load gcc/5.2.0 
module load qctool2

date


chr=${SLURM_ARRAY_TASK_ID}

printf "\n%s\n" "my inputs are in ${input_snps}"
printf "\n%s\n" "my output will go in ${output_dir}"
printf "\n%s\n" "my array number is ${chr}"

date

