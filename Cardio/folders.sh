#!/bin/bash
#SBATCH --time=0:5:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mem_bind=verbose,local
#SBATCH --mem=1G
#SBATCH --job-name=folders
#SBATCH --error=folders_%A.err
#SBATCH --output=folders_%A.out
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

mkdir -p "${main_dir}/input_all"
mkdir -p "${main_dir}/code"
mkdir -p "${main_dir}"/projects/"${proj_name}"/{input,bgen,reports,output,errors}



printf "\n%s\n" "I am working on a new project called ${proj_name}, I want my output in the standard folders set-up which now exists at ${main_dir}/projects/${proj_name}"

cp /scratch/am2609/Gwas/Read_me "${main_dir}"/Folder_Structure_Read_me
cp /scratch/am2609/Gwas/code_current/*.sh  "${main_dir}"/code/