######
# Author: Amy Mason ; based heavily on code by Stephen Burgess
# Date: Oct 2017
# Goal: Perform GWAS on biobank stroke outcomes 
# Inputs: "AdiposityID_sampleset" sample file, empty of outcomes  - in same order as 500k HRConly_EurQCp .bgen file
# 
# Outputs: sample file with outcomes from stroke_pad_all_comb
######


# libraries and setup
rm(list=ls())
setwd("C://Users/am2609/Programs/GWAS")



# load adiposity ID sample list

sample_all = read.table("AdiposityID_sampleset")


# load outcomes file and list outcomes of interest
outcomes = read.csv("./HES outcomes/stroke_pad_all_comb.csv")
wanted_outcomes<-c("isch_comb","pad_comb","sah_comb","ich_comb","haem_comb")
#########

# create function that returns column of binary outcomes for whole 
#empty results files
outcome_output<-rep(NA, nrow(sample_all))

extract_outcome<-function(sample_list, sample_id="ID_1", outcomefile, outcome, outcome_id, output){
  events<-which(outcomefile[,outcome]==1)
  cases = unique(outcomefile[events,outcome_id])
  output_temp = ifelse(as.numeric(sample_list[2:nrow(sample_list),sample_id])%in%cases, 1, 0)
  return(output_temp)
}

diag="I26"

events = which(startsWith(as.character(outcomes[,4]),diag))
cases = unique(outcomes[events,1])

diagcode_1 = rep(0, dim(sampleEUunrel)[1])
diagcode_1 = ifelse(as.numeric(sampleEUunrel_pheno[2:dim(sampleEUunrel)[1],1])%in%cases, 1, 0)
all_diagcode = c("B", diagcode_1)
sink("logfile_diagtotals.log", append=TRUE, split=TRUE)
cat("\nThis is number of people with outcome starting", diag, "\n" )
table(all_diagcode)
sink()

sample_all = cbind(sampleEUunrel_pheno[,1:3], pcs_2, all_diagcode)

write.table(sample_all, paste(diag,"lpa.sample",sep="."), row.names=FALSE, quote=FALSE)



# sense check - is the sample file same length as original file


# for extracting variants - Savita did this already
#
# module load qctool
# qctool -g /scratch/curated_genetic_data/uk_biobank/imputed/interim_release/chr6impv1.bgen \
#       -og /scratch/sb452/lpa/chr6impv_lpa.bgen -incl-positions /scratch/sb452/lpa/lpa_vars_ukbb.txt

# for running GWAS on cardio

module load snptest
snptest -data /scratch/am2609/GWAS/ukb_imp_chr6LPA_HRConly_EurQCp.bgen /scratch/am2609/GWAS/lpa_I739.sample -o /scratch/am2609/GWAS/ex.out -frequentist 1 -method score -pheno all_diagcode


###########
# 150k interim release

#pcs2 = read.table("ukbb_inclusion_list_50PCs_pcs.txt") # 150k
#sample = read.table("impv1.sample", stringsAsFactors=FALSE, header=TRUE) # 150k
#sample_case = rep(0, 152249)
#sample_case = ifelse(as.numeric(sample[2:152250,1])%in%cases, 1, 0)
#I739 = c("B", sample_case)


# read in snp output

try<- read.delim("exI60.out", header=TRUE, allowEscapes=FALSE, sep=" ",  quote="", na.strings="", comment.char="#")

