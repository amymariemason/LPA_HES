######
# Author: Amy Mason ; based heavily on code by James Staley  
# Date: Oct 2017
# Goal: Perform Mendelian randomization  
# Inputs: list of variants, data to search for those variants
# Outputs: log file reporting missing variants, dataframe containing the subset of larger data matching the variants given 
######

#####################################################
##### Set-up #####
#####################################################
#Run for


rm(list=ls())


library(ggplot2)
library(MendelianRandomization)
library(plotly)
library(htmlwidgets)
library(assertthat)
library(stringr)


#set working directory
#setwd("C://Users/am2609/Dropbox (Personal)/lpamaster/")
setwd("//me-filer1/home$/am2609/My Documents/Blood Cell Traits Data/")

# load LPA data
data_rho_master <- read.table("LPA/LPA_master_dataset_pcs_EUwinsor_withoutM.txt", header=T, sep="\t", colClasses="character")

# load lpa effects

lpa <- read.table("LPA/LPA_Variants_EUwinsor_withoutM.txt", sep="\t", header=T, colClasses="character")
fsteps <- read.table("LPA/fstep_snps_0.4_EUwinsor.txt", sep="\t", header=T, colClasses="character")
lpa <- lpa[(lpa$variantID %in% fsteps$snp),]
lpa$chr.pos <- paste0("chr", lpa$chr, ":", lpa$pos)
lpa$snp <- lpa$chr.pos; lpa$a1 <- lpa$allele1; lpa$a2 <- lpa$allele2
lpa <- lpa[, c("variantID", "snp", "chr.pos", "chr", "pos" , "a1", "a2", "beta", "se")]
#Steve says: A1 is the effect allele


# setup testing loop
#1 import data
#2 setup mr comparision
#3 output report

setwd("//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress/500k")

# import list of traits to include
TL<-scan("./HESTraitlist.txt", what="", sep="\n")
TLframe<-as.data.frame(TL)
names(TLframe)<-"trait"
TLframe$traitmin<- substring(TLframe$trait,1,10)

# load names for traits 
traitnames<-read.csv("./Traitnames.csv", header=FALSE)
traitnames<-traitnames[,1:2]
names(traitnames)<-c("trait", "name")
TLall<-merge(TLframe, traitnames, by="trait", all.x=TRUE)




  


## External





# list of fields in output file
# : ref allele: alternate allele
#VARIANT (hg19) [CHROM:POS:REF:ALT]	Variant position as [Chromosome : hg19 Position : reference allele : alternate allele]
#rsid	SNP rsID as provided by UK Biobank
#nCompleteSamples	Number of samples analyzed with non-missing phenotypes
#AC	Dosage allele count from all samples with non-missing phenotypes
#ytx	Dosage alternate allele count (ytx means "y * x" where y=phenotype and x=alternate allele dosage) 
#beta	linear regression beta coefficient
#se	linear regression standard error
#pval	linear regression p-value

#steve says: reference allele is the effect allele