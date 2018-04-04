######
# Author: Amy Mason ; based heavily on code by James Staley  and Stephen Burgess
# Date: Oct 2017
# Goal: Import data for MR 
# Inputs: various LP(a) files, and correlation matrix files
# Outputs: log file reporting missing variants, dataframe containing the subset of larger data matching the variants given 
######

#####################################################
##### Set-up #####
#####################################################

#set working directory
# this needs to be the directory where the LP(a) files below are kept
setwd("//me-filer1/home$/am2609/My Documents/Programs/MR Projects/LPA/Inputs")
#setwd("D:/CurrentWork")


##################################
# Load LPA data
##################################

# correlation matrix
data_rho_master <- read.table("LPA_master_dataset_pcs_EUwinsor_withoutM.txt", header=T, sep="\t", colClasses="character")
# all snps associations with lpa
lpa <- read.table("LPA_Variants_EUwinsor_withoutM.txt", sep="\t", header=T, colClasses="character")
# list of variants
fsteps <- read.table("fstep_UKBB_C4D_lpa.txt", sep="\t", header=T, colClasses="character")


# subset associations to only required snps
lpa <- lpa[(lpa$variantID %in% fsteps$snp),]
lpa$chr.pos <- paste0("chr", lpa$chr, ":", lpa$pos)
lpa$snp <- lpa$chr.pos; lpa$a1 <- lpa$allele1; lpa$a2 <- lpa$allele2
X_associations <- lpa[, c("variantID", "snp", "chr.pos", "chr", "pos" , "a1", "a2", "beta", "se")]
#Steve says: A1 is the effect allele

# keep correlation matrix to only required snps 
  data_rho <- data_rho_master[(data_rho_master$prev_chd==0), ]
  data_rho <- data_rho[, match(X_associations$variantID, names(data_rho))]
  data_rho <- as.matrix(data_rho)
  class(data_rho) <- "numeric"
  corr <- cor(data_rho, use="complete.obs")
  
 