######
# Author: Amy Mason ; based heavily on code by Stephen Burgess
# Date: Oct 2017
# Goal: Perform GWAS on biobank stroke outcomes 
# Inputs: list of European unrelated files (HRConly_EurQCp),
#principle components files, link adiposity to biobank file
# Outputs: "AdiposityID_sampleset" sample file, empty of outcomes 
######


rm(list=ls())
#setwd("C:/Dropbox/Dropbox Current/Dropbox/lpamaster/gwas/")
#setwd("C://Users/am2609/Dropbox (Personal)/lpamaster/gwas/")
#warninsetwd("C://Users/am2609/Programs/GWAS_inprogress")
setwd("//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress")


########
# 500k full release

# INPUT BIOBANK DATA

# principle components files
pcs  = read.table("Eur_QCp_PCs.txt", header=TRUE) 
# list of samples which are unrelated to each other (2 columns = copies of IDS, missing =0, toIncl=0)
sampleHRC = read.table("./500k/ukb_imp_genID.sample", stringsAsFactors=FALSE, header=TRUE) 
excludefiles<- read.table("./500k/toExcl_Imp_nonEur_or_QCf.txt", stringsAsFactors=FALSE, header=FALSE, sep="") 



# create ordered list of the principle components
whichlink2 = which(pcs[,1]%in%sampleHRC[,1])
pcs_1 = pcs[whichlink2,c(1,3:12)]

# add principle components to sampleHRC
sampleHRC_PC<-merge(sampleHRC[,1:3], pcs_1, by.x="ID_1", by.y="FID", all.x = TRUE)

# add col describtors 
sampleHRC_PC[1, c("PC1","PC2","PC3","PC4","PC5", "PC6", "PC7", "PC8", "PC9", "PC10")] <- rep("C", 10)


# link adiposity to biobank
samplelink  = read.table("sampleID_map.txt", stringsAsFactors=FALSE, header=TRUE) # link adiposity to biobank


#creates a list of the adiposity sample IDS
# in the same sample order as the sampleHRC file

whichlink = which(samplelink[,1]%in%sampleHRC[,1])
samplepheno = samplelink[whichlink,2]

# replace the UKbiobank IDs with Adiposity IDS
sampleHRC_pheno = sampleHRC_PC
sampleHRC_pheno[2:(dim(sampleHRC_PC)[1]), 1] = samplepheno


# check that only excluded IDs lack PCs
whichexclude=which(sampleHRC_PC[,1]%in%excludefiles[,1])

assertthat::assert_that(sum(is.na(sampleHRC_pheno[whichexclude,"PC1"]))==length(sampleHRC_pheno[whichexclude,"PC1"]),
                        msg= "Incorrect number of null values - check merge have worked correctly")

#

write.table(sampleHRC_pheno, "AdiposityID_sampleset", row.names=FALSE, quote=FALSE)
