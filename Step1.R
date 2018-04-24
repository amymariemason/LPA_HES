######
# Author: Amy Mason ; based heavily on code by Stephen Burgess with additional comments from Jess
# Date: Oct 2017
# Goal: Perform GWAS on biobank stroke outcomes 
# Inputs: list of European unrelated files (QCed_Eur_unrelated),
# principle components files (Eur_QCp_PCs.txt),
# link adiposity to biobank file (sampleID_map)
# sample file of everyone in biobank (ukb_imp_genID.sample)
# Outputs: "AdiposityID_sampleset" sample file, empty of outcomes 
######

########################################################################
# Desktop setup
#######################################################################

# this should be the folder containing the biobank data or symbolic links to that data
setwd(BB_dir)

#######################################################################
# Load Biobank files 
######################################################################
# Note this is for the full biobank release

# principle components files
#Use the first 10 PCs. The file has the following headings: ‘FID, IID, PC1, PC2,…’. Where FID is the UKBiobank ID 
pcs  = read.table("Eur_QCp_PCs.txt", header=TRUE) 

# list of all samples  (2 columns = copies of IDS, missing =0, toIncl=0)
#List of all samples with their UK Biobank ID  (2 columns = copies of IDS, missing =0, toIncl=0). Note first row  
#consist of zeros. 
sampleHRC = read.table("./ukb_imp_genID.sample", stringsAsFactors=FALSE, header=TRUE) 

# this file shows who to include to remove related europeans
#Only want to include participants who are European, unrelated and passed QC. A list of participants to include is  
#contained the text file below.  
excludefiles2<-read.table("./QCed_Eur_unrelated.txt") 

# link adiposity to biobank
#This text file contains the link between the adiposity ID and genetic ID, the first column contains the genetic ID 
#and the second contains the adiposity ID 
samplelink  = read.table("sampleID_map.txt", stringsAsFactors=FALSE, header=TRUE) 


######################################################################
# Create blank sample file
######################################################################

# create list of the principle components
#Extract the first 10 PCs by creating an  list of the PCs by using the ID numbers in the PC doc and the sample doc 

whichlink2 = which(pcs[,1]%in%sampleHRC[,1])
pcs_1 = pcs[whichlink2,c(1,3:12)]

# add principle components to sampleHRC
#Merge the PCs with the sample data using the UK Biobank ID. By having the option all.x keep values from x that do not  
#match any value from y.
sampleHRC$id<-1:nrow(sampleHRC)
sampleHRC_PC<-merge(sampleHRC, pcs_1, by.x="ID_1", by.y="FID", all.x = TRUE)
sampleHRC_PC<-sampleHRC_PC[order(sampleHRC_PC$id),]


# add col describtors 
#We need to provide column describtors for QCTOOL and SNPTEST. For all of the PC columns replace the first entry with  
#a 'C' so it knows it is a continuous variable. Recall that the first row in sampleHRC consists of zeros 

sampleHRC_PC[1, c("PC1","PC2","PC3","PC4","PC5", "PC6", "PC7", "PC8", "PC9", "PC10")] <- rep("C", 10)

#creates a list of the adiposity sample IDS in the same sample order as the sampleHRC file
#Identifies which participants are in the link file and sampleHRC and then creates a list of the adiposity sample IDS  
#in the same sample order as the sampleHRC file 

whichlink = which(samplelink[,1]%in%sampleHRC[,1])
samplepheno = samplelink[whichlink,2]

# replace the UKbiobank IDs with Adiposity IDS
#sampleHRC contains two ID ('ID_1' and 'ID_2') fields which are identical. Replace the first one with the adiposity IDs 
sampleHRC_pheno<-merge(sampleHRC_PC,samplelink , by.x="ID_1", by.y="UKB_sample_ID", all.x = TRUE)
#return to original order
sampleHRC_pheno<-sampleHRC_pheno[order(sampleHRC_pheno$id),]
# replace ID_1
sampleHRC_pheno[,"ID_1"]<-sampleHRC_pheno$Adiposity_sample_ID
#remove unneeded columns
sampleHRC_pheno<-sampleHRC_pheno[, !names(sampleHRC_pheno) %in% c("Adiposity_sample_ID", "BP_sample_ID")]


# create exclusion list from inclusion inputmer
#At the moment we only have an inclusion list, but some of the genetic tools need an exclusion list.  

sampleHRC_PC$exclude<-!(sampleHRC_PC[,1]%in%excludefiles2[,1])
exclusionlist<-sampleHRC_PC[sampleHRC_PC$exclude==TRUE,]$ID_2
write.table(exclusionlist[2:nrow(exclusionlist)], paste0(output_dir,"exclusion_list.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)

# check that only excluded IDs lack PCs
#Check only the excluded participants do not have PC data. Use the assert function to see whether there are participants  
#that should be included but don't have PC values, and print message if this is the case.   

sampleHRC_PC$missing<-is.na(sampleHRC_PC$PC1)
sampleHRC_PC$error_check<- ifelse(sampleHRC_PC$exclude==FALSE&sampleHRC_PC$missing==TRUE,1,0)
assertthat::assert_that(nrow(sampleHRC_PC[sampleHRC_PC$error_check==1,])==0, msg= "PC data missing for some samples")
#ID values of participants that should be included in the analysis but do not have PC data 
write.table(sampleHRC_PC[sampleHRC_PC$error_check==1, 1], row.names=FALSE, quote=FALSE, col.names=FALSE)

# export blank sample file
write.table(sampleHRC_pheno, paste0(output_dir,"AdiposityID_sampleset"), row.names=FALSE, quote=FALSE)



# return desktop to working directory

setwd(home_dir)
