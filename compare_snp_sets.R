##################################################################
# Compare best LP(A) snps
# Author: Amy Mason
# Date: Jan 2018
##################################################################
# Only relevant to LP(a) project, will not generalise?

##############################################################
# Load Data

# sets of snps chosen from:


## all lpa data snps
fsteps <- read.table("//me-filer1/home$/am2609/My Documents/Blood Cell Traits Data/LPA/fstep_snps_0.4_EUwinsor.txt", sep="\t", header=T, colClasses="character")


## UKBB snps

UKBB<-read.table("//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress/Outputs/fstep_UKBB.txt", header=TRUE,  colClasses="character")

#correlation table
data_rho_master <- read.table("//me-filer1/home$/am2609/My Documents/Blood Cell Traits Data/LPA/LPA_master_dataset_pcs_EUwinsor_withoutM.txt", header=T, sep="\t", colClasses="character")


# compare lpa set with UKBB set

# CHECK 1: how similar are the sets

## Step 1: which snps match

whichmatch<-which(fstep$snps %in% UKBB$snps)
match<-fsteps[whichmatch<-which(fsteps$snps %in% UKBB$snps),"snps"]

UKBB_not_lpa<-UKBB[!(UKBB$snps%in%match),]
lpa_not_UKBB<-fsteps[!(fsteps$snps%in%match),]

data_rho <- data_rho_master[(data_rho_master$prev_chd==0), ]
data_rho <- data_rho[, match(c(UKBB_not_lpa$snps, lpa_not_UKBB$snps ), names(data_rho))]
data_rho <- as.matrix(data_rho)
class(data_rho) <- "numeric"
corr <- cor(data_rho, use="complete.obs")

colmatch<-which(colnames(corr)%in%UKBB_not_lpa$snps)
rowmatch<-which(colnames(corr)%in%lpa_not_UKBB$snps)

# plot to visualise
library(gplots)
heatmap.2(corr[colmatch,rowmatch],breaks=c(-1,-0.9,0.9,1), col=c("red","white","yellow"), dendrogram = "none", trace="none", scale="none")

# no cov>0.9 either way between remaining snps

### Step 2: which gives the better score

# load data

dataset<-read.table( "//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress/Outputs/LPA_cov_set")
cov<-c("ln_lpa", "ages", "sex", "cohort", "assay_method2", "assay_method3", "gpc_1", "gpc_2", "gpc_3", "gpc_4", "gpc_5")
data_lpa<-dataset[,c(cov,fsteps$snps)]
SCORE_lpa<-lm(ln_lpa~.,data=data_lpa)

data_UKBB<-dataset[,c(cov,UKBB$snps)]
SCORE_UKBB<-lm(ln_lpa~.,data=data_UKBB)


##############
# subset LPA data to covariates only

require(data.table)
data <- fread("//me-filer1/home$/am2609/My Documents/Blood Cell Traits Data/LPA/LPA_master_dataset_pcs.txt", sep="\t", header=T, colClasses="character")
data <- data.frame(data)
# remove non-european sets
data1 <- data[data$study!="PROMIS" &data$study!="BRAVE",]
message(paste0(nrow(data)-nrow(data1), " records in PROMIS or BRAVE (removed as non-European)"))
data<-data1
rm(data1)
# remove morgam due to publication restrictions
data1 <- data[data$study!="MORGAM",]
message(paste0(nrow(data)-nrow(data1), " records in MORGAM (removed as not permitted to publish)"))
data<-data1
rm(data1)
# truncate lpa value
data$ln_lpa <- ifelse(as.numeric(data$lpa)>130, "130", data$lpa)

# filter missing data
data1 <- data[!is.na(data$chd_status),]
message(paste0(nrow(data)-nrow(data1), " records missing CHD status (removed)"))
data<-data1
rm(data1)
data1 <- data[!is.na(data$ln_lpa),]
message(paste0(nrow(data)-nrow(data1), " records missing lpa value (removed)"))
data<-data1
rm(data1)

# remove people who had CHD before start of study
data1 <- data[!(data$prev_chd==1),]
message(paste0(nrow(data)-nrow(data1), " participants CHD positive before start of study (removed)"))
data<-data1
rm(data1)

# identify which columns relate to snps in UKBBset, or LPA set
whichcols1<-which(names(data) %in% UKBB$snps)
whichcols2<-which(names(data) %in% fsteps$snps)
whichcols<-c(whichcols1, whichcols2)

# identify which covariates
#NOTE: ln_lpa is not a log transformed variable
cov<-c("ln_lpa", "ages", "sex", "cohort", "assay_method2", "assay_method3", "gpc_1", "gpc_2", "gpc_3", "gpc_4", "gpc_5")
whichcov<-which(names(data) %in% cov)

data_test<-data[,c(whichcov,whichcols)]

# subset to that data (1:82 <- covariates)
data_sub<-data[,c(whichcov, whichcols)]

# format data

data_sub$ln_lpa <- as.numeric(data_sub$ln_lpa)
data_sub$ages <- as.numeric(data_sub$ages)
data_sub$sex <- as.numeric(data_sub$sex)
data_sub$assay_method2 <- as.numeric(data_sub$assay_method2)
data_sub$assay_method3 <- as.numeric(data_sub$assay_method3)
data_sub$gpc_1 <- as.numeric(data_sub$gpc_1)
data_sub$gpc_2 <- as.numeric(data_sub$gpc_2)
data_sub$gpc_3 <- as.numeric(data_sub$gpc_3)
data_sub$gpc_4 <- as.numeric(data_sub$gpc_4)
data_sub$gpc_5 <- as.numeric(data_sub$gpc_5)

write.table(data_sub, file = "//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress/Outputs/LPA_cov_set")
