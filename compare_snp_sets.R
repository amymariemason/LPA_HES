##################################################################
# Compare best LP(A) snps from 3 sets generated
# Author: Amy Mason
# Date: Jan 2018
##################################################################
# Only relevant to LP(a) project, will not generalise?
# as all sets extremely similar use the overlap with cardiogram C4D set 

##############################################################
# Load Data

# sets of snps chosen from:


## all lpa data snps
fsteps <- read.table("//me-filer1/home$/am2609/My Documents/Blood Cell Traits Data/LPA/fstep_snps_0.4_EUwinsor.txt", sep="\t", header=T, colClasses="character")


## UKBB snps

UKBB<-read.table("//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress/Outputs/fstep_UKBB_LPA_snps.txt", header=TRUE,  colClasses="character")

## UKBB & C4D snps


UKBBC4D<-read.table("//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress/Outputs/fstep_UKBB_C4D_lpa.txt", header=TRUE,  colClasses="character")



#correlation table
data_rho_master <- read.table("//me-filer1/home$/am2609/My Documents/Blood Cell Traits Data/LPA/LPA_master_dataset_pcs_EUwinsor_withoutM.txt", header=T, sep="\t", colClasses="character")


############################ compare lpa set with UKBB set

# CHECK 1: how similar are the sets

## Step 1: which snps match

whichmatch<-which(fsteps$snps %in% UKBB$snps)
match<-fsteps[whichmatch,"snps"]
length(match)

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
heatmap.2(corr[colmatch,rowmatch],breaks=c(-1,-0.8,0.8,1), symkey=F, col=c("red","white","yellow"), dendrogram = "none", trace="none", scale="none")

# no cov>0.9 either way between remaining snps

#################### compare UKBB subset with LPA subset

## Step 1: which snps match

whichmatch<-which(UKBB$snps %in% UKBBC4D$snps)
match<-UKBB[whichmatch,"snps"]
length(match)

UKBBC4D_not_UKBB<-UKBBC4D[!(UKBBC4D$snps%in%match),]
UKBB_not_UKBBC4D<-UKBB[!(UKBB$snps%in%match),]

data_rho <- data_rho_master[(data_rho_master$prev_chd==0), ]
data_rho <- data_rho[, match(c(UKBBC4D_not_UKBB$snps, UKBB_not_UKBBC4D$snps ), names(data_rho))]
data_rho <- as.matrix(data_rho)
class(data_rho) <- "numeric"
corr <- cor(data_rho, use="complete.obs")

colmatch<-which(colnames(corr)%in%UKBBC4D_not_UKBB$snps)
rowmatch<-which(colnames(corr)%in%UKBB_not_UKBBC4D$snps)

# plot to visualise
library(gplots)
heatmap.2(corr[colmatch,rowmatch],breaks=c(-1,-0.8,0.8,1), symkey=F, col=c("red","white","yellow"), dendrogram = "none", trace="none", scale="none")

# no cov>0.9 either way between remaining snps

##########################compare UKBB with UKBB+c4D

## Step 1: which snps match

whichmatch<-which(fsteps$snps %in% UKBBC4D$snps)
match<-fsteps[whichmatch,"snps"]
length(match)

UKBBC4D_not_lpa<-UKBBC4D[!(UKBBC4D$snps%in%match),]
lpa_not_UKBBC4D<-fsteps[!(fsteps$snps%in%match),]

data_rho <- data_rho_master[(data_rho_master$prev_chd==0), ]
data_rho <- data_rho[, match(c(UKBBC4D_not_lpa$snps, lpa_not_UKBBC4D$snps ), names(data_rho))]
data_rho <- as.matrix(data_rho)
class(data_rho) <- "numeric"
corr <- cor(data_rho, use="complete.obs")

colmatch<-which(colnames(corr)%in%UKBBC4D_not_lpa$snps)
rowmatch<-which(colnames(corr)%in%lpa_not_UKBBC4D$snps)

# plot to visualise
library(gplots)
heatmap.2(corr[colmatch,rowmatch],breaks=c(-1,-0.8,0.8,1), symkey=F, col=c("red","white","yellow"), dendrogram = "none", trace="none", scale="none")

# no cov>0.9 either way between remaining snps


### Step 2: which gives the better score

# load data

dataset<-read.table( "//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress/Outputs/LPA_cov_set")
cov<-c("ln_lpa", "ages", "sex", "cohort", "assay_method2", "assay_method3", "gpc_1", "gpc_2", "gpc_3", "gpc_4", "gpc_5")
data_lpa<-dataset[,c(cov,fsteps$snps)]
SCORE_lpa<-lm(ln_lpa~.,data=data_lpa)
summary(SCORE_lpa)


data_UKBB<-dataset[,c(cov,UKBB$snps)]
SCORE_UKBB<-lm(ln_lpa~.,data=data_UKBB)
summary(SCORE_UKBB)


data_UKBBC4D<-dataset[,c(cov,UKBBC4D$snps)]
SCORE_UKBBC4D<-lm(ln_lpa~.,data=data_UKBBC4D)
summary(SCORE_UKBBC4D)

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
whichcols3<-which(names(data) %in% UKBBC4D$snps)
whichcols<-c(whichcols1, whichcols2, whichcols3)

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

####### 
# having chosen the intersection of c4d cardiogram and UKBB snps, output in a snptest suitable form 


UKBBC4D<-read.table("//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress/Outputs/fstep_UKBB_C4D_lpa.txt", header=TRUE,  colClasses="character")
lpa <- read.table("//me-filer1/home$/am2609/My Documents/Blood Cell Traits Data/LPA/LPA_Variants_EUwinsor_withoutM.txt", sep="\t", header=T, colClasses="character")
keep<-merge(UKBBC4D, lpa, by.x="snps", by.y="variantID", all.x=TRUE, all.y=FALSE)
keep$outputcol<-paste0("0",keep$chr,":",keep$pos)
write.table(t(keep$outputcol),"//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress/Outputs/lpa_snps_amy.txt",sep=" ",row.names=FALSE, quote =FALSE, col.names=FALSE)