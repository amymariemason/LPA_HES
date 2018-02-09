##################################################################
# Determine best LP(A) snps SCORE with snps present in UKBB and CardiogramC4D
# Author: Amy Mason
# Date: Jan 2018
##################################################################
# Only relevant to LP(a) project, will not generalise?

#################################################################
# programs
#################################################################

######################################
# setup desktop
######################################
# libraries
require(data.table)

create_blank_file<- function(inputfile){
  
  # create empty data frame for variant data; includes duplicate reporting
  output<-fread(file=inputfile,nrows=1)
  output$duplicate<-0
  output<- output[0,]
  
  return(output)
}

create_missing_report<-function(varlist, varname){
  missinglist<-as.data.frame(varlist[,varname])
  names(missinglist)<-varname
  missinglist$missing<-rep(0, nrow(missinglist)) 
  return(missinglist)
}


# imports data into blank file, reports missing data in list
import_data<-function(inputfile, inputname, varlist, varname=names(varlist)[[1]], maxdup=5){
  #inputfile = file to subset
  #inputname = variable in input file to search for variable
  #varlist = file containing variables to search for
  # (varname) = variable to search; defaults to first column of varlist
  # (maxdup) = maximum number of duplications to search for; default = 5
  
  # create blank report files
  missingreport<-create_missing_report(varlist,varname)
  output<-create_blank_file(inputfile)  
  num_er=0 # number of errors reported
  
  # loop to read in data from outside file; assigned as missing added to output file
  num_er =0 # error counter
  for (j in varlist[,varname]){
    # to watch progress
    print(j)
    # clear test dummy var
    #  if (exists("test")){rm(test)}
    # test if input exists for this value; reports error if missing
    test <- tryCatch(
      # this is what I want it to do:  
      fread(file=inputfile,nrows=maxdup, skip=as.character(j))
      ,
      # if error occurs
      error=function(error_message) {
        message(paste0("Error AT VALUE: ", j))
        message(error_message)
        return("ERROR")
      }
    )
    # if error, report variant as not found
    if(is.character(test)){
      missingreport[grep(j, missingreport[, varname]),"missing"]<-1;
      num_er <- num_er +1;
    }
    # if no error, add collected lines to file, reporting duplicate matches
    if(!(is.character(test))) {
      test$duplicate<-rep(0, nrow(test));
      names(test)<-names(output);
      test<-as.data.frame(test)
      testSubset <- test[grep(j, test[,inputname]), ];
      #add warning for multiple matching rows			
      if(nrow(testSubset)!=1) testSubset$duplicate<-rep(1, nrow(testSubset));
      output <- rbind(output, testSubset)
    }
  }  
  outlist<-list(output, num_er, missingreport)
  return(outlist)
}

#########################################################################
# Look up snps that are in 935 lp(a) list, c4d list and UK biobank
#########################################################################

#load data
CD4<-"//me-filer1/home$/am2609/My Documents/mi.add.030315.website.txt"
statsfile<-"//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress/Inputs/ukb_imp_chr6_HRCvars_EURsamples_snpstats.txt"
lpa <- read.table("//me-filer1/home$/am2609/My Documents/Blood Cell Traits Data/LPA/LPA_Variants_EUwinsor_withoutM.txt", sep="\t", header=T, colClasses="character")


# subset C4D data to that which overlaps with lpa

try<-import_data(CD4, inputname="bp_hg19", lpa, varname="pos") #to subset by cardio
keep<-as.data.frame(try[[1]])
message(nrow(keep), " variants in Cardiogram+C4D and lpa var list")

#subset again by what is in UKBB
try2<- import_data(statsfile, inputname="position", keep, varname="bp_hg19")  # to do UKBB and Cardio
keep2<-try2[[1]]
message(nrow(keep2), " variants in UKBioBank & Cardiogram+C4D and lpa var list")

# check for duplicates by position/ label with variant ID

lpa_label<-lpa[,c("variantID", "chr", "pos", "allele1", "allele2")] # these are all characters
keep2<- data.frame(lapply(keep2, as.character), stringsAsFactors=FALSE)
keepA<- merge(keep2, lpa_label, by.y=c("pos", "allele1", "allele2"), by.x=c("position","alleleA","alleleB"), all.x=TRUE)
keepA<-keepA[!is.na(keepA$variantID),]
keepB<- merge(keep2, lpa_label, by.y=c("pos", "allele1", "allele2"), by.x=c("position","alleleB","alleleA"), all.x=TRUE)
keepB<-keepB[!is.na(keepB$variantID),]
keep3<-rbind(keepA,keepB)
rm(keepA, keepB)

write.table(keep3,file="//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress/Outputs/LPA&UKBB&C4D_vars.txt", row.names = FALSE)
rm(keep, keep2, keep3)


########################################################################
# check imputation quality in UKBB & CD4 & LPA vars infor
########################################################################

Vars<-read.table("//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress/Outputs/LPA&UKBB&C4D_vars.txt", header=T)


# set acceptable thresholds for imputation quality
rare<-0.9
lowfreq<-0.9
common<-0.6

# print out lists of problem snps
message(paste0("these rare snps have impute values <", rare, " & Allele B frequencies <1%"))
Vars$threshold1<-ifelse(Vars$impute_info<rare & Vars$alleleB_frequency<0.01, 1,0)
Vars[Vars$threshold1==1, c("rsid","chromosome","position", "alleleA","alleleB")]

message(paste0("these low frequency snps have impute values <", lowfreq, " & Allele B frequencies 1%<= x <5%"))
Vars$threshold2<-ifelse(Vars$impute_info<rare & Vars$alleleB_frequency>=0.01 & Vars$alleleB_frequency<0.05 , 1,0)
Vars[Vars$threshold2==1,  c("rsid","chromosome","position", "alleleA","alleleB")]

message(paste0("these common snps have impute values <", common,"  & Allele B frequencies >5%"))
Vars$threshold3<-ifelse(Vars$impute_info<common & Vars$alleleB_frequency>=0.05, 1,0)
Vars[Vars$threshold3==1,  c("rsid","chromosome","position", "alleleA","alleleB")]

Vars$QCfail<- Vars$threshold3 + Vars$threshold2 +Vars$threshold1

message(nrow(Vars[Vars$QCfail==1,]), " variants have insufficiently high imputation score (removed)")

write.table(Vars[Vars$QCfail==0,],file="//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress/Outputs/QCed_LPA_BB_C4D_vars.txt", row.names = FALSE)


# 616 variants remain

########################################################################
# Do stepwise regression: setup
########################################################################

Vars_ID<-Vars[Vars$QCfail==0, c("variantID", "chromosome", "position", "alleleA", "alleleB")]

# take LP(a) data

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

# identify which columns relate to high quality imputation snps in UKbiobank
whichcols<-which(names(data) %in% Vars_ID$variantID)
# check found all variants
assertthat::assert_that(length(whichcols)==nrow(Vars_ID))

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


#


#####################################################
##### Forwards stepwise functions #####
#####################################################

fstepreg <- function(y=y, x=x, c=c, p_ent=5e-8){
  snps <- NULL
  p <- 0 
  k<- 1
  while(p<p_ent){
    steps <- top_snp(y=y,x=x,c=c)
    p <- steps$p
    if(p<p_ent){
      c <- cbind(c, x[,(names(x)==steps$snp)])
      names(c)[names(c)=="x[, (names(x) == steps$snp)]"] <- steps$snp
      x <- x[,!(names(x)==steps$snp)]
      snps <- c(snps, steps$snp)
      cat(steps$snp,"with p =",steps$p,"selected in step",k,"\n")
      k <- k + 1
    }
  }
  model <- lm(y ~ ., data=c)
  return(list(model=model, snps=snps))
}

fstepwisereg <- function(y=y, x=x, c=c, p_ent=5e-8, p_exit=1e-5, r2=0.2){
  snps <- NULL
  p <- 0 
  pb <- 1
  k<- 1
  nc <- ncol(c)+1
  nvars <- length(lm(y~.,data=c)$coef)+1
  while(p<p_ent){
    if(k==1){x1 <- x}else{x1 <- ld(x=x, c=c, nc=nc, r2=r2)}
    steps <- top_snp(y=y,x=x1,c=c)
    p <- steps$p
    if(p<p_ent){
      c <- cbind(c, x[,(names(x)==steps$snp)])
      names(c)[names(c)=="x[, (names(x) == steps$snp)]"] <- steps$snp
      x <- x[,!(names(x)==steps$snp)]
      snps <- c(snps, steps$snp)
      cat(steps$snp,"with p =",steps$p,"selected in step",k,"\n")
    }
    while(pb>p_exit){
      bsteps <- bottom_snp(y=y,c=c,nvars=nvars)
      pb <- bsteps$p
      if(pb>p_exit){
        x <- cbind(x, c[,(names(c)==bsteps$snp)])
        names(x)[names(x)=="c[, (names(c) == bsteps$snp)]"] <- bsteps$snp
        c <- c[,!(names(c)==bsteps$snp)]
        snps <- snps[snps!=bsteps$snp]
        cat(bsteps$snp,"with p =",bsteps$p,"removed in step",k,"\n")
      }
    } 
    k <- k + 1
    pb <- 1
  }
  model <- lm(y ~ ., data=c)
  p_overall <- summary(model)$coefficients[,4]
  p_overall <- p_overall[nvars:length(p_overall)]
  return(list(model=model, snps=snps, p=p_overall))
}

top_snp <- function(y=y, x=x, c=c){
  snp <- ""
  z <- 0
  p <- 1
  for(i in 1:ncol(x)){
    mod <- lm(y ~ x[,i] + ., data=c)
    z_x <- abs(summary(mod)$coefficients[2,3])
    p_x <- summary(mod)$coefficients[2,4]
    snp_x <- names(x)[i]
    if(z_x>z){z <- z_x; p <- p_x; snp <- snp_x}
    if((i/100)%%1==0){cat(i,"--DONE-- ")}
  }
  return(list(snp=snp, p=p))
}

bottom_snp <- function(y=y, x=x, c=c, nvars=nvars){
  mod <- lm(y ~ ., data=c)
  z <- abs(summary(mod)$coefficients[,3])
  z <- z[nvars:length(z)]
  zi <- 1:length(z)
  zi <- zi[z==min(z)][1]
  snp <- names(z)[z==min(z)][1]
  p <- 2*pnorm(-min(z))
  return(list(snp=snp, p=p))
}

ld <- function(x=x, c=c, nc=nc, r2=r2){
  c1 <- as.matrix(c[,nc:ncol(c)])
  x1 <- x
  for(j in 1:ncol(c1)){ld <- (cor(c1[,j],x1,use="pairwise.complete.obs"))^2; x1 <- x1[,ld<r2]}
  return(x1) 
}

##########################################################################
# run functions
###########################################################################

# find variant columns in data_sub

# identify which columns relate to high quality imputation snps in UKbiobank
whichcols2<-which(names(data_sub) %in% Vars_ID$variantID)

# turn those columns into a numeric matrix
g <- as.matrix(data_sub[,whichcols2])
class(g) <- "numeric"
g <- as.data.frame(g)
miss <- sapply(g, function(x)sum(is.na(x)))
mono <- sapply(g, function(x)length(unique(x))==1)
g <- g[,!mono]

y <- data_sub$ln_lpa
x <- g
c <- data_sub[, c("ages", "sex", "cohort", "assay_method2", "assay_method3", "gpc_1", "gpc_2", "gpc_3", "gpc_4", "gpc_5")]
c$cohort <- factor(c$cohort)
fstep <- fstepwisereg(y=y, x=x, c=c, r2=0.4)

# save selected snps
snps_selected <- fstep$snps
p <- fstep$p
results <- data.frame(snps=snps_selected, p=p) 


write.table(results, "//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress/Outputs/fstep_UKBB_C4D_lpa.txt", row.names=F, quote=F, sep="\t")


