######
# Author: Amy Mason ; based heavily on code by Stephen Burgess
# Date: Oct 2017
# Goal: Perform GWAS on biobank stroke outcomes 
# Inputs: "AdiposityID_sampleset" sample file, empty of outcomes  - in same order as 500k HRConly_EurQCp .bgen file
# 
# Outputs: sample file with outcomes from cancer subsets; also adds ages and smoking
######


# libraries and setup
#setwd("C://Users/am2609/Programs/GWAS_inprogress")
#setwd("//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress/Outputs")

# this is the directory containing inputs/output/biobank input folders
home_dir<- "//me-filer1/home$/am2609/My Documents/Programs/MR Projects/LPA_test_run/"

# output file
output_dir<-paste0(home_dir, "Outputs/")

# load adiposity ID sample list

sample_all = read.table(paste0(output_dir,"AdiposityID_sampleset"), header = TRUE)


# load outcomes file and list outcomes of interest
outcomes = read.csv("//me-filer1/home$/am2609/My Documents/Stata_outcomes/cancer_subsets_v3.csv", header=TRUE)
wanted_outcomes<-c("ca_breast", "ca_prost", "ca_lung", "lung_squa", "lung_adeno",
                   "lung_NSC", "ca_bowel", "ca_mel", "ca_nhl", "NHL_DLBC", "NHL_fol", "ca_kid", "ca_headneck",
                   "head_oral", "head_lar", "ca_brain", "ca_blad", "ca_panc", "ca_uter", "ca_leuk", "leuk_aml",
                   "leuk_cll", "ca_oeso", "ca_ov", "ca_gastric", "ca_hep", "ca_myeloma", "ca_thyroid", "ca_bil", "ca_cerv", "ca_test")

  
###############  
# create function that returns column of binary outcomes for whole 
#empty results files
output<-rep(NA, nrow(sample_all))
output<-as.data.frame(output)

extract_outcome<-function(sample_list, sample_id="ID_1", outcomefile, outcome, outcome_id="eid"){
  stopifnot(!is.null(outcomefile[,outcome]))
  stopifnot(!is.null(outcomefile[,outcome_id]))
  stopifnot(!is.null(sample_list[,sample_id]))
  events<-which(outcomefile[,outcome]==1)
  cases = unique(outcomefile[events,outcome_id])
  output_temp = ifelse(as.numeric(sample_list[2:nrow(sample_list),sample_id])%in%cases, 1, 0)
  output_temp<-as.data.frame(output_temp)
  output_temp<-rbind("B", output_temp)
  names(output_temp)<-outcome
  return(output_temp)
}

# Apply this accross all 

allframes = lapply(wanted_outcomes, function(x) extract_outcome(sample_list=sample_all, 
                                              sample_id="ID_1", 
                                              outcomefile=outcomes, 
                                              outcome=x, 
                                              outcome_id="n_eid"))
answer = do.call(cbind,allframes)

# check same row number and bind
assertthat::are_equal(nrow(sample_all),nrow(answer))

sample_out<- cbind(sample_all, answer)

# add sex, age, age squared, bmi

library(readr)
age_sex <- read_csv("~/Stata_outcomes/age_sex.csv")
sample_both<-merge(sample_out, age_sex, by.x="ID_1", by.y="n_eid", all.x = TRUE, all.y=FALSE)
assertthat::are_equal(nrow(sample_both),nrow(sample_out))
sample_both<-sample_both[order(sample_both$id),]

# add smoking indicator
smokers <- read_csv("~/Stata_outcomes/smokers.csv")
sample_both2<-merge(sample_both, smokers, by.x="ID_1", by.y="n_eid", all.x = TRUE, all.y=FALSE)
assertthat::are_equal(nrow(sample_both),nrow(sample_both2))
sample_both2<-sample_both2[order(sample_both2$id),]


# add label for line

sample_both2[1,c("ages","sex", "bmi","agesq")] <-c("C", "D", "C","C")
sample_both2[1,c("neversmoker")]<-c("D")
sample_out<-sample_both2

#return to using UKBB ID's
sample_out[,1]<-sample_out[,2]


#check row order is same as original sample file
sample_out$id2<-1:nrow(sample_out)
assertthat::are_equal(sample_out$id, sample_out$id2)

#remove unneeded columns
sample_out<-sample_out[, !names(sample_out) %in% c("id", "exclude", "error_check","id2")]

# fix error with snptest not reading missing correctly
sample_out[1, "missing"]<-0



# save sample file
write.table(sample_out,"~/Stata_outcomes/cancer_v3.sample", row.names=FALSE, quote=FALSE)

