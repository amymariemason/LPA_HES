######
# Author: Amy Mason ; based heavily on code by Stephen Burgess
# Date: Oct 2017
# Goal: Perform GWAS on biobank stroke outcomes 
# Inputs: "AdiposityID_sampleset" sample file, empty of outcomes  - in same order as 500k HRConly_EurQCp .bgen file
# 
# Outputs: sample file with outcomes from stroke_pad_all_comb
######


# libraries and setup
#setwd("C://Users/am2609/Programs/GWAS_inprogress")
#setwd("//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress/Outputs")



# load adiposity ID sample list

sample_all = read.table(paste0(output_dir,"AdiposityID_sampleset"), header = TRUE)


# load outcomes file and list outcomes of interest
outcomes = read.csv("~/Stata_outcomes/Hannah.csv", header=TRUE)
wanted_outcomes<-c("ukb_mi", "ukb_stroke", "cad_soft", "tia", "pvd_simple", "vte",
                   "hf", "hf_cmp", "af", "aa", "dmcp", "ast", "han_ckd", "han_mc",
                   "han_edc", "han_rhd", "han_sleepap", "han_rtp", "han_fgr", "tin", "park",
                    "ukb_acd", "ukb_alz", "ukb_vd", "ukb_ftd")

  
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

# add label for line

sample_both[1,c("ages","sex", "bmi","agesq")] <-c("C", "D", "C","C")
sample_out<-sample_both

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
write.table(sample_out,"~/Stata_outcomes/Hannah.sample", row.names=FALSE, quote=FALSE)

