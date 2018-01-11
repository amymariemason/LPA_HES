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
#setwd("C://Users/am2609/Programs/GWAS_inprogress")
setwd("D:/CurrentWork")
#setwd("//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress")



# load adiposity ID sample list from Step 1

sample_all = read.table("AdiposityID_sampleset", header = TRUE)


# load outcomes file and list outcomes of interest
outcomes = read.csv("./HES outcomes/stroke_pad_all_comb.csv", header=TRUE)
wanted_outcomes<-c("isch_comb","pad_comb","sah_comb","ich_comb","haem_comb")
#########

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

allframes = lapply(wanted_outcomes,
                   function(x)extract_outcome(sample_list=sample_all, 
                                              sample_id="ID_1", 
                                              outcomefile=outcomes, 
                                              outcome=x, 
                                              outcome_id="eid"))
answer = do.call(cbind,allframes)

# check same row number and bind
assertthat::are_equal(nrow(sample_all),nrow(answer))

sample_out<- cbind(sample_all, answer)

#return to using UKBB ID's
sample_out[,1]<-sample_out[,2]


# save sample file
write.table(sample_out, "./lpa_HESoutcomes.sample", row.names=FALSE, quote=FALSE)

