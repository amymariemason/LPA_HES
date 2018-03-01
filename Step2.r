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
setwd("//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress/Outputs")



# load adiposity ID sample list

sample_all = read.table("AdiposityID_sampleset", header = TRUE)


# load outcomes file and list outcomes of interest
outcomes = read.csv("C:/UKbiobank/Stata_output/All_outcomes2.csv", header=TRUE)
wanted_outcomes<-c("ukb_mi", "dvt", "pe", "ast", "af", "aa", "taa")	
#ukb_mi_nSR	ukb_stemi	ukb_stemi_nSR	ukb_nstemi	ukb_nstemi_nSR	ukb_stroke	ukb_stroke_nSR	ukb_stri	ukb_stri_nSR	ukb_ich	ukb_ich_nSR	ukb_sah	ukb_sah_nSR	cad_soft	cad_soft_nSR	cad_hard	cad_hard_nSR	cad_int	cad_int_nSR	tia	tia_nSR	pvd	pvd_nSR	dvt	dvt_nSR	pe	pe_nSR	aaa	aaa_nSR	hf	hf_nSR	dmcp	dmcp_nSR	hf_dmcp	hf_dmcp_nSR	ast	ast_nSR	af	af_nSR	ckd	ckd_nSR	hpt	hpt_nSR	taa	taa_nSR	aa	aa_nSR	hf_cmp	hf_cmp_nSR	vte	vte_nSR")
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
                                              outcome_id="n_eid"))
answer = do.call(cbind,allframes)

# check same row number and bind
assertthat::are_equal(nrow(sample_all),nrow(answer))

sample_out<- cbind(sample_all, answer)

#return to using UKBB ID's
sample_out[,1]<-sample_out[,2]


# save sample file
write.table(sample_out, "./UKBB_outcomes.sample", row.names=FALSE, quote=FALSE)

