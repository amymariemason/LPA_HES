######
# Author: Amy Mason ; based heavily on code by James Staley  and Stephen Burgess
# Date: Oct 2017
# Goal: Perform Mendelian randomization  
# Inputs: list of variants (HEStraitlist), 
#     a name/units (if applicable) file for those variants, (Traitnames)
#     association data from Step 3 to search for those variants (.out files)
#     X-associations/ data_rho file from Step 4
# Outputs: log file reporting missing variants for each .out file (./Logs/*_mr.log )
      # dataframe containing the subset of larger data matching the variants given (./Outputs/results.RData) 
######

# future goals: offer rescaling by units

#####################################################
##### Set-up #####
#####################################################

library(ggplot2)
library(MendelianRandomization)
library(plotly)
library(htmlwidgets)
library(assertthat)
library(stringr)

####################################################
# load data from Step 4 or other sources
####################################################

X_associations<-X_associations
# format should be:  c("variantID", "snp", "chr.pos", "chr", "pos" , "a1", "a2", "beta", "se")]
# program will reorder 
corr<- corr
# format should be: correlation matrix with variantID values as row.names

#######################################
# load basic trait data
#######################################

# set working directory
setwd(inputs_dir)
# this directory should contain your input files as detailled above
#setwd("D:/CurrentWork/500k")


# import list of traits to include
TL<-scan("./HESTraitlist.txt", what="", sep="\n")
TLframe<-as.data.frame(TL)
names(TLframe)<-"trait"
TLframe$traitmin<- substring(TLframe$trait,1,10)

# load names for traits 
traitnames<-read.csv("./Traitnames.csv", header=FALSE)
traitnames<-traitnames[,1:2]
names(traitnames)<-c("trait", "name")
TLall<-merge(TLframe, traitnames, by="trait", all.x=TRUE)

###########################################################################
# import data from .outfiles
###########################################################################

# create blank save file
list_save<-vector('list', length(TL))

# loop over each trait outcome
for (i in 1:length(TL)){
  # print outcome name
  outcomespec<-as.character(TLall[i,"trait"])
  outcomename<-as.character(TLall[i,"name"])
  cat("\n", as.character(outcomename),"\n")
  
  # load data as file called output 
  output<-read.table(paste(inputs_dir, outcomespec, ".out", sep = ""), header=T)
  outcome<-output[,c("chromosome", "position", "alleleA", "alleleB", "frequentist_add_beta_1", "frequentist_add_se_1")]
  # identity alignment of data -> beta is measured as effect of additional copies of A1. Coming from snptest :
  # "“SNPTEST codes allele_A as 0 and allele_B as 1 and this defines the meaning of the beta's and there se's. 
  # For example, when using the additive model the beta estimates the increase in log-odds that can be attributed 
  # to each copy of allele_B. "
  names(outcome)<-c("chr","pos","A2","A1", "outcomeBeta", "outcomeSE")
  outcome$chr<-as.character(outcome$chr)
  outcome$pos<-as.character(outcome$pos)
  
  # create reversed data set of outcome for directional realignment
  #because if I realign the X_associations, I'll need to realign the correlation matrix as well
  outcome_rev<-outcome
  outcome_rev$outcomeBeta<-(-1)*as.numeric(outcome_rev$outcomeBeta)
  
  # subset X-outcome data to those variants in the Y-outcome dataset
  outcomeA<-merge(outcome, X_associations, by.x=c("chr","pos","A1","A2"), by.y=c("chr", "pos", "a1", "a2"))
  outcomeB<-merge(outcome_rev, X_associations, by.x=c("chr","pos","A2","A1"), by.y=c("chr", "pos", "a1", "a2"))
  outcome2<-rbind(outcomeA, outcomeB)
  outcome3<-na.omit(outcome2)
  
  
  # check not lost any variants from X-outcome data; report to log if so
  #
  cat("\nThis shows an error if variants missing from Y-outcome file compared to X outcome file \n")
  tryCatch({stopifnot(nrow(outcome3) == nrow(X_associations))
  }, error = function(err.msg){
    # Add error message to the error log file
    cat("ERROR", "\n", nrow(X_associations[!(X_associations$variantID%in%outcome3$variantID),]), "variants missing from outcome data\n")
    write.table(X_associations[!(X_associations$variantID%in%outcome3$variantID),c("variantID")], row.names=FALSE, quote=FALSE, col.names=FALSE)
  }
  )
  
  cat("\nThis shows an error if variants missing from X-outcome file compared to Y outcome file \n")
  tryCatch({stopifnot(nrow(outcome3) == nrow(outcome))
  }, error = function(err.msg){
    # Add error message to the error log file
    outcome_missing<-merge(outcome, outcome3, by=c("chr","pos","A1","A2"), all.x=TRUE)
    outcome_missing<-outcome_missing[is.na(outcome_missing$variantID),]
    cat("ERROR", "\n", nrow(outcome_missing), "variants missing from outcome data\n")
    write.table(outcome_missing[,c("chr","pos","A1","A2")], row.names=FALSE, quote=FALSE, col.names=FALSE)
  }
  )
  
  # subset correlation matrix further, to match only variants in both association sets
  wantedcells<- match(outcome3$variantID, row.names(corr))
  rho <- corr[wantedcells,wantedcells]
  
  # Set-up for mr  
  bx <- as.numeric(outcome3$beta)
  by <- as.numeric(outcome3$outcomeBeta)
  byse <- as.numeric(outcome3$outcomeSE)
  bxse <- as.numeric(outcome3$se)
  rho <- as.matrix(rho)
  
  # Analysis
  MRdata_input <- mr_input(bx, bxse, by, byse, corr=rho, outcome=outcomename, exposure="Lp(a)", snps=outcome3$snp)
  
  
  tryCatch({mr = mr_ivw(mr_input(bx, bxse, by, byse, corr=rho))
  }, 
  error = function(error_message){
    cat("ERROR", "\n", outcomespec, " did not run correctly in mr_ivw\n")
    print(error_message)
  },
  warning=function(warn){
    cat("WARNING\n", "outcomespec gave warning in mr_ivw: \n")
    print(warn)
  }, 
  finally ={
    mr = mr_ivw(mr_input(bx, bxse, by, byse, corr=rho))
  }
  )
  cat("MR analysis results \n")
  print(mr)
  summary_temp<-list(outcomespec, outcomename,MRdata_input,mr)
  list_save[[i]]<-summary_temp
  # end log file

} 

#save and add to report
save(list_save, file = paste0(output_dir, "test_results.RData"))

# make report

rmarkdown::render('C:/Users/am2609/Code/Reports/Main_report.Rmd',
                  output_file = paste0(output_dir,'test_results.html'), 
                  params=list(dataset=paste0(output_dir, "test_results.RData")))