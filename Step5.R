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
corr<- corr
# format should be: correlation matrix with variantID values as row.names

#######################################
# load basic trait data
#######################################

# set working directory
setwd("//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress/500k")
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
  
  #add report to new logfile
  # create log file 
  logfile<-paste ("./Logs/", outcomespec,"_mr.log", sep="")
  sink(logfile, append=FALSE, split=TRUE)
  cat("\n This log file is showing working on the merge and mr calculations of MR_calc_loop_step_4.R with ", outcomespec, " file \n")
 
  # load data as file called output 
  output<-read.table(paste("./Inputs/", outcomespec, ".out", sep = ""), header=T)
  outcome<-output[,c("chromosome", "position", "alleleA", "alleleB", "frequentist_add_beta_1", "frequentist_add_se_1")]
  names(outcome)<-c("chr","pos","A1","A2", "outcomeBeta", "outcomeSE")
  outcome$chr<-as.character(outcome$chr)
  outcome$pos<-as.character(outcome$pos)

  # subset X-outcome data to those variants in the Y-outcome dataset
  outcomeA<-merge(outcome, X_associations, by.x=c("chr","pos","A1","A2"), by.y=c("chr", "pos", "a2", "a1"))
  outcomeB<-merge(outcome, X_associations, by.x=c("chr","pos","A1","A2"), by.y=c("chr", "pos", "a1", "a2"))
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
  sink()
  sink.number()==0
} 

#save and add to report
save(list_save, file = "./Outputs/results.RData")