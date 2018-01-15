#!/usr/bin/Rscript
###### Extract function
# Author: Amy Mason
# Date: Oct 2017
# Goal: Search large data sets of blood cell data for relevant variants, create subset if present
# Inputs: list of variants, data to search for those variants
# Outputs: log file reporting missing variants, dataframe containing the subset of larger data matching results for the variants given
######

# add: check HWE; check imputations vs genotype?

######################################
# setup desktop
######################################
# libraries
require(data.table)

####################################
# Set folders 
####################################

#working directory
setwd(home_dir)

##################################
#data input
##################################

# stats file to be scanned
statsfile<-paste0(inputs_dir, "ukb_imp_chr6_HRCvars_EURsamples_snpstats.txt")

# snps to scan for 
snplist<-read.table(paste0(inputs_dir, "lpa_snps.txt"), header=TRUE)

# NOTE TO AMY: this could be improved to search through relevant statsfiles based on snp input file


####################################
#variant input
####################################

#check format
var43<-snplist$variantID

##########################################
# Functions
#########################################

# create function to loop round variables


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

####################################################################
# create subset of snps of interest
####################################################################

# subset the data
try<-import_data(statsfile, inputname="position", snplist, varname="pos")

stats_sub<-try[[1]]

missing<- try[[3]]
# create report
message(paste0("there were ", try[[2]], "variants not found"))

message("These snps are missing")
missing[missing$missing==1, !(colnames(missing)%in%c("missing")),drop=FALSE]

message("these snps are duplicated ")
within(stats_sub[stats_sub$duplicate==1,], rm("duplicate") )
       

###########################################################
# investigate impute_info
###########################################################

write.table(stats_sub, "500k/Outputs/impute_stats", row.names=FALSE, quote=FALSE)

summary(stats_sub$impute_info)
        
message("these snps have impute values <0.9")
stats_sub[stats_sub$impute_info<0.9, c("rsid","chr","pos")]

################################################################
# random sampling of impute info


# imports data into blank file, reports missing data in list
random_data<-function(inputfile, iterate, max_start=1000000000, countalong=FALSE){
  #inputfile = file to subset
  # iterate = how many values to sample
  # max start = set longer than you think the file is
  # countalong = if you want to watch the numbers sampled
output<-create_blank_file(inputfile)  
output<-within(output,rm("duplicate"))
current_max<-max_start
i=0
while (i < iterate){
i<-i+1  
  # to watch progress
j<-sample(1:current_max,1)  
if(countalong=TRUE){print(paste0(i,":", j))}
  # test if input exists for this value; reports error if missing
  test <- tryCatch(
    # this is what I want it to do:  
    fread(file=inputfile,skip=j, nrow=1)
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
   current_max<-j-1;
  }
  # if no error, add collected lines to file, reporting duplicate matches
  if(!(is.character(test))) {
    names(test)<-names(output)
   output <- rbind(output, test, fill=TRUE)
  }
}  
  return(list(output,current_max))
}

# create random sample
sample<-random_data(inputfile, iterate=1000, max_start= 2460183)

#summarise and graph
summary(sample[[1]]$impute_info)
hist(sample[[1]]$impute_info, nclass=20)
boxplot(sample[[1]]$impute_info)

# save graph 



###########################################################################
#
###########################################################################

 