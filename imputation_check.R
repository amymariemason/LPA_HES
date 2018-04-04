#!/usr/bin/Rscript
###### Extract function
# Author: Amy Mason with comments and loop from Jess
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
#NOTE - statsfile is the document we want to subset and snplist represent the variants we want to select. Note that the  
#statsfiles are saved for each chromosome, therefore we only want to search through the statsfiles that are relevant  
#to that chromosome. 
#Hence we need to loop through the datasets.   

# stats file to be scanned: stats folder should be the folder containing all the statsfiles
statsfolder <-paste0(inputs_dir, "stats_folder/")

#Read in the variants we want to have information on and find out which statsfiles we need 
snps<-read.table(paste0(inputs_dir,"testsnps.txt"), header=TRUE) 
chr<-snps$chr 
chr<-chr[!duplicated(chr)] 


##########################################
# Functions
#########################################

# create function to loop round variables


create_blank_file<- function(inputfile){

 #create empty data frame for variant data; includes the column headings from inputfile and an additional column 
 #called 'duplicate' 
output<-fread(file=inputfile,nrows=1)
output$duplicate<-0
output<- output[0,]

return(output)
}

#This function creates a data frame with the variable 'varname' from varlist - this will be 'pos' from snplist  
#in our example. It also adds an additional column called 'missing' where each variant has a value of zero.  
#Hence, will provide information on variants that are in snplist but do not appear in the statsfile.   

create_missing_report<-function(varlist, varname=names(varlist)[[1]]){
  missinglist<-as.data.frame(varlist[,varname])
  names(missinglist)<-varname
  missinglist$missing<-rep(0, nrow(missinglist)) 
  return(missinglist)
}


# this changes the proposed column types by fread by overriding them
# varlist = empty file created for the variant data
# varname = vector of variables for which you wish to change class; selects all if not specifies
# coltype = vector of new classes; by default this overrides into character

change_col_character<-function(varlist, varchange=names(varlist), coltype=rep("character", length(varchange))) {
  current<- sapply(varlist, class)
  whichcols<-which(colnames(varlist)%in%varchange)
  current[whichcols]<-coltype
  return(current)
}


# imports data into blank file, reports missing data in list
import_data<-function(inputfile, inputname, varlist, varname=names(varlist)[[1]], maxdup=5, col_change=FALSE, varchange=NULL, coltype=NULL){
  #inputfile = file to subset
  #inputname = variable in input file to search for variable
  #varlist = file containing variables to search for
  # (varname) = variable to search; defaults to first column of varlist
  # col_change = TRUE -> changes variable types to character
    #varchange = vector of variable names to restrict col_change to; default is NULL
    #coltype = vector of classes to change the column to; default is NULL
  
        # adding coltype will impose those classes as long as they are not more restrictive than the default
      # see fread for more details on colClasses
  # (maxdup) = maximum number of duplications to search for; default = 5; if some allele identifying as logical instead of character
    # try increase max dup or use col_change to hard set  
  # Duplications of the position number may  
  #occur when there has been more than one mutation i.e. G to T and G to C etc. 
  
# create blank report files
missingreport<-create_missing_report(varlist,varname)
output<-create_blank_file(inputfile)  
num_er=0 # number of errors reported

# changes character types if needed
char_read<-fread(file=inputfile,nrows=1)
class_change <- sapply(char_read,class)
if(col_change==TRUE){
  class_change<-change_col_character(char_read,varchange, coltype)}
names(class_change)<-NULL
  
# loop to read in data from outside file. The fread() function is used to check whether the variant in snplist is  
#contained in the statsfile. To allow the code to run when the variant is missing, the tryCatch() function is used.  
#One of two things will happen: a) it doesn't find "j", in which case it will report an error and the variant will be  
#reported as missing; or b) if the variant is found then the fread() function will be executed and the data subseted 
#and saved as test. 

num_er =0 # error counter
for (j in varlist[,varname]){
  # to watch progress
  print(j)
  #add spaces to ensure exact match: Note this may still pull wrong line IF value matches another var
  j1<- paste0(" ",j," ")
  # clear test dummy var
#  if (exists("test")){rm(test)}
  # test if input exists for this value; reports error if missing
  test <- tryCatch(
    # this is what I want it to do, i.e. search for position "j" and then extract maxdup rows (e.g. 5). fread() will 
    #only extract data for the first instance of this, hence we are assuming the position is in some numerical order 
    #as otherwise we may miss duplicate entires of the same variant.  
    fread(file=inputfile,nrows=maxdup, skip=as.character(j1), colClasses=class_change)
    ,
    # if error occurs
    error=function(error_message) {
      message(paste0("Error AT VALUE: ", j))
      message(error_message)
      return("ERROR")
    }
  )
  # if error, report variant as not found
  #Hence, if the variant is not found in the inputfile then report as missing by replacing with 1  
  
  if(is.character(test)){
    #Look for j in the varname and replacing missing column with 1 if not there 
    
    missingreport[grep(j, missingreport[, varname]),"missing"]<-1;
    #If there the variant is not found increase the error counting variable num_er
    num_er <- num_er +1;
  }
  # if no error, add collected lines to file, reporting duplicate matches
  if(!(is.character(test))) {
    test$duplicate<-rep(0, nrow(test));
    #Takes the col names from output (i.e. the input file) and puts the names onto test. This will only work if the 
    #variant has been found and the relevant row(s) has been extracted by the fread() function. 
      names(test)<-names(output);
    test<-as.data.frame(test)
    testSubset<-test[which(test[,inputname] == j),]
    #If there is more than one row for the position then replace the duplicate entry with a 1       
    if(nrow(testSubset)!=1) testSubset$duplicate<-rep(1, nrow(testSubset));
    output <- rbind(output, testSubset)
  }
}  
#Return the output, number of variants not found, and information on missing variants. 
outlist<-list(output, num_er, missingreport)
return(outlist)
}

################################## 
#create subset of snps of interest 
################################## 

#Create blank data frame for the outcome data to be saved into:  
stats_sub<-create_blank_file(paste0(statsfolder,"ukb_imp_chr4_HRCvars_EURsamples_snpstats.txt")) 

#Only import the relevant statfiles and loop round. We are selecting variables based on the 'position' variable  
#in the statsfile.  
for(i in chr){ 
  statsfile<-paste0(statsfolder,"ukb_imp_chr", i, "_HRCvars_EURsamples_snpstats.txt") 
  #Only want to consider the variants that are on the ith chromosome.   
  snplist<-snps[snps$chr==i,] 
  try<-import_data(statsfile, inputname="position", snplist, varname="pos") 
  
  #Since the output is stored as a list we need to extract the stats information as seperate data frames.   
  stats_sub<-rbind(stats_sub,try[[1]]) 
  #Information on variants that were not present in the stats file.   
  missing<- try[[3]] 
  # create report on the number of missing variants  
  message(paste0("CHR ", i, " :there were ", try[[2]], " variants not found")) 
  message("These snps are missing") 
  missing[missing$missing==1, !(colnames(missing)%in%c("missing")),drop=FALSE] 
  
  #Information on the number of duplicated variants.   
  message("these snps are duplicated ") 
  within(stats_sub[stats_sub$duplicate==1,], rm("duplicate") ) 
} 


###########################################################
# investigate impute_info
###########################################################

write.table(stats_sub, paste0(output_dir, "/impute_stats"), row.names=FALSE, quote=FALSE)
summary(stats_sub$impute_info)

# set acceptable thresholds for imputation quality
rare<-0.9
lowfreq<-0.9
common<-0.6
        
message(paste0("these rare snps have impute values <", rare, " & Allele B frequencies <1%"))
stats_sub$threshold1<-ifelse(stats_sub$impute_info<rare & stats_sub$alleleB_frequency<0.01, 1,0)
stats_sub[stats_sub$threshold1==1, c("rsid","chromosome","position", "alleleA","alleleB")]

message(paste0("these low frequency snps have impute values <", lowfreq, " & Allele B frequencies 1%<= x <5%"))
stats_sub$threshold2<-ifelse(stats_sub$impute_info<rare & stats_sub$alleleB_frequency>=0.01 & stats_sub$alleleB_frequency<0.05 , 1,0)
stats_sub[stats_sub$threshold2==1,  c("rsid","chromosome","position", "alleleA","alleleB")]

message(paste0("these common snps have impute values <", common,"  & Allele B frequencies >5%"))
stats_sub$threshold3<-ifelse(stats_sub$impute_info<common & stats_sub$alleleB_frequency>=0.05, 1,0)
stats_sub[stats_sub$threshold3==1,  c("rsid","chromosome","position", "alleleA","alleleB")]

stats_sub$QCfail<- stats_sub$threshold3 + stats_sub$threshold2 +stats_sub$threshold1

# save table into folder
write.table(stats_sub, paste0(output_dir, "/impute_stats"), row.names=FALSE, quote=FALSE)


#####################################################################################################
#####################################################################################################
#####################################################################################################
# This was a custom bit of script to allow Jess and I to check what the average imputation levels were in a chr. 
# It does not need to be run to investigate quality of snps.


################################################################
# generates a random data file sampled from a single chr
################################################################

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
if(countalong==TRUE){print(paste0(i,":", j))}
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

random <-0
if(random>1){
# create random sample
sample<-random_data(inputfile, iterate=1000, max_start= 2460183)

#summarise and graph
summary(sample[[1]]$impute_info)
hist(sample[[1]]$impute_info, nclass=20)
boxplot(sample[[1]]$impute_info)

}

