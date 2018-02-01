######
# Author: Amy Mason ; Lp(a) data input (Step 4) and MR analysis (end of Step 5) based heavily on code by Stephen Burgess and James Staley
# Start Date: Oct 2017
# Goal: MR analysis from Adiposity phenotype files/LP(a) snps
# Inputs: Source files fopr Step 1, Step 2
# Outputs: Report
######
#
# The purpose of this set of programs is facilitate MR analysis with genetypic data from Biobank genetic data and
# phenotypic data from the Biobank Adiposity study 
# This occurs in several steps
#
# Firstly a file is created to link the phenotypic outcomes to the ordered sample list in UK biobank (Steps 1/2)
# Then a seperate piece of code is run to subset the Biobank data to the snps of interest; a GWAS is then run on them (Step3/snptest.sh)
# LP(a) data is loaded (this could be replaced with a repeat of the first 3 steps with a different phenotype, or 
#     your own pre-loaded data). A correlation matrix for the snps is needed as input here  (Step 4)
# Then an MR analysis is run on the results (currently IVW) (Step 5)
# Then a report is generated with a table showing each outcome, and graphs of the datapoints with lines of best fit (Report)
#
#
#Currently all the steps run individually - all this file does is ensure that their run has been logged
#In future, files will be adjusted to take their inputs from this file so that lower level editing is less nessisary
#
#
######

######################################################################
# prepare workspace: 
###########################################################################

# clear workspace
rm(list=ls())

# set home
# this should be the director containing the step1 etc
setwd("//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress/Code/")

# set log file 
#errors and warning produced by each step will be logged to this file
error_log <- "./errors.log" # alter this to change logfile name

# wipe current log
log_error <- file(error_log, open="wt")
sink(log_error, type=c("message"), append = FALSE)
message(paste("MR analysis started on: ", format(Sys.time(), "%a %b %d %X %Y")), "\n\n")
sink(type=c("message"))
close(log_error)

# check log has correctly closed
assertthat::assert_that(sink.number(type="message")==2,msg="warning: sink not closed")

##########################################################################
# set input directories
##########################################################################

#
code_dir<-getwd()

# this is the directory containing inputs/output/biobank input folders
home_dir<- "//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress"

# this is the directory containing the Biobank source files for step1
### if running on cardio, symbolic links to these files can be found in 
BB_dir<-paste0(home_dir, "/BB_inputs/")

# output file
output_dir<-paste0(home_dir, "/Outputs/")

#inputs file
inputs_dir<-paste0(home_dir, "/Inputs/")





##########################################################################
# Prep Step 1: Check imputation quality of snps
##########################################################################

#sink errors to the log file 
log_error <- file(error_log, open="a")
sink(log_error, type=c("message"), append = TRUE)

message(paste("Imputation check run on ", format(Sys.time(), "%a %b %d %X %Y")), "\n\n")


# run the file line by line to ensure all errors are caught
ll <- parse(file = "imputation_check.R")

for (i in seq_along(ll)) {
  tryCatch(eval(ll[[i]]), 
           error = function(e) message("\nOops!  ", as.character(e), "\n"),
           warning=function(warn) message("\nwarning: ", as.character(warn), "\n"))
}

message(paste("Imputation check ", format(Sys.time(), "%a %b %d %X %Y")), "\n\n")

#stop logging errors
sink(type="message")
close(log_error)

# check log has correctly closed
assertthat::assert_that(sink.number(type="message")==2,msg="warning: sink not closed")

# might be nice to add report to this


###########################################################################
# Step 1: create blank sample file 
###########################################################################

#sink errors to the log file 
log_error <- file(error_log, open="a")
sink(log_error, type=c("message"), append = TRUE)

message(paste("Step 1 run on ", format(Sys.time(), "%a %b %d %X %Y")), "\n\n")


# run the file line by line to ensure all errors are caught
ll <- parse(file = "gwas_create_sample_file_step1.R")

for (i in seq_along(ll)) {
  tryCatch(eval(ll[[i]]), 
           error = function(e) message("\nOops!  ", as.character(e), "\n"),
           warning=function(warn) message("\nwarning: ", as.character(warn), "\n"))
                }
  
message(paste("Step 1 finished run on ", format(Sys.time(), "%a %b %d %X %Y")), "\n\n")

#stop logging errors
sink(type="message")
close(log_error)

# check log has correctly closed
assertthat::assert_that(sink.number(type="message")==2,msg="warning: sink not closed")


############################################################################
# Step 2: add outcomes
############################################################################

# sink errors to log file (from here logs will append not overwrite)
log_error <- file(error_log, open="a")
sink(log_error, type=c("message"), append=TRUE)

message(paste("Step 2 run on ", format(Sys.time(), "%a %b %d %X %Y")), "\n\n")

# run the file line by line to ensure all errors are caught
ll <- parse(file = "gwas_add_outcomes_step2_2.R")

for (i in seq_along(ll)) {
  tryCatch(eval(ll[[i]]), 
           error = function(e) message("\nOops!  ", as.character(e), "\n"),
           warning=function(warn) message("\nwarning: ", as.character(warn), "\n"))
}

#stop logging
message(paste("Step 2 finished run on ", format(Sys.time(), "%a %b %d %X %Y")), "\n\n")
sink(type=c("message"))
close(log_error)
assertthat::assert_that(sink.number(type="message")==2,msg="warning: sink not closed")


###########################################################################
# Step 3
###########################################################################
# This has to run in cardio. You will need to move the .sample file output from Step 2, 
#   a list of snps of interest and the appropiete script below to cardio.
# There are two possible slurm scripts to use for this:
# 1) qc_snprun.sh: This will create a subset of the bgen file containing only the required snps and proceed as snprun.sh
# 2) snprun.sh:  This will check the association between the snps and the outcomes in your sample file
# Copy the various .out files back to your computer (there will be one for each outcome)
#
##########################################################################

##########################################################################
# Step 4: Load LP(a) data 
##########################################################################
# this step is specific to my current project and needs to be replaced with a more general option!

# open log file
log_error <- file(error_log, open="a")
sink(log_error, type=c("message"), append=TRUE)

message(paste("Step 4 run on ", format(Sys.time(), "%a %b %d %X %Y")), "\n\n")

# run the file line by line to ensure all errors are caught
ll <- parse(file = "MR_calc_loop_step_4.R")

for (i in seq_along(ll)) {
  tryCatch(eval(ll[[i]]), 
           error = function(e) message("\nOops!  ", as.character(e), "\n"),
           warning=function(warn) message("\nwarning: ", as.character(warn), "\n"))
}

# close log
message(paste("Step 4 finished run on ", format(Sys.time(), "%a %b %d %X %Y")), "\n\n")
sink(type=c("message"))
close(log_error)
assertthat::assert_that(sink.number(type="message")==2,msg="warning: sink not closed")

##########################################################################
# Step 5: MR Analysis
##########################################################################
# this step is specific to my current project and needs to be replaced with a more general option!

# open log file
log_error <- file(error_log, open="a")
sink(log_error, type=c("message"), append=TRUE)

message(paste("Step 4 run on ", format(Sys.time(), "%a %b %d %X %Y")), "\n\n")

# run the file line by line to ensure all errors are caught
ll <- parse(file = "MR_calc_loop_step_4.R")

for (i in seq_along(ll)) {
  tryCatch(eval(ll[[i]]), 
           error = function(e) message("\nOops!  ", as.character(e), "\n"),
           warning=function(warn) message("\nwarning: ", as.character(warn), "\n"))
}

# close log
message(paste("Step 4 finished run on ", format(Sys.time(), "%a %b %d %X %Y")), "\n\n")
sink(type=c("message"))
close(log_error)
assertthat::assert_that(sink.number(type="message")==2,msg="warning: sink not closed")


################################################################################
# Final: Create report
###############################################################################

require(rmarkdown) # required for md to html 
require(knitr) # required for md to html 
setwd("//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress/500k/LPA_HES/Reports")
rmarkdown::render("./Main_report.Rmd")

