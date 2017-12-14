######
# Author: Amy Mason ; based heavily on code by Stephen Burgess
# Date: Oct 2017
# Goal: GWAS - main file
# Inputs: Source files fopr Step 1, Step 2
# Outputs: Report
######

rm(list=ls())
# set working directory
setwd("//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress")
#errors and warning logged to errors.log



# Step 1: create blank sample file ----------------------------------------
log_error <- file("./errors.log", open="wt")
sink(log_error, type=c("message"), append = FALSE)

message(paste("Step 1 run on ", format(Sys.time(), "%a %b %d %X %Y")), "\n\n")
# run step 1: 

ll <- parse(file = "gwas_create_sample_file_step1.R")

for (i in seq_along(ll)) {
  tryCatch(eval(ll[[i]]), 
           error = function(e) message("\nOops!  ", as.character(e), "\n"),
           warning=function(warn) message("\nwarning: ", as.character(warn), "\n"))
                }
  
message(paste("Step 1 finished run on ", format(Sys.time(), "%a %b %d %X %Y")), "\n\n")

## reset message sink and close the file connection

sink(type="message")
close(log_error)
assertthat::assert_that(sink.number(type="message")==2,msg="warning: sink not closed")


# Step 2: add requested outcomes ----------------------------------------
log_error <- file("./errors.log", open="a")
sink(log_error, type=c("message"), append=TRUE)

message(paste("Step 2 run on ", format(Sys.time(), "%a %b %d %X %Y")), "\n\n")
# run step 2: 

ll <- parse(file = "gwas_add_outcomes_step2_2.R")

for (i in seq_along(ll)) {
  tryCatch(eval(ll[[i]]), 
           error = function(e) message("\nOops!  ", as.character(e), "\n"),
           warning=function(warn) message("\nwarning: ", as.character(warn), "\n"))
}

## reset message sink and close the file connection
message(paste("Step 2 finished run on ", format(Sys.time(), "%a %b %d %X %Y")), "\n\n")
sink(type=c("message"))
close(log_error)
assertthat::assert_that(sink.number(type="message")==2,msg="warning: sink not closed")
