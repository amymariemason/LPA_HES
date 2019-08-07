######
# Author: Amy Mason 
# Date: Jun 2019
# Goal: Add age, sex to a samplefile
# Inputs: sample file, agesex.csv
# 
# Outputs: new sample file including age and sex
######

# inputs
library(readr)
samplefile <- april_2019_outcomes_v2
samplelink <- read_table2("~/Programs/GWAS_inprogress/BB_input/sampleID_map.txt")
age_sex <- read_csv("~/Stata_outcomes/age_sex.csv")


# replace the UKbiobank IDs with Adiposity IDS
#sampleHRC contains two ID ('ID_1' and 'ID_2') fields which are identical. Replace the first one with the adiposity IDs 
samplefile$id <- 1:nrow(samplefile)
sample2<-merge(samplefile,samplelink , by.x="ID_1", by.y="UKB_sample_ID", all.x = TRUE)
#return to original order
sample2<-sample2[order(sample2$id),]
# replace ID_1
sample2[,"ID_1"]<-sample2$Adiposity_sample_ID
#remove unneeded columns
sample2<-sample2[, !names(sample2) %in% c("Adiposity_sample_ID", "BP_sample_ID")]


#### add age/sex/bmi
sample3<-merge(sample2, age_sex, by.x="ID_1", by.y="n_eid", all.x = TRUE, all.y=FALSE)
assertthat::are_equal(nrow(sample3),nrow(samplefile))
sample3<-sample3[order(sample3$id),]

# add label for line

sample3[1,c("ages","sex", "bmi","agesq")] <-c("C", "D", "C","C")
sample_out<-sample3

#return to using UKBB ID's
sample_out[,1]<-sample_out[,2]

#check row order is same as original sample file
sample_out$id2<-1:nrow(sample_out)
assertthat::are_equal(sample_out$id, sample_out$id2)

#remove unneeded columns
sample_out<-sample_out[, !names(sample_out) %in% c("id", "exclude", "error_check","id2")]

# save sample file
write.table(sample_out,"~/Stata_outcomes/april_2019_ou.sample", row.names=FALSE, quote=FALSE)


