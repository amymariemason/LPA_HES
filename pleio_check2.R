######
# Author: Jess Rees
# Date: Feb 2018
# Goal: Detecting pleiotropic variants by using phenoscanner 
######

#Using a cut off of GWS (5*10^-8) prodcues 224 in the all, 155 no proxy and 122 proxy 

rm(list=ls())
library(data.table)
library(dplyr)
setwd("/Users/jessica/Dropbox/BMI_asthma_MR/variants/phenoscanner")

#Note - the variants have already been pruned with a p-value<0.05 
vars_locke1<-read.csv(file="locke_oldtxt_33517_PhenoScanner_GWAS.csv", na.strings=c("","NA"), header=TRUE)
vars_locke2<-read.csv(file="locke_noveltxt_34040_PhenoScanner_GWAS.csv", na.strings=c("","NA"), header=TRUE)
vars_original<-rbind(vars_locke1,vars_locke2)

#Group of variants 
group_vars<-read.csv(file="/Users/jessica/Dropbox/BMI_asthma_MR/variants/adults.csv", na.strings=c("","NA"), header=TRUE)

#Remove the variants (either the 77 variants or those in high LD) that are associated with BMI from the GIANT 
#consortium from any of their publications. This includes any stratified analyses 
for (i in grep("^BMI",levels(vars_original$Trait))){
  vars_original<-vars_original[!(vars_original$Study=="GIANT" & vars_original$Trait==levels(vars_original$Trait)[i]),]
  print(dim(vars_original)[1])
}

#Remove all variants that are associated with BMI from any of the studies since we already know they are associated 
#with BMI and therefore it is not useful information.
for (i in c(54,105,
            grep('Body mass index',levels(vars_original$Trait),ignore.case=TRUE))){
  vars_original<-vars_original[vars_original$Trait!=levels(vars_original$Trait)[i],]
  print(dim(vars_original)[1])
}

#Re-label some of the traits as it is case sensitive 
old<-c("Alzheimer Disease", 	"Alzheimers disease biomarkers", 	"Blood Pressure", 	"C reactive protein CRP", 	"Coronary artery disease CAD", 	"HDL", 	"LDL", 	"LDL cholesterol in serum", 	"Type 2 diabetes combined control dataset", 	"Type 2 diabetes end stage renal disease", 	"Type II diabetes", 	"Diabetes Mellitus Type 2","Longevity exceptional longevity") 
new<-c('Alzheimers disease', 	'Alzheimers disease', 	'Blood pressure', 	'C reactive protein', 	'Coronary artery disease', 	'HDL cholesterol', 	'LDL cholesterol', 	'LDL cholesterol', 	'Type 2 diabetes', 	'Type 2 diabetes', 	'Type 2 diabetes', 	'Type 2 diabetes','Longevity')
for (i in 1:length(old)){
  levels(vars_original$Trait)[which(levels(vars_original$Trait)==old[i])]<-new[i]
}

###################################################
#Function to loop through variants and studies 
###################################################

#Explanation of function: 
#1) Proxy.Indicator==0 - when we are not considering proxy variants, then we want to loop through each SNP to 
#identify which studies have signifciant associations with this SNP. Since some studies have multiple traits (e.g.
#the GIANT study has lots of measures of adiposity), there may be multiple occasions when the variant is significantly
#associated with the trait in the same study. If this is the case we just take the entry from the first trait from the 
#most recent publication, and record the number of traits there are for each SNP in no_vars.

#2) Proxy.Indicator==1 - now considering cases where the variant is a proxy for the original SNP. Note that some of 
#of the original SNPs will not have proxies since we have only extracted information on variants that are in LD with an 
#r2>0.6. We can be more strict on this cut-off within the function. Now we have the possibility that for each study there
#will be multiple proxy SNPs and multiple traits. We want to take the entry from the first trait from the most recent 
#publication that is in the highest LD with the original SNP, and record the number of proxy SNPs there are for the 
#trait(s) in no_vars. 

#3) Alternatively we consider both sets together, hence c(0,1). In this case we search through each study for the 77
#variants we then extract the variant from the most recent publication with the highest r2 value. This will either 
#be the actual variant or a proxy to the variant. Probably better to use this appraoch as otherwise there will be a 
#lot of overlap between the results from 1) and 2).  

#NOTE - the value of no_vars may not be accurate. Some of the genetic associations were extracted twice. Once from the 
#original paper and then again from the NCBI database. Also we may have results from both the European and Mixed samples
#As a result of this, we want to create an additional variable which takes all of the unique traits from each variant 
#and study and recrods them altogether. For example, this will be handy if we are considering which other adiposity 
#traits the variants are associated with apart from BMI.  

#When we are considering a less stringent threshold we are more interested in the p-value rather than the r2 value,
#therefore it would be better to sort the variants in order of the ascending p-value to select the most strongly 
#associated variant. See the argument 'descending' refering to whether we are considering the r2 value or p-value.  

pleio<-function(proxy, LD, p.value, descending){
  comp<-NA
  if (length(proxy)==1){
    vars<-vars_original[which(vars_original$Proxy.Indicator==proxy & vars_original$r2>LD & vars_original$P<p.value),]
  }
  else {
    vars<-vars_original[which(vars_original$r2>LD & vars_original$P<p.value),]
  }
  study_u<-unique(as.numeric(vars$Study))
  var_u<-unique(as.numeric(vars$rsID))
  for (j in 1:length(var_u)){
    for (i in 1:length(study_u)){
      x<-vars[which(vars$Study==levels(vars$Study)[study_u[i]] & vars$rsID==levels(vars$rsID)[var_u[j]]),]
      if (descending==0){
        x<-setorder(x,-Year.of.Publication,-r2)
        if (dim(x)[1]!=0){ #Condition needed otherwise error will occur when there are no variants
          x$no_vars<-nrow(x) #May be misleading - see above 
          x$Traits<-paste(unique(x$Trait), collapse=", ") #All of the traits within the study for a specific variant that are significant
          comp<-rbind(comp,x[1,])
        }
      }
      else {
        x<-setorder(x,-Year.of.Publication, P)
        if (dim(x)[1]!=0){ #Condition needed otherwise error will occur when there are no variants
          x$no_vars<-nrow(x) #May be misleading - see above 
          x$Traits<-paste(unique(x$Trait), collapse=", ") #All of the traits within the study for a specific variant that are significant
          comp<-rbind(comp,x[1,])
        }
        
      }
    }
  }
  comp<-comp[-1,]
  comp<-comp[comp$Study!="EGGC",]
  return(comp)
}

######################################
#Function for summarising associations
######################################

#This function creates a table with all of the unique values of var_group and then provides a list of the factors
#in var_list that are associated with each unique value of var_group.  
summary_traits<-function(var_list, var_group){
  all_traits<-data.frame(
    var1 = rep(names(lapply(split(var_list,var_group),unique)), lapply(lapply(split(var_list,var_group),unique), length)),
    var2 = unlist(lapply(split(var_list,var_group),unique)),row.names = NULL)
  d<-paste(unique(factor(all_traits$var1)))
  gws_traits<-data.frame(matrix(ncol = 2, nrow = length(d)))
  for (i in 1:length(d)){
    gws_traits[i,1]<-d[i]
    gws_traits[i,2]<-paste(all_traits[which(all_traits$var1==d[i]),2],collapse=", ")
  }
  return(gws_traits)
}



######################################
#GWS threshold
######################################

#All of the variants together 
all<-pleio(proxy=c(0,1), LD=0.61, p.value=5*10^(-8),descending=0)

#Remove other adiposity measures 
other_adi<-c("Body fat percentage",	"Obesity",	"Waist circumference",	"Height adjusted BMI",	"Pericardial fat",	"Waist hip ratio",	"Extreme obesity with early age of onset",	"Hip circumference",	"Obesity in children and adolescents with early age of onset ",	"Fat body mass",	"Height",	"Visceral adipose tissue female",	"Percent body fat",	"Childhood weight",	"Weight",	"Obesity early onset extreme",	"Body Weight",	"Percent fat bioelectric impedance measure",	"Established obesity childhood retrospective perfectionism  ",	"Subcutaneous adipose tissue male",	"Waist circumference in Type 2 diabetes")
pleio_summary<-function(data.in){
  print(paste0(dim(data.in[data.in$Study=="GIANT",])[1]," variants are associated with at least one other adiposity measure in the GIANT consortium"))
  x<-data.in[data.in$Study!="GIANT",]
  #Remove variants that are associated with other adiposity measurements across all of the other studies 
  for (i in 1:length(other_adi)){
    level<-which(levels(vars_original$Trait)==other_adi[i])
    x<-x[x$Trait!=levels(vars_original$Trait)[level],]
  }
  print(paste0(length(unique(as.numeric(x$rsID))), " are associated with another trait not related to adiposity"))
  return(x)
}
all_reduced<-pleio_summary(data.in=all)
write.csv(all_reduced,file="summary_gws.csv",row.names=FALSE)

#Create a table with the non-adiposity traits with the variants they were associated with: 
gws_traits<-summary_traits(var_list=all_reduced$rsID, var_group=all_reduced$Trait)
write.csv(gws_traits,file="non_adiposity_traits_gws.csv",row.names=FALSE)

#Consider the other way - i.e. what traits are the variants associated with. Note here we are using the Traits variable 
#we created in the pleio function since this will create a complete list of all the traits associated with the variant.
gws_rsid<-summary_traits(var_list=all_reduced$Traits, var_group=all_reduced$rsID)

#Extract the data where the variants are associated with another adiposity measure either in GIANT or another study. 
#Then create a table with the variants are associated with the adiposity measures and combine with above 
all_adiposity<-all[which(!as.numeric(row.names(all))%in%as.numeric(row.names(all_reduced))),]
gws_rsid_adiposity<-summary_traits(var_group=all_adiposity$rsID, var_list=all_adiposity$Trait)
gws_rsid_combined<-merge(gws_rsid,gws_rsid_adiposity,by.x="X1",by.y="X1",all.x=TRUE,all.y=TRUE)
colnames(gws_rsid_combined)<-c("rsID","non_adiposity","adiposity")
gws_rsid_combined<-merge(group_vars[,c(1,2,6)],gws_rsid_combined,by.x="rsid",by.y="rsID",all.x=TRUE)
write.csv(gws_rsid_combined,file="non_bmi_variants_gws.csv",row.names=FALSE)
