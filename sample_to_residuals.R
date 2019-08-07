######
# Author: Amy Mason 
# Date: Jun 2019
# Goal: Change a sample file of outcomes to one of residuals of regression on coefficients and those missing from PC replaced with NA in outcome
# (this should speed up gwas)
# Inputs: sample file, agesex.csv
# 
# Outputs: new sample file if residuals
######

#################### inputs
library(readr)
april_2019_outcomes_v3 <- read_table2("~/Stata_outcomes/april_2019_outcomes_v3.sample")
QCed_Eur_unrelated <- read_csv("~/Programs/GWAS_inprogress/BB_input/QCed_Eur_unrelated.txt",col_names = FALSE)

################ outcomes of interest set to NA where not in non-related european
# drop variables not of intereste
outcomes<-april_2019_outcomes_v3[,c("ID_1", "ID_2", "missing", "PC1", "PC2","PC3","PC4", "PC5", "PC6", "PC7", "PC8", "PC9","PC10", "ages", "sex", "vte", "dvt")]

# create list of wanted IDs as character
QCed_Eur_unrelated$char <-as.character(QCed_Eur_unrelated$X1)

# create marker of people not in those lists
outcomes$exclude<- !(outcomes$ID_1 %in% QCed_Eur_unrelated$char)

# replace outcomes with NA for those people
outcomes$vte<-ifelse(outcomes$exclude==TRUE, NA, outcomes$vte)
outcomes$dvt<-ifelse(outcomes$exclude==TRUE, NA, outcomes$dvt)

################ regression of covariates on outcome

#outcomes_regress<-outcomes[2:nrow(outcomes),]
#outcomes_regress


#vtelogit <- glm(vte ~ ages + sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10, data = outcomes[outcomes$exclude==FALSE,], family = "binomial")
#dvtlogit <- glm(dvt ~ ages + sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10, data = outcomes[], family = "binomial")
# DO REGREssion in stata COZ R is being a butt

residuals_dvt_vte <- read_csv("~/Stata_outcomes/residuals_dvt_vte.csv", col_types = cols(id_2 = col_skip()))

outcomes$count<-1:nrow(outcomes)
outcomes2<- merge(outcomes, residuals_dvt_vte, by.x="ID_1", by.y="id_1", all=TRUE)
outcomes2<-outcomes2[order(outcomes2$count),]

# label needed columns

outcomes2[1,"vte"]<-"B"
outcomes2[1,"dvt"]<-"B"
outcomes2[1,"dvt_res"]<-"D"
outcomes2[1,"vte_res"]<-"D"

#remove unneeded columns
outcomes3<-outcomes2[, !names(outcomes2) %in% c("missing", "exclude", "count")]

# save sample file
write.table(outcomes3,"~/Stata_outcomes/vte_dvt_NA.sample", row.names=FALSE, quote=FALSE)


