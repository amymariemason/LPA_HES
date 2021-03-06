

###################################################
# Turn table of snps into qctool ready form 

qctoolify <- function(table, chr=chr, pos=pos, filename){
  #table is the table of chromosomes and positions you want to extract
  #chr is the name of the column containing chromosone information (integrer)
  #pos is the name of the column containing position information
  #filename is where the file should be saved
  ###################
  # this function creates a space seperated list of snps ready to use in qctool
  table$output<-ifelse(as.numeric(table$chr)<10, paste0("0",table$chr,":", table$pos), paste0(table$chr,":", table$pos))
  write.table(t(table$output),filename,sep=" ",row.names=FALSE, quote =FALSE, col.names=FALSE)
}
