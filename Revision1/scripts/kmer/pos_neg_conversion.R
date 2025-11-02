################
library(dplyr)
library(stringr)
################

args <- commandArgs(trailingOnly = TRUE)

###### SCRIPT ######
# a loop code to change Class column from [1,0] to [pos,neg] after pCRE_Finding_FET.py
dir <- args[1]
ends <- args[3]
prefix <- 'neg'

folder_name <- args[2]
for (i in 1:10){
    number <- i
    df <- read.delim(paste(dir, folder_name, 
                            paste0(paste(prefix, folder_name, number, sep = '_'), ends), sep = '/'))
    df$Class <- str_replace_all(df$Class, '1', 'pos') # class 1 replace with postive
    df$Class <- str_replace_all(df$Class, '0', 'neg')# class 0 replace with negative
    
    # export
    write.table(df, paste(dir, folder_name,
                            paste0(paste(prefix, folder_name, number, sep = '_'), ends), sep = '/'),
    row.names = F,
    quote = F,
    sep = '\t')
}
####################