##################
library(dplyr)
library(tidyverse)
##################


#####################
###### NOTES ########
#####################
# PLEAESE IMPLEMENT THIS SCRIPT ON HPC SYSTEM
# INPUT1: DIRECTORY TO THE BINDetect FOLDER

args <- commandArgs(trailingOnly = TRUE)
input1 <- args[1]


###### retrieving differential kmers ######
path <- input1

# make directory for promoter or distal
dir.create(paste(path, 'result', sep = '/'))

location <- list.dirs(path,
                      full.names = F,
                      recursive = F)

for (loc in location) {
  
  if (loc != 'result'){
    print(paste0('Processing: ', loc))
    data <- read.delim(paste(path, loc, 'bindetect_results.txt', sep = '/'))
    # filter(HS_CK_highlighted == 'True') is not suitable here.
    
    # make directory for results
    dir.create(paste(path, 'result', loc, sep = '/'))
    
    export_dir <- paste(path, 'result', loc, sep = '/')
    
    write.table(data, paste(export_dir, paste0(loc, '_de_kmer.txt'), sep = '/'),
                row.names = F, quote = F, sep = '\t')
    
    # all kmer
    motifs <- data$output_prefix
    
    df <- data.frame()
    for (motif in motifs) {
      print(paste0('Processing: ', motif))
      file <- read.delim(paste(path, loc, motif, paste0(motif, '_overview.txt'), sep = '/')) %>% 
        mutate(kmer_motif_name = motif)
      
      # append data
      df <- rbind(df, file)
    }
    
    
    # export bindetect files and running records
    df %>%
      distinct() %>%
      write.table(., paste(export_dir, 'bindetect_results.txt', sep = '/'),
                  row.names = F, quote = F, sep = '\t')
    save.image(paste(export_dir, 'overview.RData', sep = '/'))
    
    print(paste0('Finished processing: ', loc))
    gc()
    
  }else{
    next
  }
}
#######################################