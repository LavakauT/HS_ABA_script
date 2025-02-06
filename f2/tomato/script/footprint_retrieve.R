##################
library(dplyr)
library(tidyverse)
##################


###### retrieving differential kmers ######
#####################
###### NOTES ########
#####################
# PLEAESE IMPLEMENT THIS SCRIPT ON HPC SYSTEM
# IT NEEDS HUGE MEMORY TO PROCESS MOTIF FROM ALOTS OF INPUT
# INPUT1: DIRECTORY TO THE BINDetect FOLDER
# INPUT2: ARABIDOPSIS DAP-SEP MOTIF DATA BASE (absolute file path)

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]


# promoter or distal group
location <- list.dirs(path,
                 full.names = F,
                 recursive = F)
head(location)

for (loc in location) {
  # genotype groups
  pairs <- list.dirs(paste(path, loc, sep = '/'),
                     full.names = F,
                     recursive = F)
  head(pairs)
  
  # make directory for promoter or distal
  dir.create(paste(path, loc, 'result', sep = '/'))
  
  for (pair in pairs) {
    
    if (pair != 'result'){
      print(paste0('Processing: ', pair))
      data <- read.delim(paste(path, loc, pair, 'bindetect_results.txt', sep = '/')) %>% 
        filter(HS_CK_highlighted == 'True')
      
      # make directory for results
      dir.create(paste(path, loc, 'result', pair, sep = '/'))
      
      export_dir <- paste(path, loc, 'result', pair, sep = '/')
      
      write.table(data, paste(export_dir, paste0(pair, '_de_kmer.txt'), sep = '/'),
                  row.names = F, quote = F, sep = '\t')
      
      # all kmer
      motifs <- data$output_prefix

      df <- data.frame()
      for (motif in motifs) {
        print(paste0('Processing: ', motif))
        file <- read.delim(paste(path, loc, pair, motif, paste0(motif, '_overview.txt'), sep = '/')) %>% 
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

      print(paste0('Finished processing: ', pair))
      gc()
      
    }else{
      next
    }
  }
}
##########################################









