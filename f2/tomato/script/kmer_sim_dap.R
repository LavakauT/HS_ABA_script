#####################
library(dplyr)
library(tibble)
library(universalmotif)
library(stringr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(gridtext)
#####################


#####################
###### NOTES ########
#####################
# PLEAESE IMPLEMENT THIS SCRIPT ON HPC SYSTEM
# IT NEEDS HUGE MEMORY TO PROCESS MOTIF FROM ALOTS OF INPUT
# INPUT1: DIRECTORY TO THE "C#_kmer.txt" FILE
# INPUT2: ARABIDOPSIS DAP-SEP MOTIF DATA BASE (absolute file path)

args <- commandArgs(trailingOnly = TRUE)

input1 <- args[1]
input2 <- args[2]


dir <- paste0(input1)
filenames <- list.files(dir,
                        pattern="*kmer.txt",
                        full.names=FALSE)
head(filenames)

dap <- read_meme(input2)

for (file in filenames) {
  
    # READ FILE
    df_file <- read.delim(paste(dir, file, sep = '/'), header = T)
    df <- colnames(df_file[,-1])
    if (length(df) == 0){
      next
    }
    file_name <- str_remove(file, "_kmer.txt")
    
    # create motifs
    list.motif <- list()
    for (i in 1:length(df)) {
      m <- create_motif(df[i], name = df[i])
      list.motif <- c(list.motif, m)
      list.pwm <- merge_motifs(list.motif, method = "PCC")
      list.pwm@name <- file_name
      }
    
    
      #####################
      ### motif compare ####
      #####################
      
      comparisons <- compare_motifs(c(list.motif, dap),
                                    method = "PCC",
                                    min.mean.ic = 0,
                                    score.strat = "a.mean")
      
      # remove pCRE to pCRE PCC
      comparisons2 <- comparisons[-c(seq_along(df)), ]
      
      
      for (i in seq_along(df)) {
        if(i == 1){
          
          df <- data.frame(comparisons2[,i])
          cname <- list.motif[[i]]@name
          colnames(df) <- cname
          rname <- row.names(df)
          df$motif <- rname
          df <- df[order(-df[,1]), ] # top similar motif
          df <- df[1,c('motif', cname)]
          row.names(df) <- 1
          
        } else {
          
          df2 <- data.frame(comparisons2[,i])
          cname2 <- list.motif[[i]]@name
          colnames(df2) <- cname2
          rname2 <- row.names(df2)
          df2$motif <- rname2
          df2 <- df2[order(-df2[,1]), ] # top similar motif
          df2 <- df2[1,c('motif', cname2)]
          row.names(df2) <- 1
          
          df <- full_join(df, df2, by = 'motif')
        }
      }
    
      df3 <- gather(df, key = 'kmers', value = 'PCC', -motif) %>% 
        rename(match = motif)
      df3 <- df3 %>% filter(PCC != '')
      df3$cluster <- file_name
    
      #####################
      ### result export ##
      #####################
      # DAP motif family
      motif.family <- df3
      
      string <- motif.family$match
      # head(string)
      
      motif.family$motif <- string
      
      motif.family <- motif.family %>%
        group_by(kmers, motif) %>%
        distinct() %>% 
        summarise_all(median)
      
      motif.family <- motif.family %>% 
        column_to_rownames(var = 'kmers')
      
      
      # DAP only
      tnt <- str_locate(motif.family$motif, '\\_') %>% 
        data.frame()
      
      tnt$end <- (tnt$end)-1
      tnt$start <- 1
      
      motif.family$motif <- str_sub(motif.family$motif,
                                    start = tnt$start,
                                    end = tnt$end)
      
      output <- data.frame(motif.family)
      output <- rownames_to_column(output, var = 'kmers')
      
      # make new name
      motif_range <- unique(output$motif)
      output2 <- data.frame()
      for (mr in motif_range) {
        output.sub <- output %>% 
          filter(motif %in% mr)
        
        # rank with PCC (decreasing = T)
        output.sub <- output.sub[order(output.sub$PCC, decreasing = T),]
        
        n_row <- nrow(output.sub)
        
        # name new_name
        output.sub <- output.sub %>% 
          mutate(new_name = paste(cluster, motif, 1:n_row, sep = '_'))
        
        # append
        output2 <- rbind(output2, output.sub)
      }
      
      # check no data loss after make new_name
      if(nrow(output2) != nrow(output)){
        warning('Error: Data loss after exporting new kmer names!')
        stop()
      }
      
      # make directory for output
      dir.create(paste(dir, file_name, sep = '/'))
      
      write.table(output2,
                  paste(dir, file_name, paste0(file_name, '_top_sim_dap.txt'), sep = '/'),
                  row.names = F, quote = F, sep = '\t')
      
      save.image(file = paste(dir, file_name, paste0(file_name, '_dap.RData'),
                              sep = '/'))
}
#####################