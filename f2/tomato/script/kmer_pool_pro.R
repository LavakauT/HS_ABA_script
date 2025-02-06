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
# YOU CAN IMPLEMENT THIS SCRIPT ON HPC SYSTEM
# IT NEEDS HUGE MEMORY TO PROCESS MOTIF FROM ALOTS OF INPUT
# INPUT1: DIRECTORY TO THE UPPER DIRECTORY OF EACH "C#_top_sim_dap.txt" FILE
# INPUT2: OUTPUT MEME TYPE MOTIF NAME

args <- commandArgs(trailingOnly = TRUE)

input1 <- args[1]
input2 <- args[2]


dir <- paste0(input1)
dirs <- list.dirs(dir,
                  full.names=F,
                  recursive = F)
# dirs should be each cluster name (C1 to C15)
head(dirs)

# append all cluster motifs
df <- data.frame()
for (d in dirs) {
  
  # check whether file exist
  file <- paste0(d, '_top_sim_dap.txt')
  if(!file.exists(paste(dir, d, file, sep = '/'))){
    warning(paste0('File: ', file, ' not found in directory!'))
    next
  }
  
  # cluster name
  # head(d)
  
  # read "C#_top_sim_dap.txt" in each cluster
  df.sub <- read.delim(paste(dir, d, file, sep = '/')) %>% 
    dplyr::select(kmers, new_name)
  
  # append
  df <- rbind(df, df.sub)
}


# length of all kmers
len_kmer <- length(unique(df$kmers))
print(paste0('Total kmers found in your dataset: ', len_kmer))

# perfect match kmers among clusters,
# we give multiple name (i.g, CTTTTT name is C1_MYB_1, C2_MYB_1)

df_wide <- df %>% 
  distinct() %>% 
  mutate(presence = 1) %>% 
  spread(., key = 'new_name', value = 'presence')

df_adj <- data.frame()
for (i in 1:len_kmer) {
  
  # select kmer
  k.sub <- df_wide[i,] %>% 
    gather(., key = 'new_name', value = 'presence', -kmers) %>% 
    filter(!is.na(presence))
  
  # how many names match the kmer
  nrow_ksub <- nrow(k.sub)
  
  if(nrow_ksub == 1){
    k.sub.new <- k.sub
    
    # append kmers with multiple name
    df_adj <- rbind(df_adj, k.sub.new)
  }else{
    for (z in 1:nrow_ksub) {
      if(z == 1){
        new_name2 <- k.sub[z,]$new_name
      }else{
        tem_name <- k.sub[z,]$new_name
        new_name2 <- paste(new_name2, tem_name)
      }
      
      k.sub.new <- data.frame(kmers = k.sub[z,]$kmers,
                              new_name = new_name2,
                              presence = nrow_ksub)
    }
    
    # append kmers with multiple name
    df_adj <- rbind(df_adj, k.sub.new)
  }
}

write.table(df_adj, paste(dir, paste0(input2, '.txt'), sep = '/'),
            row.names = F, quote = F, sep = '\t')


# make MEME motif file for ATAC-seq BINDetect
# create motifs
list.motif <- list()
list.pwm <- list()
for (i in 1:len_kmer) {
  m <- create_motif(df_adj[i,]$kmers, name = df_adj[i,]$new_name)
  list.motif <- c(list.motif, m)
  m.pwm <- convert_type(m, 'PWM')
  list.pwm <- c(list.pwm, m.pwm)
}

# export PFM
# spaces in motif name will be substituted with dashed
write_meme(list.motif, paste(dir, paste0(input2, '_pfm.meme'), sep = '/'))
#export PWM
# spaces in motif name will be substituted with dashed
write_meme(list.pwm, paste(dir, paste0(input2, '_pwm.meme'), sep = '/'))