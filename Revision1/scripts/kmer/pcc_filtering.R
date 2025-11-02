#######################
library(dplyr)
library(tibble)
library(universalmotif)
library(stringr)
library(tidyverse)
library(circlize)
library(gridtext)
#######################

args <- commandArgs(trailingOnly = TRUE)

###### SCRIPT ######
# Filter definition: one kmer match perfectly to the other one or pcc > 0.9
# consensus motif comparison
dir <- args[1]
folder_name <- args[2]
dir2 <- args[3]

filenames <- list.files(paste(dir, folder_name, sep = '/'),
                        pattern="*_FETresults.txt",
                        full.names=FALSE)

filenames2 <- list.files(paste(dir, folder_name, sep = '/'),
                          pattern="*_df_p0.01.txt",
                          full.names=FALSE)

for (i in 1:10) {
  N <- ncol(read.delim(paste(dir, folder_name, filenames2[i], sep = '/')))
  if (N != 2){
    file_name <- str_remove(filenames[i], '.pcre_FETresults.txt')
    df <- read.delim(paste(dir, folder_name, filenames[i], sep = '/')) %>%
    select(1,4)
    names(df) <- c('motif', 'pvalue')
    df$times <- paste(file_name)
    
    list.motif <- list()
    for (z in 1:nrow(df)) {
      m <- create_motif(df$motif[z], name = df$motif[z], family = df$times[z])
      list.motif <- c(list.motif, m)
    }
    
    comparisons <- compare_motifs(list.motif,
                                  method = "PCC",
                                  min.mean.ic = 0,
                                  score.strat = "a.mean")
    
    kmers <- colnames(comparisons)
    sub.com = data.frame()
    for (q in 1:nrow(df)) {
      number <- q
      kmer <- kmers[q]
      sub.com2 <- data.frame(comparisons) %>%
        select(all_of(kmer)) %>% 
        rownames_to_column(var = 'kmer')
      
      colnames(sub.com2) <- c(kmer, 'pcc')
      sub.com2 <- sub.com2 %>% 
        filter(pcc > 0.9)
      colnames(sub.com2) <- c('motif', 'pcc')
      sub.com2 <- inner_join(sub.com2, df, by = 'motif')
      sub.com2$group <- kmer
      sub.com <- rbind(sub.com, sub.com2)
    }
    
    sub.com3 <- sub.com %>% 
      group_by(group) %>% 
      summarise_all(min)
    
    sub.com3 <- sub.com3[,c('motif', 'pvalue')] %>% 
      group_by(motif) %>% 
      summarise(pvalue = min(pvalue))
    
    write.table(sub.com3,
                paste(dir2,
                      paste0(file_name,'_distinct_pcc_enriched_kmer', '.txt'),
                      sep = '/'),
                row.names = F,
                quote = F,
                sep = '\t')
    
    df2 <- read.delim(paste(dir, folder_name, filenames2[i], sep = '/'))
    df2 <- cbind(df2[,1:2], df2[,sub.com3$motif])
    
    write.table(df2,
                paste(dir2, filenames2[i], sep = '/'),
                row.names = F,
                quote = F,
                sep = '\t')
  }
}
####################