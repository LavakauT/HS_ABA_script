####loading package####
library(dplyr) 
library(tibble)
library(universalmotif)
library(stringr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(gridtext)
require(gridExtra)
library(yaml)
####K-mer similarity####
setwd('/RAID1/working/R425/lavakau/kmer/mp')
dir <- getwd()
sub.dir <- list.dirs(dir,
                     full.names = F,
                     recursive = F)
dir2 <- '/RAID1/working/R425/lavakau/kmer/mp_similarity'
sub.dir2<- list.dirs(dir2,
                     full.names = F,
                     recursive = F)
# DAP-seq database from MEME-Suit (Arabidopsis)
dap <- read_meme('/RAID1/working/R425/lavakau/kmer/ArabidopsisDAPv1.meme')

for (y in seq_along(sub.dir)) {
  # meme file
  file_name <- sub.dir[y]
  df <- read_meme(paste(dir, sub.dir[y], paste0(sub.dir[y], '_consensus_motifs.meme'), sep = '/'))
  # function compare_motifs() included multiple parameters.
  # please check the manuals for comprehensive description.  
  comparisons <- compare_motifs(c(df, dap),
                                method = "PCC", # Pearson Correlation Coefficient
                                min.mean.ic = 0,
                                score.strat = "a.mean")
  # remove K-mer to K-mer PCC
  comparisons2 <- comparisons[-c(1:length(df)), ]
  
  for (z in 1:length(df)) {
    if(z == 1){
      name <- df[[z]]@name
    }else{
      name2 <- df[[z]]@name
      name <- c(name, name2)
    }
  }
  
  for (i in 1:length(df)) {
    if(i == 1){
      data <- data.frame(comparisons2[,i])
      cname <- df[[i]]@name
      colnames(data) <- cname
      rname <- row.names(data)
      data$motif <- rname
      names(data) <- c('mer', 'motif')
      data <- data[order(-data$mer),][1,] # top similar motifs
      names(data) <- c(cname, 'motif')
      data <- data[,c('motif', cname)]
      row.names(data) <- 1
      
    } else {
      
      data2 <- data.frame(comparisons2[,i])
      cname2 <- df[[i]]@name
      colnames(data2) <- cname2
      rname2 <- row.names(data2)
      data2$motif <- rname2
      names(data2) <- c('mer', 'motif')
      data2 <- data2[order(-data2$mer),][1,] # top similar motifs
      names(data2) <- c(cname2, 'motif')
      data2 <- data2[,c('motif', cname2)]
      row.names(data2) <- 1
      
      data <- full_join(data, data2, by = 'motif')
    }
    
  }
  
  df3 <- gather(data, key = 'kmers', value = 'PCC', -motif)
  df3 <- df3 %>% filter(PCC != '')
  df3$cluster <- file_name
  
  # DAP motif family strings manipulation
  motif.family <- df3
  
  
  string <- motif.family$motif
  # remove unwanted dash-line and any digit behind it as family name
  # string <- str_replace_all(string, '\\_m1', '')
  # head(string)
  
  motif.family$motif <- string
  
  motif.family <- motif.family %>%
    group_by(kmers, motif) %>%
    distinct() %>% 
    summarise_all(median)
  
  motif.family <- motif.family %>% 
    column_to_rownames(var = 'kmers')
  
  
  # final DAP-seq motif name processing
  tnt <- data.frame(str_locate(motif.family$motif, '\\_tnt\\.'))
  
  tnt$start2 <- 1
  tnt$end2 <- (tnt$start)-1
  
  
  motif.family$family <- str_sub(motif.family$motif,
                                 start = tnt$start2,
                                 end = tnt$end2)
  
  # find gene Arabidopsis gene name
  tnt2 <- data.frame(str_locate(motif.family$motif, '\\_col'))
  tnt2$start2 <- tnt$end + 1
  tnt2$end2 <- tnt2$start - 1
  
  motif.family$gene <- str_sub(motif.family$motif,
                               start = tnt2$start2,
                               end = tnt2$end2)
  
  output <- data.frame(motif.family)
  output <- output[name,] %>% 
    rownames_to_column(., var = 'name')
  fam <- data.frame(t(table(output$family))) %>% 
    pull(Var2) %>% 
    as.character()
  
  for (x in seq_along(fam)) {
    if(x == 1){
      output %>% 
        filter(family %in% fam[x]) %>% 
        mutate(new_name = paste(cluster, family, 1:nrow(.), sep = '_')) -> out
    }else{
      output %>% 
        filter(family %in% fam[x]) %>% 
        mutate(new_name = paste(cluster, family, 1:nrow(.), sep = '_')) -> out2
      rbind(out, out2) -> out
    }
  }
  out %>% 
    column_to_rownames(., var = 'name') -> out
  out[name,] -> out
  
  # rename motif name in format of 'cluster_Tf-family_num'
  for (x in seq_along(df)) {
    out[x,]$new_name -> df[[x]]@name
  }
  
  # rank consensus motif and export top5
  # apply the minimum rank from each motif's altname
  alt <- read.delim(paste(dir2, sub.dir2[y], paste0(sub.dir2[y], '_top_sim_dap.txt'), sep = '/')) %>% 
    pull(K.mers)
  yml <- read_yaml(paste(dir, sub.dir[y], paste0(sub.dir[y], '_clusters.yml'), sep = '/'))
  
  for (x in seq_along(yml)) {
    if(x ==1){
      bg <- data.frame(name = names(yml[x]),
                       kmer = yml[[x]])
    }else{
      bg2 <- data.frame(name = names(yml[x]),
                        kmer = yml[[x]])
      bg <- rbind(bg, bg2)
    }
  }
  
  # str_split(bg$motif, '\\s')
  # str_extract(bg$motif, '([^ ]+)')
  bg$kmer<- str_extract(bg$kmer, '([^ ]+)')
  bg[which(bg$kmer %in% alt),] -> bg
  rownames(bg) <- 1:nrow(bg) # remove yml format
  bg %>% 
    column_to_rownames(., var = 'kmer') -> bg
  data.frame(name = bg[alt,]) -> bg
  inner_join(bg, out %>% 
               rownames_to_column(., var = 'name'), by = 'name') -> output2
  
  # export
  dir3 <- paste(dir, file_name, sep = '/')
  write.table(output2,
              paste(dir3, paste0(file_name, '_top5_consensus_dap.txt'), sep = '/'),
              row.names = F,
              quote = F,
              sep = '\t') # export top5
  write.table(out %>% 
                rownames_to_column(., var = 'name'),
              paste(dir3, paste0(file_name, '_all_consensus_dap.txt'), sep = '/'),
              row.names = F,
              quote = F,
              sep = '\t') # export all
  
  save.image(file = paste(dir3, paste0(file_name, '_consensus_dap.RData'), sep = '/'))
  
}