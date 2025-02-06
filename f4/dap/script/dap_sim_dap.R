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

dap <- read_meme('f2/ArabidopsisDAPv1.meme')
motif <- read_meme('f4/dap/motif_1000/meme.txt_0.01_DAP-seq')

#####################
### motif compare ####
#####################

comparisons <- compare_motifs(c(motif, dap),
                              method = "PCC",
                              min.mean.ic = 0,
                              score.strat = "a.mean")

# remove pCRE to pCRE PCC
comparisons2 <- comparisons[-c(1:length(motif)), ]


for (i in seq_along(1:length(motif))) {
  if(i == 1){
    
    df <- data.frame(comparisons2[,i])
    cname <- motif[[i]]@name
    colnames(df) <- cname
    rname <- row.names(df)
    df$motif <- rname
    df <- df[order(-df[,1]), ] # top similar motif
    df <- df[1,c('motif', cname)]
    row.names(df) <- 1
    
  } else {
    
    df2 <- data.frame(comparisons2[,i])
    cname2 <- motif[[i]]@name
    colnames(df2) <- cname2
    rname2 <- row.names(df2)
    df2$motif <- rname2
    df2 <- df2[order(-df2[,1]), ] # top similar motif
    df2 <- df2[1,c('motif', cname2)]
    row.names(df2) <- 1
    
    df <- full_join(df, df2, by = 'motif')
  }
}

df3 <- gather(df, key = 'dap-seq_motifs', value = 'PCC', -motif) %>% 
  dplyr::rename(match = motif)
df3 <- df3 %>% filter(PCC != '')



#####################
### result export ##
#####################
# DAP motif family
motif.family <- df3

string <- motif.family$match
# head(string)

motif.family$motif <- string

motif.family <- motif.family %>%
  group_by(`dap-seq_motifs`, motif) %>%
  distinct() %>% 
  summarise_all(median)

motif.family <- motif.family %>% 
  column_to_rownames(var = 'dap-seq_motifs')


# DAP only
tnt <- str_locate(motif.family$motif, '\\_') %>% 
  data.frame()

tnt$end <- (tnt$end)-1
tnt$start <- 1

motif.family$motif <- str_sub(motif.family$motif,
                              start = tnt$start,
                              end = tnt$end)

output <- data.frame(motif.family)
output <- rownames_to_column(output, var = 'dap-seq_motifs')
output <- output %>% 
  dplyr::rename(family = motif)

write.table(output,
            'f4/dap/motif_1000/meme.txt_0.01_DAP-seq_top_sim_dap.txt',
            row.names = F, quote = F, sep = '\t')