############
library(dplyr)
library(tidyverse)
library(universalmotif)
library(ggplot2)
library(circlize)
library(grid)
############



input1 <- '/Users/user/Desktop/tomato_atac/kmer/sim/sly_acr_pro/proximal_pfm.meme'
input2 <- '/Users/user/Desktop/tomato_atac/kmer/sim/sly_acr_dis/distal_pfm.meme'
input3 <- '/Users/user/Documents/ATAC-seq/marchantia/tobias/ArabidopsisDAPv1.meme'

proximal <- read_meme(input1)
distal <- read_meme(input2)
dap <- read_meme(input3)


#####################
### motif compare ####
#####################

comparisons <- compare_motifs(c(proximal, dap),
                              method = "PCC",
                              min.mean.ic = 0,
                              score.strat = "a.mean")

# remove pCRE to pCRE PCC
comparisons2 <- comparisons[-c(1:length(proximal)), ]


for (i in 1:length(proximal)) {
  if(i == 1){
    
    df <- data.frame(comparisons2[,i])
    cname <- proximal[[i]]@name
    colnames(df) <- cname
    rname <- row.names(df)
    df$motif <- rname
    df <- df[order(-df[,1]), ] # top similar motif
    df <- df[1,c('motif', cname)]
    row.names(df) <- 1
    
  } else {
    
    df2 <- data.frame(comparisons2[,i])
    cname2 <- proximal[[i]]@name
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
  dplyr::rename(match = motif)
df3 <- df3 %>% filter(PCC != '')
df3$region <- 'proximal'

sequences <- data.frame()
for (i in 1:length(proximal)) {
  seqq <- proximal[[i]]@consensus
  cname <- proximal[[i]]@name
  sub <- data.frame(kmers = cname,
                    sequences = seqq)
  
  sequences <- rbind(sequences, sub)
}

df3 <- merge(df3, sequences)
df3$type <- 'PCC >= 0.9'
df3[which(df3$PCC < 0.9),]$type <- 'PCC < 0.9'
df3$remarks <- 'Similar to Arabidopsis DAP-seq'
df3[which(df3$type == 'PCC < 0.9'),]$remarks <- 'Novel to Tomato'
df3 <- df3 %>% 
  relocate(kmers, sequences, region, match, PCC, type, remarks)

write.table(df3, '/Users/user/Desktop/tomato_atac/kmer/sim/sly_acr_pro/proximal_pcc.txt',
            row.names = F, quote = F, sep = '\t')

gc()
#####################
### motif compare ####
#####################

comparisons <- compare_motifs(c(distal, dap),
                              method = "PCC",
                              min.mean.ic = 0,
                              score.strat = "a.mean")

# remove pCRE to pCRE PCC
comparisons2 <- comparisons[-c(1:length(distal)), ]


for (i in 1:length(distal)) {
  if(i == 1){
    
    df <- data.frame(comparisons2[,i])
    cname <- distal[[i]]@name
    colnames(df) <- cname
    rname <- row.names(df)
    df$motif <- rname
    df <- df[order(-df[,1]), ] # top similar motif
    df <- df[1,c('motif', cname)]
    row.names(df) <- 1
    
  } else {
    
    df2 <- data.frame(comparisons2[,i])
    cname2 <- distal[[i]]@name
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
  dplyr::rename(match = motif)
df3 <- df3 %>% filter(PCC != '')
df3$region <- 'distal'

sequences <- data.frame()
for (i in 1:length(distal)) {
  seqq <- distal[[i]]@consensus
  cname <- distal[[i]]@name
  sub <- data.frame(kmers = cname,
                    sequences = seqq)
  
  sequences <- rbind(sequences, sub)
}

df3 <- merge(df3, sequences)
df3$type <- 'PCC >= 0.9'
df3[which(df3$PCC < 0.9),]$type <- 'PCC < 0.9'
df3$remarks <- 'Similar to Arabidopsis DAP-seq'
df3[which(df3$type == 'PCC < 0.9'),]$remarks <- 'Novel to Tomato'
df3 <- df3 %>% 
  relocate(kmers, sequences, region, match, PCC, type, remarks)

write.table(df3, '/Users/user/Desktop/tomato_atac/kmer/sim/sly_acr_dis/distal_pcc.txt',
            row.names = F, quote = F, sep = '\t')
gc()
###########################################