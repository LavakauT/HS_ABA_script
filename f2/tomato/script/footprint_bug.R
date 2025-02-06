library(dplyr)
library(tidyverse)


# dir <- '/RAID1/working/R425/lavakau/atac-seq/marchantia/20230626/ATACorrect_DEP/BINDetect/promoter/Tak1_HS_vs_Tak1_CK'
# dir <- '/RAID1/working/R425/lavakau/atac-seq/marchantia/20230626/ATACorrect_DEP/BINDetect_pd/promoter/hsfa_HS_vs_hsfa_CK'
dir <- '/RAID1/working/R425/lavakau/atac-seq/marchantia/20230626/ATACorrect_DEP/BINDetect_pd/promoter/dko_HS_vs_dko_CK'
motifs <- list.dirs(dir,
                    full.names = F,
                    recursive = F)

df <- data.frame()
for (motif in motifs) {
  print(paste0('Processing: ', motif))
  file <- read.delim(paste(dir, motif, paste0(motif, '_overview.txt'), sep = '/'))
  df <- rbind(df, file)
}

df %>% 
  distinct() %>% 
  write.table(., paste(dir, 'dko_overview.txt', sep = '/'),
              row.names = F, quote = F, sep = '\t')
save.image(paste(dir, 'dko_overview.RData', sep = '/'))

data <- read.delim('/Users/user/Desktop/overview.txt')
data <- read.delim('/Users/user/Desktop/dko_overview.txt')
data <- read.delim('/Users/user/Desktop/hsfa_overview.txt')

inf <- data[is.infinite(data$HS_CK_log2fc),]
# chr7.17191792.17191908 C12 showed INF because no bound in 
# No bound in HS and CK (especially in CK no score therefore INF called)

# chr3.27116535.27116973 C7 showed INF because no bound in 
# No bound in HS and CK (especially in CK no score therefore INF called)


# chr2.10458106.10458869 C9 showed INF because no bound in 
# No bound in HS and CK (especially in CK no score therefore INF called)

bed <- read.delim('/Users/user/Desktop/f1/counts/split/DE_peaks_hclust_available.bed',
                  header = F)
removal <- data.frame(V1 = c('chr7', 'chr3', 'chr2'),
                      V2 = c(17191792, 27116535, 10458106),
                      V3 = c(17191908, 27116973, 10458869),
                      V4 = c('C12', 'C7', 'C9'))
bed2 <- anti_join(bed, removal)
write.table(bed2, '/Users/user/Desktop/f1/counts/split/DE_peaks_hclust_available2.bed',
            row.names = F, col.names = F, quote = F, sep = '\t')  
