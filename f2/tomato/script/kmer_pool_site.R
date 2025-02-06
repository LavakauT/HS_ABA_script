##################
library(dplyr)
library(ggplot2)
library(ggforce)
library(scales)
library(zoo)
library(stringr)
library(tidyverse)
library(Biostrings)
library(ggforce)
library(universalmotif)
##################



###### kmer distribution: all ######
proximal_motif <- read.delim('/Users/user/Desktop/tomato_atac/kmer/sim/sly_acr_pro/proximal_pcc.txt')
distal_motif <- read.delim('/Users/user/Desktop/tomato_atac/kmer/sim/sly_acr_dis/distal_pcc.txt')


# if you want bound_quantile_supported, please apply supported == 'Yes'
input_list_1 <- c('ATAC_M82_1h_vs_ATAC_M82_0h', 'ATAC_M82_6h_vs_ATAC_M82_0h')
input_list_2 <- c('distal')
input_list_3 <- c('strict', 'soft')
input_list_4 <- c('up')
path <- '/Users/user/Desktop/tomato_atac/tfcomb'
p_path <- '/Users/user/Desktop/tomato_atac/tfcomb/distribution'
p_path2 <- '/Users/user/Desktop/tomato_atac/tfcomb/distribution/plot'



for (z in input_list_1) {
  for (j in input_list_2) {
    for (h in input_list_3) {
      for (v in input_list_4) {
        input1 <- z
        input2 <- j
        input3 <- h
        input4 <- v
        
        if(input2 == 'proximal'){
          motifs <- proximal_motif
        }else {
          motifs <- distal_motif
        }
        
        
        if(input1 %in% input_list_1 & input2 %in% input_list_2){
          print(paste0('Processing: ', input1, ' ',input4, ' ', input2, ' + ', input3))
          if(input3 == 'soft'){
            subdir <- 'bed'
          }else{
            subdir <- 'bed_strict'
          }
          if(!file.exists(paste(path,
                               subdir,
                               paste0(input1, '_', input2, '_', input4),
                               'selected_sig',
                               paste0(input1, '_', input2, '_', input4, '_sites.txt'), sep = '/'))){
            next
          }
          bound_data <- read.delim(paste(path,
                                         subdir,
                                         paste0(input1, '_', input2, '_', input4),
                                         'selected_sig',
                                         paste0(input1, '_', input2, '_', input4, '_sites.txt'), sep = '/'))
          
          ## select novel or similar motif
          bound_data2 <- bound_data %>% 
            filter(site1_name %in% motifs$kmers | site2_name %in% motifs$kmers)
          
          
          # add family name
          meme <- read_meme('/Users/user/ArkEasePro/pCRE/atac/ArabidopsisDAPv1.meme')
          names <- data.frame()
          for (i in 1:length(meme)) {
            sub <- data.frame(name = meme[[i]]@name)
            # append
            names <- rbind(names, sub)
          }
          # DAP only
          tnt <- str_locate(names$name, '\\_') %>% 
            data.frame()
          tnt$end <- (tnt$end)-1
          tnt$start <- 1
          
          names$name<- str_sub(names$name,
                               start = tnt$start,
                               end = tnt$end)
          family <- unique(names$name)
          
          bound_data2$site1_TF <- ''
          bound_data2$site2_TF <- ''
          bound_data3 <- data.frame()
          for (i in 1:length(family)) {
            sub <- bound_data2 %>% 
              filter(., grepl(paste0('_',family[i],'_'), site1_name)) %>% 
              mutate(site1_TF = family[i])
            
            # append 1
            bound_data3 <- rbind(bound_data3, sub)
          }
          
          bound_data4 <- data.frame()
          for (i in 1:length(family)) {
            sub <- bound_data3 %>% 
              filter(., grepl(paste0('_',family[i],'_'), site2_name)) %>% 
              mutate(site2_TF = family[i])
            
            # append 2
            bound_data4 <- rbind(bound_data4, sub)
          }
          
          
          # assign distance +/-
          bound_data4[which(bound_data4$feat_strand == '+' & bound_data4$relative_location == 'Upstream'),]$distance <- -bound_data4[which(bound_data4$feat_strand == '+' & bound_data4$relative_location == 'Upstream'),]$distance
          bound_data4[which(bound_data4$feat_strand == '-' & bound_data4$relative_location == 'Upstream'),]$distance <- -bound_data4[which(bound_data4$feat_strand == '-' & bound_data4$relative_location == 'Upstream'),]$distance
          
          
          write.table(bound_data4, paste0(p_path, '/',
                                          input3, '_', input2, '_', input1, '_', input4, '_pairs_distribution.txt'),
                      row.names = F, quote = F, sep = '\t')
          
          
          # plot
          plot_data <- bound_data4 %>% 
            dplyr::select(site1_TF,
                          site2_TF,
                          distance)
          
          if (input3 == 'strict'){
            load('/Users/user/Desktop/tomato_atac/bind/strict_mode_all.RData')
          }else{
            load('/Users/user/Desktop/tomato_atac/bind/soft_mode_all.RData')
          }
          
          
          family <- top10 %>% 
            filter(region == input2,
                   genotype == input1)
          kmer_family <- family %>% 
            dplyr::select(kmer_family) %>% 
            distinct() %>% 
            pull()
          
          
          
          
          dis.sum <- data.frame()
          for (m in kmer_family) {
            df <- plot_data %>% filter(site1_TF %in% m |
                                         site2_TF %in% m)
            
            # add 1000bp to make -1kb to 0.5kb in range(0,1500)
            df$Preferential_Position <- as.numeric(df$distance) + 1000
            df <- df %>%
              group_by(Preferential_Position) %>% 
              summarise(count = n())
            names(df) <- c('location', 'tar')
            
            
            window.tar <- data.frame(location = 1:1500)
            window.tar <- left_join(window.tar, df, by = 'location')
            window.tar[is.na(window.tar$tar),]$tar <- 0
            window.tar$tar <- window.tar$tar + 1 
            
            
            window <- window.tar
            
            
            ## culculation
            # all
            # bin:100; sliding window:25; median
            tar.win <- rollapply(window$tar, width = 100, mean, by = 25, partial = FALSE)
            beg.win <- rollapply(window$beg, width = 100, mean, by = 25, partial = FALSE)
            
            
            
            dis <- data.frame(bin = seq(from = 100, to = 1500, by = 25))
            dis$tar <- tar.win
            dis$beg <- beg.win
            
            
            dis <- dis %>% 
              mutate(zscore_tar = (tar - mean(tar))/sd(tar))
            dis <- dis %>%
              gather(., key = 'group', value = 'zscore', -bin, -tar) %>% 
              mutate(kmer_family = m)
            
            # append dis
            dis.sum <- rbind(dis.sum, dis)
          }
          
          write.table(dis.sum, paste0(p_path, '/',
                                      input3, '_', input2, '_', input1, '_', input4, '_pairs_zscore.txt'),
                      row.names = F, quote = F, sep = '\t')
          
          
          plot <- dis.sum %>% 
            ggplot(.,
                   aes(x = bin,
                       y = kmer_family,
                       fill = zscore)) +
            geom_tile() +
            scale_fill_gradient2(name = 'zscore',
                                 low = 'blue',
                                 mid = 'white',
                                 high = 'red',
                                 midpoint = 0,
                                 limits = c(-2,2),
                                 oob=squish) +
            theme_classic() +
            scale_x_continuous('Distance to TSS (kb)',
                               limits = c(0,1500),
                               breaks = seq(0, 1500, 500),
                               labels = c('-1.0','-0.5', 'TSS', '0.5')) +
            theme(strip.text.x = element_text(size = 10,
                                              face = 'bold.italic'),
                  strip.text.y = element_text(size = 10,
                                              face = 'bold',
                                              angle = 0),
                  strip.background = element_blank(),
                  axis.title = element_text(size = 10, face = 'bold'),
                  axis.text = element_text(size = 10, face = 'bold'))
          pdf(paste0(p_path2, '/',
                     input3, '_', input2, '_', input1, '_', input4, '_peak_zscore.pdf'),
              width = 4, height = 2)
          print(plot)
          dev.off()
        }
      }
    }
  }
}
####################################