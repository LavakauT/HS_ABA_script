############
# load packages
library(dplyr)
library(tidyverse)
library(ggplot2)
library(universalmotif)
library(scales)
library(circlize)
library(grid)
library(ComplexHeatmap)
library(Biostrings)
library(stringr)
library(GenomicFeatures)
library(ChIPseeker)
library(WGCNA)
library(DESeq2)
library(zoo)
library(vctrs)

allowWGCNAThreads()
############


############################################# Figure 3: A
############################################# Figure 7: A,B
############################################# Supporting Figure 3: A,B
############################################# Supporting Figure 4: A,B
############################################# Supporting Figure 5: F,G
###### ATAC-SEQ PEAKS: MARCHANTIA ######
# make TP bed
# Strand of each peak need to follow its annotated gene
# Peak name, Chr, Start, End, Strand
# setwd("~/folder_to_all_the_data")

p <- 'f1/counts/split'
p_b_pro <- 'f2/kmer/bed_promoter'
p_b_dis <- 'f2/kmer/bed_distal'
p_pro <- 'f2/kmer/peak_promoter'
p_dis <- 'f2/kmer/peak_distal'
# dir.create(p_b_pro)
# dir.create(p_b_dis)
# dir.create(p_pro)
# dir.create(p_dis)

de_peak <- read.table(paste(p, 'DE_peaks_hclust.txt', sep = '/'), header = T)
head(de_peak)

Txdb_gtf <- makeTxDbFromGFF('f1/MpTak_v6.1r1.gff')

for (i in paste0('C', 1:15)) {
  # extract Chr, Start, End, Peak name
  data <- de_peak %>% 
    filter(cluster == i) %>%
    dplyr::select(peak) %>%
    mutate(Peak = peak) %>% 
    separate(., col = 'peak', into = c('Chr', 'Start', 'End'), sep = '\\.') %>%
    mutate(Start = as.numeric(Start),
           End = as.numeric(End)) %>% 
    relocate(Chr, Start, End)
  
  # make GRanges data
  gr <- makeGRangesFromDataFrame(data, keep.extra.columns=T)
  
  # annotate
  peakAnno <- annotatePeak(gr,
                           tssRegion=c(-1500, 1500),
                           genomicAnnotationPriority = c("Exon", "Intron", "5UTR", "3UTR", "Promoter", "Downstream", "Intergenic"),
                           TxDb=Txdb_gtf,
                           level ='gene')
  
  # export promoter peaks annotation
  peakAnno_tb <- as_tibble(peakAnno@anno) %>%
    dplyr::select(Peak, seqnames, start, end, geneStrand) %>% 
    relocate(Peak, seqnames, start, end, geneStrand)
  
  # export proximal bed and peak
  setwd(p_b_pro)
  write.table(peakAnno_tb, paste0(i, '.bed'), row.names = F, col.names = F, quote = F, sep = '\t')
  setwd(p_pro)
  write.table(peakAnno_tb$Peak, paste0(i, '.txt'), row.names = F, col.names = F, quote = F, sep = '\t')
  
  
  # export distal peaks annotation
  peakAnno_tb <- as_tibble(peakAnno@anno) %>%
    filter(annotation == 'Distal Intergenic') %>% 
    dplyr::select(Peak, seqnames, start, end, geneStrand) %>% 
    relocate(Peak, seqnames, start, end, geneStrand)
  
  # export distal bed and peak
  setwd(p_b_dis)
  write.table(peakAnno_tb, paste0(i, '.bed'), row.names = F, col.names = F, quote = F, sep = '\t')
  setwd(p_dis)
  write.table(peakAnno_tb$Peak, paste0(i, '.txt'), row.names = F, col.names = F, quote = F, sep = '\t')
}


# make TN bed
# Peak name, Chr, Start, End, Strand
range <- data.frame(dir = 'f1/counts/split',
                    treatment = c(paste0(c('Tak1', 'hsfa', 'hsfb', 'dko'), '_HS'),
                                  paste0(c('hsfa', 'hsfb', 'dko'), '_HS')),
                    base = c(rep('Tak1_CK', 4),
                             paste0(c('hsfa', 'hsfb', 'dko'), '_CK')))

p <- 'f1/counts/split'
p_neg <-  paste(p, 'neg', sep = '/')
sum <- c()
for (i in 1:4) {
  if(i == 1){
    # neg
    neg <- read.table(paste(p_neg, paste0(range$treatment[i], '_', range$base[i], '_n.txt'), sep = '/'), header = T)
    
    sum.sub <- distinct(neg)
    sum <- c(sum, sum.sub) 
  }else{
    z <- i + 3
    # combine /Tak1 and /genotype
    # /Tak1
    # neg
    neg <- distinct(merge(read.table(paste(p_neg, paste0(range$treatment[i], '_', range$base[i], '_n.txt'), sep = '/'), header = T),
                          read.table(paste(p_neg, paste0(range$treatment[z], '_', range$base[z], '_n.txt'), sep = '/'), header = T)))
    
    sum.sub <- distinct(neg)
    sum <- merge(sum, sum.sub)
  }
}

TN_peaks <- distinct(sum) %>% 
  dplyr::rename(peak = x) %>% 
  mutate(cluster = 'TN')
# There are 10092 TN peaks originally.
# To prevent error in downstream analysis (not found peak in GTF file), we eliminated 3 peaks.
# Therefore, there are 10089 TN peaks.


head(TN_peaks)
setwd(p_bed)
data <- TN_peaks %>% 
  dplyr::select(peak) %>%
  mutate(Peak = peak) %>% 
  separate(., col = 'peak', into = c('Chr', 'Start', 'End'), sep = '\\.') %>%
  mutate(Start = as.numeric(Start),
         End = as.numeric(End)) %>% 
  relocate(Chr, Start, End)

# export
write.table(data, 'TN_peaks.txt', row.names = F, quote = F, sep = '\t')


# make balanced negative files
range <- paste0('C', 1:15)
p_peak <- 'f2/kmer/peak'
p_pro <- 'f2/kmer/peak_promoter'
p_dis <- 'f2/kmer/peak_distal'

setwd(p_peak)
df <- read.table('TN.txt') 

dirs <- c(p_pro, p_dis)
for (dir in dirs) {
  setwd(dir)
  for (z in range) {
    cluster <- z
    neg <- df
    n.neg <- nrow(neg)
    pos <- read.table(paste0(cluster, '.txt'))
    n.pos <- nrow(pos)
    
    if (!dir.exists(cluster)){
      dir.create(cluster)
    }
    
    for (i in 1:10){
      number <- i
      set.seed(number)
      random.neg <- sample(seq_len(n.neg), size = n.pos)
      data <- neg[random.neg,]
      
      # output result
      write.table(data, 
                  paste(cluster, paste0(paste('neg',
                                              cluster,
                                              number,
                                              sep = '_'), '.txt'), sep = '/'),
                  row.names = F, col.names = F, quote = F, sep = '\t')
    }
  }
}
##################


###### ATAC-SEQ GENES: MARCHANTIA ######
p_gene <- 'f3'
dir.create(paste(p_gene, 'peak2gene', sep = '/'))

# TP
tp <- read.delim('f1/counts/split/peak2gene.txt') %>% 
  distinct()
length(unique(tp$geneId)) # 4727
clusters <- unique(tp$cluster) # 15 clusters

for (cl in clusters) {
  # make directory for TP
  dir.create(paste(p_gene, 'peak2gene', cl, sep = '/'))
  
  # export tp gene file in each cluster
  tp %>% 
    filter(cluster == cl) %>% 
    dplyr::select(geneId) %>% 
    write.table(., paste(p_gene, 'peak2gene', paste0(cl, '.txt'), sep = '/'),
                row.names = F, col.names = F, quote = F, sep = '\t')
}



# TN
# extract Chr, Start, End, Peak name
Txdb_gtf <- makeTxDbFromGFF('f1/MpTak_v6.1r1.gff')
tn <- read.delim('f1/counts/split/TN_peaks.txt')
data <- tn %>% 
  filter(cluster == 'TN') %>% 
  dplyr::select(peak) %>%
  mutate(Peak = peak) %>% 
  separate(., col = 'peak', into = c('Chr', 'Start', 'End'), sep = '\\.') %>%
  mutate(Start = as.numeric(Start),
         End = as.numeric(End)) %>% 
  relocate(Chr, Start, End)

# make GRanges data
gr <- makeGRangesFromDataFrame(data, keep.extra.columns=T)

# annotate
peakAnno <- annotatePeak(gr,
                         tssRegion=c(-1500, 1500),
                         genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
                         TxDb=Txdb_gtf,
                         level ='gene')

# export peak genes
# remove genes in TP
peakAnno_tb <- as_tibble(peakAnno@anno) %>%
  dplyr::select(geneId) %>%
  filter(!(geneId %in% unique(tp$geneId))) %>% 
  distinct() # 5207
# export
write.table(peakAnno_tb, 'f1/counts/split/peak2gene_TN.txt',
            row.names = F, quote = F, sep = '\t')

tn <- read.delim('f1/counts/split/peak2gene_TN.txt')

# make balanced negative files
range <- clusters
setwd(p_gene)
df <- tn
setwd(paste(p_gene, 'peak2gene', sep = '/'))
for (z in range) {
  cluster <- z
  neg <- df
  n.neg <- nrow(neg)
  pos <- read.table(paste0(cluster, '.txt'))
  n.pos <- nrow(pos)
  
  if (!dir.exists(cluster)){
    dir.create(cluster)
  }
  
  for (i in 1:10){
    number <- i
    set.seed(number)
    random.neg <- sample(seq_len(n.neg), size = n.pos)
    data <- neg[random.neg,]
    
    # output result
    write.table(data, 
                paste(cluster, paste0(paste('neg',
                                            cluster,
                                            number,
                                            sep = '_'), '.txt'), sep = '/'),
                row.names = F, col.names = F, quote = F, sep = '\t')
  }
}
# setwd("~/folder_to_all_the_data")
##################


###### ATAC-SEQ PEAKS: SUMMARY FEATURES COUNTS ######
# setwd("~/folder_to_all_the_data")
p_fea <- 'f2/kmer/feature'
dirs <- list.dirs(p_fea,
                  full.names = F,
                  recursive = F) # proximal (promoter) or distal
# 'SVM', 'RF', 'LogReg': algorithms of machine learning
# alg <- c('SVM', 'RF', 'LogReg') 
clusters <- paste0('C', 1:15) # all clusters in ATAC-seq data

sum <- data.frame()
for (dir in dirs) {
  # for (a in alg) {
    for (c in clusters) {
      df <- read.table(paste(p_fea, dir, #a,
                             paste(c, #a,
                                   'kmer.txt', sep = '_'), sep = '/'),
                       header = T) %>% 
        dplyr::select(-file)
      
      colname <- colnames(df)
      
      # head(colname)
      
      len_colname <- length(colname)
      
      # head(len_colname)
      
      sum.sub <- data.frame(type = dir,
                            # algorithm = a,
                            cluster = c,
                            kmer_count = len_colname)
      sum <- rbind(sum, sum.sub)
    }
  # }
}
head(sum)

sum$type <- factor(sum$type, levels = c('mp_acr_promoter', 'mp_acr_distal'),
                   labels = c('PROXIMAL', 'DISTAL'))
sum$cluster <- factor(sum$cluster, levels = paste0('C', 1:15))

plt <- sum %>%
  filter(type == "PROXIMAL") %>% 
  ggplot(.,
         aes(x = cluster,
             y = kmer_count,
             fill = type)) +
  geom_bar(stat = 'identity',
           position = 'dodge2',
           width = 0.8) +
  coord_flip() +
  labs(x = 'Cluster',
       y = 'Numbers of kmers',
       fill = 'Type') +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 10, face = 'bold')) # Figure 7, k-mers [counts] in peaks


## Figure 7-A, k-mers [counts] in peaks
pdf('f2/plots/feature_counts.pdf',
    width = 5, height = 4)
plt 
dev.off()
##################


###### ATAC-SEQ GENES: SUMMARY FEATURES COUNTS ######
# features
p_feat <- "f3/feature/mp_acr_gene"

feat <- data.frame()
files <- list.files(p_feat,
                    full.names = F)
for (file in files) {
  file_name <- str_remove(file, '_kmer.txt')
  df <- read.table(paste(p_feat, file, sep = '/'), header = T)
  
  if(colnames(df)[1] == 'file'){
    kmer_n <- ncol(df)-1
  }
  
  # append numbers of feature
  feat.sub <- data.frame(type = 'mp_rna',
                         cluster = file_name,
                         feature = kmer_n)
  
  # append numbers of feature
  feat <- rbind(feat, feat.sub)
}

# visualize feature counts
p_plot <- 'f2/plot/'
feat$cluster <- factor(feat$cluster, levels = paste0('C', 15:1))


## Figure 7-A, k-mers [counts] in genes
pdf('f2/plots/feature_counts2.pdf',
    width = 5, height = 4)
feat %>% 
  ggplot(.,
         aes(x = cluster,
             y = feature,
             fill = factor(type, levels = 'type', label = 'DEGs'))) +
  geom_bar(stat = 'identity',
           position = 'dodge',
           width = 0.5) +
  geom_text(aes(x = cluster,
                y = as.numeric(feature)*1.2,
                label = feature),
            color = 'red') +
  theme_classic() +
  coord_flip() +
  labs(x = 'Cluster',
       y = 'Feature',
       fill = 'Type') +
  scale_fill_manual(values = c('grey','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e')) +
  facet_grid(cluster~., scales = 'free_y') +
  theme(strip.background = element_blank(),
        strip.text = element_blank())
dev.off()
##################


###### ATAC-SEQ PEAKS: SUMMARY TRAIN RESULTS ######
# TRAINING PARAMETER: BalancedRuns 10 with CVfold 10
p_train <- 'f2'
results <- read.delim('f2/RESULTS_1.txt')
head(results)

# remove prefix and tail in column ID
results$ID <- str_remove(results$ID, '/RAID1/working/R425/lavakau/ML_master/')
results$ID <- str_remove(results$ID, '\\.pcre_df_p0.01_mod.txt_SVM|\\.pcre_df_p0.01_mod.txt_RF|\\.pcre_df_p0.01_mod.txt_LogReg')
head(results)

# separate ID columns as train type, cluster, and numbers of negative files
results <- results %>% 
  separate(., col = 'ID', into = c('type', 'cluster', 'negative_data'), sep = '/')
results$cluster <- factor(results$cluster, levels = paste0('C', 1:15))
results$type <- factor(results$type, levels = c('mp_acr', 'mp_acr_strict', 'mp_acr_promoter', 'mp_acr_distal'),
                       labels = c('ACR w/ TN1 (deprecated)', 'ACR w/ TN2 (deprecated)', 'PROXIMAL', 'DISTAL'))

plt <- ggplot(results %>% 
                group_by(type, cluster, Alg) %>% 
                summarise(F1_test_mean = mean(F1_test, na.rm = T),
                          AUCROC_test_mean = mean(AUCROC_test, na.rm = T)) %>% 
                filter(type == "PROXIMAL",
                       Alg == "RF",
                       ),
              aes(x = cluster,
                  y = type,
                  fill = F1_test_mean)) +
  geom_tile() +
  labs(title = expression(bold('Machine Learning Test Results')),
       x = 'Cluster',
       y = 'Type',
       fill = 'F1') + # change y axis F1/AUC
  coord_flip() +
  theme_classic() +
  scale_fill_gradientn(
    colours = c("white", "red"),         # Colors for different ranges
    limits = c(0.5, 1)                             # Set data range
  ) +
  theme(axis.text = element_text(size = 10, face = 'bold'),
        axis.title = element_text(size = 12, face = 'bold'))


## Figure 7-A, AUCROC of peaks
# results of cluster 8, 12, 14, 15 were eliminated due to too few kmers
# pdf('f2/plots/f1.pdf',
#     width = 3, height = 5)
# plt
# dev.off()


plt <- ggplot(results %>% 
                group_by(type, cluster, Alg) %>% 
                summarise(F1_test_mean = mean(F1_test, na.rm = T),
                          AUCROC_test_mean = mean(AUCROC_test, na.rm = T)) %>% 
                filter(type == "PROXIMAL",
                       Alg == "RF"),
              aes(x = cluster,
                  y = type,
                  fill = AUCROC_test_mean)) +
  geom_tile() +
  labs(title = expression(bold('Machine Learning Test Results')),
       x = 'Cluster',
       y = 'Type',
       fill = 'AUCROC') + # change y axis F1/AUC
  coord_flip() +
  theme_classic() +
  scale_fill_gradientn(
    colours = c("white", "red"),         # Colors for different ranges
    limits = c(0.5, 1)                             # Set data range
  ) +
  theme(axis.text = element_text(size = 10, face = 'bold'),
        axis.title = element_text(size = 12, face = 'bold'))


## Figure 7-A, AUCROC of peaks
# results of cluster 8, 12, 14, 15 were eliminated due to too few kmers
pdf('f2/plots/aucroc.pdf',
    width = 3, height = 5)
plt
dev.off()
##################


###### ATAC-SEQ GENES: SUMMARY TRAIN RESULTS ######
# BalancedRuns 10 CVfold 10
p_train <- 'f3'
setwd(p_train)
results <- read.delim('RESULTS.txt')
head(results)

# remove prefix and tail in column ID
results$ID <- str_remove(results$ID, '/RAID1/working/R425/lavakau/ML_master/')
results$ID <- str_remove(results$ID, '\\.fa.pcre_df_p0.01_mod.txt_SVM_lasso|\\.fa.pcre_df_p0.01_mod.txt_RF_lasso|\\.fa.pcre_df_p0.01_mod.txt_LogReg_lasso')
head(results)

# separate ID columns as train type, cluster, and numbers of negative files
results <- results %>% 
  separate(., col = 'ID', into = c('type', 'cluster', 'negative_data'), sep = '/')
results$cluster <- factor(results$cluster, levels = paste0('C', 1:15))

auc <- results %>% 
  dplyr::select(cluster, negative_data, Alg, AUCROC_test, Tag, type) %>% 
  group_by(cluster, Alg, Tag, type) %>% 
  summarise(mean_AUCROC_test = mean(AUCROC_test),
            sd_AUCROC_test = sd(AUCROC_test)) %>% 
  mutate(lab = paste0(round(mean_AUCROC_test, 2), '±', round(sd_AUCROC_test, 3)))

auc$type <- factor(auc$type, levels = c('mp_acr_gene', 'mp_rna_ovl'),
                   labels = c('Genes', 'OCR overlapping RNA genes'))
auc$Tag <- factor(auc$Tag, levels = c('lasso', 'no_lasso'),
                  labels = c('w/ Lasso', 'w/o Lasso'))


plt <- ggplot(auc %>% 
                filter(Alg == "RF",
                       type == "Genes",
                       Tag == "w/o Lasso"),
              aes(x = cluster,
                  y = type,
                  fill = mean_AUCROC_test)) +
  geom_tile() +
  labs(title = expression(bold('Machine Learning Test Results')),
       x = 'Cluster',
       y = 'type',
       fill = 'AUCROC') + # change y axis F1/AUC
  coord_flip() +
  theme_classic() +
  scale_fill_gradientn(
    colours = c("white", "red"),         # Colors for different ranges
    limits = c(0.5, 1)                             # Set data range
  ) +
  theme(axis.text = element_text(size = 10, face = 'bold'),
        axis.title = element_text(size = 12, face = 'bold'))
plt


## Figure 7-A, AUCROC of genes
pdf('f2/plots/aucroc2.pdf',
    width = 3, height = 5)
plt
dev.off()
################################


###### SUMMARY RESULTS FROM KMER BINDETECT BY TOBIAS ######
# We only discussed differential binding between HS and Cntr (CK) under the same genotype.
# Therefore, the input1 can be 'Tak1_HS_vs_Tak1_CK', 'hsfa_HS_vs_hsfa_CK', 'hsfb_HS_vs_hsfb_CK', 'dko_HS_vs_dko_CK'.
# And the input2 can be promoter (proximal) or distal.
# differential binding events in 95% or 5% quantile
input_list_1 <- c('Tak1_HS_vs_Tak1_CK', 'hsfa_HS_vs_hsfa_CK', 'hsfb_HS_vs_hsfb_CK', 'dko_HS_vs_dko_CK')
input_list_2 <- c('promoter', 'distal')


# results vectors
full <- data.frame() # all events
top10 <- data.frame() # top 10 kmer family in each genotype


for (z in input_list_1) {
  for (j in input_list_2) {
    input1 <- z
    input2 <- j
    
    if(input1 %in% input_list_1 & input2 %in% input_list_2){
      print(paste0('Processing: ', input1, ' under ', input2))
      path <- paste('f2/bind/', input2, input1, sep = '/')
      bound_data <- read.delim(paste(path, 'bindetect_results.txt', sep = '/'))
      
      # remove no bound under HS & CK
      bound_data2 <- bound_data[which(bound_data$HS_bound != 0 | bound_data$CK_bound != 0),]
      
      # remove prefix of TFBS_name
      bound_data2$kmer_motif_name <- str_remove(bound_data2$kmer_motif_name, '_')
      
      
      
      # add family name
      meme <- read_meme('f2/ArabidopsisDAPv1.meme')
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
      
      bound_data2$kmer_family <- ''
      bound_data3 <- data.frame()
      for (i in 1:length(family)) {
        sub <- bound_data2 %>% 
          filter(., grepl(paste0('_',family[i],'_'), kmer_motif_name)) %>% 
          mutate(kmer_family = family[i])
        
        # append
        bound_data3 <- rbind(bound_data3, sub)
      }
      
      # check no data loss
      # calculate all differential binding kmers to family
      # 282 kmers in 54 families
      sig_kmer <- read.delim(paste(path, paste0(input1, '_de_kmer.txt'), sep = '/')) %>% 
        dplyr::select(output_prefix, motif_id, total_tfbs, HS_bound, CK_bound, HS_CK_change) %>% 
        dplyr::rename(TFBS_name = output_prefix,
                      kmer_motif_name = motif_id,
                      total_HS_bound = HS_bound,
                      total_CK_bound = CK_bound)
      
      sum_sig_kmer <- bound_data3 %>% 
        dplyr::select(kmer_motif_name, kmer_family) %>% 
        distinct() %>% 
        group_by(kmer_family) %>%
        summarise(total_kmer = n())
      
      
      # calculation: which cluster conquer the differential of a specific kmer
      # up bound: HS_CK_log2fc > 0
      # total bound: sum(HS_bound); sum(CK_bound)
      bound_counts <- bound_data3 %>% 
        group_by(TFBS_name, peak_cluster) %>%
        filter(HS_CK_log2fc > 0) %>%  
        mutate(HS_bound_sum = sum(HS_bound),
               CK_bound_sum = sum(CK_bound)) %>% 
        dplyr::select(TFBS_name, peak_cluster, kmer_motif_name, kmer_family, HS_bound_sum, CK_bound_sum) %>% 
        distinct()
      
      # bound median log2fc
      bound_median_fc <- bound_data3 %>% 
        group_by(TFBS_name, peak_cluster) %>%
        filter(HS_CK_log2fc > 0) %>%
        mutate(median_HS_CK_log2fc = median(HS_CK_log2fc)) %>% 
        dplyr::select(TFBS_name, peak_cluster, kmer_motif_name, kmer_family, median_HS_CK_log2fc) %>% 
        distinct()
      
      # bound profile
      bound_profile <- merge(bound_counts, bound_median_fc)%>% 
        merge(., sig_kmer %>% dplyr::select(-TFBS_name)) %>% 
        mutate(hg = phyper(HS_bound_sum-1, total_HS_bound, total_tfbs-total_HS_bound, HS_bound_sum + CK_bound_sum, lower.tail = F)) %>% 
        mutate(supported = 'No')
      
      bound_profile[which(bound_profile$hg < 0.05),]$supported <- 'Yes'
      
      
      # we chose quantile 95, 85, 75, 50 and left
      qun <- data.frame(log2fc = quantile(bound_profile$median_HS_CK_log2fc,
                                          probs = c(0.5, 0.75, 0.85, 0.95))) %>% 
        rownames_to_column(var = 'quantile')
      
      # add quantile information
      bound_profile$quantile <- '<50%'
      for (i in 1:4) {
        cut_off <- qun[i,]$log2fc
        mark <- qun[i,]$quantile
        
        bound_profile[which(bound_profile$median_HS_CK_log2fc >= cut_off),]$quantile <- mark
      }
      
      # check which cut-off i want
      table(bound_profile$supported, bound_profile$quantile)
      
      # i chose 95% first
      bound_profile_95 <- bound_profile %>% 
        filter(quantile == '95%')
      
      # summary bound under 95%
      num_95 <- bound_profile_95 %>% 
        group_by(peak_cluster, kmer_family) %>% 
        summarise(number_of_bound_kmer = n(),
                  median_bound_log2fc = median(median_HS_CK_log2fc),
                  kmer_HS_bound_sum = sum(HS_bound_sum),
                  kmer_CK_bound_sum = sum(CK_bound_sum),
                  peak_total_HS_bound = sum(total_HS_bound),
                  peak_total_CK_bound = sum(total_CK_bound)) %>% 
        mutate(bound_type = 'positive',
               genotype = input1,
               region = input2,
               ratio_kmer_HS_bound_sum = kmer_HS_bound_sum/peak_total_HS_bound,
               ratio_kmer_CK_bound_sum = kmer_CK_bound_sum/peak_total_CK_bound)
      
      
      # down bound: HS_CK_log2fc > 0
      # total bound: sum(HS_bound); sum(CK_bound)
      bound_counts2 <- bound_data3 %>% 
        group_by(TFBS_name, peak_cluster) %>%
        filter(HS_CK_log2fc < 0) %>%  
        mutate(HS_bound_sum = sum(HS_bound),
               CK_bound_sum = sum(CK_bound)) %>% 
        dplyr::select(TFBS_name, peak_cluster, kmer_motif_name, kmer_family, HS_bound_sum, CK_bound_sum) %>% 
        distinct()
      
      # bound median log2fc
      bound_median_fc2 <- bound_data3 %>% 
        group_by(TFBS_name, peak_cluster) %>%
        filter(HS_CK_log2fc < 0) %>%
        mutate(median_HS_CK_log2fc = median(HS_CK_log2fc)) %>% 
        dplyr::select(TFBS_name, peak_cluster, kmer_motif_name, kmer_family, median_HS_CK_log2fc) %>% 
        distinct()
      
      # bound profile
      bound_profile2 <- merge(bound_counts2, bound_median_fc2) %>%
        merge(., sig_kmer %>% dplyr::select(-TFBS_name)) %>%
        mutate(hg = phyper(CK_bound_sum-1, total_CK_bound, total_tfbs-total_CK_bound, HS_bound_sum + CK_bound_sum, lower.tail = F)) %>% 
        mutate(supported = 'No')
      
      bound_profile2[which(bound_profile2$hg < 0.05),]$supported <- 'Yes'
      
      # we chose quantile 5, 15, 25, 50 and left
      qun2 <- data.frame(log2fc = quantile(bound_profile2$median_HS_CK_log2fc,
                                           probs = c(0.05, 0.15, 0.25, 0.5))) %>% 
        rownames_to_column(var = 'quantile')
      
      # add quantile information
      bound_profile2$quantile <- '>50%'
      for (i in 4:1) {
        cut_off <- qun2[i,]$log2fc
        mark <- qun2[i,]$quantile
        
        bound_profile2[which(bound_profile2$median_HS_CK_log2fc <= cut_off),]$quantile <- mark
      }
      
      # check which cut-off i want
      # table(bound_profile2$supported, bound_profile2$quantile)
      
      # i chose 5% first
      bound_profile_5 <- bound_profile2 %>% 
        filter(quantile == '5%')
      
      # summary bound under 5%
      num_5 <- bound_profile_5 %>% 
        group_by(peak_cluster, kmer_family) %>% 
        summarise(number_of_bound_kmer = n(),
                  median_bound_log2fc = median(median_HS_CK_log2fc),
                  kmer_HS_bound_sum = sum(HS_bound_sum),
                  kmer_CK_bound_sum = sum(CK_bound_sum),
                  peak_total_HS_bound = sum(total_HS_bound),
                  peak_total_CK_bound = sum(total_CK_bound)) %>% 
        mutate(bound_type = 'negative',
               genotype = input1,
               region = input2,
               ratio_kmer_HS_bound_sum = kmer_HS_bound_sum/peak_total_HS_bound,
               ratio_kmer_CK_bound_sum = kmer_CK_bound_sum/peak_total_CK_bound)
      
      
      # combine positive and negative bound type
      # make new TF Family name
      bound_profile_combine <- rbind(num_5, num_95) %>% 
        distinct() %>% 
        merge(., sum_sig_kmer) %>% 
        mutate(bound_kmer_ratio = number_of_bound_kmer/total_kmer,
               kmer_name = paste0(kmer_family, ' (', total_kmer, ')'))
      
      # add missing cluster
      model <- data.frame(peak_cluster = paste0('C', 1:15),
                          kmer_name = rep(unique(bound_profile_combine$kmer_name), each = 30),
                          bound_type = rep(c('positive', 'negative'), each = 15))
      bound_profile_combine <- full_join(model, bound_profile_combine)
      
      
      # add factor to cluster
      bound_profile_combine$peak_cluster <- factor(bound_profile_combine$peak_cluster, levels = paste0('C', 1:15))
      bound_profile_combine$bound_type <- factor(bound_profile_combine$bound_type, levels = c('positive', 'negative'),
                                                 labels = c('P', 'N'))
      
      
      # select top10
      kmer_top <- bound_profile_combine %>%
        filter(!is.na(kmer_family)) %>% 
        arrange(desc(ratio_kmer_HS_bound_sum), desc(median_bound_log2fc)) %>%  # Arrange by descending kmer_bound_value
        slice_head(n = 10) %>% 
        pull(kmer_family)
      
      bound_profile_combine_top <- bound_profile_combine %>% 
        filter(kmer_family %in% kmer_top)
      
      # append the data
      full <- rbind(full, bound_profile_combine)
      top10 <- rbind(top10, bound_profile_combine_top)
      
      
      # pdf(paste0('plots/', input1, '_', input2, '_bound.pdf'),
      #     width = 10, height = 7)
      # plt <- bound_profile_combine %>%
      #   ggplot(.,
      #          aes(x = bound_type,
      #              y = kmer_name,
      #              color = median_bound_log2fc,
      #              size = bound_kmer_ratio)) +
      #   geom_point() +
      #   scale_color_gradient(low = 'blue',
      #                        high = 'red') +
      #   labs(x = 'Cluster',
      #        y = 'TF Family (kmer)',
      #        color = 'Median HS/CK log2FC',
      #        size = 'DE bound kmer/total kmer') +
      #   scale_x_discrete(position = "top") +
      #   facet_grid(.~peak_cluster) +
      #   theme_classic() +
      #   theme(axis.title = element_text(size = 12, face = 'bold'),
      #         axis.text = element_text(size = 10, face = 'bold'),
      #         strip.background = element_blank(),
      #         strip.text = element_text(size = 12, face = 'bold'),
      #         strip.placement = 'outside')
      # 
      # grid.draw(plt)
      # dev.off()
      
      
      # pdf(paste0('plots/bound_quantile_support/', input1, '_', input2, '_bound.pdf'),
      #     width = 10, height = 7)
      # plt <- bound_profile_combine %>%
      #   ggplot(.,
      #          aes(x = bound_type,
      #              y = kmer_name,
      #              color = median_bound_log2fc,
      #              size = bound_kmer_ratio)) +
      #   geom_point() +
      #   scale_color_gradient(low = 'blue',
      #                        high = 'red') +
      #   labs(x = 'Cluster',
      #        y = 'TF Family (kmer)',
      #        color = 'Median HS/CK log2FC',
      #        size = 'DE bound kmer/total kmer') +
      #   scale_x_discrete(position = "top") +
      #   facet_grid(.~peak_cluster) +
      #   theme_classic() +
      #   theme(axis.title = element_text(size = 12, face = 'bold'),
      #         axis.text = element_text(size = 10, face = 'bold'),
      #         strip.background = element_blank(),
      #         strip.text = element_text(size = 12, face = 'bold'),
      #         strip.placement = 'outside')
      # 
      # grid.draw(plt)
      # dev.off()
      
      
      
      # export 95% and 5% bound bed
      # you can check whether the output is the same with bound_profile_95 or bound_profile_5
       bound_data3 %>%
         filter(HS_CK_log2fc > 0) %>%
         merge(bound_profile_95 %>% dplyr::select(kmer_motif_name, kmer_family, peak_cluster), .) %>%
         dplyr::select(TFBS_chr, TFBS_start, TFBS_end, kmer_motif_name, HS_CK_log2fc, TFBS_strand) %>%
         write.table(., paste0('f2/tfcomb/bed/', input1, '_', input2, '_up.bed'),
                     row.names = F, col.names = F, quote = F, sep = '\t')


       bound_data3 %>%
         filter(HS_CK_log2fc < 0) %>%
         merge(bound_profile_5 %>% dplyr::select(kmer_motif_name, kmer_family, peak_cluster), .) %>%
         dplyr::select(TFBS_chr, TFBS_start, TFBS_end, kmer_motif_name, HS_CK_log2fc, TFBS_strand) %>%
         write.table(., paste0('f2/tfcomb/bed/', input1, '_', input2, '_down.bed'),
                     row.names = F, col.names = F, quote = F, sep = '\t')
    }
  }
}

gc()
save(top10, full, file = 'f2/strict_mode_all.RData')



# visualization
top10$genotype <- factor(top10$genotype, levels = input_list_1,
                         labels = c('Tak-1', 'hsfa1', 'hsfb1', 'dko'))


## Supporting Figure 3-A
pdf(paste0('f2/plots/', 'proximal', '_bound_soft_mode.pdf'),
    width = 12, height = 8)
plt <- top10 %>%
  filter(region == 'promoter') %>%
  ggplot(.,
         aes(x = bound_type,
             y = kmer_family,
             color = median_bound_log2fc,
             size = ratio_kmer_HS_bound_sum)) +
  geom_point() +
  scale_color_gradient2(low = 'blue',
                        mid = 'white',
                        high = 'red',
                        midpoint = 0,
                        limits = c(-2,2),
                        oob = scales::squish) +
  labs(x = 'Cluster',
       y = 'TF Family (kmer)',
       color = 'Median HS/CK log2FC',
       size = 'DE bound under HS/total bound') +
  scale_x_discrete(position = "top") +
  facet_grid(genotype~peak_cluster, scales = 'free_y', switch = 'y') +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 10, face = 'bold'),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12, face = 'bold'),
        strip.text.y.left = element_text(size = 12, face = 'bold', angle = 0),
        strip.placement = 'outside')
grid.draw(plt)
dev.off()


top10$peak_cluster_num <- str_remove(top10$peak_cluster, 'C')
data <- top10 %>% 
  filter(peak_cluster %in% paste0('C', 1:15))
data$peak_cluster_rank <- vec_rank(as.numeric(data$peak_cluster_num),
                                   ties = "dense")

plt <- data %>%
  filter(region == 'promoter',
         !is.na(kmer_family)) %>%
  ggplot(.,
         aes(x = genotype,
             y = peak_cluster_rank,
             color = median_bound_log2fc,
             size = ratio_kmer_HS_bound_sum)) +
  # Make custom panel grid
  geom_hline(
    aes(yintercept = peak_cluster_rank), 
    color = "lightgrey"
  ) +
  geom_point() +
  scale_color_gradient2(low = 'blue',
                        mid = 'white',
                        high = 'red',
                        midpoint = 0,
                        limits = c(-1,1),
                        oob = scales::squish) +
  labs(x = 'Genotype',
       y = 'Peak cluster',
       color = 'Median HS/CK log2FC',
       size = 'DE bound under HS/total bound') +
  scale_x_discrete(position = "top") +
  facet_wrap(.~kmer_family, ncol = 6) +
  coord_polar() +
  # theme_classic() +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 10, face = 'bold'),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12, face = 'bold'),
        strip.text.y.left = element_text(size = 12, face = 'bold', angle = 0),
        strip.placement = 'outside')

plt <- plt +
  # Annotate custom scale inside plot
  annotate(
    x = Inf, 
    y = 1, 
    label = "2", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 2, 
    label = "3", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 3, 
    label = "4", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 4, 
    label = "6", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 5, 
    label = "7", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 6, 
    label = "9", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 7, 
    label = "10", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  # Scale y axis so bars don't start in the center
  scale_y_continuous(
    limits = c(-1,10),
    expand = c(0, 0),
    breaks = c(1,2,3,4,5,6,7)
  ) +
  theme(
    # Make the background white and remove extra grid lines
    panel.background = element_rect(fill = "white", color = "white"),
    # panel.grid = element_blank(),
    panel.grid.major.x = element_line(color = 'grey', linetype = 'dashed', linewidth = 0.5),
    legend.position = "bottom",
    legend.direction = 'horizontal'
  )

plt

# manipulating kmer families shared among genotypes
# kmer families specific to tak1, hsfa, hsfb
# share
sh <- c('AP2EREBP', 'BES1', 'bZIP', 'C2C2gata', 'C2C2COlike', 'C3H', 'HMG', 'HSF', 'LOBAS2', 'NAC', 'Trihelix', 'WRKY')
# specific to a genotype
t <- c('SBP', 'bHLH', 'DBP', 'E2FDP')
a <- c('ABI3VP1', 'C2C2YABBY')
b <- c('MYBrelated')
d <- c('GRF', 'Homeobox')

data$sp <- ''
data[which(data$kmer_family %in% sh),]$sp <- 'Shared'
data[which(data$kmer_family %in% t),]$sp <- 'Tak1'
data[which(data$kmer_family %in% a),]$sp <- 'hsfa'
data[which(data$kmer_family %in% b),]$sp <- 'hsfb'
data[which(data$kmer_family %in% d),]$sp <- 'dko'
data$sp <- factor(data$sp, levels = c('Shared', 'Tak1', 'hsfa', 'hsfb', 'dko'))

plt_sub <- data %>%
  filter(region == 'promoter',
         !is.na(kmer_family)) %>%
  ggplot(.,
         aes(x = genotype,
             y = peak_cluster_rank,
             color = median_bound_log2fc,
             size = ratio_kmer_HS_bound_sum)) +
  # Make custom panel grid
  geom_hline(
    aes(yintercept = peak_cluster_rank), 
    color = "lightgrey"
  ) +
  geom_point() +
  scale_color_gradient2(low = 'blue',
                        mid = 'white',
                        high = 'red',
                        midpoint = 0,
                        limits = c(-1,1),
                        oob = scales::squish) +
  labs(x = 'Genotype',
       y = 'Peak cluster',
       color = 'Median HS/CK log2FC',
       size = 'DE bound under HS/total bound') +
  scale_x_discrete(position = "top") +
  facet_wrap(~sp + kmer_family, ncol = 6) +
  coord_polar() +
  # theme_classic() +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 10, face = 'bold'),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12, face = 'bold'),
        strip.text.y.left = element_text(size = 12, face = 'bold', angle = 0),
        # strip.placement = 'outside'
  )

plt_sub <- plt_sub +
  # Annotate custom scale inside plot
  annotate(
    x = Inf, 
    y = 1, 
    label = "", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 2, 
    label = "2", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 3, 
    label = "3", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 4, 
    label = "4", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 5, 
    label = "5", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 6, 
    label = "6", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 7, 
    label = "7", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 8, 
    label = "8", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 9, 
    label = "9", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 10, 
    label = "10", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 11, 
    label = "11", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 12, 
    label = "12", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 13, 
    label = "13", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 14, 
    label = "14", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 15, 
    label = "15", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  # Scale y axis so bars don't start in the center
  scale_y_continuous(
    limits = c(-1,15),
    expand = c(0, 0),
    breaks = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
  ) +
  theme(
    # Make the background white and remove extra grid lines
    panel.background = element_rect(fill = "white", color = "white"),
    # panel.grid = element_blank(),
    panel.grid.major.x = element_line(color = 'grey', linetype = 'dashed', linewidth = 0.5),
    legend.position = "bottom",
    legend.direction = 'horizontal'
  )

plt_sub


## Supporting Figure 3-B
pdf(paste0('f2/plots/proximal', '_bound2_marked_soft_mode.pdf'),
    width = 15, height = 15)
plt_sub
dev.off()


# distal kmer binding events
## Supporting Figure 4-A
pdf(paste0('f2/plots/', 'distal', '_bound_soft_mode.pdf'),
    width = 12, height = 8)
plt <- top10 %>%
  filter(region == 'distal') %>%
  ggplot(.,
         aes(x = bound_type,
             y = kmer_family,
             color = median_bound_log2fc,
             size = ratio_kmer_HS_bound_sum)) +
  geom_point() +
  scale_color_gradient2(low = 'blue',
                        mid = 'white',
                        high = 'red',
                        midpoint = 0,
                        limits = c(-2,2),
                        oob = scales::squish) +
  labs(x = 'Cluster',
       y = 'TF Family (kmer)',
       color = 'Median HS/CK log2FC',
       size = 'DE bound under HS/total bound') +
  scale_x_discrete(position = "top") +
  facet_grid(genotype~peak_cluster, scales = 'free_y', switch = 'y') +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 10, face = 'bold'),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12, face = 'bold'),
        strip.text.y.left = element_text(size = 12, face = 'bold', angle = 0),
        strip.placement = 'outside')
grid.draw(plt)
dev.off()


data2 <- top10 %>%
  filter(peak_cluster %in% paste0('C', c(1,2,3,4,7,8,15)))
data2$peak_cluster_rank <- vec_rank(as.numeric(data2$peak_cluster_num),
                                    ties = "dense")
plt <- data2 %>%
  filter(region == 'distal',
         !is.na(kmer_family)) %>%
  ggplot(.,
         aes(x = genotype,
             y = peak_cluster_rank,
             color = median_bound_log2fc,
             size = ratio_kmer_HS_bound_sum)) +
  # Make custom panel grid
  geom_hline(
    aes(yintercept = as.numeric(peak_cluster_rank)),
    color = "lightgrey"
  ) +
  geom_point() +
  scale_color_gradient2(low = 'blue',
                        mid = 'white',
                        high = 'red',
                        midpoint = 0,
                        limits = c(-1,1),
                        oob = scales::squish) +
  labs(x = 'Genotype',
       y = 'Peak cluster',
       color = 'Median HS/CK log2FC',
       size = 'DE bound under HS/total bound') +
  scale_x_discrete(position = "top") +
  facet_wrap(.~kmer_family, ncol = 6) +
  coord_polar() +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 10, face = 'bold'),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12, face = 'bold'),
        strip.text.y.left = element_text(size = 12, face = 'bold', angle = 0),
        strip.placement = 'outside')

plt <- plt +
  # Annotate custom scale inside plot
  annotate(
    x = Inf,
    y = 1,
    label = "1",
    geom = "text",
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf,
    y = 2,
    label = "2",
    geom = "text",
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf,
    y = 3,
    label = "3",
    geom = "text",
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf,
    y = 4,
    label = "4",
    geom = "text",
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf,
    y = 5,
    label = "7",
    geom = "text",
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf,
    y = 6,
    label = "8",
    geom = "text",
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf,
    y = 7,
    label = "15",
    geom = "text",
    color = "gray12",
    size = 2
  ) +
  # Scale y axis so bars don't start in the center
  scale_y_continuous(
    limits = c(-1,10),
    expand = c(0, 0),
    breaks = c(1,2,3,4,5,6,7)
  ) +
  theme(
    # Make the background white and remove extra grid lines
    panel.background = element_rect(fill = "white", color = "white"),
    # panel.grid = element_blank(),
    panel.grid.major.x = element_line(color = 'grey', linetype = 'dashed', linewidth = 0.5),
    legend.position = "bottom",
    legend.direction = 'horizontal'
  )

plt


# manipulating kmer families shared among genotypes
# kmer families specific to tak1, hsfa, hsfb
# shared
sh <- c('MYB', 'AP2EREBP', 'TCP', 'C2C2gata', 'bZIP', 'MYBrelated', 'WRKY', 'C2H2', 'FAR1', 'NAC', 'Trihelix')
# specific to a genotype
t <- c('BES1', 'G2like', 'SBP', 'RWPRK')
a <- c('bHLH', 'C2C2dof', 'C2C2YABBY')
b <- c('ABI3VP1', 'NLP', 'BBRBPC')
d <- c('BZR', 'LOBAS2')

data2$sp <- ''
data2[which(data2$kmer_family %in% sh),]$sp <- 'Shared'
data2[which(data2$kmer_family %in% t),]$sp <- 'Tak1'
data2[which(data2$kmer_family %in% a),]$sp <- 'hsfa'
data2[which(data2$kmer_family %in% b),]$sp <- 'hsfb'
data2[which(data2$kmer_family %in% d),]$sp <- 'dko'
data2$sp <- factor(data2$sp, levels = c('Shared', 'Tak1', 'hsfa', 'hsfb', 'dko'))

plt <- data2 %>%
  filter(region == 'distal',
         !is.na(kmer_family)) %>%
  ggplot(.,
         aes(x = genotype,
             y = peak_cluster_rank,
             color = median_bound_log2fc,
             size = ratio_kmer_HS_bound_sum)) +
  # Make custom panel grid
  geom_hline(
    aes(yintercept = as.numeric(peak_cluster_rank)), 
    color = "lightgrey"
  ) +
  geom_point() +
  scale_color_gradient2(low = 'blue',
                        mid = 'white',
                        high = 'red',
                        midpoint = 0,
                        limits = c(-1,1),
                        oob = scales::squish) +
  labs(x = 'Genotype',
       y = 'Peak cluster',
       color = 'Median HS/CK log2FC',
       size = 'DE bound under HS/total bound') +
  scale_x_discrete(position = "top") +
  facet_wrap(.~sp + kmer_family, ncol = 6) +
  coord_polar() +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 10, face = 'bold'),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12, face = 'bold'),
        strip.text.y.left = element_text(size = 12, face = 'bold', angle = 0),
        strip.placement = 'outside')

plt <- plt +
  # Annotate custom scale inside plot
  annotate(
    x = Inf, 
    y = 1, 
    label = "1", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 2, 
    label = "2", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 3, 
    label = "3", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 4, 
    label = "4", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 5, 
    label = "7", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 6, 
    label = "8", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  annotate(
    x = Inf, 
    y = 7, 
    label = "15", 
    geom = "text", 
    color = "gray12",
    size = 2
  ) +
  # Scale y axis so bars don't start in the center
  scale_y_continuous(
    limits = c(-1,10),
    expand = c(0, 0),
    breaks = c(1,2,3,4,5,6,7)
  ) +
  theme(
    # Make the background white and remove extra grid lines
    panel.background = element_rect(fill = "white", color = "white"),
    # panel.grid = element_blank(),
    panel.grid.major.x = element_line(color = 'grey', linetype = 'dashed', linewidth = 0.5),
    legend.position = "bottom",
    legend.direction = 'horizontal'
  )

pdf(paste0('f2/plots/', 'distal', '_bound2_marked_soft_mode.pdf'),
    width = 15, height = 15)
plt
dev.off()
###########################


###### CO-OCCURRENCE KMERS ######
# We only discussed differential binding between HS and Cntr (CK) under the same genotype.
# Therefore, the input1 can be 'Tak1_HS_vs_Tak1_CK', 'hsfa_HS_vs_hsfa_CK', 'hsfb_HS_vs_hsfb_CK', 'dko_HS_vs_dko_CK'.
# And the input2 can be promoter (proximal) or distal.
# differential binding events in 95% or 5% quantile
input_list_1 <- c('Tak1_HS_vs_Tak1_CK', 'hsfa_HS_vs_hsfa_CK', 'hsfb_HS_vs_hsfb_CK', 'dko_HS_vs_dko_CK')
input_list_2 <- c('promoter', 'distal')
input_list_3 <- c('up', 'down')

file_path <- 'f2/tfcomb/bed'

full2 <- data.frame()

for (i in input_list_1) {
  for (z in input_list_2) {
    for (x in input_list_3) {
      
      # read file: significant rules
      print(paste0('Processing: ', i, '_', z, '_', x))
      df <- read.delim(paste(file_path, paste0(i, '_', z, '_', x), 'selected_sig',
                             paste0(i, '_', z, '_', x, '.txt'), sep = '/'))
      df2 <- read.delim(paste(file_path, paste0(i, '_', z, '_', x), 'total',
                              paste0(i, '_', z, '_', x, '.txt'), sep = '/'))
      
      # plot 2 (it wont be affected by how many selected rules we have in this genotype.)
      cos_line <- min(df$cosine)
      zsc_line <- min(df$zscore)
      
      # length of rules
      len_rules <- nrow(df)
      if(len_rules != 0){
        # find the TF family
        prefix <- data.frame(str_locate(df$TF1, '_'))
        prefix$start <- prefix$start + 1
        prefix2 <- str_locate_all(df$TF1, '_')
        for (g in 1:len_rules) {
          sub <- data.frame(prefix2[[g]])
          
          if(nrow(sub >= 2)){
            end <- sub[2,]$end
            
            prefix[g,]$end <- end-1
          }
        }
        
        df$TF_family <- str_sub(df$TF1,
                                start = prefix$start,
                                end = prefix$end)
        
        prefix <- data.frame(str_locate(df$TF2, '_'))
        prefix$start <- prefix$start + 1
        prefix2 <- str_locate_all(df$TF2, '_')
        for (h in 1:len_rules) {
          sub <- data.frame(prefix2[[h]])
          
          if(nrow(sub >= 2)){
            end <- sub[2,]$end
            
            prefix[h,]$end <- end-1
          }
        }
        
        df$TF_family2 <- str_sub(df$TF2,
                                 start = prefix$start,
                                 end = prefix$end)
        # create new name
        df$rules <- paste(pmin(df$TF_family, df$TF_family2),
                          pmax(df$TF_family, df$TF_family2), sep = "-")
        
        # motif class
        motif <- unique(c(df$TF_family, df$TF_family2))
        
        
        # extract cosine and zscore with each rules
        data <- df %>% 
          group_by(rules) %>% 
          summarise(sum_TF1_TF2_count = sum(TF1_TF2_count),
                    median_cosine = median(cosine, na.rm = T),
                    median_zsocre = median(zscore, na.rm = T))
        
        
        # add genotype, region, and expected
        data <- data %>% 
          mutate(genotype = i,
                 region = z,
                 expected = x)
        
        # append data
        full2 <- rbind(full2, data)
        
      }
    }
  }
}

# load('f2/soft_mode_all.RData')

full2_2 <- data.frame()
for (i in input_list_1) {
  for (z in input_list_2) {
    
    # kmer_family which is important
    rus <- top10 %>% 
      filter(genotype == i,
             region == z) %>% 
      pull(kmer_family) %>% 
      unique()
    
    subs <- data.frame()
    for (ru in rus) {
      sub <- full2 %>% 
        filter(genotype == i,
               expected == 'up') %>% 
        filter(., grepl(ru, rules))
      
      # append
      subs <- rbind(subs, sub)
    }
    subs <- distinct(subs)
    
    # append again
    full2_2 <- rbind(full2_2, subs)
  }
}
full2_2 <- distinct(full2_2)

rus2 <- full2_2 %>%
  filter(region == 'promoter') %>% 
  group_by(rules) %>% 
  summarise(count = sum(sum_TF1_TF2_count)) %>% 
  arrange(desc(count)) %>% 
  slice_head(n = 20) %>% 
  pull(rules)


full2_2$genotype <- factor(full2_2$genotype, levels = input_list_1,
                           labels = c('Tak1', 'hsfa', 'hsfb', 'dko'))
full2_2 <- full2_2 %>% 
  filter(rules %in% rus2)
full2_2$rules <- factor(full2_2$rules, levels = rus2)


full2_2_filter <- full2_2 %>% 
  filter(rules %in% (full2_2 %>% 
                       group_by(rules) %>% 
                       summarise(fre = n()) %>% 
                       filter(fre != 1) %>% 
                       pull(rules))) %>% 
  filter(region == 'promoter')

full2_2 <- full2_2 %>% 
  filter(region == 'promoter')


pdf(paste0('f2/rules_distribution_all/plot/',
           paste0('promoter_top_rules', '.pdf')),
    width = 8, height = 5)
plt <- full2_2 %>% 
  ggplot(data = .) +
  geom_boxplot(data = full2_2_filter,
               aes(x = rules,
                   y = median_cosine,
                   group = rules),
               width = 0.5) +
  geom_point(data = full2_2,
             aes(x = rules,
                 y = median_cosine,
                 fill = genotype,
                 size = sum_TF1_TF2_count),
             shape = 21,
             alpha = 0.5) +
  theme_classic() +
  labs(x = 'Rules',
       y = 'Median cosine',
       fill = 'Genotype',
       size = 'Co-occurence counts') +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text.x = element_text(size = 8, face = 'bold', angle = 40, hjust = 1),
        axis.text.y = element_text(size = 8, face = 'bold'))
plt
dev.off()


ct <- data.frame()
for (i in input_list_1) {
  for (z in input_list_2) {
    for (x in input_list_3) {
      
      # read file: significant rules
      print(paste0('Processing: ', i, '_', z, '_', x))
      df <- read.delim(paste(file_path, paste0(i, '_', z, '_', x), 'selected_sig',
                             paste0(i, '_', z, '_', x, '.txt'), sep = '/'))
      df2 <- read.delim(paste(file_path, paste0(i, '_', z, '_', x), 'total',
                              paste0(i, '_', z, '_', x, '.txt'), sep = '/'))
      
      
      total <- nrow(df2)
      selected <- nrow(df)
      
      sub <- data.frame(genotype = i,
                        total = total,
                        selected = selected,
                        region = z,
                        expected = x)
      # append
      ct <- rbind(ct, sub)
    }
  }
}


ct$genotype <- factor(ct$genotype, levels = input_list_1,
                      labels = c('Tak1', 'hsfa', 'hsfb', 'dko'))
ct <- ct %>%
  filter(region == 'promoter',
         expected == 'up') %>%
  dplyr::select(-region, -expected) %>% 
  gather(key = 'all_or_selected', value = 'counts', -genotype)
ct$all_or_selected <- factor(ct$all_or_selected, levels = c('total', 'selected'))

pdf(paste0('f2/rules_distribution_all/plot/',
           paste0('promoter_rules', '.pdf')),
    width = 6, height = 4)
ct %>% 
  ggplot(.,
         aes(x = genotype,
             y = counts,
             fill = all_or_selected,
             label = counts)) +
  geom_bar(stat = 'identity') +
  geom_text(data = ct,
            aes(x = genotype,
                y = counts + 100,
                label = counts),
            color = ifelse(ct$all_or_selected == 'selected', 'red', 'black')) +
  scale_fill_manual(values = c('grey', 'red')) +
  coord_flip() +
  theme_classic() +
  labs(x = 'Genotype',
       y = 'Counts',
       fill = 'Rules') +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 8, face = 'bold'))
dev.off()
###########################


###### KMER DISTRIBUTION ######
promoter_motif <- read.delim('f2/sim/mp_acr_promoter/promoter_pcc.txt')
distal_motif <- read.delim('f2/sim/mp_acr_distal/distal_pcc.txt')
input_list_1 <- c('Tak1_HS_vs_Tak1_CK', 'hsfa_HS_vs_hsfa_CK', 'hsfb_HS_vs_hsfb_CK', 'dko_HS_vs_dko_CK')
input_list_2 <- c('promoter', 'distal')
input_list_3 <- c('soft')
p_path <- 'f2/distribution_all'
p_path2 <- 'f2/distribution_all/plot'
Txdb_gtf <- makeTxDbFromGFF('f1/MpTak_v6.1r1.gff')


for (z in input_list_1) {
  for (j in input_list_2) {
    for (h in input_list_3) {
      input1 <- z
      input2 <- j
      input3 <- h
      
      # load top10 file in strict or soft mode
      # you will get top10 and full dataframe from kmer bindetect/2
      load(paste0('f2/', input3, '_mode_all.RData'))
      
      
      if(input1 %in% input_list_1 & input2 %in% input_list_2){
        print(paste0('Processing: ', input1, ' under ', input2, ' (', input3, ' mode)'))
        path <- paste('f2/bind/', input2, input1, sep = '/')
        bound_data <- read.delim(paste(path, 'bindetect_results.txt', sep = '/'))
        
        # remove no bound under HS & CK
        bound_data2 <- bound_data[which(bound_data$HS_bound != 0 | bound_data$CK_bound != 0),]
        
        # remove prefix of TFBS_name
        bound_data2$kmer_motif_name <- str_remove(bound_data2$kmer_motif_name, '_')
        
        # check all motifs are included in corresponding kmer pool
        if (input2 == 'promoter'){
          kmer_pool <- promoter_motif
          if(!all(bound_data2$kmer_motif_name %in% kmer_pool$kmers)){
            stop('Not all motifs in bound data exist in user-provided kmer pool.\n')
          }
        } else{
          kmer_pool <- distal_motif
          all(bound_data2$kmer_motif_name %in% kmer_pool$kmers)
          if(!all(bound_data2$kmer_motif_name %in% kmer_pool$kmers)){
            stop('Not all motifs in bound data exist in user-provided kmer pool.\n')
          }
        }
        
        
        # add family name
        meme <- read_meme('f2/ArabidopsisDAPv1.meme')
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
        
        bound_data2$kmer_family <- ''
        bound_data3 <- data.frame()
        for (i in 1:length(family)) {
          sub <- bound_data2 %>% 
            filter(., grepl(paste0('_',family[i],'_'), kmer_motif_name)) %>% 
            mutate(kmer_family = family[i])
          
          # append
          bound_data3 <- rbind(bound_data3, sub)
        }
        
        
        # annotation TFBS distance
        # extract Chr, Start, End, Peak name
        data <- bound_data3[,c('TFBS_chr', 'TFBS_start', 'TFBS_end', 'kmer_motif_name')]
        names(data) <- c('Chr', 'Start', 'End', 'Motif_name')
        head(data)
        
        # make GRanges data
        gr <- makeGRangesFromDataFrame(data, keep.extra.columns=T)
        
        # annotate
        peakAnno <- annotatePeak(gr,
                                 tssRegion=c(-1500, 1500),
                                 genomicAnnotationPriority = c("Exon", "Intron", "5UTR", "3UTR", "Promoter", "Downstream", "Intergenic"),
                                 TxDb=Txdb_gtf,
                                 level ='gene')
        
        # export peaks annotation
        peakAnno
        peakAnno@anno
        peakAnno_tb <- as_tibble(peakAnno@anno) %>% 
          mutate(`distanceToTSS (log10)` = log10((abs(distanceToTSS)/1000) + 1))
        
        peakAnno_tb[which(peakAnno_tb$distanceToTSS < 0),]$`distanceToTSS (log10)` <- -peakAnno_tb[which(peakAnno_tb$distanceToTSS < 0),]$`distanceToTSS (log10)`
        
        
        bound_data3_dis <- peakAnno_tb %>% 
          dplyr::select(seqnames, start, end, Motif_name, geneId, distanceToTSS, `distanceToTSS (log10)`) %>% 
          dplyr::rename(TFBS_chr = seqnames,
                        TFBS_start = start,
                        TFBS_end = end,
                        kmer_motif_name = Motif_name) %>% 
          merge(., bound_data3)
        
        
        geno <- str_remove(input1, '_.*')
        write.table(bound_data3_dis, paste0(p_path, '/', input3, '/', input2, '/',
                                            input3, '_', input2, '_', geno, '_peak_distribution.txt'),
                    row.names = F, quote = F, sep = '\t')
        
        
        # plot
        plot_data <- rbind(bound_data3_dis %>% 
                             dplyr::select(distanceToTSS,
                                           `distanceToTSS (log10)`,
                                           peak_cluster,
                                           kmer_family))
        
        
        family <- top10 %>% 
          filter(region == input2,
                 genotype == input1)
        kmer_family <- family %>% 
          dplyr::select(kmer_family) %>% 
          distinct() %>% 
          pull()
        peak_cluster <- family %>% 
          dplyr::select(peak_cluster) %>% 
          distinct()
        peak_c <- as.character(peak_cluster$peak_cluster)
        
        
        
        
        dis.sum <- data.frame()
        for (m in kmer_family) {
          df <- plot_data %>% filter(kmer_family %in% kmer_family &
                                       peak_cluster %in% peak_c)
          
          # add 1000bp to make -1kb to 0.5kb in range(0,1500)
          df$Preferential_Position <- as.numeric(df$distanceToTSS) + 1000
          df <- df %>%
            group_by(Preferential_Position, kmer_family) %>% 
            summarise(count = n())
          names(df) <- c('location', 'kmer_family', 'tar')
          
          
          window.tar <- data.frame(location = 1:1500)
          window.tar <- left_join(window.tar, df %>% 
                                    filter(kmer_family == m) %>% 
                                    dplyr::select(-kmer_family), by = 'location')
          window.tar[is.na(window.tar$tar),]$tar <- 0
          window.tar$tar <- window.tar$tar + 1 
          
          
          
          df2 <- plot_data %>% filter(kmer_family %in% kmer_family)
          # add 1000bp to make -1kb to 0.5kb in range(0,1500)
          df2$Preferential_Position <- as.numeric(df2$distanceToTSS) + 1000
          df2 <- df2 %>%
            group_by(Preferential_Position, kmer_family) %>% 
            summarise(count = n())
          names(df2) <- c('location', 'kmer_family', 'beg')
          
          
          window.beg <- data.frame(location = 1:1500)
          window.beg <- full_join(window.beg, df2 %>% 
                                    filter(kmer_family == m) %>% 
                                    dplyr::select(-kmer_family), by = 'location')
          window.beg[is.na(window.beg$beg),]$beg <- 0
          window.beg$beg <- window.beg$beg + 1 
          
          window <- merge(window.tar, window.beg, by = 'location')
          
          
          ## culculation
          # all
          # bin:100; sliding window:25; median
          tar.win <- rollapply(window$tar, width = 100, mean, by = 25, partial = FALSE)
          beg.win <- rollapply(window$beg, width = 100, mean, by = 25, partial = FALSE)
          
          
          
          dis <- data.frame(bin = seq(from = 100, to = 1500, by = 25))
          dis$tar <- tar.win
          dis$beg <- beg.win
          
          
          dis <- dis %>% 
            mutate(zscore_tar = (tar - mean(tar))/sd(tar),
                   zscore_beg = (beg - mean(beg))/sd(beg))
          dis <- dis %>% 
            gather(., key = 'group', value = 'zscore', -bin, -tar, -beg) %>% 
            mutate(kmer_family = m)
          
          # append dis
          dis.sum <- rbind(dis.sum, dis)
        }
        
        write.table(dis.sum, paste0(p_path, '/', input3, '/', input2, '/',
                                    input3, '_', input2, '_', geno, '_peak_zscore.txt'),
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
          facet_grid(.~factor(group, levels = c('zscore_tar', 'zscore_beg'), labels = c('Target', 'background'))) +
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
        
        
        ## Figure 3-A
        pdf(paste0(p_path2, '/',
                   input3, '_', input2, '_', geno, '_peak_zscore.pdf'),
            width = 7, height = 2)
        print(plot)
        dev.off()
      }
    }
  }
}
######################################


###### RNA-SEQ DEG: MARCHANTIA ######
p_gene <- 'f3'
dir.create(paste(p_gene, 'gene', sep = '/'))

# TP
tp <- read.delim('f3/all_genes_hclust.txt')

# TN
range <- rbind(data.frame(dir = rep(c('f3/DEG'),
                                    each = 4),
                          group1 = c('Tak1_heat', 'hsfa_heat',
                                     'hsfb_heat', 'dko_heat'),
                          group2 = rep(c('Tak1_NHS'), each = 4),
                          expr = c('n.txt'),
                          ex = c('N'),
                          label = paste(rep(c('HS'), each = 4),
                                        c('T', 'a', 'b', 'd'))),
               data.frame(dir = rep(c('f3/DEG'), 6),
                          group1 = rep(c('hsfa_heat', 'hsfb_heat', 'dko_heat'), each = 2),
                          group2 = rep(c('hsfa_NHS', 'hsfb_NHS', 'dko_NHS'), each = 2),
                          expr = c('n.txt'),
                          ex = c('N'),
                          label = paste(rep(c('HS'), each = 6), rep(c('a', 'b', 'd'), each = 2)))) %>% 
  distinct()

colname <- c('gene', 'group')
df <- data.frame(matrix(nrow = 0, ncol = length(colname)))
names(df) <- colname
for (i in 1:nrow(range)) {
  df.sub <- read.delim(paste(range$dir[i],
                             paste(range$group1[i], range$group2[i], range$expr[i], sep = '_'),
                             sep = '/')) %>%
    dplyr::rename(gene = x) %>% 
    mutate(group = paste(range$ex[i], range$label[i], sep = '_'))
  df <- rbind(df, df.sub)
}

tn <- df %>% 
  dplyr::select(gene) %>% 
  distinct() %>% 
  filter(!(gene %in% tp$gene))

write.table(df %>% 
              dplyr::select(gene) %>% 
              distinct() %>% 
              filter(!(gene %in% tp$gene)), paste(p_gene, 'gene', 'TN.txt', sep = '/'),
            row.names = F, quote = F, col.names = F)



# length of negative genes
len_TN <- length(unique(tn$gene))
# make balanced negative files
range <- paste0('C', 1:8)
for (r in range) {
  write.table(tp %>% 
                filter(cluster == r) %>% 
                dplyr::select(gene) %>% 
                distinct() %>% 
                pull(), paste(p_gene, 'gene', paste0(r, '.txt'), sep = '/'), row.names = F, col.names = F, quote = F, sep = '\t')
}


df <- read.table(paste(p_gene, 'gene/TN.txt', sep = '/'))
setwd(paste(p_gene, 'gene', sep = '/'))
for (z in range) {
  cluster <- z
  neg <- df
  n.neg <- nrow(neg)
  pos <- read.table(paste0(cluster, '.txt'))
  n.pos <- nrow(pos)
  
  if (!dir.exists(cluster)){
    dir.create(cluster)
  }
  
  for (i in 1:10){
    number <- i
    set.seed(number)
    random.neg <- sample(seq_len(n.neg), size = n.pos)
    data <- neg[random.neg,]
    
    # output result
    write.table(data, 
                paste(cluster, paste0(paste('neg',
                                            cluster,
                                            number,
                                            sep = '_'), '.txt'), sep = '/'),
                row.names = F, col.names = F, quote = F, sep = '\t')
  }
}
# setwd("~/folder_to_all_the_data")
###############################


###### RNA-SEQ DEG + ATAC: MARCHANTIA ######
p_gene <- 'f3'
dir.create(paste(p_gene, 'pg2gene', sep = '/'))

# TP
range <- list.dirs('f3/cluster_ovl',
                   full.names = F,
                   recursive = F)
tp <- data.frame()
for (cl in range) {
  files <- list.files(paste('f3/cluster_ovl', cl, sep = '/'),
                      full.names = F)
  tp.sub <- data.frame()
  for (file in files) {
    sub <- read.delim(paste('f3/cluster_ovl', cl, file, sep = '/')) %>% 
      dplyr::select(gene)
    
    tp.sub <- rbind(tp.sub, sub) %>% 
      distinct()
  }
  
  tp.sub <- tp.sub %>% 
    mutate(cluster = cl)
  tp <- rbind(tp, tp.sub) %>% 
    distinct()
}

# length of tp
length(unique(tp$gene)) # 933

for (cl in range) {
  # make directory for TP
  dir.create(paste(p_gene, 'pg2gene', cl, sep = '/'))
  
  # export tp gene file in each cluster
  tp %>% 
    filter(cluster == cl) %>% 
    dplyr::select(gene) %>% 
    write.table(., paste(p_gene, 'pg2gene', paste0(cl, '.txt'), sep = '/'),
                row.names = F, col.names = F, quote = F, sep = '\t')
}


# TN
# merge of tn from RNA-seq (gene) and ATAC-seq(peak2gene)
tn1 <- read.delim('f1/counts/split/peak2gene_TN.txt')
tn2 <- read.delim('f3/gene/TN.txt', header = F) %>% 
  dplyr::rename(geneId = V1)
head(tn1)
head(tn2)

tn <- inner_join(tn1, tn2) %>% 
  distinct() # 2323


# make balanced negative files
setwd(p_gene)
df <- tn
setwd(paste(p_gene, 'pg2gene', sep = '/'))
for (z in range) {
  cluster <- z
  neg <- df
  n.neg <- nrow(neg)
  pos <- read.table(paste0(cluster, '.txt'), 
                    header = F)
  n.pos <- nrow(pos)
  
  if (!dir.exists(cluster)){
    dir.create(cluster)
  }
  
  for (i in 1:10){
    number <- i
    set.seed(number)
    random.neg <- sample(seq_len(n.neg), size = n.pos)
    data <- neg[random.neg,]
    
    # output result
    write.table(data, 
                paste(cluster, paste0(paste('neg',
                                            cluster,
                                            number,
                                            sep = '_'), '.txt'), sep = '/'),
                row.names = F, col.names = F, quote = F, sep = '\t')
  }
}
###############################


###### RNA-SEQ DEG: SUMMARY FEATURES COUNTS ######
# features
p_feat <- 'f3/feature/mp_rna'

feat <- data.frame()
files <- list.files(p_feat,
                    full.names = F)
for (file in files) {
  file_name <- str_remove(file, '_kmer.txt')
  df <- read.table(paste(p_feat, file, sep = '/'), header = T)
  
  if(colnames(df)[1] == 'file'){
    kmer_n <- ncol(df)-1
  }
  
  # append numbers of feature
  feat.sub <- data.frame(type = 'mp_acr_gene',
                         cluster = file_name,
                         feature = kmer_n)
  
  # append numbers of feature
  feat <- rbind(feat, feat.sub)
}

# visualize feature counts
p_plot <- 'f3/plots/'
feat$cluster <- factor(feat$cluster, levels = paste0('C', 15:1))


## Figure 7-B, AUCROC of DEG only
pdf(paste0(p_plot,'feature_rna.pdf'),
    width = 5, height = 4)
feat %>% 
  ggplot(.,
         aes(x = cluster,
             y = feature,
             fill = factor(type, levels = 'type', label = 'DEGs'))) +
  geom_bar(stat = 'identity',
           position = 'dodge',
           width = 0.5) +
  geom_text(aes(x = cluster,
                y = as.numeric(feature) + 100,
                label = feature),
            color = 'red') +
  theme_classic() +
  coord_flip() +
  labs(x = 'Cluster',
       y = 'Feature',
       fill = 'Type') +
  facet_grid(cluster~., scales = 'free_y') +
  theme(strip.background = element_blank(),
        strip.text = element_blank())
dev.off()
###############################


###### RNA-SEQ DEG + ATAC: SUMMARY FEATURES COUNTS ######
# features
p_feat <- 'f3/feature/mp_rna_ovl'

feat <- data.frame()
files <- list.files(p_feat,
                    full.names = F)
for (file in files) {
  file_name <- str_remove(file, '_kmer.txt')
  df <- read.table(paste(p_feat, file, sep = '/'), header = T)
  
  if(colnames(df)[1] == 'file'){
    kmer_n <- ncol(df)-1
  }
  
  # append numbers of feature
  feat.sub <- data.frame(type = 'mp_acr_gene',
                         cluster = file_name,
                         feature = kmer_n)
  
  # append numbers of feature
  feat <- rbind(feat, feat.sub)
}

# visualize feature counts
p_plot <- 'f3/plots/'
feat$cluster <- factor(feat$cluster, levels = paste0('C', 15:1))


## Figure 7-B, AUCROC of DEG + ATAC
pdf(paste0(p_plot,'feature_rna_atac.pdf'),
    width = 5, height = 4)
feat %>% 
  ggplot(.,
         aes(x = cluster,
             y = feature,
             fill = factor(type, levels = 'type', label = 'DEGs'))) +
  geom_bar(stat = 'identity',
           position = 'dodge',
           width = 0.5) +
  geom_text(aes(x = cluster,
                y = as.numeric(feature) + 100,
                label = feature),
            color = 'red') +
  theme_classic() +
  coord_flip() +
  labs(x = 'Cluster',
       y = 'Feature',
       fill = 'Type') +
  facet_grid(cluster~., scales = 'free_y') +
  theme(strip.background = element_blank(),
        strip.text = element_blank())
dev.off()
###############################


###### RNA-SEQ DEG: SUMMARY TRAIN RESULTS ######
results <- read.delim("f3/RESULTS_RNA.txt")
head(results)

# remove prefix and tail in column ID
results$ID <- str_remove(results$ID, '/RAID1/working/R425/lavakau/ML_master/')
results$ID <- str_remove(results$ID, '\\.fa.pcre_df_p0.01_mod.txt_SVM|\\.fa.pcre_df_p0.01_mod.txt_RF|\\.fa.pcre_df_p0.01_mod.txt_LogReg')
head(results)

results$cluster <- factor(results$cluster, levels = paste0('C', 1:15))

table(results$type)


auc_no_lasso <- results %>% 
  group_by(cluster, Alg, type) %>% 
  summarise(mean_AUCROC_test = mean(AUCROC_test),
            sd_AUCROC_test = sd(AUCROC_test)) %>% 
  mutate(lab = paste0(round(mean_AUCROC_test, 2), '±', round(sd_AUCROC_test, 3)))


plt <- ggplot(auc_no_lasso %>% 
                filter(Alg == "RF"),
              aes(x = cluster,
                  y = factor(type, levels = "mp_rna", labels = "DEG only"),
                  fill = mean_AUCROC_test)) +
  geom_tile() +
  labs(title = expression(bold('Machine Learning Test Results')),
       x = 'Cluster',
       y = 'type',
       fill = 'AUCROC') + # change y axis F1/AUC
  coord_flip() +
  theme_classic() +
  scale_fill_gradientn(
    colours = c("white", "red"),         # Colors for different ranges
    limits = c(0.5, 1)                             # Set data range
  ) +
  theme(axis.text = element_text(size = 10, face = 'bold'),
        axis.title = element_text(size = 12, face = 'bold'))


## Figure 7-B, AUCROC of DEG only
pdf('f3/plots/AUCROC.pdf',
    width = 3, height = 5)
plt
dev.off()
###############################


###### RNA-SEQ DEG + ATAC: SUMMARY TRAIN RESULTS ######
# BalancedRuns 10 CVfold 10
p_train <- 'f3'
results <- read.delim(paste(p_train, 'RESULTS_1.txt', sep = '/'))
head(results)

# remove prefix and tail in column ID
results$ID <- str_remove(results$ID, '/RAID1/working/R425/lavakau/ML_master/')
results$ID <- str_remove(results$ID, '\\.fa.pcre_df_p0.01_mod.txt_SVM|\\.fa.pcre_df_p0.01_mod.txt_RF|\\.fa.pcre_df_p0.01_mod.txt_LogReg')
head(results)

# separate ID columns as train type, cluster, and numbers of negative files
results <- results %>% 
  separate(., col = 'ID', into = c('type', 'cluster', 'negative_data'), sep = '/')
results$cluster <- factor(results$cluster, levels = paste0('C', 1:15))

# clean again
results$negative_data <- str_remove(results$negative_data, '_lasso|_no_lasso')


result_lasso <- data.frame()
result_no_lasso <- data.frame()
for (i in 1:8) {
  clu <- paste0('C', i)
  sub <- results %>% 
    filter(cluster == clu)
  position <- which(sub$negative_data  == paste('neg', clu, 2, sep = '_'))
  
  if(length(position) > 4){
    sub_lasso <- sub[1:(position[4]-1),]
    sub_no_lasso <- sub[position[4]:nrow(sub),]
    
    result_lasso <- rbind(result_lasso, sub_lasso)
    result_no_lasso <- rbind(result_no_lasso, sub_no_lasso)
  }
  
}


write.table(results, 'f3/table/results.txt',
            row.names = F, quote = F, sep = '\t')
write.table(result_lasso, 'f3/table/result_lasso.txt',
            row.names = F, quote = F, sep = '\t')
write.table(result_no_lasso, 'f3/table/result_no_lasso.txt',
            row.names = F, quote = F, sep = '\t')

auc_no_lasso <- result_no_lasso %>% 
  dplyr::select(cluster, negative_data, Alg, AUCROC_test, type) %>% 
  group_by(cluster, Alg, type) %>% 
  summarise(mean_AUCROC_test = mean(AUCROC_test),
            sd_AUCROC_test = sd(AUCROC_test)) %>% 
  mutate(lab = paste0(round(mean_AUCROC_test, 2), '±', round(sd_AUCROC_test, 3)))


plt <- ggplot(auc_no_lasso %>% 
                filter(Alg == "RF",
                       type == "mp_rna_ovl"),
              aes(x = cluster,
                  y = factor(type, levels = "mp_rna_ovl", labels = "DEG + ATAC"),
                  fill = mean_AUCROC_test)) +
  geom_tile() +
  labs(title = expression(bold('Machine Learning Test Results')),
       x = 'Cluster',
       y = 'type',
       fill = 'AUCROC') + # change y axis F1/AUC
  coord_flip() +
  theme_classic() +
  scale_fill_gradientn(
    colours = c("white", "red"),         # Colors for different ranges
    limits = c(0.5, 1)                             # Set data range
  ) +
  theme(axis.text = element_text(size = 10, face = 'bold'),
        axis.title = element_text(size = 12, face = 'bold'))
plt


## Figure 7-B, AUCROC of DEG + ATAC
pdf('f3/plots/AUCROC2.pdf',
    width = 3, height = 5)
plt
dev.off()
###############################


###### LOAD ATAC-SEQ COUNT MATRIX: TOMATO ######
p <- 'f2/tomato'
df <- read.table(paste(p, 'counts.txt', sep = '/'), header=T)

# # generate mata files
meta <- data.frame(names = names(df[,-1])) %>%
  separate(., col = names, into = c('experiment', 'cultivar', 'treatment', 'replicate')) %>%
  dplyr::select(-experiment)
meta$treatment <- factor(meta$treatment, levels = c('0h', '1h', '6h'))

meta <- meta[order(meta$treatment),] %>%
  mutate(sample = paste0('s', 1:nrow(.)))
head(meta)
write.table(meta, paste(p, 'meta.txt', sep = '/'), row.names = F, quote = F, sep = '\t')

meta <- read.table(paste(p, 'meta.txt', sep = '/'), header=T)
df <- df %>% 
  column_to_rownames(var = 'Peaks')


names(df) <- meta$sample
meta <- meta %>% 
  column_to_rownames(var = 'sample')
#######################


####### REMOVE MISSING PEAKS: TOMATO #######
# remove all the rows where not a single sample has more than 50 reads
df <- df[apply(df, 1, max) > 50,]
# peaks from 17800 to 17794
#######################


###### DATA NORMALIZATION: TOMATO ######
dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = meta,
                              design = ~ treatment)
ddsMF <- dds
ddsMF$treatment <- relevel(ddsMF$treatment, ref = '0h')
ddsMF <- DESeq(ddsMF)

vsdMF <- vst(ddsMF, blind = FALSE)
head(assay(vsdMF), 5)
normalized.count <- data.frame(assay(vsdMF)) %>% 
  rownames_to_column(., var = 'peak')
write.table(normalized.count, paste(p, 'normalized_counts.txt', sep = '/'),
            row.names = F, quote = F, sep = '\t')
#######################


###### PCA: TOMATO ######
pcaData <- plotPCA(vsdMF, intgroup = c("treatment", "cultivar"), returnData = TRUE)
percenVar <- round(100 * attr(pcaData, "percentVar"))

plt <- ggplot(pcaData, aes(PC1, PC2,
                          shape = cultivar,
                          color = treatment)) + 
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percenVar[1], "% variance")) +
  ylab(paste0("PC2: ", percenVar[2], "% variance")) +
  labs(shape = 'cultivar', color = 'treatment') +
  scale_color_manual(values = c('#bdc9e1', '#67a9cf', '#02818a')) +
  theme_classic() +
  coord_fixed() +
  theme(axis.text = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 12, face = 'bold'),
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        axis.line = element_blank())

pdf(paste(p, 'plot/PCA.pdf', sep = '/'),
    width = 5, height = 5)
plt
dev.off()

save.image(file = paste(p, 'HS_peak.RData', sep = '/'))
#######################


###### ACR WITH PROXIMAL AND DISTAL: TOMATO ######
# make TP bed
# Strand of each peak need to follow its annotated gene
# Peak name, Chr, Start, End, Strand
# setwd("~/folder_to_all_the_data")

p <- 'f2/tomato'
p_b_pro <- paste0(p, '/bed_proximal')
p_b_dis <- paste0(p, '/bed_distal')
p_pro <- paste0(p, '/peak_proximal')
p_dis <- paste0(p, '/peak_distal')
# dir.create(p_b_pro)
# dir.create(p_b_dis)
# dir.create(p_pro)
# dir.create(p_dis)

de_peak <- read.table(paste(p, 'DE_peaks_hclust.txt', sep = '/'), header = T)
head(de_peak)

Txdb_gtf <- makeTxDbFromGFF(paste(p, 'Slycopersicum_796_ITAG5.0.gene.gff3', sep = '/'))

for (i in paste0('C', 1:5)) {
  
  # extract Chr, Start, End, Peak name
  data <- de_peak %>% 
    filter(cluster == i) %>% 
    dplyr::select(peak) %>%
    mutate(Peak = peak) %>% 
    separate(., col = 'peak', into = c('Chr', 'Start', 'End'), sep = '\\.') %>%
    mutate(Start = as.numeric(Start),
           End = as.numeric(End)) %>% 
    relocate(Chr, Start, End)
  
  # make GRanges data
  gr <- makeGRangesFromDataFrame(data, keep.extra.columns=T)
  
  # annotate
  peakAnno <- annotatePeak(gr,
                           tssRegion=c(-1500, 1500),
                           genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
                           TxDb=Txdb_gtf,
                           level ='gene')
  
  # export promoter peaks annotation
  peakAnno_tb <- as_tibble(peakAnno@anno) %>%
    filter(annotation == 'Promoter') %>% 
    dplyr::select(Peak, seqnames, start, end, geneStrand) %>% 
    relocate(Peak, seqnames, start, end, geneStrand)
  
  # export promoter bed and peak
  setwd(p_b_pro)
  write.table(peakAnno_tb, paste0(i, '.bed'), row.names = F, col.names = F, quote = F, sep = '\t')
  setwd(p_pro)
  write.table(peakAnno_tb$Peak, paste0(i, '.txt'), row.names = F, col.names = F, quote = F, sep = '\t')
  
  
  # export distal peaks annotation
  peakAnno_tb <- as_tibble(peakAnno@anno) %>%
    filter(annotation != 'Promoter') %>% 
    dplyr::select(Peak, seqnames, start, end, geneStrand) %>% 
    relocate(Peak, seqnames, start, end, geneStrand)
  
  # export distal bed and peak
  setwd(p_b_dis)
  write.table(peakAnno_tb, paste0(i, '.bed'), row.names = F, col.names = F, quote = F, sep = '\t')
  setwd(p_dis)
  write.table(peakAnno_tb$Peak, paste0(i, '.txt'), row.names = F, col.names = F, quote = F, sep = '\t')
}


# make balanced negative files
range <- paste0('C', 1:5)
p_neg <- paste0(p, '/negative')
dir.create(p_neg)
head(p_pro)
head(p_dis)


## make TN.txt
## overlapping negative peaks from 1H and 6H
# setwd("~/folder_to_all_the_data")
neg_files <- list.files(path = paste0(p, '/neg'),
                        pattern = '*_n.txt',
                        full.names = T)
neg_list <- list()
for (file in neg_files) {
  sub_list <- list(read.delim(file) %>% pull())
  neg_list <- c(neg_list, sub_list)
}
neg_intersect <- intersect(neg_list[[1]], neg_list[[2]])
setwd(p_neg)
write.table(data.frame(peaks = neg_intersect), 'TN.txt',
            row.names = F, quote = F, sep = '\t')
## we got 9491 negative peaks


df <- read.table(paste(p_neg, 'TN.txt', sep = '/'), header = T) 

dirs <- c(p_pro, p_dis)
for (dir in dirs) {
  setwd(dir)
  for (z in range) {
    cluster <- z
    neg <- df
    n.neg <- nrow(neg)
    pos <- read.table(paste0(cluster, '.txt'))
    n.pos <- nrow(pos)
    
    if (!dir.exists(cluster)){
      dir.create(cluster)
    }
    
    # balance sets
    for (i in 1:10){
      number <- i
      set.seed(number)
      random.neg <- sample(seq_len(n.neg), size = n.pos)
      data <- neg[random.neg,]
      
      # output result
      write.table(data, 
                  paste(cluster, paste0(paste('neg',
                                              cluster,
                                              number,
                                              sep = '_'), '.txt'), sep = '/'),
                  row.names = F, col.names = F, quote = F, sep = '\t')
    }
  }
}
######################################


###### SUMMARY TRAIN RESULTS: TOMATO ######
p <- 'f2/tomato'
results <- read.delim(paste(p, 'RESULTS.txt', sep = '/'))
head(results)

# remove prefix and tail in column ID
results$ID <- str_remove(results$ID, '/RAID1/working/R425/lavakau/ML_master/')
results$ID <- str_remove(results$ID, '.pcre_df_p0.01_mod.*')
head(results)

# separate ID columns as train type, cluster, and numbers of negative files
results <- results %>% 
  separate(., col = 'ID', into = c('type', 'cluster', 'negative_data'), sep = '/')
results$cluster <- factor(results$cluster, levels = paste0('C', 1:5))
results$type <- factor(results$type, levels = c("sly_acr_pro", "sly_acr_dis"),
                       labels = c("PROXIMAL", "DISTAL"))

f1 <- results %>% 
  dplyr::select(type, cluster, negative_data, Alg, F1_test, Tag) %>% 
  group_by(type, cluster, Alg, Tag) %>% 
  summarise(mean_F1_test = mean(F1_test),
            sd_F1_test = sd(F1_test))

auc <- results %>% 
  dplyr::select(type, cluster, negative_data, Alg, AUCROC_test, Tag) %>% 
  group_by(type, cluster, Alg, Tag) %>% 
  summarise(mean_AUCROC_test = mean(AUCROC_test),
            sd_AUCROC_test = sd(AUCROC_test))



plt <- ggplot(results %>% 
                group_by(type, cluster, Alg) %>% 
                summarise(F1_test_mean = mean(F1_test, na.rm = T),
                          AUCROC_test_mean = mean(AUCROC_test, na.rm = T)) %>% 
                filter(type %in% c("PROMOTER", "DISTAL"),
                       Alg == "RF"),
              aes(x = cluster,
                  y = type,
                  fill = F1_test_mean)) +
  geom_tile() +
  labs(title = expression(bold('Machine Learning Test Results')),
       x = 'Cluster',
       y = 'Type',
       fill = 'F1') + # change y axis F1/AUC
  coord_flip() +
  theme_classic() +
  scale_fill_gradientn(
    colours = c("white", "red"),         # Colors for different ranges
    limits = c(0.5, 1)                             # Set data range
  ) +
  theme(axis.text = element_text(size = 10, face = 'bold'),
        axis.title = element_text(size = 12, face = 'bold'))


plt <- ggplot(results %>% 
                group_by(type, cluster, Alg) %>% 
                summarise(F1_test_mean = mean(F1_test, na.rm = T),
                          AUCROC_test_mean = mean(AUCROC_test, na.rm = T)) %>% 
                filter(type %in% c("PROMOTER", "DISTAL"),
                       Alg == "RF"),
              aes(x = cluster,
                  y = type,
                  fill = AUCROC_test_mean)) +
  geom_tile() +
  labs(title = expression(bold('Machine Learning Test Results')),
       x = 'Cluster',
       y = 'Type',
       fill = 'AUCROC') + # change y axis F1/AUC
  coord_flip() +
  theme_classic() +
  scale_fill_gradientn(
    colours = c("white", "red"),         # Colors for different ranges
    limits = c(0.5, 1)                             # Set data range
  ) +
  theme(axis.text = element_text(size = 10, face = 'bold'),
        axis.title = element_text(size = 12, face = 'bold'))
#############################################


###### SUMMARY FEATURES COUNTS: TOMATO ######
# features
p_feat <- 'f2/tomato/feature'
dirs <- list.dirs(p_feat,
                  full.names = F,
                  recursive = F)

setwd(p_feat)

feat <- data.frame()
for (dir in dirs) {
  feat.sub <- data.frame()
  files <- list.files(dir,
                      full.names = F)
  for (file in files) {
    file_name <- str_remove(file, '_kmer.txt')
    df <- read.table(paste(dir, file, sep = '/'), header = T)
    
    if(colnames(df)[1] == 'file'){
      kmer_n <- ncol(df)-1
    }
    
    # append numbers of feature
    sub <- data.frame(type = dir,
                      cluster = file_name,
                      feature = kmer_n)
    feat.sub <- rbind(feat.sub, sub)
  }
  
  # append numbers of feature
  feat <- rbind(feat, feat.sub)
}

# visualize feature counts
p_plot <- 'f2/tomato/plot/'
feat$cluster <- factor(feat$cluster, levels = paste0('C', 5:1))

pdf(paste0(p_plot,'feature.pdf'),
    width = 4, height = 4)
feat %>% 
  ggplot(.,
         aes(x = cluster,
             y = feature,
             fill = factor(type, levels = c('sly_acr_pro', 'sly_acr_dis'),
                           labels = c('Proximal', 'Distal')))) +
  geom_bar(stat = 'identity',
           position = 'dodge',
           width = 1.5) +
  theme_classic() +
  coord_flip() +
  labs(x = 'Cluster',
       y = 'Feature',
       fill = 'Type') +
  scale_fill_manual(values = c('#d8b365','#01665e')) +
  facet_grid(cluster~., scales = 'free_y') +
  theme(strip.background = element_blank(),
        strip.text = element_blank())
dev.off()
#############################################


###### DESEQ2 DE PEAKS: TOMATO ######
# load ATAC-seq file
p <- 'f2/tomato'
load(paste(p, 'HS_peak.RData', sep = '/'))

resultsNames(ddsMF)

# export path
p_csv <- paste(p, 'table', sep = '/')
p_up <-  paste(p, 'up', sep = '/')
p_down <-  paste(p, 'down', sep = '/')
p_neg <-  paste(p, 'neg', sep = '/')

# dir.create(p_csv)
# dir.create(p_up)
# dir.create(p_down)
# dir.create(p_neg)

# export DE peaks
range <- data.frame(group = 'treatment',
                    treatment = c('1h', '6h'),
                    base = '0h')

for (i in 1:nrow(range)) {
  res <- results(ddsMF, contrast = c(range$group[i], range$treatment[i], range$base[i])) 
  res0.05 <- results(ddsMF, alpha = 0.05, contrast = c(range$group[i], range$treatment[i], range$base[i]))
  summary(res)
  summary(res0.05)
  
  write.csv(as.data.frame(res0.05), file = paste0(p_csv, '/', range$treatment[i], '_', range$base[i], '.csv'))
  res0.05 <- read.csv(paste0(p_csv, '/', range$treatment[i], '_', range$base[i], '.csv'))
  res0.05$Up_regulated <- 'NO' 
  res0.05$Up_regulated[res0.05$log2FoldChange >= 1.0 & res0.05$padj < 0.05] <- 'YES'
  res0.05$Down_regulated <- 'NO'
  res0.05$Down_regulated[res0.05$log2FoldChange <= -1.0 & res0.05$padj < 0.05] <- 'YES'
  res0.05$Negative <- 'NO' 
  res0.05$Negative[res0.05$log2FoldChange > -0.8 & res0.05$padj > 0.05] <- 'YES'
  res0.05$Negative[res0.05$log2FoldChange < 0.8 & res0.05$padj > 0.05] <- 'YES'
  write.csv(res0.05, paste0(p_csv, '/', range$treatment[i], '_', range$base[i], '.csv'), row.names = FALSE)
  
  
  # output gene list
  library(pgirmess)
  res0.05 <- read.csv(paste0(p_csv, '/', range$treatment[i], '_', range$base[i], '.csv'))
  
  up <- res0.05 %>% filter(Up_regulated =='YES') %>% dplyr::select(X)
  write.delim(up$X, file = paste0(p_up, '/', range$treatment[i], '_', range$base[i], '_u', '.txt'))
  
  dwn <- res0.05 %>% filter(Down_regulated =='YES') %>% dplyr::select(X)
  write.delim(dwn$X, file = paste0(p_down, '/', range$treatment[i], '_', range$base[i], '_d', '.txt'))
  
  neg <- res0.05 %>% filter(Negative =='YES') %>% dplyr::select(X)
  write.delim(neg$X, file = paste0(p_neg, '/', range$treatment[i], '_', range$base[i], '_n', '.txt'))
}


# summary DE peaks
p_plot <- 'f2/tomato/plot'
# dir.create(p_plot)

range <- data.frame(group = 'treatment',
                    treatment = c('1h', '6h'),
                    base = '0h')

sum <- data.frame()
for (i in 1:nrow(range)) {
  # up
  up <- nrow(read.table(paste(p_up, paste0(range$treatment[i], '_', range$base[i], '_u.txt'), sep = '/'), header = T))
  # down
  down <- nrow(read.table(paste(p_down, paste0(range$treatment[i], '_', range$base[i], '_d.txt'), sep = '/'), header = T))
  # neg
  neg <- nrow(read.table(paste(p_neg, paste0(range$treatment[i], '_', range$base[i], '_n.txt'), sep = '/'), header = T))
  
  sum.sub <- data.frame(group = range$treatment[i],
                        up_peaks = up,
                        down_peaks = down,
                        negative_peaks = neg)
  sum <- rbind(sum, sum.sub) 
}

sum <- sum %>% 
  gather(type, genes, -group) %>% 
  filter(type != 'negative_peaks') %>%
  mutate(genes = as.numeric(genes)) %>% 
  mutate(`genes (K)` = genes/1000)
sum$type <- factor(sum$type, levels = c('up_peaks', 'down_peaks'))
sum$group <- factor(sum$group, levels = c('1h', '6h'))

plt <- sum %>% 
  ggplot(.,
         aes(x = group,
             y = `genes (K)`,
             fill = type)) +
  geom_bar(stat = 'identity',
           position = 'dodge') +
  geom_text(aes(label = ifelse(type == 'up_peaks', genes, '')),
            vjust = -1, hjust = 1.1, size = 3) +
  geom_text(aes(label = ifelse(type == 'down_peaks', genes, '')),
            vjust = -1, hjust = -0.1, size = 3) +
  scale_fill_manual(values = c('#ca0025', '#91bfdb')) +
  labs(x = 'Genotype',
       y = '# of DE peaks (HS/CK)',
       fill = 'Type') +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text.x = element_text(size = 10, face = 'bold'),
        axis.text.y = element_text(size = 10, face = 'bold'),
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        axis.line = element_blank())

pdf(paste(p_plot, 'couts_DE_peaks.pdf', sep = '/'),
    width = 5, height = 5)
plt
dev.off()
#######################





###### CLUSTER DE PEAKS ######
# Hierarchical clustering
# all peaks
p <- 'f2/tomato'
p_plot <- 'f2/tomato/plot'
# dir.create(p_plot)
p_csv <- paste(p, 'table', sep = '/')
p_up <-  paste(p, 'up', sep = '/')
p_down <-  paste(p, 'down', sep = '/')
p_neg <-  paste(p, 'neg', sep = '/')

range <- data.frame(group = 'treatment',
                    treatment = c('1h', '6h'),
                    base = '0h')

sum <- c()
for (i in 1:nrow(range)) {
  # up
  up <- read.table(paste(p_up, paste0(range$treatment[i], '_', range$base[i], '_u.txt'), sep = '/'), header = T) %>% 
    pull()
  # down
  down <- read.table(paste(p_down, paste0(range$treatment[i], '_', range$base[i], '_d.txt'), sep = '/'), header = T) %>% 
    pull()
  
  sum.sub <- unique(c(up, down))
  sum <- c(sum, sum.sub)
}

all_peaks <- unique(sum) # 4427 DE peaks

norm_counts <- read.table(paste(p, 'normalized_counts.txt', sep = '/'), header = T) %>% 
  filter(peak %in% all_peaks)

# export
getwd()
write.table(norm_counts, 'DE_peaks_normalized_counts.txt', row.names = F, quote = F, sep = '\t')

# Create a matrix
df <- norm_counts %>% 
  column_to_rownames(var = 'peak')
hclust_matrix <- df %>% 
  as.matrix()

# z-score
sd <- apply(hclust_matrix, 1, sd, na.rm=TRUE)
hclust_matrix <- (hclust_matrix-rowMeans(hclust_matrix))/sd

write.table(hclust_matrix, 'DE_peaks_zscore.txt', quote = F, sep = '\t')


# if your data has been transformed or normalized
# please skip this step
# hclust_matrix <- hclust_matrix %>% 
# transpose the matrix so genes are as columns
# t() %>% 
# apply scaling to each column of the matrix (genes)
# scale() %>% 
# transpose back so genes are as rows again
# t()

peak_dist <- dist(hclust_matrix)


require(fastcluster)
peak_hclust <- fastcluster::hclust(peak_dist, method = "complete")

# The default `plot()` function can be used to produce a simple dendrogram
plot(peak_hclust, labels = FALSE)
abline(h = 4, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram
# we chose 5 clusters after screening the dendrogram plot

cutree(peak_hclust, k = 5)
peak_cluster <- cutree(peak_hclust, k = 5) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  dplyr::rename(peak = name, cluster = value) %>% 
  mutate(cluster = paste0('C', cluster))

head(peak_cluster)
table(peak_cluster$cluster)
# write.table(peak_cluster, 'DE_peaks_hclust.txt',
#             row.names = FALSE, quote = FALSE, sep = '\t')


# plot: scatter plot
meta <- read.table(paste(p, 'meta.txt', sep = '/'), header = T) %>% 
  mutate(name = paste(cultivar, treatment, replicate, sep = '_'))
peak_cluster <- read.table(paste(p, 'DE_peaks_hclust.txt', sep = '/'), header = T)
pd <- data.frame(hclust_matrix)

all(names(pd) == meta$sample) # if True, then pass

colnames(pd) <- meta$name
head(pd)

pd <- pd %>% 
  rownames_to_column(var = 'peak') %>% 
  merge(., peak_cluster) %>% 
  relocate(peak, cluster) %>% 
  gather(condition, zscore, -peak, -cluster) %>% 
  separate(., col = 'condition', into = c('cultivar', 'treatment', 'replicate'), sep = '_') %>% 
  group_by(cluster, cultivar, treatment, replicate) %>% 
  summarise(median_zscore = median(zscore))

pd$cluster <- factor(pd$cluster, levels = paste0('C', 1:5))

cluster <- data.frame(table(peak_cluster$cluster)) %>% 
  dplyr::rename(cluster = Var1,
                count = Freq)
pd <- merge(pd, cluster) %>% 
  mutate(name = paste0(cluster, ' (n=', count, ')'))
pd$name <- factor(pd$name, levels = unique(pd[order(pd$cluster),]$name))


# preliminary plot
plt <- ggplot(pd,
             aes(x = treatment,
                 y = median_zscore,
                 group = cultivar)) +
  geom_smooth(se = T,
              size = 0.5,
              color = 'grey') +
  geom_point(aes(color = cultivar),
             shape = 1,
             size = 2,
             alpha = 0.5) +
  labs(x = 'Time',
       y = 'Z-score',
       color = 'Cultivar') +
  scale_color_manual(values = c('black', '#f1a340', 'red')) +
  ylim(-2.0, 2.0) +
  facet_wrap(.~name, ncol = 5) +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text.x = element_text(size = 10, face = 'bold', angle = 40, hjust = 1),
        axis.text.y = element_text(size = 10, face = 'bold'),
        strip.background = element_rect(fill = alpha('#91bfdb', 0.5)),
        strip.text = element_text(size = 8, face = 'bold'),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = 'black'))

pdf(paste(p_plot, 'cluster.pdf', sep = '/'),
    width = 6, height = 1.75)
plt
dev.off()


# peak location/area
Txdb_gtf <- makeTxDbFromGFF('f2/tomato/Slycopersicum_796_ITAG5.0.gene.gff3')
head(peak_cluster)

data <- peak_cluster %>% 
  dplyr::select(peak) %>% 
  mutate(c = peak) %>% 
  separate(., col = 'c', into = c('Chr', 'Start', 'End'), sep = '\\.')

## make GRanges data
gr <- makeGRangesFromDataFrame(data, keep.extra.columns=T)

## annotate
peakAnno <- annotatePeak(gr,
                         tssRegion=c(-1500, 1500),
                         genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
                         TxDb=Txdb_gtf,
                         level ='gene')
## export peaks annotation
## only loss 6 peaks
peakAnno
peakAnno@anno
peakAnno_tb <- as_tibble(peakAnno@anno) %>% 
  mutate(`distanceToTSS (log10)` = log10((abs(distanceToTSS)/1000) + 1))

peakAnno_tb[which(peakAnno_tb$`distanceToTSS (log10)` < 0),]$`distanceToTSS (log10)` <- -peakAnno_tb[which(peakAnno_tb$`distanceToTSS (log10)` < 0),]$`distanceToTSS (log10)`


# plot1
plot_data <- peakAnno_tb %>% 
  dplyr::select(`distanceToTSS (log10)`)


plt <- ggplot(plot_data,
             aes(x = `distanceToTSS (log10)`)) + 
  geom_density(linewidth = 1) +
  labs(x = 'log10((Distance To TSS/1000) + 1)',
       y = 'Relative density') +
  geom_vline(aes(xintercept = log10(1)),
             linetype = 'dashed',
             color = 'black',
             linewidth = 0.25) +
  geom_vline(aes(xintercept = log10(1.5)),
             linetype = 'dashed',
             color = '#E69F00',
             linewidth = 1) +
  geom_vline(aes(xintercept = -log10(1.5)),
             linetype = 'dashed',
             color = '#E69F00',
             linewidth = 1) +
  # geom_vline(aes(xintercept = -log10(5)),
  #            linetype = 'dashed',
  #            color = 'black') +
  annotate('text',
           label = '±1.5KB',
           x = -log10(2),
           y = 2,
           size = 4,
           color = '#E69F00') +
  theme_classic() +
  scale_color_manual(values = c("#999999", "#E69F00")) +
  scale_x_continuous(limits = c(-log10(3), log10(3)),
                     breaks = c(-0.4, -0.2, 0, 0.2, 0.4),
                     labels = c(-0.4, -0.2, 0, 0.2, 0.4)) +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 10, face = 'bold'),
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        axis.line = element_blank())


pdf('merge_cluster_peak_distribution.pdf',
    width = 5, height = 3)
plt
dev.off()


# plot2
plot_data2 <- data.frame(peakAnno@annoStat)

plt <- ggplot(plot_data2,
             aes(x = Feature,
                 y = Frequency)) +
  geom_bar(stat = 'identity',
           position = 'dodge',
           fill = c("#E69F00", rep("#999999",6))
  ) +
  labs(x = 'Feature',
       y = 'Percentage') +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text.x = element_text(size = 10, face = 'bold', angle = 40, hjust = 1),
        axis.text.y = element_text(size = 10, face = 'bold'),
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        axis.line = element_blank())


pdf(paste(p_plot, 'merge_cluster_peak_distribution_ratio.pdf', sep = '/'),
    width = 3, height = 5)
plt
dev.off()
######################


###### SUMMARY RESULTS FROM KMER BINDETECT BY TOBIAS: TOMATO ######
# We only discussed differential binding between HS and Cntr (CK) under the same genotype.
# Therefore, the input1 can be 'ATAC_M82_1h_vs_ATAC_M82_0h', 'ATAC_M82_6h_vs_ATAC_M82_0h'.
# And the input2 can be promoter (proximal) or distal.
# differential binding events in 95% or 5% quantile
input_list_1 <- c('ATAC_M82_1h_vs_ATAC_M82_0h', 'ATAC_M82_6h_vs_ATAC_M82_0h')
input_list_2 <- c('proximal', 'distal')


# results vectors
full <- data.frame() # all events
top10 <- data.frame() # top 10 kmer family in each genotype


for (z in input_list_1) {
  for (j in input_list_2) {
    peak_path <- 'f2/tomato'
    peak_hclust <- read.delim(paste(peak_path, 'DE_peaks_hclust.txt', sep = '/'))
    input1 <- z
    input2 <- j
    
    if(input1 %in% input_list_1 & input2 %in% input_list_2){
      print(paste0('Processing: ', input1, ' under ', input2))
      path <- paste('f2/tomato/bind/', input2, input1, sep = '/')
      bound_data <- read.delim(paste(path, 'bindetect_results.txt', sep = '/'))
      
      # remove no bound under HS & CK
      bound_data2 <- bound_data[which(bound_data$HS_bound != 0 | bound_data$CK_bound != 0),]
      
      # remove prefix of TFBS_name
      bound_data2$kmer_motif_name <- str_remove(bound_data2$kmer_motif_name, '_')
      
      
      
      # add family name
      meme <- read_meme('f2/ArabidopsisDAPv1.meme')
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
      
      bound_data2$kmer_family <- ''
      bound_data3 <- data.frame()
      for (i in 1:length(family)) {
        sub <- bound_data2 %>% 
          filter(., grepl(paste0('_',family[i],'_'), kmer_motif_name)) %>% 
          mutate(kmer_family = family[i])
        
        # append
        bound_data3 <- rbind(bound_data3, sub)
      }
      
      # check no data loss
      # calculate all kmers to family
      # 282 kmers in 54 families
      sig_kmer <- read.delim(paste(path, paste0(input1, '_de_kmer.txt'), sep = '/')) %>% 
        dplyr::select(output_prefix, motif_id, total_tfbs, HS_bound, CK_bound, HS_CK_change) %>% 
        dplyr::rename(TFBS_name = output_prefix,
                      kmer_motif_name = motif_id,
                      total_HS_bound = HS_bound,
                      total_CK_bound = CK_bound)
      
      sum_sig_kmer <- bound_data3 %>% 
        dplyr::select(kmer_motif_name, kmer_family) %>% 
        distinct() %>% 
        group_by(kmer_family) %>%
        summarise(total_kmer = n())
      
      
      # add cluster
      bound_data3 <- bound_data3 %>%
        mutate(peak = paste(peak_chr, peak_start, peak_end, sep = '.')) %>% 
        merge(., peak_hclust, by = 'peak') %>%
        dplyr::rename(peak_cluster = cluster)
      
      # calculation: which cluster conquer the differential of a specific kmer
      # up bound: HS_CK_log2fc > 0
      # total bound: sum(HS_bound); sum(CK_bound)
      bound_counts <- bound_data3 %>%
        group_by(TFBS_name, peak_cluster) %>%
        filter(HS_CK_log2fc > 0) %>%  
        mutate(HS_bound_sum = sum(HS_bound),
               CK_bound_sum = sum(CK_bound)) %>% 
        dplyr::select(TFBS_name, peak_cluster, kmer_motif_name, kmer_family, HS_bound_sum, CK_bound_sum) %>% 
        distinct()
      
      # bound median log2fc
      bound_median_fc <- bound_data3 %>% 
        group_by(TFBS_name, peak_cluster) %>%
        filter(HS_CK_log2fc > 0) %>%
        mutate(median_HS_CK_log2fc = median(HS_CK_log2fc)) %>% 
        dplyr::select(TFBS_name, peak_cluster, kmer_motif_name, kmer_family, median_HS_CK_log2fc) %>% 
        distinct()
      
      # bound profile
      bound_profile <- merge(bound_counts, bound_median_fc)%>% 
        merge(., sig_kmer %>% dplyr::select(-TFBS_name)) %>% 
        mutate(hg = phyper(HS_bound_sum-1, total_HS_bound, total_tfbs-total_HS_bound, HS_bound_sum + CK_bound_sum, lower.tail = F)) %>% 
        mutate(supported = 'No')
      
      bound_profile[which(bound_profile$hg < 0.05),]$supported <- 'Yes'
      
      
      # we chose quantile 95, 85, 75, 50 and left
      qun <- data.frame(log2fc = quantile(bound_profile$median_HS_CK_log2fc,
                                          probs = c(0.5, 0.75, 0.85, 0.95))) %>% 
        rownames_to_column(var = 'quantile')
      
      # add quantile information
      bound_profile$quantile <- '<50%'
      for (i in 1:4) {
        cut_off <- qun[i,]$log2fc
        mark <- qun[i,]$quantile
        
        bound_profile[which(bound_profile$median_HS_CK_log2fc >= cut_off),]$quantile <- mark
      }
      
      # check which cut-off i want
      table(bound_profile$supported, bound_profile$quantile)
      
      # i chose 95% first
      bound_profile_95 <- bound_profile %>% 
        filter(quantile == '95%')
      
      # summary bound under 95%
      num_95 <- bound_profile_95 %>% 
        group_by(peak_cluster, kmer_family) %>% 
        summarise(number_of_bound_kmer = n(),
                  median_bound_log2fc = median(median_HS_CK_log2fc),
                  kmer_HS_bound_sum = sum(HS_bound_sum),
                  kmer_CK_bound_sum = sum(CK_bound_sum),
                  peak_total_HS_bound = sum(total_HS_bound),
                  peak_total_CK_bound = sum(total_CK_bound)) %>% 
        mutate(bound_type = 'positive',
               genotype = input1,
               region = input2,
               ratio_kmer_HS_bound_sum = kmer_HS_bound_sum/peak_total_HS_bound,
               ratio_kmer_CK_bound_sum = kmer_CK_bound_sum/peak_total_CK_bound)
      
      
      # down bound: HS_CK_log2fc > 0
      # total bound: sum(HS_bound); sum(CK_bound)
      bound_counts2 <- bound_data3 %>% 
        group_by(TFBS_name, peak_cluster) %>%
        filter(HS_CK_log2fc < 0) %>%  
        mutate(HS_bound_sum = sum(HS_bound),
               CK_bound_sum = sum(CK_bound)) %>% 
        dplyr::select(TFBS_name, peak_cluster, kmer_motif_name, kmer_family, HS_bound_sum, CK_bound_sum) %>% 
        distinct()
      
      # bound median log2fc
      bound_median_fc2 <- bound_data3 %>% 
        group_by(TFBS_name, peak_cluster) %>%
        filter(HS_CK_log2fc < 0) %>%
        mutate(median_HS_CK_log2fc = median(HS_CK_log2fc)) %>% 
        dplyr::select(TFBS_name, peak_cluster, kmer_motif_name, kmer_family, median_HS_CK_log2fc) %>% 
        distinct()
      
      # bound profile
      bound_profile2 <- merge(bound_counts2, bound_median_fc2) %>%
        merge(., sig_kmer %>% dplyr::select(-TFBS_name)) %>%
        mutate(hg = phyper(CK_bound_sum-1, total_CK_bound, total_tfbs-total_CK_bound, HS_bound_sum + CK_bound_sum, lower.tail = F)) %>% 
        mutate(supported = 'No')
      
      bound_profile2[which(bound_profile2$hg < 0.05),]$supported <- 'Yes'
      
      # we chose quantile 5, 15, 25, 50 and left
      qun2 <- data.frame(log2fc = quantile(bound_profile2$median_HS_CK_log2fc,
                                           probs = c(0.05, 0.15, 0.25, 0.5))) %>% 
        rownames_to_column(var = 'quantile')
      
      # add quantile information
      bound_profile2$quantile <- '>50%'
      for (i in 4:1) {
        cut_off <- qun2[i,]$log2fc
        mark <- qun2[i,]$quantile
        
        bound_profile2[which(bound_profile2$median_HS_CK_log2fc <= cut_off),]$quantile <- mark
      }
      
      # check which cut-off i want
      table(bound_profile2$supported, bound_profile2$quantile)
      
      # i chose 5% first
      bound_profile_5 <- bound_profile2 %>% 
        filter(quantile == '5%')
      
      # summary bound under 5%
      num_5 <- bound_profile_5 %>% 
        group_by(peak_cluster, kmer_family) %>% 
        summarise(number_of_bound_kmer = n(),
                  median_bound_log2fc = median(median_HS_CK_log2fc),
                  kmer_HS_bound_sum = sum(HS_bound_sum),
                  kmer_CK_bound_sum = sum(CK_bound_sum),
                  peak_total_HS_bound = sum(total_HS_bound),
                  peak_total_CK_bound = sum(total_CK_bound)) %>% 
        mutate(bound_type = 'negative',
               genotype = input1,
               region = input2,
               ratio_kmer_HS_bound_sum = kmer_HS_bound_sum/peak_total_HS_bound,
               ratio_kmer_CK_bound_sum = kmer_CK_bound_sum/peak_total_CK_bound)
      
      
      # combine positive and negative bound type
      # make new TF Family name
      bound_profile_combine <- rbind(num_5, num_95) %>% 
        distinct() %>% 
        merge(., sum_sig_kmer) %>% 
        mutate(bound_kmer_ratio = number_of_bound_kmer/total_kmer,
               kmer_name = paste0(kmer_family, ' (', total_kmer, ')'))
      
      # add missing cluster
      model <- data.frame(peak_cluster = paste0('C', 1:5),
                          kmer_name = rep(unique(bound_profile_combine$kmer_name), each = 10),
                          bound_type = rep(c('positive', 'negative'), each = 5))
      bound_profile_combine <- full_join(model, bound_profile_combine)
      
      
      # add factor to cluster
      bound_profile_combine$peak_cluster <- factor(bound_profile_combine$peak_cluster, levels = paste0('C', 1:5))
      bound_profile_combine$bound_type <- factor(bound_profile_combine$bound_type, levels = c('positive', 'negative'),
                                                 labels = c('P', 'N'))
      
      
      # select top10
      kmer_top <- bound_profile_combine %>%
        filter(!is.na(kmer_family)) %>% 
        arrange(desc(ratio_kmer_HS_bound_sum), desc(median_bound_log2fc)) %>%  # Arrange by descending kmer_bound_value
        slice_head(n = 10) %>% 
        pull(kmer_family)
      
      bound_profile_combine_top <- bound_profile_combine %>% 
        filter(kmer_family %in% kmer_top)
      
      # append the data
      full <- rbind(full, bound_profile_combine)
      top10 <- rbind(top10, bound_profile_combine_top)
      
      # 
      # pdf(paste0('f2/tomato/plot/', input1, '_', input2, '_bound.pdf'),
      #     width = 10, height = 6)
      # plot <- bound_profile_combine %>%
      #   ggplot(.,
      #          aes(x = bound_type,
      #              y = kmer_name,
      #              color = median_bound_log2fc,
      #              size = bound_kmer_ratio)) +
      #   geom_point() +
      #   scale_color_gradient2(low = 'blue',
      #                         mid = 'white',
      #                         high = 'red',
      #                         midpoint = 0) +
      #   labs(x = 'Cluster',
      #        y = 'TF Family (kmer)',
      #        color = 'Median HS/CK log2FC',
      #        size = 'DE bound kmer/total kmer') +
      #   scale_x_discrete(position = "top") +
      #   facet_grid(.~peak_cluster) +
      #   theme_classic() +
      #   theme(axis.title = element_text(size = 12, face = 'bold'),
      #         axis.text = element_text(size = 10, face = 'bold'),
      #         strip.background = element_blank(),
      #         strip.text = element_text(size = 12, face = 'bold'),
      #         strip.placement = 'outside')
      # 
      # grid.draw(plot)
      # dev.off()
      
      
      # setwd('f2/tomato')
      # pdf(paste0('plot/bound_quantile_support/', input1, '_', input2, '_bound.pdf'),
      #     width = 10, height = 7)
      # plot <- bound_profile_combine %>%
      #   ggplot(.,
      #          aes(x = bound_type,
      #              y = kmer_name,
      #              color = median_bound_log2fc,
      #              size = bound_kmer_ratio)) +
      #   geom_point() +
      #   scale_color_gradient(low = 'blue',
      #                        high = 'red') +
      #   labs(x = 'Cluster',
      #        y = 'TF Family (kmer)',
      #        color = 'Median HS/CK log2FC',
      #        size = 'DE bound kmer/total kmer') +
      #   scale_x_discrete(position = "top") +
      #   facet_grid(.~peak_cluster) +
      #   theme_classic() +
      #   theme(axis.title = element_text(size = 12, face = 'bold'),
      #         axis.text = element_text(size = 10, face = 'bold'),
      #         strip.background = element_blank(),
      #         strip.text = element_text(size = 12, face = 'bold'),
      #         strip.placement = 'outside')
      # 
      # grid.draw(plot)
      # dev.off()
      
      
      
      # export 95% and 5% bound bed
      # you can check whether the output is the same with bound_profile_95/5
      bound_data3 %>%
        filter(HS_CK_log2fc > 0) %>%
        merge(bound_profile_95 %>% dplyr::select(kmer_motif_name, kmer_family, peak_cluster), .) %>%
        dplyr::select(TFBS_chr, TFBS_start, TFBS_end, kmer_motif_name, HS_CK_log2fc, TFBS_strand) %>%
        write.table(., paste0('f2/tomato/tfcomb/bed/', input1, '_', input2, '_up.bed'),
                    row.names = F, col.names = F, quote = F, sep = '\t')
      
      
      bound_data3 %>%
        filter(HS_CK_log2fc < 0) %>%
        merge(bound_profile_5 %>% dplyr::select(kmer_motif_name, kmer_family, peak_cluster), .) %>%
        dplyr::select(TFBS_chr, TFBS_start, TFBS_end, kmer_motif_name, HS_CK_log2fc, TFBS_strand) %>%
        write.table(., paste0('f2/tomato/tfcomb/bed/', input1, '_', input2, '_down.bed'),
                    row.names = F, col.names = F, quote = F, sep = '\t')
      
      gc()
    }
  }
}

save(top10, full, file = 'f2/tomato/bind/soft_mode_all.RData')



# visualization:proximal
top10$genotype <- factor(top10$genotype, levels = input_list_1,
                         labels = c('ATAC_M82_1h_vs_ATAC_M82_0h', 'ATAC_M82_6h_vs_ATAC_M82_0h'))


## Supporting Figure 5-F
pdf(paste0('f2/tomato/plot/', 'proximal', '_bound_soft_mode.pdf'),
    width = 12, height = 8)
plt <- top10 %>%
  filter(region == 'proximal') %>%
  ggplot(.,
         aes(x = bound_type,
             y = kmer_family,
             color = median_bound_log2fc,
             size = ratio_kmer_HS_bound_sum)) +
  geom_point() +
  scale_color_gradient2(low = 'blue',
                        mid = 'white',
                        high = 'red',
                        midpoint = 0,
                        limits = c(-2,2),
                        oob = scales::squish) +
  labs(x = 'Cluster',
       y = 'TF Family (kmer)',
       color = 'Median HS/CK log2FC',
       size = 'DE bound under HS/total bound') +
  scale_x_discrete(position = "top") +
  facet_grid(genotype~peak_cluster, scales = 'free_y', switch = 'y') +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 10, face = 'bold'),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12, face = 'bold'),
        strip.text.y.left = element_text(size = 12, face = 'bold', angle = 0),
        strip.placement = 'outside')
grid.draw(plt)
dev.off()


# visualization:distal
top10$genotype <- factor(top10$genotype, levels = input_list_1,
                         labels = c('ATAC_M82_1h_vs_ATAC_M82_0h', 'ATAC_M82_6h_vs_ATAC_M82_0h'))
pdf(paste0('f2/tomato/plot/', 'distal', '_bound_soft_mode.pdf'),
    width = 12, height = 8)
plt <- top10 %>%
  filter(region == 'distal') %>%
  ggplot(.,
         aes(x = bound_type,
             y = kmer_family,
             color = median_bound_log2fc,
             size = ratio_kmer_HS_bound_sum)) +
  geom_point() +
  scale_color_gradient2(low = 'blue',
                        mid = 'white',
                        high = 'red',
                        midpoint = 0,
                        limits = c(-2,2),
                        oob = scales::squish) +
  labs(x = 'Cluster',
       y = 'TF Family (kmer)',
       color = 'Median HS/CK log2FC',
       size = 'DE bound under HS/total bound') +
  scale_x_discrete(position = "top") +
  facet_grid(genotype~peak_cluster, scales = 'free_y', switch = 'y') +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 10, face = 'bold'),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12, face = 'bold'),
        strip.text.y.left = element_text(size = 12, face = 'bold', angle = 0),
        strip.placement = 'outside')
grid.draw(plt)
dev.off()
#####################################


###### KMER DISTRIBUTION: TOMATO ######
proximal_motif <- read.delim('f2/tomato/kmer/sim/sly_acr_pro/proximal_pcc.txt')
distal_motif <- read.delim('f2/tomato/kmer/sim/sly_acr_dis/distal_pcc.txt')
input_list_1 <- c('ATAC_M82_1h_vs_ATAC_M82_0h', 'ATAC_M82_6h_vs_ATAC_M82_0h')
input_list_2 <- c('proximal', 'distal')
input_list_3 <- c('soft')
p_path <- 'f2/tomato/distribution_all'
p_path2 <- 'f2/tomato/distribution_all/plot'
Txdb_gtf <- makeTxDbFromGFF('f2/tomato/Slycopersicum_796_ITAG5.0.gene.gff3')

# dir.create(p_path)
# dir.create(p_path2)

for (z in input_list_1) {
  for (j in input_list_2) {
    for (h in input_list_3) {
      input1 <- z
      input2 <- j
      input3 <- h
      
      # load top10 file in strict or soft mode
      # you will get top10 and full dataframe from kmer bindetect/2
      load(paste0('f2/tomato/bind/', input3, '_mode_all.RData'))
      
      
      if(input1 %in% input_list_1 & input2 %in% input_list_2){
        print(paste0('Processing: ', input1, ' under ', input2, ' (', input3, ' mode)'))
        path <- paste('f2/tomato/bind/', input2, input1, sep = '/')
        bound_data <- read.delim(paste(path, 'bindetect_results.txt', sep = '/'))
        
        # remove no bound under HS & CK
        bound_data2 <- bound_data[which(bound_data$HS_bound != 0 | bound_data$CK_bound != 0),]
        
        # remove prefix of TFBS_name
        bound_data2$kmer_motif_name <- str_remove(bound_data2$kmer_motif_name, '_')
        
        # check all motifs are included in corresponding kmer pool
        if (input2 == 'proximal'){
          kmer_pool <- proximal_motif
          if(!all(bound_data2$kmer_motif_name %in% kmer_pool$kmers)){
            stop('Not all motifs in bound data exist in user-provided kmer pool.\n')
          }
        } else{
          kmer_pool <- distal_motif
          all(bound_data2$kmer_motif_name %in% kmer_pool$kmers)
          if(!all(bound_data2$kmer_motif_name %in% kmer_pool$kmers)){
            stop('Not all motifs in bound data exist in user-provided kmer pool.\n')
          }
        }
        
        
        # add family name
        meme <- read_meme('f2/ArabidopsisDAPv1.meme')
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
        
        bound_data2$kmer_family <- ''
        bound_data3 <- data.frame()
        for (i in 1:length(family)) {
          sub <- bound_data2 %>% 
            filter(., grepl(paste0('_',family[i],'_'), kmer_motif_name)) %>% 
            mutate(kmer_family = family[i])
          
          # append
          bound_data3 <- rbind(bound_data3, sub)
        }
        
        
        # annotation TFBS distance
        # extract Chr, Start, End, Peak name
        data <- bound_data3[,c('TFBS_chr', 'TFBS_start', 'TFBS_end', 'kmer_motif_name')]
        names(data) <- c('Chr', 'Start', 'End', 'Motif_name')
        head(data)
        
        # make GRanges data
        gr <- makeGRangesFromDataFrame(data, keep.extra.columns=T)
        
        # annotate
        peakAnno <- annotatePeak(gr,
                                 tssRegion=c(-1500, 1500),
                                 genomicAnnotationPriority = c("Exon", "Intron", "5UTR", "3UTR", "Promoter", "Downstream", "Intergenic"),
                                 TxDb=Txdb_gtf,
                                 level ='gene')
        
        # export peaks annotation
        peakAnno
        peakAnno@anno
        peakAnno_tb <- as_tibble(peakAnno@anno) %>% 
          mutate(`distanceToTSS (log10)` = log10((abs(distanceToTSS)/1000) + 1))
        
        peakAnno_tb[which(peakAnno_tb$distanceToTSS < 0),]$`distanceToTSS (log10)` <- -peakAnno_tb[which(peakAnno_tb$distanceToTSS < 0),]$`distanceToTSS (log10)`
        
        
        bound_data3_dis <- peakAnno_tb %>% 
          dplyr::select(seqnames, start, end, Motif_name, geneId, distanceToTSS, `distanceToTSS (log10)`) %>% 
          dplyr::rename(TFBS_chr = seqnames,
                        TFBS_start = start,
                        TFBS_end = end,
                        kmer_motif_name = Motif_name) %>% 
          merge(., bound_data3)
        
        
        write.table(bound_data3_dis, paste0(p_path, '/', input3, '/', input2, '/',
                                            input3, '_', input2, '_', input1, '_peak_distribution.txt'),
                    row.names = F, quote = F, sep = '\t')
        
        
        # plot
        peak_path <- 'f2/tomato'
        peak_hclust <- read.delim(paste(peak_path, 'DE_peaks_hclust.txt', sep = '/'))
        
        
        bound_data3_dis <- bound_data3_dis %>%
          mutate(peak = paste(peak_chr, peak_start, peak_end, sep = '.')) %>% 
          merge(., peak_hclust, by = 'peak') %>% 
          dplyr::rename(peak_cluster = cluster)
        
        plot_data <- rbind(bound_data3_dis %>% 
                             dplyr::select(distanceToTSS,
                                           `distanceToTSS (log10)`,
                                           peak_cluster,
                                           kmer_family))
        
        
        family <- top10 %>% 
          filter(region == input2,
                 genotype == input1)
        kmer_family <- family %>% 
          dplyr::select(kmer_family) %>% 
          distinct() %>% 
          pull()
        peak_cluster <- family %>% 
          dplyr::select(peak_cluster) %>% 
          distinct()
        peak_c <- as.character(peak_cluster$peak_cluster)
        
        
        
        
        dis.sum <- data.frame()
        for (m in kmer_family) {
          df <- plot_data %>% filter(kmer_family %in% kmer_family &
                                       peak_cluster %in% peak_c)
          
          # add 1000bp to make -1kb to 0.5kb in range(0,1500)
          df$Preferential_Position <- as.numeric(df$distanceToTSS) + 1000
          df <- df %>%
            group_by(Preferential_Position, kmer_family) %>% 
            summarise(count = n())
          names(df) <- c('location', 'kmer_family', 'tar')
          
          
          window.tar <- data.frame(location = 1:1500)
          window.tar <- left_join(window.tar, df %>% 
                                    filter(kmer_family == m) %>% 
                                    dplyr::select(-kmer_family), by = 'location')
          window.tar[is.na(window.tar$tar),]$tar <- 0
          window.tar$tar <- window.tar$tar + 1 
          
          
          
          df2 <- plot_data %>% filter(kmer_family %in% kmer_family)
          # add 1000bp to make -1kb to 0.5kb in range(0,1500)
          df2$Preferential_Position <- as.numeric(df2$distanceToTSS) + 1000
          df2 <- df2 %>%
            group_by(Preferential_Position, kmer_family) %>% 
            summarise(count = n())
          names(df2) <- c('location', 'kmer_family', 'beg')
          
          
          window.beg <- data.frame(location = 1:1500)
          window.beg <- full_join(window.beg, df2 %>% 
                                    filter(kmer_family == m) %>% 
                                    dplyr::select(-kmer_family), by = 'location')
          window.beg[is.na(window.beg$beg),]$beg <- 0
          window.beg$beg <- window.beg$beg + 1 
          
          window <- merge(window.tar, window.beg, by = 'location')
          
          
          ## culculation
          # all
          # bin:100; sliding window:25; median
          tar.win <- rollapply(window$tar, width = 100, mean, by = 25, partial = FALSE)
          beg.win <- rollapply(window$beg, width = 100, mean, by = 25, partial = FALSE)
          
          
          
          dis <- data.frame(bin = seq(from = 100, to = 1500, by = 25))
          dis$tar <- tar.win
          dis$beg <- beg.win
          
          
          dis <- dis %>% 
            mutate(zscore_tar = (tar - mean(tar))/sd(tar),
                   zscore_beg = (beg - mean(beg))/sd(beg))
          dis <- dis %>% 
            gather(., key = 'group', value = 'zscore', -bin, -tar, -beg) %>% 
            mutate(kmer_family = m)
          
          # append dis
          dis.sum <- rbind(dis.sum, dis)
        }
        
        write.table(dis.sum, paste0(p_path, '/', input3, '/', input2, '/',
                                    input3, '_', input2, '_', input1, '_peak_zscore.txt'),
                    row.names = F, quote = F, sep = '\t')
        
        
        plt <- dis.sum %>% 
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
          facet_grid(.~factor(group, levels = c('zscore_tar', 'zscore_beg'), labels = c('Target', 'background'))) +
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
        
        
        ## Supporting Figure 5-G
        pdf(paste0(p_path2, '/',
                   input3, '_', input2, '_', input1, '_peak_zscore.pdf'),
            width = 7, height = 2)
        print(plt)
        dev.off()
      }
    }
  }
}
######################################
