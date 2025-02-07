############
# load packages
library(dplyr)
library(tidyverse)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(Biostrings)
library(stringr)
library(GenomicFeatures)
library(ChIPseeker)
library(DiffBind)
library(BiocParallel)
library(DESeq2)
############


############################################# Figure 1: F,G,H,I
############################################# Supporting Figure 2: A,B,C,D,E
###### LOAD ATAC-SEQ COUNT MATRIX ######
# cat bed1.bed bed2.bed bedN.bed (...) | sort -k1,1 -k2,2n | bedtools merge -i - > merged.bed
# cp merged.bed to the folder of all bam files
# FeatureCounts with the bed file
# You can get the counts files to peak

# converting bed file to saf for featureCounts
# awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' merged.bed > merged.saf

# featureCounts -T 10 -a merged.saf -F SAF -o counts.txt -p -B -C *.bam

# setwd("~/folder_to_all_the_data")
p <- 'f1/counts/split' # split replicates
df <- read.table(paste(p, 'counts.txt', sep = '/'), header=T)

# generate mata files
meta <- data.frame(names = names(df[,-1])) %>%
  separate(., col = names, into = c('genotype', 'treatment', 'type', 'file')) %>%
  dplyr::select(-type, -file) %>%
  mutate(replicate = rep(c(1,2), 8))
meta$treatment <- str_remove(meta$treatment, '\\d')
meta$genotype <- factor(meta$genotype, levels = c('Tak1', 'hsfa', 'hsfb', 'dko'))
meta$treatment <- factor(meta$treatment, levels = c('CK', 'HS'))

meta <- meta[order(meta$genotype),] %>%
  mutate(sample = paste0('s', 1:nrow(.)))
head(meta)
write.table(meta, paste(p, 'meta.txt', sep = '/'), row.names = F, quote = F, sep = '\t')

meta <- read.table(paste(p, 'meta.txt', sep = '/'), header=T)
df <- df %>% 
  column_to_rownames(var = 'Peak')


names(df) <- meta$sample
dim(df) # 40873 total peaks and 16 samples
meta <- meta %>% 
  column_to_rownames(var = 'sample')
#######################


####### REMOVE MISSING PEAKS #######
# remove rows where not a single sample has more than 50 reads
df <- df[apply(df, 1, max) > 50,]
dim(df) # 40872 total peaks and 16 samples
#######################


###### DESeq2 and DATA NORMALIZATION ######
dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = meta,
                              design = ~ genotype + treatment + genotype:treatment)
ddsMF <- dds
ddsMF$group <- factor(paste(ddsMF$genotype, ddsMF$treatment, sep = '_'))
design(ddsMF) <- ~ group
head(ddsMF)
ddsMF$group <- relevel(ddsMF$group, ref = "Tak1_CK")
ddsMF <- DESeq(ddsMF)

vsdMF <- vst(ddsMF, blind = FALSE) # sample numbers < 30, VST normalization method
head(assay(vsdMF), 5)
normalized.count <- data.frame(assay(vsdMF)) %>% 
  rownames_to_column(., var = 'peak')
write.table(normalized.count, paste(p, 'normalized_counts.txt', sep = '/'),
            row.names = F, quote = F, sep = '\t')
#######################


###### PCA ######
pcaData <- plotPCA(vsdMF, intgroup = c("treatment", "genotype"), returnData = TRUE)
percenVar <- round(100 * attr(pcaData, "percentVar"))


## Figure 1-F
pdf('f1/plots/PCA.pdf',
    width = 4, height = 4)
ggplot(pcaData, aes(PC1, PC2,
                    color = factor(genotype, levels = c('Tak1', 'hsfa', 'hsfb', 'dko'), labels = c('Tak-1', 'hsfa1', 'hsfb1', 'dko')),
                    shape = factor(treatment, levels = c('CK', 'HS'), labels = c('Cntr', 'HS')))) + 
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percenVar[1], "% variance")) +
  ylab(paste0("PC2: ", percenVar[2], "% variance")) +
  labs(color = 'genotype', shape = 'treatment') +
  scale_color_manual(values = c('#f6efe9', '#bdc9e1', '#67a9cf', '#02818a')) +
  theme_classic() +
  coord_fixed() +
  theme(axis.text = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 12, face = 'bold'),
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        axis.line = element_blank())
dev.off()
#######################


###### GAIN LOSS OVERLAP PEAKS ######
# follow peak_intersect.sh to get the gain/loss/overlap peaks from each condition
# genotype between HS and CK
# gain and loss
p2 <- 'narrowpeak'
all <- read.table(paste(p2, 'all_peaks.txt', sep = '/'), header = F) %>% 
  dplyr::rename(counts = V1,
                condition = V2) %>% 
  filter(condition != 'total') %>% 
  separate(., col = 'condition', into = c('genotype_treatment', 'end'), sep = '.peaks') %>% 
  dplyr::select(-end) %>%
  separate(., col = 'genotype_treatment', into = c('genotype', 'treatment')) %>% 
  mutate(counts = as.numeric(counts),
         type = 'all') %>% 
  mutate(`counts (K)` = counts/1000)

gain <- read.table(paste(p2, 'gain_peaks.txt', sep = '/'), header = F) %>% 
  dplyr::rename(counts = V1,
                condition = V2) %>% 
  filter(condition != 'total') %>% 
  separate(., col = 'condition',
           into = c('type', 'genotype', 'treatment', 'vs', 'end'), sep = '_') %>% 
  dplyr::select(-vs, -end) %>%
  mutate(counts = as.numeric(counts)) %>% 
  mutate(`counts (K)` = counts/1000)

loss <- read.table(paste(p2, 'loss_peaks.txt', sep = '/'), header = F) %>% 
  dplyr::rename(counts = V1,
                condition = V2) %>% 
  filter(condition != 'total') %>% 
  separate(., col = 'condition',
           into = c('type', 'genotype', 'treatment', 'vs', 'end'), sep = '_') %>% 
  dplyr::select(-vs, -end) %>%
  mutate(counts = as.numeric(counts)) %>% 
  mutate(`counts (K)` = counts/1000)

intersect <- read.table(paste(p2, 'intersect_peaks.txt', sep = '/'), header = F) %>% 
  dplyr::rename(counts = V1,
                condition = V2) %>% 
  filter(condition != 'total') %>% 
  separate(., col = 'condition',
           into = c('type', 'genotype', 'treatment', 'vs', 'end'), sep = '_') %>% 
  dplyr::select(-vs, -end) %>%
  mutate(counts = as.numeric(counts)) %>% 
  mutate(`counts (K)` = counts/1000)


# combine data
glt <- rbind(all %>% 
               dplyr::select(genotype, treatment, type, counts),
             gain %>% 
               dplyr::select(genotype, treatment, type, counts),
             loss %>% 
               dplyr::select(genotype, treatment, type, counts),
             intersect %>% 
               dplyr::select(genotype, treatment, type, counts)) %>% 
  mutate(counts = as.numeric(counts)) %>% 
  mutate(`counts (K)` = counts/1000,
         name = paste0(genotype, '(', treatment, '/CK', ')'))

glt <- glt %>% 
  merge(., all %>% 
          dplyr::select(-type, -`counts (K)`) %>% 
          dplyr::rename(total_counts = counts)) %>% 
  mutate(ratio = paste0(round((counts/total_counts), 2)))


glt[which(glt$type == 'loss'),]$`counts (K)` <- -glt[which(glt$type == 'loss'),]$`counts (K)`
glt

glt$genotype <- factor(glt$genotype, levels = c('Tak1', 'hsfa', 'hsfb', 'dko'), labels = c('Tak-1', 'hsfa1', 'hsfb1', 'dko'))
glt$name <- factor(glt$name, levels = unique(glt[order(glt$genotype),]$name))
glt$counts_ratio <- paste0(glt$counts, ' (', glt$ratio, ')')


plt <- ggplot(glt %>% 
               filter(type %in% c('gain', 'loss')),
             aes(x = name,
                 y = `counts (K)`,
                 fill = type)) +
  geom_bar(stat = 'identity',
           width = 0.5) +
  geom_text(aes(label = ifelse(`counts (K)` > 0, counts_ratio, '')),
            hjust = -0.1, size = 2) +
  geom_text(aes(label = ifelse(`counts (K)` < 0 & name != 'dko(HS/CK)', counts_ratio, '')),
            hjust = 1.1, size = 2) +
  geom_text(aes(label = ifelse(`counts (K)` < 0 & name == 'dko(HS/CK)', counts_ratio, '')),
            hjust = 0.1, size = 2) +
  scale_fill_manual(values = c('#fc8d59', '#91bfdb')) +
  scale_y_continuous(limits = c(-10, 10)) +
  labs(x = '',
       y = 'counts (K)',
       fill = 'Type') +
  coord_flip() +
  theme_classic() +
  theme(axis.title = element_text(size = 10, face = 'bold'),
        axis.text = element_text(size = 8, face = 'bold'),
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        axis.line = element_blank())


## Figure 1-G
pdf(paste('f1/plots/genotype_HS_CK_peak_gain_loss.pdf', sep = '/'),
    width = 5, height = 1.5)
plt
dev.off()
#######################


Txdb_gtf <- makeTxDbFromGFF('f1/MpTak_v6.1r1.gff')


###### DESeq2 DE PEAKS ######
# load ATAC-seq file
p <- 'f1/counts/split'

resultsNames(ddsMF)

# create path to save results
p_csv <- paste(p, 'table', sep = '/')
p_up <-  paste(p, 'up', sep = '/')
p_down <-  paste(p, 'down', sep = '/')
p_neg <-  paste(p, 'neg', sep = '/')

# dir.create(p_csv)
# dir.create(p_up)
# dir.create(p_down)
# dir.create(p_neg)

# export DE peaks
range <- data.frame(group = 'group',
                    treatment = rep(paste0(c('Tak1', 'hsfa', 'hsfb', 'dko'), '_HS'), each = 4),
                    base = paste0(c('Tak1', 'hsfa', 'hsfb', 'dko'), '_CK'))

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
range <- data.frame(treatment = c(paste0(c('Tak1', 'hsfa', 'hsfb', 'dko'), '_HS'),
                                  paste0(c('hsfa', 'hsfb', 'dko'), '_HS')),
                    base = c(rep('Tak1_CK', 4),
                             paste0(c('hsfa', 'hsfb', 'dko'), '_CK')))

sum <- data.frame()
for (i in 1:4) {
  if(i == 1){
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
  }else{
    z <- i + 3
    # combine /Tak1 and /genotype
    # /Tak1
    # up
    up <- length(unique(c(read.table(paste(p_up, paste0(range$treatment[i], '_', range$base[i], '_u.txt'), sep = '/'), header = T) %>% pull(),
                          read.table(paste(p_up, paste0(range$treatment[z], '_', range$base[z], '_u.txt'), sep = '/'), header = T) %>% pull())))
    # down
    down <- length(unique(c(read.table(paste(p_down, paste0(range$treatment[i], '_', range$base[i], '_d.txt'), sep = '/'), header = T) %>% pull(),
                            read.table(paste(p_down, paste0(range$treatment[z], '_', range$base[z], '_d.txt'), sep = '/'), header = T) %>% pull())))
    # neg
    neg <- length(unique(c(read.table(paste(p_neg, paste0(range$treatment[i], '_', range$base[i], '_n.txt'), sep = '/'), header = T) %>% pull(),
                           read.table(paste(p_neg, paste0(range$treatment[z], '_', range$base[z], '_n.txt'), sep = '/'), header = T) %>% pull())))
    
    sum.sub <- data.frame(group = range$treatment[i],
                          up_peaks = up,
                          down_peaks = down,
                          negative_peaks = neg)
    sum <- rbind(sum, sum.sub)
  }
}

sum <- sum %>% 
  gather(type, genes, -group) %>% 
  filter(type != 'negative_peaks') %>%
  mutate(genes = as.numeric(genes)) %>% 
  mutate(`genes (K)` = genes/1000)
sum$type <- factor(sum$type, levels = c('up_peaks', 'down_peaks'))

sum$group <- factor(sum$group, levels = paste0(c('Tak1', 'hsfa', 'hsfb', 'dko'), '_HS'),
                    labels = c('Tak-1', 'hsfa1', 'hsfb1', 'dko'))

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


## Figure 1-H
pdf(paste('f1/plots/couts_DE_peaks.pdf', sep = '/'),
    width = 5, height = 5)
plt
dev.off()
#######################


###### CLUSTER DE PEAKS ######
# Hierarchical clustering
# all DE peaks
p <- 'f1/counts/split'
p_plot <- 'f1/plots'
p_csv <- paste(p, 'table', sep = '/')
p_up <-  paste(p, 'up', sep = '/')
p_down <-  paste(p, 'down', sep = '/')
p_neg <-  paste(p, 'neg', sep = '/')
range <- data.frame(dir = p,
                    treatment = c(paste0(c('Tak1', 'hsfa', 'hsfb', 'dko'), '_HS'),
                                  paste0(c('hsfa', 'hsfb', 'dko'), '_HS')),
                    base = c(rep('Tak1_CK', 4),
                             paste0(c('hsfa', 'hsfb', 'dko'), '_CK')))

sum <- c()
for (i in 1:4) {
  if(i == 1){
    # up
    up <- read.table(paste(p_up, paste0(range$treatment[i], '_', range$base[i], '_u.txt'), sep = '/'), header = T) %>% 
      pull()
    # down
    down <- read.table(paste(p_down, paste0(range$treatment[i], '_', range$base[i], '_d.txt'), sep = '/'), header = T) %>% 
      pull()
    
    sum.sub <- unique(c(up, down))
    sum <- c(sum, sum.sub) 
  }else{
    z <- i + 3
    # combine /Tak1 and /genotype
    # /Tak1
    # up
    up <- unique(c(read.table(paste(p_up, paste0(range$treatment[i], '_', range$base[i], '_u.txt'), sep = '/'), header = T) %>% pull(),
                   read.table(paste(p_up, paste0(range$treatment[z], '_', range$base[z], '_u.txt'), sep = '/'), header = T) %>% pull()))
    # down
    down <- unique(c(read.table(paste(p_down, paste0(range$treatment[i], '_', range$base[i], '_d.txt'), sep = '/'), header = T) %>% pull(),
                     read.table(paste(p_down, paste0(range$treatment[z], '_', range$base[z], '_d.txt'), sep = '/'), header = T) %>% pull()))
    
    sum.sub <- unique(c(up, down))
    sum <- c(sum, sum.sub)
  }
}

all_peaks <- unique(sum) 
length(all_peaks) # 5946 DE peaks

# use normalized peak read counts to do clustering
norm_counts <- read.table(paste(p, 'normalized_counts.txt', sep = '/'), header = T) %>% 
  filter(peak %in% all_peaks)

# export
getwd()
write.table(norm_counts, paste(p, 'DE_peaks_normalized_counts.txt', sep = '/'), row.names = F, quote = F, sep = '\t')

# create a matrix
df <- norm_counts %>% 
  column_to_rownames(var = 'peak')
hclust_matrix <- df %>% 
  as.matrix()

# z-score
sd <- apply(hclust_matrix, 1, sd, na.rm=TRUE)
hclust_matrix <- (hclust_matrix-rowMeans(hclust_matrix))/sd

write.table(hclust_matrix, paste(p, 'DE_peaks_zscore.txt', sep = '/'), quote = F, sep = '\t')



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
abline(h = 6.7, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram
# we chose 15 clusters after screening the dendrogram plot

cutree(peak_hclust, k = 15)
peak_cluster <- cutree(peak_hclust, k = 15) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  dplyr::rename(peak = name, cluster = value) %>% 
  mutate(cluster = paste0('C', cluster))

# write.table(peak_cluster, 'DE_peaks_hclust.txt',
#             row.names = FALSE, quote = FALSE, sep = '\t')


# plot: scatter plot
meta <- read.table(paste(p, 'meta.txt', sep = '/'), header = T) %>% 
  mutate(name = paste(genotype, treatment, replicate, sep = '_'))
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
  separate(., col = 'condition', into = c('genotype', 'treatment', 'replicate'), sep = '_') %>% 
  group_by(cluster, genotype, treatment, replicate) %>% 
  summarise(median_zscore = median(zscore))

pd$cluster <- factor(pd$cluster, levels = paste0('C', 1:15))
pd$genotype <- factor(pd$genotype, levels = c('Tak1', 'hsfa', 'hsfb', 'dko'), labels = c('Tak-1', 'hsfa', 'hsfb', 'dko'))
pd$treatment <- factor(pd$treatment, levels = c('CK', 'HS'), labels = c('Cntr', 'HS'))

cluster <- data.frame(table(peak_cluster$cluster)) %>% 
  dplyr::rename(cluster = Var1,
                count = Freq)
pd <- merge(pd, cluster) %>% 
  mutate(name = paste0(cluster, ' (n=', count, ')'))
pd$name <- factor(pd$name, levels = unique(pd[order(pd$cluster),]$name))


# preliminary plot
plt <- ggplot(pd,
             aes(x = genotype,
                 y = median_zscore,
                 group = treatment)) +
  geom_smooth(aes(color = treatment),
              se = F,
              size = 0.5) +
  geom_point(aes(color = treatment),
             shape = 1,
             size = 1,
             alpha = 0.5) +
  labs(x = 'Genotype',
       y = 'Z-score',
       color = 'Treatment') +
  scale_color_manual(values = c('#999999', '#f1a340', 'black', 'purple')) +
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


## Figure 1-I
pdf('f1/plots/cluster.pdf',
    width = 6, height = 5)
plt
dev.off()
#######################


###### PEAK DISTRIBUTION ######
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


## Supporting Figure 2-A
pdf('f1/plots/merge_cluster_peak_distribution2_revised2.pdf',
    width = 5, height = 3)
plt
dev.off()


plot_data2 <- data.frame(peakAnno@annoStat)

plt <- ggplot(plot_data2,
             aes(x = Feature,
                 y = Frequency)) +
  geom_bar(stat = 'identity',
           position = 'dodge',
           fill = "#999999") +
  labs(x = 'Feature',
       y = 'Percentage') +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text.x = element_text(size = 10, face = 'bold', angle = 40, hjust = 1),
        axis.text.y = element_text(size = 10, face = 'bold'),
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        axis.line = element_blank())


## Supporting Figure 2-B
pdf('plots/merge_cluster_peak_distribution_ratio_revised.pdf',
    width = 10, height = 5)
plt
dev.off()
#######################


###### ASSIGN PEAK TO GENE ######
# prepare each cluster to gene file
head(peakAnno_tb)
head(peak_cluster)
peakAnno_tb <- merge(peakAnno_tb, peak_cluster)
peak2gene <- peakAnno_tb %>% 
  dplyr::select(geneId, cluster) %>% 
  distinct()
peak2gene$cluster <- factor(peak2gene$cluster, levels = paste0('C', 1:15))
peak2gene <- peak2gene[order(peak2gene$cluster),]

p <- 'counts/split'
write.table(peak2gene, paste(p, 'peak2gene.txt', sep = '/'), row.names = F, quote = F, sep = '\t')

# create folder to save results
# dir.create(paste(p, 'peak2gene', sep = '/'))

for (i in paste0('C', 1:15)) {
  gene <- peak2gene %>% 
    filter(cluster == i) %>% 
    dplyr::select(geneId)
  write.table(gene, paste(p, 'peak2gene', paste0(i,'.txt'), sep = '/'),
              row.names = F, quote = F, sep = '\t')
}


# summary peaks to genes
head(peak_cluster)
head(peak2gene)

c1 <- peak_cluster %>% 
  group_by(cluster) %>% 
  summarise(peak_count = n())

c2 <- peak2gene %>% 
  group_by(cluster) %>% 
  summarise(gene_count = n())


cc <- merge(c1, c2) %>% 
  gather(type, count, -cluster)
head(cc)


cc$cluster <- factor(cc$cluster, paste0('C', 1:15))
cc$type <- factor(cc$type, c('peak_count', 'gene_count'), c('Peaks', 'Genes'))

plt <- cc %>% 
  ggplot(.,
         aes(x = cluster,
             y = count,
             fill = type)) +
  geom_bar(stat = 'identity',
           position = 'dodge') +
  labs(x = 'Cluster',
       y = '# of DE peaks and closest genes',
       fill = 'Type') +
  theme_classic() +
  scale_fill_manual(values = c('#fc8d59', '#91bfdb')) +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text.x = element_text(size = 10, face = 'bold'),
        axis.text.y = element_text(size = 10, face = 'bold'),
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        axis.line = element_blank())


## Supporting Figure 2-C
pdf('f1/plots/peak_gene_counts.pdf',
    width = 7, height = 3)
plt
dev.off()
#######################


###### CHROMATIN STATE BY ChromHMM ######
p_file <- 'f1/epi'

data <- read.table(paste(p_file, 'histone_merge_peak_counts.txt', sep = '/'), header = T) %>%
  rowwise() %>% 
  mutate(H3_Tak1 = round(median(c(H3_Tak1, Upp5_H3), na.rm = T)),
         H3K27me3_Tak1 = round(median(c(H3K27me3_Tak1, Upp5_H3K27me3), nar.rm = T)),
         H2Aub_Tak1 = Upp5_H2Aub) %>% 
  dplyr::rename(Peak = Geneid,
                ATAC_Tak1 = ATAC.seq_Tak1) %>% 
  dplyr::select(-Upp5_H2Aub, -Upp5_H3K27me3, -Upp5_H3, -ATAC_Tak1) %>% 
  column_to_rownames(var = 'Peak')

write.table(data, paste(p_file, 'modified_histone_merge_peak_counts.txt', sep = '/'),
            row.names = T, col.names = T, quote = F, sep = '\t')


# Remove rows where any sample has reads exceeding 10,000
counts_filtered <- data[apply(data, 1, function(row) all(row <= 10000)), ]


# calculate the CPM from the counts matrix
# the following command works because 
# R calculates everything by columns
cpm <- t(t(counts_filtered)*(1000000/colSums(counts_filtered)))


# remove all tiles which do not contain reads
cpm <- cpm[rowSums(cpm) > 0,]
dim(cpm) # 40847 11
#######################


###### OBSERVED/EXPECTED ######
path <- 'f1/chromHMM/peak_histone_mark_presense_table'
dirs <- list.dirs(path = path,
                  full.names = F,
                  recursive = F)

data <- data.frame()
for (dir in dirs) {
  files <- list.files(path = paste(path, dir, sep = '/'),
                      full.names = T)
  for (i in 1:length(files)) {
    
    df <- read.delim(files[i]) %>%
      filter(Class == 'pos') %>% 
      group_by(Class) %>% 
      summarise(count = n(),
                H2Aub = sum(H2Aub),
                H2AZ = sum(H2AZ),
                H3 = sum(H3),
                H3K4me1 = sum(H3K4me1),
                H3K4me3 = sum(H3K4me3),
                H3K9ac = sum(H3K9ac),
                H3K9me1 = sum(H3K9me1),
                H3K14ac = sum(H3K14ac),
                H3K27me1 = sum(H3K27me1),
                H3K27me3 = sum(H3K27me3),
                H3K36me3 = sum(H3K36me3)) %>% 
      gather(key = 'histone', value = 'overlapping', -Class, -count) %>% 
      spread(key = 'Class', value = 'overlapping') %>% 
      relocate(histone) %>% 
      mutate(time = i,
             cluster = dir)
    
    data <- rbind(data, df)
  }
}

data <- data %>% 
  group_by(cluster, histone) %>% 
  summarise(total_pos = median(count),
            count_pos = median(pos))


path2 <- 'f1/chromHMM/negative_peak_histone_mark_intersect'
files2 <- list.files(path = path2,
                     full.names = F)
data2 <- data.frame()
for (file in files2) {
  df <- read.delim(paste(path2, file, sep = '/'))
  sub <- data.frame(histone = file,
                    total_neg = 10092,
                    count_neg = nrow(df))
  data2 <- rbind(data2, sub)
  
}
data2$histone <- str_remove(data2$histone, '_Tak1_TN')


export <- merge(data,data2) %>% 
  rowwise() %>% 
  mutate(percent_neg = (count_neg/total_neg)*100,
         percent_pos = (count_pos/total_pos)*100) %>% 
  mutate(ob_ex = percent_pos/percent_neg)


# chisq.test
ch <- data.frame()
for (i in 1:nrow(export)) {
  sub <- export[i,]
  compare <- as.matrix(rbind(sub %>% 
                               dplyr::select(total_pos, count_pos) %>% 
                               mutate(class = 'pos') %>% 
                               dplyr::rename(total = total_pos,
                                             marked = count_pos),
                             sub %>% 
                               dplyr::select(total_neg, count_neg) %>% 
                               mutate(class = 'neg') %>% 
                               dplyr::rename(total = total_neg,
                                             marked = count_neg)) %>%
                         rowwise() %>% 
                         mutate(no_mark = total-marked) %>% 
                         column_to_rownames(var = 'class') %>% 
                         dplyr::select(-total))
  
  S1 <- prop.table(compare,margin = 1) 
  S2 <- chisq.test(compare,correct = TRUE)
  pvalue <- data.frame(pvalue = S2$p.value)
  
  ch <- rbind(ch, pvalue)
}

export$chisq_test <- ch$pvalue
export$sig <- ''
export[which(export$chisq_test < 0.05),]$sig <- '*'
export[which(export$chisq_test < 0.01),]$sig <- '**'
export[which(export$chisq_test < 0.001),]$sig <- '***'
export[which(export$chisq_test < 0.0001),]$sig <- '****'

export$cluster <- factor(export$cluster, levels = paste0('C', 1:15))

# percentage-odd ratio
export$histone <- factor(export$histone,
                         levels = c('H3K4me3', 'H3K36me3', 'H3K9ac', 'H3K14ac', 'H2AZ',
                                    'H2Aub', 'H3K9me1', 'H3K4me1', 'H3K27me1', 'H3K27me3'))

require(ggpubr)
plt <- export %>%
  filter(histone != 'H3') %>% 
  ggplot(.,
         aes(x = percent_pos,
             y = ob_ex,
             color = histone,
             group = cluster,
             size = count_pos)) +
  geom_smooth(method = lm,
              se = FALSE,
              color = 'grey',
              alpha = 0.1,
              key_glyph = "point") +
  geom_point(alpha = 0.5) +
  facet_wrap(.~cluster, ncol = 8) +
  coord_flip() +
  labs(x = 'Percentage (%)',
       y = 'Observed/expected',
       fill = '# of overlapiing regions',
       color = 'Histone marks and variants',
       size = 'Overlapping regions') +
  scale_color_manual(values = c('#40004b','#762a83','#9970ab','#c2a5cf','#e7d4e8','#d9f0d3','#a6dba0','#5aae61','#1b7837','#00441b')) +
  scale_y_continuous(breaks = c(0, 1, 3, 5)) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text.x = element_text(size = 8, face = 'bold'),
        axis.text.y = element_text(size = 8, face = 'bold'),
        strip.background = element_rect(fill = alpha('#91bfdb', 0.5)),
        strip.text = element_text(size = 8, face = 'bold'))
plt <- plt + stat_cor(method = "pearson",
                    label.x = 10,
                    label.y = 2.5,
                    size = 2,
                    key_glyph = "point")


## Supporting Figure 2-D
pdf('f1/plots/chromatin_state_ob_ex_correlation.pdf',
    width = 15, height = 5)
plt
dev.off()

write.table(export, 'f1/epi/expect_observe.txt',
            row.names = F, quote = F, sep = '\t')
#######################


###### EMISSION STATES TO CLUSTER PEAKS IN ATAC-SEQ #######
# file path
p <- 'f1/counts/split'
p_so <- 'f1/chromHMM/state_cluster_ovl'
so_files <- list.files(p_so,
                       pattern = '*.txt',
                       full.names = F)

peaks <- read.delim(paste(p, 'DE_peaks_hclust.txt', sep = '/'))
head(peaks)

p_emi <- 'f1/chromHMM'
emi <- read.delim(paste(p_emi, 'Tak1_12_segments.bed', sep = '/'), header = F) %>% 
  dplyr::rename(Chr = V1,
                Start = V2,
                End = V3,
                state = V4)

total_s <- emi %>% 
  group_by(state) %>% 
  summarise(total_counts_state = n())

# counts peaks in each cluster
total <- peaks %>% 
  group_by(cluster) %>% 
  summarise(total_counts = n())


sm <- data.frame()
for (file in so_files) {
  
  # cluster name
  cluster <- str_remove(file, '_state_count.txt')
  
  # read summary files
  df <- read.table(paste(p_so, file, sep = '/'), header = F) %>% 
    dplyr::rename(counts = V1,
                  state = V2) %>% 
    filter(state != 'total')
  
  # retrieve state names
  df$state <- str_remove(df$state, 'C\\d\\d_|C\\d_')
  df$state <- str_remove(df$state, '\\.bed')
  
  # pre file
  df <- df %>% 
    mutate(cluster = cluster) %>% 
    dplyr::rename(overlapping = counts)
  
  # append df to sm
  sm <- rbind(sm, df)
}

# add total counts in cluster and state
# add ratio of peak in cluster
# add hypergeometric tests from two groups
# phyper(Overlap-1, group2, Total-group2, group1,lower.tail= FALSE)
# Total is total counts in state
Total <- sum(total_s$total_counts_state)

sm2 <- sm %>% 
  merge(., total) %>% 
  merge(., total_s) %>%
  group_by(state, cluster) %>% 
  mutate(cluster_ratio = round(as.numeric(overlapping)/as.numeric(total_counts), 2)) %>% 
  mutate(hg = phyper(overlapping-1, total_counts, Total-total_counts, total_counts_state,lower.tail= FALSE)) %>% 
  mutate(`log10(hg)` = -log10(hg)) %>% 
  mutate(mark = '')

sm2[which(sm2$`hg` < 0.05),]$mark <- '*'
# sm2[which(sm2$`hg` < 0.01),]$mark <- '**'
# sm2[which(sm2$`hg` < 0.001),]$mark <- '***'

# order of state and cluster
sm2$state <- factor(sm2$state, levels = paste0('E', 1:12))
sm2$cluster <- factor(sm2$cluster, levels = paste0('C', 1:15))


plt <- sm2 %>% 
  ggplot(.,
         aes(x = cluster,
             y = state,
             fill = cluster_ratio)) +
  geom_tile() +
  geom_text(aes(label = overlapping),
            size = 4,
            fontface = 'bold') +
  geom_text(aes(label = mark),
            size = 5,
            fontface = 'bold',
            hjust = -1.5) +
  scale_fill_gradient2(low = '#5ab4ac',
                       mid = '#f5f5f5',
                       high = '#d8b365',
                       midpoint = 0.25) +
  labs(x = 'Clusters of ATAC-seq DE peaks',
       y = 'ChromHMM States',
       fill = 'Cluster Ratio') +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 10, face = 'bold'))


# clustering
df.p <- sm2 %>% 
  dplyr::select(state, cluster, cluster_ratio) %>% 
  spread(., key = 'cluster', value = 'cluster_ratio') %>% 
  column_to_rownames(var = 'state')

df.c <- sm2 %>% 
  dplyr::select(state, cluster, overlapping) %>% 
  spread(., key = 'cluster', value = 'overlapping') %>% 
  column_to_rownames(var = 'state')

df.s <- sm2 %>% 
  dplyr::select(state, cluster, mark) %>% 
  spread(., key = 'cluster', value = 'mark') %>% 
  column_to_rownames(var = 'state')


col_fun = colorRamp2(c(0, 0.25, 0.5), c('#5ab4ac', 'white', '#d8b365'))

plt <- Heatmap(df.p,
             name = 'Cluster Ratio',
             rect_gp = gpar(col = 'white',
                            lwd = 1),
             col = col_fun,
             show_row_dend = F,
             show_column_dend = T,
             cluster_rows = F,
             cluster_columns = T,
             column_names_side = 'bottom',
             column_names_gp = gpar(fontsize = 8,
                                    fontface = 'bold',
                                    hjust = 1),
             row_names_side = 'left',
             row_title = 'Chromatin state',
             column_title = 'Genome ratio',
             row_order = paste0('E', 12:1),
             row_names_gp = gpar(fontsize = 8, fontface = 'bold'),
             row_title_gp = gpar(fontsize = 10, fontface = 'bold'),
             column_title_gp = gpar(fontsize = 10, fontface = 'bold'),
             column_title_side = 'top',
             column_names_rot = 0,
             heatmap_legend_param = list(direction = 'horizontal'),
             layer_fun = function(j, i, x, y, width, height, fill, slice_r, slice_c){
               
               grid.rect(gp = gpar(lwd = 1.5, fill = "transparent"))
               
             },
             cell_fun = function(j, i, x, y, width, height, fill) {
               if(df.c[i, j] > 0)
                 grid.text(df.c[i, j], x, y, gp = gpar(fontsize = 8, fontface = 'bold'))
               if(!is.na(df.s[i, j]))
                 grid.text(df.s[i, j], x, y, gp = gpar(fontsize = 8, fontface = 'bold'),
                           hjust = -1, vjust = -0.5)
             }
)


## Supporting Figure 2-E
pdf('f1/plots/state_ovl.pdf',
    width = 9, height = 7)  
plt
dev.off()
#######################

