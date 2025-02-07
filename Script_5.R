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


############################################# Figure 6: D,E,F,G
############################################# Supporting Figure 7: A,B,C
###### LOAD ATAC-SEQ COUNT MATRIX ######
p <- 'f5/counts/split'
df <- read.table(paste(p, 'atac_split_peak_counts.txt', sep = '/'), header=T)


# # generate mata files
meta <- data.frame(names = names(df[,-1])) %>%
  separate(., col = names, into = c('genotype', 'treatment', 'type', 'file')) %>%
  dplyr::select(-type, -file) %>%
  mutate(replicate = c(rep(c(1,2), 8)))
meta$treatment <- str_remove(meta$treatment, '\\d')
meta$genotype <- factor(meta$genotype, levels = c('Tak', 'hsfa', 'hsfb', 'dko'),
                        labels = c('Tak1', 'hsfa', 'hsfb', 'dko'))
meta$treatment <- factor(meta$treatment, levels = c('CK', 'ABA'), labels = c('Mock', 'ABA'))

meta <- meta %>%
  mutate(sample = paste0('s', 1:nrow(.)),
         batch  = c(rep(c(1,2), 8)))
meta
write.table(meta, paste(p, 'meta.txt', sep = '/'), row.names = F, quote = F, sep = '\t')

meta <- read.table(paste(p, 'meta.txt', sep = '/'), header=T)
df <- df %>% 
  column_to_rownames(var = 'Geneid')


names(df) <- meta$sample
pander(dim(df), "Data dimensions") # 53390 total events and 16 samples
meta <- meta %>% 
  column_to_rownames(var = 'sample')

# We only want to discuss Tak1 and hsfa1
meta <- meta %>% 
  filter(genotype %in% c('Tak1', 'hsfa'))
df <- df[,rownames(meta)]
#######################


####### REMOVE MISSING PEAKS #######
# remove all the rows where not a single sample has more than 50 reads
dim(df)
df <- df[apply(df, 1, max) > 50,]
dim(df)
#######################


###### DATA NORMALIZATION ######
batch <- meta[, 4]
# batch effect adjustment
adjusted.count.data <- ComBat_seq(as.matrix(df), batch = batch, covar_mod = meta[,1:2])


dds <- DESeqDataSetFromMatrix(countData = adjusted.count.data,
                              colData = meta,
                              design = ~ genotype + treatment + genotype:treatment)
ddsMF <- dds
ddsMF$group <- factor(paste(ddsMF$genotype, ddsMF$treatment, sep = '_'))
design(ddsMF) <- ~ group
head(ddsMF)
ddsMF$group <- relevel(ddsMF$group, ref = "Tak1_Mock")
ddsMF <- DESeq(ddsMF)

vsdMF <- vst(ddsMF, blind = FALSE)
head(assay(vsdMF), 5)
normalized.count <- data.frame(assay(vsdMF)) %>% 
  rownames_to_column(., var = 'peak')
write.table(normalized.count, paste(p, 'normalized_counts.txt', sep = '/'),
            row.names = F, quote = F, sep = '\t')
#######################



###### PCA ######
plotPCA(vsdMF, intgroup = c("treatment", "genotype"))  

pcaData <- plotPCA(vsdMF, intgroup = c("treatment", "genotype"), returnData = TRUE)
percenVar <- round(100 * attr(pcaData, "percentVar"))

plt <- ggplot(pcaData, aes(PC1, PC2,
                          color = factor(genotype, levels = c('Tak1', 'hsfa', 'hsfb', 'dko')),
                          shape = treatment)) + 
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percenVar[1], "% variance")) +
  ylab(paste0("PC2: ", percenVar[2], "% variance")) +
  labs(color = 'genotype', shape = 'treatment') +
  scale_color_manual(values = c('#f6efe9', '#bdc9e1', '#67a9cf', '#02818a')) +
  theme_classic() +
  coord_fixed() +
  theme(axis.text = element_text(size = 10, face = 'bold'),
        axis.title = element_text(size = 12, face = 'bold'))


## Figure 6-D, ATAC-seq
pdf('f5/plot/atac_PCA.pdf',
    width = 4, height = 4)
plt
dev.off()


save.image(file = paste(p, 'ABA_peak.RData', sep = '/'))
#######################


p <- 'f5/counts/split'
load(paste(p, 'ABA_peak.RData', sep = '/'))
Txdb_gtf <- makeTxDbFromGFF('f1/MpTak_v6.1r1.gff')


###### ANNOTATE PEAK REGION ######
p_peak <- 'f5/narrowpeak/merge'
p_plot <- 'f5/plot'
p_path <- 'f5/summary'
range <- data.frame(group1 = paste(c('Tak1', 'hsfa', 'hsfb', 'dko'), 'CK', 'merge', sep = '-'),
                    group2 = paste(c('Tak1', 'hsfa', 'hsfb', 'dko'), 'ABA', 'merge', sep = '-'))

for (i in 1:nrow(range)) {
  file_name1 <- paste0(range[i,]$group1, '_peaks.narrowPeak')
  file_name2 <- paste0(range[i,]$group2, '_peaks.narrowPeak')
  
  data1 <- read.table(paste(p_peak, file_name1, sep = '/'), header = F)
  data2 <- read.table(paste(p_peak, file_name2, sep = '/'), header = F)
  
  # extract Chr, Start, End, Peak name
  data1 <- data1[,1:4]
  names(data1) <- c('Chr', 'Start', 'End', 'Peak')
  head(data1)
  
  data2 <- data2[,1:4]
  names(data2) <- c('Chr', 'Start', 'End', 'Peak')
  head(data2)
  
  # make GRanges data
  gr1 <- makeGRangesFromDataFrame(data1, keep.extra.columns=T)
  gr2 <- makeGRangesFromDataFrame(data2, keep.extra.columns=T)
  
  # annotate
  peakAnno1 <- annotatePeak(gr1,
                            tssRegion=c(-1500, 1500),
                            genomicAnnotationPriority = c("Exon", "Intron", "5UTR", "3UTR", "Promoter", "Downstream", "Intergenic"),
                            TxDb=Txdb_gtf,
                            level ='gene')
  peakAnno2 <- annotatePeak(gr2,
                            tssRegion=c(-1500, 1500),
                            genomicAnnotationPriority = c("Exon", "Intron", "5UTR", "3UTR", "Promoter", "Downstream", "Intergenic"),
                            TxDb=Txdb_gtf,
                            level ='gene')
  
  # export peaks annotation
  # to manipulate TSS = 0 bp, i added 1000 bp to all the distanceToTSS
  peakAnno1
  peakAnno1@anno
  peakAnno_tb1 <- as_tibble(peakAnno1@anno) %>% 
    mutate(`distanceToTSS (log10)` = log10((abs(distanceToTSS)/1000) + 1))
  
  peakAnno_tb1[which(peakAnno_tb1$distanceToTSS < 0),]$`distanceToTSS (log10)` <- -peakAnno_tb1[which(peakAnno_tb1$distanceToTSS < 0),]$`distanceToTSS (log10)`
  
  peakAnno2
  peakAnno2@anno
  peakAnno_tb2 <- as_tibble(peakAnno2@anno) %>% 
    mutate(`distanceToTSS (log10)` = log10((abs(distanceToTSS)/1000) + 1))
  
  peakAnno_tb2[which(peakAnno_tb2$distanceToTSS < 0),]$`distanceToTSS (log10)` <- -peakAnno_tb2[which(peakAnno_tb2$distanceToTSS < 0),]$`distanceToTSS (log10)`
  
  
  write.table(peakAnno_tb1, paste0(p_path, '/',
                                   range[i,]$group1, '_peak_annotation.txt'),
              row.names = F, quote = F, sep = '\t')
  write.table(peakAnno_tb2, paste0(p_path, '/',
                                   range[i,]$group2, '_peak_annotation.txt'),
              row.names = F, quote = F, sep = '\t')
  
  
  # plot1
  # combine two files
  plot_data <- rbind(peakAnno_tb1 %>% 
                       dplyr::select(`distanceToTSS (log10)`) %>% 
                       mutate(group = range[i,]$group1),
                     peakAnno_tb2 %>% 
                       dplyr::select(`distanceToTSS (log10)`) %>% 
                       mutate(group = range[i,]$group2))
  plot_data$group <- str_replace(plot_data$group, '-', ' ')
  
  
  # plot2
  plot_data2 <- rbind(data.frame(peakAnno1@annoStat) %>%
                        mutate(group = range[i,]$group1),
                      data.frame(peakAnno2@annoStat) %>%
                        mutate(group = range[i,]$group2))
  plot_data2$group <- str_replace(plot_data2$group, '-', ' ')
  
  write.table(plot_data, paste0(p_path, '/',
                                range[i,]$group1, '_vs_', range[i,]$group2, '_peak_distribution.txt'),
              row.names = F, quote = F, sep = '\t')
  write.table(plot_data2, paste0(p_path, '/',
                                 range[i,]$group1, '_vs_', range[i,]$group2, '_peak_ratio.txt'),
              row.names = F, quote = F, sep = '\t')
}
#######################


###### LOAD RNA_SEQ COUNT MATRIX ######
dir <- "f3/DEG"
coldata_aba <- read.csv('f3/coldata_ABA.csv')
count_data_aba <- read.csv('f3/counts_ABA.csv')
count_data_aba <- column_to_rownames(count_data_aba, var = 'Geneid')

coldata_aba$sample <- paste0('s', 121:144)
colnames(count_data_aba) <- coldata_aba$sample
coldata_aba <- coldata_aba[,c('sample', 'genotype', 'treatment')]
coldata_aba <- column_to_rownames(coldata_aba, var = 'sample')
coldata_aba$treatment <- str_replace_all(coldata_aba$treatment, 'CK', 'mock')


coldata_aba$batch <- 3
coldata_aba$time <- '2hr'


## combind data--------------
df2 <- count_data_aba
metadata2 <- coldata_aba


# QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(df2))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data2 <- df2[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree2 <- hclust(dist(t(data2)), method = "average")
plot(htree2)


# outlier: s126, s128, s135
### NOTE: If there are batch effects observed, correct for them before moving ahead


# exclude outlier samples
samples.to.be.excluded <- c('s126', 's128', 's135')
df2.subset <- df2[,!(colnames(df2) %in% samples.to.be.excluded)]
metadata2.subset <- metadata2[!(row.names(metadata2) %in% samples.to.be.excluded),]
#######################


###### NORMALIZATION ######
# create a deseq2 dataset
#ComBat-seq batch adjusting: no apply
require(airway)
require(sva)
require(devtools)
# devtools::install_github("zhangyuqing/sva-devel", force = TRUE)

gsg <- goodSamplesGenes(t(df2.subset))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data2 <- df2.subset[gsg$goodGenes == TRUE,]
df <- data.frame(t(data2))
indx <- data.frame(sapply(df, function(x) any(is.na(x) | is.infinite(x))))
colnames(indx) <- 'is_inf_na'
indx <- indx %>% 
  filter(is_inf_na == 'FALSE')
indx <- row.names(indx)
data2 <- data2[indx,]
count_data_mat2 <- as.matrix(data2)


# making the rownames and column names identical
all(rownames(metadata2.subset) %in% colnames(count_data_mat2))
all(rownames(metadata2.subset) == colnames(count_data_mat2))


# create dds
count_data_mat2 <- data.frame(count_data_mat2)
dds <- DESeqDataSetFromMatrix(countData = count_data_mat2,
                              colData = metadata2.subset,
                              design = ~ genotype + treatment + genotype:treatment) # not spcifying model

ddsMF <- dds
ddsMF$group <- factor(paste(ddsMF$genotype, ddsMF$treatment, sep = '_'))
design(ddsMF) <- ~ group
head(ddsMF)
ddsMF$group <- relevel(ddsMF$group, ref = "Tak1_mock")
ddsMF <- DESeq(ddsMF)


## remove all genes with counts < 15 in more than 75% of samples (21*0.75=16)
## remember to change this
## suggested by WGCNA on RNA-seq FAQ

dds75 <- ddsMF[rowSums(counts(ddsMF) >= 15) >= 16,]
nrow(dds75) # 13692/19472 genes


#transform counts data into newly DESeqDataSets
vsdMF <- vst(dds75, blind = FALSE)
head(assay(vsdMF), 5) #assay is to extract the matrux of normalized values
write.csv(data.frame(assay(vsdMF)), paste0('f3/vst_norm_counts/vst_norm_counts_aba.csv'))


pcaData <- plotPCA(vsdMF, intgroup = c("treatment", "genotype"), returnData = TRUE)
percenVar <- round(100 * attr(pcaData, "percentVar"))

g2 <- ggplot(pcaData, aes(PC1, PC2, color = genotype, shape = treatment)) + 
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percenVar[1], "% variance")) +
  ylab(paste0("PC2: ", percenVar[2], "% variance")) +
  labs(color = 'treatment', shape = 'genotype') +
  theme_classic() +
  coord_fixed()

res <- results(dds75)
res

# Tak1 and hsfa only
meta <- metadata2.subset %>% 
  filter(genotype %in% c('Tak1', 'hsfa'))
count_data_mat2_2 <- count_data_mat2[,rownames(meta)]
dds <- DESeqDataSetFromMatrix(countData = count_data_mat2_2,
                              colData = meta,
                              design = ~ genotype + treatment + genotype:treatment) # not spcifying model

ddsMF <- dds
ddsMF$group <- factor(paste(ddsMF$genotype, ddsMF$treatment, sep = '_'))
design(ddsMF) <- ~ group
head(ddsMF)
ddsMF$group <- relevel(ddsMF$group, ref = "Tak1_mock")
ddsMF <- DESeq(ddsMF)

## remove all genes with counts < 15 in more than 75% of samples (12*0.75=9)
## remember to change this
## suggested by WGCNA on RNAseq FAQ

dds75 <- ddsMF[rowSums(counts(ddsMF) >= 15) >= 9,]
nrow(dds75) # 138222/19472 genes


#transform counts data into newly DESeqDataSets
vsdMF <- vst(dds75, blind = FALSE)
head(assay(vsdMF), 5) 

pcaData <- plotPCA(vsdMF, intgroup = c("treatment", "genotype"), returnData = TRUE)
percenVar <- round(100 * attr(pcaData, "percentVar"))

g2 <- ggplot(pcaData, aes(PC1, PC2,
                          color = factor(genotype, levels = c('Tak1', 'hsfa')),
                          shape = treatment)) + 
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percenVar[1], "% variance")) +
  ylab(paste0("PC2: ", percenVar[2], "% variance")) +
  labs(color = 'genotype', shape = 'treatment') +
  scale_color_manual(values = c('#f6efe9', '#bdc9e1')) +
  theme_classic() +
  coord_fixed() +
  theme(axis.text = element_text(size = 10, face = 'bold'),
        axis.title = element_text(size = 12, face = 'bold'))


pdf('f5/plot/rna_pca.pdf',
    width = 4, height = 4)
g2
dev.off()
#######################


####### OUTPUT DEG ######
resultsNames(dds75)
range <- data.frame(group = 'group',
                    treatment = rep(paste0(c('Tak1', 'hsfa'), '_ABA'), each = 4),
                    base = paste0(c('Tak1', 'hsfa'), '_mock'))
dir <- 'f3/DEG/'


for (i in 1:nrow(range)) {
  res <- results(dds75, contrast = c(range$group[i], range$treatment[i], range$base[i])) 
  res0.05 <- results(dds75, alpha = 0.05, contrast = c(range$group[i], range$treatment[i], range$base[i]))
  summary(res)
  summary(res0.05)
  
  write.csv(as.data.frame(res0.05), file = paste0(dir, range$treatment[i], '_', range$base[i], '.csv'))
  res0.05 <- read.csv(paste0(dir, range$treatment[i], '_', range$base[i], '.csv'))
  res0.05$Up_regulated <- 'NO' 
  res0.05$Up_regulated[res0.05$log2FoldChange >= 1.0 & res0.05$padj < 0.05] <- 'YES'
  res0.05$Down_regulated <- 'NO'
  res0.05$Down_regulated[res0.05$log2FoldChange <= -1.0 & res0.05$padj < 0.05] <- 'YES'
  res0.05$Negative <- 'NO' 
  res0.05$Negative[res0.05$log2FoldChange > -0.8 & res0.05$padj > 0.05] <- 'YES'
  res0.05$Negative[res0.05$log2FoldChange < 0.8 & res0.05$padj > 0.05] <- 'YES'
  write.csv(res0.05, paste0(dir, range$treatment[i], '_', range$base[i], '.csv'), row.names = FALSE)
  
  
  # output gene list
  library(pgirmess)
  res0.05 <- read.csv(paste0(dir, range$treatment[i], '_', range$base[i], '.csv'))
  
  up <- res0.05 %>% filter(Up_regulated =='YES') %>% select(X)
  write.delim(up$X, file = paste0(dir, range$treatment[i], '_', range$base[i], '_u', '.txt'))
  
  dwn <- res0.05 %>% filter(Down_regulated =='YES') %>% select(X)
  write.delim(dwn$X, file = paste0(dir, range$treatment[i], '_', range$base[i], '_d', '.txt'))
  
  neg <- res0.05 %>% filter(Negative =='YES') %>% select(X)
  write.delim(neg$X, file = paste0(dir, range$treatment[i], '_', range$base[i], '_n', '.txt'))
  
}
# DEGs can be sent to ABA
#######################


###### ATAC-SEQ BETWEEN HS AND ABA OF ABA MARKER GENES  ######
# The bed is from: awk '$3=="gene"' MpTak_v6.1r1.gff | awk 'BEGIN {OFS="\t"} {print $1,$4,$5,$9,$6,$7}' | sed -E 's/ID=(AT.{7});.*\t\./\1\t\./' > MpTak_v6.1_whole_gene_for_deeptools.bed
# This bed script is supported by https://github.com/WangLab-CEMPS/2020STARProtocols_pipeline/tree/main/UpStreamAnalysis
# The aba_marker_genes.txt is from our lab
bed <- read.delim('f5/MpTak_v6.1_whole_gene_for_deeptools.bed', header = F)
marker <- read.delim('f5/aba_marker_genes.txt', header = F)

bed$V4 <- str_remove(bed$V4, 'ID=')
bed %>% 
  filter(V4 %in% marker$V1) %>% 
  dplyr::select(-V5) %>% 
  write.table('f5/aba_marker_genes.bed',
              row.names = F, col.names = F, quote = F, sep = '\t')

# ABA
p <- 'f5/counts/split'
range <- data.frame(dir = 'f5/counts/split/table',
                    treatment = c(paste0(c('Tak1', 'hsfa', 'hsfb', 'dko'), '_ABA'),
                                  paste0(c('hsfa', 'hsfb', 'dko'), '_ABA')),
                    base = c(rep('Tak1_Mock', 4),
                             paste0(c('hsfa', 'hsfb', 'dko'), '_Mock')))

exp <- data.frame()
for (i in 1:4) {
  if(i == 1){
    # expression
    exp_sub <- read.csv(paste(range$dir[i], paste0(range$treatment[i], '_', range$base[i], '.csv'), sep = '/'), header = T) %>% 
      dplyr::rename(peak = X)
    
    data <- data.frame(peak = exp_sub$peak) %>% 
      dplyr::select(peak) %>% 
      mutate(c = peak) %>% 
      separate(., col = 'c', into = c('Chr', 'Start', 'End'), sep = '\\.') %>% 
      relocate(Chr, Start, End)
    
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
    
    peakAnno_tb2 <- peakAnno_tb %>% 
      filter(geneId %in% marker$V1,
             annotation == 'Promoter') %>% 
      dplyr::select(peak, geneId)
    
    exp.sub <- inner_join(peakAnno_tb2, exp_sub %>% dplyr::select(peak, log2FoldChange)) %>% 
      group_by(geneId) %>% 
      summarise(peak_expr = median(log2FoldChange)) %>% 
      mutate(target = range$treatment[i]) %>% 
      separate(., target, c('genotype', 'treatment'))
    exp <- rbind(exp, exp.sub)
    
  }else{
    
    # expression
    exp_sub <- read.csv(paste(range$dir[i], paste0(range$treatment[i], '_', range$base[i], '.csv'), sep = '/'), header = T) %>% 
      dplyr::rename(peak = X)
    
    data <- data.frame(peak = exp_sub$peak) %>% 
      dplyr::select(peak) %>% 
      mutate(c = peak) %>% 
      separate(., col = 'c', into = c('Chr', 'Start', 'End'), sep = '\\.') %>% 
      relocate(Chr, Start, End)
    
    ## make GRanges data
    gr <- makeGRangesFromDataFrame(data, keep.extra.columns=T)
    
    ## annotate
    peakAnno <- annotatePeak(gr,
                             tssRegion=c(-1500, 1500),
                             genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
                             TxDb=Txdb_gtf,
                             level ='gene')
    ## export peaks annotation
    ## only loss six peaks
    peakAnno
    peakAnno@anno
    peakAnno_tb <- as_tibble(peakAnno@anno) %>% 
      mutate(`distanceToTSS (log10)` = log10((abs(distanceToTSS)/1000) + 1))
    
    peakAnno_tb2 <- peakAnno_tb %>% 
      filter(geneId %in% marker$V1,
             annotation == 'Promoter') %>% 
      dplyr::select(peak, geneId)
    
    exp.sub <- inner_join(peakAnno_tb2, exp_sub %>% dplyr::select(peak, log2FoldChange)) %>% 
      group_by(geneId) %>% 
      summarise(peak_expr = median(log2FoldChange)) %>% 
      mutate(target = range$treatment[i]) %>% 
      separate(., target, c('genotype', 'treatment'))
    exp <- rbind(exp, exp.sub)
  }
}


# HS
p <- 'f1/counts/split'
range <- data.frame(dir = 'f1/counts/split/table',
                    treatment = c(paste0(c('Tak1', 'hsfa', 'hsfb', 'dko'), '_HS'),
                                  paste0(c('hsfa', 'hsfb', 'dko'), '_HS')),
                    base = c(rep('Tak1_CK', 4),
                             paste0(c('hsfa', 'hsfb', 'dko'), '_CK')))

for (i in 1:4) {
  if(i == 1){
    # expression
    exp_sub <- read.csv(paste(range$dir[i], paste0(range$treatment[i], '_', range$base[i], '.csv'), sep = '/'), header = T) %>% 
      dplyr::rename(peak = X)
    
    data <- data.frame(peak = exp_sub$peak) %>% 
      dplyr::select(peak) %>% 
      mutate(c = peak) %>% 
      separate(., col = 'c', into = c('Chr', 'Start', 'End'), sep = '\\.') %>% 
      relocate(Chr, Start, End)
    
    ## make GRanges data
    gr <- makeGRangesFromDataFrame(data, keep.extra.columns=T)
    
    ## annotate
    peakAnno <- annotatePeak(gr,
                             tssRegion=c(-1500, 1500),
                             genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
                             TxDb=Txdb_gtf,
                             level ='gene')
    ## export peaks annotation
    ## only loss six peaks
    peakAnno
    peakAnno@anno
    peakAnno_tb <- as_tibble(peakAnno@anno) %>% 
      mutate(`distanceToTSS (log10)` = log10((abs(distanceToTSS)/1000) + 1))
    
    peakAnno_tb2 <- peakAnno_tb %>% 
      filter(geneId %in% marker$V1,
             annotation == 'Promoter') %>% 
      dplyr::select(peak, geneId)
    
    exp.sub <- inner_join(peakAnno_tb2, exp_sub %>% dplyr::select(peak, log2FoldChange)) %>% 
      group_by(geneId) %>% 
      summarise(peak_expr = median(log2FoldChange)) %>% 
      mutate(target = range$treatment[i]) %>% 
      separate(., target, c('genotype', 'treatment'))
    exp <- rbind(exp, exp.sub)
    
  }else{
    
    # expression
    exp_sub <- read.csv(paste(range$dir[i], paste0(range$treatment[i], '_', range$base[i], '.csv'), sep = '/'), header = T) %>% 
      dplyr::rename(peak = X)
    
    data <- data.frame(peak = exp_sub$peak) %>% 
      dplyr::select(peak) %>% 
      mutate(c = peak) %>% 
      separate(., col = 'c', into = c('Chr', 'Start', 'End'), sep = '\\.') %>% 
      relocate(Chr, Start, End)
    
    ## make GRanges data
    gr <- makeGRangesFromDataFrame(data, keep.extra.columns=T)
    
    ## annotate
    peakAnno <- annotatePeak(gr,
                             tssRegion=c(-1500, 1500),
                             genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
                             TxDb=Txdb_gtf,
                             level ='gene')
    ## export peaks annotation
    ## only loss six peaks
    peakAnno
    peakAnno@anno
    peakAnno_tb <- as_tibble(peakAnno@anno) %>% 
      mutate(`distanceToTSS (log10)` = log10((abs(distanceToTSS)/1000) + 1))
    
    peakAnno_tb2 <- peakAnno_tb %>% 
      filter(geneId %in% marker$V1,
             annotation == 'Promoter') %>% 
      dplyr::select(peak, geneId)
    
    exp.sub <- inner_join(peakAnno_tb2, exp_sub %>% dplyr::select(peak, log2FoldChange)) %>% 
      group_by(geneId) %>% 
      summarise(peak_expr = median(log2FoldChange)) %>% 
      mutate(target = range$treatment[i]) %>% 
      separate(., target, c('genotype', 'treatment'))
    exp <- rbind(exp, exp.sub)
  }
}


xp <- full_join(exp %>% filter(genotype == 'Tak1', treatment == 'HS') %>% dplyr::select(geneId, peak_expr) %>% dplyr::rename(HS = peak_expr),
                exp %>% filter(genotype == 'Tak1', treatment == 'ABA') %>% dplyr::select(geneId, peak_expr) %>%  dplyr::rename(ABA = peak_expr))

xp <- full_join(data.frame(geneId = marker$V1),
                xp)
xp[is.na(xp)] <- 0


## Figure 6-F, Tak1
pdf('f5/plot/atac_hs_aba_cor_tak1.pdf',
    width = 4, height = 4)
p <- xp %>% 
  ggplot(.,
         aes(x = HS,
             y = ABA)) +
  geom_point() +
  geom_hline(yintercept = 0,
             color = 'red',
             alpha = 0.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 0,
             color = 'red',
             alpha = 0.5,
             linetype = 'dashed') +
  theme_classic() +
  labs(x = 'ATAC-seq Log2FC (HS)',
       y = 'ATAC-seq Log2FC (ABA)') +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 8, face = 'bold'),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        axis.line = element_blank())
p
dev.off()


xp <- full_join(exp %>% filter(genotype == 'hsfa', treatment == 'HS') %>% dplyr::select(geneId, peak_expr) %>% dplyr::rename(HS = peak_expr),
                exp %>% filter(genotype == 'hsfa', treatment == 'ABA') %>% dplyr::select(geneId, peak_expr) %>%  dplyr::rename(ABA = peak_expr))

xp <- full_join(data.frame(geneId = marker$V1),
                xp)
xp[is.na(xp)] <- 0


## Figure 6-F, hsfa
pdf('f5/plot/atac_hs_aba_cor_hsfa.pdf',
    width = 4, height = 4)
p <- xp %>% 
  ggplot(.,
         aes(x = HS,
             y = ABA)) +
  geom_point() +
  geom_hline(yintercept = 0,
             color = 'red',
             alpha = 0.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 0,
             color = 'red',
             alpha = 0.5,
             linetype = 'dashed') +
  theme_classic() +
  labs(x = 'ATAC-seq Log2FC (HS)',
       y = 'ATAC-seq Log2FC (ABA)') +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 8, face = 'bold'),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        axis.line = element_blank())
p
dev.off()


# export
geno <- c('Tak1', 'hsfa', 'hsfb', 'dko')
expt <- data.frame()
for (gg in geno) {
  xp2 <- full_join(exp %>% filter(genotype == gg, treatment == 'ABA') %>% dplyr::select(geneId, peak_expr, genotype) %>% dplyr::rename(HS_ATAC = peak_expr),
                   exp %>% filter(genotype == gg, treatment == 'HS') %>% dplyr::select(geneId, peak_expr, genotype) %>%  dplyr::rename(ABA_ATAC = peak_expr))
  
  xp2 <- full_join(data.frame(geneId = marker$V1),
                   xp2)
  xp2[is.na(xp2)] <- 0
  xp2 <- xp2 %>% 
    mutate(genotype = gg)
  
  expt <- rbind(expt, xp2)
}

write.table(expt, 'f5/summary/atac_atac.txt',
            row.names = F, quote = F, sep = '\t')
#######################


###### RNA-SEQ EXPRESSION BETWEEN HS AND ABA OF ABA MARKER GENES ######
# ABA
p <- 'f3/ABA'
range <- data.frame(dir = 'f3/ABA',
                    treatment = c(paste0(c('Tak1', 'hsfa', 'hsfb', 'dko'), '_ABA'),
                                  paste0(c('hsfa', 'hsfb', 'dko'), '_ABA')),
                    base = c(rep('Tak1_mock', 4),
                             paste0(c('hsfa', 'hsfb', 'dko'), '_mock')))

exp_rna <- data.frame()
for (i in 1:4) {
  if(i == 1){
    # expression
    exp_sub <- read.csv(paste(range$dir[i], paste0(range$treatment[i], '_', range$base[i], '.csv'), sep = '/'), header = T) %>% 
      dplyr::rename(peak = X)
    
    exp.sub <- exp_sub %>%
      filter(peak %in% marker$V1) %>% 
      dplyr::select(peak, log2FoldChange) %>% 
      dplyr::rename(geneId = peak) %>% 
      group_by(geneId) %>% 
      summarise(rna_expr = median(log2FoldChange)) %>% 
      mutate(target = range$treatment[i]) %>% 
      separate(., target, c('genotype', 'treatment'))
    
    exp_rna <- rbind(exp_rna, exp.sub)
    
  }else{
    
    # expression
    exp_sub <- read.csv(paste(range$dir[i], paste0(range$treatment[i], '_', range$base[i], '.csv'), sep = '/'), header = T) %>% 
      dplyr::rename(peak = X)
    
    exp.sub <- exp_sub %>%
      filter(peak %in% marker$V1) %>% 
      dplyr::select(peak, log2FoldChange) %>% 
      dplyr::rename(geneId = peak) %>% 
      group_by(geneId) %>% 
      summarise(rna_expr = median(log2FoldChange)) %>% 
      mutate(target = range$treatment[i]) %>% 
      separate(., target, c('genotype', 'treatment'))
    
    exp_rna <- rbind(exp_rna, exp.sub)
  }
}

aba <- exp_rna %>% 
  dplyr::rename(ABA = rna_expr)

# HS
p <- 'f3/DEG'
range <- data.frame(dir = 'f3/DEG',
                    treatment = c(paste0(c('Tak1', 'hsfa', 'hsfb', 'dko'), '_heat'),
                                  paste0(c('hsfa', 'hsfb', 'dko'), '_heat')),
                    base = c(rep('Tak1_NHS', 4),
                             paste0(c('hsfa', 'hsfb', 'dko'), '_NHS')))

exp_rna <- data.frame()
for (i in 1:4) {
  if(i == 1){
    # expression
    exp_sub <- read.csv(paste(range$dir[i], paste0(range$treatment[i], '_', range$base[i], '.csv'), sep = '/'), header = T) %>% 
      dplyr::rename(peak = X)
    
    exp.sub <- exp_sub %>%
      filter(peak %in% marker$V1) %>% 
      dplyr::select(peak, log2FoldChange) %>% 
      dplyr::rename(geneId = peak) %>% 
      group_by(geneId) %>% 
      summarise(rna_expr = median(log2FoldChange)) %>% 
      mutate(target = range$treatment[i]) %>% 
      separate(., target, c('genotype', 'treatment'))
    
    exp_rna <- rbind(exp_rna, exp.sub)
    
  }else{
    
    # expression
    exp_sub <- read.csv(paste(range$dir[i], paste0(range$treatment[i], '_', range$base[i], '.csv'), sep = '/'), header = T) %>% 
      dplyr::rename(peak = X)
    
    exp.sub <- exp_sub %>%
      filter(peak %in% marker$V1) %>% 
      dplyr::select(peak, log2FoldChange) %>% 
      dplyr::rename(geneId = peak) %>% 
      group_by(geneId) %>% 
      summarise(rna_expr = median(log2FoldChange)) %>% 
      mutate(target = range$treatment[i]) %>% 
      separate(., target, c('genotype', 'treatment'))
    
    exp_rna <- rbind(exp_rna, exp.sub)
  }
}

exp_rna$treatment <- str_replace(exp_rna$treatment, 'heat', 'HS')

hs <- exp_rna %>% 
  dplyr::rename(HS = rna_expr)

xp2 <- full_join(aba %>% filter(genotype == 'Tak1', treatment == 'ABA') %>% dplyr::select(-treatment),
                 hs %>% filter(genotype == 'Tak1', treatment == 'HS') %>% dplyr::select(-treatment),
                 by = c('geneId', 'genotype'))

xp2 <- full_join(data.frame(geneId = marker$V1),
                 xp2)
xp2[is.na(xp2)] <- 0


## Figure 6-G, Tak1
pdf('plot/rna_hs_aba_cor_tak1.pdf',
    width = 4, height = 4)
p <- xp2 %>% 
  ggplot(.,
         aes(x = HS,
             y = ABA)) +
  geom_point() +
  geom_hline(yintercept = 0,
             color = 'red',
             alpha = 0.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 0,
             color = 'red',
             alpha = 0.5,
             linetype = 'dashed') +
  theme_classic() +
  labs(x = 'RNA-seq Log2FC (HS)',
       y = 'RNA-seq Log2FC (ABA)') +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 8, face = 'bold'),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        axis.line = element_blank())
p
dev.off()


xp2 <- full_join(aba %>% filter(genotype == 'hsfa', treatment == 'ABA') %>% dplyr::select(-treatment),
                 hs %>% filter(genotype == 'hsfa', treatment == 'HS') %>% dplyr::select(-treatment),
                 by = c('geneId', 'genotype'))

xp2 <- full_join(data.frame(geneId = marker$V1),
                 xp2)
xp2[is.na(xp2)] <- 0


## Figure 6-G, hsfa
pdf('plot/rna_hs_aba_cor_hsfa.pdf',
    width = 4, height = 4)
p <- xp2 %>% 
  ggplot(.,
         aes(x = HS,
             y = ABA)) +
  geom_point() +
  geom_hline(yintercept = 0,
             color = 'red',
             alpha = 0.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 0,
             color = 'red',
             alpha = 0.5,
             linetype = 'dashed') +
  theme_classic() +
  labs(x = 'RNA-seq Log2FC (HS)',
       y = 'RNA-seq Log2FC (ABA)') +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 8, face = 'bold'),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        axis.line = element_blank())
p
dev.off()


# export
geno <- c('Tak1', 'hsfa', 'hsfb', 'dko')
expt <- data.frame()
for (gg in geno) {
  xp2 <- full_join(aba %>% filter(genotype == gg, treatment == 'ABA') %>% dplyr::select(-treatment),
                   hs %>% filter(genotype == gg, treatment == 'HS') %>% dplyr::select(-treatment),
                   by = c('geneId', 'genotype'))
  
  xp2 <- full_join(data.frame(geneId = marker$V1),
                   xp2)
  xp2[is.na(xp2)] <- 0
  
  expt <- rbind(expt, xp2)
}

write.table(expt, 'f5/summary/rna_rna.txt',
            row.names = F, quote = F, sep = '\t')
#######################


###### RNA-SEQ VS ATAC-SEQ OF ABA MARKER GENES ######
# ABA RNA-seq
marker <- read.delim('f5/aba_marker_genes.txt', header = F)
Txdb_gtf <- makeTxDbFromGFF('f1/MpTak_v6.1r1.gff')
p <- 'f3/ABA'
range <- data.frame(dir = 'f3/ABA',
                    treatment = c(paste0(c('Tak1', 'hsfa', 'hsfb', 'dko'), '_ABA'),
                                  paste0(c('hsfa', 'hsfb', 'dko'), '_ABA')),
                    base = c(rep('Tak1_mock', 4),
                             paste0(c('hsfa', 'hsfb', 'dko'), '_mock')))

exp_rna <- data.frame()
for (i in 1:4) {
  if(i == 1){
    # expression
    exp_sub <- read.csv(paste(range$dir[i], paste0(range$treatment[i], '_', range$base[i], '.csv'), sep = '/'), header = T) %>% 
      dplyr::rename(peak = X)
    
    exp.sub <- exp_sub %>%
      filter(peak %in% marker$V1) %>% 
      dplyr::select(peak, log2FoldChange) %>% 
      dplyr::rename(geneId = peak) %>% 
      group_by(geneId) %>% 
      summarise(rna_expr = median(log2FoldChange)) %>% 
      mutate(target = range$treatment[i]) %>% 
      separate(., target, c('genotype', 'treatment'))
    
    exp_rna <- rbind(exp_rna, exp.sub)
    
  }else{
    
    # expression
    exp_sub <- read.csv(paste(range$dir[i], paste0(range$treatment[i], '_', range$base[i], '.csv'), sep = '/'), header = T) %>% 
      dplyr::rename(peak = X)
    
    exp.sub <- exp_sub %>%
      filter(peak %in% marker$V1) %>% 
      dplyr::select(peak, log2FoldChange) %>% 
      dplyr::rename(geneId = peak) %>% 
      group_by(geneId) %>% 
      summarise(rna_expr = median(log2FoldChange)) %>% 
      mutate(target = range$treatment[i]) %>% 
      separate(., target, c('genotype', 'treatment'))
    
    exp_rna <- rbind(exp_rna, exp.sub)
  }
}

rna <- exp_rna %>% 
  dplyr::rename(RNA = rna_expr)

# ABA
p <- 'f5/counts/split'
range <- data.frame(dir = 'f5/counts/split/table',
                    treatment = c(paste0(c('Tak1', 'hsfa', 'hsfb', 'dko'), '_ABA'),
                                  paste0(c('hsfa', 'hsfb', 'dko'), '_ABA')),
                    base = c(rep('Tak1_Mock', 4),
                             paste0(c('hsfa', 'hsfb', 'dko'), '_Mock')))

exp <- data.frame()
for (i in 1:4) {
  if(i == 1){
    # expression
    exp_sub <- read.csv(paste(range$dir[i], paste0(range$treatment[i], '_', range$base[i], '.csv'), sep = '/'), header = T) %>% 
      dplyr::rename(peak = X)
    
    data <- data.frame(peak = exp_sub$peak) %>% 
      dplyr::select(peak) %>% 
      mutate(c = peak) %>% 
      separate(., col = 'c', into = c('Chr', 'Start', 'End'), sep = '\\.') %>% 
      relocate(Chr, Start, End)
    
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
    
    peakAnno_tb2 <- peakAnno_tb %>% 
      filter(geneId %in% marker$V1,
             annotation == 'Promoter') %>% 
      dplyr::select(peak, geneId)
    
    exp.sub <- inner_join(peakAnno_tb2, exp_sub %>% dplyr::select(peak, log2FoldChange)) %>% 
      group_by(geneId) %>% 
      summarise(peak_expr = median(log2FoldChange)) %>% 
      mutate(target = range$treatment[i]) %>% 
      separate(., target, c('genotype', 'treatment'))
    exp <- rbind(exp, exp.sub)
    
  }else{
    
    # expression
    exp_sub <- read.csv(paste(range$dir[i], paste0(range$treatment[i], '_', range$base[i], '.csv'), sep = '/'), header = T) %>% 
      dplyr::rename(peak = X)
    
    data <- data.frame(peak = exp_sub$peak) %>% 
      dplyr::select(peak) %>% 
      mutate(c = peak) %>% 
      separate(., col = 'c', into = c('Chr', 'Start', 'End'), sep = '\\.') %>% 
      relocate(Chr, Start, End)
    
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
    
    peakAnno_tb2 <- peakAnno_tb %>% 
      filter(geneId %in% marker$V1,
             annotation == 'Promoter') %>% 
      dplyr::select(peak, geneId)
    
    exp.sub <- inner_join(peakAnno_tb2, exp_sub %>% dplyr::select(peak, log2FoldChange)) %>% 
      group_by(geneId) %>% 
      summarise(peak_expr = median(log2FoldChange)) %>% 
      mutate(target = range$treatment[i]) %>% 
      separate(., target, c('genotype', 'treatment'))
    exp <- rbind(exp, exp.sub)
  }
}

atac <- exp %>% 
  dplyr::rename(ATAC = peak_expr)

xp2 <- full_join(rna %>% filter(genotype == 'Tak1', treatment == 'ABA') %>% dplyr::select(-treatment),
                 atac %>% filter(genotype == 'Tak1', treatment == 'ABA') %>% dplyr::select(-treatment),
                 by = c('geneId', 'genotype'))

xp2 <- full_join(data.frame(geneId = marker$V1),
                 xp2)
xp2[is.na(xp2)] <- 0


## Figure 6-E, Tak1
pdf('g5/plot/atac_rna_aba_cor_tak1.pdf',
    width = 4, height = 4)
p <- xp2 %>% 
  ggplot(.,
         aes(x = RNA,
             y = ATAC)) +
  geom_point() +
  geom_hline(yintercept = 0,
             color = 'red',
             alpha = 0.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 0,
             color = 'red',
             alpha = 0.5,
             linetype = 'dashed') +
  theme_classic() +
  labs(x = 'RNA-seq Log2FC (ABA)',
       y = 'ATAC-seq Log2FC (ABA)') +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 8, face = 'bold'),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        axis.line = element_blank())
p
dev.off()


xp2 <- full_join(rna %>% filter(genotype == 'hsfa', treatment == 'ABA') %>% dplyr::select(-treatment),
                 atac %>% filter(genotype == 'hsfa', treatment == 'ABA') %>% dplyr::select(-treatment),
                 by = c('geneId', 'genotype'))

xp2 <- full_join(data.frame(geneId = marker$V1),
                 xp2)
xp2[is.na(xp2)] <- 0


## Figure 6-E, hsfa
pdf('g5/plot/atac_rna_aba_cor_hsfa.pdf',
    width = 4, height = 4)
p <- xp2 %>% 
  ggplot(.,
         aes(x = RNA,
             y = ATAC)) +
  geom_point() +
  geom_hline(yintercept = 0,
             color = 'red',
             alpha = 0.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 0,
             color = 'red',
             alpha = 0.5,
             linetype = 'dashed') +
  theme_classic() +
  labs(x = 'RNA-seq Log2FC (ABA)',
       y = 'ATAC-seq Log2FC (ABA)') +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 8, face = 'bold'),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        axis.line = element_blank())
p
dev.off()


# export
geno <- c('Tak1', 'hsfa', 'hsfb', 'dko')
expt <- data.frame()
for (gg in geno) {
  xp2 <- full_join(atac %>% filter(genotype == gg, treatment == 'ABA') %>% dplyr::select(-treatment) %>% dplyr::rename(ABA_ATAC = peak_expr),
                   rna %>% filter(genotype == gg, treatment == 'ABA') %>% dplyr::select(-treatment) %>% dplyr::rename(ABA_RNA = rna_expr),
                   by = c('geneId', 'genotype'))

  xp2 <- full_join(data.frame(geneId = marker$V1),
                   xp2)
  xp2[is.na(xp2)] <- 0

  expt <- rbind(expt, xp2)
}

write.table(expt, 'f5/summary/aba_atac_rna.txt',
            row.names = F, quote = F, sep = '\t')
#######################


###### RNA-SEQ VS ATAC-SEQ OF ABA MARKER GENES UNDER HS ######
# RNA
p <- 'f3/DEG'
range <- data.frame(dir = 'f3/DEG',
                    treatment = c(paste0(c('Tak1', 'hsfa', 'hsfb', 'dko'), '_heat'),
                                  paste0(c('hsfa', 'hsfb', 'dko'), '_heat')),
                    base = c(rep('Tak1_NHS', 4),
                             paste0(c('hsfa', 'hsfb', 'dko'), '_NHS')))

exp_rna <- data.frame()
for (i in 1:4) {
  if(i == 1){
    # expression
    exp_sub <- read.csv(paste(range$dir[i], paste0(range$treatment[i], '_', range$base[i], '.csv'), sep = '/'), header = T) %>% 
      dplyr::rename(peak = X)
    
    exp.sub <- exp_sub %>%
      filter(peak %in% marker$V1) %>% 
      dplyr::select(peak, log2FoldChange) %>% 
      dplyr::rename(geneId = peak) %>% 
      group_by(geneId) %>% 
      summarise(rna_expr = median(log2FoldChange)) %>% 
      mutate(target = range$treatment[i]) %>% 
      separate(., target, c('genotype', 'treatment'))
    
    exp_rna <- rbind(exp_rna, exp.sub)
    
  }else{
    
    # expression
    exp_sub <- read.csv(paste(range$dir[i], paste0(range$treatment[i], '_', range$base[i], '.csv'), sep = '/'), header = T) %>% 
      dplyr::rename(peak = X)
    
    exp.sub <- exp_sub %>%
      filter(peak %in% marker$V1) %>% 
      dplyr::select(peak, log2FoldChange) %>% 
      dplyr::rename(geneId = peak) %>% 
      group_by(geneId) %>% 
      summarise(rna_expr = median(log2FoldChange)) %>% 
      mutate(target = range$treatment[i]) %>% 
      separate(., target, c('genotype', 'treatment'))
    
    exp_rna <- rbind(exp_rna, exp.sub)
  }
}

exp_rna$treatment <- str_replace(exp_rna$treatment, 'heat', 'HS')

rna <- exp_rna %>% 
  dplyr::rename(RNA = rna_expr)


# ATAC
p <- 'f1/counts/split'
range <- data.frame(dir = 'f1/counts/split/table',
                    treatment = c(paste0(c('Tak1', 'hsfa', 'hsfb', 'dko'), '_HS'),
                                  paste0(c('hsfa', 'hsfb', 'dko'), '_HS')),
                    base = c(rep('Tak1_CK', 4),
                             paste0(c('hsfa', 'hsfb', 'dko'), '_CK')))

exp <- data.frame()
for (i in 1:4) {
  if(i == 1){
    # expression
    exp_sub <- read.csv(paste(range$dir[i], paste0(range$treatment[i], '_', range$base[i], '.csv'), sep = '/'), header = T) %>% 
      dplyr::rename(peak = X)
    
    data <- data.frame(peak = exp_sub$peak) %>% 
      dplyr::select(peak) %>% 
      mutate(c = peak) %>% 
      separate(., col = 'c', into = c('Chr', 'Start', 'End'), sep = '\\.') %>% 
      relocate(Chr, Start, End)
    
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
    
    peakAnno_tb2 <- peakAnno_tb %>% 
      filter(geneId %in% marker$V1,
             annotation == 'Promoter') %>% 
      dplyr::select(peak, geneId)
    
    exp.sub <- inner_join(peakAnno_tb2, exp_sub %>% dplyr::select(peak, log2FoldChange)) %>% 
      group_by(geneId) %>% 
      summarise(peak_expr = median(log2FoldChange)) %>% 
      mutate(target = range$treatment[i]) %>% 
      separate(., target, c('genotype', 'treatment'))
    exp <- rbind(exp, exp.sub)
    
  }else{
    
    # expression
    exp_sub <- read.csv(paste(range$dir[i], paste0(range$treatment[i], '_', range$base[i], '.csv'), sep = '/'), header = T) %>% 
      dplyr::rename(peak = X)
    
    data <- data.frame(peak = exp_sub$peak) %>% 
      dplyr::select(peak) %>% 
      mutate(c = peak) %>% 
      separate(., col = 'c', into = c('Chr', 'Start', 'End'), sep = '\\.') %>% 
      relocate(Chr, Start, End)
    
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
    
    peakAnno_tb2 <- peakAnno_tb %>% 
      filter(geneId %in% marker$V1,
             annotation == 'Promoter') %>% 
      dplyr::select(peak, geneId)
    
    exp.sub <- inner_join(peakAnno_tb2, exp_sub %>% dplyr::select(peak, log2FoldChange)) %>% 
      group_by(geneId) %>% 
      summarise(peak_expr = median(log2FoldChange)) %>% 
      mutate(target = range$treatment[i]) %>% 
      separate(., target, c('genotype', 'treatment'))
    exp <- rbind(exp, exp.sub)
  }
}

atac <- exp %>% 
  dplyr::rename(ATAC = peak_expr)

xp2 <- full_join(rna %>% filter(genotype == 'Tak1', treatment == 'HS') %>% dplyr::select(-treatment),
                 atac %>% filter(genotype == 'Tak1', treatment == 'HS') %>% dplyr::select(-treatment),
                 by = c('geneId', 'genotype'))

xp2 <- full_join(data.frame(geneId = marker$V1),
                 xp2)
xp2[is.na(xp2)] <- 0


pdf('g5/plot/atac_rna_hs_cor_tak1.pdf',
    width = 4, height = 4)
p <- xp2 %>% 
  ggplot(.,
         aes(x = RNA,
             y = ATAC)) +
  geom_point() +
  geom_hline(yintercept = 0,
             color = 'red',
             alpha = 0.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 0,
             color = 'red',
             alpha = 0.5,
             linetype = 'dashed') +
  theme_classic() +
  labs(x = 'RNA-seq Log2FC (ABA)',
       y = 'ATAC-seq Log2FC (ABA)') +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 8, face = 'bold'),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        axis.line = element_blank())
p
dev.off()


xp2 <- full_join(rna %>% filter(genotype == 'hsfa', treatment == 'ABA') %>% dplyr::select(-treatment),
                 atac %>% filter(genotype == 'hsfa', treatment == 'ABA') %>% dplyr::select(-treatment),
                 by = c('geneId', 'genotype'))

xp2 <- full_join(data.frame(geneId = marker$V1),
                 xp2)
xp2[is.na(xp2)] <- 0


pdf('g5/plot/atac_rna_hs_cor_hsfa.pdf',
    width = 4, height = 4)
p <- xp2 %>% 
  ggplot(.,
         aes(x = RNA,
             y = ATAC)) +
  geom_point() +
  geom_hline(yintercept = 0,
             color = 'red',
             alpha = 0.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 0,
             color = 'red',
             alpha = 0.5,
             linetype = 'dashed') +
  theme_classic() +
  labs(x = 'RNA-seq Log2FC (ABA)',
       y = 'ATAC-seq Log2FC (ABA)') +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 8, face = 'bold'),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        axis.line = element_blank())
p
dev.off()


# export
geno <- c('Tak1', 'hsfa', 'hsfb', 'dko')
expt <- data.frame()
for (gg in geno) {
  xp2 <- full_join(exp %>% filter(genotype == gg, treatment == 'HS') %>% dplyr::select(-treatment) %>% dplyr::rename(HS_ATAC = peak_expr),
                   exp_rna %>% filter(genotype == gg, treatment == 'HS') %>% dplyr::select(-treatment) %>% dplyr::rename(HS_RNA = rna_expr),
                   by = c('geneId', 'genotype'))

  xp2 <- full_join(data.frame(geneId = marker$V1),
                   xp2)
  xp2[is.na(xp2)] <- 0

  expt <- rbind(expt, xp2)
}

write.table(expt, 'f5/summary/hs_atac_rna.txt',
            row.names = F, quote = F, sep = '\t')
#######################


###### TARGET ANALYSIS SUMMARY ######
rna_rna <- read.delim('f5/summaryrna_rna.txt') %>% 
  filter(ABA > 0 & HS > 0) %>% 
  dplyr::select(geneId, genotype) %>% 
  mutate(RNA_RNA = 'V')
atac_atac <- read.delim('f5/summary/atac_atac.txt') %>% 
  filter(ABA_ATAC > 0 & HS_ATAC > 0) %>% 
  dplyr::select(geneId, genotype) %>% 
  mutate(ATAC_ATAC = 'V')
hs_atac_rna <- read.delim('f5/summary/hs_atac_rna.txt') %>% 
  filter(HS_ATAC > 0 & HS_ATAC > 0) %>% 
  dplyr::select(geneId, genotype) %>% 
  mutate(HS_ATAC_ATAC = 'V')
aba_atac_rna <- read.delim('f5/summary/aba_atac_rna.txt') %>% 
  filter(ABA_ATAC > 0 & ABA_ATAC > 0) %>% 
  dplyr::select(geneId, genotype) %>% 
  mutate(ABA_ATAC_ATAC = 'V')

data <-full_join(rna_rna, atac_atac, by = c('geneId', 'genotype')) %>% 
  full_join(., hs_atac_rna, by = c('geneId', 'genotype')) %>% 
  full_join(., aba_atac_rna, by = c('geneId', 'genotype'))


write.table(data, 'f5/summary/summary_table_52_aba_marker.txt',
            row.names = F, quote = F, sep = '\t')s
##############################
