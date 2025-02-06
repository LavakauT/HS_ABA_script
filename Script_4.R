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


############################################# Figure 4: A,B
############################################# Supporting Figure 7: A,B,C
###### GENIE3 ######
# setwd("~/folder_to_all_the_data")
# load package
source('f4/GENIE3.R')

# load library
require(randomForest)
require(GENIE3)
require(pgirmess)

vst <- read.delim('f3/vst_norm_counts/vst_norm_counts_hs.txt') %>% 
  dplyr::rename(geneId = X)
meta <- read.delim('f3/vst_norm_counts/vst_norm_counts_hs_metadata.txt')

# gene
# cluster genes
clu <- read.delim('f3/all_genes_hclust.txt')
deg <- read.delim('f3/All_DEGs.txt') %>% 
  mutate(expression = str_remove(group, '_.*')) %>% 
  spread(., key = 'group', value = 'expression') %>% 
  gather(., key = 'group', value = 'expression', -gene) %>% 
  mutate(group = str_remove(group, '.*(?=\\s)')) %>% 
  distinct() %>% 
  filter(!is.na(expression))

# plantTFDB database for transcription factors annotation----------------
# before mapping the TF send the all new version genes to MP.database for old version
# load the converted ID file
mr.tf <- read.delim('f4/MpTak_v6.1_TF.txt') %>%
  distinct(.keep_all = TRUE)

id.conv <- deg
id.conv <- id.conv %>% 
  dplyr::rename(ID = gene)

id.conv <- left_join(id.conv, mr.tf, by = 'ID') %>% distinct()
# combine the msg.2 with transcription annotation
attr <- id.conv
attr[is.na(attr)] <- ''

# RESTRICT THE CANDIDATE REGULATORS (according to the msg.3, we know some of the genes are TFS)
input.genes <- attr %>%
  filter(TF != '') %>%
  dplyr::select('ID') %>%
  distinct()

# out of bounds: solv: inner_join
all.genes <- attr %>% 
  dplyr::select(ID) %>% 
  distinct()

expr.matrix <- vst %>% 
  dplyr::rename(ID = geneId)

expr.matrix <- inner_join(all.genes, expr.matrix, by = 'ID')
input.genes <- inner_join(expr.matrix, input.genes, by = 'ID') %>% 
  dplyr::select('ID') %>% 
  pull()
expr.matrix <- expr.matrix %>% 
  column_to_rownames(var = 'ID') %>% 
  as.matrix()

weight.matrix.tf <- GENIE3(expr.matrix, regulators = input.genes)

link.list <- getLinkList(weight.matrix.tf, threshold = 0.01)
qua <- data.frame(quantile = quantile(link.list$weight, c(0.25, 0.5, 0.75)))


link.list <- getLinkList(weight.matrix.tf, threshold = qua[1,])
write.delim(link.list, 'f4/hclust/genie3_tf_low.txt')

link.list <- getLinkList(weight.matrix.tf, threshold = qua[2,])
head(link.list)
write.delim(link.list, 'f4/hclust/genie3_tf_medium.txt')

link.list <- getLinkList(weight.matrix.tf, threshold = qua[3,])
head(link.list)
write.delim(link.list, 'f4/hclust/genie3_tf_high.txt')


write.table(attr, 'f4/hclust/genie3_node.txt',
            row.names = F, quote = F, sep = '\t')
##################


###### ARACNe ######
matrix <- rownames_to_column(as.data.frame(expr.matrix), var = 'gene')
write.delim(matrix,
            'f4/hclust/ARACNe-AP-master/hclust/expr.txt',
            row.names = FALSE)

write.delim(input.genes, 'f4/hclust/ARACNe-AP-master/hclust/tfs.txt')
##################


###### RETRIEVE KMER ON ATAC-SEQ CLUSTER ######
require(GenomicFeatures)
require(GenomicRanges)
require(ChIPseeker)

peak <- read.delim('f1/counts/split/DE_peaks_hclust.txt')
data <- peak %>% 
  dplyr::select(peak) %>% 
  mutate(c = peak) %>% 
  separate(., col = 'c', into = c('Chr', 'Start', 'End'), sep = '\\.') %>% 
  relocate(Chr, Start, End)

## make GRanges data
Txdb_gtf <- makeTxDbFromGFF('f1/MpTak_v6.1r1.gff')
gr <- makeGRangesFromDataFrame(data, keep.extra.columns=T)

## annotate
peakAnno <- annotatePeak(gr,
                         tssRegion=c(-1500, 1500),
                         genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"),
                         TxDb=Txdb_gtf,
                         level ='gene')
## export peaks annotation
## only loss 6 peaks
peakAnno
peakAnno@anno
peakAnno_tb <- as_tibble(peakAnno@anno) %>% 
  dplyr::select(peak, geneId, annotation)
peakAnno_tb[which(peakAnno_tb$annotation != 'Distal Intergenic'),]$annotation <- 'Proximal'

mat <- data.frame()
for (i in 1:15) {
  if(file.exists(paste0('f4/atac_kmer/C', i, '/',
                        'C', i, '_top_sim_dap.txt')) == TRUE){
    kmer <- read.delim(paste0('f4/atac_kmer/C', i, '/',
                              'C', i, '_top_sim_dap.txt'))
    files <- list.files(path = paste0('f4/atac/C', i))
    
    one_hot_sub <- data.frame()
    for (file in files) {
      one_hot <- read.delim(paste0('f4/atac/C', i, '/', file)) %>%
        dplyr::rename(peak = X) %>% 
        left_join(., peakAnno_tb) %>% 
        relocate(peak, geneId, Class, annotation) %>% 
        filter(annotation == 'Proximal') # critial steps to limit the TSS binding events
      kmer_family <- kmer %>% 
        filter(kmers %in% names(one_hot[,5:ncol(one_hot)]))
      
      one_hot_sub2 <- one_hot %>% 
        gather(., key = 'kmers', value = 'presence', -peak, -geneId, -Class, -annotation) %>% 
        merge(., kmer_family %>% dplyr::select(kmers, motif)) %>% 
        group_by(geneId, motif) %>% 
        summarise(presence = max(presence))
      
      one_hot_sub <- rbind(one_hot_sub, one_hot_sub2)
      
    }
    
    one_hot_sub <- one_hot_sub %>% 
      group_by(geneId, motif) %>% 
      summarise(presence = max(presence)) %>% 
      mutate(cluster = paste0('C', i))
    
    
    mat <- rbind(mat, one_hot_sub)
  }
}

mat2 <- mat %>% 
  group_by(geneId, motif) %>% 
  summarise(presence = max(presence, na.rm = T))
table(mat2$presence)


# source to target 
# The edge table was generated by Intersection GRN from genie3_tf_high.txt (GENIE3) and network.txt (ARACNe) by Cytoscape 
edge1 <- read.csv('f4/hclust/edge_table.csv') %>% 
  separate(., col = 'shared.name', into = c('source', 'target'), sep = '\\s\\(interacts with\\)\\s') %>% 
  left_join(., mr.tf %>% 
              dplyr::select(ID, family) %>% 
              dplyr::rename(source = ID,
                            source_TF = family)) %>% 
  left_join(., mr.tf %>% 
              dplyr::select(ID, family) %>% 
              dplyr::rename(target = ID,
                            target_TF = family)) %>% 
  left_join(., mat2 %>% 
              dplyr::rename(target = geneId,
                            source_TF = motif))
edge1[which(edge1$presence == 1),]$presence <- 'source'
edge1[which(edge1$presence != 'source'),]$presence <- ''
edge1[is.na(edge1$presence),]$presence <- ''


# reverse, source to target
edge2 <- read.csv('f4/hclust/edge_table.csv') %>% 
  separate(., col = 'shared.name', into = c('target', 'source'), sep = '\\s\\(interacts with\\)\\s') %>% 
  left_join(., mr.tf %>% 
              dplyr::select(ID, family) %>% 
              dplyr::rename(source = ID,
                            source_TF = family)) %>% 
  left_join(., mr.tf %>% 
              dplyr::select(ID, family) %>% 
              dplyr::rename(target = ID,
                            target_TF = family)) %>% 
  left_join(., mat2 %>% 
              dplyr::rename(target = geneId,
                            source_TF = motif))
edge2[which(edge2$presence == 1),]$presence <- 'target'
edge2[which(edge2$presence != 'target'),]$presence <- ''
edge2[is.na(edge2$presence),]$presence <- ''


# merge two edge files
left_names <- edge2[which(edge2$presence != ''),]
edge2 <- edge1 %>% 
  filter(name %in% left_names$name)
edge2$presence <- str_replace(edge2$presence, 'source', 'source_target')
edge2[which(edge2$presence == ''),]$presence <- 'target'
edge <- rbind(edge1 %>% 
                filter(!(name %in% left_names$name)),edge2) %>% 
  distinct()
write.csv(edge, 'f4/hclust/edge_table_atac.csv',
          row.names = F)
##################


###### RETRIEVE ATAC-SEQ PEAKS ######
# all peaks
p <- 'f1/counts/split'
p_plot <- 'f1/plot'
p_csv <- paste(p, 'table', sep = '/')
p_up <-  paste(p, 'up', sep = '/')
p_down <-  paste(p, 'down', sep = '/')
p_neg <-  paste(p, 'neg', sep = '/')
range <- data.frame(dir = 'f1/counts/split',
                    treatment = c(paste0(c('Tak1', 'hsfa', 'hsfb', 'dko'), '_HS'),
                                  paste0(c('hsfa', 'hsfb', 'dko'), '_HS')),
                    base = c(rep('Tak1_CK', 4),
                             paste0(c('hsfa', 'hsfb', 'dko'), '_CK')))

sum <- data.frame()
for (i in 1:4) {
  if(i == 1){
    # up
    up <- read.table(paste(p_up, paste0(range$treatment[i], '_', range$base[i], '_u.txt'), sep = '/'), header = T) %>% 
      pull()
    sum.sub_1 <- data.frame(name = unique(up),
                            genotype = range$treatment[i],
                            ocr = 'up')
    # down
    down <- read.table(paste(p_down, paste0(range$treatment[i], '_', range$base[i], '_d.txt'), sep = '/'), header = T) %>% 
      pull()
    sum.sub_2 <- data.frame(name = unique(down),
                            genotype = range$treatment[i],
                            ocr = 'down')
    
    sum.sub <- rbind(sum.sub_1, sum.sub_2) %>% 
      dplyr::select(-genotype)
    names(sum.sub) <- c('name', range$treatment[i])
    sum <- rbind(sum, sum.sub) 
  }else{
    z <- i + 3
    # combine /Tak1 and /genotype
    # /Tak1
    # up
    up <- unique(c(read.table(paste(p_up, paste0(range$treatment[i], '_', range$base[i], '_u.txt'), sep = '/'), header = T) %>% pull(),
                   read.table(paste(p_up, paste0(range$treatment[z], '_', range$base[z], '_u.txt'), sep = '/'), header = T) %>% pull()))
    sum.sub_1 <- data.frame(name = unique(up),
                            genotype = range$treatment[i],
                            ocr = 'up')
    # down
    down <- unique(c(read.table(paste(p_down, paste0(range$treatment[i], '_', range$base[i], '_d.txt'), sep = '/'), header = T) %>% pull(),
                     read.table(paste(p_down, paste0(range$treatment[z], '_', range$base[z], '_d.txt'), sep = '/'), header = T) %>% pull()))
    sum.sub_2 <- data.frame(name = unique(down),
                            genotype = range$treatment[i],
                            ocr = 'down')
    
    sum.sub <- rbind(sum.sub_1, sum.sub_2) %>% 
      dplyr::select(-genotype)
    names(sum.sub) <- c('name', range$treatment[i])
    sum <- full_join(sum, sum.sub) 
  }
}

sum$dup <- duplicated(sum$name)
# remove duplicated row which compare to its background as neg
sum <- sum %>% 
  filter(dup == FALSE) %>% 
  dplyr::select(-dup)

Txdb_gtf <- makeTxDbFromGFF('f1/MpTak_v6.1r1.gff')


# extract Chr, Start, End, Peak name
data <- sum %>% 
  mutate(Peak = name) %>%
  dplyr::select(name, Peak) %>% 
  separate(., col = 'name', into = c('Chr', 'Start', 'End'), sep = '\\.')
names(data) <- c('Chr', 'Start', 'End', 'Peak')
head(data)

# make GRanges data
gr <- makeGRangesFromDataFrame(data, keep.extra.columns=T)

# annotate
peakAnno1 <- annotatePeak(gr,
                          tssRegion=c(-1500, 1500),
                          genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
                          TxDb=Txdb_gtf,
                          level ='gene')

# export peaks annotation
peakAnno1
peakAnno1@anno
peakAnno_tb1 <- as_tibble(peakAnno1@anno) %>% 
  mutate(`distanceToTSS (log10)` = log10((abs(distanceToTSS)/1000) + 1))

peakAnno_tb1[which(peakAnno_tb1$distanceToTSS < 0),]$`distanceToTSS (log10)` <- -peakAnno_tb1[which(peakAnno_tb1$distanceToTSS < 0),]$`distanceToTSS (log10)`

sum <- left_join(sum, peakAnno_tb1 %>% dplyr::rename(name = Peak) %>% dplyr::select(name, geneId, annotation))
sum_pro <- sum %>% filter(annotation != 'Distal Intergenic')

sum_pro$dup <- duplicated(sum_pro$geneId)
sum_pro <- sum_pro %>% 
  filter(dup == FALSE) %>% 
  relocate(geneId, Tak1_HS, hsfa_HS, hsfb_HS, dko_HS, name, annotation) %>% 
  dplyr::rename(ID = geneId,
                T_OCR = Tak1_HS,
                a_OCR = hsfa_HS,
                b_OCR = hsfb_HS,
                d_OCR = dko_HS,
                Peak = name,
                OCR_to_TSS = annotation) %>% 
  dplyr::select(-dup)


node_original <- read.delim('f4/hclust/genie3_node.txt')
node_new <- left_join(node_original, sum_pro)
node_new[is.na(node_new)] <- ''

write.table(node_new, 'f4/hclust/genie3_node_ocr.txt',
            row.names = F, quote = F, sep = '\t')
##################



###### RETRIEVE VALUES OF GENE EXPRESSION AND ACCESSIBLE######
gene_expression <- read.delim('f3/vst_norm_counts/hs_zscore.txt')
meta <- read.delim('f3/vst_norm_counts/vst_norm_counts_hs_metadata.txt')
all(colnames(gene_expression) == meta$sample)
e_zscore <- gene_expression %>% 
  rowwise() %>% 
  mutate(Tak1_HS = median(s103, s104, s105, rm.na =T),
         hsfa_HS = median(s112, s116, rm.na = T),
         hsfb_HS = median(s119, s120, rm.na = T),
         dko_HS = median(s106, s107, s108, rm.na =T)) %>% 
  dplyr::select(Tak1_HS, hsfa_HS, hsfb_HS, dko_HS)
rownames(e_zscore) <- rownames(gene_expression)
write.table(e_zscore, 'f4/hclust/expression_zscore.txt',
            quote = F, sep = '\t')


atac_oc <- read.delim('f1/counts/split/DE_peaks_zscore.txt')
meta_oc <- read.delim('f1/counts/split/meta.txt')
all(colnames(atac_oc) == meta_oc$sample)
a_zscore <- atac_oc %>% 
  rowwise() %>% 
  mutate(Tak1_HS_oc = median(s3, s4, rm.na =T),
         hsfa_HS_oc = median(s7, s8, rm.na = T),
         hsfb_HS_oc = median(s11, s12, rm.na = T),
         dko_HS_oc = median(s15, s16, rm.na =T)) %>% 
  dplyr::select(Tak1_HS_oc, hsfa_HS_oc, hsfb_HS_oc, dko_HS_oc)
rownames(a_zscore) <- rownames(atac_oc)
write.table(a_zscore, 'f4/hclust/atac_zscore.txt',
            quote = F, sep = '\t')

node_new2 <- node_new %>% 
  left_join(., e_zscore %>% 
              rownames_to_column(., var = 'ID')) %>% 
  left_join(., a_zscore %>% 
              rownames_to_column(., var = 'Peak'))
write.table(node_new2, 'f4/hclust/genie3_node_ocr2.txt',
            row.names = F, quote = F, sep = '\t')
##################


###### RETRIEVE KMER ON RNA-SEQ CLUSTER ######
require(GenomicFeatures)
require(GenomicRanges)
require(ChIPseeker)

gene <- read.delim('f3/all_genes_hclust.txt')
mr.tf <- read.delim('f4/MpTak_v6.1_TF.txt') %>%
  distinct(.keep_all = TRUE)


mat <- data.frame()
for (i in 1:8) {
  if(file.exists(paste0('f3/kmer/sim/mp_rna/C', i, '/',
                        'C', i, '_top_sim_dap.txt')) == TRUE){
    kmer <- read.delim(paste0('f3/kmer/sim/mp_rna/C', i, '/',
                              'C', i, '_top_sim_dap.txt'))
    files <- list.files(path = paste0('f3/one_hot/C', i))
    
    one_hot_sub <- data.frame()
    for (file in files) {
      one_hot <- read.delim(paste0('f3/one_hot/C', i, '/', file)) %>%
        dplyr::rename(geneId = X) %>%
        relocate(geneId, Class)
      kmer_family <- kmer %>% 
        filter(kmers %in% names(one_hot[,5:ncol(one_hot)]))
      
      one_hot_sub2 <- one_hot %>% 
        gather(., key = 'kmers', value = 'presence', -geneId, -Class) %>% 
        merge(., kmer_family %>% dplyr::select(kmers, motif)) %>% 
        group_by(geneId, motif) %>% 
        summarise(presence_rna = max(presence))
      
      one_hot_sub <- rbind(one_hot_sub, one_hot_sub2)
      
    }
    
    one_hot_sub <- one_hot_sub %>% 
      group_by(geneId, motif) %>% 
      summarise(presence_rna = max(presence_rna)) %>% 
      mutate(cluster = paste0('C', i))
    
    
    mat <- rbind(mat, one_hot_sub)
  }
}

mat2 <- mat %>% 
  group_by(geneId, motif) %>% 
  summarise(presence_rna = max(presence_rna, na.rm = T))
table(mat2$presence_rna)

edge1 <- read.csv('f4/hclust/edge_table_atac.csv') %>% 
  mutate(shared.name = name) %>% 
  separate(., col = 'name', into = c('source', 'target'), sep = '\\s\\(interacts with\\)\\s') %>%
  dplyr::select(-source_TF, -target_TF) %>% 
  left_join(., mr.tf %>% 
              dplyr::select(ID, family) %>% 
              dplyr::rename(source = ID,
                            source_TF = family)) %>% 
  left_join(., mr.tf %>% 
              dplyr::select(ID, family) %>% 
              dplyr::rename(target = ID,
                            target_TF = family)) %>% 
  left_join(., mat2 %>% 
              dplyr::rename(target = geneId,
                            source_TF = motif))
edge1[which(edge1$presence_rna == 1),]$presence_rna <- 'source'
edge1[which(edge1$presence_rna != 'source'),]$presence_rna <- ''
edge1[is.na(edge1$presence_rna),]$presence_rna <- ''

edge2 <- read.csv('f4/hclust/edge_table_atac.csv') %>% 
  mutate(shared.name = name) %>%
  separate(., col = 'name', into = c('target', 'source'), sep = '\\s\\(interacts with\\)\\s') %>%
  dplyr::select(-source_TF, -target_TF) %>% 
  left_join(., mr.tf %>% 
              dplyr::select(ID, family) %>% 
              dplyr::rename(source = ID,
                            source_TF = family)) %>% 
  left_join(., mr.tf %>% 
              dplyr::select(ID, family) %>% 
              dplyr::rename(target = ID,
                            target_TF = family)) %>% 
  left_join(., mat2 %>% 
              dplyr::rename(target = geneId,
                            source_TF = motif))
edge2[which(edge2$presence_rna == 1),]$presence_rna <- 'target'
edge2[which(edge2$presence_rna != 'target'),]$presence_rna <- ''
edge2[is.na(edge2$presence_rna),]$presence_rna <- ''


# merge two edge files
left_names <- edge2[which(edge2$presence_rna != ''),]
edge2 <- edge1 %>% 
  filter(shared.name %in% left_names$shared.name)
edge2$presence_rna <- str_replace(edge2$presence_rna, 'source', 'source_target')
edge2[which(edge2$presence_rna == ''),]$presence_rna <- 'target'
edge <- rbind(edge1 %>% 
                filter(!(shared.name %in% left_names$shared.name)),edge2) %>% 
  distinct()

# annotate edge with RNA-seq only as RNA-seq
# annotate edge without any prediction as ''
edge[which(edge$presence == '' & edge$presence_rna != ''),]$presence <- 'RNA-seq'
edge[which(edge$presence == '' & edge$presence_rna == ''),]$presence <- ''

write.csv(edge, 'f4/hclust/edge_table_atac_rna.csv',
          row.names = F)
##################


###### RETRIEVE ATAC-SEQ + DAP-SEQ OF HSFA1 ######
# Tak1
file_path <- 'f4/dap/con1/bindetect_results.txt'
file <- read.delim(file_path)
motif_name <- 'MEME-1_CBRGAWDBTTCTAGA_DAP-seq_MEME-1_CBRGAWDBTTCTAGA_DAP-seq'

# remove no bound sites
file <- file %>% 
  filter(HS_bound !=0 | CK_bound != 0,
         HS_CK_log2fc > 0,
         TFBS_name %in% motif_name)


require(GenomicFeatures)
require(ChIPseeker)
Txdb_gtf <- makeTxDbFromGFF('f1/MpTak_v6.1r1.gff')

# extract Chr, Start, End, Peak name
data <- file %>% 
  dplyr::select(peak_chr, peak_start, peak_end) %>% 
  mutate(peaks = paste(peak_chr, peak_start, peak_end, sep =  '.')) %>%
  dplyr::rename(Chr = peak_chr,
                Start = peak_start,
                End = peak_end) %>% 
  relocate(Chr, Start, End, peaks)

# make GRanges data
gr <- makeGRangesFromDataFrame(data, keep.extra.columns=T)

# annotate
peakAnno <- annotatePeak(gr,
                         tssRegion=c(-1500, 1500),
                         genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
                         TxDb=Txdb_gtf,
                         level ='gene')

# export peaks annotation
peakAnno
peakAnno@anno
peakAnno_tb <- as_tibble(peakAnno@anno) %>% 
  mutate(`distanceToTSS (log10)` = log10((abs(distanceToTSS)/1000) + 1))


file2 <- file %>% 
  mutate(peaks = paste(peak_chr, peak_start, peak_end, sep =  '.')) %>% 
  merge(., peakAnno_tb %>% 
          dplyr::select(peaks, annotation, geneId)) %>% 
  distinct()

write.table(file, 'f4/hclust/HSF3_Tak1_bind.txt',
            row.names = F, quote = F, sep = '\t')



# export edge files
file2 %>% 
  dplyr::select(geneId, HS_CK_log2fc) %>% 
  mutate(source = 'MpVg00470') %>% 
  dplyr::rename(target = geneId) %>% 
  relocate(source, target) %>% 
  write.csv(., 'f4/hclust/edge_dap.csv',
            row.names = F, quote = F)
##################


###### RETRIEVE ATAC-SEQ + DAP-SEQ OF HSFA1 GENOMEWIDE ######
# Tak1
file_path <- 'f4/dap/con2/bindetect_results.txt'
file <- read.delim(file_path)
motif_name <- 'MEME-1_CBRGAWDBTTCTAGA_DAP-seq_MEME-1_CBRGAWDBTTCTAGA_DAP-seq'

# remove no bound sites
file <- file %>% 
  filter(HS_bound !=0, #CK_bound != 0,
         #HS_CK_log2fc > 0,
         TFBS_name %in% motif_name)

Txdb_gtf <- makeTxDbFromGFF('f1/MpTak_v6.1r1.gff')

# extract Chr, Start, End, Peak name
data <- file %>% 
  dplyr::select(TFBS_chr, TFBS_start, TFBS_end) %>% 
  mutate(TFBS = paste(TFBS_chr, TFBS_start, TFBS_end, sep =  '.')) %>%
  dplyr::rename(Chr = TFBS_chr,
                Start = TFBS_start,
                End = TFBS_end) %>% 
  relocate(Chr, Start, End, TFBS)

# make GRanges data
gr <- makeGRangesFromDataFrame(data, keep.extra.columns=T)

# annotate
peakAnno <- annotatePeak(gr,
                         tssRegion=c(-1500, 1500),
                         genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
                         TxDb=Txdb_gtf,
                         level ='gene')

# export peaks annotation
peakAnno
peakAnno@anno
peakAnno_tb <- as_tibble(peakAnno@anno) %>% 
  mutate(`distanceToTSS (log10)` = log10((abs(distanceToTSS)/1000) + 1))

peakAnno_tb <- peakAnno_tb %>% 
  filter(annotation != 'Distal Intergenic')

file2 <- file %>% 
  mutate(TFBS = paste(TFBS_chr, TFBS_start, TFBS_end, sep =  '.')) %>% 
  merge(., peakAnno_tb %>% 
          dplyr::select(TFBS, annotation, geneId)) %>% 
  distinct() %>% 
  filter(annotation != 'Distal Intergenic')

mr.tf <- read.delim('f4/MpTak_v6.1_TF.txt') %>%
  distinct(.keep_all = TRUE)

file3 <- file2 %>% 
  mutate(ID = geneId) %>% 
  merge(., mr.tf) %>% 
  dplyr::select(-ID)

write.table(file, 'f4/hclust/HSF3_Tak1_bind_genomewide.txt',
            row.names = F, quote = F, sep = '\t')

# export edge files
# merge with network
# The genie3_node_ocr2_clean.csv was generated from genie3_node_ocr.txt (removal of duplicated rows) by Cytoscape
nodes <- read.csv('f4/hclust/genie3_node_ocr2_clean.csv')
edges <- read.csv('f4/hclust/edge_table_atac_rna.csv')

file2 <- file2 %>% 
  dplyr::select(geneId) %>% 
  mutate(source = 'MpVg00470') %>% 
  dplyr::rename(target = geneId) %>% 
  relocate(source, target) %>% 
  distinct() %>% 
  filter(target %in% nodes$name) %>%
  mutate(weight = 0.12,
         presence = 'DAP-seq',
         presence_rna = 'source',
         source_TF = 'HSF',
         target_TF = '')
for (gene in file3$geneId) {
  family <- file3[which(file3$geneId == gene),]$family
  file2[which(file2$target %in% file3$geneId),]$target_TF <- family
}

write.csv(file2, 'f4/hclust/edge_dap_genomewide.csv',
          row.names = F, quote = F)

# revise edges
edges <- read.csv('f4/hclust/edge_table_atac_rna.csv')
edges$presence_rna_revised <- edges$presence_rna
edges[which(edges$presence != edges$presence_rna & edges$presence != '' & edges$presence != 'RNA-seq'),]$presence_rna_revised <- edges[which(edges$presence != edges$presence_rna & edges$presence != '' & edges$presence != 'RNA-seq'),]$presence

write.csv(edges,'f4/hclust/edge_table_atac_rna.csv', row.names = F, quote = F)
##################


###### DAP-SEQ PEAKS GO ANALYSIS ######
## orthologous
ortho <- read.csv('f3/MpTak_v6.1r1.protein__v__Athaliana_447_Araport11.protein.csv')
ortho2 <- str_split(ortho$MpTak_v6.1r1.protein, ", ", n = Inf, simplify =  TRUE)
ortho3 <- cbind(ortho, ortho2)
ortho3 <- ortho3[,3:ncol(ortho3[])]
ortho4 <- gather(ortho3, key = "Mpolymorpha",
                 value = "ID", -Athaliana_447_Araport11.protein)
ortho4 <- ortho4[,-2]
ortho4 <- ortho4 %>% filter(ID != "")
rep.id <- as.data.frame(str_replace_all(ortho4$ID, '\\.\\d', ''))
ortho4 <- cbind(ortho4, rep.id)
ortho4 <- ortho4[,-2]
colnames(ortho4) <- c("Athaliana", "Mpolymorpha")


order.id <- read.csv('f4/hclust/edge_dap.csv')
order.sub.ara <- ortho4 %>% filter(Mpolymorpha %in% order.id$target)
order.sub.ara <- as.data.frame(str_split(order.sub.ara$Athaliana, ", ", n = Inf, simplify = TRUE))
order.sub.ara <- gather(order.sub.ara, key = "X", value = "geneId")
order.sub.ara <- order.sub.ara %>%
  dplyr::select(geneId) %>% 
  distinct()

write.table(order.sub.ara,
            'f4/hclust/edge_dap_ara_node.txt',
            row.names = F,
            quote = F,
            sep = '\t')


## go reduced input
read.delim('f4/hclust/dap_ara_go.txt') %>% 
  filter(q.value < 0.05) %>% 
  dplyr::select(GO.ID, q.value) %>% 
  write.table(., 'f4/hclust/dap_ara_go_input.txt',
              row.names = F, quote = F, sep = '\t')


## loading go ontology
colname <- c('GO.ID', 'Term', 'Value', 'Count', 'Group')
df <- data.frame(matrix(nrow = 0, ncol = length(colname)))
names(df) <- colname
df.ori <- read.delim('f4/hclust/dap_ara_go.txt')
df.rvi <- read.delim('f4/hclust/dap_ara_go_revi.tsv') %>% 
  filter(Uniqueness < 0.9) %>% 
  dplyr::rename(GO.ID = TermID,
         Term = Name)
df.rvi <- df.rvi[order(df.rvi$Value),][1:10,] %>% 
  dplyr::select(GO.ID, Term, Value)
df.ori <- merge(df.rvi, df.ori, by = c('GO.ID', 'Term')) %>% 
  dplyr::select(GO.ID, Term, Value, Count) %>% 
  mutate(Group = 'HSF3_DAP_footprint')
df <- rbind(df, df.ori)

# export
write.table(df, 'f4/hclust/dap_go.txt', row.names = F, quote = F, sep = '\t')

# heatmap
require(circlize)
require(ComplexHeatmap)
col_fun = colorRamp2(c(-0.1, 0, 3, 7.5), c('grey60', 'white', 'darkturquoise', 'gold'))

df %>% 
  dplyr::select(Term, Group, Value) %>%
  mutate(Value = -Value) %>% 
  spread(., key = 'Group', value = 'Value') %>% 
  column_to_rownames(var = 'Term') -> df.p
df.p[is.na(df.p)] <- -0.1

df %>% 
  dplyr::select(Term, Group, Count) %>%
  spread(., key = 'Group', value = 'Count') %>% 
  column_to_rownames(var = 'Term') -> df.c

# split row
df.c[is.na(df.c)] <- 0

plt <- Heatmap(df.p,
                name = '-log10(q.value)',
                rect_gp = gpar(col = 'white',
                               lwd = 0.5),
                col = col_fun,
                show_row_dend = F,
                show_column_dend = F,
                cluster_rows = F,
                cluster_columns = F,
                column_names_side = 'top',
                column_names_gp = gpar(fontsize = 8,
                                       fontface = 'bold',
                                       hjust = 1),
                column_names_rot = 40,
                row_title = 'Representative GO terms',
                row_names_gp = gpar(fontsize = 6, hjust = 1),
                row_title_gp = gpar(fontsize = 10, fontface = 'bold'),
                column_title_gp = gpar(fontsize = 10, fontface = 'bold'),
                column_title_side = 'top',
                heatmap_legend_param = list(direction = 'horizontal'),
                layer_fun = function(j, i, x, y, width, height, fill, slice_r, slice_c){
                  
                  grid.rect(gp = gpar(lwd = 1.5, fill = "transparent"))
                  
                },
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if(df.c[i, j] > 0)
                    grid.text(df.c[i, j], x, y, gp = gpar(fontsize = 6))
                }
)
##################


###### DOWNSTREAM ANALYSIS AFTER GRN CONSTURCTION ######
# edge/node_merged_dap_footprint_target.csv were generated from the genie3_node_ocr2_clean.csv merged with DAP-seq targets
nodes <- read.csv('f4/hclust/node_merged_dap_footprint_target.csv')
edges <- read.csv('f4/hclust/edge_merged_dap_footprint_target.csv')

dap <- read.delim('f4/dap/hsfa_targets.txt')
nodes$dap <- 'NO'
nodes[which(nodes$Matching.Attribute %in% dap$gene),]$dap <- 'YES'
table(nodes$dap)


path <- 'f3/DEG'
range <- data.frame(genotype_1 = c('Tak1', 'hsfa', 'hsfb', 'dko'),
                    condition_1 = 'heat',
                    genotype_2 = c('Tak1', 'hsfa', 'hsfb', 'dko'),
                    condition_2 = 'NHS')

for (i in 1:nrow(range)) {
  temp <- read.csv(paste0(path, '/', range[i,]$genotype_1, '_', range[i,]$condition_1, '_',
                          range[i,]$genotype_2, '_', range[i,]$condition_2, '.csv')) %>% 
    dplyr::select(X, log2FoldChange)
  names(temp) <- c('gene', range[i,]$genotype_1)
  
  if(i == 1){
    expression <- temp
  }else{
    expression <- merge(expression, temp)
  }
}

a_bind <- nodes %>% 
  filter(dap == 'YES')

# a-independent
a_ind <- nodes[which(nodes$`T` == nodes$b & nodes$`T` == nodes$a),] %>%
  distinct()

a_bind_ind <- intersect(a_ind$Matching.Attribute, a_bind$Matching.Attribute)
a_no_bind_ind <- setdiff(a_ind$Matching.Attribute, a_bind$Matching.Attribute)

# a-dependent
a_dep <- nodes[which(nodes$`T` != nodes$a & nodes$`T` == nodes$b),] %>%
  distinct()

a_bind_dep <- intersect(a_dep$Matching.Attribute, a_bind$Matching.Attribute)
a_no_bind_dep <- setdiff(a_dep$Matching.Attribute, a_bind$Matching.Attribute)


require(ggpubr)


## Supporting Figure 7-B, middle
pdf('f4/hclust/plot/a_bind_vs_no_bind_dep_up.pdf',
    width = 5, height = 3)
rbind(expression %>% 
        filter(gene %in% a_no_bind_dep) %>% 
        dplyr::select(gene, Tak1, hsfa) %>%
        filter(Tak1 >= 0) %>% 
        mutate(group = 'a_no_bind_dep'),
      expression %>% 
        filter(gene %in% a_bind_dep) %>% 
        dplyr::select(gene, Tak1, hsfa) %>%
        filter(Tak1 >= 0) %>% 
        mutate(group = 'a_bind_dep')) %>% 
  ggplot(.,
         aes(x = Tak1,
             y = hsfa,
             color = factor(group, levels = c('a_bind_dep', 'a_no_bind_dep'),
                            labels = c('HSF3-binding',
                                       'No HSF3-binding')))) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(method = "pearson",
           label.x.npc = 'middle',
           size = 3.5,
           key_glyph = "point") +
  theme_classic() +
  labs(x = 'Tak1 (RNA-seq)',
       y = 'hsfa (RNA-seq)',
       color = 'HSF3-Dependent\nsub-group') +
  theme(axis.title = element_text(size = 12, face = 'bold'))
dev.off()


## Supporting Figure 7-B, right
pdf('f4/hclust/plot/a_bind_vs_no_bind_dep_down.pdf',
    width = 5, height = 3)
rbind(expression %>% 
        filter(gene %in% a_no_bind_dep) %>% 
        dplyr::select(gene, Tak1, hsfa) %>%
        filter(Tak1 < 0) %>% 
        mutate(group = 'a_no_bind_dep'),
      expression %>% 
        filter(gene %in% a_bind_dep) %>% 
        dplyr::select(gene, Tak1, hsfa) %>%
        filter(Tak1 < 0) %>% 
        mutate(group = 'a_bind_dep')) %>% 
  ggplot(.,
         aes(x = Tak1,
             y = hsfa,
             color = factor(group, levels = c('a_bind_dep', 'a_no_bind_dep'),
                            labels = c('HSF3-binding',
                                       'No HSF3-binding')))) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(method = "pearson",
           label.x.npc = 'left',
           size = 3.5,
           key_glyph = "point") +
  theme_classic() +
  labs(x = 'Tak1 (RNA-seq)',
       y = 'hsfa (RNA-seq)',
       color = 'HSF3-Dependent\nsub-group') +
  theme(axis.title = element_text(size = 12, face = 'bold'))
dev.off()


## Supporting Figure 7-B, left
pdf('f4/hclust/plot/a_bind_vs_no_bind_ind_up.pdf',
    width = 5, height = 3)
rbind(expression %>% 
        filter(gene %in% a_no_bind_ind) %>% 
        dplyr::select(gene, Tak1, hsfa) %>%
        mutate(group = 'a_no_bind_ind'),
      expression %>% 
        filter(gene %in% a_bind_ind) %>% 
        dplyr::select(gene, Tak1, hsfa) %>%
        mutate(group = 'a_bind_ind')) %>% 
  ggplot(.,
         aes(x = Tak1,
             y = hsfa,
             color = factor(group, levels = c('a_bind_ind', 'a_no_bind_ind'),
                            labels = c('HSF3-binding',
                                       'No HSF3-binding')))) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(method = "pearson",
           label.x.npc = 'middle',
           label.y.npc = 'bottom',
           size = 3.5,
           key_glyph = "point") +
  theme_classic() +
  labs(x = 'Tak1 (RNA-seq)',
       y = 'hsfa (RNA-seq)',
       color = 'HSF3-Independent\nsub-group') +
  theme(axis.title = element_text(size = 12, face = 'bold'))
dev.off()


# regulators change among genotypes
# TF only
nodes <- read.csv('f4/hclust/node_merged_dap_footprint_target.csv')
Total <- length(nodes %>% 
                  dplyr::select(Matching.Attribute) %>% 
                  distinct() %>% 
                  pull(Matching.Attribute))

tab <- nodes %>%
  filter(TF != '') %>% 
  dplyr::select(`T`, a, b, d, Matching.Attribute) %>% 
  gather(., key = 'genotype', value = 'expression', -Matching.Attribute)
Total <- nrow(nodes %>%
                filter(TF != '') %>% 
                dplyr::select(Matching.Attribute) %>% 
                distinct())

enrich <- data.frame()
for (ge in c('a', 'b', 'd')) {
  for (ex in c('U', '', 'D')) {
    # all nodes
    bac_u <- tab %>% 
      filter(genotype == 'T',
             expression == ex) %>% 
      pull(Matching.Attribute)
    tar_u <- tab %>% 
      filter(genotype == ge,
             expression == 'U') %>%  
      pull(Matching.Attribute)
    bt_u_ovl <- intersect(bac_u, tar_u)
    u_sta <- phyper(length(bt_u_ovl)-1, length(bac_u), Total-length(bac_u), length(tar_u),lower.tail= FALSE)
    
    tar_n <- tab %>% 
      filter(genotype == ge,
             expression == '') %>%  
      pull(Matching.Attribute)
    bt_un_ovl <- intersect(bac_u, tar_n)
    n_sta <- phyper(length(bt_un_ovl)-1, length(bac_u), Total-length(bac_u), length(tar_n),lower.tail= FALSE)
    
    tar_d <- tab %>% 
      filter(genotype == ge,
             expression == 'D') %>% 
      pull(Matching.Attribute)
    bt_ud_ovl <- intersect(bac_u, tar_d)
    d_sta <- phyper(length(bt_ud_ovl)-1, length(bac_u), Total-length(bac_u), length(tar_d),lower.tail= FALSE)
    
    temp <- data.frame(genotype = rep(ge, 3),
                       Tak1_type = rep(ex, 3),
                       Tak1_number = length(bac_u),
                       type = c('U', 'N', 'D'),
                       ovl = c(length(bt_u_ovl), length(bt_un_ovl), length(bt_ud_ovl)),
                       enrichemnt = c(u_sta, n_sta, d_sta))
    enrich <- rbind(enrich, temp) 
  }
}

enrich[which(enrich$Tak1_type == ''),]$Tak1_type <- 'N'
head(enrich)
enrich$genotype <- factor(enrich$genotype, levels = c('a', 'b', 'd'))
enrich$Tak1_type <- factor(enrich$Tak1_type, levels = c('U', 'N', 'D'))
enrich$type <- factor(enrich$type, levels = c('U', 'N', 'D'))
enrich$sign <- ''
enrich[which(enrich$enrichemnt < 0.05),]$sign <- '*'
enrich[which(enrich$enrichemnt < 0.01),]$sign <- '**'
enrich[which(enrich$enrichemnt < 0.001),]$sign <- '***'
enrich$type_name <- paste0(enrich$Tak1_type, ' (', enrich$Tak1_number, ')')
enrich$type_name <- factor(enrich$type_name, levels = unique(enrich[order(enrich$Tak1_type),]$type_name))


## Supporting Figure 7-A
pdf('f4/hclust/plot/enrichment_TF.pdf',
    width = 4, height = 3)
enrich %>% 
  ggplot(.,
         aes(x = genotype,
             y = round((ovl/Tak1_number)*100, 2),
             fill = type,
             label = ifelse(enrichemnt < 0.05, sign, ''))) +
  geom_bar(stat="identity") +
  geom_text(position=position_stack(vjust=0.5)) +
  facet_grid(.~type_name) +
  scale_fill_manual(values = c('red', 'grey', 'blue')) +
  labs(y = 'percentage(%)',
       fill = 'expression') +
  theme_classic()
dev.off()
##################


