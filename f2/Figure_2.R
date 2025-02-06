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


############################################# Figure 2: A,B,C,D,E
###### LOAD RNA-SEQ DATASET 1 ######
count_data <- read.csv("f3/raw_reads_counts_hs.csv")

rname <- count_data[,1]
count_data <- count_data[,2:ncol(count_data)]
row.names(count_data) <- rname

coldata <- read.csv('f3/condition_hs.csv')
head(coldata)

# replace the sample names of count data by the column 1 in coldata
colname <- coldata[,1]
colnames(count_data) <- colname
head(count_data)

rname2 <- coldata[,1]
coldata <- coldata[,2:ncol(coldata)]
row.names(coldata) <- rname2
coldata$genotype <- str_replace_all(coldata$genotype, 'tak', 'Tak1')
coldata$treatment <- str_replace_all(coldata$treatment, 'hs', 'heat')
coldata$treatment <- str_replace_all(coldata$treatment, 'cntr', 'NHS')
head(coldata)
##################


###### LOAD RNA-SEQ DATASET 2 ######
count_data_hsf <- read.csv('f3/GSE178776_1_raw_reads_new_version.csv')
coldata_hsf <- read.csv('f3/sample_information.csv')

select <- c('tak1', 'tak1_yko', 'tak1_bko')
coldata_hsf <- coldata_hsf %>%
  filter(genotype %in% select)

select.id <- coldata_hsf[,1]
count_data_hsf <- count_data_hsf %>% 
  column_to_rownames(var = 'Geneid') %>% 
  dplyr::select(all_of(select.id))

coldata$treatment <- gsub('heat', 'hs', coldata$treatment)
coldata_hsf$treatment <- str_replace_all(coldata_hsf$treatment, 'hs', 'heat')
coldata_hsf$genotype <- gsub('tak1_yko', 'hsfa', coldata_hsf$genotype)
coldata_hsf$genotype <- gsub('tak1_bko', 'hsfb', coldata_hsf$genotype)
coldata_hsf$samples <- paste0('s', 109:120)
colnames(count_data_hsf) <- coldata_hsf$samples
coldata_hsf$treatment <- str_replace_all(coldata_hsf$treatment, 'cntr', 'NHS')
coldata_hsf$genotype <- str_replace_all(coldata_hsf$genotype, 'tak1', 'Tak1')
head(coldata_hsf)
##################


###### DATA PROCESSING ######
coldata <- coldata[,1:2]
treatment <- as.data.frame(gsub(pattern = "1h", coldata[,2], replacement = "heat"))
row.names(treatment) <- row.names(coldata)
colnames(treatment) <- "treatment"
coldata <- cbind(coldata, treatment)
coldata <- coldata[,-2]
coldata$batch <- '1' # adding batch


coldata_hsf <- coldata_hsf %>% 
  column_to_rownames(var = 'samples') %>% 
  dplyr::select('genotype', 'treatment')
coldata_hsf$batch <- '2'


coldata$time <- '0hr'
coldata$time[coldata$treatment == 'heat'] <- '1hr'


coldata_hsf$time <- '0hr'
coldata_hsf$time[coldata_hsf$treatment == 'heat'] <- '2hr' 


coldata.hs <- rbind(coldata, coldata_hsf)

id <- row.names(count_data)
count_data_hsf <- count_data_hsf[ id, ]
count.data.hs <- cbind(count_data, count_data_hsf)
head(count.data.hs)
##################


###### COMBINE DATA ######
df <- count.data.hs
metadata <- coldata.hs


# QC - outlier detection
gsg <- goodSamplesGenes(t(df))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detected as outliers
data <- df[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
sample.name <- paste(metadata$treatment, metadata$batch, sep = '_')
sample.name2 <- colnames(data)
colnames(data) <- sample.name
htree <- hclust(dist(t(data)), method = "average")
plot(htree)

colnames(data) <- sample.name2
# no sample seems to be outlier in two dataset


# pca - method 2
# filter out inf row with bellow function
df2 <- data.frame(t(data))
indx <- data.frame(sapply(df2, function(x) any(is.na(x) | is.infinite(x))))
colnames(indx) <- 'is_inf_na'
indx <- indx %>% 
  filter(is_inf_na == 'FALSE')
indx <- row.names(indx)
data <- data[indx,]

pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

### Warning: If there are batch effects observed, correct for them before moving ahead
# exclude outlier samples
# samples.to.be.excluded <- c('GSM4615000', 'GSM4614993', 'GSM4614995')
# data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]
############################


###### NORMALIZATION ######
# exclude outlier samples
# colData <- phenoData %>% 
# filter(!row.names(.) %in% samples.to.be.excluded)

# create a deseq2 dataset

# ComBat-seq batch adjustment
# devtools::install_github("zhangyuqing/sva-devel", force = TRUE)
require(airway)
require(sva)
require(devtools)

batch <- metadata[, 3]

gsg <- goodSamplesGenes(t(df))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detected as outliers
data <- df[gsg$goodGenes == TRUE,]
df <- data.frame(t(data))
indx <- data.frame(sapply(df, function(x) any(is.na(x) | is.infinite(x))))
colnames(indx) <- 'is_inf_na'
indx <- indx %>% 
  filter(is_inf_na == 'FALSE')
indx <- row.names(indx)
data <- data[indx,]
count_data_mat <- as.matrix(data)

# set group in genotype
adjusted.count.data <- ComBat_seq(count_data_mat, batch = batch, covar_mod = metadata[,1:2])

# making the row names and column names identical
all(rownames(metadata) %in% colnames(adjusted.count.data))
all(rownames(metadata) == colnames(adjusted.count.data))

# create dds
adjusted.count.data <- data.frame(adjusted.count.data)
dds <- DESeqDataSetFromMatrix(countData = adjusted.count.data,
                              colData = metadata,
                              design = ~ genotype + treatment + genotype:treatment) # not specifying mode

ddsMF <- dds
ddsMF$group <- factor(paste(ddsMF$genotype, ddsMF$treatment, sep = '_'))
design(ddsMF) <- ~ group
head(ddsMF)
ddsMF$group <- relevel(ddsMF$group, ref = "Tak1_NHS")
ddsMF <- DESeq(ddsMF)


## remove all genes with counts < 15 in more than 75% of samples (24*0.75=18)
## remember to change this
## suggested by WGCNA on RNAseq FAQ
dds75 <- ddsMF[rowSums(counts(ddsMF) >= 15) >= 18,]
nrow(dds75) # 11878/19472 genes

#transform counts data into newly DESeqDataSets
vsdMF <- vst(dds75, blind = FALSE)
head(assay(vsdMF), 5) #assay is to extract the matrux of normalized values

dir.create("vst_norm_counts")
write.csv(data.frame(assay(vsdMF)), "vst_norm_counts/vst_norm_counts_hs.txt")


# pcaData <- plotPCA(vsdMF, intgroup = c("treatment", "genotype"), returnData = TRUE)
# percenVar <- round(100 * attr(pcaData, "percentVar"))
# plt <- ggplot(pcaData, aes(PC1, PC2, color = genotype, shape = treatment)) + 
#   geom_point(size = 3) +
#   xlab(paste0("PC1: ", percenVar[1], "% variance")) +
#   ylab(paste0("PC2: ", percenVar[2], "% variance")) +
#   labs(color = 'treatment', shape = 'genotype') +
#   theme_classic() +
#   coord_fixed()
# plt

res <- results(dds75)
res
########################




###### output HS_NHS ######
resultsNames(dds75)
range <- data.frame(group = 'group',
                    treatment = rep(paste0(c('Tak1', 'hsfa', 'hsfb', 'dko'), '_heat'), each = 4),
                    base = paste0(c('Tak1', 'hsfa', 'hsfb', 'dko'), '_NHS'))

dir.create("f3/DEG")
dir <- 'f3/DEG'

for (i in 1:nrow(range)) {
  res <- results(dds75, contrast = c(range$group[i], range$treatment[i], range$base[i])) 
  res0.05 <- results(dds75, alpha = 0.05, contrast = c(range$group[i], range$treatment[i], range$base[i]))
  summary(res)
  summary(res0.05)
  
  write.csv(as.data.frame(res0.05), file = paste0(dir, '/', range$treatment[i], '_', range$base[i], '.csv'))
  res0.05 <- read.csv(paste0(dir, '/', range$treatment[i], '_', range$base[i], '.csv'))
  res0.05$Up_regulated <- 'NO' 
  res0.05$Up_regulated[res0.05$log2FoldChange >= 1.0 & res0.05$padj < 0.05] <- 'YES'
  res0.05$Down_regulated <- 'NO'
  res0.05$Down_regulated[res0.05$log2FoldChange <= -1.0 & res0.05$padj < 0.05] <- 'YES'
  res0.05$Negative <- 'NO' 
  res0.05$Negative[res0.05$log2FoldChange > -0.8 & res0.05$padj > 0.05] <- 'YES'
  res0.05$Negative[res0.05$log2FoldChange < 0.8 & res0.05$padj > 0.05] <- 'YES'
  write.csv(res0.05, paste0(dir, '/', range$treatment[i], '_', range$base[i], '.csv'), row.names = FALSE)
  
  
  # export gene list
  library(pgirmess)
  res0.05 <- read.csv(paste0(dir, '/', range$treatment[i], '_', range$base[i], '.csv'))

  up <- res0.05 %>% filter(Up_regulated =='YES') %>% dplyr::select(X)
  write.delim(up$X, file = paste0(dir, '/', range$treatment[i], '_', range$base[i], '_u', '.txt'))

  dwn <- res0.05 %>% filter(Down_regulated =='YES') %>% dplyr::select(X)
  write.delim(dwn$X, file = paste0(dir, '/', range$treatment[i], '_', range$base[i], '_d', '.txt'))

  neg <- res0.05 %>% filter(Negative =='YES') %>% dplyr::select(X)
  write.delim(neg$X, file = paste0(dir, '/', range$treatment[i], '_', range$base[i], '_n', '.txt'))
  
}
##########################


###### ASSIGNED GENES EXPRESSION FROM ATAC-SEQ PEAKS ######
p_vst <- 'f3/vst_norm_counts'
peak2gene <- read.delim('f1/counts/split/peak2gene.txt')


# plot: scatter plot
meta <- read.delim(paste(p_vst, 'vst_norm_counts_hs_metadata.txt', sep = '/')) %>%
  mutate(replicate = c(rep(c(1,2,3), 4), rep(c(4,1), 2), rep(c(5,2), 2), c(1,2,1,2))) %>% 
  mutate(name = paste(genotype, treatment, replicate, sep = '_'))
vst_counts <- read.delim(paste(p_vst, 'vst_norm_counts_hs.txt', sep = '/')) %>%
  filter(X %in% unique(peak2gene$geneId)) %>% 
  column_to_rownames(var = 'X')

# Create a matrix
hclust_matrix <- vst_counts %>% 
  as.matrix()

# z-score
sd <- apply(hclust_matrix, 1, sd, na.rm=TRUE)
hclust_matrix <- (hclust_matrix-rowMeans(hclust_matrix))/sd

pd <- data.frame(hclust_matrix)

all(names(pd) == meta$sample) # if True, then pass

colnames(pd) <- meta$name
head(pd)

pd2 <- data.frame()
for (i in unique(peak2gene$cluster)) {
  chosen <- peak2gene %>% 
    filter(cluster == i) %>% 
    dplyr::select(geneId) %>% 
    distinct() %>% 
    pull()
  
  sub <- pd %>% 
    rownames_to_column(var = 'gene') %>% 
    filter(gene %in% chosen) %>%
    mutate(cluster = i) %>% 
    relocate(gene, cluster) %>% 
    gather(condition, zscore, -gene, -cluster) %>% 
    separate(., col = 'condition', into = c('genotype', 'treatment', 'replicate'), sep = '_') %>% 
    group_by(cluster, genotype, treatment, replicate) %>% 
    summarise(median_zscore = median(zscore))
  
  # append
  pd2 <- rbind(pd2, sub)
}



pd2$cluster <- factor(pd2$cluster, levels = paste0('C', 1:15))
pd2$genotype <- factor(pd2$genotype, levels = c('Tak1', 'hsfa', 'hsfb', 'dko'), labels = c('Tak-1', 'hsfa', 'hsfb', 'dko'))

cluster <- data.frame(table(peak2gene$cluster)) %>% 
  dplyr::rename(cluster = Var1,
                count = Freq)
pd2 <- merge(pd2, cluster) %>% 
  mutate(name = paste0(cluster, ' (n=', count, ')'))
pd2$name <- factor(pd2$name, levels = unique(pd2[order(pd2$cluster),]$name))
pd2$treatment <- factor(pd2$treatment, levels = c('NHS', 'heat'),
                        labels = c('Cntr', 'HS'))


# preliminary plot
plt <- ggplot(pd2,
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
        strip.background = element_rect(fill = alpha('#756bb1', 0.5)),
        strip.text = element_text(size = 8, face = 'bold'),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = 'black'))


## Figure 2-A
pdf('f3/plots/cluster_rna.pdf',
    width = 6, height = 5)
plt
dev.off()
##########################


###### GO terms ######
Txdb_gtf <- makeTxDbFromGFF('f1/MpTak_v6.1r1.gff')
peak_cluster <- read.table('f1/counts/split/DE_peaks_hclust.txt', header = T)
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

# GO terms
## orthologous
ortho <- read.csv('MpTak_v6.1r1.protein__v__Athaliana_447_Araport11.protein.csv')
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


dir <- 'f1/counts/split/peak2gene'
dir2 <- 'f1/counts/split/peak2gene_ara'
for (i in paste0('C', 1:15)) {
  order.id <- read.table(paste(dir, paste0(i, '.txt'), sep = '/'), header = T)
  order.sub.ara <- ortho4 %>% filter(Mpolymorpha %in% order.id$geneId)
  order.sub.ara <- as.data.frame(str_split(order.sub.ara$Athaliana, ", ", n = Inf, simplify = TRUE))
  order.sub.ara <- gather(order.sub.ara, key = "X", value = "geneId")
  order.sub.ara <- order.sub.ara %>%
    dplyr::select(geneId) %>% 
    distinct()
  
  write.table(order.sub.ara,
              paste(dir2, paste0(i, '_ara.txt'), sep = '/'),
              row.names = F,
              quote = F,
              sep = '\t')
}


## go reduced input
dir <- 'f1/counts/split/peak2gene_go_input'
dir2 <- 'f1/counts/split/peak2gene_go'
range <- list.files(dir2,
                    pattern = '*go.txt',
                    full.names = F)
for (i in seq_along(range)) {
  read.delim(paste(dir2, range[i], sep = '/')) %>% 
    filter(q.value < 0.05) %>% 
    dplyr::select(GO.ID, q.value) %>% 
    write.table(., paste(dir, range[i], sep = '/'),
                row.names = F, quote = F, sep = '\t')
}


## loading go ontology
range <- paste0('C', 1:15)
dir <- 'f1/counts/split/peak2gene_go'
dir2 <- 'f1/counts/split/peak2gene_go_revi'

colname <- c('GO.ID', 'Term', 'Value', 'Count', 'Group')
df <- data.frame(matrix(nrow = 0, ncol = length(colname)))
names(df) <- colname
for (i in seq_along(range)) {
  clu <- range[i]
  df.ori <- read.delim(paste(dir, paste(clu, 'ara_go.txt', sep = '_'), sep = '/'))
  
  file <- paste(dir2, paste(clu, 'ara_go_revi.tsv', sep = '_'), sep = '/')
  if (!file.exists(file)){
    next
  }
  df.rvi <- read.delim(paste(dir2, paste(clu, 'ara_go_revi.tsv', sep = '_'), sep = '/')) %>% 
    filter(Uniqueness < 0.9) %>% 
    rename(GO.ID = TermID,
           Term = Name)
  df.rvi <- df.rvi[order(df.rvi$Value),][1:10,] %>% 
    dplyr::select(GO.ID, Term, Value)
  df.ori <- merge(df.rvi, df.ori, by = c('GO.ID', 'Term')) %>% 
    dplyr::select(GO.ID, Term, Value, Count) %>% 
    mutate(Group = clu)
  
  df <- rbind(df, df.ori)
}

# export
write.table(df, 'f1/counts/split/go.txt', row.names = F, quote = F, sep = '\t')

# heatmap
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
colname <- c('term', 'group')
da <- data.frame(matrix(nrow = 0, ncol = length(colname)))
names(da) <- colname
for (i in seq_along(rownames(df.c))) {
  name <- rownames(df.c[i,])
  df.c[i,] %>%
    t() %>%
    data.frame() %>%
    filter_all(., any_vars(. != 0)) -> df.c.c
  if(nrow(df.c.c) > 1){
    df.c.c <- 'multiple Clusters'
    da.sub <- data.frame(term = name,
                         group = df.c.c)
    da <- rbind(da, da.sub)
  }else{
    df.c.c = rownames(df.c.c)
    da.sub <- data.frame(term = name,
                         group = df.c.c)
    da <- rbind(da, da.sub)
  }
  
}

da$group <- factor(da$group, levels = c(paste0('C', 1:15), 'multiple Clusters'))

plt <- Heatmap(df.p,
               name = '-log10(q.value)',
               rect_gp = gpar(col = 'white',
                              lwd = 0.5),
               col = col_fun,
               show_row_dend = F,
               show_column_dend = F,
               cluster_rows = F,
               cluster_columns = F,
               column_names_side = 'bottom',
               column_names_gp = gpar(fontsize = 8,
                                      fontface = 'bold',
                                      hjust = 1),
               column_names_rot = 50,
               row_title = 'Representative GO terms',
               column_title = 'Clusters',
               clustering_distance_rows = 'pearson',
               row_names_gp = gpar(fontsize = 6, hjust = 1),
               row_title_gp = gpar(fontsize = 10, fontface = 'bold'),
               column_title_gp = gpar(fontsize = 10, fontface = 'bold'),
               column_title_side = 'top',
               row_split = da$group,
               column_order = paste0('C', c(1:11, 13:15)),
               heatmap_legend_param = list(direction = 'horizontal'),
               layer_fun = function(j, i, x, y, width, height, fill, slice_r, slice_c){
                 
                 grid.rect(gp = gpar(lwd = 1.5, fill = "transparent"))
                 
               },
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(df.c[i, j] > 0)
                   grid.text(df.c[i, j], x, y, gp = gpar(fontsize = 6))
               }
)


# revised version
df2 <- read.delim('f1/counts/split/go_revised.txt')
df2.p <- df2 %>% 
  dplyr::select(Term, Group, Value) %>%
  mutate(Value = -Value) %>% 
  spread(., key = 'Group', value = 'Value') %>% 
  column_to_rownames(var = 'Term')
df2.p[is.na(df2.p)] <- -0.1

df2.c <- df2 %>% 
  dplyr::select(Term, Group, Count) %>%
  spread(., key = 'Group', value = 'Count') %>% 
  column_to_rownames(var = 'Term')

# split row
df2.c[is.na(df2.c)] <- 0
colname2 <- c('term', 'group')
da2 <- data.frame(matrix(nrow = 0, ncol = length(colname2)))
names(da2) <- colname2
for (i in seq_along(rownames(df2.c))) {
  name <- rownames(df2.c[i,])
  df2.c[i,] %>%
    t() %>%
    data.frame() %>%
    filter_all(., any_vars(. != 0)) -> df2.c.c
  if(nrow(df2.c.c) > 1){
    df2.c.c <- 'multiple Clusters'
    da2.sub <- data.frame(term = name,
                          group = df2.c.c)
    da2 <- rbind(da2, da2.sub)
  }else{
    df2.c.c = rownames(df2.c.c)
    da2.sub <- data.frame(term = name,
                          group = df2.c.c)
    da2 <- rbind(da2, da2.sub)
  }
  
}

da2$group <- factor(da2$group, levels = c(paste0('C', 1:15), 'multiple Clusters'))

col_fun = colorRamp2(c(-0.1, 0, 3, 7.5), c('grey', 'white', 'darkturquoise', 'gold'))
plt <- Heatmap(df2.p,
                 name = '-log10(q.value)',
                 rect_gp = gpar(col = 'white',
                                lwd = 0.5),
                 col = col_fun,
                 show_row_dend = F,
                 show_column_dend = F,
                 cluster_rows = F,
                 cluster_columns = F,
                 column_names_side = 'bottom',
                 column_names_gp = gpar(fontsize = 8,
                                        fontface = 'bold',
                                        hjust = 1),
                 column_names_rot = 40,
                 row_title = 'Representative GO terms',
                 column_title = 'Clusters',
                 clustering_distance_rows = 'pearson',
                 row_names_gp = gpar(fontsize = 8, fontface = 'bold', hjust = 1),
                 row_title_gp = gpar(fontsize = 10, fontface = 'bold'),
                 column_title_gp = gpar(fontsize = 10, fontface = 'bold'),
                 column_title_side = 'top',
                 row_split = da2$group,
                 column_order = paste0('C', c(1:11, 13:15)),
                 heatmap_legend_param = list(direction = 'horizontal'),
                 layer_fun = function(j, i, x, y, width, height, fill, slice_r, slice_c){
                   
                   grid.rect(gp = gpar(lwd = 1.5, fill = "transparent"))
                   
                 },
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   if(df2.c[i, j] > 0)
                     grid.text(df2.c[i, j], x, y, gp = gpar(fontsize = 8, fontface = 'bold'))
                 }
)


## Figure 2-B
pdf('f1/plots/go.pdf',
    width = 5, height = 5)
draw(plt, heatmap_legend_side = 'bottom')
dev.off()
##########################



###### Hierarchical clustering ######
# all  DEGs
p_file <- 'f3'
range <- rbind(data.frame(dir = rep('DEG',
                                    each = 4),
                          group1 = c('Tak1_heat', 'hsfa_heat',
                                     'hsfb_heat', 'dko_heat'),
                          group2 = rep(c('Tak1_NHS'), each = 4),
                          expr = rep(c('u.txt', 'd.txt'), each = 4),
                          ex = rep(c('U', 'D'), each = 4),
                          label = paste(rep(c('HS'), each = 4),
                                        c('T', 'a', 'b', 'd'))),
               data.frame(dir = rep('DEG', 6),
                          group1 = rep(c('hsfa_heat', 'hsfb_heat', 'dko_heat'), each = 2),
                          group2 = rep(c('hsfa_NHS', 'hsfb_NHS', 'dko_NHS'), each = 2),
                          expr = rep(c('u.txt', 'd.txt'), 3),
                          ex = rep(c('U', 'D'), 3),
                          label = paste(rep(c('HS'), each = 6), rep(c('a', 'b', 'd'), each = 2))))

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

df <- distinct(df) # 5947 (genes + condition)
length(unique(df$gene)) # 3166 genes 

write.table(df, paste(p_file, 'All_DEGs.txt', sep = '/'),
            row.names = FALSE, quote = FALSE , sep = '\t')


# export DEG
# deg_count <- df %>% 
#   group_by(group) %>% 
#   summarise(count = n()) %>% 
#   separate(group, c('expression', 'genotype'), sep = '\\s') %>% 
#   separate(expression, c('expression', 'treatment'), sep = '_')
# 
# p <- deg_count %>% 
#   ggplot(.,
#          aes(x = genotype,
#              y = count,
#              fill = expression)) +
#   geom_bar(stat = 'identity',
#            position = 'dodge') +
#   geom_text(aes(label = ifelse(expression == 'U', count, '')), vjust = -1, hjust = 0) +
#   geom_text(aes(label = ifelse(expression == 'D', count, '')), vjust = -1, hjust = 1) +
#   scale_fill_manual(values = c('#ca0025', '#91bfdb')) +
#   labs(x = 'Genotype (HS)',
#        y = 'counts',
#        fill = 'Expression') +
#   theme_classic() +
#   theme(axis.title = element_text(size = 12, face = 'bold'),
#         axis.text.x = element_text(size = 10, face = 'bold'),
#         axis.text.y = element_text(size = 10, face = 'bold'))

p_vst <- 'vst_norm_counts'
meta <- read.delim(paste(p_vst, 'vst_norm_counts_hs_metadata.txt', sep = '/'))
vst_counts <- read.delim(paste(p_vst, 'vst_norm_counts_hs.txt', sep = '/')) %>%
  filter(X %in% unique(df$gene)) %>% 
  column_to_rownames(var = 'X')

# Create a matrix
hclust_matrix <- vst_counts %>% 
  as.matrix()

# z-score
sd <- apply(hclust_matrix, 1, sd, na.rm=TRUE)
hclust_matrix <- (hclust_matrix-rowMeans(hclust_matrix))/sd

# if your data has been transformed or normalized
# please skip this step
# hclust_matrix <- hclust_matrix %>% 
# transpose the matrix so genes are as columns
# t() %>% 
# apply scaling to each column of the matrix (genes)
# scale() %>% 
# transpose back so genes are as rows again
# t()

gene_dist <- dist(hclust_matrix)

require(fastcluster)
gene_hclust <- fastcluster::hclust(gene_dist, method = "complete")

# The default `plot()` function can be used to produce a simple dendrogram
plot(gene_hclust, labels = FALSE)
abline(h = 9, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram
# we chose 8 clusters after screening the dendrogram plot

cutree(gene_hclust, k = 8)
gene_cluster <- cutree(gene_hclust, k = 8) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  dplyr::rename(gene = name, cluster = value)

gene_cluster$cluster <- paste0('C', gene_cluster$cluster)
head(gene_cluster)
# write.table(gene_cluster, paste(p_file, 'all_genes_hclust.txt', sep = '/'),
#             row.names = FALSE, quote = FALSE, sep = '\t')


# plot: scatter plot
meta <- read.delim(paste(p_vst, 'vst_norm_counts_hs_metadata.txt', sep = '/')) %>%
  mutate(replicate = c(rep(c(1,2,3), 4), rep(c(4,1), 2), rep(c(5,2), 2), c(1,2,1,2))) %>% 
  mutate(name = paste(genotype, treatment, replicate, sep = '_'))
gene_cluster <- read.table('all_genes_hclust.txt', header = T)
pd <- data.frame(hclust_matrix)

all(names(pd) == meta$sample) # if True, then pass

colnames(pd) <- meta$name
head(pd)

pd <- pd %>% 
  rownames_to_column(var = 'gene') %>% 
  merge(., gene_cluster) %>% 
  relocate(gene, cluster) %>% 
  gather(condition, zscore, -gene, -cluster) %>% 
  separate(., col = 'condition', into = c('genotype', 'treatment', 'replicate'), sep = '_') %>% 
  group_by(cluster, genotype, treatment, replicate) %>% 
  summarise(median_zscore = median(zscore))

pd$cluster <- factor(pd$cluster, levels = paste0('C', 1:14))
pd$genotype <- factor(pd$genotype, levels = c('Tak1', 'hsfa', 'hsfb', 'dko'), labels = c('Tak-1', 'hsfa', 'hsfb', 'dko'))

cluster <- data.frame(table(gene_cluster$cluster)) %>% 
  dplyr::rename(cluster = Var1,
                count = Freq)
pd <- merge(pd, cluster) %>% 
  mutate(name = paste0(cluster, ' (n=', count, ')'))
pd$name <- factor(pd$name, levels = unique(pd[order(pd$cluster),]$name))
pd$treatment <- factor(pd$treatment, levels = c('NHS', 'heat'),
                       labels = c('Cntr', 'HS'))


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
       y = 'Z-score (median)',
       color = 'Treatment') +
  scale_color_manual(values = c('#999999', '#f1a340', 'black', 'purple')) +
  ylim(-2.0, 2.0) +
  facet_wrap(.~name, ncol = 3) +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text.x = element_text(size = 10, face = 'bold', angle = 40, hjust = 1),
        axis.text.y = element_text(size = 10, face = 'bold'),
        strip.background = element_rect(fill = alpha('#756bb1', 0.5)),
        strip.text = element_text(size = 8, face = 'bold'),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = 'black'))


## Figure 2-C
pdf('plots/cluster.pdf',
    width = 4, height = 5)
plt
dev.off()
###############################


####### OVERLAPPING GENES BETWEEN CLUSTERS IN ATAC-SEQ AND RNA-SEQ ######
rna_cluster <- read.delim('~/f3/all_genes_hclust.txt')
ocr_cluster <- read.delim('~/f1/counts/split/peak2gene.txt') %>% 
  dplyr::rename(gene = geneId)

head(rna_cluster)
head(ocr_cluster)

# we need to get the peak on promoter region
# peak location/area
Txdb_gtf <- makeTxDbFromGFF('~/f1/MpTak_v6.1r1.gff')
peak_cluster <- read.table('~/f1/counts/split/DE_peaks_hclust.txt', header = T)
head(peak_cluster)

data <- peak_cluster %>% 
  dplyr::select(peak) %>% 
  mutate(c = peak) %>% 
  separate(., col = 'c', into = c('Chr', 'Start', 'End'), sep = '\\.')

data <- merge(data, peak_cluster)

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

ocr_cluster <- peakAnno_tb %>% 
  dplyr::select(geneId, peak, cluster) %>% 
  dplyr::rename(gene = geneId)


p_export <- 'cluster_ovl'
# dir.create(p_export)

sum_cluster <- data.frame()
sum_cluster2 <- data.frame()
for (i in 1:8) {
  c_rna <- paste0('C', i)
  dir.create(paste(p_export, c_rna, sep = '/'))
  
  
  rna <- rna_cluster %>% 
    filter(cluster == c_rna) %>% 
    dplyr::rename(RNA_cluster = cluster)
  
  for (z in 1:15) {
    c_ocr <- paste0('C', z)
    ocr <- ocr_cluster %>% 
      filter(cluster == c_ocr) %>% 
      dplyr::rename(ATAC_cluster = cluster)
    
    ovl <- inner_join(rna, ocr) %>% 
      distinct()
    ovl_n <- length(unique(ovl$gene))
    
    # write.table(ovl, paste(p_export, c_rna, paste0('RNA_', c_rna, '_ATAC_', c_ocr, '.txt'), sep = '/'),
    #             row.names = F, quote = F, sep = '\t')
    sub <- data.frame(RNA_cluster = c_rna,
                      ATAC_cluster = c_ocr,
                      overlapping = ovl_n)
    
    # append data
    sum_cluster <- rbind(sum_cluster, sub)
    
    if(ovl_n != 0){
      sub2 <- data.frame(RNA_cluster = c_rna,
                         ATAC_cluster = c_ocr,
                         overlapping_gene = ovl$gene,
                         peak = ovl$peak)
      # append data2
      sum_cluster2 <- rbind(sum_cluster2, sub2)
    }
  }
}


# hypergeometric test
Total <- length(unique(ocr_cluster$gene, rna_cluster$gene))

ocr_gene <- data.frame(table(ocr_cluster$cluster)) %>% 
  dplyr::rename(ATAC_cluster = Var1,
                ATAC_cluster_genes = Freq)

rna_gene <- data.frame(table(rna_cluster$cluster)) %>% 
  dplyr::rename(RNA_cluster = Var1,
                RNA_cluster_genes = Freq)

df <- sum_cluster %>% 
  merge(., rna_gene) %>% 
  merge(., ocr_gene) %>% 
  mutate(hg = phyper(overlapping-1, ATAC_cluster_genes, Total-ATAC_cluster_genes, RNA_cluster_genes, lower.tail= FALSE)) %>% 
  mutate(`-log10(hg)` = -log10(hg)) %>% 
  mutate(ratio = overlapping/ATAC_cluster_genes)

df$sig <- ''
df[which(df$hg < 0.05),]$sig <- '*'

df$RNA_cluster <- factor(df$RNA_cluster, paste0('C', 1:8))
df$ATAC_cluster <- factor(df$ATAC_cluster, paste0('C', 1:15))

df.tem <- df
df.tem[which(df.tem$hg > 0.05),]$`-log10(hg)` <- -0.1

df.tem <- df.tem %>% 
  mutate(ATAC_name = paste0(ATAC_cluster, ' (',ATAC_cluster_genes, ')'),
         RNA_name = paste0(RNA_cluster, ' (', RNA_cluster_genes, ')'))

df.p <- df.tem %>% 
  dplyr::select(RNA_name, ATAC_name, ratio) %>% 
  spread(key = 'ATAC_name', value = ratio) %>% 
  column_to_rownames(var = 'RNA_name')


df.c <- df.tem %>% 
  dplyr::select(RNA_name, ATAC_name, overlapping) %>% 
  spread(key = 'ATAC_name', value = overlapping) %>% 
  column_to_rownames(var = 'RNA_name')

df.m <- df.tem %>% 
  dplyr::select(RNA_name, ATAC_name, sig) %>% 
  spread(key = 'ATAC_name', value = sig) %>% 
  column_to_rownames(var = 'RNA_name')


# export
write.table(df, 'hg.txt', row.names = F, quote = F, sep = '\t')
write.table(df.p, 'percentage.txt', quote = F, sep = '\t')
write.table(df.c, 'counts.txt', quote = F, sep = '\t')

col_fun = colorRamp2(c(0, 0.1, 0.2), c('darkturquoise', 'white', 'gold'))

column_order <- unique(df.tem[order(df.tem$ATAC_cluster),]$ATAC_name)
row_order <- unique(df.tem[order(df.tem$RNA_cluster, decreasing = T),]$RNA_name)

plt <- Heatmap(df.p,
             name = 'Fraction',
             rect_gp = gpar(col = 'white',
                            lwd = 0.5),
             col = col_fun,
             show_row_dend = F,
             show_column_dend = F,
             cluster_rows = F,
             cluster_columns = F,
             column_names_side = 'bottom',
             column_names_gp = gpar(fontsize = 8,
                                    fontface = 'bold',
                                    hjust = 1),
             column_names_rot = 40,
             row_title = 'Clusters of RNA-seq',
             column_title = 'Clusters of  ATAC-seq',
             clustering_distance_rows = 'pearson',
             row_names_side = 'left',
             row_names_gp = gpar(fontsize = 8, fontface = 'bold'),
             row_title_gp = gpar(fontsize = 10, fontface = 'bold'),
             column_title_gp = gpar(fontsize = 10, fontface = 'bold', hjust = 1),
             column_title_side = 'top',
             row_order = row_order,
             column_order = column_order,
             heatmap_legend_param = list(direction = 'horizontal'),
             layer_fun = function(j, i, x, y, width, height, fill, slice_r, slice_c){
               
               grid.rect(gp = gpar(lwd = 1.5, fill = "transparent"))
               
             },
             cell_fun = function(j, i, x, y, width, height, fill) {
               if(df.c[i, j] > 0)
                 grid.text(df.c[i, j], x, y, gp = gpar(fontsize = 6, fontface = 'bold'))
               grid.text(df.m[i, j], x*0.98, y*1.1, gp = gpar(fontsize = 8))
               
             }
)


## Figure 2-D
pdf('plots/cluster_ovl.pdf',
    width = 5, height = 4)  
draw(plt, heatmap_legend_side = 'bottom')
dev.off()
###############################


###### CORRELATION OF OVERLAPPING GENES #########
sum_cluster <- data.frame()
sum_cluster2 <- data.frame()
for (i in 1:8) {
  c_rna <- paste0('C', i)
  dir.create(paste(p_export, c_rna, sep = '/'))
  
  
  rna <- rna_cluster %>% 
    filter(cluster == c_rna) %>% 
    dplyr::rename(RNA_cluster = cluster,
                  geneId = gene)
  
  for (z in 1:15) {
    c_ocr <- paste0('C', z)
    ocr <- peakAnno_tb %>% 
      filter(cluster == c_ocr,
             annotation == 'Promoter') %>% 
      dplyr::rename(ATAC_cluster = cluster)
    
    ovl <- inner_join(rna, ocr, by = 'geneId') %>% 
      distinct()
    ovl_n <- nrow(ovl)
    
    write.table(ovl %>% dplyr::select(geneId, RNA_cluster, ATAC_cluster, peak, seqnames, start, end, width, distanceToTSS),
                paste(p_export, c_rna, paste0('RNA_', c_rna, '_ATAC_', c_ocr, '.txt'), sep = '/'),
                row.names = F, quote = F, sep = '\t')
    sub <- data.frame(RNA_cluster = c_rna,
                      ATAC_cluster = c_ocr,
                      overlapping = ovl_n)
    # append data
    sum_cluster <- rbind(sum_cluster, sub)
    
    
    if(ovl_n != 0){
      sub2 <- data.frame(RNA_cluster = c_rna,
                         ATAC_cluster = c_ocr,
                         overlapping_gene = ovl$geneId,
                         overlapping_peak = ovl$peak)
      
      # append data2
      sum_cluster2 <- rbind(sum_cluster2, sub2)
    }
    
  }
}


# expression matrix rna-seq
df <- read.delim('all_genes_hclust.txt') %>% 
  dplyr::select(gene) %>% 
  distinct()
p_vst <- 'vst_norm_counts'
meta <- read.delim(paste(p_vst, 'vst_norm_counts_hs_metadata.txt', sep = '/')) %>%
  mutate(replicate = c(rep(c(1,2,3), 4), rep(c(4,1), 2), rep(c(5,2), 2), c(1,2,1,2))) %>% 
  mutate(name = paste(genotype, treatment, replicate, sep = '_'))
vst_counts <- read.delim(paste(p_vst, 'vst_norm_counts_hs.txt', sep = '/')) %>%
  filter(X %in% unique(df$gene)) %>% 
  column_to_rownames(var = 'X')


# Create a matrix
hclust_matrix <- vst_counts %>% 
  as.matrix()

# z-score
sd <- apply(hclust_matrix, 1, sd, na.rm=TRUE)
hclust_matrix <- (hclust_matrix-rowMeans(hclust_matrix))/sd

# pd rna
pd <- data.frame(hclust_matrix)

all(names(pd) == meta$sample) # if True, then pass

colnames(pd) <- meta$name
head(pd)


# accessible matrix atac-seq
p2 <- '～/f1/counts/split'
norm_counts <- read.table(paste(p2, 'normalized_counts.txt', sep = '/'), header = T) %>% 
  filter(peak %in% peak_cluster$peak)

# Create a matrix
df2 <- norm_counts %>% 
  column_to_rownames(var = 'peak')
hclust_matrix2 <- df2 %>% 
  as.matrix()

# z-score
sd2 <- apply(hclust_matrix2, 1, sd, na.rm=TRUE)
hclust_matrix2 <- (hclust_matrix2-rowMeans(hclust_matrix2))/sd2

meta2 <- read.table(paste(p2, 'meta.txt', sep = '/'), header = T) %>% 
  mutate(name = paste(genotype, treatment, replicate, sep = '_'))
pd2 <- data.frame(hclust_matrix2)

all(names(pd2) == meta2$sample) # if True, then pass

colnames(pd2) <- meta2$name
head(pd2)

# calculate PCC
require(Hmisc)
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

rna_range <- paste0('C', 1:8)
atac_range <- paste0('C', 1:15)
pcc <- data.frame()
for (g in 1:length(rna_range)) {
  for (z in 1:length(atac_range)) {
    rz <- atac_range[z]
    rg <- rna_range[g]
    if(nrow(sum_cluster2 %>% filter(ATAC_cluster == rz,
                                    RNA_cluster == rg)) != 0){
      # rna
      chose <- sum_cluster2 %>% 
        filter(RNA_cluster == rg,
               ATAC_cluster == rz) %>%
        dplyr::rename(gene = overlapping_gene) %>% 
        dplyr::select(gene)
      
      pd_r <- pd %>% 
        rownames_to_column(var = 'gene') %>% 
        merge(., chose) %>% 
        relocate(gene) %>% 
        gather(condition, zscore, -gene) %>% 
        separate(., col = 'condition', into = c('genotype', 'treatment', 'replicate'), sep = '_') %>% 
        group_by(genotype, treatment) %>% 
        summarise(median_zscore = median(zscore))
      
      pd_r$genotype <- factor(pd_r$genotype, levels = c('Tak1', 'hsfa', 'hsfb', 'dko'))
      pd_r$treatment <- factor(pd_r$treatment, levels = c('NHS', 'heat'),
                               labels = c('NHS', 'HEAT'))
      
      # atac
      chose2 <- sum_cluster2 %>% 
        filter(RNA_cluster == rg,
               ATAC_cluster == rz) %>%
        dplyr::select(peak)
      
      pd_a <- pd2 %>% 
        rownames_to_column(var = 'peak') %>% 
        merge(., chose2) %>% 
        relocate(peak) %>% 
        gather(condition, zscore, -peak) %>% 
        separate(., col = 'condition', into = c('genotype', 'treatment', 'replicate'), sep = '_') %>% 
        group_by(genotype, treatment) %>% 
        summarise(median_zscore = median(zscore))
      
      pd_a$genotype <- factor(pd_a$genotype, levels = c('Tak1', 'hsfa', 'hsfb', 'dko'))
      pd_a$treatment <- factor(pd_a$treatment, levels = c('CK', 'HS'),
                               labels = c('NHS', 'HEAT'))
      
      # merge pd_r and pd_a
      pd_m <- merge(pd_r %>% dplyr::rename(rna = median_zscore),
                    pd_a %>% dplyr::rename(atac = median_zscore)) %>% 
        mutate(name = paste(genotype, treatment, sep = '_')) %>% 
        dplyr::select(-genotype, -treatment) %>% 
        column_to_rownames(., var = 'name')
      
      # correlation
      res_sub <- rcorr(as.matrix(pd_m), type = 'pearson')
      res2_sub <- flattenCorrMatrix(res_sub$r, res_sub$P) %>% 
        mutate(RNA_cluster = rg,
               ATAC_cluster = rz) %>% 
        dplyr::select(-row, -column) %>% 
        relocate(RNA_cluster, ATAC_cluster)
      
      # append data
      pcc <- rbind(pcc, res2_sub)
    }else{
      res2_sub <- data.frame(RNA_cluster = rg,
                             ATAC_cluster = rz,
                             cor = NA,
                             p = NA)
      # append data
      pcc <- rbind(pcc, res2_sub)
    }
    
  }
}

pcc$sig <- ''
pcc[which(pcc$p < 0.05),]$sig <- '*'
pcc[which(pcc$p < 0.01),]$sig <- '**'
pcc[which(pcc$p < 0.001),]$sig <- '***'

write.table(pcc, 'correlation_matrix.txt',
            row.names = F, quote = F, sep = '\t')

pcc$ATAC_cluster <- factor(pcc$ATAC_cluster, levels = paste0('C', 1:15))

plt <- pcc %>% 
  ggplot(.,
         aes(x = ATAC_cluster,
             y = RNA_cluster,
             fill = cor)) +
  geom_tile() +
  geom_text(data = pcc,
            aes(x = ATAC_cluster,
                y = RNA_cluster,
                label = sig)) +
  scale_fill_gradient2(low = 'darkturquoise',
                       mid = 'white',
                       high = 'gold',
                       midpoint = 0) +
  labs(x = 'Clusters of ATAC-seq',
       y = 'Clusters of RNA-seq', 
       fill = 'PCC of RNA-seq and ATAC-seq') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 8, face = 'bold', angle = 40, hjust = 1),
        axis.text.y = element_text(size = 8, face = 'bold'),
        axis.title  = element_text(size = 12, face = 'bold'),
        legend.direction = 'horizontal',
        legend.position = 'bottom',
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        axis.line = element_blank())


## Figure 2-E
pdf('plots/cor_p.pdf',
    width = 5, height = 4)
plt
dev.off()
#############################
