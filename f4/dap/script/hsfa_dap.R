###################
library(dplyr)
library(tibble)
library(tidyverse)
library(stringr)
library(Hmisc)
library(GenomicFeatures)
library(GenomicRanges)
library(ChIPseeker)
###################


###################
path <- 'f4/tftargetcaller/R'
files <- list.files(path, full.names = T, recursive = F)
for (file in files) {
  source(file)
  print(paste0("Load R source: ", file))
}

setwd('f4/tftargetcaller/data')
load('peakPosition.RData')
load('genePosition.RData')
load('chromSizes.RData')

# # Binary method
# TF_target_score <- TFTargetCaller(peakPosition, genePosition, method="Binary", n=5000)
# 
# # Linear method
# TF_target_score <- TFTargetCaller(peakPosition, genePosition, method="Linear", n=50000)
# 
# # Ouyang method
# TF_target_score <- TFTargetCaller(peakPosition, genePosition, method="Ouyang", n=5000)
# 
# # Cheng method
# TF_target_score <- TFTargetCaller(paste(.libPaths()[1],"/TFTargetCaller/data/wigfile.wig",sep=""), genePosition, method="Cheng", n=10000, chr.len=chromSizes, ChengScore="pvalue")
# 
# # Chen method
# TF_target_score <- TFTargetCaller(peakPosition, genePosition, method="Chen", chr.len=chromSizes)
# 
# # ClosestGene method returning score
# TF_target_score <- TFTargetCaller(peakPosition, genePosition, method="ClosestGene", ClosestGeneScore="score")
# 
# # ClosestGene method returning q-value
# TF_target_score <- TFTargetCaller(peakPosition, genePosition, method="ClosestGene", ClosestGeneScore="qvalue")


Txdb_gtf <- makeTxDbFromGFF('f1/MpTak_v6.1r1.gff')
tss <- transcripts(Txdb_gtf)
tss_matrix <- data.frame(tss[1:length(tss)])
tss_matrix$gene_id <- str_remove(tss_matrix$tx_name, '\\.\\d')
tss_matrix_first <- tss_matrix %>% 
  filter(grepl("\\.1", tx_name))

gene <- genes(Txdb_gtf)
gene_matrix <- data.frame(gene[1:length(gene)])

all(nrow(tss_matrix_first) == length(gene))

genePosition_mar <- tss_matrix_first %>% 
  mutate(geneID = gene_id) %>% 
  column_to_rownames(., var = 'gene_id') %>% 
  dplyr::rename(chromosome = seqnames,
                tss = start) %>% 
  dplyr::select(geneID, chromosome, strand, tss)
genePosition_mar$strand <- as.character(genePosition_mar$strand)
genePosition_mar[which(genePosition_mar$strand == "+"),]$strand <- 1
genePosition_mar[which(genePosition_mar$strand == "-"),]$strand <- -1

setwd('f4/dap/macs2_narroPeak_copy')
# ClosestGene method
chrs <- unique(genePosition_mar$chromosome)
genePosition_mar$chromosome <- as.character(genePosition_mar$chromosome)
genePosition_mar$strand <- as.numeric(genePosition_mar$strand)

peakfile_1 <- read.delim('HSFA_22CY5FLT4.txt') %>% 
  dplyr::rename(chromosome = chr,
                center = abs_summit,
                intensity = pileup) %>% 
  filter(chromosome %in% chrs)
target_score_1 <- cbind(data.frame(TFTargetCaller(peakfile_1, genePosition_mar, method="ClosestGene", ClosestGeneScore="score")),
                        data.frame(TFTargetCaller(peakfile_1, genePosition_mar, method="ClosestGene", ClosestGeneScore="qvalue")))
gc()
names(target_score_1) <- c('score', 'qvalue')
target_score_1_q <- target_score_1 %>%
  arrange(qvalue) %>% 
  dplyr::slice(1:2000)

peakfile_2 <- read.delim('A_22CV7GLT4.txt') %>% 
  dplyr::rename(chromosome = chr,
                center = abs_summit,
                intensity = pileup) %>% 
  filter(chromosome %in% chrs)
target_score_2 <- cbind(data.frame(TFTargetCaller(peakfile_2, genePosition_mar, method="ClosestGene", ClosestGeneScore="score")),
                        data.frame(TFTargetCaller(peakfile_2, genePosition_mar, method="ClosestGene", ClosestGeneScore="qvalue")))
gc()
names(target_score_2) <- c('score', 'qvalue')
target_score_2_q <- target_score_2%>%
  arrange(qvalue) %>% 
  dplyr::slice(1:2000)

peakfile_3 <- read.delim('50_22HML7LT3.txt') %>% 
  dplyr::rename(chromosome = chr,
                center = abs_summit,
                intensity = pileup) %>% 
  filter(chromosome %in% chrs)
target_score_3 <- cbind(data.frame(TFTargetCaller(peakfile_3, genePosition_mar, method="ClosestGene", ClosestGeneScore="score")),
                        data.frame(TFTargetCaller(peakfile_3, genePosition_mar, method="ClosestGene", ClosestGeneScore="qvalue")))
gc()
names(target_score_3) <- c('score', 'qvalue')
target_score_3_q <- target_score_3 %>%
  arrange(qvalue) %>% 
  dplyr::slice(1:2000)


save(target_score_1, target_score_2, target_score_3, target_score_1_q, target_score_2_q, target_score_3_q,
     file = 'f4/dap/macs2_narroPeak_copy/hsfa_target.RData')


load('f4/dap/macs2_narroPeak_copy/hsfa_target.RData')

data <- rbind(data.frame(gene = row.names(target_score_1_q),
                         group = 'HSFA_22CY5FLT4'),
              data.frame(gene = row.names(target_score_2_q),
                         group = 'A_22CV7GLT4'),
              data.frame(gene = row.names(target_score_3_q),
                         group = '50_22HML7LT3')) %>% 
  mutate(count = 'V') %>% 
  spread(., key = 'group', value = 'count')

data_intersect <- data[which(!is.na(data$`50_22HML7LT3`) & !is.na(data$A_22CV7GLT4) & !is.na(data$HSFA_22CY5FLT4)),]

# export
write.table(data, 'f4/dap/macs2_narroPeak_copy/hsfa_targets.txt',
            row.names = F, quote = F, sep = '\t')
write.table(data_intersect, 'f4/dap/macs2_narroPeak_copy/hsfa_targets_intersect.txt',
            row.names = F, quote = F, sep = '\t')




# DAP-seq nodes functional analysis----------------
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


order.id <- data_intersect
order.sub.ara <- ortho4 %>% filter(Mpolymorpha %in% order.id$gene)
order.sub.ara <- as.data.frame(str_split(order.sub.ara$Athaliana, ", ", n = Inf, simplify = TRUE))
order.sub.ara <- gather(order.sub.ara, key = "X", value = "geneId")
order.sub.ara <- order.sub.ara %>%
  dplyr::select(geneId) %>% 
  distinct()

write.table(order.sub.ara,
            'f4/dap/macs2_narroPeak_copy/edge_dap_ara.txt',
            row.names = F,
            quote = F,
            sep = '\t')


## go reduced input
read.delim('f4/dap/macs2_narroPeak_copy/dap_ara_go.txt') %>% 
  filter(q.value < 0.05) %>% 
  dplyr::select(GO.ID, q.value) %>% 
  write.table(., 'f4/dap/macs2_narroPeak_copy/dap_ara_go_input.txt',
              row.names = F, quote = F, sep = '\t')


## loading go ontology
colname <- c('GO.ID', 'Term', 'Value', 'Count', 'Group')
df <- data.frame(matrix(nrow = 0, ncol = length(colname)))
names(df) <- colname
df.ori <- read.delim('f4/dap/macs2_narroPeak_copy/dap_ara_go.txt')
df.rvi <- read.delim('f4/dap/macs2_narroPeak_copy/dap_ara_go_revi.tsv') %>% 
  filter(Uniqueness < 0.9) %>% 
  dplyr::rename(GO.ID = TermID,
         Term = Name)
df.rvi <- df.rvi[order(df.rvi$Value),]%>% 
  dplyr::select(GO.ID, Term, Value)
df.ori <- merge(df.rvi, df.ori, by = c('GO.ID', 'Term')) %>% 
  dplyr::select(GO.ID, Term, Value, Count) %>% 
  mutate(Group = 'HSF3_DAP_target')
df <- rbind(df, df.ori)

# export
write.table(df, 'f4/dap/macs2_narroPeak_copy/dap_go.txt', row.names = F, quote = F, sep = '\t')

# heatmap
library(circlize)
library(ComplexHeatmap)
col_fun = colorRamp2(c(-0.1, 0, 3, 7.5), c('grey60', 'white', 'darkturquoise', 'gold'))


df <- read.delim('f4/dap/macs2_narroPeak_copy/dap_go_revise.txt')
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

p_go <- Heatmap(df.p,
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

pdf('f4/dap/macs2_narroPeak_copy/dap_go.pdf',
    width = 3.75, height = 3)
p_go
dev.off()






# export dap-seq target to TF-GRN------------------------
dap_target <- read.delim('f4/dap/macs2_narroPeak_copy/hsfa_targets.txt') %>% 
  dplyr::select(gene) %>% 
  distinct()

mr.tf <- read.delim('f4/MpTak_v6.1_TF.txt') %>%
  distinct(.keep_all = TRUE)

dap_target_tf <- dap_target %>% 
  mutate(ID = gene) %>% 
  merge(., mr.tf) %>% 
  dplyr::select(-ID)

# export edge files
# merge with network
nodes <- read.csv('f4/hclust/genie3_node_ocr2_clean.csv')
edges <- read.csv('f4/hclust/edge_table_atac_rna.csv')



dap_target_tf <- dap_target_tf %>%
  dplyr::rename(geneId = gene) %>% 
  dplyr::select(geneId, family) %>% 
  mutate(source = 'MpVg00470') %>% 
  dplyr::rename(target = geneId,
                target_TF = family) %>% 
  relocate(source, target) %>% 
  distinct() %>% 
  filter(target %in% nodes$name) %>%
  mutate(weight = 0.12,
         presence = 'DAP-seq_target',
         presence_rna = 'source',
         source_TF = 'HSF')

write.csv(dap_target_tf, 'f4/hclust/dap_target.csv',
          row.names = F, quote = F)


# revise edges
edges <- read.csv('f4/hclust/edge_table_atac_rna.csv')
edges$presence_rna_revised <- edges$presence_rna
edges[which(edges$presence != edges$presence_rna & edges$presence != '' & edges$presence != 'RNA-seq'),]$presence_rna_revised <- edges[which(edges$presence != edges$presence_rna & edges$presence != '' & edges$presence != 'RNA-seq'),]$presence

write.csv(edges,'f4/hclust/edge_table_atac_rna.csv', row.names = F, quote = F)
##################
