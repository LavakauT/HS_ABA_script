getGenePosition <- function(species="Mus musculus", type=c("protein_coding", "lincRNA", "miRNA"), archiveHost=NULL, identifiers=NULL) {
  
  #if archiveHost is set, it allows access to a prior version of the ensembl database 
  #e.g. for hg18 the paramenter would be "may2009.archive.ensembl.org"
  if (!is.null(archiveHost)){
    mart <- useMart("ENSEMBL_MART_ENSEMBL",
                    host=archiveHost,
                    path="/biomart/martservice",
                    archive=FALSE)
  } else {
    mart <- useMart("ensembl") 
  }
  
  tmp <- listDatasets(mart)
  dataset <- as.character(tmp[grep(species, tmp[,2]),1])
  mart <- useDataset(dataset, mart=mart)
  
  info <- getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id", "chromosome_name", "strand", "transcript_start", "transcript_end", "gene_biotype",identifiers), mart=mart)
  
  if (! is.null(type)) {
    info <- info[which(info$gene_biotype %in% type),]
  }
  
  info <- cbind(info, tss=findTSS(info$strand, info$transcript_start, info$transcript_end))
  
  geneID2tss <- split(info$tss, info$ensembl_gene_id)
  geneID2strand <- split(info$strand, info$ensembl_gene_id)
  geneID2strand <- lapply(geneID2strand, unique)
  geneID2strand <- unlist(geneID2strand)
  geneID2tss <- mapply(function(a,b) ifelse(b==1, min(a), max(a)), geneID2tss, geneID2strand[names(geneID2tss)])
  
  tmpNames <- colnames(info)[8:(length(colnames(info))-1)]
  info <- data.frame(geneID=names(geneID2tss), chromosome=info[match(names(geneID2tss), info$ensembl_gene_id), "chromosome_name"], strand=info[match(names(geneID2tss), info$ensembl_gene_id), "strand"], tss=geneID2tss, info[match(names(geneID2tss), info$ensembl_gene_id),8:(length(colnames(info))-1)], stringsAsFactors=FALSE)
  
  colnames(info)[5:length(colnames(info))] <- tmpNames
  
  info$chromosome <- sapply(info$chromosome, function(x) paste(c("chr", x), collapse=""))
  
  return(info)
}