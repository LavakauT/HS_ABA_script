generateDistanceMatrix <- function(peak, gene) {

    chromosomes <- intersect(peak$chromosome, gene$chromosome)
    
    distanceMatrixPerChromosome <- lapply(chromosomes, function(chr) {
        
        geneChr <- gene[which(gene$chromosome == chr),]
        peakChr <- peak[which(peak$chromosome == chr),]
            
        distance <- matrix(sapply(peakChr$center, function(x) x-geneChr$tss), nrow=nrow(geneChr), ncol=nrow(peakChr))
        distance <- distance*geneChr$strand
        rownames(distance) <- geneChr$geneID
            
        return(distance)
    })

    return(distanceMatrixPerChromosome)
}