TFTargetCaller <- function(peak, gene, method="ClosestGene", chromosomes=c(1:19, "X"), n=5000, chr.len=NULL, smooth=TRUE, ChengScore="pvalue", ClosestGeneScore="score", randomizations=500) {
    
    if ((! all(gene$strand %in% c(1,-1))) & (! all(gene$strand %in% c("+", "-")))) stop("Gene strand must be given as (1, -1) or (+, -)")
    gene$strand[which(gene$strand == "+")] <- 1
    gene$strand[which(gene$strand == "-")] <- -1
    gene$strand <- as.numeric(gene$strand)
    
    if (method %in% c("Binary","Linear","Ouyang","Chen","ClosestGene")) if (length(grep("chr", peak$chromosome)) == 0) peak$chromosome <- sapply(peak$chromosome, function(x) paste(c("chr", x), collapse=""))
    
    if (length(grep("chr", gene$chromosome)) == 0) gene$chromosome <- sapply(gene$chromosome, function(x) paste(c("chr", x), collapse=""))

    chromosomes <- unique(gene$chromosome)
    
    if (method %in% c("Binary","Linear","Ouyang","Chen")) distanceMatrix <- generateDistanceMatrix(peak, gene)

    if (method == "Cheng") {
    
        if (is.null(chr.len)) stop("Argument chr.len must be specified!")
        if (! all(names(chr.len) %in% chromosomes)) warning(paste(c("Chromosome length was not provided or chromosomes: ", paste(setdiff(chromosomes, names(chr.len)), collapse=", ")), collapse=""))
        chr.len <- chr.len[intersect(names(chr.len), chromosomes)]
        chr.nam <- names(chr.len)
        
        mygene <- data.frame(name=gene$geneID, chr=gene$chromosome, str=gene$strand, sta=gene$tss, end=gene$tss, stringsAsFactors=FALSE)
        mygene[which(mygene$str == 1),"sta"] <- mygene[which(mygene$str == 1),"sta"] - n
        mygene[which(mygene$str == 1),"end"] <- mygene[which(mygene$str == 1),"end"] + n
        mygene[which(mygene$str == -1),"sta"] <- mygene[which(mygene$str == -1),"sta"] + n
        mygene[which(mygene$str == -1),"end"] <- mygene[which(mygene$str == -1),"end"] - n
        
        conIn = file(peak, "r")
        data = readLines(conIn, -1)
        close(conIn)
        
        pos = grep("fixedStep", data)
        info = data[pos]
        pos = c(pos, length(data)+1)
        info = unlist(strsplit(info, " "))
        cnum = length(info)/4
        mychr = info[(1:cnum)*4-2]
        mysta = info[(1:cnum)*4-1]
        mychr = gsub("chrom=", "", mychr)
        mysta = as.numeric(gsub("start=", "", mysta))
        pos.sta = pos[1:(length(pos)-1)]+1
        pos.end = pos[2:length(pos)]-1
        sig.tk = data.frame(mychr, mysta, pos.sta, pos.end)
        
        myw = rep(0, n*2+1)
        
        for(k in 1:length(chr.nam))
        {
            # cat("\r", chr.nam[k])
            # signal vector
            read.cov = rep(0, chr.len[k])
            tag = sig.tk[,1]==chr.nam[k]
            if(sum(tag)==0) next
            mysig = sig.tk[sig.tk[,1]==chr.nam[k], ]
            for(i in 1:nrow(mysig))
            {
                tmp1 = mysig[i,3]   # start line
                tmp2 = mysig[i,4]   # end line
                tmp.sta = mysig[i,2]  # start position in chr
                read.cov[(tmp.sta+1):(tmp.sta+tmp2-tmp1+1)] = as.numeric(data[tmp1:tmp2])
            }

            curgene = mygene[mygene[,2]==chr.nam[k],]
            for(i in 1:nrow(curgene))
            {
                tmp.myw = read.cov[curgene[i,4]:curgene[i,5]]
                tmp.myw[which(is.na(tmp.myw))] = 0
                myw = myw+tmp.myw
            }
        }
        
        myw = myw/nrow(mygene)
        tmp = myw
        if(smooth==T)
        {
            for(i in 1:length(myw))
            {
                myw[i] = mean(tmp[max(i-250,1):min(i+250, length(myw))])	
            }
        }
        myw = myw/sum(myw)
        
        mysco = rep(0, nrow(mygene))
        gname = rep("",nrow(mygene))
        count = 0
        for(k in 1:length(chr.nam))
        {
            # cat("\r", chr.nam[k])
            # signal vector
            read.cov = rep(0, chr.len[k])
            tag = sig.tk[,1]==chr.nam[k]
            if(sum(tag)==0) next
            mysig = sig.tk[sig.tk[,1]==chr.nam[k], ]
            for(i in 1:nrow(mysig))
            {
                tmp1 = mysig[i,3]   # start line
                tmp2 = mysig[i,4]   # end line
                tmp.sta = mysig[i,2]  # start position in chr
                read.cov[(tmp.sta+1):(tmp.sta+tmp2-tmp1+1)] = as.numeric(data[tmp1:tmp2])
            }

            curgene = mygene[mygene[,2]==chr.nam[k],]
            for(i in 1:nrow(curgene))
            {
                count = count + 1
                tmp = read.cov[curgene[i,4]:curgene[i,5]]
                tmp[which(is.na(tmp))] = 0
                mysco[count] = sum(tmp*myw)
                gname[count] = curgene[i,1]
            }
        }
        if (ChengScore == "score") score <- mysco
        if (ChengScore == "zscore") score <- (mysco-mean(mysco))/sd(mysco)
        if (ChengScore == "pvalue") score <- pnorm(-((mysco-mean(mysco))/sd(mysco)))
        names(score) <- gname
        score <- score[gene$geneID]
    }
    
    if (method == "Chen") {
    
        boundaries <- c(-10^12, -100000, -50000, -20000, -10000, -5000, -2000, -1000, 0, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 10^12)
    
        if (is.null(chr.len)) stop("Argument chr.len must be specified!")
        if (! all(names(chr.len) %in% chromosomes)) warning(paste(c("Chromosome length was not provided for chromosomes: ", paste(setdiff(chromosomes, names(chr.len)), collapse=", ")), collapse=""))
        chr.len <- chr.len[intersect(names(chr.len), chromosomes)]
        chr.fraction <- chr.len/sum(as.numeric(chr.len))
        
        distanceVectorRandom <- c()
        
        for (i in 1:3) {
            
            fractionOfPeaks <- round(chr.fraction*nrow(peak))
            
            chromosome2positions <- mapply(function(u,v,w) igraph.sample(low=u, high=v, length=w), structure(.Data=rep(1,length(fractionOfPeaks)),names=names(fractionOfPeaks)), chr.len[names(fractionOfPeaks)], fractionOfPeaks)
            
            peakRandom <- data.frame(chromosome=rep(names(chromosome2positions), times=sapply(chromosome2positions, length)), center=unlist(chromosome2positions), stringsAsFactors=FALSE)
            
            distanceMatrixRandom <- generateDistanceMatrix(peakRandom, gene)
            
            distanceVectorRandom <- c(distanceVectorRandom, unlist(lapply(distanceMatrixRandom, function(mat) apply(mat, 1, function(x) x[which(abs(x) == min(abs(x)))][1]))))
        }
        
        histogramRandom <- hist(distanceVectorRandom, breaks=boundaries, plot=FALSE)$counts/3
        
        distanceVector <- unlist(lapply(distanceMatrix, function(mat) apply(mat, 1, function(x) x[which(abs(x) == min(abs(x)))][1])))
        
        histogram <- hist(distanceVector, breaks=boundaries, plot=FALSE)$counts
        
        histFunction <- (histogram-histogramRandom)/histogram
        histFunction[histFunction < 0] <- 0
        histFunction[1] <- 0
        histFunction[length(histFunction)] <- 0
        
        scoreCalculated <- sapply(distanceVector, function(x) max(hist(x, breaks=boundaries, plot=FALSE)$counts*histFunction))
        
        score <- rep(0,nrow(gene))
        names(score) <- gene$geneID
        score[names(scoreCalculated)] <- scoreCalculated 
    }
    
    if (method == "Ouyang") {
    
        if (! "intensity" %in% colnames(peak)) stop("Peak intensity is missing!")
        
        peakIntensity <- lapply(intersect(peak$chromosome, gene$chromosome), function(chr) {
            if (chr %in% peak$chromosome && chr %in% gene$chromosome) return(peak[which(peak$chromosome == chr),"intensity"])
        })
        
        scoreCalculated <- unlist(mapply(function(a,b) {
            a <- abs(a)
            a[a > 10^6] <- NA
            result <- apply(a, 1, function(x) {
                y <- b*exp(-x/n)
                return(sum(y, na.rm=TRUE))
            })
            return(result)
        }, distanceMatrix, peakIntensity))
        
        score <- rep(0,nrow(gene))
        names(score) <- gene$geneID
        score[names(scoreCalculated)] <- scoreCalculated
    }
    
    if (method == "Binary") {
    
        scoreCalculated <- unlist(lapply(distanceMatrix, function(x) {
            result <- apply(x, 1, function(y) ifelse(any(abs(y) < n), 1, 0))
            return(result)
        }))
        
        score <- rep(0,nrow(gene))
        names(score) <- gene$geneID
        score[names(scoreCalculated)] <- scoreCalculated
    }
    
    if (method == "Linear") {
    
        scoreCalculated <- unlist(lapply(distanceMatrix, function(x) {
            result <- apply(x, 1, function(y) {
                y <- 1-abs(y)/n
                z <- sum(y[y > 0])
                return(z)
            })
            return(result)
        }))
        
        score <- rep(0,nrow(gene))
        names(score) <- gene$geneID
        score[names(scoreCalculated)] <- scoreCalculated
    }
    
    if (method == "ClosestGene") {
    
        score <- NULL
        if (ClosestGeneScore == "qvalue") scoreRandom <- lapply(1:randomizations, function(i) as.numeric())      
        distribution <- NULL
        
        for (chr in chromosomes) {
        
            if (chr %in% peak$chromosome && chr %in% gene$chromosome) {
            
                geneChr <- gene[which(gene$chromosome == chr),]
                peakChr <- peak[which(peak$chromosome == chr),]
                
                distance <- matrix(sapply(peakChr$center, function(x) x-geneChr$tss), nrow=nrow(geneChr), ncol=nrow(peakChr))
                distance <- distance*geneChr$strand
                rownames(distance) <- geneChr$geneID
                
                distribution <- c(distribution, as.numeric(distance[abs(distance) <= 10^6]))
                
                distance <- apply(distance, 2, function(x) {
                    x[which(abs(x) != min(abs(x)))] <- NA
                    return(x)
                })
                
                # 20241016; modified by Lavakau for only one gene bug
                if(!is.matrix(distance)) {
                  scoreTemp <- list(gene = distance)
                  names(scoreTemp) <- geneChr$geneID
                  score <- c(score, scoreTemp)
                } else{
                  scoreTemp <- apply(distance, 1, function(x) x[which(! is.na(x))])
                  scoreTemp <- scoreTemp[which(sapply(scoreTemp, length) > 0)]
                  score <- c(score, scoreTemp)
                }
                
                rm(distance)
                gc()
                
                if (ClosestGeneScore == "qvalue") {
                    
                    scoreRandomTemp <- lapply(1:randomizations, function(i) split(unlist(scoreTemp, use.names=FALSE), sample(geneChr$geneID, length(unlist(scoreTemp)), replace=TRUE)))
                        
                    scoreRandom <- mapply(function(a,b) c(a,b), scoreRandom, scoreRandomTemp, SIMPLIFY=FALSE)
                }
            }
        }
        
        upstream <- distribution[distribution <= 0]
        downstream <- distribution[distribution >= 0]
    
        check <- (0 %in% round(distribution, 2))
    
        rm(distribution)
        gc()
    
        upstreamEcdf <- Ecdf(-upstream, pl=FALSE)
        downstreamEcdf <- Ecdf(downstream, pl=FALSE)
    
        rm(upstream, downstream)
        gc()
    
        x <- c(-upstreamEcdf$x[2:length(upstreamEcdf$x)], 0, downstreamEcdf$x[2:length(downstreamEcdf$x)])
        y <- c(upstreamEcdf$y[2:length(upstreamEcdf$y)], ifelse(check, NA, 0), downstreamEcdf$y[2:length(downstreamEcdf$y)])
    
        rm(upstreamEcdf, downstreamEcdf)
        gc()
    
        scoringFunction <- approxfun(x, y, yleft=1, yright=1)
        
        tmp <- unlist(score)
        tmp <- scoringFunction(tmp)
        tmp <- relist(tmp, score)
        scoreFinal <- sapply(tmp, prod)
        scoreComplete <- rep(1, nrow(gene))
        names(scoreComplete) <- gene$geneID
        scoreComplete[names(scoreFinal)] <- scoreFinal
        
        if (ClosestGeneScore == "qvalue") {
        
            scoreRandomComplete <- lapply(scoreRandom, function(item) {
                tmp <- unlist(item)
                tmp <- scoringFunction(tmp)
                tmp <- relist(tmp, item)
                itemFinal <- sapply(tmp, prod)
                itemComplete <- rep(1, nrow(gene))
                names(itemComplete) <- gene$geneID
                itemComplete[names(itemFinal)] <- itemFinal
                return(itemComplete)
            })
            
            #calculation of q-values. Here the empirical cumulative distribution function 
            #is calculated for the score calculated by randomizations
            randomDistribution <- ecdf(as.numeric(unlist(scoreRandomComplete)))
            
            #because calculatd p-values can be very small, which can result in 0 values,
            #the 0 p-values are replaced with 1/randomizations
            pvalues <- randomDistribution(scoreComplete)
            pvalues[pvalues==0] <- 1/randomizations
            qvalue <- structure(.Data=p.adjust(pvalues, method="fdr"), names=names(scoreComplete))
        }
        
        if (ClosestGeneScore == "score") {
            scoreComplete[which(scoreComplete == 0)] <- min(scoreComplete[which(scoreComplete > 0)])
            scoreComplete <- -log10(scoreComplete)
            score <- scoreComplete
        }
        
        if (ClosestGeneScore == "qvalue") score <- qvalue
    }
    
    return(score)
}