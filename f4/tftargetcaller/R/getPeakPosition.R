getPeakPosition <- function(filename, skip=0, chr=1, start=2, end=3, summit=NULL, intensity=NULL) {

    input <- read.table(filename, skip=skip, header=FALSE, as.is=TRUE)
    
    if (! is.null(summit)) {
        output <- data.frame(chromosome=input[,chr], center=input[,summit], stringsAsFactors=FALSE)
    } else {
        output <- data.frame(chromosome=input[,chr], center=rowMeans(input[,c(start,end)]), stringsAsFactors=FALSE)
    }
    
    if (! is.null(intensity)) output <- cbind(output, intensity=input[,intensity])
    
    return(output)
}