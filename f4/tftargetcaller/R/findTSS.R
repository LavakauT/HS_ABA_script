findTSS <- function(strand, start, end) {

    tss <- mapply(function(o, s, e) ifelse(o == 1, s, e), strand, start, end)
    
    return(tss)
}