\name{findTSS}
\alias{findTSS}
\title{Find the position of a gene's TSS}
\description{Find the position of a gene's TSS given its strand, start and end position. Function can be applied both to gene and transcript positions.}
\usage{findTss(strand, start, end)}
\arguments{
	\item{strand}{A vector containing strand information for genes or transcripts. Strand has to be encoded as 1 (forward) or -1 (reverse).}
	\item{start}{A vector containing start positions of genes or transcripts.}
	\item{end}{A vector containing end positions of genes or transcripts.}}
\details{\code{strand}, \code{start} and \code{end} must be of the same length and the entries should correspond to the same genes or transcripts.}
\value{A vector containing the TSS positions of genes or transcripts.}
% \references{}
% \note{}
% \author{}
\seealso{\code{\link{getGenePosition}}}
\examples{
exampleGenes <- vector(mode="list")
exampleGenes$name <- c("Ldb1","Gfi1b","Tal1","Gata1","Lmo2","Rag2","Lyl1","Cebpa","Cbfa2t","Erg","Slc22a","Foxo1","Cebpb","Ebf1","Stat4")
exampleGenes$strand <- c(-1,-1,1,-1,1,1,1,1,1,-1,1,1,1,1,1) 
exampleGenes$start <- c(46032593,28609450,115056426,7959260,103957986,101624718,84701457,35119293,154436481,95359169,12584189,52268337,167688915,44617317,51987148)
exampleGenes$end <- c(46045214,28621982,115071755,7978071,103981878,101632529,84704715,35121926,154539356,95586593,12628465,52350109,167690418,45008091,52107188)
exampleGenes$tss <- findTSS(exampleGenes$strand, exampleGenes$start, exampleGenes$end)
}
% \keyword{}
