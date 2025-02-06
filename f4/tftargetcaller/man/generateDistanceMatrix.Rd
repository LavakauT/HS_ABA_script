\name{generateDistanceMatrix}
\alias{generateDistanceMatrix}
\title{Generate a gene-to-peak distance matrix for each chromosome}
\description{Generate a matrix of distances between all genes and all peaks for each chromosome.}
\usage{generateDistanceMatrix(peak, gene)}
\arguments{
	\item{peak}{A data.frame with at least two columns (chromosome, center) containing information on peak positions.}
	\item{gene}{A data.frame with four columns (geneID, chromosome, strand, tss) containing information on gene positions.}}
% \details{}
\value{A list with elements corresponding to chromosomes. Each element is a matrix of gene-to-peak distances.}
% \references{}
% \note{}
% \author{}
\seealso{\code{getPeakPosition}, \code{getGenePosition}}
\examples{
data(peakPosition)
data(genePosition)
distMat <- generateDistanceMatrix(peakPosition, genePosition)
}
% \keyword{}
