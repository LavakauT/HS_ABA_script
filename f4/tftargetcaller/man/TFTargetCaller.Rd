\name{TFTargetCaller}
\alias{TFTargetCaller}
\title{TF-target scoring}
\description{Assign a TF-target score using selected method.}
\usage{TFTargetCaller(peak, gene, method="ClosestGene", n=5000, chr.len=NULL, smooth=TRUE, ChengScore="pvalue", ClosestGeneScore="score", randomizations=500)}
\arguments{
	\item{peak}{For \code{Binary}, \code{Linear}, \code{Chen} and \code{ClosestGene}: a data.frame with at least two columns (chromosome and center) containing information on peak positions. Fore \code{Ouyang} a data.frame with three columns (chromosome, center and intensity) containing information on peak positions and intensities. For \code{Cheng} a character string specifying filename of the .wig file. Package contains an exemplary .wig file for transcription factor CEBPa from HemoChIP compendium.}
	\item{gene}{A data.frame with four columns (geneID, chromosome, strand, tss) containing information on gene positions.}
	\item{method}{Character string specifying which method to use. The available methods are: \code{Binary}, \code{Linear}, \code{Ouyang}, \code{Cheng}, \code{Chen}, \code{ClosestGene}. The default is \code{ClosestGene}.}
	\item{n}{An integer specifying the value of the parameter required by some of the methods: window size around the TSS [-n, n] for \code{Binary}, \code{Linear} and \code{Cheng} methods and constant for \code{Ouyang} method. The default is \code{5000}. Original \code{Cheng} and \code{Ouyang} methods were used with n=10000 and n=5000, respectively.}
	\item{chr.len}{A named integer vector specifying lengths of the chromosomes. An exemplary vector containing the lengths of mouse chromosomes extracted from UCSC is stored in the object \code{chromSizes}. Used for \code{Cheng} and \code{Chen} methods.}
	\item{smooth}{A boolean specifying if smoothing should be applied. Used for \code{Cheng} method.}
	\item{ChengScore}{Character string specifying which score should be extracted for \code{Cheng} method. \code{score}, \code{zscore} and \code{pvalue} correspond to the original score, z-score and p-value respectively. The default is \code{pvalue}.}
	\item{ClosestGeneScore}{Character string specifying which score should be extracted for \code{ClosestGene} method. \code{score} and \code{qvalue} correspond to the original score and q-value respectively. The default is \code{score}.}
	\item{randomizations}{An integer specifying number of randomizations to be performed for q-value calculation. The default is \code{500}.}}
% \details{}
\value{A vector containing the score for each TF-target pair.}
\references{Cheng C. et al.
    \emph{TIP: a probabilistic method for identifying transcription factor target genes from ChIP-seq binding profiles.} Bioinformatics. 2011 Dec 1;27(23):3221-7,
    Chen X. et al.
    \emph{Integration of external signaling pathways with the core transcriptional network in embryonic stem cells.} Cell. 2008 Jun 13;133(6):1106-17,
    Ouyang Z. et al.
    \emph{ChIP-Seq of transcription factors predicts absolute and differential gene expression in embryonic stem cells.} Proc Natl Acad Sci U S A. 2009 Dec 22;106(51):21521-6,
    Hannah R. et al.
    \emph{A compendium of genome-wide hematopoietic transcription factor maps supports the identification of gene regulatory control mechanisms.} Exp Hematol. 2011 May;39(5):531-41}
% \note{}
% \author{}
% \seealso{}
\examples{
	data(peakPosition)
	data(genePosition)
	data(chromSizes)
	
	# Binary method
	TF_target_score <- TFTargetCaller(peakPosition, genePosition, method="Binary", n=5000)
	
	# Linear method
	TF_target_score <- TFTargetCaller(peakPosition, genePosition, method="Linear", n=50000)
	
	# Ouyang method
	TF_target_score <- TFTargetCaller(peakPosition, genePosition, method="Ouyang", n=5000)
	
	# Cheng method
	TF_target_score <- TFTargetCaller(paste(.libPaths()[1],"/TFTargetCaller/data/wigfile.wig",sep=""), genePosition, method="Cheng", n=10000, chr.len=chromSizes, ChengScore="pvalue")
	
	# Chen method
	TF_target_score <- TFTargetCaller(peakPosition, genePosition, method="Chen", chr.len=chromSizes)
	
	# ClosestGene method returning score
	TF_target_score <- TFTargetCaller(peakPosition, genePosition, method="ClosestGene", ClosestGeneScore="score")
	
	# ClosestGene method returning q-value
	TF_target_score <- TFTargetCaller(peakPosition, genePosition, method="ClosestGene", ClosestGeneScore="qvalue")
}
% \keyword{}
