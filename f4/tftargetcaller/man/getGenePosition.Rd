\name{getGenePosition}
\alias{getGenePosition}
\title{Get information on gene positions from BioMart}
\description{Uses the getBM function from the biomaRt package to get the information on gene positions from the Ensembl database. Function allows to restrict transcript types to be used in the analysis. For genes with several transcripts the most 5' position of their TSSs is considered as the gene's TSS.}
\usage{getGenePosition(species="Mus musculus", type=c("protein_coding", "lincRNA", "miRNA"))}
\arguments{
    \item{species}{String representing the organism in which the target calling is performed, e.g. for mouse \code{Mus musculus} and for human \code{Homo sapiens}. The default value is \code{Mus musculus}. List of all supported organisms can be obtained with \code{biomaRt} function \code{listDatasets(mart="ensembl")}.}
    \item{type}{Types of transcripts included in the analysis. By default: \code{protein_coding}, \code{lincRNA} and \code{miRNA}. If \code{type=NULL} no filtering is performed.}}
% \details{}
\value{A data.frame containing four columns: geneID, chromosome, strand and tss.}
\references{Steffen Durinck and Wolfgang Huber
    \emph{biomaRt: Interface to BioMart databases (e.g. Ensembl, COSMIC, Wormbase and Gramene).}}
% \note{}
% \author{}
\seealso{\code{getBM}, \code{\link{findTSS}}, \code{\link{genePosition}}}
\examples{
genePosition <- getGenePosition()
genePosition <- getGenePosition("Homo sapiens", type="protein_coding")
}
% \keyword{}
