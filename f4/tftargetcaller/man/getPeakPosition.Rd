\name{getPeakPosition}
\alias{getPeakPosition}
\title{Process a BED file to get the information on the positions of peaks centers and intensities}
\description{Process a file in BED format to extract peak locations and find the centers and intensities of the peaks.}
\usage{getPeakPosition(filename, skip=0, chr=1, start=2, end=3, summit=NULL, intensity=NULL)}
\arguments{
	\item{filename}{Name of the BED file. BED files should contain at least three columns with chromosome, peak start position and peak end position respectively. Optionally there might be the fourth column containing the position of the peak summit and the fifth column containing the peak intensity. There should be no header and no comment lines. If they occur they can be omitted using the \code{skip} argument.}
	\item{skip}{Integer specifying number of lines to be skipped from the top of the BED file. The default value is \code{0}.}
	\item{chr}{Integer specifying which column contains chromosome location. The default value is \code{1}.}
	\item{start}{Integer specifying which column contains peak start. The default value is \code{2}.}
	\item{end}{Integer specifying which column contains peak end. The default value is \code{3}.}
	\item{summit}{Integer specifying which column contains peak summit. If specified, the summit will be used as a peak location parameter. If \code{NULL}, the peak location parameter will be calculated as the center of the peak. The default value is \code{NULL}.}
	\item{intensity}{Integer specifying which column contains peak intensity. The default value is \code{NULL}.}}
% \details{}
\value{A data.frame containing two (chromosome and center) or three (chromosome, center and intensity) columns.}
% \references{}
\note{If intensity is not provided, some of the methods cannot be used.}
% \author{}
\seealso{\code{peakPosition}}
\examples{
# library path might have to be adapted
peakPosition <- getPeakPosition(filename=paste(.libPaths()[1], "/TFTargetCaller/data/bedfile.bed", sep=""), intensity=4)
}
% \keyword{}
