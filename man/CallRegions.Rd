\name{CallRegions}
\alias{CallRegions}

\title{Call and merge regions using joint posterior probabilities calculated by BAC.}

\description{
Call and merge regions using joint posterior probabilities calculated by BAC.
}

\usage{
CallRegions(position,jointPP,cutoff=0.5,maxGap=500)
}

\arguments{
\item{position}{A vector containing the probe genomic positions}
\item{jointPP}{A vector containing the joint posterior probabilities as returned by BAC.}
\item{cutoff}{The cutoff used to call regions.}
\item{maxGap}{The maximum gap allowed between regions. Regions that are less than maxGap bps away will be merged.}
}

\details{

}

\value{
A vector containing the region index for each probe. Probes with the same positive index belong to the same region, whereas probe with index zero are 
background probes (not part of a bound region). These indices can be used to form a BED file, see example below.
}

\author{
Raphael Gottardo, \email{raph@stat.ubc.ca}
}

\note{
}

\seealso{
BAC
}

\examples{
# Load the data
data(ER)
# Only select the first 5000 probes for speed-up
ER<-ER[1:5000,]
# Calculate the joint posterior probabilities
#Only use 100 iterations for speed up (You should use more! See default value) 
BAConER<-BAC(ER[,5:7], ER[,2:4], B=100,verbose=FALSE,w=5)
# For Regions using 0.5 cut-off for the joint posterior probabilities
ERregions<-CallRegions(ER[,1],BAConER$jointPP,cutoff=0.5,maxGap=500)
# Create the BED file
nRegions<-max(ERregions)
BED<-matrix(0,nRegions,4)
for(i in 1:nRegions)
{
BED[i,2:3]<-range(ER[ERregions==i,1])
#The score should be between 0 and 1000
BED[i,4]<-max(BAConER$jointPP[ERregions==i])*1000
}
BED<-data.frame(BED)
# The ER data is a subset of chr 21
BED[,1]<-"chr21"
names(BED)<-c("chrom","chromStart","chromEnd","Score")
# print it
print(BED)
}

\keyword{models}

