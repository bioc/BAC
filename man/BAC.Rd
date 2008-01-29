\name{BAC}
\alias{BAC}

\title{Bayesian Analysis of ChIP-chip tiling arrays}

\description{
Bayesian Analysis of ChIP-chip tiling arrays 
}

\usage{
BAC(C, I, B=15000,verbose=FALSE,w=5)
}

\arguments{
\item{C}{The matrix of control measurements. Rows correspond to probes and columns to samples.}
\item{I}{The matrix of IP measurements. Rows correspond to probes and columns to samples.}
\item{B}{Number of iterations used the MCMC. Default to 15000.}
\item{verbose}{Logical parameter. If TRUE, some progression }
\item{w}{The window size. Default to 5. See details below for more about this parameter.}
}

\details{
The window size should be calculated in function of the resolution and the shearing resolution. 
For example, for Affymetrix human tiling arrays, the shearing resolution is 500-1000bps, the tiling resolution is 35bps and the probe length is 25bps. 
Then one would expect a bound region to contain 500-1000/(35+25)~8-16 probes. Thus we decided to set w to 5. Note that the exact value of w is not crucial. 
}

\value{
The marginal posterior probabilities and the joint posterior probabilities computed from the Bayesian hierarchical model. 
We recommend using the joint posterior probabilities to call enriched regions.
}

\author{
Raphael Gottardo, \email{raph@stat.ubc.ca}
}

\note{
}

\seealso{
CallRegions
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

