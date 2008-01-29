CallRegions<-function(position,jointPP,cutoff=0.5,maxGap=500)
{

	obj<-.C("callRegions",
	as.integer(position),
	as.integer(length(position)),
	as.double(maxGap),
	as.double(jointPP),
	as.double(cutoff),
	regions=as.integer(rep(0,length(position))),
	package="BAC")

	obj$regions
}
