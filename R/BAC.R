"BAC"<-function(C, I, B=15000,verbose=FALSE,w=5)
{

	mu<-apply(C,1,"mean")
	gamma<-apply(I,1,"mean")-mu
	G<-dim(C)[1]
	R1<-dim(C)[2]
	R2<-dim(I)[2]

	model<-as.integer(gamma>1)
	gamma[gamma<=1]<-0

	lambda1<-1/(apply(C,1,"var")+0.001)
	lambda2<-1/(apply(I,1,"var")+0.001)

	a.eps1<-median(lambda1)
	b.eps1<-mad(lambda1)
	a.eps2<-median(lambda2)
	b.eps2<-mad(lambda2)

	obj<-.C("BAC",
	as.double(t(C)),
	as.double(t(I)),
	as.integer(R1),
	as.integer(R2),
	as.integer(G),
	as.double(mu),
	as.double(gamma),
	as.integer(model),
	margPP=as.double(rep(0,G)),
	jointPP=as.double(rep(0,G)),
	as.double(lambda1),
	as.double(a.eps1),
	as.double(b.eps1),
	as.double(lambda2),
	as.double(a.eps2),
	as.double(b.eps2),
	as.integer(B),
	as.integer(verbose),
	as.integer(w),
  package="BAC")

	list(margPP=obj$margPP/B,jointPP=obj$jointPP/B)
}
