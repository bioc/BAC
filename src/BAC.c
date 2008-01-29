#include "util.h" 


double up_date_lambda(double sum_w_yy, double sum_w_y, double sum_w, int R, double mu, double gamma, double a_eps, double b_eps);

double up_date_lambda_baseline(double sum_w_yy, double sum_w_y, double sum_w, int R, double mu, double a_eps, double b_eps);

double up_date_mu(double sum_w_y1, double sum_w_y2, double sum_w1, double sum_w2, double gamma, double lambda1, double lambda2, double xi, double psi);

void up_date_gamma_metropolis(double sum_w_y, double sum_w, double mu, double *gamma, double lambda, double xi_plus, double psi_plus, double w, int *z);

double log_target_theta_tiling(double x, double gamma, int g, double *theta, int G, double lambda_theta, int nb_neighbors);

double slice_sampling_theta_tiling(double theta0, double w, int p, double gamma, int g, double *theta, int G, double lambda_theta, int nb_neighbors);

double up_date_lambda_theta_tiling(double *theta, int G, int nb_neighbors, double* joint_post_p, double* gamma, int k);

void up_date_normal_hyperpriors(double sum_xx, double sum_x, int length, double *xi, double *lambda);

void up_date_normal_hyperpriors_gamma(double sum_xx, double sum_x, int length, double *xi, double *lambda);

void BAC(double *y1_vec, double *y2_vec, int *R1, int *R2, int *G, double *mu,	double *gamma, int *model, double *post_p, double *joint_post_p, double *lambda1, double *a_eps1, double *b_eps1, 
  double *lambda2, double *a_eps2, double *b_eps2, int *B, int *verbose, int *nb_neighbors);

void callRegions(int *position, int *nProbes, double *dMerge, double *Score, double *cutoff, int *regions);


static const R_CMethodDef CEntries[] = {
  {"BAC", (DL_FUNC)&BAC, 20},
  {"callRegions", (DL_FUNC)&callRegions, 7},  
  {NULL, NULL, 0}
};


void BAC(double *y1_vec, double *y2_vec, int *R1, int *R2, int *G, double *mu,	double *gamma, int *model, double *post_p, double *joint_post_p, double *lambda1, double *a_eps1, double *b_eps1, 
double *lambda2, double *a_eps2, double *b_eps2, int *B, int *verbose, int *nb_neighbors)
{
	int k;
	double **y1;
	double **y2;
	int i,g;
	double *theta,*w;
	/* Initialize some values */
	double lambda_theta=0.01;
	double xi=0,psi=1.,xi_plus=4,psi_plus=0.15;
	double sum_lambda1=0;
	double sum_log_lambda1=0;
	double sum_lambda2=0;
	double sum_log_lambda2=0;
	double *sum_w_y1;
	double *sum_w_y2;
	double *sum_w_yy1;
	double *sum_w_yy2;  
	double *sum_w1;
	double *sum_w2;
	double sum_mu=0;
	double sum_mumu=0;
	double sum_gamma=0;
	double sum_gammagamma=0;
	int nb_gamma_plus=0;

	y1=dmatrix(*G, *R1);
	y2=dmatrix(*G, *R2);
	theta=dvector(*G,-10);
	w=dvector(*G,0.01);
	sum_w1=dvector(*G,0);
	sum_w2=dvector(*G,0);
	sum_w_y1=dvector(*G,0);
	sum_w_y2=dvector(*G,0);
	sum_w_yy1=dvector(*G,0);
	sum_w_yy2=dvector(*G,0);

	vec_mat(y1_vec,G,R1,y1);
	vec_mat(y2_vec,G,R2,y2);


	GetRNGstate();


	for(g=0;g<*G;g++)
	{
		for(i=0;i<*R1;i++)
		{
			sum_w1[g]+=1;
			sum_w_y1[g]+=y1[g][i];
			sum_w_yy1[g]+=y1[g][i]*y1[g][i];
		}

		for(i=0;i<*R2;i++)
		{
			sum_w2[g]+=1;
			sum_w_y2[g]+=y2[g][i];
			sum_w_yy2[g]+=y2[g][i]*y2[g][i];
		}
	}
	for(k=0;k<(*B+1000);k++)
	{ 


		if((((k+1)*100)%(10**B)==0) & (*verbose==1))
		{
			printf("%d percent completed \n",(((k+1)*100)/(*B)));
		}

		theta[0]=-10;
		for(g=1;g<*G;g++)
		{
			theta[g]=slice_sampling_theta_tiling(theta[g], 0.1, 100, gamma[g], g, theta, *G, lambda_theta, *nb_neighbors);
		}

		sum_lambda1=0;
		sum_log_lambda1=0;
		sum_lambda2=0;
		sum_log_lambda2=0;
		sum_mu=0;
		sum_mumu=0;
		sum_gamma=0;
		sum_gammagamma=0;
		nb_gamma_plus=0;
		for(g=0;g<*G;g++)
		{


			/*	w[g]=1./(1.+exp(-theta[g]));*/
			w[g]=exp(-theta[g])/(1.+exp(-theta[g]));

			if(k>=1000)
				post_p[g]+=model[g];
			mu[g]=up_date_mu(sum_w_y1[g], sum_w_y2[g], sum_w1[g], sum_w2[g],gamma[g], lambda1[g], lambda2[g], xi, psi);	      
			up_date_gamma_metropolis(sum_w_y2[g], sum_w2[g], mu[g], gamma+g,lambda2[g], xi_plus, psi_plus, w[g], model+g);

			lambda1[g]=up_date_lambda_baseline(sum_w_yy1[g],sum_w_y1[g], sum_w1[g], *R1, mu[g], *a_eps1, *b_eps1);
			lambda2[g]=up_date_lambda(sum_w_yy2[g], sum_w_y2[g], sum_w2[g], *R2, mu[g], gamma[g], *a_eps2, *b_eps2);

			sum_mu+=mu[g];
			sum_mumu+=mu[g]*mu[g];

			if(gamma[g]!=0)
			{
				sum_gamma+=gamma[g];
				sum_gammagamma+=gamma[g]*gamma[g];	      
				nb_gamma_plus++;
			}

			sum_lambda1+=lambda1[g];
			sum_log_lambda1+=log(lambda1[g]);
			sum_lambda2+=lambda2[g];
			sum_log_lambda2+=log(lambda2[g]);
		}

		up_date_normal_hyperpriors(sum_mumu, sum_mu, *G, &xi, &psi);

		lambda_theta=up_date_lambda_theta_tiling(theta, *G, *nb_neighbors, joint_post_p, gamma, k);
		up_date_normal_hyperpriors_gamma(sum_gammagamma, sum_gamma, nb_gamma_plus, &xi_plus, &psi_plus);
/*		The mean of the background is set to zero*/
		xi=0;

		*a_eps1=slice_sampling_a(*a_eps1, 2, 50, sum_log_lambda1, sum_lambda1, *b_eps1, *G);
		*b_eps1=slice_sampling_b(*b_eps1, 2, 50, sum_log_lambda1, sum_lambda1, *a_eps1, *G);
		*a_eps2=slice_sampling_a(*a_eps2, 2, 50, sum_log_lambda2, sum_lambda2, *b_eps2, *G);
		*b_eps2=slice_sampling_b(*b_eps2, 2, 50, sum_log_lambda2, sum_lambda2, *a_eps2, *G);

	}

	PutRNGstate();

  // free(w);
  // free_dmatrix(y1, *G);
  // free_dmatrix(y2, *G);
  // free(theta);
  // free(sum_w1);
  // free(sum_w2);
  // free(sum_w_y1);
  // free(sum_w_y2);
  // free(sum_w_yy1);
  // free(sum_w_yy2);
}

void callRegions(int *position, int *nProbes, double *dMerge, double *Score, double *cutoff, int *regions)
{
	int i=0,p=0;
	int w=0,ww=0,max=0;
	int nRegions=0;

	while(p<*nProbes)
	{
		if(Score[p]>*cutoff)
		{
			nRegions++;
			regions[p]=nRegions;
			w=p+1;
			max=p;
			ww=p;

			while((w<*nProbes) && ((position[w]-position[ww])<*dMerge))
			{
				if(Score[w]>*cutoff)
				{
					max=w;
					ww=max;
				}
				w++;
			}
			for(i=p;i<=max;i++)
			{
				regions[i]=nRegions;
			}
			p=w;
		}
		else
		{
			regions[p]=0;
			p++;
		}
	}
}


double up_date_mu(double sum_w_y1, double sum_w_y2, double sum_w1, double sum_w2, double gamma, double lambda1, double lambda2, double xi, double psi)
{

	double sum_tmp=0;
	double mean=0;

	mean=sum_w_y1*lambda1+sum_w_y2*lambda2-sum_w2*lambda2*gamma+psi*xi;
	sum_tmp=sum_w1*lambda1+sum_w2*lambda2+psi;

	return(rnorm(mean/sum_tmp,1./sqrt(sum_tmp)));
}


double up_date_lambda_baseline(double sum_w_yy, double sum_w_y, double sum_w, int R, double mu, double a_eps, double b_eps)
{

	double SSR=0;

	SSR=sum_w_yy-2*sum_w_y*mu+sum_w*mu*mu;

	return(rgamma(R/2.+a_eps*a_eps/b_eps,1./(SSR/2.+a_eps/b_eps)));
}


double up_date_lambda(double sum_w_yy, double sum_w_y, double sum_w, int R, double mu, double gamma, double a_eps, double b_eps)
{

	double SSR=0;

	SSR=sum_w_yy-2*sum_w_y*(mu+gamma)+sum_w*(mu+gamma)*(mu+gamma);

	return(rgamma(R/2.+a_eps*a_eps/b_eps,1./(SSR/2.+a_eps/b_eps)));
}


void up_date_gamma_metropolis(double sum_w_y, double sum_w, double mu, double *gamma, double lambda, double xi_plus, double psi_plus, double w, int *z)
{

	double sum_tmp=0;
	double mean=0;
	double one_minus_w_star=0;
	double gamma_new=0;
	int z_new;
	double log_ratio=0;


/** Compute the normalizing constant **/

	sum_tmp=sum_w*lambda+psi_plus;
	mean=lambda*(sum_w_y-sum_w*mu)+psi_plus*xi_plus;

	one_minus_w_star=(1.-w)/(1.-w+w*pnorm(0,mean/sum_tmp,1./sqrt(sum_tmp),0,0)/pnorm(0,xi_plus,1./sqrt(psi_plus),0,0)*sqrt(psi_plus)/sqrt(sum_tmp)*exp(0.5*mean*mean/sum_tmp));


	if(runif(0,1)<one_minus_w_star)
	{
		z_new=0;
		gamma_new=0;
	}
	else
	{
		gamma_new=rnorm(mean/sum_tmp,1./sqrt(sum_tmp));
		z_new=1;
	}

	if(gamma_new>=0)
	{
		if(gamma_new>0)
		{
			log_ratio=0;
		}
		else if(gamma_new==0)
		{
			if(*gamma>0)
				log_ratio=pnorm(0,mean/sum_tmp,1./sqrt(sum_tmp),0,1);	  
		}
		if(log(runif(0,1))<log_ratio)
		{
			*z=z_new;
			*gamma=gamma_new;
		}
	}
}


void up_date_normal_hyperpriors(double sum_xx, double sum_x, int length, double *xi, double *lambda)
{
	double sum_residuals=0;

	sum_residuals=sum_xx-2*sum_x**xi+length**xi**xi;

	*lambda=rgamma(length/2.+1.,1./(sum_residuals/2.+0.0001));

	*xi=rnorm((sum_x**lambda+0.0001*0)/(length**lambda+0.0001),1./sqrt(length**lambda+0.0001));
}

void up_date_normal_hyperpriors_gamma(double sum_xx, double sum_x, int length, double *xi, double *lambda)
{
	double new_value=0;
	double log_ratio=0;
	double trunc=0;

	new_value=*xi+rnorm(0,0.25);

	log_ratio+=-length*pnorm(trunc,new_value,1./sqrt(*lambda),0,1)+0.5*length*log(*lambda)-0.5**lambda*(sum_xx-2*sum_x*new_value+length*new_value*new_value);
	log_ratio-=-length*pnorm(trunc,*xi,1./sqrt(*lambda),0,1)+0.5*length*log(*lambda)-0.5**lambda*(sum_xx-2*sum_x**xi+length**xi**xi);

	log_ratio+=dunif(new_value,0,10,1)-dunif(*xi,0,10,1);

	if((log_ratio>log(runif(0,1))))
		*xi=new_value;

	log_ratio=0;
	new_value=*lambda+rnorm(0,.5);

	log_ratio+=-length*pnorm(trunc,*xi,1./sqrt(new_value),0,1)+0.5*length*log(new_value)-0.5*new_value*(sum_xx-2*sum_x**xi+length**xi**xi);
	log_ratio-=-length*pnorm(trunc,*xi,1./sqrt(*lambda),0,1)+0.5*length*log(*lambda)-0.5**lambda*(sum_xx-2*sum_x**xi+length**xi**xi);

	log_ratio+=dgamma(new_value,1,1./1.,1)-dgamma(*lambda,1,1./1.,1);

	if((log_ratio>log(runif(0,1))) & (new_value>0) & (new_value<1.0))
		*lambda=new_value;
}


double log_target_theta_tiling(double x, double gamma, int g, double *theta, int G, double lambda_theta, int nb_neighbors)
{
	double log_target=0;
	int i,nb=0;
	double mean=0;

	int start=fmax2(0,g-nb_neighbors);
	int end=fmin2(G-1,g+nb_neighbors);

	for(i=start;i<=end;i++)
	{
		if(i!=g)
		{
			mean+=theta[i];
			nb++;
		}
	}

	log_target=dnorm(x,mean/nb,1./sqrt(lambda_theta*nb/(2.*nb_neighbors)),1.);

	if(gamma>0)
		log_target+=-log(1.+exp(x));
	else
		log_target+=x-log(1.+exp(x));

	return(log_target);

}

double slice_sampling_theta_tiling(double theta0, double w, int p, double gamma, int g, double *theta, int G, double lambda_theta, int nb_neighbors)
{

	double L, R, z, theta1;
	int J,K;
	double log_f_R=0,log_f_L=0;


	z=log_target_theta_tiling(theta0, gamma, g, theta, G, lambda_theta, nb_neighbors)-rgamma(1.,1.);

/** lower bound **/
	L=theta0-w*runif(0,1);
	R=L+w;
	J=(int)(p*runif(0,1));
	K=(p-1)-J;


/** Find the interval **/
	log_f_R=log_target_theta_tiling(R, gamma, g, theta, G, lambda_theta, nb_neighbors);
	log_f_L=log_target_theta_tiling(L, gamma, g, theta, G, lambda_theta, nb_neighbors);

	while((J>0) & (z<log_f_L))
	{
		L=L-w;
		log_f_L=log_target_theta_tiling(L, gamma, g, theta, G, lambda_theta, nb_neighbors);

		J--;
	}

	while((K>0) & (z<log_f_R))
	{
		R=R+w;
		log_f_R=log_target_theta_tiling(R, gamma, g, theta, G, lambda_theta, nb_neighbors);

		K--;
	}

	theta1=runif(L,R);
	while(log_target_theta_tiling(theta1, gamma, g, theta, G, lambda_theta, nb_neighbors)<z)
	{
	/** Shrink the interval **/
		if(theta1<theta0)
			L=theta1;
		else
			R=theta1;
		theta1=runif(L,R);
	}

	return(theta1);
}

double up_date_lambda_theta_tiling(double *theta, int G, int nb_neighbors, double* joint_post_p, double* gamma, int k)
{
	int g=0,i=0;
	double sum=0;
	int nb=0;
	int start,end;
	double mean=0;
	double sum_xy=0;
	double sum_model=0;

	for(g=0;g<G;g++)
	{
		nb=0;
		mean=0;
		sum_xy=0;
		sum_model=0;
		start=fmax2(0,g-nb_neighbors);
		end=fmin2(G-1,g+nb_neighbors);

		for(i=start;i<=end;i++)
		{
			if(i!=g)
			{
				mean+=theta[i];
				nb++;
				sum_xy+=theta[g]*theta[i];
			}
			if(gamma[i]>0)
			{
				sum_model++;
			}
		}

		if((k>=1000) && (sum_model>nb_neighbors))
		{
			joint_post_p[g]++;
		}
/*	sum+=theta[g]*(nb/nb_neighbors*theta[g]-mean/nb_neighbors);*/
		sum+=nb/(2.*nb_neighbors)*theta[g]*theta[g]-sum_xy/(2.*nb_neighbors);
	}
	return(rgamma(G/2.+0.001,1./(sum/2.+0.001)));
}


void R_init_BAC(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    // R_useDynamicSymbols(dll, FALSE);
}
