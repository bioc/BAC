#include "util.h"



double code_miss=-9999999;
double c=9;


double stdd(double *vector,int *length, int *is_finite)
{
/* Compute the standard deviation of a vector */

	int i,count=0;
	double sum=0;
	double x_bar;
	double result;


	x_bar=mean_vec(vector,length);
	if(x_bar==code_miss)
		return(code_miss);
	else
	{
		for(i=0;i<*length;i++)
		{
			if(vector[i]!=code_miss)
			{
				count=count+1;
				sum=sum+(vector[i]-x_bar)*(vector[i]-x_bar);
			}
		}
		*is_finite=count;
		if(count>1)
		{
			result=(sqrt(sum/((double)count-1)));
			return(result);
		}
		else
		{
			return(code_miss);
		}
	}
}

void vec_mat(double *array_vec,int* nb_row,int *nb_col,double **array)
{
	int i,j;

	for(i=0;i<*nb_row;i++)
		for(j=0;j<*nb_col;j++)
		array[i][j]=array_vec[i**nb_col+j];
}

void mat_vec(double *array_vec,int* nb_row,int *nb_col,double **array)
{
	int i,j;

	for(i=0;i<*nb_row;i++)
		for(j=0;j<*nb_col;j++)
		array_vec[i**nb_col+j]=array[i][j];
}

void init_ivector(int *vector,int *length, int  value)
{
	int i;

	for(i=0;i<*length;i++)
		vector[i]=value;
}

double **dmatrix(int nb_row,int nb_col)
{
	double **array;
	int i,j;

/* Allocate the memory */
	array=(double **)R_alloc(nb_row, sizeof(double));

	for(i=0;i<nb_row;i++)
		array[i]=(double *)R_alloc(nb_col, sizeof(double));

/* Initialize to zero*/
	for(i=0;i<nb_row;i++)
		for(j=0;j<nb_col;j++)
		array[i][j]=0;

	return(array);
}

double *dvector(int length, int init)
{
	int i;
	double *vector;

/* Allocate the memory */
	vector=(double *)R_alloc(length, sizeof(double));

/* Initialize the memory */
	for(i=0;i<length;i++)
		vector[i]=init;

	return(vector);
}



int *ivector(int length, int init)
{
	int i;
	int *vector;

/* Allocate the memory */
	vector=(int *)R_alloc(length, sizeof(int));

/* Initialize the memory */
	for(i=0;i<length;i++)
		vector[i]=init;

	return(vector);
}



double  mean_vec(double *vector,int *length)
{
	int i,count=0;
	double sum=0;

	for(i=0;i<*length;i++)
	{
		if(vector[i]!=code_miss)
		{
			count=count+1;
			sum=sum+vector[i];
		}
	}
	if (count>0)
	{
		return(sum/(double)count);
	}
	else
	{
		return(code_miss);
	}

}


void free_dmatrix(double **array, int nb_row)
{

	int i;
	for(i=0;i<nb_row;i++)
		free(array[i]);

	free(array);
}


void quicksort2(double *a, int *b, int *p, int *r)
{
	int q;
	int q_p;

	if (*p<*r)
	{
		q=rand_part2(a,b, *p, *r);
		quicksort2(a,b,p,&q);
		q_p=q+1;
		quicksort2(a,b,&q_p,r);
	}
}
int partition2(double *a, int *b, int p, int r)
{
	double x=a[p];
	int i=p-1;
	int j=r+1;
	double temp;
	int itemp;

	for(;;)
	{
		do
		{
			j--;
		}while(a[j]>x);
		do
		{
			i++;
		}while(a[i]<x);
		if(i<j)
		{
			temp=a[i];
			a[i]=a[j];
			a[j]=temp;
			itemp=b[i];
			b[i]=b[j];
			b[j]=itemp;
		}
		else
		{
			return(j);
		}
	}

}

int rand_part2(double *a, int *b,  int p, int r)
{
	int i;
	double temp;
	int itemp;

	i=uni_rand(p,r);
	temp=a[p];
	a[p]=a[i];
	a[i]=temp;
	itemp=b[p];
	b[p]=b[i];
	b[i]=itemp;

	return(partition2(a,b,p,r));

}


void idquicksort2(int *a, double *b, int *p, int *r)
{
	int q;
	int q_p;

	if (*p<*r)
	{
		q=idrand_part2(a,b, *p, *r);
		idquicksort2(a,b,p,&q);
		q_p=q+1;
		idquicksort2(a,b,&q_p,r);
	}
}
int idpartition2(int *a, double *b, int p, int r)
{
	int x=a[p];
	int i=p-1;
	int j=r+1;
	double temp;
	int itemp;

	for(;;)
	{
		do
		{
			j--;
		}while(a[j]>x);
		do
		{
			i++;
		}while(a[i]<x);
		if(i<j)
		{
			itemp=a[i];
			a[i]=a[j];
			a[j]=itemp;
			temp=b[i];
			b[i]=b[j];
			b[j]=temp;
		}
		else
		{
			return(j);
		}
	}

}

int idrand_part2(int *a, double *b,  int p, int r)
{
	int i;
	double temp;
	int itemp;

	i=uni_rand(p,r);
	itemp=a[p];
	a[p]=a[i];
	a[i]=itemp;
	temp=b[p];
	b[p]=b[i];
	b[i]=temp;

	return(idpartition2(a,b,p,r));

}

int uni_rand(int min,int max)
{
	int rand_nb;


	GetRNGstate();
	rand_nb=(int)(unif_rand()*(max+1-min)+min);
	PutRNGstate();

	return(rand_nb);

}


void init_dvector(double *vector, int *length, int value)
{
	int i;


/* Initialize to value */
	for(i=0;i<*length;i++)
		vector[i]=value;
}

double  sum_vec(double *vector,int *length)
{
	int i,count=0;
	double sum=0;

	for(i=0;i<*length;i++)
	{
		if(vector[i]!=code_miss)
		{
			count=count+1;
			sum=sum+vector[i];
		}
	}
	if (count>0)
	{
		return(sum);
	}
	else
	{
		return(code_miss);
	}
}


void quicksort(double *a, int *p, int *r)
{
	int q;
	int q_p;

	if (*p<*r)
	{
		q=rand_part(a, *p, *r);
		quicksort(a,p,&q);
		q_p=q+1;
		quicksort(a,&q_p,r);
	}
}


int partition(double *a, int p, int r)
{
	double x=a[p];
	int i=p-1;
	int j=r+1;
	double temp;

	for(;;)
	{
		do
		{
			j--;
		}while(a[j]>x);
		do
		{
			i++;
		}while(a[i]<x);
		if(i<j)
		{
			temp=a[i];
			a[i]=a[j];
			a[j]=temp;
		}
		else
			return(j);
	}

}

int rand_part(double *a, int p, int r)
{
	int i;
	double temp;

	i=uni_rand(p,r);
	temp=a[p];
	a[p]=a[i];
	a[i]=temp;
	return(partition(a,p,r));

}

double dabs(double a)
{
	if (a<0)
		return(-a);
	else
		return(a);
}

double dmax(double a, double b)
{
	if(a<b)
		return(b);
	else
		return(a);
}

void product_matrix(double **mat1, int *n1, int *n2, double **mat2, int *m1, int *m2, double **res)
{

/** res has to be n1*n2 **/
	int i,j,k;
	double sum;

	for(i=0;i<*n1;i++)
		for(j=0;j<*m2;j++)
	{
		sum=0;
		for(k=0;k<*n2;k++)
			sum+=mat1[i][k]*mat2[k][j];
		res[i][j]=sum;
	}
}

double product_vec_vec(double *vec1, int *n1, double *vec2)
{

	int i;
	double sum=0;

	for(i=0;i<*n1;i++)
		sum=vec1[i]*vec2[i]+sum;

	return(sum);
}

void product_mat_vec(double **mat, int *n1, int *n2, double *vec, double *res)
{

/** res has to be n1*n2 **/
	int i,j;
	double sum;

	for(i=0;i<*n1;i++)
	{
		sum=0;
		for(j=0;j<*n2;j++)
			sum+=mat[i][j]*vec[j];
		res[i]=sum;
	}
}

double rexp_trunc(double lambda, double min, double max)
	/* generate a random deviate from a truncated exponential */
{

	return(-1/lambda*log(runif(exp(-lambda*max),exp(-lambda*min))));

}

double dexp_trunc(double x, double lambda, double min, double max)
	/* generate a random deviate from a truncated exponential */
{
	if(x<min || x>max)
		return(0);
	else
		return(lambda*exp(-lambda*x)/(exp(-lambda*min)-exp(-lambda*max)));
}

double log2(double x)
{
	return(log(x)/log(2.));
}

double log_f_ab(double sum_log_lambda, double sum_lambda, double a, double b, int n)
{

	double log_f;

	log_f=(a*a/b-1.)*sum_log_lambda-a/b*sum_lambda+n*a*a/b*log(a/b)-n*lgammafn(a*a/b)+dunif(a,0,10000.,1.)+dunif(b,0,10000.,1.);

	return(log_f);

}

double slice_sampling_b(double b0, double w, int p, double sum_log_lambda, double sum_lambda, double a, int n)
{

	double L, R, z, b1;
	int J,K;
	double log_f_R=0,log_f_L=0;


	z=log_f_ab(sum_log_lambda, sum_lambda,a,b0,n)-rgamma(1.,1.);
/** lower bound **/
	L=b0-w*runif(0,1);
	R=L+w;
	J=(int)(p*runif(0,1));
	K=(p-1)-J;


/** Find the interval **/
	log_f_R=log_f_ab(sum_log_lambda, sum_lambda,a,R,n);
	log_f_L=log_f_ab(sum_log_lambda, sum_lambda,a,L,n);

	while((J>0) & (z<log_f_L))
	{
		L=L-w;
		log_f_L=log_f_ab(sum_log_lambda, sum_lambda,a,L,n);

		J--;
	}

	while((K>0) & (z<log_f_R))
	{
		R=R+w;
		log_f_R=log_f_ab(sum_log_lambda, sum_lambda,a,R,n);

		K--;
	}

/** Make sure it is in the range **/
	L=fmax2(0.,L);
	R=fmin2(10000.,R);

	b1=runif(L,R);
	while(log_f_ab(sum_log_lambda, sum_lambda, a, b1, n)<z)
	{
	/** Shrink the interval **/
		if(b1<b0)
			L=b1;
		else
			R=b1;
		b1=runif(L,R);
	}

	return(b1);

}

double slice_sampling_a(double a0, double w, int p, double sum_log_lambda, double sum_lambda, double b, int n)
{

	double L, R, z, a1;
	int J,K;
	double log_f_R=0,log_f_L=0;


	z=log_f_ab(sum_log_lambda, sum_lambda,a0,b,n)-rgamma(1.,1.);
/** lower bound **/
	L=a0-w*runif(0,1);
	R=L+w;
	J=(int)(p*runif(0,1));
	K=(p-1)-J;


/** Find the interval **/
	log_f_R=log_f_ab(sum_log_lambda, sum_lambda,R,b,n);
	log_f_L=log_f_ab(sum_log_lambda, sum_lambda,L,b,n);

	while((J>0) & (z<log_f_L))
	{
		L=L-w;
		log_f_L=log_f_ab(sum_log_lambda, sum_lambda,L,b,n);

		J--;
	}

	while((K>0) & (z<log_f_R))
	{
		R=R+w;
		log_f_R=log_f_ab(sum_log_lambda, sum_lambda,R,b,n);

		K--;
	}

/** Make sure it is in the range **/
	L=fmax2(0.,L);
	R=fmin2(10000.,R);

	a1=runif(L,R);
	while(log_f_ab(sum_log_lambda, sum_lambda, a1, b, n)<z)
	{
	/** Shrink the interval **/
		if(a1<a0)
			L=a1;
		else
			R=a1;
		a1=runif(L,R);
	}

	return(a1);

}

double slice_sampling_rho(double rho, double w, int p, double SSR1, double SSR2, double SS12, int n)
{

	double L, R, z, rho_new;
	int J,K;
	double log_f_R=0,log_f_L=0;



	z=log_f_rho(SSR1, SSR2, SS12, rho, n)-rgamma(1,1);

/** lower bound **/
	L=rho-w*runif(0,1);
	R=L+w;
	J=(int)(p*runif(0,1));
	K=(p-1)-J;


/** Find the interval **/
	log_f_L=log_f_rho(SSR1, SSR2, SS12, L, n);
	log_f_R=log_f_rho(SSR1, SSR2, SS12, R, n);

	while((J>0) & (z<log_f_L))
	{
		L=L-w;
		log_f_L=log_f_rho(SSR1, SSR2, SS12,L,n);

		J--;

	}

	while((K>0) & (z<log_f_R))
	{
		R=R+w;
		log_f_R=log_f_rho(SSR1, SSR2, SS12,R,n);

		K--;
	}


/** Make sure it is in the range **/
	L=fmax2(-1.,L);
	R=fmin2(1.,R);


	rho_new=runif(L,R);

	while(log_f_rho(SSR1, SSR2, SS12, rho_new, n)<z)
	{
	/** Shrink the interval **/
		if(rho_new<rho)
			L=rho_new;
		else
			R=rho_new;
		rho_new=runif(L,R);
	}

	return(rho_new);

}


double log_f_rho(double SSR1, double SSR2, double SS12, double rho, int n)
{

	double log_f;

	log_f=-n/2.*log(1-rho*rho)-1./(2.*(1-rho*rho))*(SSR1-2*rho*SS12+SSR2)+dunif(rho,-1,1,1);

	return(log_f);

}


double log_f_p(int sum_gamma_equ, int sum_gamma_diff, double prob)
{

	double log_f;

	log_f=sum_gamma_equ*log(1.-prob)+sum_gamma_diff*log(prob)+dunif(prob,0.,1.,1);

	return(log_f);

}


double slice_sampling_p(double prob0, double w, int p, int sum_gamma_equ, int sum_gamma_diff)
{

	double L, R, z, prob1;
	int J,K;
	double log_f_R=0,log_f_L=0;


	z=log_f_p(sum_gamma_equ, sum_gamma_diff, prob0)-rgamma(1,1);

/** lower bound **/
	L=prob0-w*runif(0,1);
	R=L+w;
	J=(int)(p*runif(0,1));
	K=(p-1)-J;


/** Find the interval **/
	log_f_R=log_f_p(sum_gamma_equ, sum_gamma_diff, R);
	log_f_L=log_f_p(sum_gamma_equ, sum_gamma_diff, L);

/** Find the interval **/
	while((J>0) & (z<log_f_L))
	{
		L=L-w;
		log_f_L=log_f_p(sum_gamma_equ, sum_gamma_diff, L);

		J--;
	}

	while((K>0) & (z<log_f_R))
	{
		R=R+w;
		log_f_R=log_f_p(sum_gamma_equ, sum_gamma_diff, R);

		K--;
	}


/** Make sure it is in the range **/
	L=fmax2(0.,L);
	R=fmin2(1.,R);


	prob1=runif(L,R);
	while(log_f_p(sum_gamma_equ, sum_gamma_diff, prob1)<z)
	{
		if(prob1<prob0)
			L=prob1;
		else
			R=prob1;
		prob1=runif(L,R);
	}

	return(prob1);

}

