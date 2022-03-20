/*****************************************************************************/
/* lmfit.c 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Standalone library for weighted linear and nonlinear fit of real	     */
/* data with a function defined on an arbitrary domain.			     */
/* (c) 2003, 2004, 2006; Pal, A. (apal@szofi.elte.hu). Some functions are    */
/* based on the appropriate chapters of Numerical Recipes (in C)             */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in lmfit.h         */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "lmfit.h"

/*****************************************************************************/

double **matrix_alloc_gen(int x,int y)
{
 int	i;
 double **m;
 m=(double **)malloc(sizeof(double *)*(y+1));
 if ( m==NULL )	return(NULL);
 for ( i=0 ; i<y ; i++ )
  {	m[i]=(double *)malloc(sizeof(double)*x);
	if ( m[i]==NULL )
	 {	matrix_free(m);
		return(NULL);
	 }
  }
 m[y]=NULL;
 return(m);
}

double **matrix_alloc(int d)
{
 double **m;
 m=matrix_alloc_gen(d,d);
 return(m);
}

void matrix_free(double **m)
{
 int	i;
 for ( i=0 ; m[i] != NULL ; i++ )
  {	free(m[i]);		}
 free(m);
}
double *vector_alloc(int d)
{
 double *v;
 v=(double *)malloc(sizeof(double)*d);
 return(v);
}
void vector_free(double *v)
{
 free(v);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static void swap_matrix_rows(double **a,int r,int s)
{
 double *t;
 if ( r==s )	return;
 t=a[r],a[r]=a[s],a[s]=t;
}
static void swap_vector_elements(double *v,int r,int s)
{
 double t;
 if ( r==s )	return;
 t=v[r],v[r]=v[s],v[s]=t;
}
static double normalize_vector(double *v,int n)
{
 int	i;
 double	mm;

 mm=0;
 for ( i=0 ; i<n ; i++ )
  {	if ( -v[i] > mm )	mm=-v[i];
	else if ( v[i] > mm )	mm=v[i];
  }
 if ( mm==0.0 )	return(mm);

 for ( i=0 ; i<n ; i++ )
  {	v[i]/=mm;				}
 return(mm);
}

int solve_gauss(double **a,double *b,int n)
{
 int	i,j,k;
 double	t;

 for ( i=0 ; i<n ; i++ )
  {	t=normalize_vector(a[i],n);
	if ( t==0 )	return(-1);
	else		b[i]/=t;
  }
 for ( i=0 ; i<n ; i++ )
  {	for ( j=i ; j<n ; j++ )
	 {	if ( a[j][i] != 0 )	break;		}
	if ( j==n )	return(-1);
	swap_matrix_rows(a,i,j);
	swap_vector_elements(b,i,j);
	for ( j=i+1 ; j<n ; j++ )
	 {	t=a[j][i]/a[i][i];
		for ( k=0 ; k<n ; k++ )
		 {	a[j][k]-=t*a[i][k];		}
		b[j]-=t*b[i];
	 }
  }
 for ( i=n-1 ; i>=0 ; i-- )
  {	t=b[i];
	for ( j=i+1 ; j<n ; j++ )
	 {	t-=b[j]*a[i][j];		}
	t/=a[i][i];
	b[i]=t;
  }
 return(0);
}

int invert_gauss(double **a,int n)
{
 int	i,j,k;
 double	t,**e;

 e=matrix_alloc(n);

 for ( i=0 ; i<n ; i++ )
  {	for ( j=0 ; j<n ; j++ )
	 {	e[i][j]=0.0;		}
	e[i][i]=1.0;
  }
 for ( i=0 ; i<n ; i++ )
  {	t=normalize_vector(a[i],n);
	if ( t==0 )	return(-1);
	else 
	 {	for ( j=0 ; j<n ; j++ )
		 {	e[i][j]/=t;		}
	 }
  }
 for ( i=0 ; i<n ; i++ )
  {	for ( j=i ; j<n ; j++ )
	 {	if ( a[j][i] != 0 )	break;		}
	if ( j==n )	return(-1);
	swap_matrix_rows(a,i,j);
	swap_matrix_rows(e,i,j);
	for ( j=i+1 ; j<n ; j++ )
	 {	t=a[j][i]/a[i][i];
		for ( k=0 ; k<n ; k++ )
		 {	a[j][k]-=t*a[i][k];
			e[j][k]-=t*e[i][k];
		 }
	 }
  }
 for ( i=n-1 ; i>=0 ; i-- )
  {	for ( k=0 ; k<n ; k++ )
	 {	t=e[i][k];
		for ( j=i+1 ; j<n ; j++ )
		 {	t-=e[j][k]*a[i][j];		}
		t/=a[i][i];
		e[i][k]=t;
	 }
  }

 for ( i=0 ; i<n ; i++ )
  {	for ( j=0 ; j<n ; j++ )
	 {	a[i][j]=e[i][j];		}
  }

 matrix_free(e);
 
 return(0);
}

/*****************************************************************************/

double nlm_fit_base_con(void **x,double *y,double *a,double *weight,
	void (*funct)(void *,double *,double *,double *,void *),
	int n,int ndata,void *param,constraint *cc,
	double lambda,double lambda_mpy)
{
 int	i,j,k,nc;
 double *da,**amatrix,*bvector,yy,chi2old,chi2new,ww;

 if ( cc==NULL || cc->nc <= 0 )	nc=0;
 else				nc=cc->nc;

 da=bvector=NULL;amatrix=NULL;
 da     =vector_alloc(n);
 bvector=vector_alloc(n+nc);
 amatrix=matrix_alloc(n+nc);
 if ( da==NULL || bvector==NULL || amatrix==NULL )
  {	if ( da      != NULL )	vector_free(da);
	if ( bvector != NULL )	vector_free(bvector);
	if ( amatrix != NULL )	matrix_free(amatrix);
	return(-1.0);
  }
 
 chi2old=0.0;
 for ( i=0 ; i<n+nc ; i++ )
  {	for ( j=0 ; j<n+nc ; j++ )
	 {	amatrix[i][j]=0.0;		}
	bvector[i]=0.0;
  }

 for ( k=0 ; k<ndata ; k++ )
  {	funct(x[k],a,&yy,da,param);
	if ( weight != NULL )	ww=weight[k];
	else			ww=1.0;
	yy=y[k]-yy;
	chi2old+=yy*yy*ww;
	for ( i=0 ; i<n ; i++ )
	 {	for ( j=0 ; j<=i ; j++ )
		 {	amatrix[i][j]+=da[i]*da[j]*ww;	}
		bvector[i]+=yy*da[i]*ww;
	 }
  }

 for ( i=0 ; i<n ; i++ )
  {	amatrix[i][i]*=(1.0+lambda);	}

 for ( i=0 ; i<n+nc ;  i++ )
  {	for ( j=i+1 ; j<n+nc ; j++ )
	 {	amatrix[i][j]=amatrix[j][i];		}
  }

 if ( nc>0 && cc != NULL )
  {	for ( i=0 ; i<nc ; i++ )
	 {	bvector[i+n]=cc->cvector[i];
		for ( j=0 ; j<n ; j++ )
		 {	amatrix[i+n][j]=cc->cmatrix[i][j];
			amatrix[j][i+n]=cc->cmatrix[i][j];
			bvector[i+n]-=cc->cmatrix[i][j]*a[j];
		 }
	 }
  }


 i=solve_gauss(amatrix,bvector,n+nc);
 if ( i )
  {	vector_free(da);
	vector_free(bvector);
	matrix_free(amatrix);
	return(0.0);
  }

 for ( i=0 ; i<n ; i++ )
  {	bvector[i]+=a[i];		}

 chi2new=0.0;
 for ( k=0 ; k<ndata ; k++ )
  {	funct(x[k],bvector,&yy,NULL,param);
	if ( weight != NULL )	ww=weight[k];
	else			ww=1.0;
	yy=y[k]-yy;
	chi2new+=yy*yy*ww;
  }

 if ( chi2new>=chi2old )
  {	lambda=lambda*lambda_mpy;	}
 else
  {	for ( i=0 ; i<n ; i++ )
	 {	a[i]=bvector[i];		}
	lambda=lambda/lambda_mpy;
  }

 vector_free(da);
 vector_free(bvector);
 matrix_free(amatrix);

 return(lambda);
}

double nlm_fit_base(void **x,double *y,double *a,double *weight,
	void (*funct)(void *,double *,double *,double *,void *),
	int n,int ndata,void *param,double lambda,double lambda_mpy)
{
 double	ret;

 ret=nlm_fit_base_con(x,y,a,weight,funct,n,ndata,param,NULL,lambda,lambda_mpy);
 
 return(ret);
 
}

double nlm_fit_nmdf_con(void **x,double *y,double *a,double *weight,
	void (*funct)(void *,double *,double *,double *,void *),
	int n,int ndata,void *param,constraint *cc,
	double lambda,double lambda_mpy,double *dda)
{
 int	i,j,k,nc;
 double *da,*wa,**amatrix,*bvector,yy,y1,y2,chi2old,chi2new,ww,dd;

 if ( cc==NULL || cc->nc <= 0 )	nc=0;
 else				nc=cc->nc;

 da=bvector=NULL;amatrix=NULL;
 da     =vector_alloc(n);
 wa     =vector_alloc(n);
 bvector=vector_alloc(n+nc);
 amatrix=matrix_alloc(n+nc);
 if ( da==NULL || bvector==NULL || amatrix==NULL )
  {	if ( da      != NULL )	vector_free(da);
	if ( bvector != NULL )	vector_free(bvector);
	if ( amatrix != NULL )	matrix_free(amatrix);
	return(-1.0);
  }
 
 chi2old=0.0;
 for ( i=0 ; i<n+nc ; i++ )
  {	for ( j=0 ; j<n+nc ; j++ )
	 {	amatrix[i][j]=0.0;		}
	bvector[i]=0.0;
  }

 for ( k=0 ; k<ndata ; k++ )
  {	funct(x[k],a,&yy,NULL,param);

	for ( i=0 ; i<n ; i++ )
	 {	wa[i]=a[i];		}
	for ( i=0 ; i<n ; i++ )
	 {	dd=0.5*dda[i];
		wa[i]=a[i]+dd;
		funct(x[k],wa,&y2,NULL,param);
		wa[i]=a[i]-dd;
		funct(x[k],wa,&y1,NULL,param);
		wa[i]=a[i];
		da[i]=(y2-y1)/dda[i];
	 }

	if ( weight != NULL )	ww=weight[k];
	else			ww=1.0;
	yy=y[k]-yy;
	chi2old+=yy*yy*ww;
	for ( i=0 ; i<n ; i++ )
	 {	for ( j=0 ; j<=i ; j++ )
		 {	amatrix[i][j]+=da[i]*da[j]*ww;	}
		bvector[i]+=yy*da[i]*ww;
	 }
  }

 for ( i=0 ; i<n ; i++ )
  {	amatrix[i][i]*=(1.0+lambda);	}

 for ( i=0 ; i<n+nc ;  i++ )
  {	for ( j=i+1 ; j<n+nc ; j++ )
	 {	amatrix[i][j]=amatrix[j][i];		}
  }

 if ( nc>0 && cc != NULL )
  {	for ( i=0 ; i<nc ; i++ )
	 {	bvector[i+n]=cc->cvector[i];
		for ( j=0 ; j<n ; j++ )
		 {	amatrix[i+n][j]=cc->cmatrix[i][j];
			amatrix[j][i+n]=cc->cmatrix[i][j];
			bvector[i+n]-=cc->cmatrix[i][j]*a[j];
		 }
	 }
  }


 i=solve_gauss(amatrix,bvector,n+nc);
 if ( i )
  {	vector_free(da);
	vector_free(bvector);
	matrix_free(amatrix);
	return(0.0);
  }

 for ( i=0 ; i<n ; i++ )
  {	bvector[i]+=a[i];		}

 chi2new=0.0;
 for ( k=0 ; k<ndata ; k++ )
  {	funct(x[k],bvector,&yy,NULL,param);
	if ( weight != NULL )	ww=weight[k];
	else			ww=1.0;
	yy=y[k]-yy;
	chi2new+=yy*yy*ww;
  }

 if ( chi2new>=chi2old )
  {	lambda=lambda*lambda_mpy;	}
 else
  {	for ( i=0 ; i<n ; i++ )
	 {	a[i]=bvector[i];		}
	lambda=lambda/lambda_mpy;
  }

 matrix_free(amatrix);
 vector_free(bvector);
 vector_free(wa);
 vector_free(da);

 return(lambda);
}

double nlm_fit_nmdf(void **x,double *y,double *a,double *weight,
	void (*funct)(void *,double *,double *,double *,void *),
	int n,int ndata,void *param,
	double lambda,double lambda_mpy,double *dda)
{
 double	ret;

 ret=nlm_fit_nmdf_con(x,y,a,weight,funct,n,ndata,param,NULL,lambda,lambda_mpy,dda);
 
 return(ret);
}

/*****************************************************************************/

int lin_fit_con(void **x,double *y,double *a,double *weight,
	void (*funct)(void *,double *,double *,double *,void *),
	int n,int ndata,void *param,constraint *cc,double *err)
{
 int	i,j,k,nc;
 double *da,**amatrix,**ematrix,*bvector,yy,ww,ch2;
 
 if ( cc==NULL || cc->nc <= 0 )	nc=0;
 else				nc=cc->nc;

 da=bvector=NULL;
 amatrix=ematrix=NULL;
 da     =vector_alloc(n);
 bvector=vector_alloc(n+nc);
 amatrix=matrix_alloc(n+nc);
 ematrix=matrix_alloc(n);
 if ( da==NULL || bvector==NULL || amatrix==NULL || ematrix==NULL )
  {	if ( da != NULL )	vector_free(da);
	if ( bvector != NULL )	vector_free(bvector);
	if ( amatrix != NULL )	matrix_free(amatrix);
	if ( ematrix != NULL )	matrix_free(ematrix);
	return(-1);
  }
 
 for ( i=0 ; i<n+nc ; i++ )
  {	for ( j=0 ; j<n+nc ; j++ )
	 {	amatrix[i][j]=0.0;		}
	bvector[i]=0.0;
  }

 for ( k=0 ; k<ndata ; k++ )
  {	funct(x[k],a,&yy,da,param);
	if ( weight != NULL )	ww=weight[k];
	else			ww=1.0;
	for ( i=0 ; i<n ; i++ )
	 {	for ( j=0 ; j<=i ; j++ )
		 {	amatrix[i][j]+=da[i]*da[j]*ww;	}
		bvector[i]+=y[k]*da[i]*ww;
	 }
  }

 for ( i=0 ; i<n+nc ;  i++ )
  {	for ( j=i+1 ; j<n+nc ; j++ )
	 {	amatrix[i][j]=amatrix[j][i];		}
  }

 if ( nc>0 && cc != NULL )
  {	for ( i=0 ; i<nc ; i++ )
	 {	for ( j=0 ; j<n ; j++ )
		 {	amatrix[i+n][j]=cc->cmatrix[i][j];
			amatrix[j][i+n]=cc->cmatrix[i][j];
		 }
		bvector[i+n]=cc->cvector[i];
	 }
  }

 for ( i=0 ; i<n ; i++ )
  {	for ( j=0 ; j<n && err != NULL ; j++ )
	 {	ematrix[i][j]=amatrix[i][j];		}
  }

 j=solve_gauss(amatrix,bvector,n+nc);

 if ( ! j )
  {	for ( i=0 ; i<n ; i++ )
	 {	a[i]=bvector[i];		}
  }
 if ( ! j && err != NULL )
  {	ch2=0.0;
	for ( k=0 ; k<ndata ; k++ )	
	 {	funct(x[k],a,&yy,da,param);
		if ( weight != NULL )	ww=weight[k];
		else			ww=1.0;
		ch2+=(yy-y[k])*(yy-y[k])*ww;
	 }
	invert_gauss(ematrix,n);
	/*if ( weight != NULL )
	 {	for ( i=0 ; i<n ; i++ )
		 {	err[i]=sqrt(ematrix[i][i]);	}
	 }
	else*/ if ( ndata <= n )
	 {	for ( i=0 ; i<n ; i++ )
		 {	err[i]=0.0;			}
	 }
	else
	 {	ch2=sqrt(ch2/(double)(ndata-n));
		for ( i=0 ; i<n ; i++ )
		 {	err[i]=ch2*sqrt(ematrix[i][i]);	}
	 }
  }

 matrix_free(ematrix);
 matrix_free(amatrix);
 vector_free(bvector);
 vector_free(da);

 if ( j )	return(1);
 else		return(0);
}

int lin_fit(void **x,double *y,double *a,double *weight,
	void (*funct)(void *,double *,double *,double *,void *),
	int n,int ndata,void *param,double *err)
{
 int	ret;

 ret=lin_fit_con(x,y,a,weight,funct,n,ndata,param,NULL,err);

 return(ret);
}

/*****************************************************************************/
                   

