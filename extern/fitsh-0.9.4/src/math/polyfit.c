/*****************************************************************************/
/* polyfit.c [-> poly.c,lmfit.c,point.h]				     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* A library for fitting polynomials in R^2 (in a plane).		     */
/* The library is NOT standalone, it needs functions from lmfit.c and poly.c */
/* (c) 2001, 2004, 2005; Pal, A. (apal@szofi.elte.hu). 			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in polyfit.h       */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "fit/lmfit.h"
#include "poly.h"
#include "point.h"
#include "polyfit.h"

/*****************************************************************************/

static int fit_poly_alloc_arrays(int nvar,double **rfvars,double **rbvector,double ***ramatrix)
{
 int	i,j;

 *rfvars=vector_alloc(nvar);
 if ( *rfvars==NULL )
  {	return(1);		}
 *rbvector=vector_alloc(nvar);
 if ( *rbvector==NULL )
  {	vector_free(*rfvars);
	return(1); 
  }
 *ramatrix=matrix_alloc(nvar);
 if ( *ramatrix==NULL )
  {	vector_free(*rfvars);
	vector_free(*rbvector);
	return(1);
  }

 for ( i=0 ; i<nvar ; i++ )
  {	for ( j=0 ; j<nvar ; j++ )
	 {	(*ramatrix)[i][j]=0.0;		}
	(*rbvector)[i]=0.0;
  }

 return(0);
}

static int fit_poly_free_arrays(double *fvars,double *bvector,double **amatrix)
{
 if ( fvars != NULL )	vector_free(fvars);
 if ( bvector != NULL )	vector_free(bvector);
 if ( amatrix != NULL )	matrix_free(amatrix);
 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int fit_2d_poly(point *data,int nd,int order,double *fit,double ox,double oy,double scale)
{
 int	i,j,m,n,nvar;
 double	**amatrix,*bvector,*fvars;
 double	weight,f;

 nvar=(order+1)*(order+2)/2;

 if ( nd < nvar ) 	return(1);

 if ( fit_poly_alloc_arrays(nvar,&fvars,&bvector,&amatrix) )	return(-1);

 for ( i=0 ; i<nd ; i++ )
  {	eval_2d_monoms(data[i].x,data[i].y,order,fvars,ox,oy,scale);
	f=data[i].value,weight=data[i].weight;
	for ( m=0 ; m<nvar ; m++ )
	 {	bvector[m]+=f*fvars[m]*weight;		}
	for ( m=0 ; m<nvar ; m++ )
	 {	for ( n=0 ; n<nvar ; n++ )
		 {	amatrix[m][n]+=fvars[m]*fvars[n]*weight;	}
 	 }
  }
 i=solve_gauss(amatrix,bvector,nvar);
 for ( j=0 ; j<nvar && ! i ; j++ )
  {	fit[j]=bvector[j];		}
 fit_poly_free_arrays(fvars,bvector,amatrix);
 if ( i )	return(1);
 else		return(0);
}

int fit_2d_leg_poly(point *data,int nd,int order,double *fit,double ox,double oy,double scale)
{
 int	i,j,m,n,nvar;
 double	**amatrix,*bvector,*fvars;
 double	weight,f;

 nvar=(order+1)*(order+2)/2;

 if ( nd < nvar ) 	return(1);

 if ( fit_poly_alloc_arrays(nvar,&fvars,&bvector,&amatrix) )	return(-1);

 for ( i=0 ; i<nd ; i++ )
  {	eval_2d_leg_monoms(data[i].x,data[i].y,order,fvars,ox,oy,scale);
	f=data[i].value,weight=data[i].weight;
	for ( m=0 ; m<nvar ; m++ )
	 {	bvector[m]+=f*fvars[m]*weight;		}
	for ( m=0 ; m<nvar ; m++ )
	 {	for ( n=0 ; n<nvar ; n++ )
		 {	amatrix[m][n]+=fvars[m]*fvars[n]*weight;	}
 	 }
  }
 i=solve_gauss(amatrix,bvector,nvar);
 for ( j=0 ; j<nvar && ! i ; j++ )
  {	fit[j]=bvector[j];		}
 fit_poly_free_arrays(fvars,bvector,amatrix);
 if ( i )	return(1);
 else		return(0);
}

/*****************************************************************************/

int fit_1d_poly(int order,double *x,double *y,int nd,double *wv,double *coeff)
{
 int	i,j,m,n,nvar,o;
 double	**amatrix,*bvector,*fvars;
 double	weight,f;

 nvar=order+1;

 if ( nd < nvar ) 	return(1);

 if ( fit_poly_alloc_arrays(nvar,&fvars,&bvector,&amatrix) )	return(-1);

 for ( i=0 ; i<nd ; i++ )
  {	f=1.0;
	for ( o=0 ; o<=order ; o++ )
	 {	fvars[o]=f;
		f=f*x[i];
	 }
	f=y[i];
	if ( wv==NULL )
		weight=1.0;
	else
		weight=wv[i];

	for ( m=0 ; m<nvar ; m++ )
	 {	bvector[m]+=f*fvars[m]*weight;		}
	for ( m=0 ; m<nvar ; m++ )
	 {	for ( n=0 ; n<nvar ; n++ )
		 {	amatrix[m][n]+=fvars[m]*fvars[n]*weight;	}
 	 }
  }
 i=solve_gauss(amatrix,bvector,nvar);
 for ( j=0 ; j<nvar && ! i ; j++ )
  {	coeff[j]=bvector[j];		}
 fit_poly_free_arrays(fvars,bvector,amatrix);
 if ( i )	return(1);
 else		return(0);
}

/*****************************************************************************/
                                                                      
