/*****************************************************************************/
/* downhill.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Simple minimalization function based on the Nedler-Mead downhill simplex  */
/* algorithm and Numerical Recipes. The algorithm has been improved by the   */
/* author, namely the simplex can propagate among an 'nvar' dimensional      */
/* hypersurface of the 'ndim' dimensional vector space.			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2007; Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "downhill.h"

/*****************************************************************************/

static double downhill_simplex_amotry(double **p,double *y,double *psum,
	int ndim,int nvar,double (*funct)(void *param,double *arr),void *param,
	int ihi,double fac,double *ptry)
{
 int	j;
 double	fac1,fac2,ytry;

 fac1=(1.0-fac)/(double)nvar;
 fac2=fac1-fac;

 for ( j=0 ; j<ndim ; j++ )
  {	ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;	}

 ytry=funct(param,ptry);

 if ( ytry < y[ihi] )
  {	y[ihi]=ytry;
	for ( j=0 ; j<ndim ; j++ )
	 {	psum[j] += ptry[j]-p[ihi][j];
		p[ihi][j]=ptry[j];
	 }
  }

 return(ytry);
}

int downhill_simplex(double **p,int ndim,int nvar,
	double (*funct)(void *param,double *arr),void *param,
	double precision,int maxstep)
{
 double	*y,*psum,*ptry;
 int	i,ihi,ilo,inhi,j,mpts,nfunk;
 double	rtol,sum,swap,ysave,ytry;

 y   =(double *)malloc(sizeof(double)*(nvar+1));
 psum=(double *)malloc(sizeof(double)*ndim);
 ptry=(double *)malloc(sizeof(double)*ndim);

 for ( i=0 ; i<nvar+1 ; i++ )
	y[i]=funct(param,p[i]);

 mpts=nvar+1;

 for ( j=0 ; j<ndim ; j++ )
  {	for ( sum=0.0,i=0 ; i<mpts ; i++ )
	 {	sum += p[i][j];		}
	psum[j]=sum;
  }

 nfunk=0;

 if ( precision <= 0.0 )	precision=1e-11;

 while ( 1 )
  {	ilo=0;
	ihi=( y[0]>y[1] ? (inhi=1,0) : (inhi=0,1) );
	for ( i=0 ; i<mpts ; i++ )
	 {	if ( y[i] <= y[ilo] )
			ilo=i;
		if ( y[i] >  y[ihi] )
		 {	inhi=ihi;
			ihi=i;
		 }
		else if ( y[i] > y[inhi] && i != ihi )
			inhi=i;
	 }
	rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+1e-15);
	if ( rtol<precision )
	 {	swap=y[0],y[0]=y[ilo],y[ilo]=swap;
		for ( i=0 ; i<ndim ; i++ )
			swap=p[0][i],p[0][i]=p[ilo][i],p[ilo][i]=swap;
		break;
	 }
	if ( maxstep>0 && nfunk >= maxstep )
	 {	nfunk=-1;
		break;
	 }

	nfunk+=2;
	ytry=downhill_simplex_amotry(p,y,psum,ndim,nvar,funct,param,ihi,-1.0,ptry);
	if ( ytry <= y[ilo] )
	 {	ytry=downhill_simplex_amotry(p,y,psum,ndim,nvar,
			funct,param,ihi,2.0,ptry);
	 }
	else if ( ytry >= y[inhi] )
	 {	ysave=y[ihi];
		ytry=downhill_simplex_amotry(p,y,psum,ndim,nvar,
			funct,param,ihi,0.5,ptry);
		if ( ytry>=ysave )
		 {	for ( i=0 ; i<mpts ; i++ )
			 {	if ( i != ilo )
				 {	for ( j=0 ; j<ndim ; j++ )
					 {	p[i][j]=psum[j]=0.5*
						(p[i][j]+p[ilo][j]);
					 }
				 	y[i]=funct(param,psum);
				 }
		 	 }
			nfunk+=ndim;
			for ( j=0 ; j<ndim ; j++ )
			 {	for ( sum=0.0,i=0 ; i<mpts ; i++ )
				 {	sum += p[i][j];		}
				psum[j]=sum;
			 }
		 }
	 }
	else
		nfunk--;
  }

 free(ptry);
 free(psum);
 free(y);

 return(nfunk);
}

/*****************************************************************************/
