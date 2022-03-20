/*****************************************************************************/
/* splinefit.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Functions related to cubic spline fitting in 1 and 2 dimension.	     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2006; Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "point.h"
#include "fit/lmfit.h"
#include "splinefit.h"

#define		RECIP6		(1.0/6.0)

/*****************************************************************************/

static int basespline_coeff(int n,int c,double *y2,double *tmp)
{
 int	i;
 double	p,qn,un,*u;

 u=tmp;

 y2[0]=u[0]=0.0;
 for ( i=1 ; i<n-1 ; i++ )
  {	p=0.5*y2[i-1]+2.0;
	y2[i]=-0.5/p;
	switch ( i-c )
	 {	case -1: u[i]=+1.0 ;break;
	     	case  0: u[i]=-2.0;break;
		case +1: u[i]=+1.0 ;break;
		default: u[i]= 0.0 ;break;
	 }	
	u[i]=(3.0*u[i]-0.5*u[i-1])/p;
  }
 qn=un=0.0;
 y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
 for ( i=n-2 ; i>=0 ; i-- )
  {	y2[i]=y2[i]*y2[i+1]+u[i];					}

 return(0);
}

/*****************************************************************************/

int fit_1d_spline(double x0,double x1,int nx,
	double *xv,double *yv,int npoint,double *wv,double *coeff)
{
 int	nvar,i,j,k,t;
 double	**splinebase,*tmp;
 double	**amatrix,*bvector,*fvars,x,y,a,b,w;
 double	yl,yr,y2l,y2r,c2l,c2r;

 nvar=nx+1;

 if ( x0==x1 )		return(-1);
 else if ( x0>x1 )	x=x0,x0=x1,x1=x;

 if ( xv==NULL || yv==NULL )
	return(-1);
 if ( npoint<nvar )
	return(1);

 splinebase=matrix_alloc(nvar);
 tmp=vector_alloc(nvar);
 for ( i=0 ; i<=nx ; i++ )
  {	basespline_coeff(nx+1,i,splinebase[i],tmp);	}
 vector_free(tmp);

 amatrix=matrix_alloc(nvar);
 bvector=vector_alloc(nvar);
 fvars=vector_alloc(nvar);

 for ( i=0 ; i<nvar ; i++ )
  {	for ( j=0 ; j<nvar ; j++ )
	 {	amatrix[i][j]=0.0;		}
	bvector[i]=0.0;
  }

 for ( k=0 ; k<npoint ; k++ )
  {	x=xv[k],
	y=yv[k];
	if ( x<x0 || x>=x1 )	continue;
	t=(int)floor((double)nx*(x-x0)/(x1-x0));
	b=(double)nx*(x-x0)/(x1-x0)-(double)t;
	a=1.0-b;
	c2l=(a*a*a-a)*RECIP6;
	c2r=(b*b*b-b)*RECIP6;
	for ( i=0 ; i<=nx ; i++ )
	 {	if ( i==t+0 )		yl=1.0,yr=0.0;
		else if ( i==t+1 )	yl=0.0,yr=1.0;
		else			yl=0.0,yr=0.0;
		y2l=splinebase[i][t+0];
		y2r=splinebase[i][t+1];
		fvars[i]=a*yl+b*yr+c2l*y2l+c2r*y2r;
	 }
	if ( wv == NULL )
		w=1.0;
	else
		w=wv[k];

	for ( i=0 ;i<nvar ; i++ )
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	amatrix[i][j]+=w*fvars[i]*fvars[j];		}
		bvector[i]+=w*y*fvars[i];
	 }
  }

 t=solve_gauss(amatrix,bvector,nvar);
 if ( ! t )
  {	for ( i=0 ; i<nvar ; i++ )
	 {	coeff[i]=bvector[i];		}
  }

 vector_free(fvars);
 vector_free(bvector);
 matrix_free(amatrix);
 matrix_free(splinebase);

 if ( t )	return(1); 
 else		return(0);
}

/*****************************************************************************/

int fit_2d_spline(double x0,double x1,int nx,double y0,double y1,int ny,
	point *points,int npoint,double **coeff)
{
 int	nvar,mx,t,tx,ty,i,j,k;
 double	x,y,z,w;
 double	**splinexbase,**splineybase,*tmp;
 double	**amatrix,*bvector,*fvars,*zvsx,*zvsy;
 double	ax,ay,bx,by,rx,ry;
 double	zlx,zrx,z2lx,z2rx,
	zly,zry,z2ly,z2ry,
	zix,ziy;
 double	c2lx,c2rx,c2ly,c2ry;

 if ( x0==x1 )		return(-1);
 else if ( x0>x1 )	x=x0,x0=x1,x1=x;
 if ( y0==y1 )		return(-1);
 else if ( y0>y1 )	y=y0,y0=y1,y1=y;
 if ( points==NULL )	return(-1);

 nvar=(nx+1)*(ny+1);
 if ( npoint<nvar )	return(1);
 mx=(nx>ny?nx:ny)+1;

 splinexbase=matrix_alloc(nx+1);
 splineybase=matrix_alloc(ny+1);
 tmp=vector_alloc(mx);
 for ( i=0 ; i<=nx ; i++ )
  {	basespline_coeff(nx+1,i,splinexbase[i],tmp);	}
 for ( i=0 ; i<=ny ; i++ )
  {	basespline_coeff(ny+1,i,splineybase[i],tmp);	}
 vector_free(tmp);

 amatrix=matrix_alloc(nvar);
 bvector=vector_alloc(nvar);
 fvars  =vector_alloc(nvar);
 zvsx=(double *)malloc(sizeof(double)*(nx+1));
 zvsy=(double *)malloc(sizeof(double)*(ny+1));

 for ( i=0 ; i<nvar ; i++ )
  {	for ( j=0 ; j<nvar ; j++ )
	 {	amatrix[i][j]=0.0;	 }
	bvector[i]=0.0;
  }

 for ( k=0 ; k<npoint ; k++ )
  {	x=points[k].x,
	y=points[k].y,
	z=points[k].value;
	w=points[k].weight;

	rx=(double)nx*(x-x0)/(x1-x0);
	tx=(int)floor(rx),bx=rx-(double)tx,ax=1.0-bx;
	ry=(double)ny*(y-y0)/(y1-y0);
	ty=(int)floor(ry),by=ry-(double)ty,ay=1.0-by;

	c2lx=(ax*ax*ax-ax)*RECIP6,c2rx=(bx*bx*bx-bx)*RECIP6;
	c2ly=(ay*ay*ay-ay)*RECIP6,c2ry=(by*by*by-by)*RECIP6;

	for ( i=0 ; i<=ny ; i++ )
	 {	if ( i==ty+0 )		zly=1.0,zry=0.0;
		else if ( i==ty+1 )	zly=0.0,zry=1.0;
		else			zly=0.0,zry=0.0;
		z2ly=splineybase[i][ty+0];
		z2ry=splineybase[i][ty+1];
		zvsy[i]=ay*zly+by*zry+c2ly*z2ly+c2ry*z2ry;
	 }
	for ( j=0 ; j<=nx ; j++ )
	 {	if ( j==tx+0 )		zlx=1.0,zrx=0.0;
		else if ( j==tx+1 )	zlx=0.0,zrx=1.0;
		else			zlx=0.0,zrx=0.0;
		z2lx=splinexbase[j][tx+0];
		z2rx=splinexbase[j][tx+1];
		zvsx[j]=ax*zlx+bx*zrx+c2lx*z2lx+c2rx*z2rx;
	 }
	for ( i=0 ; i<=ny ; i++ )
	 {	for ( j=0 ; j<=nx ; j++ )
		 {	zix=zvsx[j];
			ziy=zvsy[i];
			fvars[i*(nx+1)+j]=zix*ziy;
		 }
	 }
	for ( i=0 ; i<nvar ; i++ )
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	amatrix[i][j]+=w*fvars[i]*fvars[j];	}
		bvector[i]+=w*fvars[i]*z;
	 }
  }
 t=solve_gauss(amatrix,bvector,nvar);
 if ( ! t )
  {	for ( i=0,k=0 ; i<=ny ; i++ )
	 {	for ( j=0 ; j<=nx ; j++,k++ )
		 {	coeff[i][j]=bvector[k];		}
	 }
  }

 free(zvsy);
 free(zvsx);
 vector_free(fvars);
 vector_free(bvector);
 matrix_free(amatrix);
 
 matrix_free(splineybase);
 matrix_free(splinexbase);
 
 if ( t )	return(1);
 else		return(0);
}

/*****************************************************************************/
                      
