/*****************************************************************************/
/* bicubic.c 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Standalone library for 2D cubic spline interpolation.		     */
/* (c) 2004, Pal, A. (apal@szofi.elte.hu)				     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in bicubic.h	     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "bicubic.h"

/*****************************************************************************/

static int cubic_spline_coeff(double **y,int px,int py,int dx,int dy,int n,double **y2,int lx,int ly,double *u)
{
 int	i,k;
 double	p,qn,un,yprev,ycurr,ynext;

 if ( n<=0 )
	return(-1);
 else if ( n<2 )
  { 	y2[ly][lx]=0.0;
	return(0);
  }

 y2[ly][lx]=u[0]=0.0;
 yprev=y[py][px];
 py+=dy,px+=dx;
 ycurr=y[py][px];
 for ( i=1 ; i<n-1 ; i++ )
  {	p=0.5*y2[ly][lx]+2.0;
	ly+=dy,lx+=dx;
	y2[ly][lx]=-0.5/p;
	py+=dy,px+=dx;
	ynext=y[py][px];
	u[i]=ynext-2.0*ycurr+yprev;
	u[i]=(3.0*u[i]-0.5*u[i-1])/p;
	yprev=ycurr,
	ycurr=ynext;
  }
 qn=un=0.0;
 ynext=y2[ly+dy][lx+dx]=(un-qn*u[n-2])/(qn*y2[ly][lx]+1.0);
 for ( k=n-2 ; k>=0 ; k-- )
  {	ycurr=y2[ly][lx]=y2[ly][lx]*ynext+u[k];
	ly-=dy,lx-=dx;
	ynext=ycurr;
  }

 return(0);
}

int bicubic_coeff(double **y,int sx,int sy,double **c,char **mask)
{
 int	i,j,ms;
 double	*tmp;

 ms=(sx>sy?sx:sy);
 tmp=(double *)malloc(sizeof(double)*ms);

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	c[2*i][2*j]=y[i][j];		}
  }

 if ( mask==NULL )
  {	for ( i=0 ; i<sy ; i++ )
	 {	cubic_spline_coeff(c,0,2*i,2,0,sx,c,1,2*i,tmp);		}
	for ( i=0 ; i<2*sx ; i++ )
	 {	cubic_spline_coeff(c,i,0  ,0,2,sy,c,i,1  ,tmp);		}
  }
 else
  {	int	i0,j0;
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; )
		 {	while ( j<sx &&   mask[i][j] )	j++;
			j0=j;
			while ( j<sx && ! mask[i][j] )	j++;
			if ( j>j0 )
			 {	
				cubic_spline_coeff(c,2*j0,2*i,2,0,j-j0,c,2*j0+1,2*i,tmp);
		 	 }
		 }
	 }
	for ( j=0 ; j<sx ; j++ )
	 {	for ( i=0 ; i<sy ; )
		 {	while ( i<sy &&   mask[i][j] )	i++;
			i0=i;
			while ( i<sy && ! mask[i][j] )	i++;
			if ( i>i0 )
			 {	cubic_spline_coeff(c,2*j+0,2*i0,0,2,i-i0,c,2*j+0,2*i0+1,tmp);
				cubic_spline_coeff(c,2*j+1,2*i0,0,2,i-i0,c,2*j+1,2*i0+1,tmp);
			 }
		 }
	 }
  }

 free(tmp);

 return(0);
}

double bicubic_inter(double **c,double x,double y)
{
 int	ix,iy;
 double	ax,bx,ay,by,ax3,bx3;
 double	zb,zt,z2b,z2t,z;

 ix=(int)x,x-=(double)ix;ix=ix*2;
 iy=(int)y,y-=(double)iy;iy=iy*2;

 ay=(1.0-y),by=y;
 ax=(1.0-x),bx=x;
 ax3=(ax*ax*ax-ax)/6.0,
 bx3=(bx*bx*bx-bx)/6.0;
 
 zb =ax*c[iy  ][ix]+bx*c[iy  ][ix+2]+ax3*c[iy  ][ix+1]+bx3*c[iy  ][ix+3];
 z2b=ax*c[iy+1][ix]+bx*c[iy+1][ix+2]+ax3*c[iy+1][ix+1]+bx3*c[iy+1][ix+3];
 zt =ax*c[iy+2][ix]+bx*c[iy+2][ix+2]+ax3*c[iy+2][ix+1]+bx3*c[iy+2][ix+3];
 z2t=ax*c[iy+3][ix]+bx*c[iy+3][ix+2]+ax3*c[iy+3][ix+1]+bx3*c[iy+3][ix+3];
 z=ay*zb+by*zt+((ay*ay*ay-ay)*z2b+(by*by*by-by)*z2t)/6.0;

 return(z);
}

/*****************************************************************************/
                                                                    
                 
