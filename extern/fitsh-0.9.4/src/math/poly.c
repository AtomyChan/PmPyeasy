/*****************************************************************************/
/* poly.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Standalone library for evaluating polynomials (standard and Legendre)     */
/* defined in R^2 (in a plane).						     */
/* (c) 2001, 2004, 2005, 2006; Pal, A. (apal@szofi.elte.hu). 		     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in poly.h	     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "poly.h"

/*****************************************************************************/

int eval_2d_monoms(double x,double y,int order,double *ret,double ox,double oy,double scale)
{
 int	i,j,k,n;

 x=(x-ox)/scale;
 y=(y-oy)/scale;

 switch ( order )
  {  case 0:
	ret[0]=1.0;
	break;
     case 1:
	ret[0]=1.0;
	ret[1]=x;ret[2]=y;
	break;
     default:
	ret[0]=1.0;
	for ( i=0,k=0 ; k<order ; k++ )
	 {	for ( j=0,n=k+1 ; j<=k ; j++,n-- )
			ret[i+k+1+j]=ret[i+j]*x/(double)n;
		ret[i+2*k+2]=ret[i+k]*y/(double)(k+1);
		i+=k+1;
	 };
	break;
  }
 return(0);
}

double eval_2d_poly(double x,double y,int order,double *coeff,double ox,double oy,double scale)
{
 int	m,n,i,k,p,yf;
 double	ym,ret,w;

 x=(x-ox)/scale;
 y=(y-oy)/scale;

 switch ( order )
  {  case 0:
	return(coeff[0]);
     case 1:
	return(coeff[0]+x*coeff[1]+y*coeff[2]);
     default:
	ret=0.0;
	m=(order)*(order+1)/2;n=order;ym=1.0;yf=1;
	while ( n>=0 )
	 {	i=n;p=m;w=coeff[m];k=order;
		while ( i>=1 )
		 {	w=w*x/(double)i;
			p-=k;
			w+=coeff[p];
			i--;k--;
		 };
		ret+=w*ym;
		ym=ym*y/(double)yf;
		yf++;
		n--,m++;
	 };
	return(ret);
  }
}

/*****************************************************************************/

int eval_2d_leg_monoms(double x,double y,int order,double *ret,double ox,double oy,double scale)
{
 int	i,j,k,n,p1,p2;

 x=(x-ox)/scale;
 y=(y-oy)/scale;

 switch ( order )
  {  case 0:
	ret[0]=1.0;
	break;
     case 1:
	ret[0]=1.0;
	ret[1]=x,ret[2]=y;
	break;
     default:
	ret[0]=1.0;
	ret[1]=x,ret[2]=y;
	for ( k=1,p1=0,p2=1,n=3 ; k<order ; k++,p1=p2,p2=n,n+=k+1 )
	 {	ret[n]=((2*k+1)*x*ret[p2]-k*ret[p1])/(double)(k+1);	}
	for ( k=1,p1=0,p2=2,n=5 ; k<order ; k++,p1=p2,p2=n,n+=k+2 )
	 {	ret[n]=((2*k+1)*y*ret[p2]-k*ret[p1])/(double)(k+1);	}
	for ( k=2,p1=1,n=3 ; k<=order ; k++ )
	 {	n++;
		for ( i=p1,p2=2,j=1 ; j<k ; j++ )
		 {	ret[n]=ret[i]*ret[p2];
			n++;p2+=j+2;i-=k-j;
		 }
		n++;
		p1+=k;
	 }
	break;
  }	
 return(0);
}

double eval_2d_leg_poly(double x,double y,int order,double *coeff,double ox,double oy,double scale)
{
 int	i,n,j,m;
 double	ret,xlp,xl,xln,ylp,yl,yln;

 x=(x-ox)/scale;
 y=(y-oy)/scale;

 switch ( order )
  {  case 0:
	return(coeff[0]);
     case 1:
	return(coeff[0]+x*coeff[1]+y*coeff[2]);
     default:
	ret=0.0;
	ylp=1.0,yl=y;m=0;
	for ( j=0 ; j<=order ; )
	 {	ret+=coeff[m]*ylp;
		n=m+j+1;
		xlp=ylp,xl=x*ylp;
		for ( i=1 ; i<=order-j ; i++ )
		 {	ret+=coeff[n]*xl;
			if ( i<order-j )
			 {	xln=((2*i+1)*x*xl-i*xlp)/(double)(i+1);
				xlp=xl,xl=xln;
				n+=j+i+1;
			 }
		 }
		j++;
		m+=j+1;
		yln=((2*j+1)*y*yl-j*ylp)/(double)(j+1);
		ylp=yl,yl=yln;
	 }
	return(ret);
  }
}

/*****************************************************************************/

int diff_2d_poly(double *coeff,int order,double *dxcoeff,double *dycoeff)
{
 int	o,i,k,kx,ky;
 k=kx=ky=0;
 for ( o=0 ; o<=order ; o++ )
  {	for ( i=0 ; i<=o ; i++,k++ )
	 {	if ( dxcoeff != NULL && i<o )	dxcoeff[kx]=coeff[k],kx++;
		if ( dycoeff != NULL && 0<i )	dycoeff[ky]=coeff[k],ky++;
	 }
  }
 return(0);
}

/*****************************************************************************/

double calc_2d_unitarity(double *xfit,double *yfit,int order)
{
 double ma,mb,mc,md;
 double n1,n2,nn,dd,det,dx,dy;
 double unitarity;
 int	nvar;
 double	*dxfit,*dyfit,*tarr;

 if ( order>=2 )
  {	ma=xfit[1],mb=xfit[2];
	mc=yfit[1],md=yfit[2];
	det=ma*md-mb*mc;
	dx=(-xfit[0]*md+yfit[0]*mb)/det;
	dy=(-ma*yfit[0]+mc*xfit[0])/det;
	if ( det==0.0 )
	 {	return(-1.0);		}
	nvar=(order)*(order+1)/2;
	tarr=(double *)malloc(sizeof(double)*nvar*2);
	dxfit=tarr,dyfit=tarr+nvar;
	diff_2d_poly(xfit,order,dxfit,dyfit);
	ma=eval_2d_poly(dx,dy,order-1,dxfit,0,0,1);
	mb=eval_2d_poly(dx,dy,order-1,dyfit,0,0,1);
	diff_2d_poly(yfit,order,dxfit,dyfit);
	mc=eval_2d_poly(dx,dy,order-1,dxfit,0,0,1);
	md=eval_2d_poly(dx,dy,order-1,dyfit,0,0,1);
	free(tarr);
  }
 else
  {	ma=xfit[1],mb=xfit[2];
	mc=yfit[1],md=yfit[2];
  }

 n1=(ma-md)*(ma-md)+(mb+mc)*(mb+mc);
 n2=(ma+md)*(ma+md)+(mb-mc)*(mb-mc);
 nn=(n1<n2?n1:n2);
 dd=ma*ma+mb*mb+mc*mc+md*md;

 if ( dd<=0.0 )
  {	return(-1.0);		}
 else
  {	if ( nn<=0.0 )	unitarity=0.0;
	else		unitarity=sqrt(nn/dd);
	return(unitarity);
  }
}

/*****************************************************************************/

static int multiply_2d_poly_with_affine(double *coeff,int order,double *lin)
{
 int	o,n,i;

 n=(order)*(order+1)/2;
 for ( o=order ; o>=0 ; o-- )
  {	for ( i=0 ; i<=o ; i++ )
	 {	coeff[n+o+1+i]+=lin[1]*coeff[n+i];
		coeff[n+o+2+i]+=lin[2]*coeff[n+i];
		coeff[n+i]*=lin[0];
	 }
	n-=o;
  }

 return(0);
}

int compose_2d_poly_with_affine(double *coeff,int order,double *xlin,double *ylin,double *rc)
{
 double	*xm,*ym,c,m,f;
 int	nvar,i,o,j,n;

 if ( order<0 || coeff==NULL || xlin==NULL || ylin==NULL || rc==NULL )
	return(-1);
 else if ( order==0 )
  {	rc[0]=coeff[0];
	return(0);
  }

 nvar=(order+1)*(order+2)/2;
 xm=(double *)malloc(sizeof(double)*nvar);
 ym=(double *)malloc(sizeof(double)*nvar);

 xm[0]=ym[0]=1.0,rc[0]=0.0;
 for ( i=1 ; i<nvar ; i++ )
  {	xm[i]=0.0;
	ym[i]=0.0;
	rc[i]=0.0;
  }

 f=1.0;
 for ( o=0 ; o<=order ; o++ )
  {	for ( i=0 ; i<(o+1)*(o+2)/2 ; i++ )
	 {	xm[i]=ym[i];			}
	for ( ; i<nvar ; i++ )
	 {	xm[i]=0;			}

	m=f;
	for ( j=0 ; o+j<=order ; j++ )
	 {	n=(o+j+1)*(o+j+2)/2;
		c=coeff[(o+j)*(o+j+1)/2+o]/m;
		for ( i=0 ; i<n ; i++ )
		 {	rc[i]+=xm[i]*c;		}

		if ( o+j<order )
			multiply_2d_poly_with_affine(xm,o+j,xlin);
		m=m*(double)(j+1);
	 }
	if ( o<order )
		multiply_2d_poly_with_affine(ym,o,ylin);
	f=f*(double)(o+1);
  }

 f=1.0;
 for ( o=0,i=0 ; o<=order ; o++ )
  {	m=f;
	for ( j=0 ; o+j<=order ; j++ )
	 {	rc[(o+j)*(o+j+1)/2+o]*=m;
		m=m*(double)(j+1);
	 }
	f=f*(double)(o+1);
  }
 free(ym);
 free(xm);

 return(0);
}

/*****************************************************************************/
                    
