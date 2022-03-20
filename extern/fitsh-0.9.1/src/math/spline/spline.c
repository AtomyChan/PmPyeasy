/*****************************************************************************/
/* spline.c 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Standalone library for 1D spline interpolation.			     */
/* (c) 2004, Pal, A. (apal@szofi.elte.hu), these functions are based on      */
/* the book Numerical Reicpes (in C).		 			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in spline.h        */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>

#include "spline.h"

/*****************************************************************************/

int spline_coeff(double *x,double *y,int n,double *yp1,double *ypn,double *y2)
{
 int	i,k;
 double	p,qn,sig,un,*u;

 u=(double *)malloc(sizeof(double)*n);
 if ( yp1==NULL )
  {	y2[0]=u[0]=0.0;							}
 else
  {	y2[0]=-0.5;
	u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-(*yp1));
  }
 for ( i=1 ; i<n-1 ; i++ )
  {	sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
	p=sig*y2[i-1]+2.0;
	y2[i]=(sig-1.0)/p;
	u[i]=(y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
	u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
 if ( ypn==NULL )
  {	qn=un=0.0;							}
 else
  {	qn=0.5;
	un=(3.0/(x[n-1]-x[n-2]))*((*ypn)-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
 y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
 for ( k=n-2 ; k>=0 ; k-- )
  {	y2[k]=y2[k]*y2[k+1]+u[k];					}
 free(u);

 return(0);
}

double spline_inter(double *xa,double *ya,double *y2a,int n,double x)
{
 int	klo,khi,k;
 double	a,b,h,y;

 klo=0,khi=n-1;
 while ( khi-klo>1 )
  {	k=(khi+klo)/2;
	if ( xa[k] > x )	khi=k;
	else			klo=k;
  }
 h=xa[khi]-xa[klo];
 a=(xa[khi]-x)/h,b=(x-xa[klo])/h;
 y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
 return(y);
}

double spline_inter_der(double *xa,double *ya,double *y2a,int n,double x)
{
 int	klo,khi,k;
 double	a,b,h,y;

 klo=0,khi=n-1;
 while ( khi-klo>1 )
  {	k=(khi+klo)/2;
	if ( xa[k] > x )	khi=k;
	else			klo=k;
  }
 h=xa[khi]-xa[klo];
 a=(xa[khi]-x)/h,b=(x-xa[klo])/h;
 y=(-ya[klo]+ya[khi])/h+(-(3*a*a-1.0)*y2a[klo]+(3*b*b-1.0)*y2a[khi])*(h)/6.0;
 return(y);
}

/*****************************************************************************/

int natspline_coeff(double *y,int n,double *y2)
{
 int	i,k;
 double	p,qn,un,*u;

 u=(double *)malloc(sizeof(double)*n);
 y2[0]=u[0]=0.0;
 for ( i=1 ; i<n-1 ; i++ )
  {	p=0.5*y2[i-1]+2.0;
	y2[i]=-0.5/p;
	u[i]=(y[i+1]-y[i])-(y[i]-y[i-1]);
	u[i]=(3.0*u[i]-0.5*u[i-1])/p;
  }
 qn=un=0.0;
 y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
 for ( k=n-2 ; k>=0 ; k-- )
  {	y2[k]=y2[k]*y2[k+1]+u[k];					}
 free(u);

 return(0);
}

double natspline_inter(double *ya,double *y2a,int n,double x)
{
 int	klo,khi;
 double	a,b,y,d;
 
 klo=(int)x;
 khi=klo+1;d=x-(double)klo;
 a=(1.0-d),b=(d);
 y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])/6.0;
 return(y);
}


double natspline_inter_der(double *ya,double *y2a,int n,double x)
{
 int	klo,khi;
 double	a,b,y,d;
 
 klo=(int)x;
 khi=klo+1;d=x-(double)klo;
 a=(1.0-d),b=(d);
 y=(-ya[klo]+ya[khi])+(-(3*a*a-1.0)*y2a[klo]+(3*b*b-1.0)*y2a[khi])/6.0;
 return(y);
}
  
/*****************************************************************************/

int cyspline_coeff(double *y,int n,double *y2)
{
 if ( n <= 0 )
	return(-1); 
 else if ( n==1 )
  {	y2[0]=0.0;
	return(0);
  }
 else if ( n==2 )
  {	/*
	y2[0]=(4.0*(2*y[1]-2*y[0])-1.0*(2*y[0]-2*y[1]))/15.0;
	y2[1]=(4.0*(2*y[0]-2*y[1])-1.0*(2*y[1]-2*y[0]))/15.0;
	*/
 	y2[0]=9*(2*y[1]-2*y[0])/3.0;
 	y2[1]=9*(2*y[0]-2*y[1])/3.0;
	return(0);
  }
 else if ( n==3 )
  {	/*
	y2[0]=(5.0*(y[1]-2*y[0]+y[2])-1.0*(y[2]-2*y[1]+y[0])-1.0*(y[0]-2*y[2]+y[1]))/18.0;
	*/
	y2[0]=6*(-2*y[0]+y[1]+y[2])/3.0;
	y2[1]=6*(-2*y[1]+y[2]+y[0])/3.0;
	y2[2]=6*(-2*y[2]+y[0]+y[1])/3.0;
	return(0);
  }
 else
  {	int	*da,*db,i,j,s;
	da=(int *)malloc(sizeof(int)*(n+1));
	db=(int *)malloc(sizeof(int)*(n+1));
	da[0]=da[1]=db[0]=db[1]=0;
	da[2]=0,da[3]=1;
	db[2]=1,db[3]=0;
	for ( i=4 ; i<=n ; i++ )
	 {	if ( i%2 )	
		 {	da[i]=2*da[i-1]+da[i-2];
			db[i]=2*db[i-1]+db[i-2];
		 }
		else
		 {	da[i]=da[i-1]+da[i-2];
			db[i]=db[i-1]+db[i-2];
		 }
	 }
	for ( i=0 ; i<n ; i++ )
	 {	y2[i]=-6*(2*da[n]+db[n])*y[i];
		s=1;
		for ( j=0 ; j<(n-1)/2 ; j++,s=-s )
		 {	y2[i]+=s*6*(da[n-2*j]+db[n-2*j])*(y[(i+j+1)%n]+y[(i+n-j-1)%n]);	}
		if ( n%2==0 )
		 {	y2[i]+=s*6*(da[n-2*j]+db[n-2*j])*(y[(i+n/2)%n]);		}
		y2[i]/=(3*da[n]+db[n]);
	 }
	free(db);
	free(da);
	return(0);
  }			

}

double cyspline_inter(double *ya,double *y2a,int n,double x)
{
 int	klo,khi;
 double	a,b,y,d;
 
 klo=(int)x;
 d=x-(double)klo;
 if ( d<0.0 )	d+=1.0,klo--;
 if ( klo<0 )	klo+=n*((-klo)/n+2);
 klo=klo%n;
 khi=(klo+1)%n;
 a=(1.0-d),b=(d);
 y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])/6.0;
 return(y);
}

double cyspline_inter_der(double *ya,double *y2a,int n,double x)
{  
 int	klo,khi;
 double	a,b,y,d;
 
 klo=(int)x;
 d=x-(double)klo;
 if ( d<0.0 )	d+=1.0,klo--;
 if ( klo<0 )	klo+=n*((-klo)/n+2);
 klo=klo%n;
 khi=(klo+1)%n;
 a=(1.0-d),b=(d);
 y=(-ya[klo]+ya[khi])+(-(3*a*a-1.0)*y2a[klo]+(3*b*b-1.0)*y2a[khi])/6.0;
 return(y);
}

/*****************************************************************************/

int intspline_coeff(double *y,int p,double *l)
{
 double	*gam,bet;
 int	j;
 if ( p<=0 )	return(1);
 else if ( p==1 )
  {	l[0]=l[1]=y[0];
	return(0);
  }
 
 gam=(double *)malloc(sizeof(double)*(p+1));
 l[0]=2.0*y[0];
 bet=1.0;
 for ( j=1 ; j<p ; j++ )
  {	gam[j]=1.0/bet;
	bet=4.0-gam[j];
	l[j]=(3.0*(y[j]+y[j-1])-l[j-1])/bet;
  }
 gam[p]=1.0/bet;
 bet=1.0-gam[p];
 l[p]=(2.0*y[p-1]-l[p-1])/bet;
 for ( j=p-1 ; j>=0 ; j-- )
  {	l[j] -= gam[j+1]*l[j+1];	}
 free(gam);

 return(0);
}

double intspline_inter(double *y,double *l,int n,double x,double len)
{
 double tlen,ret,px,lx,ex;
 int	p;

 tlen=(double)(n+1);
 if ( x<0.0 || x+len>tlen )	return(0.0);

 ret=0.0;p=(int)x;px=x-(double)p;
 while ( len>0.0 )
  {	lx=len;
	ex=px+lx;
	if ( ex>1.0 )	lx=1.0-px,ex=1.0;
	
	if ( lx==1.0 )	ret+=y[p];
	else		ret+=lx*(l[p]+ 
				 (ex+px)*(-2*l[p]-l[p+1]+3*y[p])+ 
				 (ex*(ex+px)+px*px)*(l[p]+l[p+1]-2*y[p]));

	p++;
	x=(double)p;px=0.0;
	len-=lx;
  };
 return(ret);

}

/*****************************************************************************/
                         
