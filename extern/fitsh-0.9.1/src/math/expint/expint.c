/*****************************************************************************/
/* expint.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Standalone library to calculate one/two-dimensional Gaussian intergrals.  */
/* For more information (about the integrals and recurrence formulae)	     */
/* see the `FI package documentation`, chapter `Frequently Used Functions`.  */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2005-06, Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "expint.h"

/*****************************************************************************/

static double expfunct(double s,double d,double k,double x,double y)
{
 double ret;

 ret=exp(-0.5*((s+d)*(x*x)+(s-d)*(y*y))-k*x*y);
 return(ret);
}
double expint_numerical(double s,double d,double k,double x1,double x2,double y1,double y2)
{
 double	ret,x,y,dx,dy,w,fx,fy,dd;
 int	n,i,j;

 n=1000,ret=0.0;
 dx=(x2-x1)/(double)n,dy=(y2-y1)/(double)n;

 dd=dx*dy;
 for ( i=0 ; i<=n ; i++ )
  {	for ( j=0 ; j<=n ; j++ )
	 {	x=x1+dx*((double)j);
		y=y1+dy*((double)i);

		w=expfunct(s,d,k,x,y);

		if ( i==0 || i==n )	fy=1.0;
		else if ( i%2==1 )	fy=4.0;
		else			fy=2.0;
		if ( j==0 || j==n )	fx=1.0;
		else if ( j%2==1 )	fx=4.0;
		else			fx=2.0;
		
		ret+=fx*fy*w*dd;
	 }
  }	

 return(ret/9.0);
}

/*****************************************************************************/

double expint_taylor_ee(double s,double d,double k,double x1,double x2,double y1,double y2,expintee *ee)
{
 double	ssx,ssy,ret,retprev,iprevx,iprevy,icurrx,icurry,inextx,inexty,
	x1c,x2c,y1c,y2c,m;
 int	i,n;

 if ( s<=0.0 || s*s-d*d-k*k<=0.0 )	return(-1.0);

 ssx=1.0/(s+d),
 ssy=1.0/(s-d);
 
 x1c=ee->expx1,x2c=ee->expx2;
 y1c=ee->expy1,y2c=ee->expy2;

 iprevx=sqrt(0.5*M_PI*ssx)*(ee->erfx2-ee->erfx1),
 iprevy=sqrt(0.5*M_PI*ssy)*(ee->erfy2-ee->erfy1);
 icurrx=-(x2c-x1c)*ssx,
 icurry=-(y2c-y1c)*ssy;

 m=-k;
 ret=iprevx*iprevy+m*icurrx*icurry;

 n=100;
 for ( i=1 ; i<=n ; i++ )
  {	x1c*=x1,x2c*=x2,
	y1c*=y1,y2c*=y2;
	inextx=((double)i*iprevx-(x2c-x1c))*ssx,
	inexty=((double)i*iprevy-(y2c-y1c))*ssy;
	m=-m*k/(double)(i+1);
	retprev=ret;
	ret+=m*inextx*inexty;
	if ( ret==retprev )	break;
	iprevx=icurrx,icurrx=inextx;
	iprevy=icurry,icurry=inexty;
  }
/* fprintf(stderr,"iterations=%d\n",i); */

 return(ret);
}

double expint_taylor(double s,double d,double k,double x1,double x2,double y1,double y2)
{
 double		cx,cy,ret;
 expintee	ee;

 cx=s+d,cy=s-d; 
 ee.expx1=exp(-0.5*cx*x1*x1),ee.erfx1=erf(sqrt(0.5*cx)*x1);
 ee.expx2=exp(-0.5*cx*x2*x2),ee.erfx2=erf(sqrt(0.5*cx)*x2);
 ee.expy1=exp(-0.5*cy*y1*y1),ee.erfy1=erf(sqrt(0.5*cy)*y1);
 ee.expy2=exp(-0.5*cy*y2*y2),ee.erfy2=erf(sqrt(0.5*cy)*y2);

 ret=expint_taylor_ee(s,d,k,x1,x2,y1,y2,&ee);
 return(ret);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int expint_taylor_ee_diff(double s,double d,double k,double x1,double x2,double y1,double y2,double *dlist,expintee *ee)
{
 double	ssx,ssy,iprevx,iprevy,icurrx,icurry,inextx,inexty,
	i0,ix,iy,ixx,ixy,iyy,i0prev,x1c,x2c,y1c,y2c,m;
 int	i,n;

 if ( s<=0.0 || s*s-d*d-k*k<=0.0 )	return(-1);

 ssx=1.0/(s+d),
 ssy=1.0/(s-d);
 
 x1c=ee->expx1,x2c=ee->expx2;
 y1c=ee->expy1,y2c=ee->expy2;

 iprevx=sqrt(0.5*M_PI*ssx)*(ee->erfx2-ee->erfx1),
 iprevy=sqrt(0.5*M_PI*ssy)*(ee->erfy2-ee->erfy1);
 icurrx=-(x2c-x1c)*ssx,
 icurry=-(y2c-y1c)*ssy;

 m=1.0;
 i0=ix=iy=ixx=ixy=iyy=0.0;

 n=100;
 for ( i=1 ; i<=n ; i++ )
  {	x1c*=x1,x2c*=x2,
	y1c*=y1,y2c*=y2;
	inextx=((double)i*iprevx-(x2c-x1c))*ssx,
	inexty=((double)i*iprevy-(y2c-y1c))*ssy;
	i0prev=i0;
	i0+=m*iprevx*iprevy;
	ix+=m*icurrx*iprevy;
	iy+=m*iprevx*icurry;
	ixx+=m*inextx*iprevy;
	ixy+=m*icurrx*icurry;
	iyy+=m*iprevx*inexty;
	if ( i0prev==i0 )	break;
	m=-m*k/(double)i;
	iprevx=icurrx,icurrx=inextx;
	iprevy=icurry,icurry=inexty;
  }
 dlist[0]=i0;
 dlist[1]=(s+d)*ix+k*iy;
 dlist[2]=k*ix+(s-d)*iy;
 dlist[3]=-0.5*(ixx+iyy);
 dlist[4]=-0.5*(ixx-iyy);
 dlist[5]=-ixy;

 return(0);
}

int expint_taylor_ee_shift_diff(double s,double d,double k,double x0,double y0,double x1,double x2,double y1,double y2,double *dlist,expintee *ee)
{
 return(expint_taylor_ee_diff(s,d,k,x1-x0,x2-x0,y1-y0,y2-y0,dlist,ee));
}

int expint_taylor_diff(double s,double d,double k,double x1,double x2,double y1,double y2,double *dlist)
{
 double		cx,cy,mxex,mxef,myex,myef;
 expintee	ee;

 cx=s+d,cy=s-d; 
 mxex=-0.5*cx,mxef=sqrt(0.5*cx);
 myex=-0.5*cy,myef=sqrt(0.5*cy);
 ee.expx1=exp(mxex*x1*x1),ee.erfx1=erf(mxef*x1);
 ee.expx2=exp(mxex*x2*x2),ee.erfx2=erf(mxef*x2);
 ee.expy1=exp(myex*y1*y1),ee.erfy1=erf(myef*y1);
 ee.expy2=exp(myex*y2*y2),ee.erfy2=erf(myef*y2);

 return(expint_taylor_ee_diff(s,d,k,x1,x2,y1,y2,dlist,&ee));
}

/*****************************************************************************/

int expint_list_e(double s,double x1,double x2,int n,double *list,expinte *e)
{
 double	ss,iprev,icurr,inext,x1c,x2c;
 int	i;

 if ( s<=0.0 || list==NULL || n<0 )	return(-1);
 ss=1.0/(s);
 
 x1c=e->exp1,x2c=e->exp2;

 iprev=sqrt(0.5*M_PI*ss)*(e->erf2-e->erf1);
 list[0]=iprev;
 if ( n==0 )	return(0);
 icurr=-(x2c-x1c)*ss;
 list[1]=icurr;
 for ( i=1 ; i<n ; i++ )
  {	x1c*=x1,x2c*=x2;
	inext=((double)i*iprev-(x2c-x1c))*ss;
	list[i+1]=inext;
	iprev=icurr,icurr=inext;
  }
 return(0);
}

int expint_list(double s,double x1,double x2,int n,double *list)
{
 expinte	e;
 double		mex,mef;

 mex=-0.5*s;
 mef=sqrt(0.5*s);

 e.exp1=exp(mex*x1*x1),e.erf1=erf(mef*x1);
 e.exp2=exp(mex*x2*x2),e.erf2=erf(mef*x2);

 return(expint_list_e(s,x1,x2,n,list,&e));
}

/*****************************************************************************/

int expint_primitive_list(double s,double x,int n,double *list)
{
 double	ss,iprev,icurr,inext,xc;
 int	i;

 if ( s<=0.0 || list==NULL || n<0 )	return(-1);
 ss=1.0/(s);
 
 xc=exp(-0.5*s*x*x);
 iprev=sqrt(0.5*M_PI*ss)*erf(sqrt(0.5*s)*x);
 list[0]=iprev;
 if ( n==0 )	return(0);
 icurr=-xc*ss;
 list[1]=icurr;
 for ( i=1 ; i<n ; i++ )
  {	xc*=x;
	inext=((double)i*iprev-xc)*ss;
	list[i+1]=inext;
	iprev=icurr,icurr=inext;
  }
 return(0);
}
                                               
