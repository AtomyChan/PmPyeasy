/*****************************************************************************/
/* star-draw.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2004, 2005, 2006; Pal, A. (apal@szofi.elte.hu)			     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include "math/expint/expint.h"
#include "math/poly.h"
#include "math/spline/biquad.h"
#include "math/spline/biquad-isc.h"

#include "tensor.h"
#include "stars.h"
#include "psf.h"

int star_draw_gauss
	(double **iarr,int sx,int sy,
	double x0,double y0,double is,double id,double ik)
{
 int		i,j;
 double		gs,gd,gk,det,dx,dy,is_deviated;
 static double	**itensor=NULL;
 static int	afsize=0;
 double		*iverx,*ihorx,*iverf,*ihorf,*ilx,*ily;
 
 if ( iarr==NULL )		return(-1);
 if ( sx<=0 || sy<=0 )		return(-1);

 det=is*is-id*id-ik*ik;
 if ( det<=0.0 )	return(0);
 gs=(is*is+id*id+ik*ik)/(det*det),
 gd=-2*is*id/(det*det),
 gk=-2*is*ik/(det*det);

 if ( gk == 0.0 )	is_deviated=0;
 else			is_deviated=1;

 if ( sx>afsize || sy>afsize || itensor==NULL )
  {	afsize=(sx>sy?sx:sy);
	if ( itensor != NULL )	tensor_free(itensor);
	itensor=(double **)tensor_alloc_2d(double,afsize+1,6);
  }
 ihorx=itensor[0],ihorf=itensor[1],ilx=itensor[2];
 iverx=itensor[3],iverf=itensor[4],ily=itensor[5];

 if ( ! is_deviated )
  {	double	sqx,sqy;

	sqx=sqrt(0.5*(gs+gd));
	sqy=sqrt(0.5*(gs-gd));
	for ( j=0 ; j<=sx ; j++ )
	 {	dx=(double)j-x0;
		ihorf[j]=erf(sqx*dx);
	 }
	for ( i=0 ; i<=sy ; i++ )
	 {	dy=(double)i-y0;
		iverf[i]=erf(sqy*dy);
	 }
	for ( j=0 ; j<sx ; j++ )	ihorf[j]=ihorf[j+1]-ihorf[j];
	for ( i=0 ; i<sy ; i++ )	iverf[i]=iverf[i+1]-iverf[i];

	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	iarr[i][j]=iverf[i]*ihorf[j];	}
	 }
  }
 else
  {	double		x1,y1,x2,y2,sqx,sqy;
	expintee	ee;
	
	sqx=sqrt(0.5*(gs+gd));
	sqy=sqrt(0.5*(gs-gd));

	for ( j=0 ; j<=sx ; j++ )
	 {	ilx[j]=dx=(double)j-x0;
		ihorf[j]=erf(sqx*dx);
		ihorx[j]=exp(-0.5*(gs+gd)*dx*dx);
	 }
	for ( i=0 ; i<=sy ; i++ )
	 {	ily[i]=dy=(double)i-y0;
		iverf[i]=erf(sqy*dy);
		iverx[i]=exp(-0.5*(gs-gd)*dy*dy);
	 }

	for ( i=0 ; i<sy ; i++ )
	 {	y1=ily[i];
		y2=ily[i+1];
		ee.expy1=iverx[i],ee.expy2=iverx[i+1];
		ee.erfy1=iverf[i],ee.erfy2=iverf[i+1];
		for ( j=0 ; j<sx ; j++ )
		 {	x1=ilx[j];
			x2=ilx[j+1];
			ee.expx1=ihorx[j],ee.expx2=ihorx[j+1];
			ee.erfx1=ihorf[j],ee.erfx2=ihorf[j+1];
			iarr[i][j]=expint_taylor_ee(gs,gd,gk,x1,x2,y1,y2,&ee);
		 }
	 }
  }

 return(0);

}


int star_draw_deviated
	(double **iarr,int sx,int sy,
	double x0,double y0,double gs,int order,double *mom)
{
 int		i,j,o,k,l;
 double		dx,dy;
 static double	**iverl=NULL,**ihorl=NULL,*fact=NULL;
 static int	afsize=0,aorder=0;

 if ( iarr==NULL )		return(-1);
 if ( sx<=0 || sy<=0 )		return(-1);

 if ( gs<=0.0 )			return(0);

 if ( sx>afsize || sy>afsize || order>aorder)
  {	afsize=(sx>sy?sx:sy);
	aorder=order;
	if ( ihorl != NULL )	tensor_free(ihorl);
	if ( iverl != NULL )	tensor_free(iverl);
	if ( fact  != NULL )	tensor_free(fact);
	ihorl=(double **)tensor_alloc_2d(double,aorder+1,afsize+1);
	iverl=(double **)tensor_alloc_2d(double,aorder+1,afsize+1);
	fact =(double  *)tensor_alloc_1d(double,(order+1)*(order+2)/2);
	for ( o=0,l=0,i=1 ; o<=order ; o++ )
	 {	j=i;
		for ( k=0 ; k<=o ; k++,l++ )
		 {	fact[l]=1.0/(double)(j);
			if ( k<o )	j=j*(k+1)/(o-k);
		 }
		i=i*(o+1);
	 }
  }

 for ( i=0 ; i<=sy ; i++ )
  {	dy=(double)(i)-y0;
	expint_primitive_list(gs,dy,order,iverl[i]);
  }
 for ( j=0 ; j<=sx ; j++ )
  {	dx=(double)(j)-x0;
	expint_primitive_list(gs,dx,order,ihorl[i]);
  }
	
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	iarr[i][j]=(iverl[i+1][0]-iverl[i][0])*(ihorl[j+1][0]-ihorl[j][0]);
		for ( o=2,l=0 ; o<=order ; o++ )
		 {	for ( k=0 ; k<=o ; k++,l++ )
			 {	iarr[i][j]+=fact[l+3]*mom[l]*(iverl[i+1][k]-iverl[i][k])*(ihorl[j+1][o-k]-ihorl[j][o-k]);		}
		 }
	 }
  }

 return(0);
}

int star_draw_psf(double **iarr,int sx,int sy,
	double x0,double y0,psf *p,double px,double py,
	double is,double id,double ik,double il)
{
 static double	**bqc=NULL,**coeff=NULL,*cpoly=NULL;
 static int	abx=0,aby=0,anvar=0;
 int		i,j,k,nvar,bx,by;
 int		is_deviated;
 double		x1,y1,x2,y2;
 double		lmxx,lmxy,lmyx,lmyy,det;
 double		xdl,ydl,xdr,ydr,xul,yul,xur,yur,afm,afo;

 if ( iarr==NULL )		return(-1);
 if ( p==NULL )			return(-1);

 det=is*is+il*il-id*id-ik*ik;
 if ( is<=0.0 && det<=0.0 )
  {	lmxx=1.0,lmyy=1.0;
	lmxy=0.0,lmyx=0.0;
  }
 else
  {	lmxx=(+is-id)/det;
	lmxy=(-ik-il)/det;
	lmyx=(-ik+il)/det;
	lmyy=(+is+id)/det;
  }
 if ( lmxy != 0 || lmyx != 0.0 )	is_deviated=1;
 else					is_deviated=0;

 nvar=(p->order+1)*(p->order+2)/2;
 bx=p->grid*(2*p->hsize+1);
 by=p->grid*(2*p->hsize+1);

 if ( bx>abx || by>aby )
  {	if ( bqc != NULL )	tensor_free(bqc); 
	if ( coeff != NULL )	tensor_free(coeff); 
	coeff=(double **)tensor_alloc_2d(double,bx,by);
	bqc  =(double **)tensor_alloc_2d(double,2*bx+1,2*by+1);
	abx=bx,aby=by;
  }
 if ( nvar>anvar )
  {	if ( cpoly != NULL )	tensor_free(cpoly); 
	cpoly=(double *)tensor_alloc_1d(double,nvar);
	anvar=nvar;
  }

 for ( i=0 ; i<by ; i++ )
  {	for ( j=0 ; j<bx ; j++ )
	 {	for ( k=0 ; k<nvar ; k++ )
		 {	cpoly[k]=p->coeff[k][i][j];		}
		coeff[i][j]=eval_2d_poly(px,py,p->order,cpoly,p->ox,p->oy,p->scale);
	 }
  }
 biquad_coeff(coeff,bx,by,bqc,NULL);

 afm=(double)p->grid;
 afo=(double)p->grid*(0.5+(double)p->hsize);
 for ( i=0 ; i<sy ; i++ )
  {  for ( j=0 ; j<sx ; j++ )
      {	x1=(double)(j+0)-x0;
	x2=(double)(j+1)-x0;
	y1=(double)(i+0)-y0;
	y2=(double)(i+1)-y0;
	if ( ! is_deviated )
	 {	x1=afm*lmxx*x1+afo;
		x2=afm*lmxx*x2+afo;
		y1=afm*lmyy*y1+afo;
		y2=afm*lmyy*y2+afo;
		if ( x1<0.0 )	x1=0.0;
		if ( x1>=bx )	x1=bx;
		if ( x2<0.0 )	x2=0.0;
		if ( x2>=bx )	x2=bx;
		if ( y1<0.0 )	y1=0.0;
		if ( y1>=by )	y1=by;
		if ( y2<0.0 )	y2=0.0;
		if ( y2>=by )	y2=by;
		if ( x1==x2 || y1==y2 )
			iarr[i][j]=0.0;	
		else
			iarr[i][j]=biquad_isc_int_rectangle(bqc,x1,y1,x2,y2);
	 }
	else
	 {	xdl=afm*(lmxx*x1+lmxy*y1)+afo,ydl=afm*(lmyx*x1+lmyy*y1)+afo;
		xdr=afm*(lmxx*x2+lmxy*y1)+afo,ydr=afm*(lmyx*x2+lmyy*y1)+afo;
		xul=afm*(lmxx*x1+lmxy*y2)+afo,yul=afm*(lmyx*x1+lmyy*y2)+afo;
		xur=afm*(lmxx*x2+lmxy*y2)+afo,yur=afm*(lmyx*x2+lmyy*y2)+afo;
		iarr[i][j]=biquad_isc_int_triangle(bqc,1,xdl,ydl,xdr,ydr,xul,yul,bx,by)+
			   biquad_isc_int_triangle(bqc,1,xdr,ydr,xur,yur,xul,yul,bx,by);
	 }
     }
  }

 return(0); 
}

