/*****************************************************************************/
/* star-model.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Functions related to star model fitting to a FITS image. This is the main */
/* module to create a star-list from a candidate-list.			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2006, Pal, A. (apal@szofi.elte.hu), part of the 'FI' package	     */
/*****************************************************************************/

#define	__STAR_MODEL_C	1

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <fits/fits.h>

#include "fitsmask.h"
#include "math/fit/lmfit.h"
#include "math/poly.h"
#include "math/expint/expint.h"
#include "index/multiindex.h"

#include "fitsh.h"
#include "common.h"
#include "stars.h"
#include "tensor.h"
#include "statistics.h"

#include "star-model.h"

/*****************************************************************************/

typedef struct
 {	double	efx,efy;
	double	exx,exy;
 } eeval;

typedef struct
 {	int	ix0,iy0,isx,isy;
	double	gs,gd,x0,y0;
	eeval	*eevals;
	double	efc,rc2s,efm;
	int	derivative_level;	/* < 0: no derivatives,		*/
					/*   0: location /A+B/,		*/
					/*   1: location /A+B+X+Y/	*/
					/* >=2: location + shape	*/
 } gfparam;

/*****************************************************************************/

/* F=B+Aexp(-1/2*S*((x-x0)^2+(y-y0)^2))					     */ 
/* a[] = { A, B, x0, y0, S }						     */
/* xpnt = { ix, iy }							     */
void gauss_2d_nmom_exact(void *xpnt,double *a,double *yy,double *dyda,void *param)
{
 double		x,y,x1,x2,y1,y2,x0,y0,s;
 double		efx1,efx2,efy1,efy2,efmul,efc,efm;
 double		exx1,exx2,exy1,exy2,rc2s;

 x=((ipoint *)xpnt)->x,x0=a[2],x1=x-x0,x2=x1+1.0;
 y=((ipoint *)xpnt)->y,y0=a[3],y1=y-y0,y2=y1+1.0;
 s=a[4];

 efc=sqrt(s/2);
 efx1=erf(efc*x1),efx2=erf(efc*x2);
 efy1=erf(efc*y1),efy2=erf(efc*y2);
 rc2s=0.5/s;

 efmul=(efx2-efx1)*(efy2-efy1)*M_PI*rc2s;

 *yy=a[0]*efmul+a[1];

 if ( dyda != NULL )
  {	dyda[0]=efmul;				/* \partial_A	*/
	dyda[1]=1.0;				/* \partial_B	*/
	efm=sqrt(M_PI*rc2s);
	efx1*=efm,efx2*=efm;
	efy1*=efm,efy2*=efm;
	exx1=exp(-0.5*s*x1*x1),exx2=exp(-0.5*s*x2*x2);
	exy1=exp(-0.5*s*y1*y1),exy2=exp(-0.5*s*y2*y2);

	dyda[2]=a[0]*(efy2-efy1)*(exx1-exx2);	/* \partial_x0	*/
	dyda[3]=a[0]*(efx2-efx1)*(exy1-exy2);	/* \partial_y0	*/
	
	dyda[4]=a[0]*rc2s*( (efy2-efy1)*(x2*exx2-efx2-x1*exx1+efx1)+
			    (efx2-efx1)*(y2*exy2-efy2-y1*exy1+efy1) );
	
  }
 
}

/* F=B+Aexp(-1/2*S*((x-x0)^2+(y-y0)^2))					     */ 
/* a[] = { A, B, x0, y0, S }						     */
/* xpnt = { ix, iy }							     */
void gauss_2d_nmom_exact_init(void *xpnt,double *a,double *yy,double *dyda,void *param)
{
 double		x,y,x1,x2,y1,y2,x0,y0,gs;
 double		efx1,efx2,efy1,efy2,efmul,efc,efm;
 double		exx1,exx2,exy1,exy2,rc2s;
 int		ix,iy;
 gfparam	*gfp;

 if ( param==NULL )
  {	gauss_2d_nmom_exact(xpnt,a,yy,dyda,NULL);
	return;
  }

 gfp=(gfparam *)param;

 x0=a[2],y0=a[3];
 gs=a[4];

 if ( x0 != gfp->x0 || y0 != gfp->y0 || gs != gfp->gs )
  {	int	j,k;

	efc=sqrt(0.5*gs);

	if ( gfp->derivative_level>0 )
	 {	for ( k=0 ; k<gfp->isx+1 ; k++ )
		 {	x=(double)(k+gfp->ix0)-x0;
			gfp->eevals[k].efx=erf(efc*x);
			gfp->eevals[k].exx=exp(-0.5*gs*x*x);
		 }
		for ( j=0 ; j<gfp->isy+1 ; j++ )
		 {	y=(double)(j+gfp->iy0)-y0;
			gfp->eevals[j].efy=erf(efc*y);
			gfp->eevals[j].exy=exp(-0.5*gs*y*y);
		 }
	 }
	else
	 {	for ( k=0 ; k<gfp->isx+1 ; k++ )
		 {	x=(double)(k+gfp->ix0)-x0;
			gfp->eevals[k].efx=erf(efc*x);
		 }
		for ( j=0 ; j<gfp->isy+1 ; j++ )
		 {	y=(double)(j+gfp->iy0)-y0;
			gfp->eevals[j].efy=erf(efc*y);
		 }
	 }

	gfp->efc=efc;
	gfp->rc2s=0.5/gs;
	gfp->efm=sqrt(M_PI*gfp->rc2s);
	
	gfp->x0=x0;
	gfp->y0=y0;
	gfp->gs=gs;
  }
 else
	gfp->rc2s=0.0;
 
 ix=((ipoint *)xpnt)->x-gfp->ix0;
 iy=((ipoint *)xpnt)->y-gfp->iy0;

 x=((ipoint *)xpnt)->x,x1=x-x0,x2=x1+1.0;
 y=((ipoint *)xpnt)->y,y1=y-y0,y2=y1+1.0;

 efc=gfp->efc;
 rc2s=gfp->rc2s;
 efx1=gfp->eevals[ix+0].efx;
 efx2=gfp->eevals[ix+1].efx;
 efy1=gfp->eevals[iy+0].efy;
 efy2=gfp->eevals[iy+1].efy;

 efmul=(efx2-efx1)*(efy2-efy1)*M_PI*rc2s;

 *yy=a[0]*efmul+a[1];

 if ( dyda != NULL && gfp->derivative_level>=0 )
  {	
	dyda[0]=efmul;				/* \partial_A	*/
	dyda[1]=1.0;				/* \partial_B	*/
	if ( gfp->derivative_level>0 )
	 {	efm=gfp->efm;
		efx1*=efm,efx2*=efm;
		efy1*=efm,efy2*=efm;
		exx1=gfp->eevals[ix+0].exx;
		exx2=gfp->eevals[ix+1].exx;
		exy1=gfp->eevals[iy+0].exy;
		exy2=gfp->eevals[iy+1].exy;

		dyda[2]=a[0]*(efy2-efy1)*(exx1-exx2);	/* \partial_x0	*/
		dyda[3]=a[0]*(efx2-efx1)*(exy1-exy2);	/* \partial_y0	*/
	
		if ( gfp->derivative_level>=2 )
		 {	dyda[4]=a[0]*rc2s*( (efy2-efy1)*(x2*exx2-efx2-x1*exx1+efx1)+
			    (efx2-efx1)*(y2*exy2-efy2-y1*exy1+efy1) );
		 }
	 }
  }
 
}

/*****************************************************************************/

/* F=B+Aexp(-1/2*[S((x-x0)^2+(y-y0)^2)+D((x-x0)^2-(y-y0)^2)+2K(x-x0)(y-y0)]) */
/* a[] = { A, B, x0, y0, S, D, K }					     */
/* xpnt = { ix, iy }							     */
void gauss_2d_wmom_exact(void *xpnt,double *a,double *yy,double *dyda,void *p)
{
 double x,y,x1,x2,y1,y2,x0,y0,gs,gd,gk,dlist[6];

 x=((ipoint *)xpnt)->x,x0=a[2],x1=x-x0,x2=x1+1.0;
 y=((ipoint *)xpnt)->y,y0=a[3],y1=y-y0,y2=y1+1.0;
 gs=a[4],gd=a[5],gk=a[6];

 expint_taylor_diff(gs,gd,gk,x1,x2,y1,y2,dlist);

 *yy=a[0]*dlist[0]+a[1];

 if ( dyda != NULL )
  {	dyda[0]=dlist[0];		/* \partial_A	*/
	dyda[1]=1.0;			/* \partial_B	*/
	dyda[2]=a[0]*dlist[1];		/* \partial_x0	*/
	dyda[3]=a[0]*dlist[2];		/* \partial_y0	*/

	dyda[4]=a[0]*dlist[3];
	dyda[5]=a[0]*dlist[4];
	dyda[6]=a[0]*dlist[5];
  }
 
}

/* F=B+Aexp(-1/2*[S((x-x0)^2+(y-y0)^2)+D((x-x0)^2-(y-y0)^2)+2K(x-x0)(y-y0)]) */
/* a[] = { A, B, x0, y0, S, D, K }					     */
/* xpnt = { ix, iy }							     */
void gauss_2d_wmom_exact_init(void *xpnt,double *a,double *yy,double *dyda,void *param)
{
 double		x,y,x1,x2,y1,y2,x0,y0,gs,gd,gk,dlist[6],efcx,efcy;
 int		ix,iy;
 gfparam	*gfp;
 expintee	ee;

 if ( param==NULL )
  {	gauss_2d_wmom_exact(xpnt,a,yy,dyda,NULL);
	return;
  }
 
 gfp=(gfparam *)param;

 x=((ipoint *)xpnt)->x,x0=a[2],x1=x-x0,x2=x1+1.0;
 y=((ipoint *)xpnt)->y,y0=a[3],y1=y-y0,y2=y1+1.0;
 gs=a[4],gd=a[5],gk=a[6];

 if ( x0 != gfp->x0 || y0 != gfp->y0 || gs != gfp->gs || gd != gfp->gd )
  {	int	j,k;

	efcx=sqrt(0.5*(gs+gd));
	efcy=sqrt(0.5*(gs-gd));

	for ( k=0 ; k<gfp->isx+1 ; k++ )
	 {	x=(double)(k+gfp->ix0)-x0;
		gfp->eevals[k].efx=erf(efcx*x);
		gfp->eevals[k].exx=exp(-0.5*(gs+gd)*x*x);
	 }
	for ( j=0 ; j<gfp->isy+1 ; j++ )
	 {	y=(double)(j+gfp->iy0)-y0;
		gfp->eevals[j].efy=erf(efcy*y);
		gfp->eevals[j].exy=exp(-0.5*(gs-gd)*y*y);
	 }
	
	gfp->x0=x0;
	gfp->y0=y0;
	gfp->gs=gs;
	gfp->gd=gd;
  }
  
 ix=((ipoint *)xpnt)->x-gfp->ix0;
 iy=((ipoint *)xpnt)->y-gfp->iy0;

 ee.expx1=gfp->eevals[ix+0].exx,ee.expx2=gfp->eevals[ix+1].exx;
 ee.expy1=gfp->eevals[iy+0].exy,ee.expy2=gfp->eevals[iy+1].exy;
 ee.erfx1=gfp->eevals[ix+0].efx,ee.erfx2=gfp->eevals[ix+1].efx;
 ee.erfy1=gfp->eevals[iy+0].efy,ee.erfy2=gfp->eevals[iy+1].efy;

 expint_taylor_ee_diff(gs,gd,gk,x1,x2,y1,y2,dlist,&ee);

 *yy=a[0]*dlist[0]+a[1];

 if ( dyda != NULL && gfp->derivative_level>=0 )
  {	dyda[0]=dlist[0];		/* \partial_A	*/
	dyda[1]=1.0;			/* \partial_B	*/
	if ( gfp->derivative_level>0 )
	 {	dyda[2]=a[0]*dlist[1];		/* \partial_x0	*/
		dyda[3]=a[0]*dlist[2];		/* \partial_y0	*/
	 	
		if ( gfp->derivative_level>=2 )
		 {	dyda[4]=a[0]*dlist[3];
			dyda[5]=a[0]*dlist[4];
			dyda[6]=a[0]*dlist[5];
	 	 }
	 }
  }
 
}

/*****************************************************************************/

static int wdev_order;

/* F=B+Aexp(-1/2*[S((x-x0)^2+(y-y0)^2)]) *				     */
/*     ( 1 + 1/2*M_20(x-x0)^2+M_11*(x-x0)(y-y0)+1/2*M_02(y-y0)^2 + ... )     */ 
/* a[] = { A, B, x0, y0, S, M20, M11, M02, M30, M21, M12, M03, M40, ... }    */
/* xpnt = { ix, iy }							     */
void gauss_2d_wdev_exact_init(void *xpnt,double *a,double *yy,double *dyda,void *param)
{
 int		order,ndev,i,j,k,jx,jy;
 double		x,y,x1,x2,y1,y2,x0,y0,i0,s,w;
 double		cf[MXC],cfi[MXT],cfdx[MXT],cfdy[MXT],cfds[MXT],
		ix[MXO+3],iy[MXO+3];
 gfparam	*gfp;
 int		derivative_level;
 expinte	e;

 if ( param==NULL )
  {	gfp=NULL;
	derivative_level=2;
  }
 else
  {	gfp=(gfparam *)param;
	derivative_level=gfp->derivative_level;
  }

 x=((ipoint *)xpnt)->x,x0=a[2],x1=x-x0,x2=x1+1.0;
 y=((ipoint *)xpnt)->y,y0=a[3],y1=y-y0,y2=y1+1.0;
 s=a[4];

 order=wdev_order;
 ndev=(order+1)*(order+2)/2;

 cf[0]=1.0;
 cf[1]=cf[2]=0.0;
 for ( i=3 ; i<ndev ; i++ )
  {	cf[i]=a[i+2];			}

 if ( gfp==NULL )
  {	expint_list(s,x1,x2,order+2,ix);
	expint_list(s,y1,y2,order+2,iy);
  }
 else
  {	jx=((ipoint *)xpnt)->x-gfp->ix0;
	e.exp1=gfp->eevals[jx+0].exx;
	e.exp2=gfp->eevals[jx+1].exx;
	e.erf1=gfp->eevals[jx+0].efx;
	e.erf2=gfp->eevals[jx+1].efx;
	expint_list_e(s,x1,x2,order+2,ix,&e);
	jy=((ipoint *)xpnt)->y-gfp->iy0;
	e.exp1=gfp->eevals[jy+0].exy;
	e.exp2=gfp->eevals[jy+1].exy;
	e.erf1=gfp->eevals[jy+0].efy;
	e.erf2=gfp->eevals[jy+1].efy;
	expint_list_e(s,y1,y2,order+2,iy,&e);
  }

 switch ( order )
  {  case 2:
 	cfi[I00]=1.0;
	cfi[I10]=0.0;
	cfi[I01]=0.0;
	cfi[I20]=M20*cf[I20];
	cfi[I11]=M11*cf[I11];
	cfi[I02]=M02*cf[I02];

	i0=ix[0]*iy[0]+
		cfi[I20]*ix[2]*ix[0]+
		cfi[I11]*ix[1]*iy[1]+
		cfi[I02]*ix[0]*iy[2];

	*yy=a[1]+a[0]*i0;
	if ( dyda != NULL && derivative_level>=0 )
	 {	dyda[0]=i0;			/* \partial_A	*/
		dyda[1]=1.0;			/* \partial_B	*/

		if ( derivative_level>=1 )
		 {	cfdx[I00]=0.0;
			cfdx[I10]=s-cf[I20];
			cfdx[I01]=-cf[I11];
			cfdx[I20]=0.0;
			cfdx[I11]=0.0;
			cfdx[I02]=0.0;
			cfdx[I30]=M20*cf[I20]*s;
			cfdx[I21]=M11*cf[I11]*s;
			cfdx[I12]=M02*cf[I02]*s;
			cfdx[I03]=0.0;
			for ( i=0,w=0.0,k=0 ; i<=3 ; i++ )
			 {	for ( j=0 ; j<=i ; j++,k++ )
					w+=cfdx[k]*ix[i-j]*iy[j];
			 }
			dyda[2]=w*a[0];			/* \partial_x0	*/

			cfdy[I00]=0.0;
			cfdy[I10]= -cf[I11];
			cfdy[I01]=s-cf[I02];
			cfdy[I20]=0.0;
			cfdy[I11]=0.0;
			cfdy[I02]=0.0;
			cfdy[I30]=0.0;
			cfdy[I21]=M20*cf[I20]*s;
			cfdy[I12]=M11*cf[I11]*s;
			cfdy[I03]=M02*cf[I02]*s;
			for ( i=0,w=0.0,k=0 ; i<=3 ; i++ )
			 {	for ( j=0 ; j<=i ; j++,k++ )
					w+=cfdy[k]*ix[i-j]*iy[j];
			 }
			dyda[3]=w*a[0];			/* \partial_y0	*/
		 }

		if ( derivative_level>=2 )
		 {	cfds[I00]=0.0;
			cfds[I10]=0.0;
			cfds[I01]=0.0;
			cfds[I20]=1.0;
			cfds[I11]=0.0;
			cfds[I02]=1.0;
			cfds[I30]=0.0;
			cfds[I21]=0.0;
			cfds[I12]=0.0;
			cfds[I03]=0.0;
			cfds[I40]=M20*cf[I20];
			cfds[I31]=M11*cf[I11];
			cfds[I22]=M20*cf[I20]+M02*cf[I02];
			cfds[I13]=M11*cf[I11];
			cfds[I04]=M02*cf[I02];
			for ( i=0,dyda[4]=0.0,k=0 ; i<=4 ; i++ )
			 {	for ( j=0 ; j<=i ; j++,k++ )
					dyda[4]+=cfds[k]*ix[i-j]*iy[j];	
			 }
			dyda[4]=-0.5*a[0]*dyda[4];
	
			dyda[5]=a[0]*M20*ix[2]*iy[0];
			dyda[6]=a[0]*M11*ix[1]*iy[1];
			dyda[7]=a[0]*M02*ix[0]*iy[2];
		 }
	 } 
	break;

     case 3:
 	cfi[I00]=1.0;
	cfi[I10]=0.0;
	cfi[I01]=0.0;
	cfi[I20]=M20*cf[I20];
	cfi[I11]=M11*cf[I11];
	cfi[I02]=M20*cf[I02];
	cfi[I30]=M30*cf[I30];
	cfi[I21]=M21*cf[I21];
	cfi[I12]=M12*cf[I12];
	cfi[I03]=M03*cf[I03];
	
	i0=ix[0]*iy[0]+
		cfi[I20]*ix[2]*ix[0]+
		cfi[I11]*ix[1]*iy[1]+
		cfi[I02]*ix[0]*iy[2]+
		cfi[I30]*ix[3]*iy[0]+
		cfi[I21]*ix[2]*iy[1]+
		cfi[I12]*ix[1]*iy[2]+
		cfi[I03]*ix[0]*iy[3];

	*yy=a[1]+a[0]*i0;
	if ( dyda != NULL && derivative_level>=0 )
	 {	dyda[0]=i0;			/* \partial_A	*/
		dyda[1]=1.0;			/* \partial_B	*/

		if ( derivative_level>=1 )
		 {	cfdx[I00]=0.0;
			cfdx[I10]=s-M10*cf[I20];
			cfdx[I01]= -M01*cf[I11];
			cfdx[I20]= -M20*cf[I30];
			cfdx[I11]= -M11*cf[I21];
			cfdx[I02]= -M02*cf[I12];
			cfdx[I30]=M20*s*cf[I20];
			cfdx[I21]=M11*s*cf[I11];
			cfdx[I12]=M02*s*cf[I02];
			cfdx[I03]=0.0;
			cfdx[I40]=M30*s*cf[I30];
			cfdx[I31]=M21*s*cf[I21];
			cfdx[I22]=M12*s*cf[I12];
			cfdx[I13]=M03*s*cf[I03];
			cfdx[I04]=0.0;
			for ( i=0,w=0.0,k=0 ; i<=4 ; i++ )
			 {	for ( j=0 ; j<=i ; j++,k++ )
				w+=cfdx[k]*ix[i-j]*iy[j];
				 }
			dyda[2]=w*a[0];			/* \partial_x0	*/
	
			cfdy[I00]=0.0;
			cfdy[I10]= -M10*cf[I11];
			cfdy[I01]=s-M01*cf[I02];
			cfdy[I20]= -M20*cf[I21];
			cfdy[I11]= -M11*cf[I12];
			cfdy[I02]= -M02*cf[I30];
			cfdy[I30]=0.0;
			cfdy[I21]=M20*s*cf[I20];
			cfdy[I12]=M11*s*cf[I11];
			cfdy[I03]=M02*s*cf[I02];
			cfdy[I40]=0.0;
			cfdy[I31]=M30*s*cf[I30];
			cfdy[I22]=M21*s*cf[I21];
			cfdy[I13]=M12*s*cf[I12];
			cfdy[I04]=M03*s*cf[I03];
			for ( i=0,w=0.0,k=0 ; i<=4 ; i++ )
			 {	for ( j=0 ; j<=i ; j++,k++ )
					w+=cfdy[k]*ix[i-j]*iy[j];
			 }
			dyda[3]=w*a[0];			/* \partial_y0	*/
		 }
	
		if ( derivative_level>=2 )
	 	 {	cfds[I00]=0.0;
			cfds[I10]=0.0;
			cfds[I01]=0.0;
			cfds[I20]=1.0;
			cfds[I11]=0.0;
			cfds[I02]=1.0;
			cfds[I30]=0.0;
			cfds[I21]=0.0;
			cfds[I12]=0.0;
			cfds[I03]=0.0;
			cfds[I40]=	      M20*cf[I20];
			cfds[I31]=            M11*cf[I11];
			cfds[I22]=M20*cf[I20]+M02*cf[I02];
			cfds[I13]=M11*cf[I11];
			cfds[I04]=M02*cf[I02];
			cfds[I50]=            M30*cf[I30];
			cfds[I41]=            M21*cf[I21];
			cfds[I32]=M30*cf[I30]+M12*cf[I12];
			cfds[I23]=M21*cf[I21]+M03*cf[I03];
			cfds[I14]=M12*cf[I12];
			cfds[I05]=M03*cf[I03];
			for ( i=0,dyda[4]=0.0,k=0 ; i<=5 ; i++ )
			 {	for ( j=0 ; j<=i ; j++,k++ )
					dyda[4]+=cfds[k]*ix[i-j]*iy[j];	
			 }
			dyda[ 4]=-0.5*a[0]*dyda[4];
	
			dyda[ 5]=a[0]*M20*ix[2]*iy[0];
			dyda[ 6]=a[0]*M11*ix[1]*iy[1];
			dyda[ 7]=a[0]*M02*ix[0]*iy[2];
			dyda[ 8]=a[0]*M30*ix[3]*iy[0];
			dyda[ 9]=a[0]*M21*ix[2]*iy[1];
			dyda[10]=a[0]*M12*ix[1]*iy[2];
			dyda[11]=a[0]*M03*ix[0]*iy[3];
		 }
	 } 
	break;

     case 4:
 	cfi[I00]=1.0;
	cfi[I10]=0.0;
	cfi[I01]=0.0;
	cfi[I20]=M20*cf[I20];
	cfi[I11]=M11*cf[I11];
	cfi[I02]=M20*cf[I02];
	cfi[I30]=M30*cf[I30];
	cfi[I21]=M21*cf[I21];
	cfi[I12]=M12*cf[I12];
	cfi[I03]=M03*cf[I03];
	cfi[I40]=M40*cf[I40];
	cfi[I31]=M31*cf[I31];
	cfi[I22]=M22*cf[I22];
	cfi[I13]=M13*cf[I13];
	cfi[I04]=M04*cf[I04];
	
	i0=ix[0]*iy[0]+
		cfi[I20]*ix[2]*ix[0]+
		cfi[I11]*ix[1]*iy[1]+
		cfi[I02]*ix[0]*iy[2]+
		cfi[I30]*ix[3]*iy[0]+
		cfi[I21]*ix[2]*iy[1]+
		cfi[I12]*ix[1]*iy[2]+
		cfi[I03]*ix[0]*iy[3]+
		cfi[I40]*ix[4]*iy[0]+
		cfi[I31]*ix[3]*iy[1]+
		cfi[I22]*ix[2]*iy[2]+
		cfi[I13]*ix[1]*iy[3]+
		cfi[I04]*ix[0]*iy[4];

	*yy=a[1]+a[0]*i0;
	if ( dyda != NULL && derivative_level>=0 )
	 {	dyda[0]=i0;		/* \partial_A	*/
		dyda[1]=1.0;		/* \partial_B	*/

		if ( derivative_level>=1 )
	 	 {	cfdx[I00]=0.0;
			cfdx[I10]= s		-M10*cf[I20];
			cfdx[I01]=		-M01*cf[I11];
			cfdx[I20]=		-M20*cf[I30];
			cfdx[I11]=		-M11*cf[I21];
			cfdx[I02]=		-M02*cf[I12];
			cfdx[I30]=M20*s*cf[I20]	-M30*cf[I40];
			cfdx[I21]=M11*s*cf[I11]	-M21*cf[I31];
			cfdx[I12]=M02*s*cf[I02]	-M12*cf[I22];
			cfdx[I03]=		-M03*cf[I13];
			cfdx[I40]=M30*s*cf[I30];
			cfdx[I31]=M21*s*cf[I21];
			cfdx[I22]=M12*s*cf[I12];
			cfdx[I13]=M03*s*cf[I03];
			cfdx[I04]=0.0;
			cfdx[I50]=M40*s*cf[I40];
			cfdx[I41]=M31*s*cf[I31];
			cfdx[I32]=M22*s*cf[I22];
			cfdx[I23]=M13*s*cf[I13];
			cfdx[I14]=M04*s*cf[I04];
			cfdx[I05]=0.0;
			for ( i=0,w=0.0,k=0 ; i<=5 ; i++ )
			 {	for ( j=0 ; j<=i ; j++,k++ )
					w+=cfdx[k]*ix[i-j]*iy[j];
			 }
			dyda[2]=w*a[0];		/* \partial_x0	*/
	
			cfdy[I00]=0.0;
			cfdy[I10]=		-M10*cf[I11];
			cfdy[I01]= s		-M01*cf[I02];
			cfdy[I20]=		-M20*cf[I21];
			cfdy[I11]=		-M11*cf[I12];
			cfdy[I02]=		-M02*cf[I03];
			cfdy[I30]=		-M30*cf[I31];
			cfdy[I21]=M20*s*cf[I20]	-M21*cf[I22];
			cfdy[I12]=M11*s*cf[I11]	-M12*cf[I13];
			cfdy[I03]=M02*s*cf[I02]	-M03*cf[I04];
			cfdy[I40]=0.0;
			cfdy[I31]=M30*s*cf[I30];
			cfdy[I22]=M21*s*cf[I21];
			cfdy[I13]=M12*s*cf[I12];
			cfdy[I04]=M03*s*cf[I03];
			cfdy[I50]=0.0;
			cfdy[I41]=M40*s*cf[I40];
			cfdy[I32]=M31*s*cf[I31];
			cfdy[I23]=M22*s*cf[I22];
			cfdy[I14]=M13*s*cf[I13];
			cfdy[I05]=M04*s*cf[I04];
			for ( i=0,w=0.0,k=0 ; i<=5 ; i++ )
			 {	for ( j=0 ; j<=i ; j++,k++ )
					w+=cfdy[k]*ix[i-j]*iy[j];
			 }
			dyda[3]=w*a[0];		/* \partial_y0	*/
		 }
	
		if ( derivative_level>=2 )
		 {	cfds[I00]=0.0;
			cfds[I10]=0.0;
			cfds[I01]=0.0;
			cfds[I20]=1.0;
			cfds[I11]=0.0;
			cfds[I02]=1.0;
			cfds[I30]=0.0;
			cfds[I21]=0.0;
			cfds[I12]=0.0;
			cfds[I03]=0.0;
			cfds[I40]=	      M20*cf[I20];
			cfds[I31]=            M11*cf[I11];
			cfds[I22]=M20*cf[I20]+M02*cf[I02];
			cfds[I13]=M11*cf[I11];
			cfds[I04]=M02*cf[I02];
			cfds[I50]=            M30*cf[I30];
			cfds[I41]=            M21*cf[I21];
			cfds[I32]=M30*cf[I30]+M12*cf[I12];
			cfds[I23]=M21*cf[I21]+M03*cf[I03];
			cfds[I14]=M12*cf[I12];
			cfds[I05]=M03*cf[I03];
			cfds[I60]=	      M40*cf[I40];
			cfds[I51]=	      M31*cf[I31];
			cfds[I42]=M40*cf[I40]+M22*cf[I22];
			cfds[I33]=M31*cf[I31]+M13*cf[I13];
			cfds[I24]=M22*cf[I22]+M04*cf[I04];
			cfds[I15]=M13*cf[I13];
			cfds[I06]=M04*cf[I04];
			for ( i=0,dyda[4]=0.0,k=0 ; i<=6 ; i++ )
			 {	for ( j=0 ; j<=i ; j++,k++ )
					dyda[4]+=cfds[k]*ix[i-j]*iy[j];	
			 }
			dyda[ 4]=-0.5*a[0]*dyda[4];
	
			dyda[ 5]=a[0]*M20*ix[2]*iy[0];
			dyda[ 6]=a[0]*M11*ix[1]*iy[1];
			dyda[ 7]=a[0]*M02*ix[0]*iy[2];
			dyda[ 8]=a[0]*M30*ix[3]*iy[0];
			dyda[ 9]=a[0]*M21*ix[2]*iy[1];
			dyda[10]=a[0]*M12*ix[1]*iy[2];
			dyda[11]=a[0]*M03*ix[0]*iy[3];
			dyda[12]=a[0]*M40*ix[4]*iy[0];
			dyda[13]=a[0]*M31*ix[3]*iy[1];
			dyda[14]=a[0]*M22*ix[2]*iy[2];
			dyda[15]=a[0]*M13*ix[1]*iy[3];
			dyda[16]=a[0]*M04*ix[0]*iy[4];
		 } 
	 }
	break;
  }
 
}

/* F=B+Aexp(-1/2*[S((x-x0)^2+(y-y0)^2)]) *				     */
/*     ( 1 + 1/2*M_20(x-x0)^2+M_11*(x-x0)(y-y0)+1/2*M_02(y-y0)^2 + ... )     */ 
/* a[] = { A, B, x0, y0, S, M20, M11, M02, M30, M21, M12, M03, M40, ... }    */
/* xpnt = { ix, iy }							     */
void gauss_2d_wdev_exact(void *xpnt,double *a,double *yy,double *dyda,void *p)
{
 gauss_2d_wdev_exact_init(xpnt,a,yy,dyda,NULL);
}

/*****************************************************************************/

/* level: 				*/
/*    -	0: only A and B			*/
/*    -	1: A, B + position (X, Y)	*/
/*    -	2: A, B, X, Y + shape (S)	*/

int refine_fit_model_gauss(int nipoint,double *yvals,void **fitpnt,
	starlocation *loc,starshape *shp,starfit *sfp,int level)
{
 int		i;
 int		ix0,iy0,ix1,iy1,isx,isy,ix,iy,inum,nfit;
 double		lam,afp[5],gs;
 void		(*g2d_nmom)(void *,double *,double *,double *,void *);
 gfparam	gfp;
 static	eeval	*eevals=NULL;
 static int	neeval=0;

 g2d_nmom=gauss_2d_nmom_exact_init;

 if ( nipoint < 5 )	return(1);

 ix0=ix1=((ipoint *)fitpnt[0])->x,
 iy0=iy1=((ipoint *)fitpnt[0])->y;
 for ( i=1 ; i<nipoint ; i++ )
  {	ix=((ipoint *)fitpnt[i])->x,
	iy=((ipoint *)fitpnt[i])->y;
	if ( ix<ix0 )	ix0=ix;	
	if ( ix>ix1 )	ix1=ix;
	if ( iy<iy0 )	iy0=iy; 
	if ( iy>iy1 )	iy1=iy;
  }

 gfp.ix0=ix0;
 gfp.iy0=iy0;
 gfp.isx=isx=ix1-ix0+1;
 gfp.isy=isy=iy1-iy0+1;
 
 inum=(isx>isy?isx:isy)+1;
 if ( eevals==NULL || inum>neeval )
  {	neeval=inum;
	eevals=(eeval *)realloc(eevals,sizeof(eeval)*neeval);
  }
 gfp.eevals=eevals;
 
 afp[0]=loc->gamp,afp[1]=loc->gbg;
 afp[2]=loc->gcx,afp[3]=loc->gcy;
 afp[4]=shp->gs;

 gfp.x0=afp[2],gfp.y0=afp[3];
 gfp.gs=afp[4]+1.0;
 gfp.gd=0.0;

 gfp.derivative_level=level;

 if ( level <= 0 )
  {	lin_fit(fitpnt,yvals,afp,NULL,g2d_nmom,2,nipoint,&gfp,NULL);	}
 else
  {	if ( level==1 )	nfit=4;
	else		nfit=5;
	lam=0.001;
	for ( i=0 ; i<sfp->iter_symmetric && isfinite(afp[0]) ; i++ )
		lam=nlm_fit_base(fitpnt,yvals,afp,NULL,g2d_nmom,nfit,nipoint,&gfp,lam,10.0);
  }

 for ( i=0 ; i<5 ; i++ )
  {	if ( ! isfinite(afp[i]) )	return(1);	}

 loc->gamp=afp[0],loc->gbg=afp[1];
 loc->gcx=afp[2],loc->gcy=afp[3];

 gs=afp[4];
 if ( gs<=0.0 )	return(1);

 shp->model=SHAPE_GAUSS;
 shp->order=0;
 shp->gs=gs;
 shp->gd=shp->gk=shp->gl=0.0;
 
 return(0);
}

int fit_model_gauss(int nipoint,double *yvals,void **fitpnt,candidate *wc,
	starlocation *loc,starshape *shp,starfit *sfp)
{
 int	ret;

 loc->gamp=wc->amp,
 loc->gbg =wc->bg;
 loc->gcx =wc->cx;
 loc->gcy =wc->cy;

 shp->gs=(wc->sxx+wc->syy)/2.0;

 ret=refine_fit_model_gauss(nipoint,yvals,fitpnt,loc,shp,sfp,2);

 return(ret);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* level: 				*/
/*    -	0: only A and B			*/
/*    -	1: A, B + position (X, Y)	*/
/*    -	2: A, B, X, Y + shape (S, D, K)	*/

int refine_fit_model_elliptic_gauss(int nipoint,double *yvals,void **fitpnt,
	starlocation *loc,starshape *shp,starfit *sfp,int level)
{
 int		i,ix,iy;
 int		ix0,iy0,ix1,iy1,isx,isy,inum,nfit;
 double		lam,afp[7],gs,gd,gk;
 void		(*g2d_wmom)(void *,double *,double *,double *,void *);
 gfparam	gfp;
 static	eeval	*eevals=NULL;
 static int	neeval=0;

 g2d_wmom=gauss_2d_wmom_exact_init;

 if ( level<=0 )	nfit=2;
 else if ( level==1 )	nfit=4;
 else			nfit=7;
	
 if ( nipoint < nfit )	return(1);

 ix0=ix1=((ipoint *)fitpnt[0])->x,
 iy0=iy1=((ipoint *)fitpnt[0])->y;
 for ( i=1 ; i<nipoint ; i++ )
  {	ix=((ipoint *)fitpnt[i])->x,
	iy=((ipoint *)fitpnt[i])->y;
	if ( ix<ix0 )	ix0=ix;	
	if ( ix>ix1 )	ix1=ix;
	if ( iy<iy0 )	iy0=iy; 
	if ( iy>iy1 )	iy1=iy;
  }

 gfp.ix0=ix0;
 gfp.iy0=iy0;
 gfp.isx=isx=ix1-ix0+1;
 gfp.isy=isy=iy1-iy0+1;
 
 inum=(isx>isy?isx:isy)+1;
 if ( eevals==NULL || inum>neeval )
  {	neeval=inum;
	eevals=(eeval *)realloc(eevals,sizeof(eeval)*neeval);
  }
 gfp.eevals=eevals;

 afp[0]=loc->gamp,afp[1]=loc->gbg;
 afp[2]=loc->gcx, afp[3]=loc->gcy;
 afp[4]=shp->gs;
 afp[5]=shp->gd;
 afp[6]=shp->gk;

 gfp.x0=afp[2],gfp.y0=afp[3];
 gfp.gs=afp[4]+1.0;
 gfp.gd=0.0; 

 gfp.derivative_level=level;

 if ( level<=0 )
  {	lin_fit(fitpnt,yvals,afp,NULL,g2d_wmom,2,nipoint,&gfp,NULL);	}
 else
  {	lam=0.001;
/*	for ( i=0 ; i<sfp->iter_symmetric && isfinite(afp[0]); i++ )
		lam=nlm_fit_base(fitpnt,yvals,afp,NULL,g2d_nmom,(nfit<5?nfit:5),nipoint,&gfp,lam,10.0);
*/
	for ( i=0 ; i<sfp->iter_general && isfinite(afp[0]) ; i++ )
		lam=nlm_fit_base(fitpnt,yvals,afp,NULL,g2d_wmom,nfit,nipoint,&gfp,lam,10.0);
  }

 for ( i=0 ; i<7 ; i++ )
  {	if ( ! isfinite(afp[i]) )	return(1);	}

 loc->gamp=afp[0],loc->gbg=afp[1];
 loc->gcx=afp[2],loc->gcy=afp[3];

 gs=afp[4],gd=afp[5],gk=afp[6];
 if ( gs<=0.0 || gs*gs-gd*gd-gk*gk<=0.0 )	return(1);

 shp->model=SHAPE_ELLIPTIC;
 shp->order=0;
 shp->gs=gs;
 shp->gd=gd,
 shp->gk=gk;
 shp->gl=0.0;

 return(0);
}

int fit_model_elliptic_gauss(int nipoint,double *yvals,void **fitpnt,candidate *wc,
	starlocation *loc,starshape *shp,starfit *sfp)
{
 int	ret;

 loc->gamp=wc->amp,
 loc->gbg =wc->bg;
 loc->gcx =wc->cx;
 loc->gcy =wc->cy;

 shp->gs=(wc->sxx+wc->syy)/2.0;
 shp->gd=0.0;
 shp->gk=0.0;

 ret=refine_fit_model_elliptic_gauss(nipoint,yvals,fitpnt,loc,shp,sfp,2);

 return(ret);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int refine_fit_model_deviated(int nipoint,double *yvals,void **fitpnt,
	starlocation *loc,starshape *shp,starfit *sfp,int level)
{
 return(0);
}

int fit_model_deviated(int nipoint,double *yvals,void **fitpnt,candidate *wc,
	starlocation *loc,starshape *shp,starfit *sfp,int order)
{
 int	i;
 double	lam,afp[MAX_DEVIATION_COEFF+5],gs,wgs,wgd,wgk,ndev,ndev0,ntot,*dcf;
 void	(*g2d_nmom)(void *,double *,double *,double *,void *);
 void	(*g2d_wdev)(void *,double *,double *,double *,void *);

 g2d_nmom=gauss_2d_nmom_exact;
 g2d_wdev=gauss_2d_wdev_exact;

 ndev0=(order+1)*(order+2)/2;
 ndev=ndev0-3;
 ntot=ndev+5;

 if ( nipoint < ntot )	return(1);

 lam=0.001;
 afp[0]=wc->amp,afp[1]=wc->bg;
 afp[2]=wc->cx ,afp[3]=wc->cy;
 afp[4]=(wc->sxx+wc->syy)/2.0;

 for ( i=0 ; i<sfp->iter_symmetric ; i++ )
	lam=nlm_fit_base(fitpnt,yvals,afp,NULL,g2d_nmom,5,nipoint,NULL,lam,10.0);

 for ( i=0 ; i<5 ; i++ )
  {	if ( ! isfinite(afp[i]) )	return(1);	}

 dcf=&afp[5];
 for ( i=0 ; i<ndev ; i++ )
  {	dcf[i]=0.0;		}
 wdev_order=order;
 for ( i=0 ; i<sfp->iter_general ; i++ )
	lam=nlm_fit_base(fitpnt,yvals,afp,NULL,g2d_wdev,ntot,nipoint,NULL,lam,10.0);

 for ( i=0 ; i<ntot ; i++ )
  {	if ( ! isfinite(afp[i]) )	return(1);	}

 loc->gamp=afp[0],loc->gbg=afp[1];
 loc->gcx =afp[2],loc->gcy=afp[3];

 gs=afp[4];
 wgs=afp[4]-(afp[5]+afp[7])*0.5;
 wgd=-(afp[5]-afp[7])*0.5;
 wgk=-(afp[6]);
 if ( gs <= 0.0 || wgs<=0.0 || wgs*wgs-wgd*wgd-wgk*wgk<=0.0 )	return(1);

 shp->model=SHAPE_DEVIATED;
 shp->order=order;
 shp->gs=gs;
 shp->gd=wgd,
 shp->gk=wgk;
 shp->gl=0.0;
 for ( i=0 ; i<ndev ; i++ )
  {	shp->mom[i]=dcf[i];		}
 for ( ; i<MAX_DEVIATION_COEFF ; i++ )
  {	shp->mom[i]=0.0;		}

 return(0);
}

/*****************************************************************************/

int fit_model_combined(int nipoint,double *yvals,void **fitpnt,candidate *wc,
	starlocation *loc,starshape *mshps,int n,
	starmodelfit *smfs,starfit *sfp)
{
 int		i,j,nshape,nfv,is_failed;
 double		lam,*afp,bg,gs,gd,gk,gg;
 modelparam	*mps;
 void		(*funct)(void *,double *,double *,double *,void *);

 mps=(modelparam *)malloc(sizeof(modelparam)*(n+1));

 nfv=3;
 for ( i=0 ; i<n ; i++ )
  {	switch ( smfs[i].model )
	 {   case SHAPE_GAUSS:
		funct=gauss_2d_nmom_exact;
		nshape=1;
		break;
	     case SHAPE_ELLIPTIC:
		funct=gauss_2d_wmom_exact;
		nshape=3;
		break;
	     case SHAPE_DEVIATED:
		funct=gauss_2d_wdev_exact;
		nshape=1+(smfs[i].order+1)+(smfs[i].order+2)/2-3;
		break;
	     default:
		funct=NULL;
		nshape=0;
	 }
	if ( funct==NULL || nshape<=0 )
	 {	free(mps);return(1);		}

	nfv+=1+nshape;		/* aplidude/factor + shape parameters */
	
	mps[i].funct=funct;
	mps[i].nshape=nshape;
	mps[i].param=NULL;
  }


 mps[n].funct=NULL;
 mps[n].nshape=0;
 mps[n].param=NULL;

 if ( nipoint < nfv )
  {	free(mps);return(1); 		}

 funct=gauss_2d_nmom_exact;

 afp=(double *)malloc(sizeof(double)*nfv);

 lam=0.001;
 
 afp[0]=wc->amp,afp[1]=wc->bg;
 afp[2]=wc->cx,afp[3]=wc->cy;
 afp[4]=(wc->sxx+wc->syy)/2.0;

 for ( i=0 ; i<sfp->iter_symmetric ; i++ )
 	lam=nlm_fit_base(fitpnt,yvals,afp,NULL,funct,5,nipoint,NULL,lam,10.0);

 bg=afp[1];
 afp[1]=afp[2];
 afp[2]=afp[3];
 afp[3]=afp[0];
 afp[0]=bg;
 
 nfv=3;
 for ( i=0 ; i<n ; i++ )
  {	switch ( smfs[i].model )
	 {   case SHAPE_GAUSS:
		nshape=1;
		break;
	     case SHAPE_ELLIPTIC:
		nshape=3;
		break;
	     case SHAPE_DEVIATED:
		nshape=1+(smfs[i].order+1)+(smfs[i].order+2)/2-3;
		break;
	     default:
		nshape=0;
	 }

	if ( i>0 )
	 {	afp[0+nfv]=0.1*afp[3];
		afp[1+nfv]=0.5*afp[4];
	 }

	nfv+=1+nshape;
  }

/* fprintf(stderr,"nfv=%d\n",nfv);
 for ( i=0 ; i<nfv ; i++ )
  {	fprintf(stderr,"%12g\n",afp[i]);	} */

 if ( ! ( n==1 && smfs[0].model==SHAPE_GAUSS ) )
  {	for ( i=0 ; i<sfp->iter_general ; i++ )
		lam=nlm_fit_base(fitpnt,yvals,afp,NULL,model_combine,nfv,nipoint,mps,lam,10.0);
  }

/* fprintf(stderr,"nfv=%d\n",nfv);
 for ( i=0 ; i<nfv ; i++ )
  {	fprintf(stderr,"%12g\n",afp[i]);	} */

 loc->gbg=afp[0];
 loc->gcx=afp[1];
 loc->gcy=afp[2];
 loc->gamp=0.0;

 nfv=3;
 is_failed=0;

 for ( i=0 ; i<n ; i++ )
  {	switch ( smfs[i].model )
	 {   case SHAPE_GAUSS:
		mshps[i].model=SHAPE_GAUSS;
		mshps[i].order=0;
		mshps[i].gs=afp[1+nfv];
		if ( i==0 && ( mshps[i].gs <= 0.0 ) )
			is_failed=1;
		mshps[i].gd=mshps[i].gk=mshps[i].gl=0.0;
		for ( j=0 ; j<MAX_DEVIATION_COEFF ; j++ )
		 {	mshps[i].mom[j]=0.0;			}
		mshps[i].factor=afp[0+nfv];
		nshape=1;
		break;
	     case SHAPE_ELLIPTIC:
		mshps[i].model=SHAPE_ELLIPTIC;
		mshps[i].order=0;
		gs=mshps[i].gs=afp[1+nfv];
		gd=mshps[i].gd=afp[2+nfv];
		gk=mshps[i].gk=afp[3+nfv];
		if ( i==0 )
		 {	gg=gs*gs-gd*gd-gk*gk;
			if ( gs<=0.0 || gg<=0.0 )
				is_failed=1;		
		 }
		mshps[i].gl=0.0;
		for ( j=0 ; j<MAX_DEVIATION_COEFF ; j++ )
		 {	mshps[i].mom[j]=0.0;			}
		mshps[i].factor=afp[0+nfv];
		nshape=3;
		break;
	     case SHAPE_DEVIATED:
		mshps[i].model=SHAPE_DEVIATED;
		mshps[i].order=smfs[i].order;
		mshps[i].gs=afp[1+nfv];
		gs=afp[1+nfv]-(afp[2+nfv]+afp[4+nfv])*0.5;
		gd=mshps[i].gd=-(afp[2+nfv]-afp[4+nfv])*0.5;
		gk=mshps[i].gk=-(afp[3+nfv]);
		if ( i==0 )
		 {	gg=gs*gs-gd*gd-gk*gk;
			if ( mshps[i].gs<=0.0 || gs<=0.0 || gg<=0.0 )	
				is_failed=1;		
		 }
		mshps[i].gl=0.0;
		nshape=1+(smfs[i].order+1)+(smfs[i].order+2)/2-3;
		for ( j=0 ; j<nshape-1 ; j++ )
		 {	mshps[i].mom[j]=afp[2+nfv+j];			}
		for ( ; j<MAX_DEVIATION_COEFF ; j++ )
		 {	mshps[i].mom[j]=0.0;				}
		mshps[i].factor=afp[0+nfv];
		break;
	     default:
		funct=NULL;
		nshape=0;
	 }

	loc->gamp += mshps[i].factor;
	nfv+=1+nshape;
  }

 if ( ! is_failed )
  {	for ( i=0 ; i<n ; i++ )
	 {	mshps[i].factor /= loc->gamp;		}
  }

 free(afp);
 free(mps);

 return(is_failed); 
 
}

/*****************************************************************************/

int fit_star_single_model(fitsimage *img,char **mask,candidate *cands,int ncand,
	star **rstars,int *rnstar,starfit *sfp,int shapemodel,int modelorder)
{
 int		i,n,bsize;
 ipoint  	*ipoints;
 void		**fitpnt;double *yvals;

 star		*stars,*ws;
 candidate	*wc;
 int		nstar;

 starshape	shp;	
 starlocation	loc;

 if ( rstars != NULL )	*rstars=NULL;
 if ( rnstar != NULL )	*rnstar=0;
 
 if ( ncand == 0 )	return(0);
 else if ( ncand < 0 )	return(1);

 bsize=0;
 for ( n=0 ; n<ncand ; n++ )
  {	if ( cands[n].nipoint>bsize )
		bsize=cands[n].nipoint;
  }

 fitpnt=(void  **)malloc(sizeof(void *)*bsize);
 yvals =(double *)malloc(sizeof(double)*bsize);

 stars=NULL;
 nstar=0;
 for ( n=0 ; n<ncand ; n++ )
  {	wc=&cands[n];

	if ( wc->marked )				continue;
	if ( wc->nipoint==0 || wc->ipoints==NULL )	continue;

	ipoints=wc->ipoints;
	for ( i=0 ; i<wc->nipoint ; i++ )
	 {	int	ix,iy;
		ix=wc->ipoints[i].x,
		iy=wc->ipoints[i].y;
		yvals[i]=img->data[iy][ix];
		fitpnt[i]=(void *)(&ipoints[i]);
	 }

	switch ( shapemodel )
	 {   case SHAPE_GAUSS:
		i=fit_model_gauss(wc->nipoint,yvals,fitpnt,wc,&loc,&shp,sfp);
		break;
	     case SHAPE_ELLIPTIC:
		i=fit_model_elliptic_gauss(wc->nipoint,yvals,fitpnt,wc,&loc,&shp,sfp);
		break;
	     case SHAPE_DEVIATED:
		i=fit_model_deviated(wc->nipoint,yvals,fitpnt,wc,&loc,&shp,sfp,modelorder);
		break;
	     default:
		i=-1;
		break;
	 }

	if ( i )	continue;

	stars=(star *)realloc(stars,sizeof(star)*(nstar+1));
	ws=&stars[nstar];
	nstar++;

	memcpy(&ws->location,&loc,sizeof(starlocation));
	memcpy(&ws->shape,&shp,sizeof(starshape));

	#ifdef	STAR_MULTIMODEL
	ws->mshapes=NULL;
	ws->nmshape=1;
	#endif

	star_set_common_shape_params(shp.gs,shp.gd,shp.gk,ws);

	ws->flux=loc.gamp*star_get_unity_flux(&shp);

	ws->marked=0;
	ws->cand=wc;
  }

 free(yvals);
 free(fitpnt); 

 if ( rstars != NULL )	*rstars=stars;
 if ( rnstar != NULL )	*rnstar=nstar;

 return(0);
}

/*****************************************************************************/

#ifdef	STAR_MULTIMODEL
int fit_star_multi_models(fitsimage *img,char **mask,candidate *cands,int ncand,
	star **rstars,int *rnstar,starfit *sfp,starmodelfit *smfs,int nsmf)
{
 int		i,sx,sy,n,bsize,nmshp;
 ipoint  	*ipoints;
 void		**fitpnt;double *yvals;
 double		gs,gd,gk,gg,gm;

 star		*stars,*ws;
 candidate	*wc;
 int		nstar;

 starshape	*mshps,*sw;	
 starlocation	loc;

 if ( rstars != NULL )	*rstars=NULL;
 if ( rnstar != NULL )	*rnstar=0;
 
 if ( ncand == 0 )	return(0);
 else if ( ncand < 0 )	return(1);

 if ( nsmf<=0 )		return(1);

 sx=img->sx,sy=img->sy;

 bsize=0;
 for ( n=0 ; n<ncand ; n++ )
  {	if ( cands[n].nipoint>bsize )
		bsize=cands[n].nipoint;
  }

 fitpnt=(void  **)malloc(sizeof(void *)*bsize);
 yvals =(double *)malloc(sizeof(double)*bsize);

 stars=NULL;
 nstar=0;
 for ( n=0 ; n<ncand ; n++ )
  {	wc=&cands[n];

	if ( wc->marked )				continue;
	if ( wc->nipoint==0 || wc->ipoints==NULL )	continue;

	ipoints=wc->ipoints;
	for ( i=0 ; i<wc->nipoint ; i++ )
	 {	int	ix,iy;
		ix=wc->ipoints[i].x,
		iy=wc->ipoints[i].y;
		yvals[i]=img->data[iy][ix];
		fitpnt[i]=(void *)(&ipoints[i]);
	 }

	nmshp=nsmf;
	mshps=(starshape *)malloc(sizeof(starshape)*nmshp);

	i=fit_model_combined(wc->nipoint,yvals,fitpnt,wc,&loc,mshps,nmshp,smfs,sfp);

	if ( i )
	 {	free(mshps);continue;	}	/* fit failed */

	gs=mshps[0].gs,	/* mshps[0] is assumed to be the dominant model */
	gd=mshps[0].gd,
	gk=mshps[0].gk;
	
	gg=gs*gs-gd*gd-gk*gk;
	if ( gg<=0.0 )
	 {	free(mshps);continue;	}	/* fit failed */
	gm=sqrt(gg);

	stars=(star *)realloc(stars,sizeof(star)*(nstar+1));
	ws=&stars[nstar];
	nstar++;

	memcpy(&ws->location,&loc,sizeof(starlocation));

	if ( nmshp==1 )		/* only one model... */
	 {	memcpy(&ws->shape,mshps,sizeof(starshape));
		ws->mshapes=NULL;
		ws->nmshape=1;
		free(mshps);
	 }
	else			/* ...combined model */
	 {	ws->mshapes=mshps;
		ws->nmshape=nmshp;
	 }

	star_set_common_shape_params(gs,gd,gk,ws);

	ws->flux=0;
	for ( i=0 ; i<nmshp ; i++ )	/* calculate the cumulative flux */
	 {	if ( ws->mshapes != NULL )	sw=&ws->mshapes[i];
		else				sw=&ws->shape;
		ws->flux += loc.gamp*sw->factor*star_get_unity_flux(sw);
	 }

	ws->marked=0;
	ws->cand=wc;
  }

 free(yvals);
 free(fitpnt); 

 if ( rstars != NULL )	*rstars=stars;
 if ( rnstar != NULL )	*rnstar=nstar;

 return(0);
}

/*
int fit_star_model(fitsimage *img,char **mask,candidate *cands,int ncand,
	star **rstars,int *rnstar,starfit *sfp,int shapemodel,int modelorder)
{
 starmodelfit	smf;
 
 smf.model=shapemodel;
 smf.order=modelorder;

 fit_star_multi_models(img,mask,cands,ncand,rstars,rnstar,sfp,&smf,1);

 return(0);
}
*/
#endif

/*****************************************************************************/

int drawback_model_gauss(ipoint *ipoints,int nipoint,double *yvals,starlocation *loc,starshape *shp,double mul)
{
 int		i;
 int		ix0,iy0,ix1,iy1,isx,isy,ix,iy,inum;
 double		afp[5],y;
 void		(*g2d_nmom)(void *,double *,double *,double *,void *);
 gfparam	gfp;
 static	eeval	*eevals=NULL;
 static int	neeval=0;

 g2d_nmom=gauss_2d_nmom_exact_init;

 ix0=ix1=ipoints[0].x,
 iy0=iy1=ipoints[0].y;
 for ( i=1 ; i<nipoint ; i++ )
  {	ix=ipoints[i].x,
	iy=ipoints[i].y;
	if ( ix<ix0 )	ix0=ix;	
	if ( ix>ix1 )	ix1=ix;
	if ( iy<iy0 )	iy0=iy; 
	if ( iy>iy1 )	iy1=iy;
  }

 gfp.ix0=ix0;
 gfp.iy0=iy0;
 gfp.isx=isx=ix1-ix0+1;
 gfp.isy=isy=iy1-iy0+1;
 
 inum=(isx>isy?isx:isy)+1;
 if ( eevals==NULL || inum>neeval )
  {	neeval=inum;
	eevals=(eeval *)realloc(eevals,sizeof(eeval)*neeval);
  }
 gfp.eevals=eevals;
 
 afp[0]=loc->gamp,afp[1]=0.0;
 afp[2]=loc->gcx,afp[3]=loc->gcy;
 afp[4]=shp->gs;

 gfp.x0=afp[2],gfp.y0=afp[3];
 gfp.gs=afp[4]+1.0;
 gfp.gd=0.0;
 gfp.derivative_level=-1;	/* no derivatives required for drawback */

 for ( i=0 ; i<nipoint ; i++ )
  {	/* call callback directly... */
	g2d_nmom((void *)(&ipoints[i]),afp,&y,NULL,&gfp);
	yvals[i] += mul * y;
  }
 return(0);
}
int drawback_model_elliptic_gauss(ipoint *ipoints,int nipoint,double *yvals,starlocation *loc,starshape *shp,double mul)
{
 int		i;
 int		ix0,iy0,ix1,iy1,isx,isy,ix,iy,inum;
 double		afp[7],y;
 void		(*g2d_wmom)(void *,double *,double *,double *,void *);
 gfparam	gfp;
 static	eeval	*eevals=NULL;
 static int	neeval=0;

 g2d_wmom=gauss_2d_wmom_exact_init;

 ix0=ix1=ipoints[0].x,
 iy0=iy1=ipoints[0].y;
 for ( i=1 ; i<nipoint ; i++ )
  {	ix=ipoints[i].x,
	iy=ipoints[i].y;
	if ( ix<ix0 )	ix0=ix;	
	if ( ix>ix1 )	ix1=ix;
	if ( iy<iy0 )	iy0=iy; 
	if ( iy>iy1 )	iy1=iy;
  }

 gfp.ix0=ix0;
 gfp.iy0=iy0;
 gfp.isx=isx=ix1-ix0+1;
 gfp.isy=isy=iy1-iy0+1;
 
 inum=(isx>isy?isx:isy)+1;
 if ( eevals==NULL || inum>neeval )
  {	neeval=inum;
	eevals=(eeval *)realloc(eevals,sizeof(eeval)*neeval);
  }
 gfp.eevals=eevals;
 
 afp[0]=loc->gamp,afp[1]=0.0;
 afp[2]=loc->gcx,afp[3]=loc->gcy;
 afp[4]=shp->gs;
 afp[5]=shp->gd;
 afp[6]=shp->gk;

 gfp.x0=afp[2],gfp.y0=afp[3];
 gfp.gs=afp[4]+1.0;
 gfp.gd=0.0;
 gfp.derivative_level=-1;	/* no derivatives required for drawback */

 for ( i=0 ; i<nipoint ; i++ )
  {	/* call callback directly... */
	g2d_wmom((void *)(&ipoints[i]),afp,&y,NULL,&gfp);
	yvals[i] += mul * y;
  }
 return(0);
}
int drawback_model_deviated(ipoint *ipoints,int nipoint,double *yvals,starlocation *loc,starshape *shp,double mul)
{
 return(0);
}

int drawback_model(ipoint *ipoints,int nipoint,double *yvals,starlocation *loc,starshape *shp,double mul)
{
 switch ( shp->model )
  {  case SHAPE_GAUSS:
	drawback_model_gauss(ipoints,nipoint,yvals,loc,shp,mul);
	break;
     case SHAPE_ELLIPTIC:
	drawback_model_elliptic_gauss(ipoints,nipoint,yvals,loc,shp,mul);
	break;
     case SHAPE_DEVIATED:
	drawback_model_deviated(ipoints,nipoint,yvals,loc,shp,mul);
	break;
  }
 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int collective_fit_star_single_model(fitsimage *img,char **mask,
	star *stars,int nstar,ipointlist *ipl,starfit *sfp,int is_putback,int level)
{
 int		i,n,bsize,nipoint;
 ipoint  	*ipoints;
 void		**fitpnt;
 double		*yvals;

 star		*ws;

 starshape	shp;	
 starlocation	loc;

 bsize=0;
 for ( n=0 ; n<nstar ; n++ )
  {	if ( ipl[n].nipoint>bsize )
		bsize=ipl[n].nipoint;
  }

 fitpnt=(void  **)malloc(sizeof(void *)*bsize);
 yvals =(double *)malloc(sizeof(double)*bsize);

 for ( n=0 ; n<nstar ; n++ )
  {	ws=&stars[n];
	ipoints=ipl[n].ipoints;
	nipoint=ipl[n].nipoint;
	
	if ( ws->marked )		continue;

	for ( i=0 ; i<nipoint ; i++ )
	 {	int	ix,iy;
		ix=ipoints[i].x,
		iy=ipoints[i].y;
		yvals[i]=img->data[iy][ix];
		fitpnt[i]=(void *)(&ipoints[i]);
	 }
	
	if ( is_putback )
	 {	drawback_model(ipoints,nipoint,yvals,&ws->location,&ws->shape,+1);	}

	memcpy(&loc,&ws->location,sizeof(starlocation));
	memcpy(&shp,&ws->shape,sizeof(starshape));

	switch ( shp.model )
	 {   case SHAPE_GAUSS:
		i=refine_fit_model_gauss(nipoint,yvals,fitpnt,&loc,&shp,sfp,level);
		break;
	     case SHAPE_ELLIPTIC:
		i=refine_fit_model_elliptic_gauss(nipoint,yvals,fitpnt,&loc,&shp,sfp,level);
		break;
	     case SHAPE_DEVIATED:
		i=refine_fit_model_deviated(nipoint,yvals,fitpnt,&loc,&shp,sfp,level);
		break;
	     default:
		i=-1;
		break;
	 }

	/* if the fit was successful, overwrite the old star parameters
	   with the refined ones,... */
	if ( ! i && loc.gamp>0 )
	 {	memcpy(&ws->location,&loc,sizeof(starlocation));
		memcpy(&ws->shape,&shp,sizeof(starshape));
	 }
	/* ... otherwise, skip the star... */
	else
		ws->marked=1;
	 
  }

 free(yvals);
 free(fitpnt); 

 return(0);
}

int collective_fit_star_single_model_iterative(fitsimage *img,char **mask,
	star *stars,int nstar,ipointlist *ipl,starfit *sfp,int level,int niter)
{
 int		i,j,ix,iy,sx,sy,n,k,bsize,nipoint;
 ipoint  	*ipoints;
 double		*yvals;
 fitsimage	fitimg_data,*fitimg;
 star		*ws;

 sx=img->sx,sy=img->sy;

 bsize=0;
 for ( n=0 ; n<nstar ; n++ )
  {	if ( ipl[n].nipoint>bsize )
		bsize=ipl[n].nipoint;
  }

 yvals =(double *)malloc(sizeof(double)*bsize);

 collective_fit_star_single_model(img,mask,stars,nstar,ipl,sfp,0,0);

 if ( niter>0 )
  {	fitimg=&fitimg_data;
	fits_image_duplicate(fitimg,img,1);

	for ( k=0 ; k<niter ; k++ )
	 {	for ( n=0 ; n<nstar ; n++ )
		 {	ws=&stars[n];
			ipoints=ipl[n].ipoints;
			nipoint=ipl[n].nipoint;
			if ( ws->marked )		continue;
			for ( i=0 ; i<nipoint ; i++ )
			 {	ix=ipoints[i].x,
				iy=ipoints[i].y;
				yvals[i]=fitimg->data[iy][ix];
			 }
			drawback_model(ipoints,nipoint,yvals,&ws->location,&ws->shape,-1);
			for ( i=0 ; i<nipoint ; i++ )
			 {	ix=ipoints[i].x,
				iy=ipoints[i].y;
				fitimg->data[iy][ix]=yvals[i];
			 }
		 }
		collective_fit_star_single_model(fitimg,mask,stars,nstar,ipl,sfp,1,level);
		for ( i=0 ; i<sy ; i++ )
		 {	for ( j=0 ; j<sx ; j++ )
			 {	fitimg->data[i][j]=img->data[i][j];	}
		 }
	 }

	fits_image_free(fitimg);
  }

 for ( n=0 ; n<nstar ; n++ )
  {	ws=&stars[n];
	star_set_common_shape_params(ws->shape.gs,ws->shape.gd,ws->shape.gk,ws);
	ws->flux=ws->location.gamp*star_get_unity_flux(&ws->shape);
  }

 free(yvals);

 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static int index_compare_star_sort_x(int p1,int p2,void *param)
{
 star	*stars=(star *)param;
 if ( stars[p1].location.gcx < stars[p2].location.gcx )		return(-1);
 else if ( stars[p1].location.gcx > stars[p2].location.gcx )	return(1);
 else								return(0);
}
static int index_compare_star_sort_y(int p1,int p2,void *param)
{
 star	*stars=(star *)param;
 if ( stars[p1].location.gcy < stars[p2].location.gcy )		return(-1);
 else if ( stars[p1].location.gcy > stars[p2].location.gcy )	return(1);
 else								return(0);
}
static int index_compare_star_value_x(int p,double x,void *param)
{
 star	*stars=(star *)param;
 if ( stars[p].location.gcx < x )	return(-1);
 else if ( stars[p].location.gcx > x )	return(1);
 else					return(0);
}
static int index_compare_star_value_y(int p,double y,void *param)
{
 star	*stars=(star *)param;
 if ( stars[p].location.gcy < y )	return(-1);
 else if ( stars[p].location.gcy > y )	return(1);
 else					return(0);
}

typedef struct
 {	int	ndat;
	double	*bgs;
	double	*amps;
 } sfdata;

int collective_fit_star_single_model_blocked(fitsimage *img,char **mask,
	star *stars,int nstar,ipointlist *ipl,double bhsize)
{
 int		i,j,k,n,asize,nipoint,ix0,ix1,iy0,iy1,isx,isy,ix,iy;
 int		p,q,r,t,c;
 ipoint  	*ipoints;
 double		*yvaldata,**ypl,x0,x1,y0,y1,**amatrix,*bvector,*fvars;
 star		*ws;
 multiindex	mi;
 int		*ixs,***idmx,*ipx;
 int		nix,mnx,nvar;
 sfdata		*sfs;

 multiindex_create(stars,nstar,
	index_compare_star_sort_x,
	index_compare_star_sort_y,
	&mi);

 sfs=(sfdata *)malloc(sizeof(sfdata)*nstar);
 for ( n=0 ; n<nstar ; n++ )
  {	sfs[n].ndat=0;
	sfs[n].bgs=NULL;
	sfs[n].amps=NULL;
  }

 asize=0;
 for ( n=0 ; n<nstar ; n++ )
  {	asize+=ipl[n].nipoint;		}
 yvaldata=(double *)malloc(sizeof(double)*asize);
 for ( i=0 ; i<asize ; i++ )
  {	yvaldata[i]=0.0;		}
 ypl=(double **)malloc(sizeof(double *)*nstar);
 asize=0;
 for ( n=0 ; n<nstar ; n++ )
  {	ypl[n]=yvaldata+asize;
	asize+=ipl[n].nipoint;
  }

 for ( n=0 ; n<nstar ; n++ )
  {	ws=&stars[n];
	ipoints=ipl[n].ipoints;
	nipoint=ipl[n].nipoint;
	ws->location.gamp=1.0;
	drawback_model(ipoints,nipoint,ypl[n],&ws->location,&ws->shape,+1.0);
  }

 mnx=64;
 ixs=(int *)malloc(sizeof(int)*mnx);
 amatrix=matrix_alloc(mnx+1);
 bvector=vector_alloc(mnx+1);
 fvars  =vector_alloc(mnx+1);
 for ( n=0 ; n<nstar ; n++ )
  {
	ws=&stars[n];

	x0=ws->location.gcx-bhsize,x1=ws->location.gcx+bhsize;
	y0=ws->location.gcy-bhsize,y1=ws->location.gcy+bhsize;
	nix=multiindex_range_query(&mi,
		index_compare_star_value_x,
		index_compare_star_value_y,
		x0,x1,y0,y1,ixs,mnx);

	while ( nix>=mnx )
	 {	vector_free(fvars);
		vector_free(bvector);
		matrix_free(amatrix);
		mnx=2*mnx;
		ixs=(int *)realloc(ixs,sizeof(int)*mnx);
		nix=multiindex_range_query(&mi,
			index_compare_star_value_x,
			index_compare_star_value_y,
			x0,x1,y0,y1,ixs,mnx);
		amatrix=matrix_alloc(mnx+1);
		bvector=vector_alloc(mnx+1);
		fvars  =vector_alloc(mnx+1);
	 };
	nvar=nix+1;

	for ( i=0 ; i<nix ; i++ )
	 {	if ( ixs[i]==n )	break;	}
	if ( i==nix || nix==0 )	
	 {	fprintf(stderr,"fatal: n=%d\n",n);	}

	for ( i=0 ; i<nvar ; i++ )
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	amatrix[i][j]=0.0;		}
		bvector[i]=0.0;
	 }

	ix0=ix1=ipl[n].ipoints[0].x;
	iy0=iy1=ipl[n].ipoints[0].y;
	for ( k=0 ; k<nix ; k++ )
	 {	c=ixs[k];
		ipoints=ipl[c].ipoints;
		nipoint=ipl[c].nipoint;
		for ( i=0 ; i<nipoint ; i++ )
		 {	if ( ipoints[i].x < ix0 )	ix0=ipoints[i].x;
			if ( ipoints[i].x > ix1 )	ix1=ipoints[i].x;
			if ( ipoints[i].y < iy0 )	iy0=ipoints[i].y;
			if ( ipoints[i].y > iy1 )	iy1=ipoints[i].y;
		 }
	 }
	ix1++,iy1++;
	isx=ix1-ix0;
	isy=iy1-iy0;
	idmx=tensor_alloc_3d(int,2*nix+1,isx,isy);
	for ( i=0 ; i<isy ; i++ )
	 {	for ( j=0 ; j<isx ; j++ )
		 {	idmx[i][j][0]=0;		}
	 }
	for ( k=0 ; k<nix ; k++ )
	 {	c=ixs[k];
		ipoints=ipl[c].ipoints;
		nipoint=ipl[c].nipoint;
		for ( i=0 ; i<nipoint ; i++ )
		 {	ix=ipoints[i].x-ix0;
			iy=ipoints[i].y-iy0;
			ipx=idmx[iy][ix];
			j=ipx[0];
			ipx[2*j+1]=k;
			ipx[2*j+2]=i;
			ipx[0]++;
		 }
	 }
	t=0;
	for ( i=0 ; i<isy ; i++ )
	 {	for ( j=0 ; j<isx ; j++ )
		 {	ipx=idmx[i][j];
			r=ipx[0];
			if ( r<=0 )	  /* comment this if you want all    */
				continue; /* background data to be fitted... */	
			fvars[0]=1.0;
			for ( k=1 ; k<nvar ; k++ )
			 {	fvars[k]=0.0;		}
			for ( p=0 ; p<r ; p++ )
			 {	k=ipx[2*p+1];
				c=ixs[k];
				fvars[k+1]=ypl[c][ipx[2*p+2]];
			 }
			for ( p=0 ; p<nvar ; p++ )
			 {	if ( fvars[p]==0.0 )	continue;
				for ( q=0 ; q<nvar ; q++ )
				 {	if ( fvars[q]==0.0 )	continue;
					amatrix[p][q]+=fvars[p]*fvars[q];
				 }
				bvector[p]+=fvars[p]*img->data[iy0+i][ix0+j];
			 }
			t++;
		 }
	 }
	/*fprintf(stderr,"nvar=%d t=%d isx=%d isy=%d\n",nvar,t,isx,isy);*/
	solve_gauss(amatrix,bvector,nvar);
	for ( k=0 ; k<nix ; k++ )
	 {	c=ixs[k];
		r=sfs[c].ndat;
		sfs[c].bgs =(double *)realloc(sfs[c].bgs ,sizeof(double)*(r+1));
		sfs[c].amps=(double *)realloc(sfs[c].amps,sizeof(double)*(r+1));
		sfs[c].bgs [r]=bvector[0];
		sfs[c].amps[r]=bvector[k+1];
		sfs[c].ndat++;
	 }
	
	/*
	ws->location.gbg=bvector[0];
	for ( k=0 ; k<nix ; k++ )
	 {	if ( ixs[k]==n )
			ws->location.gamp=bvector[k+1];
	 }
	*/

	tensor_free(idmx);
  }
 vector_free(fvars);
 vector_free(bvector);
 matrix_free(amatrix);
 free(ixs);

 free(ypl);
 free(yvaldata);
 multiindex_free(&mi);

 for ( n=0 ; n<nstar ; n++ )
  {	ws=&stars[n];
	r=sfs[n].ndat;
	ws->location.gbg =median(sfs[n].bgs ,r);
	ws->location.gamp=median(sfs[n].amps,r);
	free(sfs[n].amps);
	free(sfs[n].bgs);
	star_set_common_shape_params(ws->shape.gs,ws->shape.gd,ws->shape.gk,ws);
	ws->flux=ws->location.gamp*star_get_unity_flux(&ws->shape);
  }

 free(sfs);

 return(0);
}

/*****************************************************************************/

