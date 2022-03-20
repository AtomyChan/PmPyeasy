/*****************************************************************************/
/* weight-star.c							     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Functions related to weight stamp handling:				     */
/* create weight list from star list.					     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2006; Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <stdarg.h>

#include <fits/fits.h>

#include "tensor.h"
#include "stars.h"
#include "psf.h"

#include "weight.h"

/*****************************************************************************/

int weight_draw(weightlist *wl,star *stars,int nstar,int hsize,int grid,psf *tpd)
{
 int		o,l,n,order,fsize,ix0,iy0,i,j;
 double		s,d,k,dflux,sum,zf,iz,gs,mom[MAX_DEVIATION_COEFF],*wmom,wz,x,y;
 star		*ws;
 double		***zdata,**iarr,x0,y0;
 weight		*ww;

 if ( stars==NULL || nstar<0 )	return(1);

 if ( grid<1 )	zf=1.0,grid=1;
 else		zf=(double)grid;
 iz=1.0/zf;

 fsize=grid*(2*hsize+1);
 zdata=(double ***)tensor_alloc_3d(double,fsize,fsize,nstar);
 iarr=(double **)tensor_alloc_2d(double,fsize,fsize);

 wl->zdata=zdata;
 wl->weights=(weight *)malloc(sizeof(weight)*nstar);
 wl->nweight=nstar;
 wl->hsize=hsize;
 wl->grid=grid;

 for ( n=0 ; n<nstar ; n++ )
  {	ws=&stars[n];

	dflux=ws->flux;
	if ( dflux<=0.0 )	continue;

	x0=ws->location.gcx;
	y0=ws->location.gcy;
	ix0=(int)floor(x0),
	iy0=(int)floor(y0);		
	x=(x0-ix0)+hsize;
	y=(y0-iy0)+hsize;

	if ( ws->shape.model==SHAPE_GAUSS || ws->shape.model==SHAPE_ELLIPTIC )
	 {	s=ws->gsig;
		if ( ws->shape.model==SHAPE_ELLIPTIC )
			d=ws->gdel,k=ws->gkap;
		else	
			d=k=0;

		star_draw_gauss(iarr,fsize,fsize,x*zf,y*zf,s*zf,d*zf,k*zf);
	 }

	else if ( ws->shape.model==SHAPE_DEVIATED )
	 {	if ( grid>1 )
	 	 {	wz=iz*iz;
			gs=ws->shape.gs*wz;
			order=ws->shape.order;
			for ( o=2,l=0 ; o<=order ; o++ )
			 {	for ( j=0 ; j<=o ; j++,l++ )
				 {	mom[l]=ws->shape.mom[l]*wz;	}
				wz=wz*iz;
			 }
			wmom=mom;
		 }
		else
		 {	gs=ws->shape.gs;
			order=ws->shape.order;
			wmom=ws->shape.mom;
		 }

		star_draw_deviated(iarr,fsize,fsize,x*zf,y*zf,gs,order,wmom);
	 }
	
	else if ( ws->shape.model == SHAPE_PSF )
	 {	star_draw_psf(iarr,fsize,fsize,x*zf,y*zf,tpd,
			ws->location.gcx,ws->location.gcy,
			ws->shape.gs*zf,ws->shape.gd*zf,ws->shape.gk*zf,ws->shape.gl*zf);
	 }

	sum=0.0;
	for ( i=0 ; i<fsize ; i++ )
	 {	for ( j=0 ; j<fsize ; j++ )
		 {	sum+=iarr[i][j];
			/*fprintf(stderr,"%12g\n",iarr[i][j]);*/
		 }
	 }

	for ( i=0 ; i<fsize ; i++ )
	 {	for ( j=0 ; j<fsize ; j++ )
		 {	zdata[n][i][j]=iarr[i][j]*dflux/sum;
		 }
	 }

	ww=&wl->weights[n];
	ww->x=x0,
	ww->y=y0;
	ww->ix=ix0,
	ww->iy=iy0;
	ww->flux=dflux;
	ww->iarr=zdata[n];
  }

 tensor_free(iarr);

 return(0);		
}

/*****************************************************************************/

