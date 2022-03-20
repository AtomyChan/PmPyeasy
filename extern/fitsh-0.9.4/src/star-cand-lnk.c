/*****************************************************************************/
/* star-cand-lnk.c							     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Star candidate searching based on the `uplink` algorithm...		     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <fits/fits.h>

#include "io/iof.h"
#include "fitsmask.h"
#include "math/fit/lmfit.h"
#include "math/poly.h"
#include "math/point.h"
#include "math/tpoint.h"
#include "math/delaunay.h"
#include "link/linkpoint.h"

#include "statistics.h"
#include "tensor.h"
#include "fitsh.h"
#include "common.h"
#include "stars.h"

/*****************************************************************************/

typedef struct
 {	double	sbg,sbg2,sfx;
	int	nbg,cindex;
 } lbackground;

int search_star_candidates_link(fitsimage *img,char **mask,
	candidate **rcands,int *rncand,
	range *srcrange,double threshold,double fluxthreshold,
	double critical_prominence)
{
 candidate	*cands,*wc;
 int		icand,ncand;
 int		i,j,k,l,m,sx,sy,mx,my;
 double		bg,amp,cx,cy,axx,axy,ayy,peak,flux,noise;

 linkpoint	**lparr;
 linkreference	**lrarr;

 linkrange	lrin,*lr;
 lbackground	*lbgs;

 sx=img->sx,sy=img->sy;

 if ( srcrange==NULL )
 	lr=NULL;
 else
  {	lrin.xmin=srcrange->xmin,lrin.xmax=srcrange->xmax;
	lrin.ymin=srcrange->ymin,lrin.ymax=srcrange->ymax;
	lr=&lrin;
  }

 lparr=linkpoint_create(img->data,sx,sy,lr,mask,1);
 linkpoint_reconnect(lparr,sx,sy);
 lrarr=linkreference_create(lparr,sx,sy);

 if ( critical_prominence >= 0.0 )
  {	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	if ( lrarr[i][j].mcnt && lrarr[i][j].nngh > lrarr[i][j].ncnt )
			 {	int	nx,ny,n,cmask;
				cmask=linkpoint_mask_same_group(lparr,sx,sy,j,i);
				n=linkpoint_local_extreme(img->data,0,0,sx,sy,mask,j,i,+1,cmask,&nx,&ny);
				if ( n<=0 )
					continue;
				else if ( ! linkpoint_is_same_endpoint(lparr,j,i,nx,ny) )
			 	 {	lparr[i][j].nx=nx;
					lparr[i][j].ny=ny;
					/* fprintf(stderr,"connect (%d,%d)->(%d,%d) [%d,%d]:[%d,%d]\n",j,i,nx,ny,lparr[i][j].mx,lparr[i][j].my,lparr[ny][nx].mx,lparr[ny][nx].my); */
				 }
			 }
		 }
	 }
	linkreference_free(lrarr);

	linkpoint_reconnect(lparr,sx,sy);
	lrarr=linkreference_create(lparr,sx,sy);
  }

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	lrarr[i][j].identifier=-1;
		lrarr[i][j].aux=0;
	 }
  }

 icand=0;

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( lparr[i][j].nx<0 || lparr[i][j].ny<0 )
			continue;
		mx=lparr[i][j].mx;
		my=lparr[i][j].my;
		if ( lrarr[my][mx].identifier < 0 )
		 {	lrarr[my][mx].identifier=icand;
			icand++;
		 }
	 }
  }

 if ( icand>0 )
  {	lbgs=(lbackground *)malloc(sizeof(lbackground)*icand);
	for ( i=0 ; i<icand ; i++ )
	 {	lbgs[i].sbg =0.0;
		lbgs[i].sbg2=0.0;
		lbgs[i].sfx =0.0;
		lbgs[i].nbg =0;
		lbgs[i].cindex=-1;
	 }
  }
 else
	lbgs=NULL;

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( lparr[i][j].nx<0 || lparr[i][j].ny<0 )
			continue;
		mx=lparr[i][j].mx;
		my=lparr[i][j].my;
		k=lrarr[my][mx].identifier;
		lbgs[k].sfx += img->data[i][j];
		if ( lrarr[i][j].nsgn<lrarr[i][j].nngh && lrarr[i][j].ncnt<=0 )
		 {	lbgs[k].sbg +=img->data[i][j];
			lbgs[k].sbg2+=img->data[i][j]*img->data[i][j];
			lbgs[k].nbg ++;
		 }
	 }
  }

 cands=NULL;
 ncand=0;
 
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( lrarr[i][j].mcnt <= 0 || lrarr[i][j].identifier<0 )
			continue;
		
		if ( lrarr[i][j].nsgn<8 )
			continue;

		k=lrarr[i][j].identifier;
		
		m=lbgs[k].nbg;
		bg   =lbgs[k].sbg/(double)m;
		noise=lbgs[k].sbg2/(double)m-bg*bg;
		if ( noise>0.0 )	noise=sqrt(noise);
		else			noise=0.0;
		flux=lbgs[k].sfx-lrarr[i][j].mcnt*bg;
		
		k=fit_small_parabola_block_param(img,j,i,&cx,&cy,&axx,&axy,&ayy,&peak);
		if ( k )
			continue;

		/* fprintf(stderr,"bg=%g flux=%g\n",bg,flux); */

		amp=peak-bg;
		if ( amp<0.0 )
			continue;
		else if ( threshold>0.0 && amp<threshold )
			continue;
		else if ( fluxthreshold>0.0 && flux<fluxthreshold )
			continue;

		cands=(candidate *)realloc(cands,sizeof(candidate)*(ncand+1));
		wc=&cands[ncand];

		wc->ix=j,wc->cx=cx+(double)j+0.5;
		wc->iy=i,wc->cy=cy+(double)i+0.5;
		wc->peak=peak;
		wc->amp=amp;
		wc->bg=bg;
		wc->noise=noise;
		wc->sxx=-axx/amp;
		wc->sxy=-axy/amp;
		wc->syy=-ayy/amp;

		wc->marked=0;

		wc->ipoints=(ipoint *)malloc(sizeof(ipoint)*lrarr[i][j].mcnt);
		wc->nipoint=0;
		wc->area=(double)lrarr[i][j].mcnt;
		wc->flux=flux;

		k=lrarr[i][j].identifier;
		lbgs[k].cindex=ncand;

		ncand++;
	 }
  }

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( lparr[i][j].mx<0 || lparr[i][j].my<0 )
			continue;
		mx=lparr[i][j].mx;
		my=lparr[i][j].my;
		k=lrarr[my][mx].identifier;
		l=lbgs[k].cindex;
		if ( l<0 )
			continue;
		wc=&cands[l];
		wc->ipoints[wc->nipoint].x=j;
		wc->ipoints[wc->nipoint].y=i;
		wc->nipoint++;
	 }
  }

 if ( rcands != NULL )	*rcands=cands;
 if ( rncand != NULL )	*rncand=ncand;

 if ( lbgs != NULL )	free(lbgs);

 linkreference_free(lrarr);

 linkpoint_free(lparr);

 return(0);
}

/*****************************************************************************/

