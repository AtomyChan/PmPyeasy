/*****************************************************************************/
/* weight-gen.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Functions related to weight stamp handling.				     */
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

#include "fi.h"

#include "tensor.h"
#include "stars.h"
#include "psf.h"

#include "weight.h"

/*****************************************************************************/

static int weight_compare(const void *v1,const void *v2)
{
 weight	*w1,*w2;
 w1=(weight *)v1,
 w2=(weight *)v2;
 if ( w1->x < w2->x )	return(-1);
 else 			return(1);
}

int weight_sort(weightlist *wl)
{
 if ( wl==NULL || wl->weights==NULL || wl->nweight<=0 )	return(1);

 if ( wl->nweight>1 )
  {	qsort(wl->weights,wl->nweight,sizeof(weight),weight_compare);	}

 return(0);
}

/*****************************************************************************/

weight *weight_get_closest(weightlist *wl,double x0,double y0)
{
 int	best,min,max,mid;
 double	dist,xc,yc,xc0,yc0,xx,yy,xdist,ydist,cdist,mindist,maxdist;
 weight	*ww;
 int	outl,outr,pl,pr;

 if ( wl==NULL || wl->weights==NULL || wl->nweight<=0 )	return(NULL);

 min=0;
 max=wl->nweight;
 while ( max>min )
  {	mid=(min+max)/2;
	dist=x0-wl->weights[mid].x;
	if ( fabs(dist)<1e-10 )	break;
	else if ( dist>0.0 )	min=mid+1;
	else			max=mid;
  }

 xc0=wl->weights[mid].x,
 yc0=wl->weights[mid].y;

 xdist=fabs(x0-xc0),
 ydist=fabs(y0-yc0);
 
 mindist=xdist*xdist+ydist*ydist;
 maxdist=xdist+ydist;

 best=mid;
 for ( outl=outr=1,pl=best-1,pr=best+1 ; outl || outr ; pl--,pr++ )
  {	if ( outl )
	 {	if ( pl>=0 )
		 {	xc=wl->weights[pl].x,
			yc=wl->weights[pl].y;
			xx=x0-xc;
			if ( xx<maxdist )
			 {	yy=fabs(yc-y0);
				if ( yy<=ydist )
			 	 {	ydist=yy;
					dist=xx*xx+yy*yy;
					if ( dist<=mindist )
					 {	best=pl,mindist=dist;
						cdist=fabs(xx)+yy;
						if ( cdist<maxdist )
							maxdist=cdist;
					 }
				 }
			 }
			else	outl=0;
		 }
		else	outl=0;
	 }
	if ( outr )
	 {	if ( pr<wl->nweight )
		 {	xc=wl->weights[pr].x,
			yc=wl->weights[pr].y;
			xx=xc-x0;
			if ( xx<maxdist )
			 {	yy=fabs(yc-y0);
				if ( yy<=ydist )
			 	 {	ydist=yy;
					dist=xx*xx+yy*yy;
					if ( dist<=mindist )
					 {	best=pr,mindist=dist;
						cdist=fabs(xx)+yy;
						if ( cdist<maxdist )
							maxdist=cdist;
					 }
				 }
			 }
			else	outr=0;
		 }
		else	outr=0;
	 }
  }

 ww=&wl->weights[best];
 return(ww);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int weight_get_intersec_list(weightlist *wl,int x0,int y0,int sx,int sy,int *list,int maxcnt)
{
 int	n,fr,fl,xr,xl,i,cnt,yu,yd;
 int	fsize;
 weight	*ww;

 if ( wl==NULL || wl->weights==NULL || wl->nweight<=0 )	return(-1);
 fsize=2*wl->hsize+1;

 xl=x0-1-wl->hsize;
 xr=x0+sx+1+wl->hsize;
 yd=y0-1-wl->hsize;
 yu=y0+sy+1+wl->hsize;

 fl=0,n=wl->nweight;
 while ( n )
  {	if ( xl>wl->weights[fl+n/2].ix )	fl+=n/2+1,n-=n/2+1;
	else					n=n/2;
  };
 fr=0,n=wl->nweight;
 while ( n )
  {	if ( xr>wl->weights[fr+n/2].ix )	fr+=n/2+1,n-=n/2+1;
	else					n=n/2;
  };
 
 cnt=0;
 for ( i=fl ; i<=fr && ( maxcnt<0 || cnt<maxcnt ) ; i++ )
  {	ww=&wl->weights[i];
	if ( yd <= ww->iy && ww->iy <= yu )
	 {	list[cnt]=i;
		cnt++;
	 }
  }

 return(cnt);
}

/*****************************************************************************/
