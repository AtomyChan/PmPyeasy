/*****************************************************************************/
/* star-cand-pp.c							     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Star candidate searching based on the (old) `parabola peak' algorithm.    */
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
#include "link/floodfill.h"

#include "fitsh.h"
#include "common.h"
#include "stars.h"

/*****************************************************************************/

int search_star_candidates(fitsimage *img,char **mask,candidate **rcands,int *rncand,range *srcrange,double threshold,spatial *spbg,double skysigma)
{
 candidate	*cands,*wc;
 int		ncand;

 int		i,j,k,l,sx,sy,imin,imax,jmin,jmax;
 double		w,is,is2,ns,ia,ia2,na;

 double		pfit[6],a,ax,ay,axx,axy,ayy,det,tr,cx,cy,peak,amp,bg;

 sx=img->sx,sy=img->sy;
 ncand=0;cands=NULL;

 if ( srcrange==NULL )
  {	imin=2,imax=sy-3,
	jmin=2,jmax=sx-3;
  }
 else
  {	imin=srcrange->ymin,imax=srcrange->ymax;
	jmin=srcrange->xmin,jmax=srcrange->xmax;
	if ( imin>imax )	i=imin,imin=imax,imax=i;
	if ( jmin>jmax )	j=jmin,jmin=jmax,jmax=j;
  }
 if ( imin<2 )		imin=2;
 if ( imax>=sy-2 )	imax=sy-3;
 if ( jmin<2 )		jmin=2;
 if ( jmax>=sx-2 )	jmax=sx-3;

 for ( i=imin ; i<=imax ; i++ )
  {	for ( j=jmin ; j<=jmax ; j++ )
	 {	double	sky;
		if ( mask != NULL )
		 {	if ( mask[i][j] || mask[i][j+1] || mask[i][j-1] ||
			     mask[i+1][j] || mask[i-1][j] )
					continue;
		 }

		w=img->data[i][j];
		if ( ! ( img->data[i][j-1]<w && img->data[i][j+1]<w ) )	continue;
		if ( ! ( img->data[i-1][j]<w && img->data[i+1][j]<w ) )	continue;
		k=0;

		w=img->data[i-1][j];l=0;
		if ( img->data[i-1][j-1]<w )	l++; 
		if ( img->data[i-1][j+1]<w )	l++;
		k+=l*2;

		w=img->data[i+1][j];l=0;
		if ( img->data[i+1][j-1]<w )	l++; 
		if ( img->data[i+1][j+1]<w )	l++;
		k+=l*2;

		w=img->data[i][j-1];l=0;
		if ( img->data[i-1][j-1]<w )	l++; 
		if ( img->data[i+1][j-1]<w )	l++;
		k+=l*2;

		w=img->data[i][j+1];l=0;
		if ( img->data[i-1][j+1]<w )	l++; 
		if ( img->data[i+1][j+1]<w )	l++;
		k+=l*2;

		if ( k<12 )	continue;
		
		is=is2=0.0;ns=0;
		for ( k=-1 ; k<=1 ; k++ )
		 {  for ( l=-1 ; l<=1 ; l++ )
		     {	w=img->data[i+k][j+l];
			is+=w,is2+=w*w;ns++;
		     }
		 }
		is/=(double)ns;is2/=(double)ns;
		
		ia=ia2=0.0;na=0;
		for ( k=-1 ; k<=1 ; k++ )
		 {	w=img->data[i-2][j+k];ia+=w,ia2+=w*w;
			w=img->data[i+2][j+k];ia+=w,ia2+=w*w;
			w=img->data[i+k][j-2];ia+=w,ia2+=w*w;	
			w=img->data[i+k][j+2];ia+=w,ia2+=w*w; na+=4;
		 }
		ia/=(double)na;ia2/=(double)na;
		
		if ( ! ( is>ia+skysigma ) )	continue;
		
		fit_small_parabola_block(img,j,i,pfit);

		a=pfit[0],ax=pfit[1],ay=pfit[2],
		axx=pfit[3],axy=pfit[4],ayy=pfit[5];
		det=axx*ayy-axy*axy; tr=axx+ayy;
		if ( det<=0.0 || tr>=0.0 )	continue;
		cx=-(+ayy*ax-axy*ay)/det,
		cy=-(-axy*ax+axx*ay)/det;
		if ( fabs(cx)>1 || fabs(cy)>1 )	continue;

		peak=a+ax*cx+ay*cy+0.5*(axx*cx*cx+2.0*axy*cx*cy+ayy*cy*cy);

		if ( spbg==NULL )
			sky=0.0;
		else
			sky=eval_2d_poly((double)j,(double)i,spbg->order,spbg->coeff,spbg->ox,spbg->oy,spbg->scale);

		if ( peak-sky < threshold && sky>0.0 && threshold>0.0 )	continue;

		bg=sky;
		amp=peak-sky;
	
		cands=(candidate *)realloc(cands,sizeof(candidate)*(ncand+1));
	
		wc=&cands[ncand];

		wc->ix=j,wc->cx=cx+(double)j+0.5;
		wc->iy=i,wc->cy=cy+(double)i+0.5;
		wc->peak=peak;
		wc->amp=amp;
		wc->bg=bg;
		wc->sxx=-axx/amp;
		wc->sxy=-axy/amp;
		wc->syy=-ayy/amp;

		wc->marked=0;

		wc->ipoints=NULL;
		wc->nipoint=0;
	
		wc->area=0.0;
		wc->noise=0.0;

		ncand++;
	 }
  }

 if ( rcands != NULL )	*rcands=cands;
 if ( rncand != NULL )	*rncand=ncand;

 return(0);
}

/*****************************************************************************/

typedef struct 
 {	fitsimage	*img;
	int		bx,by;
	int		bsx,bsy;
	int		*fillmask;
	double		cx,cy;
 } fillparam;

int getpixel_fill(void *param,int x,int y)
{
 fillparam	*fp;
 double		iix,iiy,dx,dy,w;

 fp=(fillparam *)param;
 if ( x<0 || y<0 || x>=fp->bsx       || y>=fp->bsy         )	return(1);
 if ( fp->fillmask[y*fp->bsx+x] ) 				return(1);
 if ( x==fp->bsx/2 && y==fp->bsy/2 )				return(0);
 x+=fp->bx,y+=fp->by;
 if ( x<1 || y<1 || x>=fp->img->sx-1 || y>=fp->img->sy-1   )	return(1);

 iix=(fp->img->data[y-1][x+1]-fp->img->data[y-1][x-1])+
     (fp->img->data[y  ][x+1]-fp->img->data[y  ][x-1])+
     (fp->img->data[y+1][x+1]-fp->img->data[y+1][x-1]),

 iiy=(fp->img->data[y+1][x-1]-fp->img->data[y-1][x-1])+
     (fp->img->data[y+1][x  ]-fp->img->data[y-1][x  ])+
     (fp->img->data[y+1][x+1]-fp->img->data[y-1][x+1]);

 dx=-((double)x+0.5-fp->cx);
 dy=-((double)y+0.5-fp->cy);
 w=iix*dx+iiy*dy;
 if ( w<0.0 )						return(1);
 if ( w*w <= (dx*dx+dy*dy)*(iix*iix+iiy*iiy)*0.9 )	return(1);
 return(0);
}
void setpixel_fill(void *param,int x,int y)
{
 fillparam *fp;
 fp=(fillparam *)param;
 if ( x<0 || y<0 || x>=fp->bsx    || y>=fp->bsy      )	return;
 fp->fillmask[y*fp->bsx+x]=1;
}

int markout_candidates(fitsimage *img,char **mask,candidate *cands,int ncand)
{
 int		i,j,k,sx,sy,n,bx,by;
 double		avgd,bg;
 int		*fillmask;
 int		hbs,bs;
 fillparam	fp;
 candidate	*wc;

 if ( ncand == 0 )	return(0);
 if ( ncand <  0 )	return(1);

 sx=img->sx,sy=img->sy;
 avgd=sqrt((double)sx*(double)sy/(double)ncand);

 hbs=(int)(avgd*2.0);
 hbs=10;

 bs=2*hbs;
 fillmask=(int *)malloc(sizeof(int)*bs*bs);

 for ( n=0 ; n<ncand ; n++ )
  {	
	wc=&cands[n];
	wc->marked=0;

	fp.img=img;
	fp.bsx=fp.bsy=bs;
	fp.bx=bx=wc->ix-hbs;
	fp.by=by=wc->iy-hbs;
	fp.cx=wc->cx,
	fp.cy=wc->cy;
	for ( i=0 ; i<bs*bs ; i++ )
	 {	fillmask[i]=0;			}
	fp.fillmask=fillmask;
	floodfill(hbs,hbs,getpixel_fill,setpixel_fill,(void *)&fp);

	for ( i=0,k=0 ; i<bs ; i++ )
	 {  for ( j=0 ; j<bs ; j++ )	
	     {	int	li,lj;
		li=i+by,lj=j+bx;
		if ( li<1 || lj<1 || li>=sy-1 || lj>=sx-1 )	continue;
		if ( mask != NULL && mask[li][lj] )		continue;
	    	if ( fillmask[i*bs+j] )	k++;
	     }
	 }
	if ( k<9 )
	 {	wc->marked=1;
		continue; 
	 }

	wc->nipoint=k;
	wc->ipoints=(ipoint *)malloc(sizeof(ipoint)*k);
	bg=0;
	for ( i=0,k=0 ; i<bs  ; i++ )
	 {  for ( j=0 ; j<bs ; j++ )
	     {	int	li,lj;
		li=i+by,lj=j+bx;
		if ( li<1 || lj<1 || li>=sy-1 || lj>=sx-1 )	continue;
		if ( mask != NULL && mask[li][lj] )		continue;
		if ( k==0 || img->data[li][lj]<bg )	bg=img->data[li][lj];
	    	if ( fillmask[i*bs+j] )	
			wc->ipoints[k].x=lj,
			wc->ipoints[k].y=li,
			k++;
	     }
	 }
	wc->bg=bg;
	wc->amp=wc->peak-bg;
	wc->area=(double)wc->nipoint;
	wc->noise=0.0;
  }

 free(fillmask);

 return(0);
}

