/*****************************************************************************/
/* star-cand-trb.c							     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Star candidate searching based on the triangulation algorithm...	     */
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
#include "statistics.h"

#include "fitsh.h"
#include "common.h"
#include "stars.h"

#define		FLAG_MIN	1
#define		FLAG_MAX	2

/*****************************************************************************/

int is_in_triangle(int x1,int y1,int x2,int y2,int x3,int y3,int x,int y)
{
 int    a,b,c;
 a=(x2-x1)*(y-y1)-(x-x1)*(y2-y1);
 b=(x3-x2)*(y-y2)-(x-x2)*(y3-y2);
 c=(x1-x3)*(y-y3)-(x-x3)*(y1-y3);
 if ( (a>=0 && b>=0 && c>=0) || (a<=0 && b<=0 && c<=0) )
	return(1); 
 else
	return(0);
}

static int is_neighbour(tpoint *p1,tpoint *p2,tpoint *points,tpoint ***allneig)
{
 tpoint	**cn; 
 int	i;
 cn=allneig[p2-points];
 for ( i=0 ; cn[i] != NULL ; i++ )
  {	if ( cn[i] == p1 )	return(1);	}
 return(0);
}

static int order_point_neighbours(tpoint *ps,tpoint *p,tpoint **neig,tpoint ***allneig)
{
 tpoint	*cp,*w;
 int	i,j;

 for ( i=0 ; neig[i+1] != NULL && neig[i] != NULL ; i++ )
  {	cp=neig[i];
	for ( j=i+1 ; neig[j] != NULL ; j++ )
	 {	if ( is_neighbour(cp,neig[j],ps,allneig) )	break;	}
	if ( neig[j]==NULL )	return(1);
	else if ( neig[j] != NULL && j>i+1 )	
		w=neig[i+1],neig[i+1]=neig[j],neig[j]=w;
  }
 if ( ! is_neighbour(neig[0],neig[i],ps,allneig) )	return(1); 
 else							return(0);
}

/*****************************************************************************/
/* debug stuff */
/*
 {	FILE	*fw;
	tpoint	*p1,*p2,*p3;
	fw=stdout;
	for ( i=0 ; i<ntri ; i++ )
	 {	p1=tris[i].vertices[0];
		p2=tris[i].vertices[1];
		p3=tris[i].vertices[2];
		fprintf(fw,"%4g %4g %4g %4g\n",p1->xcoord,p1->ycoord,p2->xcoord-p1->xcoord,p2->ycoord-p1->ycoord);
		fprintf(fw,"%4g %4g %4g %4g\n",p2->xcoord,p2->ycoord,p3->xcoord-p2->xcoord,p3->ycoord-p2->ycoord);
		fprintf(fw,"%4g %4g %4g %4g\n",p3->xcoord,p3->ycoord,p1->xcoord-p3->xcoord,p1->ycoord-p3->ycoord);
	 }
  }
*/
/*
 i=0;
 for ( j=0 ; ti.neigs[i][j] != NULL ; j++ )
  {	fprintf(stderr,"%g %g\n",ti.neigs[i][j]->xcoord,ti.neigs[i][j]->ycoord);	}
*/
/*
 for ( i=0 ; i<ncand ; i++ )
  {	double	area;
	area=0.0;
	for ( j=0 ; ti.neigs[i][j] != NULL ; j++ )
	 {	k=(ti.neigs[i][j+1]==NULL?0:j+1);
	 	area+=ti.neigs[i][j]->xcoord*ti.neigs[i][k]->ycoord;
	 	area-=ti.neigs[i][k]->xcoord*ti.neigs[i][j]->ycoord;
	 }
	fprintf(stderr,"%d %d : %g\n",cands[i].ix,cands[i].iy,area);
  }
*/
/*****************************************************************************/

int search_star_candidates_trb(fitsimage *img,char **mask,candidate **rcands,int *rncand,range *srcrange,double treshold)
{
 candidate	*cands,*wc,*nc,*nc1,*nc2;
 int		ncand;
 int		c,i,j,k,l,n,sx,sy,imin,imax,jmin,jmax;
 int		ismin,ismax;
 double		w,w0,area,*darr,med,sig,peak;
 tpoint		*tps,**neig;
 triinfo	ti;
 triangle	*tris;
 int		ntri;

 sx=img->sx,sy=img->sy;
 ncand=0;cands=NULL;

 if ( srcrange==NULL )
  {	imin=0,imax=sy-1,
	jmin=0,jmax=sx-1;
  }
 else
  {	imin=srcrange->ymin,imax=srcrange->ymax;
	jmin=srcrange->xmin,jmax=srcrange->xmax;
	if ( imin>imax )	i=imin,imin=imax,imax=i;
	if ( jmin>jmax )	j=jmin,jmin=jmax,jmax=j;
  }
 if ( imin<0 )		imin=0;
 if ( imax>=sy )	imax=sy-1;
 if ( jmin<0 )		jmin=0;
 if ( jmax>=sx )	jmax=sx-1;

 for ( i=imin ; i<=imax ; i++ )
  {	for ( j=jmin ; j<=jmax ; j++ )
	 {	
		if ( mask != NULL && mask[i][j] )
			continue;

		w0=img->data[i][j];
		ismin=ismax=1;n=0;
		for ( k=-1 ; k<=1 ; k++ )
		 {	if ( i+k<0 || i+k>=sy )	continue;
			for ( l=-1 ; l<=1 ; l++ )
			 {	if ( j+l<0 || j+l>=sx )	continue;
				if ( mask[i+k][j+l] )	continue;
				w=img->data[i+k][j+l];
				if ( w>w0 )		ismax=0;
				else if ( w<w0 )	ismin=0;
				n++;
			 }
		 }
		if ( n==0 || ! ( ismax || ismin ) )	continue;
	
		cands=(candidate *)realloc(cands,sizeof(candidate)*(ncand+1));

		wc=&cands[ncand];

		wc->flags=0;
		if ( ismin )	wc->flags |= FLAG_MIN;
		if ( ismax )	wc->flags |= FLAG_MAX;

		wc->ix=j,wc->cx=(double)j+0.5;
		wc->iy=i,wc->cy=(double)i+0.5;
		wc->peak=0.0;
		wc->sxx=0.0;
		wc->sxy=0.0;
		wc->syy=0.0;

		wc->marked=0;

		wc->ipoints=NULL;
		wc->nipoint=0;

		ncand++;
	 }
  }

 tps=(tpoint *)malloc(sizeof(tpoint)*ncand);
 for ( i=0 ; i<ncand ; i++ )
  {	tps[i].xcoord=cands[i].cx;
	tps[i].ycoord=cands[i].cy;
	tps[i].id=i;
  }

 delaunay_triangulation(tps,ncand,&tris,&ntri,&ti);

 darr=NULL;
 for ( i=0 ; i<ncand ; i++ )
  {	wc=&cands[i];
	neig=ti.neigs[i];
	c=order_point_neighbours(tps,&tps[i],neig,ti.neigs);
	if ( c || (wc->flags&(FLAG_MAX|FLAG_MIN)) != FLAG_MAX )
	 {	wc->marked=1;
		continue;
	 }

	for ( n=0 ; neig[n] != NULL ; )	n++;
	area=0.0;
	darr=(double *)malloc(sizeof(double)*n);
	for ( k=0 ; k<n ; k++ )
	 {	l=(k<n-1?k+1:0);
		area+=neig[k]->xcoord*neig[l]->ycoord;
		area-=neig[l]->xcoord*neig[k]->ycoord;
		nc=cands+(neig[k]-tps);
		darr[k]=img->data[nc->iy][nc->ix];
	 }
	med=median(darr,n);
	for ( k=0 ; k<n ; k++ )
	 {	darr[k]=fabs(darr[k]-med);	}
	sig=median(darr,n);
	peak=img->data[wc->iy][wc->ix]-med;
	if ( peak > 5 * sig && peak > treshold )
	 {	int	xmax,xmin,ymax,ymin,x,y;
		wc->ipoints=NULL;
		wc->nipoint=0;
		xmax=xmin=0;
		ymax=ymin=0;
		for ( k=0 ; k<n ; k++ )
		 {	nc=cands+(neig[k]-tps);
			if ( k==0 ) 
			 {	xmax=xmin=nc->ix;
				ymax=ymin=nc->iy;
			 }
			else
			 {	if ( nc->ix > xmax )	xmax=nc->ix;
				if ( nc->ix < xmin )	xmin=nc->ix;
				if ( nc->iy > ymax )	ymax=nc->iy;
				if ( nc->iy < ymin )	ymin=nc->iy;
			 }
		 }
		wc->ipoints=NULL;
		wc->nipoint=0;
		for ( y=ymin ; y<=ymax ; y++ )
		 {	for ( x=xmin ; x<=xmax ; x++ )
		 	 {	for ( k=0 ; k<n ; k++ )
				 {	l=(k<n-1?k+1:0);
					nc1=cands+(neig[k]-tps);
					nc2=cands+(neig[l]-tps);
					if ( is_in_triangle(wc->ix,wc->iy,nc1->ix,nc1->iy,nc2->ix,nc2->iy,x,y) )	break;
				 }
				if ( k==n )	continue;
				wc->ipoints=(ipoint *)realloc(wc->ipoints,sizeof(ipoint)*(wc->nipoint+1));
				wc->ipoints[wc->nipoint].x=x;
				wc->ipoints[wc->nipoint].y=y;
				wc->nipoint++;
			 }
		 }
		wc->marked=0;
	 }
	else
	 {	wc->ipoints=NULL;
		wc->nipoint=0;
		wc->marked=1;
	 }	

	wc->area=(double)wc->nipoint;
	wc->noise=0.0;

	free(darr);
  }

 free(ti.neigs);
 free(ti.neigpoints);

 cleanup_candlist(&cands,&ncand);

 if ( rcands != NULL )	*rcands=cands;
 if ( rncand != NULL )	*rncand=ncand;

 return(0);
}

