/*****************************************************************************/
/* pselect.c [-> point.h, sort.c] 					     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* A library for selecting points homogenously from a given set.	     */
/* This library is not standalone, it requires the declaration of the 	     */
/* structure 'point' defined in point.h and some functions from sort.c.	     */
/* (c) 2004, Pal, A. (apal@szofi.elte.hu). 				     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in pselect.h	     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "math/point.h"
#include "index/sort.h"
#include "pselect.h"

/*****************************************************************************/

static int point_select_count(point *points,int npoint,double x0,double y0,
			      double x1,double y1,char *ret,int val)
{
 int	n;
 n=0;
 if ( val<0 || ret==NULL )
  {	while ( npoint )
	 {	if ( x0 <= points->x && points->x < x1 &&
		y0 <= points->y && points->y < y1 )
			n++;
		npoint--,points++;
	 };
  }
 else
  {	while ( npoint )
	 {	if ( x0 <= points->x && points->x < x1 &&
		y0 <= points->y && points->y < y1 )
			n++,*ret=val;
		npoint--,points++,ret++;
	 };
  }
 return(n);
}

static int point_select_num(point *points,int npoint,double x0,double y0,
			    double x1,double y1,char *ret,int val,int nn)
{
 int	n;

 if ( val<0 || ret==NULL || nn<0 )	return(0);

 n=0;
 while ( npoint && nn>0 )
  {	if ( x0 <= points->x && points->x < x1 &&
	y0 <= points->y && points->y < y1 )
		n++,nn--,*ret=val;
	npoint--,points++,ret++;
  };
 return(n);
}

static int point_select_compare(int i1,int i2,void *param)
{
 point	*points;

 points=(point *)param;
      if ( points[i1].weight  < points[i2].weight )	return(1);
 else if ( points[i1].weight == points[i2].weight )	return(0);
 else							return(-1);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int point_select(point *points,int npoint,char *ret,int nselect,
		 double x0,double y0,double x1,double y1,int level)
{
 double	xl,xm,xh,
	yl,ym,yh,w;
 int	i,k,n,npp;
 int	n_bl,n_br,n_tl,n_tr,
	s_bl,s_br,s_tl,s_tr;

 if ( nselect==0 )	return(0);

 if ( x0>x1 )	w=x0,x0=x1,x1=w;
 if ( y0>y1 )	w=y0,y0=y1,y1=w;

 npp=0;
 if ( level==0 )
  {	int	*index;
	n=point_select_count(points,npoint,x0,y0,x1,y1,ret,2);
	if ( nselect>=n )
	 {	for ( i=0 ; i<npoint ; i++ )
		 {	if ( ret[i] )	ret[i]=1;		}
		return(n);
	 }
	else
	 {	index=(int *)malloc(sizeof(int)*n);
		for ( i=0,k=0 ; i<npoint ; i++ )
		 {	if ( ret[i]==2 )	index[k]=i,k++;		}
		index_qsort(index,n,point_select_compare,(void *)points);
		for ( i=0 ; i<nselect ; i++ )
		 {	k=index[i];
			ret[k]=1;
		 }
		for ( i=nselect ; i<n ; i++ )
		 {	k=index[i];
			ret[k]=0;
		 }
		free(index);
		return(nselect);
	 }
  }
 else
  {	xl=x0,xm=(x0+x1)*0.5,xh=x1;
	yl=y0,ym=(y0+y1)*0.5,yh=y1;

	if ( nselect==1 )
	 {	npp=point_select_num(points,npoint,xl,yl,xh,yh,ret,1,1);
		return(npp);
	 }

	npp=0;
	n_bl=point_select_count(points,npoint,xl,yl,xm,ym,NULL,0);
	n_br=point_select_count(points,npoint,xm,yl,xh,ym,NULL,0);
	n_tl=point_select_count(points,npoint,xl,ym,xm,yh,NULL,0);
	n_tr=point_select_count(points,npoint,xm,ym,xh,yh,NULL,0);

	n=n_bl+n_br+n_tl+n_tr;
	if ( n==0 )	return(0);

	k=nselect/4;
	s_bl=s_br=s_tl=s_tr=k;
	switch ( nselect-k*4 )
	 {	case 3:	s_br++;
		case 2:	s_tr++;
		case 1:	s_bl++;
	 }

	if ( 0<n_bl && n_bl<=s_bl )
		npp+=point_select_count(points,npoint,xl,yl,xm,ym,ret,1);
	else if ( 0<n_bl )
		npp+=point_select(points,npoint,ret,s_bl,xl,yl,xm,ym,level-1);

	if ( 0<n_br && n_br<=s_br )
		npp+=point_select_count(points,npoint,xm,yl,xh,ym,ret,1);
	else if ( 0<n_br )
		npp+=point_select(points,npoint,ret,s_br,xm,yl,xh,ym,level-1);

	if ( 0<n_tl && n_tl<=s_tl )
		npp+=point_select_count(points,npoint,xl,ym,xm,yh,ret,1);
	else if ( 0<n_tl )
		npp+=point_select(points,npoint,ret,s_tl,xl,ym,xm,yh,level-1);

	if ( 0<n_tr && n_tr<=s_tr )
		npp+=point_select_count(points,npoint,xm,ym,xh,yh,ret,1);
	else if ( 0<n_tr )
		npp+=point_select(points,npoint,ret,s_tr,xm,ym,xh,yh,level-1);

	return(npp);
  }
 return(0);
}
                                          
