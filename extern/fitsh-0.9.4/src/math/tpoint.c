/*****************************************************************************/
/* tpoint.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Operations on 2D points and arrays.					     */
/*****************************************************************************/

/* tpoint.c : operation on 2D points and arrays (domsa@konkoly.hu (2002)) */
/* tpoint.c,v 5.5 2003/05/13 16:09:43 domsa Exp */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tpoint.h"

/*****************************************************************************/

/* create/set/destroy a point array */
tpointarr  *tpoint_createarr(void)
{
 tpointarr *arr=malloc(sizeof(*arr));
 if ( arr != NULL )
  {	arr->length = 0;
	arr->points = NULL;
  }
 return(arr);
}

tpointarr  *tpoint_buildarr(int length,tpoint *points)
{
 tpointarr *arr=malloc(sizeof(*arr));

 if ( arr != NULL )
  {	arr->length = length;
	arr->points = points;
  }
 return(arr);
}

void  tpoint_destroyarr(tpointarr *pa)
{
 if ( pa != NULL )
  {	if ( pa->points != NULL )	free(pa->points);
	free(pa);
  }
 return;
}

/* calculate eucledian distance between two 2D points */
double  tpoint_eucdist(tpoint *p1,tpoint *p2)
{
 double	xx,yy;

 xx=p1->xcoord-p2->xcoord;
 yy=p1->ycoord-p2->ycoord;

 return(sqrt(xx*xx+yy*yy));
}

/* calculate squared eucledian distance between two 2D points */
double  tpoint_eucdist2(tpoint *p1, tpoint *p2)
{
 double  xx, yy;

 xx=p1->xcoord-p2->xcoord;
 yy=p1->ycoord-p2->ycoord;

 return(xx*xx+yy*yy);
}

/* return minimum in x coordinates or 0 */
double tpoint_minx(int nelem,tpoint *points)
{
 int	ii;
 double	 minx;

 if ( nelem < 1 || points == NULL )	return(0);
 minx=points[0].xcoord;
 for ( ii=1 ; ii<nelem ; ii++ )
  {	if ( points[ii].xcoord<minx )	minx=points[ii].xcoord;		}

 return(minx);
}

/* return minimum in y coordinates or 0 */
double tpoint_miny(int nelem,tpoint *points)
{
 int	ii;
 double	miny;

 if ( nelem < 1 || points == NULL )	return(0);
 miny=points[0].ycoord;
 for ( ii=1 ; ii<nelem ; ii++ )
  {	if ( points[ii].ycoord<miny )	miny=points[ii].ycoord;		}

 return(miny);
}

/* return maximum in x coordinates or 0 */
double tpoint_maxx(int nelem,tpoint *points)
{
 int	ii;
 double	 maxx;

 if ( nelem < 1 || points == NULL )	return(0);
 maxx=points[0].xcoord;
 for ( ii=1 ; ii<nelem ; ii++ )
  {	if ( points[ii].xcoord>maxx )	maxx=points[ii].xcoord;		}

 return(maxx);
}

/* return maximum in y coordinates or 0 */
double tpoint_maxy(int nelem,tpoint *points)
{
 int	ii;
 double	maxy;

 if ( nelem < 1 || points == NULL )	return(0);
 maxy=points[0].ycoord;
 for ( ii=1 ; ii<nelem ; ii++ )
  {	if ( points[ii].ycoord>maxy )	maxy=points[ii].ycoord;		}

 return(maxy);
}

/* return the minimum and maximum coordinates along the x-axis */
void tpoint_minmaxx(int nelem,tpoint *points,double *minx,double *maxx)
{
 int  ii;

 *minx=*maxx=points->xcoord;
 for ( ii = 0 ; ii<nelem ; ii++ )
  {	if ( *minx>points[ii].xcoord )		*minx=points[ii].xcoord;
	else if ( *maxx<points[ii].xcoord )	*maxx=points[ii].xcoord;
  }
}

/* return the minimum and maximum coordinates along the y-axis */
void tpoint_minmaxy(int nelem,tpoint *points,double *miny,double *maxy)
{
 int  ii;

 *miny=*maxy=points->ycoord;
 for ( ii = 0 ; ii<nelem ; ii++ )
  {	if ( *miny>points[ii].ycoord )		*miny=points[ii].ycoord;
	else if ( *maxy<points[ii].ycoord )	*maxy=points[ii].ycoord;
  }
}

/* get length of x side */
double tpoint_lengthx(int nelem,tpoint *points)
{
 double	minx,maxx;
 tpoint_minmaxx(nelem,points,&minx,&maxx);
 return(maxx-minx);
}

/* get length of y side */
double tpoint_lengthy(int nelem,tpoint *points)
{
 double miny,maxy;
 tpoint_minmaxy(nelem,points,&miny,&maxy);
 return(maxy-miny);
}

/* flip along x axis */
void tpoint_flipx(int nelem,tpoint *points)
{
 int	ii;
 double	minx,maxx,midx;

 tpoint_minmaxx(nelem,points,&minx,&maxx);
 midx=(minx+maxx)/2.0;
 for ( ii=0 ; ii<nelem ; ii++ )
  {	points[ii].xcoord=2.0*midx-points[ii].xcoord;	}
}

/* flip along y axis */
void tpoint_flipy(int nelem,tpoint *points)
{
 int	ii;
 double	miny,maxy,midy;

 tpoint_minmaxy(nelem,points,&miny,&maxy);
 midy=(miny+maxy)/2.0;
 for ( ii=0 ; ii<nelem ; ii++ )
  {	points[ii].ycoord=2.0*midy-points[ii].ycoord;	}
}

/* zoom in x axis */
void tpoint_magx(int nelem,tpoint *points,double xmag)
{
 int	ii;
 for ( ii=0 ; ii<nelem ; ii++ )
  {	points[ii].xcoord *= xmag;	}
}

/* zoom in y axis */
void tpoint_magy(int nelem, tpoint *points, double ymag)
{
 int	ii;
 for ( ii=0 ; ii<nelem ; ii++ )
  {	points[ii].ycoord *= ymag;	}
}

/* zoom out points */
void tpoint_zoom(int nelem, tpoint *points, double zoom)
{
 int	ii;
 for ( ii=0 ; ii<nelem ; ii++ )
  {	points[ii].xcoord *= zoom;
	points[ii].ycoord *= zoom;
  }
}

/* shift along x axis */
void tpoint_shiftx(int nelem, tpoint *points, double xshift)
{
 int	ii;
 for ( ii=0 ; ii<nelem ; ii++ )
  {	points[ii].xcoord += xshift;	}
}

/* shift along y axis */
void tpoint_shifty(int nelem, tpoint *points, double yshift)
{
 int	ii;
 for ( ii=0 ; ii<nelem ; ii++ )
  {	points[ii].ycoord += yshift;	}
}

/* rotate points clockwise around the origo (rotation angle in radians) */
void tpoint_rotate(int nelem, tpoint *points,double rotang)
{
 int	ii;
 double	xcoord, ycoord;
 double	cphi=cos(rotang);
 double	sphi=sin(rotang);

 for ( ii=0 ; ii<nelem ; ii++ )
  {	xcoord=points[ii].xcoord;
	ycoord=points[ii].ycoord;
	points[ii].xcoord =  cphi*xcoord + sphi*ycoord;
	points[ii].ycoord = -sphi*xcoord + cphi*ycoord;
  }
}

/* rotate points clockwise around the origo (rotation angle in degrees) */
void tpoint_rrotate(int nelem, tpoint *points, double rotang)
{
 tpoint_rotate(nelem, points, rotang*M_PI/180.0);
}

/* sort along x-axis (for qsort()) */
int tpoint_sortx(const void *p1,const void *p2)
{
 double	x1=((tpoint *)p1)->xcoord;
 double	x2=((tpoint *)p2)->xcoord;
 return(x2<x1?1:-1);
}

/* sort along y-axis (for qsort()) */
int tpoint_sorty(const void *p1,const void *p2)
{
 double y1=((tpoint *)p1)->ycoord;
 double y2=((tpoint *)p2)->ycoord;
 return(y2<y1?1:-1);
}

/* sort points based on properties (for qsort()) */
/*
int tpoint_sortprop(const void *p1,const void *p2)
{
 double prop1=((tpoint *)p1)->prop;
 double prop2=((tpoint *)p2)->prop;
 return(prop2<prop1?1:-1);
}
*/

/* sort points based on ids (for qsort()) */
int tpoint_sortid(const void *p1,const void *p2)
{
 int	prop1=((tpoint *)p1)->id;
 int	prop2=((tpoint *)p2)->id;
 return(prop2<prop1?1:-1);
}

/* copy a point array into another one */
tpoint *tpoint_copy(int nelem,tpoint *points)
{
 tpoint  *dest;
 size_t  size=nelem*sizeof(*dest);

 if ( nelem < 1 || points == NULL )
	return(NULL);
 if ( (dest=malloc(size)) != NULL )
	memcpy((void *)dest,(void *)points,size);
 return(dest);
}

/* remove a point from an array */
void tpoint_remove(int idx, int nelem,tpoint *points)
{
 if ( nelem > 1 && points != NULL && idx > 0 && idx < --nelem )
  {	points += idx;
	memmove((void *)points,(void *)(points+1),(nelem-idx)*sizeof(*points));
  }
}

/*****************************************************************************/
    
