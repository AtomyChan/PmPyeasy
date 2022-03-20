/*****************************************************************************/
/* convexhull.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* An implementation of the Graham scan algorithm to figure out convex hulls */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2016, 2017; Pal, Andras <apal@szofi.net>				     */
/*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "point.h"
#include "convexhull.h"

static int convexhull_qsort_cb_point_y(const void *v1,const void *v2)
{
 const	point	*p1=v1,*p2=v2;
 if ( p1->y < p2->y )		return(-1);
 else if ( p2->y < p1->y )	return(+1);
 else				return(0);
}

static int convexhull_qsort_cb_point_aux(const void *v1,const void *v2)
{
 const	point	*p1=v1,*p2=v2;
 if ( p1->weight < p2->weight )		return(-1);
 else if ( p2->weight < p1->weight )	return(+1);
 else					return(0);
}

static double convexhull_point_ccw(point *p1,point *p2,point *p3)
{
 return((p2->x-p1->x)*(p3->y-p1->y)-(p2->y-p1->y)*(p3->x-p1->x));
}

/*****************************************************************************/

int convexhull_compute(point *points,int npoint)
{
 point	pzero;
 int	i,m;

 qsort(points,npoint,sizeof(point),convexhull_qsort_cb_point_y);

 for ( i=0 ; i<npoint ; i++ )
  {	double	dx,dy,r;
	dx=points[i].x-points[0].x;
	dy=points[i].y-points[0].y;
	r=sqrt(dx*dx+dy*dy);
	points[i].weight=-dx/r;
  }
 qsort(points+1,npoint-1,sizeof(point),convexhull_qsort_cb_point_aux);
 pzero=points[npoint-1];
 
 m=1;
 for ( i=2 ; i<=npoint ; i++ )
  {	point	wp;
	while ( convexhull_point_ccw(1<m?&points[m-2]:&pzero,&points[m-1],&points[i-1]) <= 0 )
	 {	if ( 1<m )
		 {	m--;
			continue;
		 }
		else if ( i==npoint )
			break;
		else
			i++;
	 }
	m++;
	if ( 0<m )	wp=points[m-1],points[m-1]=points[i-1],points[i-1]=wp;
	else		wp=pzero,pzero=points[i-1],points[i-1]=wp;
  }

 return(m);
}

double convexhull_area(point *hullpoints,int nhullpoint)
{
 int	i;
 double	area;

 if ( nhullpoint<3 )
	return(0.0);

 area=0.0;
 for ( i=1 ; i<nhullpoint-1 ; i++ )
  {	area += hullpoints[i].x*( hullpoints[i+1].y - hullpoints[i-1].y );	}
 area += hullpoints[i].x*(hullpoints[0].y-hullpoints[i-1].y);i++;
 area += hullpoints[0].x*(hullpoints[1].y-hullpoints[i-1].y);

 return(area/2.0);
}


/*****************************************************************************/
