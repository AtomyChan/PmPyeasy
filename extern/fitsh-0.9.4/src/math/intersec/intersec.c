/*****************************************************************************/
/* intersec.c 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Standalone library for the determination of intersections of recangles    */
/* and/or circles.							     */
/* (c) 2004, Pal, A. (apal@szofi.elte.hu).	 			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in intersec.h      */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "intersec.h"

int intersec_rectangles(drectangle *r1,drectangle *r2,drectangle *ir)
{
 double	p1,l1,p2,l2,p,l;
 ir->x=ir->y=ir->sx=ir->sy=0.0;
 p1=r1->x,l1=r1->sx,
 p2=r2->x,l2=r2->sx;
 if ( fabs((2*p1+l1)-(2*p2+l2)) >= l1+l2 )	return(1);
 if ( p2>p1 )	{ p=p2;l=p1+l1-p2;if ( l2<l ) l=l2; }
 else		{ p=p1;l=p2+l2-p1;if ( l1<l ) l=l1; }
 ir->x=p,ir->sx=l;
 p1=r1->y,l1=r1->sy,
 p2=r2->y,l2=r2->sy;
 if ( fabs((2*p1+l1)-(2*p2+l2)) >= l1+l2 )	return(1);
 if ( p2>p1 )	{ p=p2;l=p1+l1-p2;if ( l2<l ) l=l2; }
 else		{ p=p1;l=p2+l2-p1;if ( l1<l ) l=l1; }
 ir->y=p,ir->sy=l;
 return(0);
}
double area_rectangle(drectangle *r)
{
 double	a;
 a=r->sx*r->sy;
 return(a);
}
double area_circle(dcircle *c)
{
 double a;
 a=M_PI*(c->radius)*(c->radius);
 return(a);
}

/* The indefinite integral of \sqrt{r^2-x^2}... */
static double area_indefinite(double r,double x)
{
 double	t,w;
 if ( x<=-r )		return(0.5*r*r*asin(-1.0));
 else if ( x>=r )	return(0.5*r*r*asin(+1.0));
 t=x/r;
 w=0.5*r*r*(asin(t)+t*sqrt(1.0-t*t));
 return(w);
}
static double area_intersec_of_rect_semic(double cx,double cy,double r,double by,double sy,double x1,double x2)
{
 double dx1,dx2,s,w,a;

 if ( x1==x2 )		return(0.0);
 else if ( x1>x2 )	w=x1,x1=x2,x2=w;
 if ( sy==0.0 )		return(0.0);
 else if ( sy<0 )	by+=sy,sy=-sy;

 if ( cy+r<=by || by+sy<=cy )	return(0.0);
 if ( by<cy )			sy-=cy-by,by=cy;
 if ( sy<=0.0 )			return(0.0);

 if ( by+sy>=cy+r )
  {	if ( by==cy )	dx1=cx-r,dx2=cx+r,s=0.0;
	else		s=by-cy,w=sqrt(r*r-s*s),dx1=cx-w,dx2=cx+w;
	if ( x2<=dx1 || x1>=dx2 )	return(0.0);
	if ( x1<dx1 )	x1=dx1;
	if ( x2>dx2 )	x2=dx2;
	a=area_indefinite(r,x2-cx)-area_indefinite(r,x1-cx)-s*(x2-x1);
	return(a);
  }
 else
  {	double	hs,hx1,hx2,hdx1,hdx2,a1,a2;
	if ( by==cy )	dx1=cx-r,dx2=cx+r,s=0.0;
	else		s=by-cy,w=sqrt(r*r-s*s),dx1=cx-w,dx2=cx+w;
	if ( x2<=dx1 || x1>=dx2 )	return(0.0);
	if ( x1<dx1 )	x1=dx1;
	if ( x2>dx2 )	x2=dx2;
	hs=by+sy-cy;w=sqrt(r*r-hs*hs);hdx1=cx-w,hdx2=cx+w;
	hx1=x1,hx2=x2;
	if ( hx1<hdx1 )	hx1=hdx1;
	if ( hx2>hdx2 )	hx2=hdx2;
	if ( hx1<hx2 )
	 {	if ( x2==hx2 )	a2=0.0;
		else		a2=area_indefinite(r,x2-cx)-area_indefinite(r,hx2-cx);
		if ( x1==hx1 )	a1=0.0;
		else		a1=area_indefinite(r,x1-cx)-area_indefinite(r,hx1-cx);
		a=(a2-a1) - (s*(x2-x1)-hs*(hx2-hx1));
	 }
	else
		a=area_indefinite(r,x2-cx)-area_indefinite(r,x1-cx)-s*(x2-x1);
	return(a);
  }
}

double area_intersec_of_rect_circ(drectangle *r,dcircle *c)
{
 double aupper,alower;
 aupper=area_intersec_of_rect_semic(c->x,+c->y,c->radius,+r->y,r->sy,r->x,r->x+r->sx);
 alower=area_intersec_of_rect_semic(c->x,-c->y,c->radius,-r->y-r->sy,r->sy,r->x,r->x+r->sx);
 return(aupper+alower);
}

/*****************************************************************************/
