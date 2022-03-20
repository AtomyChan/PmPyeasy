/*****************************************************************************/
/* intersec-cri.c 							     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Standalone library for integrating polynomials surfaces on intersections  */
/* of circles and rectangles.						     */
/* (c) 2006; Pal, A. (apal@szofi.elte.hu).	 			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in intersec-cri.h  */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "intersec-cri.h"

/*****************************************************************************/

static int intersec_cri_int_indefinite(double *pc,int order,double mul,
	double r,double x)
{
 double	t,s,c;
 double	r2,r3,r4;
 double	c2,s2,c4,sc;

 s=x/r;
 if ( s<-1.0 )	s=-1.0;
 if ( s>+1.0 )	s=+1.0;
 t=asin(s);c=cos(t);

 c2=c*c,s2=s*s,c4=c2*c2,sc=s*c;
 r2=r*r,r3=r2*r,r4=r2*r2;

 switch ( order )
  {  case 2:
	pc[3]+=mul*r4*0.125*(t+sc*(s2-c2));		/* x^2		*/
	pc[4]-=mul*r4*0.125*c4;				/* xy		*/
	pc[5]+=mul*r4*(3*t+sc*(2.0*c2+3.0))/24.0;	/* y^2		*/
     case 1:
	pc[1]-=mul*r3*(c2*c)/3.0;			/* x		*/
	pc[2]+=mul*r3*s*(0.5-s2/6.0);			/* y		*/
     case 0:
	pc[0]+=mul*r2*0.5*(t+sc);			/* 1		*/
	break;
  }

 return(0);
}

static int intersec_cri_int_semirect_definite(double *pc,int order,double mul,
	double x1,double x2,double y2)
{
 double	w,sx,sx2,sy,sy2;

 sx=x2+x1,sx2=x2*sx+x1*x1;sx=sx/2.0,sx2=sx2/3.0;
 sy=y2   ,sy2=y2*sy      ;sy=sy/2.0,sy2=sy2/3.0;
 
 w=(x2-x1)*y2*mul;

 switch ( order )
  {  case 2:
	pc[3]+=w*sx2;					/* x^2		*/
	pc[4]+=w*sx*sy;					/* xy		*/
	pc[5]+=w*sy2;					/* y^2		*/
     case 1:
	pc[1]+=w*sx;					/* x		*/
	pc[2]+=w*sy;					/* y		*/
     case 0:
	pc[0]+=w;					/* 1		*/
	break;
  }
 
 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static int intersec_cri_int_semicircle(double *pc,int order,
	double x1,double x2,double y1,double y2,double r)
{
 double	by,sy,dx1,dx2,s,w;

 by=y1;
 sy=y2-y1;

 if ( r<=by || by+sy<=0 )	return(0);
 if ( by<0.0 )			sy+=by,by=0.0;
 if ( sy<=0.0 )			return(0);
 
 if ( by+sy>=r )
  {	if ( by==0.0 )	dx1=-r,dx2=+r,s=0.0;
	else		s=by,w=sqrt(r*r-s*s),dx1=-w,dx2=+w;
	if ( x2<=dx1 || x1>=dx2 )	return(0);
	if ( x1<dx1 )	x1=dx1;
	if ( x2>dx2 )	x2=dx2;
	intersec_cri_int_indefinite(pc,order,+1,r,x2);
	intersec_cri_int_indefinite(pc,order,-1,r,x1);
	intersec_cri_int_semirect_definite(pc,order,-1,x1,x2,by);
	return(0);
  }
 else
  {	double	hs,hx1,hx2,hdx1,hdx2;
	if ( by==0.0 )	dx1=-r,dx2=+r,s=0.0;
	else		s=by,w=sqrt(r*r-s*s),dx1=-w,dx2=+w;
	if ( x2<=dx1 || x1>=dx2 )	return(0);
	if ( x1<dx1 )	x1=dx1;
	if ( x2>dx2 )	x2=dx2;
	hs=by+sy;w=sqrt(r*r-hs*hs);hdx1=-w,hdx2=+w;
	hx1=x1,hx2=x2;
	if ( hx1<hdx1 )	hx1=hdx1;
	if ( hx2>hdx2 )	hx2=hdx2;
	if ( hx1<hx2 )
	 {	if ( x2 != hx2 )
		 {	intersec_cri_int_indefinite(pc,order,+1,r,x2);
			intersec_cri_int_indefinite(pc,order,-1,r,hx2);
		 }
		if ( x1 != hx1 )
		 {	intersec_cri_int_indefinite(pc,order,-1,r,x1);
			intersec_cri_int_indefinite(pc,order,+1,r,hx1);
		 }
		intersec_cri_int_semirect_definite(pc,order,-1,x1,x2,by);
		intersec_cri_int_semirect_definite(pc,order,+1,hx1,hx2,by+sy);
	 }
	else
	 {	intersec_cri_int_indefinite(pc,order,+1,r,x2);
		intersec_cri_int_indefinite(pc,order,-1,r,x1);
		intersec_cri_int_semirect_definite(pc,order,-1,x1,x2,by);
	 }
	return(0);
  }
}

int	intersec_cri_integrate_monoms(double x0,double y0,double sx,double sy,
	double cr,double *coeff,int order)
{
 double	cdis,dx,dy,hdiam;
 double	aclower[6],acupper[6];
 int	nvar,i;

 if ( order<0 )	return(-1);
 if ( cr<0.0 )	cr=-cr;

 if ( order>2 )	order=2;
 nvar=(order+1)*(order+2)/2;

 for ( i=0 ; i<nvar ; i++ )
  {	coeff[i]=0.0;		}

 /* Obvious cases: the circle and the rectangle are disjoint... */
 if ( cr <= x0 || x0+sx <= -cr || cr <= y0 || y0+sy <= -cr || cr==0.0 )
	return(0);

/* Distance of the centers is sqrt(dx*dx+dy*dy)...: */
 dx=x0+0.5*sx;
 dy=y0+0.5*sy;
 cdis =sqrt(dx*dx+dy*dy);
 hdiam=0.5*sqrt(sx*sx+sy*sy);

/* Disjoint: */
 if ( cdis >= hdiam + cr )
	return(0);

/* The rectangle is contained by the circle: */
 if ( cr >= hdiam + cdis )
  {	switch ( order )
	 {   case 2:
		coeff[3]=sy*sx*(x0*x0+x0*sx+sx*sx/3.0);
		coeff[4]=sx*(x0+0.5*sx)*sy*(y0+0.5*sy);
		coeff[5]=sx*sy*(y0*y0+y0*sy+sy*sy/3.0);
	     case 1:
		coeff[1]=sy*sx*(x0+0.5*sx);
		coeff[2]=sx*sy*(y0+0.5*sy);
	     case 0:
		coeff[0]=sx*sy;
		break;
	 }
  }

/* Shift the polynomial coefficients to set the center of the circle to 0:   */
/* acf[0]=pcf[0]+cx*(pcf[1]+cx*pcf[2])+
	cy*((pcf[3]+cx*(pcf[4]+cx*pcf[5]))+
            cy*(pcf[6]+cx*(pcf[7]+cx*pcf[8])));
 acf[1]=pcf[1]+2*cx*pcf[2]+cy*((pcf[4]+2*cx*pcf[5])+cy*(pcf[7]+2*cx*pcf[8]));
 acf[2]=pcf[2]+cy*(pcf[5]+cy*pcf[8]);
 acf[3]=pcf[3]+cx*(pcf[4]+cx*pcf[5])+2*cy*(pcf[6]+cx*(pcf[7]+cx*pcf[8]));
 acf[4]=pcf[4]+2*cx*pcf[5]+2*cy*(pcf[7]+2*cx*pcf[8]);
 acf[5]=pcf[5]+2*cy*pcf[8];
 acf[6]=pcf[6]+cx*(pcf[7]+cx*pcf[8]);
 acf[7]=pcf[7]+2*cx*pcf[8];
 acf[8]=pcf[8]; */

 for ( i=0 ; i<nvar ; i++ )
  {	aclower[i]=0.0;
	acupper[i]=0.0;
  }

/* Calculate the integrals for the upper semicircle: */
 intersec_cri_int_semicircle(acupper,order,x0,x0+sx,+y0,+y0+sy,cr);
/* Calculate the integrals for the lower semicircle: */
 intersec_cri_int_semicircle(aclower,order,x0,x0+sx,-y0-sy,-y0,cr);

 switch ( order )
  {  case 2:
	coeff[3]=acupper[3]+aclower[3];
	coeff[4]=acupper[4]-aclower[4];
	coeff[5]=acupper[5]+aclower[5];
     case 1:
	coeff[1]=acupper[1]+aclower[1];
	coeff[2]=acupper[2]-aclower[2];
     case 0:
	coeff[0]=acupper[0]+aclower[0];
	break;
  }

 return(0);
}

/*****************************************************************************/
                                                         
