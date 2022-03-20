/*****************************************************************************/
/* biquad-isc.c 							     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Standalone library for integrating biquadratic surfaces on intersections: */
/* on circles, triangles, polygons and rectangles. Primary developed for     */
/* weighted aperture photometry. This library is standalone, however, some   */
/* related functions (to create and evaluate biquadratic interpolation       */
/* surfaces) can be found in biquad.c.					     */
/* (c) 2004-05, Pal, A. (apal@szofi.elte.hu).	 			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in biquad-isc.h    */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "biquad-isc.h"

/*****************************************************************************/

/* Nontrivial coefficients of the biquadratic polynomial; derived from the   */
/* appropriate subset of the array 'c' holding the surface coefficients      */
/* (which are calculated by biquad_coeff(), for example).		     */
static int biquad_poly_coefficients(double **cf,double *pcf)
{
 double	wa1,wa2,wb1,wb2;

 pcf[0]=cf[0][0];					/* 1		*/
 pcf[1]=-4.0*cf[0][0]+6.0*cf[0][1]-2.0*cf[0][2];	/* x		*/
 wa1   =-4.0*cf[1][0]+6.0*cf[1][1]-2.0*cf[1][2];	
 wa2   =-4.0*cf[2][0]+6.0*cf[2][1]-2.0*cf[2][2];
 pcf[2]=3.0*(cf[0][0]-2.0*cf[0][1]+cf[0][2]);		/* x^2		*/
 wb1   =3.0*(cf[1][0]-2.0*cf[1][1]+cf[1][2]);
 wb2   =3.0*(cf[2][0]-2.0*cf[2][1]+cf[2][2]);
 pcf[3]=-4.0*cf[0][0]+6.0*cf[1][0]-2.0*cf[2][0];	/* y		*/
 pcf[4]=-4.0*pcf[1]+6.0*wa1-2.0*wa2;			/* xy		*/
 pcf[5]=-4.0*pcf[2]+6.0*wb1-2.0*wb2;			/* x^2y		*/
 pcf[6]=3.0*(cf[0][0]-2.0*cf[1][0]+cf[2][0]);		/* y^2		*/
 pcf[7]=3.0*(pcf[1]-2.0*wa1+wa2);			/* xy^2		*/
 pcf[8]=3.0*(pcf[2]-2.0*wb1+wb2);			/* x^2y^2	*/

 return(0);
}

static double biquad_isc_integral_indefinite(double *pc,double r,double x)
{
 double	t,w,s,c,
	r2,r3,r4,r5,r6,
	c2,s2,c4,s4,sc,s2c2;

 s=x/r;
 if ( s<-1.0 )	s=-1.0;
 if ( s>+1.0 )	s=+1.0;
 t=asin(s);c=cos(t);

 c2=c*c,s2=s*s,c4=c2*c2,s4=s2*s2,sc=s*c,s2c2=s2*c2;
 r2=r*r,r3=r2*r,r4=r2*r2,r5=r2*r3,r6=r3*r3;

 w=+pc[0]*r2*0.5*(t+sc)					/* 1		*/
   -pc[1]*r3*(c2*c)/3.0					/* x		*/
   +pc[2]*r4*0.125*(t+sc*(s2-c2))			/* x^2		*/
   +pc[3]*r3*s*(0.5-s2/6.0)				/* y		*/
   -pc[4]*r4*c4/8.0					/* xy		*/
   +pc[5]*r5*s*(s4+2.5*s2c2)/15.0			/* x^2y		*/
   +pc[6]*r4*(3*t+sc*(2.0*c2+3.0))/24.0			/* y^2		*/
   -pc[7]*r5*c4*c/15.0					/* xy^2		*/
   +pc[8]*r6*((t+sc*(s4-c4))/48.0+s2c2*sc/18.0);	/* x^2y^2	*/

 return(w);
}

static double biquad_isc_integral_semirect_definite(double *pc,double x1,double x2,double y2)
{
 double	w,sx,sx2,sy,sy2;

 sx=x2+x1,sx2=x2*sx+x1*x1;sx=sx/2.0,sx2=sx2/3.0;
 sy=y2   ,sy2=y2*sy      ;sy=sy/2.0,sy2=sy2/3.0;
 
 w=(x2-x1)*y2* (    (pc[0]+pc[1]*sx+pc[2]*sx2) +
		sy *(pc[3]+pc[4]*sx+pc[5]*sx2) +
		sy2*(pc[6]+pc[7]*sx+pc[8]*sx2) );

 return(w);
}

static double biquad_isc_integral_rect_definite(double *pc,double x1,double x2,double y1,double y2)
{
 double	w,sx,sx2,sy,sy2;

 sx=x2+x1,sx2=x2*sx+x1*x1;sx=sx/2.0,sx2=sx2/3.0;
 sy=y2+y1,sy2=y2*sy+y1*y1;sy=sy/2.0,sy2=sy2/3.0;
 
 w=(x2-x1)*(y2-y1)*(        (pc[0]+pc[1]*sx+pc[2]*sx2) +
			sy *(pc[3]+pc[4]*sx+pc[5]*sx2) +
			sy2*(pc[6]+pc[7]*sx+pc[8]*sx2) );

 return(w);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static double biquad_isc_int_pixel_semicircle(double *pc,double x1,double x2,
					      double y1,double y2,double r)
{
 double	by,sy,dx1,dx2,s,w,a;

 by=y1;
 sy=y2-y1;

 if ( r<=by || by+sy<=0 )	return(0.0);
 if ( by<0.0 )			sy+=by,by=0.0;
 if ( sy<=0.0 )			return(0.0);
 
 if ( by+sy>=r )
  {	if ( by==0.0 )	dx1=-r,dx2=+r,s=0.0;
	else		s=by,w=sqrt(r*r-s*s),dx1=-w,dx2=+w;
	if ( x2<=dx1 || x1>=dx2 )	return(0.0);
	if ( x1<dx1 )	x1=dx1;
	if ( x2>dx2 )	x2=dx2;
	a=biquad_isc_integral_indefinite(pc,r,x2)-
	  biquad_isc_integral_indefinite(pc,r,x1)-
	  biquad_isc_integral_semirect_definite(pc,x1,x2,by);
	return(a);
  }
 else
  {	double	hs,hx1,hx2,hdx1,hdx2,a1,a2;
	if ( by==0.0 )	dx1=-r,dx2=+r,s=0.0;
	else		s=by,w=sqrt(r*r-s*s),dx1=-w,dx2=+w;
	if ( x2<=dx1 || x1>=dx2 )	return(0.0);
	if ( x1<dx1 )	x1=dx1;
	if ( x2>dx2 )	x2=dx2;
	hs=by+sy;w=sqrt(r*r-hs*hs);hdx1=-w,hdx2=+w;
	hx1=x1,hx2=x2;
	if ( hx1<hdx1 )	hx1=hdx1;
	if ( hx2>hdx2 )	hx2=hdx2;
	if ( hx1<hx2 )
	 {	if ( x2==hx2 )	a2=0.0;
		else		a2=biquad_isc_integral_indefinite(pc,r,x2)-
				   biquad_isc_integral_indefinite(pc,r,hx2);
		if ( x1==hx1 )	a1=0.0;
		else		a1=biquad_isc_integral_indefinite(pc,r,x1)-
				   biquad_isc_integral_indefinite(pc,r,hx1);
		a=(a2-a1)-
		  biquad_isc_integral_semirect_definite(pc,x1,x2,by)+
		  biquad_isc_integral_semirect_definite(pc,hx1,hx2,by+sy);
	 }
	else
	 {	a=biquad_isc_integral_indefinite(pc,r,x2)-
		  biquad_isc_integral_indefinite(pc,r,x1)-
		  biquad_isc_integral_semirect_definite(pc,x1,x2,by);
	 }
	return(a);
  }
}

double biquad_isc_int_pixel_circle(double **c,int ix,int iy,
				   double cx,double cy,double cr)
{
 double	cdis2,dx,dy;
 double	iupper,ilower,*cf[3],pcf[9],acf[9];

 if ( cx+cr <= (double)ix )	return(0.0);	/* Obvious cases: the circle */
 if ( (double)(ix+1) <= cx-cr )	return(0.0);	/* and the [ix:ix+1,iy:iy+1] */
 if ( cy+cr <= (double)iy )	return(0.0);	/* pixel are disjoint...     */
 if ( (double)(iy+1) <= cy-cr )	return(0.0);

/* Distance of the centers is sqrt(dx*dx+dy*dy)...: */
 dx=(double)ix+0.5-cx;
 dy=(double)iy+0.5-cy;
 cdis2=dx*dx+dy*dy;

/* Disjoint: */
 if ( cdis2 >= cr*(cr+1.42)+0.5 )
	return(0.0);	
/* The pixel is contained by the circle: */
 if ( cdis2 <= cr*(cr-1.42)+0.5 && cr>=0.71 )
	return(c[2*iy+1][2*ix+1]);

/* The subset of the coefficients 'c' required by the integration...:	     */
 cf[0]=&c[2*iy+0][2*ix];
 cf[1]=&c[2*iy+1][2*ix];
 cf[2]=&c[2*iy+2][2*ix];
 biquad_poly_coefficients(cf,pcf);

/* Shift the pixel and the circle to the bottom-left corner: */
 cx-=(double)ix,
 cy-=(double)iy;

/* Shift the polynomial coefficients to set the center of the circle to 0:   */
 acf[0]=pcf[0]+cx*(pcf[1]+cx*pcf[2])+
	cy*((pcf[3]+cx*(pcf[4]+cx*pcf[5]))+
            cy*(pcf[6]+cx*(pcf[7]+cx*pcf[8])));
 acf[1]=pcf[1]+2*cx*pcf[2]+cy*((pcf[4]+2*cx*pcf[5])+cy*(pcf[7]+2*cx*pcf[8]));
 acf[2]=pcf[2]+cy*(pcf[5]+cy*pcf[8]);
 acf[3]=pcf[3]+cx*(pcf[4]+cx*pcf[5])+2*cy*(pcf[6]+cx*(pcf[7]+cx*pcf[8]));
 acf[4]=pcf[4]+2*cx*pcf[5]+2*cy*(pcf[7]+2*cx*pcf[8]);
 acf[5]=pcf[5]+2*cy*pcf[8];
 acf[6]=pcf[6]+cx*(pcf[7]+cx*pcf[8]);
 acf[7]=pcf[7]+2*cx*pcf[8];
 acf[8]=pcf[8];

/* Calculate the integral for the upper semicircle: */
 iupper=biquad_isc_int_pixel_semicircle(acf,-cx,1.0-cx,-cy,1.0-cy,cr);
/* Negate the coefficients of the polynomial which has odd order in 'y': */
 acf[3]=-acf[3],
 acf[4]=-acf[4],
 acf[5]=-acf[5];
/* Calculate the integral for the lower semicircle: */
 ilower=biquad_isc_int_pixel_semicircle(acf,-cx,1.0-cx, cy-1.0,cy,cr);

 return(iupper+ilower);
}

/*****************************************************************************/

double biquad_isc_int_subpixel(double **c,int ix,int iy,double x0,double y0,double x1,double y1)
{
 double	*cf[3],pcf[9],w;

 if ( x0>x1 )	w=x0,x0=x1,x1=w;
 if ( y0>y1 )	w=y0,y0=y1,y1=w;
 if ( x0<0.0 )	x0=0.0;
 if ( y0<0.0 )	y0=0.0;
 if ( x0>1.0 )	x0=1.0;
 if ( y0>1.0 )	y0=1.0;

 cf[0]=&c[2*iy+0][2*ix];
 cf[1]=&c[2*iy+1][2*ix];
 cf[2]=&c[2*iy+2][2*ix];
 biquad_poly_coefficients(cf,pcf);

 w=biquad_isc_integral_rect_definite(pcf,x0,x1,y0,y1);

 return(w);
}

static double biquad_isc_int_subpixel_nt(double **c,int ix,int iy,double x0,double y0,double x1,double y1)
{
 double	*cf[3],pcf[9],w;

 cf[0]=&c[2*iy+0][2*ix];
 cf[1]=&c[2*iy+1][2*ix];
 cf[2]=&c[2*iy+2][2*ix];
 biquad_poly_coefficients(cf,pcf);

 w=biquad_isc_integral_rect_definite(pcf,x0,x1,y0,y1);

 return(w);
}

double biquad_isc_int_rectangle(double **c,double x0,double y0,double x1,double y1)
{
 double	w,ret;
 int	sig,ix0,ix1,iy0,iy1,x,y;
 sig=1;
 if ( x0>x1 )	sig=-sig,w=x0,x0=x1,x1=w;
 if ( y0>y1 )	sig=-sig,w=y0,y0=y1,y1=w;
 ix0=(int)x0,iy0=(int)y0;
 ix1=(int)x1,iy1=(int)y1;
 x0-=(double)ix0,y0-=(double)iy0;
 x1-=(double)ix1,y1-=(double)iy1;
 if ( ix0==ix1 )
  {	if ( iy0==iy1 )
		ret=biquad_isc_int_subpixel_nt(c,ix0,iy0,x0,y0,x1,y1);
	else
	 {	ret=biquad_isc_int_subpixel_nt(c,ix0,iy0,x0,y0,x1,1.0);
		for ( y=iy0+1 ; y<iy1 ; y++ )
		 	ret+=biquad_isc_int_subpixel_nt(c,ix0,y,x0,0.0,x1,1.0);	
		if ( y1>0.0 )
			ret+=biquad_isc_int_subpixel_nt(c,ix0,iy1,x0,0.0,x1,y1);
	 }
  }
 else
  {	if ( iy0==iy1 )
	 {	ret=biquad_isc_int_subpixel_nt(c,ix0,iy0,x0,y0,1.0,y1);
		for ( x=ix0+1 ; x<ix1 ; x++ )
			ret+=biquad_isc_int_subpixel_nt(c,x,iy0,0.0,y0,1.0,y1);	
		if ( x1>0.0 )
			ret+=biquad_isc_int_subpixel_nt(c,ix1,iy0,0.0,y0,x1,y1);
	 }
	else
	 {	ret=biquad_isc_int_subpixel_nt(c,ix0,iy0,x0,y0,1.0,1.0);
		for ( x=ix0+1 ; x<ix1 ; x++ )
			ret+=biquad_isc_int_subpixel_nt(c,x,iy0,0.0,y0,1.0,1.0);
		if ( x1>0.0 )
			ret+=biquad_isc_int_subpixel_nt(c,ix1,iy0,0.0,y0,x1,1.0);
		for ( y=iy0+1 ; y<iy1 ; y++ )
		 {	ret+=biquad_isc_int_subpixel_nt(c,ix0,y,x0,0.0,1.0,1.0);
			for ( x=ix0+1 ; x<ix1 ; x++ )
				ret+=c[2*y+1][2*x+1];
			if ( x1>0.0 )
				ret+=biquad_isc_int_subpixel_nt(c,ix1,y,0.0,0.0,x1,1.0);
		 }
		if ( y1>0.0 )
		 {	ret+=biquad_isc_int_subpixel_nt(c,ix0,iy1,x0,0.0,1.0,y1);
			for ( x=ix0+1 ; x<ix1 ; x++ )
				ret+=biquad_isc_int_subpixel_nt(c,x,iy1,0.0,0.0,1.0,y1);
			if ( x1>0.0 )
				ret+=biquad_isc_int_subpixel_nt(c,ix1,iy1,0.0,0.0,x1,y1);
		 }
	 }
  }
 if ( sig>0 )	return(ret);
 else		return(-ret);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int biquad_isc_int_block_subpixels(double **c,int ix,int iy,int nx,int ny,double **ret)
{
 double	*cf[3],pcf[9],w,x0,y0,x1,y1;
 int	i,j;

 if ( nx<=0 || ny<=0 )	return(1);

 cf[0]=&c[2*iy+0][2*ix];
 cf[1]=&c[2*iy+1][2*ix];
 cf[2]=&c[2*iy+2][2*ix];
 biquad_poly_coefficients(cf,pcf);

 y1=0.0;
 for ( i=0 ; i<ny ; i++ )
  {	y0=y1;
	y1=(double)(i+1)/(double)ny;
	x1=0.0;
	for ( j=0 ; j<nx ; j++ )
	 {	x0=x1;
		x1=(double)(j+1)/(double)nx;
		w=biquad_isc_integral_rect_definite(pcf,x0,x1,y0,y1);
		ret[i][j]=w;
	 }
  }

 return(0);
}

/*****************************************************************************/

#define		STATIC_VERTICES		64

typedef struct 
 {	int	type,side;
	double	delta;
	double	x,y;
 } vertex;

static int biquad_isc_is_inside(double x1,double y1,double x2,double y2,double x3,double y3,double x,double y)
{
 if ( ((x2-x1)*(y-y1)-(x-x1)*(y2-y1)) < 0.0 )	return(0);
 if ( ((x3-x2)*(y-y2)-(x-x2)*(y3-y2)) < 0.0 )	return(0);
 if ( ((x1-x3)*(y-y3)-(x-x3)*(y1-y3)) < 0.0 )	return(0);
 return(1);
}

static int biquad_isc_it_build_side(vertex *v,double x1,double y1,double x2,double y2,int ix1,int iy1,int ix2,int iy2,int sid)
{
 int	nv;
 int	ix,iy,n;
 double	dd;
 v[0].type=0;
 v[0].side=sid,v[0].delta=0.0;
 v[0].x=x1,v[0].y=y1;

 nv=1;

 if ( ix1<ix2 )	ix=ix1+1,n=ix2-ix1;
 else 		ix=ix2+1,n=ix1-ix2;
 for ( ; n ; ix++,n--,nv++ )
  {	dd=((double)ix-x1)/(x2-x1);
	v[nv].delta=dd,
	v[nv].side=sid,
	v[nv].type=1;		/* intersection with a vertical grid line    */
	v[nv].x=(double)ix,
	v[nv].y=y1+dd*(y2-y1);
  }

 if ( iy1<iy2 )	iy=iy1+1,n=iy2-iy1;
 else 		iy=iy2+1,n=iy1-iy2;
 for ( ; n ; iy++,n--,nv++ )
  {	dd=((double)iy-y1)/(y2-y1);
	v[nv].delta=dd,
	v[nv].side=sid,
	v[nv].type=2;		/* intersection with a horizontal grid line  */
	v[nv].x=x1+dd*(x2-x1),
	v[nv].y=(double)iy;
  }
 return(nv);
}

static int biquad_isc_it_compare(const void *c1,const void *c2)
{
 if ( ((vertex *)c1)->side<((vertex *)c2)->side )		return(-1);
 else if ( ((vertex *)c1)->side>((vertex *)c2)->side )		return(1);
 else if ( ((vertex *)c1)->delta<((vertex *)c2)->delta )	return(-1);
 else if ( ((vertex *)c1)->delta>((vertex *)c2)->delta )	return(1);
 else								return(0);
}
static double biquad_isc_int_subpolygon(double *pcf,double *xcf,double *ps,int np)
{
 double	x,y,x1,y1,x2,y2;
 double	dx1,dy1,dx2,dy2,jac;
 int	i;
 double	acf[15],bcf[15],ccf[15],ret,tret,w;

 x=ps[0],y=ps[1],x1=ps[2],y1=ps[3];ps+=4;dx1=x1-x,dy1=y1-y;
 np-=2;

 tret=0.0;
 for ( ; np>0 ; x1=x2,dx1=dx2,y1=y2,dy1=dy2,ps+=2,np-- )
  {	x2=ps[0],y2=ps[1];
	dx2=x2-x,dy2=y2-y;
	jac=dx1*dy2-dx2*dy1;
	if ( jac==0.0 )	continue;

	ret=xcf[0]*pcf[0];

	acf[0]=x,
	acf[1]=dx1,
	acf[2]=dx2;
	for ( w=0.0,i=0 ; i<3 ; i++ )	w+=acf[i]*xcf[i];
	ret+=pcf[1]*w;

	bcf[0]=x*x,
	bcf[1]=2*x*dx1,
	bcf[2]=2*x*dx2,
	bcf[3]=dx1*dx1,
	bcf[4]=2*dx1*dx2,
	bcf[5]=dx2*dx2;
	for ( w=0.0,i=0 ; i<6 ; i++ )	w+=bcf[i]*xcf[i];
	ret+=pcf[2]*w;

	acf[0]=y;
	acf[1]=dy1;
	acf[2]=dy2;
	for ( w=0.0,i=0 ; i<3 ; i++ )	w+=acf[i]*xcf[i];
	ret+=pcf[3]*w;

	bcf[0]=x*acf[0],
	bcf[1]=x*acf[1]+dx1*acf[0],
	bcf[2]=x*acf[2]+dx2*acf[0],
	bcf[3]=dx1*acf[1],
	bcf[4]=dx1*acf[2]+dx2*acf[1],
	bcf[5]=dx2*acf[2];
	for ( w=0.0,i=0 ; i<6 ; i++ )	w+=bcf[i]*xcf[i];
	ret+=pcf[4]*w;

	ccf[0]=x*bcf[0],
	ccf[1]=x*bcf[1]+dx1*bcf[0],
	ccf[2]=x*bcf[2]+dx2*bcf[0],
	ccf[3]=x*bcf[3]+dx1*bcf[1],
	ccf[4]=x*bcf[4]+dx1*bcf[2]+dx2*bcf[1],
	ccf[5]=x*bcf[5]+dx2*bcf[2],
	ccf[6]=dx1*bcf[3],
	ccf[7]=dx1*bcf[4]+dx2*bcf[3],
	ccf[8]=dx1*bcf[5]+dx2*bcf[4],
	ccf[9]=dx2*bcf[5];
	for ( w=0.0,i=0 ; i<10 ; i++ )	w+=ccf[i]*xcf[i];
	ret+=pcf[5]*w;

	acf[0]=y*y,
	acf[1]=2*y*dy1,
	acf[2]=2*y*dy2,
	acf[3]=dy1*dy1,
	acf[4]=2*dy1*dy2,
	acf[5]=dy2*dy2;
	for ( w=0.0,i=0 ; i<6 ; i++ )	w+=acf[i]*xcf[i];
	ret+=pcf[6]*w;

	bcf[0]=x*acf[0],
	bcf[1]=x*acf[1]+dx1*acf[0],
	bcf[2]=x*acf[2]+dx2*acf[0],
	bcf[3]=x*acf[3]+dx1*acf[1],
	bcf[4]=x*acf[4]+dx1*acf[2]+dx2*acf[1],
	bcf[5]=x*acf[5]+dx2*acf[2],
	bcf[6]=dx1*acf[3],
	bcf[7]=dx1*acf[4]+dx2*acf[3],
	bcf[8]=dx1*acf[5]+dx2*acf[4],
	bcf[9]=dx2*acf[5];
	for ( w=0.0,i=0 ; i<10 ; i++ )	w+=bcf[i]*xcf[i];
	ret+=pcf[7]*w;

	ccf[0]=x*bcf[0],
	ccf[1]=x*bcf[1]+dx1*bcf[0],
	ccf[2]=x*bcf[2]+dx2*bcf[0],
	ccf[3]=x*bcf[3]+dx1*bcf[1],
	ccf[4]=x*bcf[4]+dx1*bcf[2]+dx2*bcf[1],
	ccf[5]=x*bcf[5]+dx2*bcf[2],
	ccf[6]=x*bcf[6]+dx1*bcf[3],
	ccf[7]=x*bcf[7]+dx1*bcf[4]+dx2*bcf[3],
	ccf[8]=x*bcf[8]+dx1*bcf[5]+dx2*bcf[4],
	ccf[9]=x*bcf[9]+dx2*bcf[5],
	ccf[10]=dx1*bcf[6],
	ccf[11]=dx1*bcf[7]+dx2*bcf[6],
	ccf[12]=dx1*bcf[8]+dx2*bcf[7],
	ccf[13]=dx1*bcf[9]+dx2*bcf[8],
	ccf[14]=dx2*bcf[9];
	for ( w=0.0,i=0 ; i<15 ; i++ )	w+=ccf[i]*xcf[i];
	ret+=pcf[8]*w;

	tret+=ret*jac;
	
  }
 return(tret);
}

#define	int_round(x,i) do { if ( fabs(x-i)<1e-8 ) x=i;else if ( fabs(x-1.0-i)<1e-8 ) x=i+1,i++; } while (0)

double biquad_isc_int_triangle(double **c,int is_spline,double x1,double y1,double x2,double y2,double x3,double y3,int sx,int sy)
{
 int	ix1,iy1,ix2,iy2,ix3,iy3;
 int	minx,maxx,miny,maxy,ix,iy,i,j,k,ii;
 int	nv,cv,iv,np;
 vertex	vert_stat[STATIC_VERTICES],*vert_dyn,*vert,*v,vl[8];
 double	ret,x,y,points[20],pcf[9],*cf[3],det,xcf[15],fa,fb;

 if ( c==NULL )		return(0.0);

 /* Re-arrange the point into counter-clockwise order */
 det=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
 if ( det==0.0 )	return(0.0);
 else if ( det<0.0 )	x=x3,x3=x2,x2=x,y=y3,y3=y2,y2=y;

 ix1=(int)x1,ix2=(int)x2,ix3=(int)x3;
 iy1=(int)y1,iy2=(int)y2,iy3=(int)y3;
/* int_round(x1,ix1);
 int_round(x2,ix2);
 int_round(x3,ix3);
 int_round(y1,iy1);
 int_round(y2,iy2);
 int_round(y3,iy3);*/

/* fprintf(stdout,"%.16f %.16f %.16f %.16f %.16f %.16f\n",x1,y1,x2,y2,x3,y3); */

 nv=3+	abs(ix1-ix2)+abs(ix2-ix3)+abs(ix3-ix1)+
	abs(iy1-iy2)+abs(iy2-iy3)+abs(iy3-iy1);
 if ( nv<=STATIC_VERTICES )
  {	vert_dyn=NULL,
	vert=vert_stat;
  }
 else
  {	vert_dyn=(vertex *)malloc(sizeof(vertex)*nv),
	vert=vert_dyn;
  }

 cv=0;
 cv+=biquad_isc_it_build_side(vert+cv,x1,y1,x2,y2,ix1,iy1,ix2,iy2,0);
 cv+=biquad_isc_it_build_side(vert+cv,x2,y2,x3,y3,ix2,iy2,ix3,iy3,1);
 cv+=biquad_isc_it_build_side(vert+cv,x3,y3,x1,y1,ix3,iy3,ix1,iy1,2);
 if ( cv > nv )	/* inconsistency */
  {	return(0.0);		}
 else		nv=cv;
 for ( i=0 ; i<nv ; i++ )
  {	if ( vert[i].delta>=1.0 || ( vert[i].delta<=0.0 && vert[i].type>0 ) )
		vert[i].type=-1;
  }
 qsort(vert,nv,sizeof(vertex),biquad_isc_it_compare);
 for ( i=0 ; i<nv ; i++ )
  {	if ( 0<i && vert[i-1].side==vert[i].side && vert[i-1].type>0 &&
	vert[i].type>0 && vert[i-1].delta==vert[i].delta )
		vert[i].type=-1;
  }

/* fprintf(stderr,"biit(): db0: (%g,%g) (%g,%g) (%g,%g)\n",x1,y1,x2,y2,x3,y3);
 for ( i=0 ; i<nv ; i++ )
  {	fprintf(stderr,"biit(): db0: [%2d] %d %10.14g (%10g,%10g)\n",
		vert[i].type,vert[i].side,vert[i].delta,vert[i].x,vert[i].y);
  }*/  /*DB0*/
 
 minx=ix1;if ( ix2<minx ) minx=ix2;if ( ix3<minx ) minx=ix3;
 maxx=ix1;if ( ix2>maxx ) maxx=ix2;if ( ix3>maxx ) maxx=ix3;
 miny=iy1;if ( iy2<miny ) miny=iy2;if ( iy3<miny ) miny=iy3;
 maxy=iy1;if ( iy2>maxy ) maxy=iy2;if ( iy3>maxy ) maxy=iy3;

 xcf[0]=fa=0.5;k=1;
 for ( i=1 ; i<=4 ; i++ )
  {	fa=fa*(double)i/(double)(i+2);
	for ( j=0,fb=fa ; j<=i ; j++,k++ )
	 {	xcf[k]=fb;
		if ( j<i )	fb=fb*(double)(j+1)/(double)(i-j);
	 }
  }

 ret=0.0;fb=0.0;
 for ( iy=miny ; iy<=maxy ; iy++ )
  {  if ( sy>=0 && ( iy<0 || iy>=sy ) )	continue;
     for ( ix=minx ; ix<=maxx ; ix++ )
      {	if ( sx>=0 && ( ix<0 || ix>=sx ) ) continue;
	for ( i=0,iv=0,v=vert ; i<nv ; i++,v++ )
	 {	if ( v->type==0 && ix<=v->x && v->x<=ix+1 
		&& iy<=v->y && v->y<=iy+1 )
		 {	vl[iv].type=0;
			x=vl[iv].x=v->x-ix,
			y=vl[iv].y=v->y-iy;
			if ( x<=0.0 )		vl[iv].type=1;
			else if ( x>=1.0 )	vl[iv].type=3;
			else if ( y<=0.0 )	vl[iv].type=2;
			else if ( y>=1.0 )	vl[iv].type=4;
			vl[iv].side=v->side;
			vl[iv].delta=0.0;
			iv++;
		 }
		else if ( v->type==1 && iy<=v->y && v->y<=iy+1 )
		 {	if ( v->x==ix )
			 {	vl[iv].type=1,
				vl[iv].x=0.0,
				vl[iv].y=v->y-iy;
				vl[iv].side=v->side;
				vl[iv].delta=v->delta;
				iv++;
			 }
			if ( v->x==ix+1 )
			 {	vl[iv].type=3,
				vl[iv].x=1.0,
				vl[iv].y=v->y-iy;
				vl[iv].side=v->side;
				vl[iv].delta=v->delta;
				iv++;
			 }
		 }
		else if ( v->type==2 && ix<=v->x && v->x<=ix+1 )
		 {	if ( v->y==iy )
			 {	vl[iv].type=2,
				vl[iv].x=v->x-ix,
				vl[iv].y=0.0;
				vl[iv].side=v->side;
				vl[iv].delta=v->delta;
				iv++;
			 }
			if ( v->y==iy+1 )
			 {	vl[iv].type=4,
				vl[iv].x=v->x-ix,
				vl[iv].y=1.0;
				vl[iv].side=v->side;
				vl[iv].delta=v->delta;
				iv++;
			 }
		 }
	 }
	if ( iv<=1 )
	 {	x=(double)ix+0.5,
		y=(double)iy+0.5;
		/* is (x,y) inside the triangle or not...?! */
		ii=biquad_isc_is_inside(x1,y1,x2,y2,x3,y3,x,y);
		if ( ii )
		 {	if ( is_spline )	ret+=c[2*iy+1][2*ix+1];
			else			ret+=c[iy][ix];
		 }
		/*fprintf(stderr,"[%d,%d] [%d]\n",ix,iy,ii);*/ /*DB1*/
	 }
	else
	 {	vertex	*pv,*vv;
		double	*pp;
		int	k;
		pv=&vl[iv-1];
		pp=&points[0];
		pp[0]=pv->x,pp[1]=pv->y;
		np=1;pp+=2;
		for ( i=0,vv=vl ; i<iv ; i++,pv=vv,vv++ )
		 {	if ( pv->delta<=0.0 && ((pv->side+1)%3==(vv->side)%3) )
			 {	pp[0]=vv->x,pp[1]=vv->y;
				pp+=2;np++;
			 }
			else if ( vv->delta<=0.0 && ((pv->side+1)%3==(vv->side)%3) )
			 {	pp[0]=vv->x,pp[1]=vv->y;
				pp+=2;np++;
			 }
			else if ( pv->side==vv->side && pv->delta<vv->delta )
			 {	pp[0]=vv->x,pp[1]=vv->y;
				pp+=2;np++;
			 }
			else
			 {	k=pv->type;
				while ( k%4 != (vv->type)%4 )
				 {	if ( k%4==1 )
						pp[0]=0.0,pp[1]=0.0;
					if ( k%4==2 )
						pp[0]=1.0,pp[1]=0.0;
					if ( k%4==3 )
						pp[0]=1.0,pp[1]=1.0;
					if ( k%4==0 )
						pp[0]=0.0,pp[1]=1.0;
					pp+=2;np++,k++;
				 };	
				pp[0]=vv->x,pp[1]=vv->y;
				pp+=2;np++;
			 }
		 }
		if ( points[2*np-2]==points[0] && points[2*np-1]==points[1] )
			np--;

		/*fprintf(stdout,"[np=%d]\n",np);
		fprintf(stdout,"[%d,%d]",ix,iy);
		for ( k=0 ; k<np ; k++ )
		 {	fprintf(stdout," (%g,%g)",points[2*k],points[2*k+1]);	}
		fprintf(stdout,"\n");*/ /*DB2*/

		if ( np<=2 )	np=0;

		/*fa=points[2*(np-1)]*points[1]-points[2*(np-1)+1]*points[0];
		for ( k=0 ; k<2*(np-1) ; k+=2 )
			fa+=points[k]*points[k+3]-points[k+1]*points[k+2];
		fb+=fa;*/ /*DB3*/

		if ( is_spline )
		 {	cf[0]=&c[2*iy+0][2*ix];
			cf[1]=&c[2*iy+1][2*ix];
			cf[2]=&c[2*iy+2][2*ix];
			biquad_poly_coefficients(cf,pcf);
			ret+=biquad_isc_int_subpolygon(pcf,xcf,points,np);
		 }
		else
		 {	fa=points[2*(np-1)]*points[1]-points[2*(np-1)+1]*points[0];
			for ( k=0 ; k<2*(np-1) ; k+=2 )
				fa+=points[k]*points[k+3]-points[k+1]*points[k+2];
			ret+=c[iy][ix]*0.5*fa;
		 }
	 }
	
      }
  }

/* fprintf(stdout,"A=%g\n\n",fb);*/ /*DB3*/

 if ( vert_dyn != NULL )	free(vert_dyn);
 
 return(ret);
}

/*****************************************************************************/
                                       
