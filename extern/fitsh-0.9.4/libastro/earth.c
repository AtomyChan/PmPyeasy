/*****************************************************************************/
/* earth.c 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Calculates the heliocentric ecliptical Cartesian coordinates of Earth.    */
/* This part of the library is not standalone (depends on cmath.[ch]).	     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2003, 2006; Pal, A. (apal@szofi.elte.hu).			     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <limits.h>

#include <astro/cmath.h>
#include <astro/earth.h>

/*****************************************************************************/

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841968
#endif

/*****************************************************************************/

int get_earth_coords(double jdate,double *rx,double *ry,double *rz)
{
 double		jc1,jc2,jc3,x;
 double		l0,m0,e0,v0,th,em,o0;
 double		la,dis;
 static double	chx,chy,chz,cjdate;
 static int	is_valid_cache=0;

 if ( is_valid_cache && cjdate==jdate )
  {	*rx=chx;*ry=chy;*rz=chz;
	return(0);
  }

 jc1=(jdate-2415020.0)/36525.0;
 jc2=jc1*jc1;jc3=jc2*jc1;

 l0=279.69668+36000.76892*jc1+0.0003025*jc2;
 m0=358.47583+35999.04975*jc1-0.0001500*jc2-0.0000033*jc3;
 e0=0.01675101-0.0000418*jc1-0.000000126*jc2;

 em=solve_kepler_equ(m0*M_PI/180.0,e0);
 v0=sqrt((1.0+e0)/(1.0-e0))*tan(em/2.0);

 v0=2*atan(v0)*180.0/M_PI;
 th=l0+v0-m0+360.0;
 dis=1.0000002*(1-e0*cos(em));

 o0=(259.18-  1934.1420*jc1)*M_PI/180.0;
 th=th-0.00569-0.00479*sin(o0);
 
 x =(153.23+ 22518.7541*jc1)*M_PI/180.0,
 th+=0.00134*cos(x),dis+=0.00000543*sin(x);

 x =(216.57+ 45037.5082*jc1)*M_PI/180.0,
 th+=0.00154*cos(x),dis+=0.00001575*sin(x);

 x =(312.68+ 32964.3577*jc1)*M_PI/180.0,
 th+=0.00200*cos(x),dis+=0.00001627*sin(x);

 x =(350.74+445267.1142*jc1-0.00144*jc2)*M_PI/180.0;
 th+=0.00179*sin(x),dis+=0.00003076*cos(x);

 x=(231.19+   20.2000*jc1)*M_PI/180.0,th +=0.00178000*sin(x);
 x=(353.40+65928.7155*jc1)*M_PI/180.0,dis+=0.00000927*sin(x);

 la=(th-180.0)*M_PI/180.0;

 chx=*rx=dis*cos(la);
 chy=*ry=dis*sin(la);
 chz=*rz=0.0;

 is_valid_cache=1;
 cjdate=jdate;

 return(0);
}

/*****************************************************************************/
                                               
