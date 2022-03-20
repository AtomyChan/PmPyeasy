/*****************************************************************************/
/* xfunct.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Some extra functions for extending the default 'lfit' functionality:	     */
/*	- eccentric_offset_{q,p}(l,k,h) -- eop(), eoq(): eccentric offset    */
/*	  functions and their derivatives (these functions are analytic for  */
/*	  k^2 + h^2 < 1 and for all values of l).			     */
/*	- HJD and BJD functions						     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2007; Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <astro/astro.h>

#include "xfunct.h"

/*****************************************************************************/

double eccentric_offset_q(double lambda,double k,double h)
{
 double	e,E;
 e=sqrt(h*h+k*k);
 if ( e<=0.0 )
	return(0.0); 
 else
  {	E=solve_kepler_equ(lambda-atan2(h,k),e);
	return(e*cos(E));
  }
}
double eccentric_offset_p(double lambda,double k,double h)
{
 double	e,E;
 e=sqrt(h*h+k*k);
 if ( e<=0.0 )
	return(0.0); 
 else
  {	E=solve_kepler_equ(lambda-atan2(h,k),e);
	return(e*sin(E));
  }
 
}

double eccentric_trigonometric_c(double lambda,double k,double h)
{
 return(cos(lambda+eccentric_offset_p(lambda,k,h)));
}

double eccentric_trigonometric_s(double lambda,double k,double h)
{
 return(sin(lambda+eccentric_offset_p(lambda,k,h)));
}

/*****************************************************************************/

double get_hjd(double jd,double ra,double dec)
{
 double	hjd;

 hjd=get_heliocentric_julian_date(jd,ra,dec);

 return(hjd);
}

double get_bjd(double jd,double ra,double dec)
{
 double	bjd;

 bjd=get_barycentric_julian_date(jd,ra,dec);

 return(bjd);
}

/*****************************************************************************/
                          
