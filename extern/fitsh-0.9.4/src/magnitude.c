/*****************************************************************************/
/* magnitude.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Some simple functions which converts between magnitudes and fluxes.	     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2006; Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "magnitude.h"

/*****************************************************************************/

double mag_to_flux(double mag,magflux *mf)
{
 double	flux;
 
 flux=mf->intensity*exp(-0.4*M_LN10*(mag-mf->magnitude));

 return(flux);
}

double flux_to_mag(double flux,magflux *mf)
{
 double	mag;

 if ( flux<=0.0 || mf->intensity<=0 )
	return(0.0);
 else
  {	mag=(-2.5*log10(flux/mf->intensity))+mf->magnitude;
	return(mag);
  }
}

/*****************************************************************************/

int flux_to_mag_magerr(double f,double fe,magflux *mf,double *rm,double *rme)
{
 double	mg,me;
 int	ret;

 if ( f <= 0.0 || mf->intensity <= 0.0 )
  {	mg=0.0,
	me=0.0;
	ret=1;
  }
 else
  {	mg=flux_to_mag(f,mf);
	me=1.08574 * (fabs(fe)/f);	/* 1.08574... = 5/ln(100)	*/
	ret=0;
  }

 if ( rm  != NULL )	*rm =mg;
 if ( rme != NULL )	*rme=me;

 return(ret);
}

/*****************************************************************************/
      
