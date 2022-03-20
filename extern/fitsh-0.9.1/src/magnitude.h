/*****************************************************************************/
/* magnitude.h								     */
/*****************************************************************************/

#ifndef	__MAGNITUDE_H_INCLUDED
#define	__MAGNITUDE_H_INCLUDED	1

/*****************************************************************************/

typedef struct
 {	double	magnitude;
	double	intensity;
 } magflux;

double	mag_to_flux(double magn,magflux *mf);
double	flux_to_mag(double flux,magflux *mf);

int	flux_to_mag_magerr(double f,double fe,magflux *mf,double *rm,double *rme);

/*****************************************************************************/

#endif

/*****************************************************************************/
                                                      
