/*****************************************************************************/
/* xfunct.h								     */
/*****************************************************************************/

#ifndef	__XFUNCT_H_INCLUDED
#define	__XFUNCT_H_INCLUDED	1

/*****************************************************************************/

/* eccentric_offset_q(), eccentric_offset_p():
   These functions calculate the components of the eccentric offset
   vector. The eccentric offset vector (q,p) is e*(cos E,sin E) where
   e^2 = k^2 + h^2 and E is the solution of E - e*sin E = lambda-arg(k,h).   */
double	eccentric_offset_q(double lambda,double k,double h);
double	eccentric_offset_p(double lambda,double k,double h);

double	eccentric_trigonometric_c(double lambda,double k,double h);
double	eccentric_trigonometric_s(double lambda,double k,double h);

/*****************************************************************************/

double	elliptic_complete_first(double k);
double	elliptic_complete_second(double k);
double	elliptic_complete_third(double n,double k);

/*****************************************************************************/

/* wrapper to get_heiocentric_julian_date(): */
double	get_hjd(double jd,double ra,double dec);
/* wrapper to get_barycentric_julian_date(): */
double	get_bjd(double jd,double ra,double dec);

/*****************************************************************************/

#endif

/*****************************************************************************/
                                                   
                                     
