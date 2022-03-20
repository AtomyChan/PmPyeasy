/*****************************************************************************/
/* easpec.h 								     */
/*****************************************************************************/

#ifndef	__EASPEC_H_INCLUDED
#define	__EASPEC_H_INCLUDED	1

/* get_observer_coords():
   Calculates the geocentric equatorial Cartesian coordiantes of the observer
   located at (lon,lat) at the time instance 'jdate'.			     */
void	get_observer_coords(double jdate,double lon,double lat,
	double *x,double *y,double *z);

/* get_sideral_time():
  Calculates and returns the sideral time in Greenwih at the instance 'jdate'*/
double	get_sideral_time(double jdate);

/* get_earth_axis_angle():
   Calculates and returns the deviation angle of the rotational axis of the
   Earth to the norm of the ecliptical plane at the given instance 'jdate'.
   (Note that this is a non-covariant quantity, so this function might be
   replaced or changed in future releases of this library...)		     */
double	get_earth_axis_angle(double jdate);

/* get_heiocentric_julian_date():
   Calculates the heliocentric (not barycentric!) correction for the
   time instance 'jdate' at the spherical equatorial point ('ra', 'dec').
   The corrected time instance is returned.				     */
double	get_heliocentric_julian_date(double jdate,double ra,double dec);

/* get_heiocentric_julian_date_epoch():
   Same as get_heiocentric_julian_date(), but the coordinates 'ra' and 'dec'
   are for the epoch (in Julian years) 'epoch'.				     */
double	get_heliocentric_julian_date_epoch(double jdate,double ra,double dec,double y0);

/* get_barycentric_julian_date():
   Calculates the barycentric correction for the time instance 'jdate' at the 
   spherical equatorial point ('ra', 'dec'). The corrected time instance is 
   returned.				   				     */
double	get_barycentric_julian_date(double jdate,double ra,double dec);

/* get_annular_parallax_correction():
   This function calculates the annular parallax correction (i.e. as it 
   would observed from the Sun), for the celestial object located at 
   (ra,dec) (RA, DEC, in degrees) within a distance of 'dis' parsec at 
   the instance 'jdate' (JD). The corrected coordinates are stored in 'rra'
   and 'rdec'. 								     */
int	get_annular_parallax_correction(double jdate,double ra,double dec,
	double dis,double *rra,double *rdec);

#endif
                                    
