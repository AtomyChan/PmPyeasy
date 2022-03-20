/*****************************************************************************/
/* earth.h 								     */
/*****************************************************************************/

#ifndef	__EARTH_H_INCLUDED
#define	__EARTH_H_INCLUDED	1

/* get_earth_coords():
   Returns the heliocentric ecliptical Cartesian coordiantes of the barycentre 
   of the Earth-Moon system at the given instance 'jdate'. The coordinates 
   are in astronomical units.                                                */
int	get_earth_coords(double jdate,double *rx,double *ry,double *rz);

#endif
                  
