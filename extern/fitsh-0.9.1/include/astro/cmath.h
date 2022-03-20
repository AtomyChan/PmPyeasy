/*****************************************************************************/
/* cmath.h 								     */
/*****************************************************************************/

#ifndef	__CMATH_H_INCLUDED
#define	__CMATH_H_INCLUDED	1

typedef	double	vector[3];
typedef	double	matrix[3][3];

/* sina(), cosa(), asina(), acosa():
   Funcions like sin(), cos(), asin() and acos(), but all angles are in degs.*/
double	sina(double deg);
double	cosa(double deg);
double	asina(double s);
double	acosa(double c);

/* rotate():
   Rotates (*x,*y) counter-clockwise by the specified angle 'phi' (radians). */
void	rotate(double *x,double *y,double phi);

/* getpcoords(), get3dpcoords().
   Calculates the polar coordinate of the 2d point (x,y) or 3d point (x,y,z).
   All angles returned are in degrees (not in radians).			     */
double	getpcoords(double x,double y);
void	get3dpcoords(double x,double y,double z,double *longi,double *latit);

/* get_angular_distance():
   Calculates the angular separation (in degrees) between the two spherical
   points (lon1,lat1) and (lon2,lat2). These spherical coordinates should
   also be specified in degrees.					     */
double	get_angular_distance(double lon1,double lat1,double lon2,double lat2);

/* hypot3d():
   Calculates the Eucledian distance of the point (x,y,z) from the origin.   */

/* solve_kepler_equ():
   Solves the Equation of Kepler for the mean anomaly 'm' and eccentricity 
   'ex'. All angles ('m' and the returned eccentric anomaly) are in radians. */
double	solve_kepler_equ(double m,double ex);

/* matrix_rotation():
   Creates (and stores in 'm') a rotational matrix which rotates by the angle
   'phi' (it is in degrees) around the axis 'ax' ('ax' is 0 for x, 1 for y and 
   2 for z-axis).							     */
void	matrix_rotation(matrix m,int ax,double phi);

/* matrix_mul():
   Calculates 'm1'*'m2' and stores the result in 'r'.			     */
void	matrix_mul(matrix r,matrix m1,matrix m2);

#endif
                                
