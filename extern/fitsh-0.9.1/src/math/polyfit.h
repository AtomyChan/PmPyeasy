/*****************************************************************************/
/* polyfit.h								     */
/*****************************************************************************/

#ifndef __POLYFIT_H_INCLUDED
#define	__POLYFIT_H_INCLUDED	1

#include "point.h"
#include "poly.h"

/* fit_2d_poly():
   Fits a 2 dimensional polynomial up to order 'order' to the points 'dat'.
   The number of the points 'nd' should be greater or equal than the
   number of monomials to the order 'order':
	nd >= (order+1)*(order+2)/2,					(1)
   otherwise the fit cannot be performed (the least-squares matrix will be
   singular). Before the fit, all points are shifted and scaled using 
   the parameters 'ox', 'oy' and 'scale'. If the fit was successful, the 
   coefficients of the fitted polynomial is retured in the array 'fit'
   and 0 is returned. If the allocation of the temporary resources fails,
   the function returns a negative value. If the fit cannot be performed
   because there are not enough points (so, (1) is not satisfied) or the 
   least-squares matrix is singular, the function returns a positive value.  */
int	fit_2d_poly   (	point *dat,int nd,int order,double *fit,
			double ox,double oy,double scale);

/* fit_2d_leg_poly():
   Fits a 2 dimensional Legendre polynomial up to order 'order' to the 
   points 'dat'. See also fit_2d_poly().				     */
int	fit_2d_leg_poly(point *dat,int nd,int order,double *fit,
			double ox,double oy,double scale);

/* fit_1d_poly():
   Fits an 1 dimensional non-normalized polynomial up to order 'order' to
   the points (x,y)[]. The number of points is 'nd'. The fitted coefficients
   are returned in 'coeff' (order+1 values). If nd is less or equal than
   order, the fit is not performed and the function returns a non-zero value.*/
int	fit_1d_poly(int order,double *x,double *y,int nd,double *w,double *coeff);

#endif
                  
