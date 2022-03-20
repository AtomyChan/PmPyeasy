/*****************************************************************************/
/* splinefit.h								     */
/*****************************************************************************/

#ifndef	__SPLINEFIT_H_INCLUDED
#define	__SPLINEFIT_H_INCLUDED	1

#include "point.h"

/*****************************************************************************/

/* fit_1d_spline():
   This function fits a natural cubic spline on the interval [x0,x1] to the
   data points 'xv','yv' ('npoint' points, altogether). The spline will
   have 'bx' equidistant sections, therefore it has 'bx'+1 control points
   (at x=x0, x0+(x1-x0)/bx, x0+2*(x1-x0)/bx, ..., x1). The spline values
   at the control points are returned in 'coeff'. Note that in the array
   'coeff' there should be enough space for 'bx'+1 points. If the fit
   fails (e.g. 'npoint' is less than 'bx'+1 or x0=x1, or xv, yv or coeff
   is NULL) the function returns a non-zero value, otherwise it returns 0.
   The array 'wv' can be used to define weights to the data points.          */
int fit_1d_spline(double x0,double x1,int bx,
	double *xv,double *yv,int npoint,double *wv,double *coeff);

/* fit_2d_spline():
   This function fits a natural bicubic spline on the reactangle
   [x0,x1] x [y0,y1] to the data points described by the pointset 'points'
   (with 'npoint' pont). The bicubic spline will be defined on 'bx' times 'by'
   equilateral rectangles. The control points of the fitted spline 
   will be returned in 'coeff'. Note that in the array 'coeff'
   there should be enough space for 'bx'+1 times 'by'+1 points. If the fit
   fails (e.g. 'npoint' is less than 'bx'+1 times 'by'+1 or x0=x1, or y0=y1,
   or either 'points' or 'coeff' is NULL) the function returns a non-zero 
   value, otherwise it returns 0.					     */
int fit_2d_spline(double x0,double x1,int bx,double y0,double y1,int by,
	point *points,int npoint,double **coeff);

/*****************************************************************************/

#endif

/*****************************************************************************/
                   
