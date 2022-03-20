/*****************************************************************************/
/* spline.h								     */
/*****************************************************************************/

#ifndef	__SPLINE_H_INCLUDED
#define	__SPLINE_H_INCLUDED	1

/*****************************************************************************/

/* spline_coeff():
   Calculates the cubic spline coefficients for the data points y[]
   defined in the domain x[]. Both of the sets should contain 'n' points. The 
   spline coefficients are returned  in 'y2' ('n' real numbers). The 'yp1' 
   and 'ypn' pointers should point to two real numbers, to the second 
   derivatives for the fitted spline in the first (x[0]) and last (x[n-1]) 
   point. If any of these pointers are NULL, the natural spline approximation 
   is used (so the second derivative in the first and the last point is 0.   */
int	spline_coeff(double *x,double *y,int n,double *yp1,double *ypn,double *y2);
/* spline_inter(), spline_inter_der():
   Interpolates the x[]-y[] function and its derivative in the point 'xp' 
   using a cubic spline interpolation with the coefficiens 'y2' (calculated 
   by spline_coeff()).							     */
double	spline_inter(double *x,double *y,double *y2,int n,double xp);
double	spline_inter_der(double *x,double *y,double *y2,int n,double xp);

/* natspline_coeff():
   Calculated the natural cubic spline coefficiens of the data points y[] 
   defined in the domain 0,...,n-1. The spline coefficients are returned 
   in 'y2'. The coefficients calculated are equal to the coefficients returned 
   by spline_coeff() if x[i]=i (i=0,...,n-1) and yp1=ypn=NULL.		     */
int	natspline_coeff(double *y,int n,double *y2);
/* natspline_inter(), natspline_inter_der():
   Interpolates the y[i] function and its derivative in the point 'x' using 
   a (natural) cubic spline interpolation with the coefficients 'y2'.	     */
double	natspline_inter(double *ya,double *y2a,int n,double x);
double	natspline_inter_der(double *ya,double *y2a,int n,double x);

/* cyspline_coeff():
   Calculated the cyclic cubic spline coefficiens of the data points y[] 
   defined in the domain 0,...,n-1. The spline coefficients are returned in 
   the array 'y2'. 							     */
int	cyspline_coeff(double *y,int n,double *y2);
/* cyspline_inter(), cyspline_inter_der():
   Interpolates the y[i] function and its derivative in the point 'x' using 
   a cyclic cubic spline interpolation with the coefficients 'y2'.	     */
double	cyspline_inter(double *ya,double *y2a,int n,double x);
double	cyspline_inter_der(double *ya,double *y2a,int n,double x);

int	intspline_coeff(double *y,int p,double *l);
double	intspline_inter(double *y,double *l,int n,double x,double len);

/*****************************************************************************/

#endif

/*****************************************************************************/
  
