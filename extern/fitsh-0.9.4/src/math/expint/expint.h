/*****************************************************************************/
/* expint.h								     */
/*****************************************************************************/

#ifndef	__EXPINT_H_INCLUDED
#define	__EXPINT_H_INCLUDED	1

typedef struct
 {	double	exp1,exp2;
	double	erf1,erf2;
 } expinte;

typedef struct 
 {	double	expx1,expx2,expy1,expy2;
	double	erfx1,erfx2,erfy1,erfy2;
 } expintee;

/* expint_numerical(), expint_taylor(), expint_taylor_ee():
   Calculates the definite integral of the general elliptical Gaussian function
   exp[-1/2*((s+d)x^2+(s-d)y^2+k*(2xy))] on the rectangle [x1,x2]x[y1,y2]. 
   If the calculation is successful, the desired integral is returned
   (which always is positive if x1<x2 and y1<y2), otherwise the function 
   returns -1. Expint_taylor_ee() assumes the values of exp() and erf()
   function at the points -0.5*T*T*(s+d) and -0.5*T*T*(s-d) (where 
   T=x1 or x2 for `s+d` and T=y1 or y2 for `s-d`) are previously calculated
   and stored in the structure 'ee' (the function expint_taylor() just 
   calculates these values, store them in a temporary array and simply 
   calls expint_taylor_ee(), see the source code for details).	             */
double	expint_taylor_ee(double s,double d,double k,
	double x1,double x2,double y1,double y2,expintee *ee);
double	expint_taylor(double s,double d,double k,
	double x1,double x2,double y1,double y2);
/* expint_numerical():
   Calculates the definite integral of the general elliptical Gaussian function
   exp[-1/2*((s+d)x^2+(s-d)y^2+k*(2xy))] on the rectangle [x1,x2]x[y1,y2]
   using Simphson's algorithm. The build-in precision (grid size) is 
   approximately for 7-8 significant figures.				     */
double	expint_numerical(double s,double d,double k,
	double x1,double x2,double y1,double y2);

/* expint_taylor_ee_diff(), expint_taylor_diff():
   Calculates the definite integral of the general elliptical Gaussian function
   exp[-1/2*((s+d)x^2+(s-d)y^2+k*(2xy))]  and its derivatives (by the 
   parameters x0, y0, s, d and k) on the rectangle [x1,x2]x[y1,y2]. 
   The integral and the integrals of the derivatives (which are the derivatives
   of the integrals also...), totally six numbers are stored in the above
   order in the array dlist[].						     */
int	expint_taylor_ee_diff(double s,double d,double k,
	double x1,double x2,double y1,double y2,double *dlist,expintee *ee);
int	expint_taylor_diff(double s,double d,double k,
	double x1,double x2,double y1,double y2,double *dlist);

/* expint_taylor_ee_shift_diff():
   Calculates the definite integral of the general elliptical Gaussian function
   exp[-1/2*((s+d)x'^2+(s-d)y'^2+k*(2x'y'))] and its derivatives (by the 
   parameters x0, y0, s, d and k) on the rectangle [x1,x2]x[y1,y2],
   where x'=x-x0, y'=y-y0. This function simply subtract x0 and y0 from
   x1, x2 and y1, y2, respectively and calls expint_taylor_ee_diff().	     */
int	expint_taylor_ee_shift_diff(double s,double d,double k,double x0,double y0,
	double x1,double x2,double y1,double y2,double *dlist,expintee *ee);

/* expint_list(), expint_list_e():
   Stores the definite integrals of the functions (x^i)*exp(-1/2*s*x^2)
   (for i=0..n) on the interval [x1,x2] and stores the result in
   (totally n+1 numbers) in the array list[]. Note that this method
   is faster than using expint_primitve_list() and subtraction.		     */
int	expint_list_e(double s,double x1,double x2,int n,double *list,expinte *e);
int	expint_list(double s,double x1,double x2,int n,double *list);

/* expint_primitive_list():
   Calculates the primitive functions of (x^i)*exp(-1/2*s*x^2) (for i=0..n)
   at the point 'x' and stores (totally n+1 numbers) in the array list[].    */
int	expint_primitive_list(double s,double x,int n,double *list);

#endif
             
