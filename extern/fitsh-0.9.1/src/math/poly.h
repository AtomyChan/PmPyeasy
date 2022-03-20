/*****************************************************************************/
/* poly.h								     */
/*****************************************************************************/

#ifndef __POLY_H_INCLUDED
#define	__POLY_H_INCLUDED	1

/* eval_2d_monoms():
   Evaluates 2 dimensional monomials up to order 'order' in the point (x,y).
   Before the evaluation, the point (x,y) is shifted and rescaled using the 
   parameters 'ox', 'oy' and 'scale'. Not to perform such shift and scaling,
   set ox=oy=0.0 and scale=1.0. The monoms are returned in the array 'ret'
   where should be enough space for the monomials. Up to order 'order', there 
   are (order+1)*(order+2)/2 monomials. 
   For example, if 'order' is 3, the array 'ret' has the following elements:
	ret[0]=1.0
	ret[1]=x'^1/1!
	ret[2]=y'^1/1!
	ret[3]=x'^2/2!
	ret[4]=x'^1/1!*y'^1/1!
	ret[5]=y'^2/2!
	ret[6]=x'^3/3!
	ret[7]=x'^2/2!*y'^1/1!
	ret[8]=x'^1/1!*y'^2/2!
	ret[9]=y'^3/3!,
   where x'=(x-ox)/scale and y'=(y-oy)/scale. 				     */

int	eval_2d_monoms(	double x,double y,int order,double *ret  ,  
			double ox,double oy,double scale);


/* eval_2d_poly():
   Evaluates a 2 dimensional polynomial up to order 'order' in the point (x,y) 
   using the coefficients 'coeff'. Before the evaluation, the point (x,y) is 
   shifted and rescaled using the parameters 'ox', 'oy' and 'scale'.	     */

double	eval_2d_poly  (	double x,double y,int order,double *coeff,
			double ox,double oy,double scale);

/* eval_2d_leg_monoms():
   Evaluates 2 dimensional Legendre monomials up to order 'order' in the 
   point (x,y). See also eval_2d_monoms().				     */
int	eval_2d_leg_monoms(	double x,double y,int order,double *ret,  
				double ox,double oy,double scale);
/* eval_2d_leg_poly():
   Evaluates a 2 dimensional Legendre polynomial up to order 'order' in the 
   point (x,y). See also eval_2d_poly().				     */
double	eval_2d_leg_poly(double x,double y,int order,double *coeff,
			 double ox,double oy,double scale);

/* diff_2d_poly():
   Calculates the derivative polynomial of the polynom described by 
   the coefficients 'coeff' (up to order 'order'). The calculated 
   derivatives are stored in 'dxcoeff' and 'dycoeff' if they are not NULL. 
   Note that the resulted polynomials have an order 'order'-1, therefore
   they have less coefficients.						     */
int	diff_2d_poly(double *coeff,int order,double *dxcoeff,double *dycoeff);

/* compose_2d_poly_with_affine():
   Composes the polynomial described by the 'coeff' coefficients (up to
   order 'order') with the affine transformation stored in xlin[0],[1],[2]
   and ylin[0],[1],[2]. The coefficients of the new polynomial is stored 
   in 'rc'.								     */
int	compose_2d_poly_with_affine(double *coeff,int order,
		double *xlin,double *ylin,double *rc);

/* calc_2d_unitarity():
   Calculates the unitarity of the 2d->2d transformation described by 
   the coefficients 'xfit' and 'yfit' up to the order 'order'.		     */
double	calc_2d_unitarity(double *xfit,double *yfit,int order);

#endif
                                                                   
