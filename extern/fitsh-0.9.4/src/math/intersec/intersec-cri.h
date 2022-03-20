/*****************************************************************************/
/* intersec-cri.h 							     */
/*****************************************************************************/

#ifndef	__INTERSEC_CRI_H_INCLUDED
#define	__INTERSEC_CRI_H_INCLUDED	1

/* intersec_cri_integrate_monoms():
   Integrates the monomials 1, x, y, x^2, xy and y^2 (depending on the 
   polynomial order 'order' which can be 0, 1 or 2) on the intersection of
   the rectangle [x0,x0+sx] x [y0,y0+sy] and the circle which is centered to 
   the origin and has a radius of 'cr'. The returned integrals are stored
   in the array 'coeff' (1, 3 or 6 numbers for order=0, 1 or 2, rescpectively).
   All of these integrals are zero if the circle and the rectangle are 
   disjoint. The function returns 0 on success, otherwise it returns 
   a non-zero value (the only cases are when 'cr' is non-positive or 'order' is 
   not from the set {0,1,2}).						     */
int	intersec_cri_integrate_monoms(double x0,double y0,double sx,double sy,
	double cr,double *coeff,int order);

#endif
                                 
