/*****************************************************************************/
/* bicubic.c 								     */
/*****************************************************************************/

#ifndef	__BICUBIC_H_INCLUDED	
#define	__BICUBIC_H_INCLUDED	1

/* bicubic_coeff():
   Calculates the bicubic spline coefficients of the surface tabulated at
   the array 'y' (which has a size of 'sx' times 'sy'). The coefficients are
   stored at 'c', which should have a size of at least 2*'sx' times 2*'sy'.
   All information required by the bicubic interpolation is stored in 'c',
   therefore the array 'y' can be released after this call, if it is 
   not needed anymore. The optional argument 'mask' can be used for marking
   points of 'y' to be unused ones. If the element mask[j][i] is zero(!),
   the point y[j][i] will be used during the interpolation, 
   otherwise it won't be used. If 'mask' is set to be NULL, all points
   of the array 'y' are used.						     */
int	bicubic_coeff(double **y,int sx,int sy,double **c,char **mask);

/* bicubic_inter():
   Calculates the interpolation value of the bicubic surface stored in 'c' 
   at the point ('x','y'). The bicubic coefficients can be calculated with
   the function bicubic_coeff(). Note that if the optional argument 'mask'
   was used during the call of this function and the element mask[j][i]
   was marked as a point which mustn't be used, the arguments x,y in the 
   range [j:j+1,i:i+1] are not recommended to use and can cause unexpected
   results (NaNs).							     */
double	bicubic_inter(double **c,double x,double y);

#endif
 
