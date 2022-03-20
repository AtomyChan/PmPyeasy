/*****************************************************************************/
/* biquad.h								     */
/*****************************************************************************/

#ifndef	__BIQUAD_H_INCLUDED
#define	__BIQUAD_H_INCLUDED	1

/* biquad_coeff():
   Calculates the biquadratic spline coefficients of the two-dimensional
   array 'y'. The array has a size of 'sx' times 'sy'. The calculated
   coefficients are stored in the array 'c', which should have a size
   of 2*'sx'+1 times 2*'sy'+1. The array 'c' store the whole information
   required by the biquadratic spline interpolation, so, after the 
   call, the array 'y' might be released if the tabbed grid data stored in
   'y' are not needed anymore. The optional argument 'mask' indicates 
   which points of the grid data 'y' should be used in the interpolation.
   If mask[i][j] is zero(!), the point y[i][j] is used, otherwise
   it is skipped. If 'mask' is NULL, all points are used.	             */
int	biquad_coeff(double **y,int sx,int sy,double **c,char **mask);

/* biquad_smooth():
   Smooths (makes C^1) the biquadratic surface described by the coefficients
   'c' and having a size of 'sx' by 'sy' using the mask 'mask' (if it is
   non-NULL). The function always returns successfully (with 0).	     */
int	biquad_smooth(double **c,int sx,int sy,char **mask);

/* biquad_eval():
   Evaluates the biquadratic surface stored in 'c' (returned by biquad_coeff())
   at the point (x,y). 							     */
double	biquad_eval(double **c,double x,double y);

/* biquad_quad_fit_coeff():
   Calculates the coefficients of a quadratic surface fitted to the pixel
   point (x,y) of the biquadratic surface stored in 'c'. The fitted 
   coefficients (6 numbers) are stored in 'coeff'.			     */
int	biquad_quad_fit_coeff(double **c,int x,int y,double *coeff);

/* biquad_quad_fit_eval():
   Evaluates the fitted quadratic surface at the point (x,y).		     */
double	biquad_quad_fit_eval(double **c,double x,double y);


int	biquad_quad_fit_minmax(double **c,int x,int y,double *mmr,
			       double *rx,double *ry);

/* biquad_int_pixel():
   Integrates the biquadratic surface described by 'c' on the pixel 
   [x:x+1,y:y+1]. According to the definition and the storage methods of these 
   surfaces, the return value is c[2*y+1][2*x+1]...			     */
double	biquad_int_pixel(double **c,int x,int y);

/* biquad_int_square_pixel():
   Integrates the point-by-point squared function of the biquadratic surface
   desribed by the coefficients 'c' on the pixel [x:x+1,y:y+1].		     */
double	biquad_int_square_pixel(double **c,int x,int y);

/* biquad_scatter():
   Calculates the scatter of the biquadratic surface on the pixel 
   [x:x+1,y:y+1]. This function calls the previous two functions and derives
   the scatter (standard deviation) via the usual way.			     */
double	biquad_scatter(double **c,int x,int y);

/* biquad_diff(), biquad_diff_x(), biquad_diff_y():
   Calculates the biquad-coefficients of the partial derivative (by x
   or y) surface of the surface described by the coefficients 'c' (with
   the size of 'sx' and 'sy'). The coefficients of the new surface are
   going to be stored in the array 'd' (which like 'c', has to have 
   an allocated size of (2*sx+1) times (2*sy+1)). 'd' can be the same 
   as 'c', in this case the original surface is replaced by one of it's 
   (the desired) partial derivative. The optional argument 'mask' can be 
   used as a badness-mask (see biquad_coeff()), it's highly recommended to 
   use the same mask which was used for calculating the coefficients 'c'.
   If 'var' is zero, the partial derivative by x, otherwise the derivative 
   by is calculated. 
   Note that biquadratic surfaces are can only be differentiated once by
   a given variable. Otherwise this call yields unexpected results...	     */ 
int	biquad_diff(double **c,int sx,int sy,double **d,char **mask,int var);

#define	biquad_diff_x(c,sx,sy,d,mask)	biquad_diff((c),(sx),(sy),(d),(mask),0)
#define	biquad_diff_y(c,sx,sy,d,mask)	biquad_diff((c),(sx),(sy),(d),(mask),1)

#endif
                  
