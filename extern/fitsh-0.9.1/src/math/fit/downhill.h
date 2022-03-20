/*****************************************************************************/
/* downhill.h								     */
/*****************************************************************************/

#ifndef	__DOWNHILL_H_INCLUDED
#define	__DOWNHILL_H_INCLUDED	1

/*****************************************************************************/

/* downhill_simplex():
   This function minimizes the simplex p[nvar+1][ndim] using the function
   funct(). This funct() function has two arguments: an arbitrary opaque 
   parameter 'param' which is passed to funct() directly on each call
   and a vector of real numbers (doubles) with the dimension of 'dim'.
   The number of the points in the simplex is 'nvar+1', thus the minimalization
   is done on a hypersurface with the dimension of 'nvar' since the 
   consequent points of the simplex propagation are always an affine 
   combination of the previous points. 
   The function returns the number of steps which is required to reach the
   minimum within the relative precision of 'precision' unless more steps
   are required than 'maxstep'. In this case the function returns a negative
   value indicating non-convergence.					     */

int	downhill_simplex(double **p,int ndim,int nvar,
		double (*funct)(void *param,double *arr),void *param,
		double precision,int maxstep);

/*****************************************************************************/

#endif

/*****************************************************************************/
        
