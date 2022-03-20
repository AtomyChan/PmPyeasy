/*****************************************************************************/
/* imgtrans.h 								     */
/*****************************************************************************/

#ifndef	__IMGTRANS_H_INCLUDED
#define	__IMGTRANS_H_INCLUDED	1

#include <fits/fits.h>

int	laplace_of_image(fitsimage *img);
int	laplace_of_image_ign(fitsimage *img);
int	cyclic_laplace_of_image(fitsimage *img);
int	cyclic_laplace_of_image_ign(fitsimage *img);

#endif
