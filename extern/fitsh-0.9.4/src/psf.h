/*****************************************************************************/
/* psf.h								     */
/*****************************************************************************/

#ifndef	__PSF_H_INCLUDED
#define	__PSF_H_INCLUDED	1

/*****************************************************************************/

/* This structure describes a point-spread function (PSF) in a _tabulated_   */
/* way, including spatial variations. The half-size of the PSF stamp is      */
/* 'hsize', each pixel has 'grid' times 'grid' equilateral subpixels while   */
/* the order of the spatial variation is 'order' (it means polynomial order).*/
/* Therefore, the whole data describing this spatially vairating tabulated   */
/* PSF requires (2*hsize+1) x (2*hsize+1) x ((order+1)*(order+2)/2) numbers  */
/* which are stored in 'coeff'. During the polynomial evaluation of the	     */
/* tabulated data, the parameters 'ox', 'oy' and 'scale' are also used	     */
/* (see the appropriate functions in ./math/poly.c, e.g. eval_2d_poly()).    */

typedef struct
 {	int	hsize;			/* FITS header:	'PSFHSIZE'	*/
	int	grid;			/*		'PSFSGRID'	*/
	int	order;			/*		'PSFORDER'	*/
	double	ox,			/*		'PSFOFFSX'	*/
		oy,			/*		'PSFOFFSY'	*/
		scale;			/*		'PSFSCALE'	*/

	double	***coeff;		/* NAXIS=3, [NX3][NX2][NX1]	*/
					/* NAXIS1=grid*(2*hsize+1)	*/
					/* NAXIS2=grid*(2*hsize+1)	*/
					/* NAXIS3=(order+1)*(order+2)/2	*/
 } psf;

/* The whole tabulated PSF can be stored as a three-dimensional FITS image   */
/* also. The paremeters 'hsize', 'grid', 'order', 'ox', 'oy' and 'scale' are */
/* stored as special mandatory keywords, while the coefficients 'coeff' are  */
/* stored in the primary image array. Altough some of the parameters can be  */
/* calculated directly from FITS header values, the presence of _all_ of     */
/* these headers are required.						     */

/*****************************************************************************/

#ifndef __IPOINT_STRUCT_DEFINED		/* this one also can be defined in   */
#define __IPOINT_STRUCT_DEFINED		/* "stars.h", however, this library  */
typedef struct				/* is independent from 'stars.a'.    */
 {	int     x;
	int	y;
 } ipoint;
#endif

typedef struct
 {	ipoint	*ipoints;
	double	*yvals;
	int	nipoint;
	double	bg,amp;
	double	cx,cy;
 } psfcandidate;

/*****************************************************************************/

#endif

/*****************************************************************************/
                                                                 
