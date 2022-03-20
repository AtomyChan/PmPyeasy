/*****************************************************************************/
/* kernel.h								     */
/*****************************************************************************/

#ifndef	__FI_KERNEL_H_INCLUDED
#define	__FI_KERNEL_H_INCLUDED	1

/*****************************************************************************/

#define		KERNEL_UNKNOWN		0
#define		KERNEL_BACKGROUND	1
#define		KERNEL_IDENTITY		2
#define		KERNEL_GAUSSIAN		3
#define		KERNEL_DDELTA		4

/*****************************************************************************/

typedef struct 
 {	
	int	type;		/* 0: unknown (hsize, image and offset used) */
				/* 1: background			     */
				/* 2: identity				     */
				/* 3: Gaussian (hsize, image and offset are  */
				/*    used; sigma, bx, by: parameters        */
				/* 4: double-delta (hsize and image don't    */
				/*    have to be used, offset is used, para- */
				/*    meters: bx,by                          */
						
	double	sigma;		/* sigma, bx, by: the parameters of the	     */
	int	bx,by;		/* theoretical kernel (used only in real     */
				/* gaussian kernels, to describe them, not   */
				/* used during the convolution		     */
				/* if type==4, kernel is equivalent to       */
				/*   - image[hsize   ][hsize   ]=-0.5,       */
				/*   - image[hsize+by][hsize+bx]=+0.5        */

	int	hsize;		/* The half width/height of the kernel,	the  */
				/* real width/height is 2*hsize+1 pixels.    */
	double	**image;	/* Array of (2*hsize+1)*(2*hsize+1) doubles: */
				/* the kernel if type==0 or 3                */
	double	offset;		/* After convolution, the sum should be inc- */
				/* reased with this value (if type==0 or 3)  */
				/* (note that all other kernels can be desc- */
				/* ribed with image/offset also, but the co- */
				/* nvolution function convolve_point() is    */
				/* optimized for these types...              */

	int	order;		/* Order of spatial variation.		     */
	double	*coeff;		/* Polynomial coefficients of the variation, */
				/* a dynamically allocated array by 	     */
				/* fit_kernel_coefficients_*() with a size   */
				/* of (order+1)*(order+2)/2, obviously...    */

	int	flag;		/* A flag, for any kind of usage...	     */
				/* Actually, is used to flag whether the     */
				/* given kernel should be		     */
				/*  - totally ignored (<0),		     */
				/*  - taken as a fitted kernel (=0), or      */
				/*  - mark for fitting (>0). 		     */
				/* These values are determined after reading */
				/* the kernel lists (-k, -ik), and partially */
				/* exported (see the switch -ok) if desired. */

	int	target;		/* target flag for the convolution:	     */
				/* 0: reference should be convolved,	     */
				/* !0: image should be convolved.	     */

 } kernel;

typedef struct
 {	kernel	*kernels;	/* List of the kernels.		*/
	int	nkernel;	/* Number of kernels.		*/

	double	ox,oy,scale;	/* Shift & scale of the poly-	*/
				/* nomial coefficients.  	*/
	int	type;		/* The type of the kernel set: 	*/
				/*  0: not fitted		*/
				/*  !0: fitted, polynomial var. */
				/*  (coeffs stored at kernels)	*/
	
 } kernellist;

typedef struct 			/* Stamp area:			*/
 {	int	x1,y1;		/* 	x1 <= x < x2,		*/
	int	x2,y2;		/*	y1 <= y < y2;		*/
 } rectan;

/*****************************************************************************/

int	kernel_set_verbosity(int vlevel);

/*****************************************************************************/

double	hermite(int n,double x);
double	monom(int n,double x);
double	eval_gaussian(kernel *k,double u,double v);

int	kernel_image_norm(kernel *k,double sum);
int	kernel_image_calc_gaussian(kernel *k);
int	kernel_image_calc_linear(kernel *k);

int	kernel_image_subtract(kernel *k,kernel *subk);

double	convolve_point(fitsimage *img,kernel *k,int x,int y);
int	convolve_image(fitsimage *img,char **mask,kernel *k,fitsimage *out);

int	fit_kernel_poly_coefficients_block(fitsimage *refimg,fitsimage *img,
	char **mask,double **weight,int nbx,int nby,kernellist *klist,kernellist *xlist);

int	fit_kernel_poly_coefficients_stamp(fitsimage *refimg,fitsimage *img,
	char **mask,double **weight,rectan *rects,int nrect,kernellist *klist,kernellist *xlist);

int	convolve_with_kernel_set(fitsimage *refimg,
	char **mask,kernellist *klist,fitsimage *outimg);

int	convolve_to_subtracted(fitsimage *ref,fitsimage *img,char **mask,
	kernellist *klist,kernellist *xlist,fitsimage *outimg);

/*****************************************************************************/

kernel*	add_kernel_empty(kernellist *kl);
int	add_kernel_gaussian(kernellist *kl,double sigma,int hx,int hy,int hsize,int sporder);
int	add_kernel_gaussian_set(kernellist *kl,double sigma,int order,int hsize,int sporder);
int	add_kernel_background(kernellist *kl,int sporder);
int	add_kernel_identity(kernellist *kl,int sporder);
int	add_kernel_linear(kernellist *kl,int px,int py,int sporder);
int	add_kernel_linear_set(kernellist *kl,int ksize,int sporder);

int	kernel_init_images(kernellist *kl);
int	create_kernels_from_kernelarg(char *kernelarg,kernellist *kl);


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int	dump_kernel(FILE *fw,kernel *k);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int	kernel_info_dump_image(FILE *fw,kernel *k,int is_comment);
int	kernel_info_write(FILE *fw,kernellist *kl);

int	kernel_info_read(FILE *fr,kernellist *kl);

/*****************************************************************************/

#endif	
