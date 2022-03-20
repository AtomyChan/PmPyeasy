/*****************************************************************************/
/* psf-determine.h							     */
/*****************************************************************************/

#ifndef	__PSF_DETERMINE_H_INCLUDED
#define	__PSF_DETERMINE_H_INCLUDED	1

#include <fits/fits.h>	/* typedef: fits, fitsimage	*/

#include "psf.h"	/* typedef: psf			*/
#include "stars.h"	/* typedef: candidate		*/

/*****************************************************************************/

#define		PSF_DET_NATIVE		1
#define		PSF_DET_INTEGRAL	2
#define		PSF_DET_CIRCLE		3

typedef struct
 {	int	use_biquad;
 } psfdeterminenative;

typedef struct
 {	double	kappa;
 } psfdetermineintegral;

typedef struct
 {	double	width;
	int	order;
 } psfdeterminecircle;

typedef union
 {	psfdeterminenative	native;
	psfdetermineintegral	integral;
	psfdeterminecircle	circle;
 } psfdetparam;

typedef struct
 {	int		type;
	int		hsize,grid;
        int		order;
	psfdetparam	param;
        int		is_symmetrize;	
 } psfdetermine;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int	psf_determine_native  (fitsimage *img,char **mask,
		psfcandidate *cands,int ncand,int is_subtracted,
		int hsize,int grid,int order,psf *p,int use_biquad);
int	psf_determine_integral(fitsimage *img,char **mask,
		psfcandidate *cands,int ncand,int is_subtracted,
		int hsize,int grid,int order,psf *p,double kappa);
int	psf_determine_circle  (fitsimage *img,char **mask,
		psfcandidate *cands,int ncand,int is_subtracted,
		int hsize,int grid,int order,psf *p,double circwd,int circorder);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int	psf_determine(fitsimage *img,char **mask,
		psfcandidate *cands,int ncand,int is_subtracted,
		psfdetermine *pd,psf *p);

int	drawback_psf(ipoint *ipoints,int nipoint,double *yvals,
		double x0,double y0,double amp,psf *p,double mul);

int	psf_bgamp_fit(fitsimage *img,char **mask,
		psfcandidate *cands,int ncand,int is_subtracted,psf *p);

/*****************************************************************************/

#endif

/*****************************************************************************/
