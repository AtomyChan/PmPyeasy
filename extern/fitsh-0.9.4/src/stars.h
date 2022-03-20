/*****************************************************************************/
/* stars.h								     */
/*****************************************************************************/

#ifndef	__STARS_H_INCLUDED
#define	__STARS_H_INCLUDED	1

#include <fits/fits.h>
#include "background.h"

typedef struct
 {	int	xmin,xmax;
	int	ymin,ymax;
 } range;

#ifndef	__IPOINT_STRUCT_DEFINED
#define	__IPOINT_STRUCT_DEFINED
typedef struct
 {	int	x,y;
 } ipoint;
#endif

typedef struct
 {	ipoint	*ipoints;
	int	nipoint;
 } ipointlist;

typedef struct
 {	ipoint	*ipoints;
	double	*yvals;
	int	nipoint;
 } mpointlist;

typedef struct
 {	int	ix,iy;			/* Integer coordinates of a star can-*/
					/* didate, the image has a local max-*/
					/* imum in the point (ix,iy), or     */
					/* something like that...	     */

	double	cx,cy;			/* A better approximation for the    */
					/* center of a star candidate, deri- */
					/* ved by a second-order polynomial  */
					/* fit (or something else, see also  */
					/* refine_candidate_params()...	     */

	double	peak,amp,bg;			
	double	sxx,syy,sxy;		/* Other parameters derived from this*/
					/* second-order polynomial fit.	     */
					/* All of these parameters are set   */
					/* by search_star_candidates() and   */
					/* some of them are used by 	     */
					/* markout_stars() and 		     */
					/* fit_gaussians().		     */

	ipoint	*ipoints;		/* ipoints[] contains nipoint	     */
	int	nipoint;		/* elements, (x,y) pairs, which	     */
					/* possible belong to the star.      */
					/* Determined by markout_*() and     */
					/* used by fit_gaussians().	     */

	double	area,			/* area (it is always an integer)    */
		noise,			/* bg noise			     */
		flux;			/* flux				     */

	int	flags,marked;

 } candidate;

#include "psf.h"

#define		SHAPE_GAUSS		1
#define		SHAPE_ELLIPTIC		2
#define		SHAPE_DEVIATED		3
#define		SHAPE_PSF		4

#define		MAX_DEVIATION_ORDER	4
#define		MAX_DEVIATION_COEFF	15	/* (MDO+1)*(MDO+2)/2	*/

/* #define	STAR_MULTIMODEL		1 */	/* obsoleted		*/

typedef struct
 {	double	gamp,gbg;
	double	gcx,gcy;
 } starlocation;

typedef struct
 {	int	model,order;	/* `model` can be: SHAPE_GAUSS, ELLIPTIC,    */
				/* DEVIATED. `order` is only defined	     */
				/* for SHAPE_DEVIATED (act. between 2 and 4) */
				/* PSF fit results are stored in another     */
				/* form (see (starlocation)star->psf).	     */

	double	gs,gd,gk,gl;	/* coeff's for SHAPE_ELLIPTIC (gs, gd, gk),  */
				/* DEVIATED (gs) and PSF (gs, gd, gk and gl).*/

	double	mom[MAX_DEVIATION_COEFF];	/* coeff's for SHAPE_DEVIATED*/

	double	factor;		/* factor of multi-model fitting (should not */
				/* really expand out from the interval [0,1])*/

 } starshape;

typedef struct
 {	starlocation	location;	/* An approximation for the center,  */
					/* background and amplitude of the   */
					/* star	(`gamp` is model-dependent). */

	starshape	shape;		/* Shape parameters of the star.     */

	starlocation	psf;		/* PSF fitting yields only info's    */
					/* like this: centroid coordinates,  */
					/* background and amplitude (later is*/
					/* equivalent to the total flux...). */

	double		gsig,gdel,gkap;	/* Derived parameters from the Gau-  */
	double		gfwhm,		/* sian fit (FWHM, ellpticity and    */
			gellip,		/* position angle), also set by      */
			gpa;		/* fit_gaussians().		     */

	double		flux;		/* The total flux of the star.	     */
					/* Set by fit functions, and used    */
					/* by firandom also.		     */

	int		marked;		/* A flag for cleanup_starlist().    */

	candidate	*cand;		/* candidate pointer (if available). */
 } star;

/* temporary removed from typedef struct { ... } star; */
#ifdef	STAR_MULTIMODEL
	starshape	*mshapes;	/* Shape parameters of the star if   */
	int		nmshape;	/* multi-model fitting was used.     */
#endif

typedef struct
 {	int	ix,iy;
	double	value;
 } imgpoint;

typedef struct
 {	int	model;	/* can be SHAPE_{GAUSS,ELLIPTIC,DEVIATED}	     */
	int	order;	/* up to MAX_DEVIATION_ORDER, only for SHAPE_DEVIATED*/
/*	double	igs;*/	/* initial value of the 'gs' starshape parameter     */
 } starmodelfit;

/* star-base.c */ /***********************************************************/

/* fit_small_parabola_{point,block,block_param}():
   Some common parabola fitting routines used during candidate searching
   (see the appropriate star-cand-*.c modules).				     */
int	fit_small_parabola_point(fitsimage *img,int x,int y,double *fit);
int	fit_small_parabola_block(fitsimage *img,int x,int y,double *fit);
int	fit_small_parabola_block_param(fitsimage *img,int j,int i,double *rcx,double *rcy,double *raxx,double *raxy,double *rayy,double *rpeak);

/* order_candidates_by_peak():
   Orders star candidate array cands[] (with 'ncand' elements) by peak.	     */
int	order_candidates_by_peak(candidate *cands,int ncand);

/* refine_candidate_params():
   Refines the candidate parameters - the centroid coordinates and the 
   shape parameters - of the candidates in the array cands[]. This refine
   method is based on the calculation of the central momenta (up to order of 2)
   of the points belonging to the given star candidate.			     */
int	refine_candidate_params(fitsimage *img,candidate *cands,int ncand);

int	star_set_common_shape_params(double gs,double gd,double gk,star *ws);
double	star_get_unity_flux(starshape *sw);

/* convert_candidates():
   Convert candidate information into star information. The stars are modelled
   by an elliptical Gaussian function (so, star->model is SHAPE_ELLIPTIC),
   the parameters (gs, gd, gk, amplitude, ...) are calculated from the
   approximations sxx, sxy, syy (see the definition of the 'star' and 
   'candidate' structures).						     */
int	convert_candidates(candidate *cands,int ncand,star **rstars,int *rnstar);

/* cleanup_candlist(), cleanup_starlist():
   Remove marked candidates/stars from the arraies pointed by rcands and rstars
   respectively. 							     */
int	cleanup_candlist(candidate **rcands,int *rncand);
int	cleanup_starlist(star **rstars,int *rnstar);

/* free_stars(), free_candidates(): 
   Releases fully the array 'stars' and 'cands' (thus all internal dynamical
   allocations are also freed). 					     */
int	free_stars(star *stars,int nstar);
int	free_candidates(candidate *cands,int ncand);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

typedef struct
 {	void	(*funct)(void *,double *,double *,double *,void *);
	int	nshape;
	void	*param;
 } modelparam;

/* model_merge(), model_combine(): 
   Combination of star-models (including PSF). The parameter 'p' has a 
   type of 'modelparam[]', an array of a set of models described 
   by 'funct' and 'nshape' number of shape parameters.		             */
void	model_merge(void *xpnt,double *a,double *yy,double *dyda,void *p);
void	model_combine(void *xpnt,double *a,double *yy,double *dyda,void *p);

/* star-model.c */ /**********************************************************/
void	gauss_1d_prof(void *xpnt,double *a,double *yy,double *dyda);

void	gauss_2d_nmom(void *xpnt,double *a,double *yy,double *dyda);
void	gauss_2d_wmom(void *xpnt,double *a,double *yy,double *dyda);
void	gauss_2d_nmom_subps(void *xpnt,double *a,double *yy,double *dyda);
void	gauss_2d_wmom_subps(void *xpnt,double *a,double *yy,double *dyda);

#define	FIT_AB			0x01
#define	FIT_XY			0x02
#define	FIT_WIDTH		0x04
#define	FIT_DEVIATION		0x08

typedef struct
 {	int	fit_flags;
	int	iter_symmetric;
	int	iter_general;
 } starfit;

int	fit_star_single_model(fitsimage *img,char **mask,candidate *cands,int ncand,
		star **rstars,int *rnstar,starfit *fgp,int model,int modelorder);

#ifdef	STAR_MULTIMODEL
int	fit_star_multi_models(fitsimage *img,char **mask,candidate *cands,int ncand,
		star **rstars,int *rnstar,starfit *fgp,starmodelfit *smfs,int nsmf);
int	fit_star_model(fitsimage *img,char **mask,candidate *cands,int ncand,
		star **rstars,int *rnstar,starfit *fgp,int model,int modelorder);
#endif

int	collective_fit_star_single_model(fitsimage *img,char **mask,
		star *stars,int nstar,ipointlist *ipl,
		starfit *sfp,int is_putback,int level);
int	collective_fit_star_single_model_iterative(fitsimage *img,char **mask,
		star *stars,int nstar,ipointlist *ipl,
		starfit *sfp,int level,int niter);
int	collective_fit_star_single_model_blocked(fitsimage *img,char **mask,
	        star *stars,int nstar,ipointlist *ipl,double bhsize);

/* star-psf.c */ /************************************************************/
typedef struct
 {	int	iterations;
 } psffit;

int	fit_star_psf_native(fitsimage *img,char **mask,candidate *cands,int ncand,
		star **rstars,int *rnstar,psffit *pfp,psf *tpd);

int	fit_star_psf(fitsimage *img,char **mask,candidate *cands,int ncand,
		star **rstars,int *rnstar,psffit *pfp,psf *tpd);

/* star-cand-pp.c */ /********************************************************/
int	search_star_candidates(fitsimage *img,char **mask,
	candidate **rcands,int *rncand,range *srcrange,
	double threshold,spatial *bg,double skysigma);

int	markout_candidates(fitsimage *img,char **mask,candidate *cands,int ncand);

/* star-cand-biq.c */ /*******************************************************/

int	search_star_candidates_biquad(fitsimage *img,char **mask,
	candidate **rcands,int *rncand,range *srcrange);

/* star-cand-trb.c */ /*******************************************************/
int	search_star_candidates_trb(fitsimage *img,char **mask,
	candidate **rcands,int *rncand,range *srcrange,
	double treshold);

/* star-cand-lnk.c */ /*******************************************************/
int	search_star_candidates_link(fitsimage *img,char **mask,
	candidate **rcands,int *rncand,range *srcrange,
	double threshold,double fluxthreshold,double critical_prominence);

/* star-draw.c */ /***********************************************************/

/* star_draw_gauss(), star_draw_deviated(), star_draw_psf():
   Draws analytic stars with the profile models Gaussian, deviated
   or PSF-based. The star will be written to the array 'iarr' with the
   size of 'sx' time 'sy'. The center of the star profile is (x0,y0). 
   The additional parameters are model-dependent (is, id, ik: covariance
   parameters for Gaussian profile; gs, order, mom: central coefficient,
   order and moments for deviated profiles and the PSF 'p' itself with
   extra parameters /px, py, is, id, ik and il/).			     */
int star_draw_gauss
	(double **iarr,int sx,int sy,
	 double x0,double y0,double is,double id,double ik);

int star_draw_deviated
	(double **iarr,int sx,int sy,
	 double x0,double y0,double gs,int order,double *mom);

int star_draw_psf
	(double **iarr,int sx,int sy,
	 double x0,double y0,psf *p,double px,double py,
	 double is,double id,double ik,double il);

/*****************************************************************************/

int drawback_model(ipoint *ipoints,int nipoint,double *yvals,
	starlocation *loc,starshape *shp,double mul);

/*****************************************************************************/

#endif

/*****************************************************************************/
