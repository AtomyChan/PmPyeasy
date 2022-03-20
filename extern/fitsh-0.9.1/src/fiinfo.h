/*****************************************************************************/
/* fiinfo.h								     */
/*****************************************************************************/

#ifndef	__FIINFO_H_INCLUDED
#define	__FIINFO_H_INCLUDED	1

#include <fits/fits.h>

#include "math/point.h"

/*****************************************************************************/

typedef struct
 {	int	r,g,b;
 } color;

typedef struct
 {	color	beg,end;
 } gradient;

typedef struct
 {	int		is_color,is_invert;
	int		minmaxmethod,scalemethod;
	int		mmin_set,mmax_set;
	double		manmin,manmax;
	double		zcontrast,percentage;
	double		contrast,brightness;
	gradient	*palette;	
	int		ncol;
	int		is_flip,is_mirror;
 } pnmparam;

/* fiinfo-pnm.c */ /**********************************************************/

#define	MM_MINMAX	0	/* automatic min, automatic max	*/
#define	MM_ZSCALE	1	/* zscale			*/
#define	MM_ZMIN		2	/* zmin				*/
#define	MM_ZMAX		3	/* zmax				*/
#define	MM_PERCENTAGE	4	/* percentage (see also DS9)	*/
#define	MM_MANUAL	5	/* manual min, manual max	*/

#define	SCALE_LINEAR	0
#define	SCALE_HISTEQU	1
#define	SCALE_LOG	2
#define	SCALE_SQRT	3
#define	SCALE_SQUARED	4

int	fitsimage_dump_pnm(fitsimage *img,char **mask,FILE *fw,pnmparam *pp);
int	parse_palette(char *pstr,gradient **rpal,int *rncol);

/* fiinfo-image.c */ /********************************************************/

int	create_link_background(fitsimage *img,char **mask,fitsimage *bgi,int nx,int ny);
int	fits_stat_basic(fitsimage *img,char **mask,double *rmin,double *rmax,double *rmean,double *rstdd);

double	estimate_skysigma_naive(double *rawdata,int k,double sky,double smm,double spp);
int	fits_stat_raw_sky_skysigma(double *data,int ndat,double stddev,double *rsky,double *rskysigma);

int	fits_stat_sky(fitsimage *img,double stddev,double *rsky,double *rskysigma);
double	fits_stat_median(fitsimage *img);

int	fits_stat_background(fitsimage *img,int nbx,int nby,point **pdsky,point **pdsigma,double stddev);

int	fits_stat_sky_like_isis(fitsimage *img,double stddev,double *rsky,double *rskysigma);
int	fits_stat_sky_biquad(fitsimage *img,double stddev,double *rsky,double *rskysigma);


/*****************************************************************************/

#endif

/*****************************************************************************/
       
