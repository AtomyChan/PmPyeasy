/*****************************************************************************/
/* apphot.h								     */
/*****************************************************************************/

#ifndef __FITSH_APPHOT_H_INCLUDED
#define	__FITSH_APPHOT_H_INCLUDED	1

typedef struct
 {	double	area;
	double	flux;
 } fracpixel;

#define		BGTYPE_MEAN		1
#define		BGTYPE_MEDIAN		2
#define		BGTYPE_MODE		3

#define		BGSCATTER_STDDEV	0
#define		BGSCATTER_MAD		1

typedef struct
 {	int	type;
	int	scatter;
	int	rejniter;
	double	rejlower,rejupper;
 } bgmode;

typedef struct
 {	double	ac_radius;
 } aperture_circular;

typedef struct
 {	double	x;
	double	y;
 } aperture_vertex;

typedef struct
 {	int		ap_npoint;
	aperture_vertex	*ap_vertices;
	int		ap_nvertex;
 } aperture_polygon;

typedef struct
 {	double	r0,ra,da;
	double	gain;
	bgmode	bgm;
	double	*r0_poly,*ra_poly,*da_poly;
	int	nr0,nra,nda;
 } apphotpar;

typedef struct
 {	double	fw;			/* zeroth order moment: flux 	     */
	double	fwx,fwy;		/* first order moments: centroid     */
	double	fwxx,fwxy,fwyy;		/* second order moments		     */
 } apphot_out;

int	fracpixel_order(fracpixel *fp,int n);
double	fracpixel_median(fracpixel *fp,int n);
int	fracpixel_stat(fracpixel *fp,int n,double *rs,double *rmean,double *rsigma);
int	fracpixel_reject(fracpixel *fp,int n,double lower,double higher);

int	aperture_photometry_back_int(double **data,char **mask,
		int sx,int sy,double cx,double cy,apphotpar *par,
		double *rbgarea,double *rbgflux,double *rbgavg,double *rbgsigma);

int	aperture_photometry_back_ring(double **data,char **mask,
		int sx,int sy,double cx,double cy,apphotpar *par,
		double *rbgarea,double *rbgflux,double *rbgavg,double *rbgsigma,
		int *ratot,int *rabad,char **ringmask);

int	aperture_photometry_back_polygons(double **data,char **mask,
		int sx,int sy,double cx,double cy,apphotpar *par,
		double *rbgarea,double *rbgflux,double *rbgavg,double *rbgsigma,
		int *ratot,int *rabad,char **ringmask);

int	aperture_photometry_flux_circle(double **data,char **mask,int sx,int sy,
		double cx,double cy,double r0,
		double *rarea,double *rflux,
		apphot_out *out,double background,
		int *rrtot,int *rrbad,int *rrign,
		int maskignore,double **subpixeldata,int subg,char **xmask);

int	aperture_photometry_flux_polygon(double **data,char **mask,int sx,int sy,
	        double cx,double cy,double *poly,int nppoly,
		double *rarea,double *rflux,
        	apphot_out *out,double background,
	        int *rrtot,int *rrbad,int *rrign,
	        int maskignore,char **xmask);


int	aperture_photometry(fitsimage *img,char **mask,
		double cx,double cy,
		apphotpar *par,double *rflux,double *rfluxerr);

int	aperture_photometry_flux_biquad(double **c,char **mask,int sx,int sy,
		double cx,double cy,double r0,
		double *rarea,double *rflux,
		int *rrtot,int *rrbad,int *rrign,
		int maskingonre,char **xmask);

#endif
