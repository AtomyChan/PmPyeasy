/*****************************************************************************/
/* fiphot.h								     */
/*****************************************************************************/

#ifndef	__FIPHOT_H_INCLUDED
#define	__FIPHOT_H_INCLUDED	1

#include "apphot.h"		/* typedef: bgmode */

/*****************************************************************************/

#define		USE_WEIGHT_SUBTRACTED		1
#define		USE_WEIGHT_WEIGHTED		2

/*****************************************************************************/

typedef struct
 {	double	r0,
		ra,
		da;
 } apgeom;

typedef struct
 {	apgeom	ag;
	double	bgarea,bgflux,bgmedian,bgsigma;
	double	flux,fluxerr;
	double	cntr_x,cntr_x_err,
		cntr_y,cntr_y_err,
		cntr_width,cntr_w_err,cntr_w_d,cntr_w_k;
	double	mag ,magerr ;
	int	flag,rtot,rbad,rign,atot,abad;
 } photflux;

typedef struct
 {	double		x,y;

	apgeom		*inaps;
	int		ninap;

	int		n;		
	photflux	*fluxes;	/* size == n 			     */
	apgeom		optimal;

	photflux	*rfflux;	/* size == n 			     */
					/* these are read from the raw pho-  */
					/* tometry file			     */

	int		use_ref;	/* use a reference magnitude?	     */
	double		ref_mag;	/* reference magnitude for imgsub    */
	double		ref_col;	/* reference color for postmagfit    */
	double		ref_err;	/* reference magnitude error	     */

	char		*id;
 } photstar;

typedef struct
 {	int	colx,coly;
	int	colid,colmag,colcol,colerr;
	int	*colap;	
 } colread;

typedef struct
 {	bgmode	bgm;
	int	maskignore;
	int	use_biquad;
	double	wconfdist;
	double	sky;
	int	use_sky;
	int	is_disjoint_rings;
	int	is_disjoint_apertures;
	double	disjoint_radius;
	double	**subpixeldata;
	int	subg;
 } xphotpar;

typedef struct
 {	int	orders[4];	/* spatial orders		*/
	int	norder;		/* order in the color		*/
	int	niter;		/* number of iterations		*/
	double	sigma;		/* rejection limit in sigma	*/
 } magfitparams;

typedef struct
 {	int	nstar;		/* total number of stars		*/
	int	naperture;	/* total number of apertures used	*/
	int	*ninit;		/* number of initially used stars	*/
	int	*nrejs;		/* number of rejected stars		*/
 } magfitstat;

typedef struct
 {	int	order;		/* spatial order of gain variations	     */
	double	*coeff;		/* gain variation coefficients, B(order+2,2) */
	double	vmin;		/* minimal value of the gain		     */
 } spatialgain;

/*****************************************************************************/

extern	int	is_comment,is_verbose;

/*****************************************************************************/

int	read_input_star_list(FILE *fr,photstar **rps,int *rnp,colread *col,int zoom);
int	read_raw_photometry(FILE *fr,photstar **rps,int *rnp);

int	get_id_format_parameters(photstar *ps,int np,int *rid_len,char *idf,char *ndf);
int	write_output_star_list(FILE *fw,photstar *ps,int np);

int	photometry_status(char *buff,photflux *pf);
int	write_photometry(FILE *fw,photstar *ps,int np,char *ofxy,char *ofph,spatialgain *sg,int basistype,char *nanstring,char *serialstring,int sx,int sy);
int	write_raw_photometry(FILE *fw,photstar *ps,int np,int basistype,char *nanstring);

char *	output_format_xy(char *oformat);
char *	output_format_ph(char *oformat);

char *	subtracted_format_xy(char *sformat);
char *	subtracted_format_ph(char *sformat);

int *	create_col_ap_data(char *apcolpar);
int	create_input_ap_data(char *appar,apgeom **rinaps,int *rninap,int zoom);

int	check_apertures(photstar *ps,int np);

/*****************************************************************************/

#endif

/*****************************************************************************/

