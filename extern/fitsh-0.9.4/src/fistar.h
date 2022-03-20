/*****************************************************************************/
/* fistar.h								     */
/*****************************************************************************/

#ifndef	__FISTAR_H_INCLUDED
#define	__FISTAR_H_INCLUDED	1

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define		MARK_SYM_DOTS		0
#define		MARK_SYM_SQUARES	1
#define		MARK_SYM_CIRCLES	2

typedef struct
 {	int	id,x,y;
	int	s,d,k;
	int	bg,flux,mag;
 } colinfo;

typedef struct
 {	double	x,y;
	double	s,d,k;
	double	bg,flux;
 } initcand;

typedef struct
 {	char	*name;
	double	x,y;
 } iposition;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define	I_EMPTY		0	/* empty, "-" 				     */

#define	I_ID		1	/* identifier				     */

#define	I_IX		2	/* integer part of the centroid x coordinate */
#define	I_IY		3	/* integer part of the centroid y coordinate */

/* candidate parameters */
#define	I_CX		4	/* candidate centroid x coordinate	     */
#define	I_CY		5	/* candidate centroid y coordinate	     */
#define	I_CBG		6	/* candidate background			     */
#define	I_CAMP		7	/* candidate amplitude			     */
#define	I_CMAX		8	/* candidate peak value			     */
#define	I_CPIX		9	/* number of pixels belonging to candidate   */
#define	I_CDEVS		12	/* candidate deviation momentum 'S'	     */
#define	I_CDEVD		13	/* candidate deviation momentum 'D'	     */
#define	I_CDEVK		14	/* candidate deviation momentum 'K'	     */
#define	I_CFLUX		15	/* candidate flux			     */

/* fitted model parameters */
#define	I_X		16	/* centroid x coordinate		     */
#define	I_Y		17	/* centroid x coordinate		     */
#define	I_BG		18	/* background				     */
#define	I_AMP		19	/* amplitude				     */
#define	I_DEVS		20	/* deviation momentum 'S'		     */
#define	I_DEVD		21	/* deviation momentum 'D'		     */
#define	I_DEVK		22	/* deviation momentum 'K'		     */
#define	I_DEVL		23	/* not used				     */

#define	I_MOM		24	/* polynomial deviation momenta		     */
#define	I_DSIG		25	/* covariance momentum 'sigma'		     */
#define	I_DDEL		26	/* covariance momentum 'delta'		     */
#define	I_DKAP		27	/* covariance momentum 'kappa'		     */

#define	I_FWHM		28	/* full width at half magnitude		     */
#define	I_ELLIP		29	/* ellipticity				     */
#define	I_PA		30	/* position angle			     */

#define	I_FLUX		32	/* flux					     */
#define	I_MAG		33	/* magnitude (see --mag-flux)		     */
#define	I_NOISE		34	/* something which could be called "noise".  */
#define	I_SN		35	/* signal-to-noise ratio, whatever is "noise"*/

/* fitted PSF parameters */
#define	I_PX		36	/* PSF centroid x coordinate		     */
#define	I_PY		37	/* PSF centroid y coordinate		     */
#define	I_PBG		38	/* PSF background 			     */
#define	I_PAMP		39	/* PSF amplitude, equivalent to flux	     */
#define	I_PMOMS		40	/* PSF distortion momentum 's', 1.0 exported */
#define	I_PMOMD		41	/* PSF distortion momentum 'd', 0.0 exported */
#define	I_PMOMK		42	/* PSF distortion momentum 'k', 0.0 exported */
#define	I_PMOML		43	/* PSF distortion momentum 'l', 0.0 exported */

int *	output_format_stars(char *iformat);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int	read_star_candidates(FILE *fr,colinfo *col,
	initcand **ricands,int *rnicand,magflux *mf);

int	read_input_position_list(FILE *fr,colinfo *col,iposition **,int *);

int	draw_mark_star_pixels(fitsimage *img,star *stars,int nstar);
int	draw_mark_stars(fitsimage *img,star *stars,int nstar,int msym,int msize);

int	write_stars(FILE *fw,star *stars,int nstar,
	int *indx,int is_comment,int *oformat,magflux *mf);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#endif

/*****************************************************************************/

