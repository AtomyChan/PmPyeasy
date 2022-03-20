/*****************************************************************************/
/* combine.h								     */
/*****************************************************************************/

#ifndef	__COMBINE_H_INCLUDED
#define	__COMBINE_H_INCLUDED	1

/*****************************************************************************/

#include <fits/fits.h>

#include "fitsmask.h"

/*****************************************************************************/

/* classic averaging modes (REJ_AVG and REJ_MED uses `niter` and `[lh]sigma`)*/
#define		COM_MODE_MEAN		0	/* mean			     */
#define		COM_MODE_MEDIAN		1	/* median		     */
#define		COM_MODE_REJ_MEAN	2	/* sigma rejection w/ mean   */
#define		COM_MODE_REJ_MEDIAN	3	/* sigma rejection w/ median */
#define		COM_MODE_TRUNC_MEAN	4	/* truncated mean	     */
#define		COM_MODE_WINS_MEAN	5	/* winsorized mean	     */

/* some other combinations for other purposes: */
#define		COM_MODE_SUM		16	/* sum of the images/pixels  */
#define		COM_MODE_SQSUM		17	/* squared sum		     */
#define		COM_MODE_SCT		32	/* standard deviation	     */

#define		COM_MODE_MIN		128	/* minimum		     */
#define		COM_MODE_MAX		129	/* maximum		     */

#define		COM_MODE_REJ_DEPRECATED	255	/* old stuff		     */

/* ignore modes: */
#define		COM_IGNORE_NEGATIVE	0x01
#define		COM_IGNORE_ZERO		0x02
#define		COM_IGNORE_POSITIVE	0x04

typedef struct
 {	fitsimage	*img;
	double		scale;
	int		x0,y0;			
 } presubdata;

typedef struct
 {	int		mode;			/* see COM_MODE_* definition */
	int		niter;			/* see COM_MODE_REJ_*        */
	double		lower,upper;		/* see COM_MODE_REJ_*        */
	int		lowest,highest;		/* see COM_MODE_MEAN_LH_REJ  */
	int		ignore_flag;		/* see COM_IGNORE_*	     */
	int		logicalmethod;		/* (0): or (!0): and 	     */
 } combine_parameters;

typedef struct 
 {	fits		*img;
	FILE		*fr;
 } comimg;

/*****************************************************************************/

/* Combination of images: user interface specific stuff: */

int	combine_parameters_reset(combine_parameters *cp);
int	combine_parse_mode(char *modstr,combine_parameters *cp);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* Combination of images: lower and higher level combination specific stuff  */

double	combine_points(double *points,int n,combine_parameters *cp);

int	combine_lines(double **lines,int n,int sx,double *out,
		combine_parameters *cp,char **wmask,char *outmask);

int	combine_images_from_files(comimg *inputs,int n,fits *outimg,
		combine_parameters *cp,char **inmask,char **outmask,
		presubdata *pss,int nps,size_t maxmem);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int	combine_cleanup(comimg *inputs,int ninput);

/*****************************************************************************/

#endif
                                                              
/*****************************************************************************/

