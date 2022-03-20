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
#define		COM_MODE_AVG		0	/* average, simple mean	     */
#define		COM_MODE_MED		1	/* median		     */
#define		COM_MODE_REJ_AVG	2	/* sigma rejection mean	     */
#define		COM_MODE_REJ_MED	3	/* sigma rejection median    */

/* some other combinations for other purposes: */
#define		COM_MODE_SUM		4	/* sum of the images/pixels  */
#define		COM_MODE_SQSUM		5	/* squared sum		     */
#define		COM_MODE_SCT		6	/* standard deviation	     */

#define		COM_MODE_MIN		8	/* minimum		     */
#define		COM_MODE_MAX		9	/* maximum		     */

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
	int		niter;			/* used by COM_MODE_REJ_*    */
	double		lower,upper;		/* used by COM_MODE_REJ_*    */
	int		ignore_flag;		/* see COM_IGNORE_*	     */
	int		logicalmethod;		/* (0): or (!0): and 	     */
 } compar;

typedef struct 
 {	fits		*img;
	FILE		*fr;
 } comimg;

/*****************************************************************************/

/* Combination of images: user interface specific stuff (gonna be removed?!) */

int	combine_parse_mode(char *modstr,compar *cp);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* Combination of images: lower and higher level combination specific stuff  */

double	combine_points(double *points,int n,compar *cp);

int	combine_lines(double **lines,int n,int sx,double *out,
		compar *cp,char **wmask,char *outmask);

int	combine_images_from_files(comimg *inputs,int n,fits *outimg,
		compar *cp,char **inmask,char **outmask,
		presubdata *pss,int nps,size_t maxmem);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int	combine_cleanup(comimg *inputs,int ninput);

/*****************************************************************************/

#endif
                                                              
/*****************************************************************************/

