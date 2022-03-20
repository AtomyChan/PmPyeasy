/*****************************************************************************/
/* common.h								     */
/*****************************************************************************/

#ifndef	__FITSH_COMMON_H_INCLUDED
#define	__FITSH_COMMON_H_INCLUDED

/*****************************************************************************/

#include <fits/fits.h>
#include "fitsmask.h"

#define		OOSQ2PI		0.3989422804014326779
#define		OOSQ2		0.7071067811865475244

/*****************************************************************************/

typedef struct 
 {	double	x,y;
 	double	dx,dy;
 } dpoint;

/*****************************************************************************/

int	mark_saturated_pixels(fitsimage *img,char **mask,fitsimage *satimg,double param,int method);

int	join_masks_from_files(char **mask,int sx,int sy,char **inmasklist);

int	mark_integerlimited_pixels(fitsimage *img,char **mask,int bitpix,int is_corr,int mvlo,int mvhi);

/*****************************************************************************/

#endif

/*****************************************************************************/
                                                                   
           
