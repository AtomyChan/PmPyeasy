/*****************************************************************************/
/* common.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Some common stuff for the 'fi' package.				     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2004-2008; Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <fits/fits.h>

#include "fitsmask.h"
#include "math/fit/lmfit.h"
#include "math/spline/spline.h"
#include "statistics.h"
#include "math/poly.h"
#include "math/polyfit.h"
#include "io/tokenize.h"
#include "io/iof.h"
#include "tensor.h"

#include "fitsh.h"
#include "common.h"

/*****************************************************************************/

int mark_saturated_pixels(fitsimage *img,char **mask,fitsimage *satimg,double param,int method)
{
 int	i,j,sx,sy,c;

 if ( mask==NULL )	return(-1);

 sx=img->sx,
 sy=img->sy;

 if ( satimg != NULL )
  {	if ( sx != satimg->sx || sy != satimg->sy )	return(1);
	for ( i=0 ; i<sy ; i++ )
	 {  for ( j=0 ; j<sx ; j++ )
	     {	if ( img->data[i][j]>=satimg->data[i][j]*param )
			mask[i][j] |= MASK_OVERSATURATED;
	     }
	 }
  }
 else
  {	for ( i=0 ; i<sy ; i++ )
	 {  for ( j=0 ; j<sx ; j++ )
	     {	if ( img->data[i][j]>=param )
			mask[i][j] |= MASK_OVERSATURATED;
	     }
	 }
  }

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	c=mask[i][j];
		if ( ! ( c & MASK_OVERSATURATED ) )	continue;
		if ( method & 1 )
		 {	if ( i>0    )	mask[i-1][j] |= MASK_LEAKED;
			if ( i<sy-1 )	mask[i+1][j] |= MASK_LEAKED;
		 }
		if ( method & 2 )
		 {	if ( j>0    )	mask[i][j-1] |= MASK_LEAKED;
			if ( j<sx-1 )	mask[i][j+1] |= MASK_LEAKED;
		 }
	 }
  }
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	c=mask[i][j];
		if ( (c & MASK_OVERSATURATED) && (c & MASK_LEAKED) )
			mask[i][j] &= ~MASK_LEAKED;
	 }
  }

 return(0);
}

/*****************************************************************************/

int join_masks_from_files(char **mask,int sx,int sy,char **inmasklist)
{
 char		*inmaskname,**inmask;
 fitsheaderset	header;
 fitsimage	fi;
 int		i,j;
 FILE		*ft;

 if ( inmasklist != NULL )
  {
	for ( ; *inmasklist != NULL ; inmasklist++ )
	 {	inmaskname=*inmasklist;
		ft=fopenread(inmaskname);
		if ( ft==NULL )		return(1);
		fits_headerset_reset(&header);
		fits_headerset_read(ft,&header);
		fclose(ft);
		fits_image_get_params(&header,&fi);
		if ( sx != fi.sx || sy != fi.sy )	return(2);
		inmask=fits_mask_read_from_header(&header,sx,sy,NULL);
		fits_headerset_free(&header);
		for ( i=0 ; i<sy ; i++ )
		 {	for ( j=0 ; j<sx ; j++ )
			 {	mask[i][j] |= inmask[i][j];	}
		 }
		fits_mask_free(inmask);
	 }
  }

 return(0);
}

/*****************************************************************************/

int mark_integerlimited_pixels(fitsimage *i,char **mask,int bitpix,int is_corr,int mvlo,int mvhi)
{
 double	llo,lhi;
 int	k,l;

 if ( i==NULL || i->sx <=0 || i->sy <=0 || i->data==NULL ) /* invalid */
	return(-1);
 else if ( bitpix<0 )	/* real data, do nothing */
	return(0);
 else if ( bitpix==8 )
  {	llo=-128.0;
	lhi=+127.0;
  }
 else if ( bitpix==16 )
  {	llo=-32768.0;
	lhi=+32767.0;
  }
 else if ( bitpix==32 )
  {	llo=-2147483648.0;
	lhi=+2147483647.0;
  }
 else			/* invalid bitpix value	*/
	return(-1);

 
 for ( k=0 ; k<i->sy ; k++ )
  {	for ( l=0 ; l<i->sx ; l++ )
	 {	if ( i->data[k][l]<llo )
		 {	if ( mask != NULL )	mask[k][l] |= mvlo;
			if ( is_corr )		i->data[k][l]=llo;
		 }
		else if ( i->data[k][l]>lhi )
		 {	if ( mask != NULL )	mask[k][l] |= mvhi;
			if ( is_corr )		i->data[k][l]=lhi;
		 }
		else if ( is_corr )
			i->data[k][l]=floor(i->data[k][l]);
	 }
  }

 return(0);
}

/*****************************************************************************/
                                                                             
