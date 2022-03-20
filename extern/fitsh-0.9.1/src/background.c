/*****************************************************************************/
/* background.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Functions related to background determination...			     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <fits/fits.h>

#include "fitsmask.h"
#include "math/fit/lmfit.h"
#include "statistics.h"
#include "math/poly.h"
#include "math/polyfit.h"
#include "tensor.h"

#include "background.h"

/*****************************************************************************/

int determine_background(fitsimage *img,spatial *bg,int gx,int gy,int order)
{
 int	i,j,sx,sy,nvar;
 point	*dp;

 if ( img==NULL )	return(1);
 if ( img->data==NULL )	return(1);

 sx=img->sx,sy=img->sy;

 if ( gx<=0 || gy<=0 )	gx=gy=order+1;

 nvar=(order+1)*(order+2)/2;
 bg->order=order;
 bg->ox=0.5*(double)sx,
 bg->oy=0.5*(double)sy,
 bg->scale=0.5*(double)sx;
 bg->coeff=(double *)malloc(sizeof(double)*nvar);

 dp=(point *)malloc(sizeof(point)*gx*gy);
 
 for ( i=0 ; i<gy ; i++ )
  { for ( j=0 ; j<gx ; j++ )
     {	int	imin,imax,jmin,jmax,k;
	int	ii,jj,n;
	double	*rawdata,med;
	imin=i*sy/gy,imax=(i+1)*sy/gy;
	jmin=j*sx/gx,jmax=(j+1)*sx/gx;
	rawdata=(double *)malloc(sizeof(double)*(imax-imin)*(jmax-jmin));
	for ( ii=imin,n=0 ; ii<imax ; ii++ )
	 {  for ( jj=jmin ; jj<jmax ; jj++ )
	     {	if ( img->data[ii][jj]>0.0 )
			rawdata[n]=img->data[ii][jj],n++;
	     }
	 }
	med=median(rawdata,n);
	free(rawdata);
	k=i*gx+j;
	dp[k].x=(double)(jmax+jmin-1)*0.5;
	dp[k].y=(double)(imax+imin-1)*0.5;
	dp[k].value=med;
	dp[k].weight=1.0;
     }
  }	
 fit_2d_poly(dp,gx*gy,order,bg->coeff,bg->ox,bg->oy,bg->scale);
 
 free(dp);

 return(0);
}

/*****************************************************************************/
                                                                        
