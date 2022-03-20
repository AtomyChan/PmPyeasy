/*****************************************************************************/
/* star-cand-biq.c							     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* An extension to stars.c, used by `fistar` to search stars or star	     */
/* candidates on images... based on the method of biquadratic surfaces.	     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2005, Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <fits/fits.h>

#include "fitsmask.h"
#include "math/spline/biquad.h"
#include "math/fit/lmfit.h"
#include "math/poly.h"
#include "io/iof.h"
#include "tensor.h"

#include "fitsh.h"
#include "common.h"
#include "stars.h"

/*****************************************************************************/

#define		NUM_ITER		6

/*****************************************************************************/

static int biquad_poly_coefficients(double **cf,double *pcf)
{
 double wa1,wa2,wb1,wb2;

 pcf[0]=cf[0][0];
 pcf[1]=-4.0*cf[0][0]+6.0*cf[0][1]-2.0*cf[0][2];
 wa1   =-4.0*cf[1][0]+6.0*cf[1][1]-2.0*cf[1][2];
 wa2   =-4.0*cf[2][0]+6.0*cf[2][1]-2.0*cf[2][2];
 pcf[2]=3.0*(cf[0][0]-2.0*cf[0][1]+cf[0][2]);
 wb1   =3.0*(cf[1][0]-2.0*cf[1][1]+cf[1][2]);
 wb2   =3.0*(cf[2][0]-2.0*cf[2][1]+cf[2][2]);
 pcf[3]=-4.0*cf[0][0]+6.0*cf[1][0]-2.0*cf[2][0];
 pcf[4]=-4.0*pcf[1]+6.0*wa1-2.0*wa2;
 pcf[5]=-4.0*pcf[2]+6.0*wb1-2.0*wb2;
 pcf[6]=3.0*(cf[0][0]-2.0*cf[1][0]+cf[2][0]);
 pcf[7]=3.0*(pcf[1]-2.0*wa1+wa2);
 pcf[8]=3.0*(pcf[2]-2.0*wb1+wb2);

 return(0);
}
static int biquad_poly_coeff_rearrange(double *pcf,double x,double y,double *out)
{
 double	ncf[9];
 int	i;

 ncf[0]=pcf[0]+x*(pcf[1]+pcf[2]*x)+y*((pcf[3]+x*(pcf[4]+pcf[5]*x))+y*
	(pcf[6]+x*(pcf[7]+pcf[8]*x)));
 ncf[1]=pcf[1]+2*pcf[2]*x+y*((pcf[4]+2*pcf[5]*x)+y*(pcf[7]+2*pcf[8]*x));
 ncf[2]=pcf[2]+y*(pcf[5]+y*pcf[8]);
 ncf[3]=pcf[3]+2*pcf[6]*y+x*((pcf[4]+2*pcf[7]*y)+x*(pcf[5]+2*pcf[8]*y));
 ncf[4]=pcf[4]+2*pcf[5]*x+2*y*(pcf[7]+2*pcf[8]*x);
 ncf[5]=pcf[5]+2*pcf[8]*y;
 ncf[6]=pcf[6]+x*(pcf[7]+x*pcf[8]);
 ncf[7]=pcf[7]+2*pcf[8]*x;
 ncf[8]=pcf[8];
 for ( i=0 ; i<9 ; i++ )
  {	out[i]=ncf[i];		}
 return(0);
}

/*****************************************************************************/

int search_star_candidates_biquad(fitsimage *img,char **mask,candidate **rcands,int *rncand,range *srcrange)
{
 double		**bqc,*bql[3],xd,xu,yl,yr,det,x,y,pcf[9],x0,y0,peak;
 int		sx,sy,i,j,k;
 candidate	*cands,*wc;
 int		ncand;

 if ( img==NULL || img->data==NULL )	return(1);
 sx=img->sx,sy=img->sy;
 if ( sx<=0 || sy<=0 )			return(1);

 bqc=(double **)tensor_alloc_2d(double,2*sx+1,2*sy+1);

 biquad_coeff(img->data,sx,sy,bqc,mask);
 cands=NULL;
 ncand=0;

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( mask != NULL && mask[i][j] )	continue;
		bql[0]=&bqc[2*i+0][2*j];
		bql[1]=&bqc[2*i+1][2*j];
		bql[2]=&bqc[2*i+2][2*j];

		k=0;
		if ( 2*bql[0][0]+  bql[0][2]<3*bql[0][1] && 
		       bql[0][0]+2*bql[0][2]<3*bql[0][1] )	k++;
		if ( 2*bql[1][0]+  bql[1][2]<3*bql[1][1] &&
		       bql[1][0]+2*bql[1][2]<3*bql[1][1] )	k++;
		if ( 2*bql[2][0]+  bql[2][2]<3*bql[2][1] &&
		       bql[2][0]+2*bql[2][2]<3*bql[2][1] )	k++;
		if ( k<3 )	continue;

		if ( 2*bql[0][0]+  bql[2][0]<3*bql[1][0] &&
		       bql[0][0]+2*bql[2][0]<3*bql[1][0] )	k++;
		if ( 2*bql[0][1]+  bql[2][1]<3*bql[1][1] &&
		       bql[0][1]+2*bql[2][1]<3*bql[1][1] )	k++;
		if ( 2*bql[0][2]+  bql[2][2]<3*bql[1][2] &&
		       bql[0][2]+2*bql[2][2]<3*bql[1][2] )	k++;
		if ( k<6 )	continue;

		xd=(2*bql[0][0]-3*bql[0][1]+bql[0][2])/(3*(bql[0][0]-2*bql[0][1]+bql[0][2]));
		xu=(2*bql[2][0]-3*bql[2][1]+bql[2][2])/(3*(bql[2][0]-2*bql[2][1]+bql[2][2]));
		yl=(2*bql[0][0]-3*bql[1][0]+bql[2][0])/(3*(bql[0][0]-2*bql[1][0]+bql[2][0]));
		yr=(2*bql[0][2]-3*bql[1][2]+bql[2][2])/(3*(bql[0][2]-2*bql[1][2]+bql[2][2]));

		det=1.0-(xu-xd)*(yr-yl);
		if ( det<=0.0 )	continue;
		det=1.0/det;
		x=det*(xd+yl*(xu-xd));
		y=det*(yl+xd*(yr-yl));

		biquad_poly_coefficients(bql,pcf);
		x0=x,y0=y;
		for ( k=0 ; k<NUM_ITER ; k++ )
		 {	/*fprintf(stderr,"(%12g %12g)\n",x0,y0);*/
			biquad_poly_coeff_rearrange(pcf,x,y,pcf);
			det=4*pcf[2]*pcf[6]-pcf[4]*pcf[4];
			x=-(+2*pcf[6]*pcf[1]-pcf[4]*pcf[3])/det;
			y=-(-pcf[4]*pcf[1]+2*pcf[2]*pcf[3])/det;
			x0+=x,y0+=y;
		 }
		if ( x0<0.0 || y0<0.0 || x0>=1.0 || y0>=1.0 )	continue;
		/*fprintf(stderr,"-> (%12g %12g)\n",x0,y0);*/

		x=x0,x0+=(double)j;
		y=y0,y0+=(double)i;
		peak=biquad_eval(bql,x,y);
		
		cands=(candidate *)realloc(cands,sizeof(candidate)*(ncand+1));

		wc=&cands[ncand];

		wc->ix=j,wc->iy=i;
		
		wc->cx=x0,
		wc->cy=y0;
		wc->peak=peak;

		wc->sxx=wc->syy=wc->sxy=0.0;

		wc->marked=0;
		wc->ipoints=NULL;
		wc->nipoint=0;
		
		ncand++;
	 }
  }

 tensor_free(bqc);

 if ( rcands != NULL )	*rcands=cands;
 if ( rncand != NULL )	*rncand=ncand;
 
 return(0);
}
