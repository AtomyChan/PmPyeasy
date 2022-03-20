/*****************************************************************************/
/* fiinfo-image.c							     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line user interface to get some statistics from FITS data:	     */
/* Gather information about FITS images.				     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>

#include <fits/fits.h>

#include "fi.h"

#include "fitsmask.h"
#include "math/spline/biquad.h"
#include "math/spline/bicubic.h"
#include "math/spline/spline.h"
#include "math/fit/lmfit.h"
#include "statistics.h"
#include "math/poly.h"
#include "math/polyfit.h"
#include "math/splinefit.h"
#include "io/iof.h"
#include "io/scanarg.h"
#include "io/tokenize.h"
#include "link/linkpoint.h"
#include "tensor.h"
#include "common.h"

#include "fiinfo.h"

/*****************************************************************************/

int create_link_background(fitsimage *img,char **mask,fitsimage *bgi,int nx,int ny)
{
 double		**bqc;
 int		sx,sy,i,j,mi,mj /*,k,l,ct,cm */;
 linkpoint	**lparr;

 if ( img==NULL || img->data==NULL )	return(1);
 if ( bgi==NULL || bgi->data==NULL )	return(1);
 sx=img->sx,sy=img->sy;
 if ( sx<=0 || sy<=0 )			return(1);
 if ( bgi->sx != sx || bgi->sy != sy )	return(1);

 bqc=tensor_alloc_2d(double,2*sx+1,2*sy+1);
 biquad_coeff(img->data,sx,sy,bqc,mask);
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( mask != NULL && mask[i][j] )	continue;
		bgi->data[i][j]=biquad_scatter(bqc,j,i);
	 }
  }
 lparr=linkpoint_create(bgi->data,sx,sy,NULL,mask,-1);
 linkpoint_reconnect(lparr,sx,sy);

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( mask != NULL && mask[i][j] )
		 {	bgi->data[i][j]=0.0;continue;		}
		mi=lparr[i][j].my,
		mj=lparr[i][j].mx;
		/*ct=cm=0;
		for ( k=-1 ; k<=1 ; k++ )
		 {	if ( i+k<0 || i+k>=sy )
				continue;
			for ( l=-1 ; l<=1 ; l++ )
			 {	if ( j+l<0 || j+l>=sx || ( k==0 && l==0 ) )
					continue;
				if ( mask != NULL && mask[i+k][j+l] )
					continue;
				ct++;
				if ( lparr[i+k][j+l].mx==mj && lparr[i+k][j+l].my==mi )
					cm++;
			 }
		 }
		if ( cm>=1 || ct<1 )*/
			bgi->data[i][j]=img->data[mi][mj];
		/*else
		 {	bgi->data[i][j]=0.0;
			cm=ct=0;
			for ( k=-1 ; k<=1 ; k++ )
			 {	if ( i+k<0 || i+k>=sy )
					continue;
				for ( l=-1 ; l<=1 ; l++ )
				 {	if ( j+l<0 || j+l>=sx || ( k==0 && l==0 ) )
						continue;
					if ( mask != NULL && mask[i+k][j+l] )
						continue;
					ct++;
					if ( lparr[i+k][j+l].mx != mj || lparr[i+k][j+l].my != mi )
					 {	bgi->data[i][j]+=img->data[lparr[i+k][j+l].my][lparr[i+k][j+l].mx];
						cm++;
					 }
				 }
			 }
			if ( cm>0 )
				bgi->data[i][j]/=(double)cm;
		 }*/
		
	 }
  }

 /* do spline smoothing */
 if ( nx>0 && ny>0 )
  {	double	**scc,**bcc,x,y;
	point	*points;
	int	npoint;

	scc=tensor_alloc_2d(double,nx+1,ny+1);      /* spline control points */
	bcc=tensor_alloc_2d(double,2*nx+2,2*ny+2);  /* bicubic coeffs	     */
	points=(point *)malloc(sizeof(point)*sx*sy);/* points for fitting... */
	npoint=0;
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0; j<sx ; j++ )
		 {	if ( mask != NULL && mask[i][j] )
				continue;
			if ( i%16 || j%16 )	continue;
			points[npoint].x=(double)j+0.5;
			points[npoint].y=(double)i+0.5;
			points[npoint].value=bgi->data[i][j];
			points[npoint].weight=1.0;
			npoint++;
		 }
	 }
	fit_2d_spline(0,sx,nx,0,sy,ny,points,npoint,scc);
	bicubic_coeff(scc,nx+1,ny+1,bcc,NULL);
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	x=(double)nx*((double)j+0.5)/(double)sx;
			y=(double)nx*((double)i+0.5)/(double)sy;
			bgi->data[i][j]=bicubic_inter(bcc,x,y);
		 }
	 }
	free(points);
	tensor_free(bcc);
	tensor_free(scc);
  }
 tensor_free(lparr);
 return(0);

}
/*****************************************************************************/

int fits_stat_basic(fitsimage *img,char **mask,double *rmin,double *rmax,double *rmean,double *rstdd)
{
 int	i,j,n,is_first;
 double	mean,mean2,min,max,stdd,w;

 min=max=mean=mean2=0.0;n=0;
 is_first=1;
 for ( i=0 ; i<img->sy ; i++ )
  {	for ( j=0 ; j<img->sx ; j++ )
	 {	w=img->data[i][j];
		if ( mask != NULL && mask[i][j] )	continue;
		mean +=w;
		mean2+=w*w;
		if ( is_first )		min=max=w,is_first=0;
		else if ( w<min )	min=w;
		else if ( w>max )	max=w;
		n++;
	 }
  }
 mean /=(double)n;
 mean2/=(double)n;
 stdd  =sqrt(mean2-mean*mean);

 if ( rmin  != NULL )	*rmin =min;
 if ( rmax  != NULL )	*rmax =max;
 if ( rmean != NULL )	*rmean=mean;
 if ( rstdd != NULL )	*rstdd=stdd;

 return(0);
}

double estimate_skysigma_naive(double *rawdata,int k,double sky,double smm,double spp)
{
 double sp,sp2,ssigma,w;
 int	i,j;

 sp=0.0;sp2=0.0;j=0;
 for ( i=0 ; i<k ; i++ )
  {	w=rawdata[i];
	if ( sky+smm <= w && w <= sky+spp )
	 {	sp+=w;sp2+=w*w;j++;			}
  }
 sp /=(double)j;
 sp2/=(double)j;
 ssigma=sqrt(sp2-sp*sp);

 return(ssigma);
}

int fits_stat_raw_sky_skysigma(double *data,int ndat,double stddev,double *rsky,double *rskysigma)
{
 int	i,j,m,n;
 double	sky,skysigma,ss0;
 double	sky_max,sky_fit,sigr;

 double *histo,dh,dlon,h0,w,hfl;
 int	nhist;

 double	**amatrix,*bvector; 
 double	fvars[4],f; int	nvar;
 double	fp,fq,fr,fs,x1,x2,x0,xmx,hmx;

 sky=median(data,ndat);
 ss0=estimate_skysigma_naive(data,ndat,sky,-stddev,0.0);
 ss0=estimate_skysigma_naive(data,ndat,sky,-ss0*2.0,+ss0*2.0);

 dh=1.0,dlon=3.0;	/* global parameters of sky & skysigma estimation!!! */
 hfl=0.5;		/* ... why not ... */

 nhist=(int)(ss0*2.0*dlon/dh);
 h0=sky-ss0*dlon;
 histo=(double *)malloc(sizeof(double)*nhist);
 for ( i=0 ; i<nhist ; i++ )	histo[i]=0.0;

 for ( i=0 ; i<ndat ; i++ )
  {	w=data[i];
	if ( w<h0 )	continue;
	j=(w-h0)/dh;
	if ( j>=nhist )	continue;
	histo[j]+=1.0;
  }

/* for ( i=0 ; i<nhist ; i++ ) */
/*  {	fprintf(stdout,"%g %g\n",h0+((double)i+0.5)*dh,histo[i]);	} */

 for ( i=1,j=0 ; i<nhist ; i++ )
  {	if ( histo[i]>histo[j] )	j=i;		}
 sky_max=h0+((double)j+0.5)*dh;
 hmx=histo[j];

 nvar=4;
 amatrix=matrix_alloc(nvar);
 bvector=vector_alloc(nvar);
 for ( i=0 ; i<nvar ; i++ )
  {	for ( j=0 ; j<nvar ; j++ )
	 {	amatrix[i][j]=0.0;		}
	bvector[i]=0.0;
  }
 for ( i=0 ; i<nhist/2 ; i++ )
  {	if ( histo[i] < hmx * hfl )	continue;
	w=h0+((double)i+0.5)*dh;
	fvars[3]=1.0;
	fvars[2]=w;
	fvars[1]=w*w;
	fvars[0]=w*fvars[1];
	f=histo[i];
	for ( m=0 ; m<nvar ; m++ )
	 {	bvector[m]+=f*fvars[m];		}
	for ( m=0 ; m<nvar ; m++ )
	 {	for ( n=0 ; n<nvar ; n++ )
		 {	amatrix[m][n]+=fvars[m]*fvars[n];	}
 	 }
  }
 i=solve_gauss(amatrix,bvector,nvar);
 matrix_free(amatrix);
 fp=bvector[0],fq=bvector[1],fr=bvector[2],fs=bvector[3];
 vector_free(bvector);

 w=4*fq*fq-12*fp*fr;
 if ( w < 0.0 )	return(1);	/* defective histogram...  */
 w=sqrt(w);
 x1=(-2.0*fq+w)/(6*fp);
 x2=(-2.0*fq-w)/(6*fp);
 sigr=-(6*fp*x1+2*fq);
 if ( sigr > 0.0 )	x0=x1;
 else			x0=x2;
 
 sky_fit=x0;
 xmx=((fp*x0+fq)*x0+fr)*x0+fs;
 sigr=-((6*fp*x0+2*fq)/xmx);

 skysigma=1/sqrt(sigr);
 
 if ( rsky      != NULL )	*rsky     =sky_fit;
 if ( rskysigma != NULL )	*rskysigma=skysigma;

 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int fits_stat_sky(fitsimage *img,double stddev,double *rsky,double *rskysigma)
{
 int	i,j,k,sx,sy;
 double	*rawdata,sky,skysigma;

 sx=img->sx,sy=img->sy;
  
 rawdata=(double *)malloc(sizeof(double)*sx*sy);
 for ( i=0,k=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	rawdata[k]=img->data[i][j];
		k++;
	 }
  }
 fits_stat_raw_sky_skysigma(rawdata,sx*sy,stddev,&sky,&skysigma);
 free(rawdata);

 if ( rsky      != NULL )	*rsky     =sky;
 if ( rskysigma != NULL )	*rskysigma=skysigma;

 return(0);
}

double fits_stat_median(fitsimage *img)
{
 int	i,j,k,sx,sy;
 double	*rawdata,med;

 sx=img->sx,sy=img->sy;
  
 rawdata=(double *)malloc(sizeof(double)*sx*sy);
 for ( i=0,k=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	rawdata[k]=img->data[i][j];
		k++;
	 }
  }
 med=median(rawdata,sx*sy);
 free(rawdata);

 return(med);
}

int fits_stat_background(fitsimage *img,int nbx,int nby,point **pdsky,point **pdsigma,double stddev)
{
 int	i,j,k,sx,sy;
 int	ii,jj,in,jn,imin,imax,jmin,jmax,tn;
 double	*rdata,x,y,val,sig;
 
 sx=img->sx,sy=img->sy;

 for ( i=0 ; i<nby ; i++ )
  {	for ( j=0 ; j<nbx ; j++ )
	 {	imin=sy*i/nby;imax=sy*(i+1)/nby;
		jmin=sx*j/nbx;jmax=sx*(j+1)/nbx;

		in=imax-imin,jn=jmax-jmin;
		tn=in*jn;
		rdata=(double *)malloc(sizeof(double)*tn);
		for ( ii=imin,k=0 ; ii<imax ; ii++ )
		 {	for ( jj=jmin ; jj<jmax ; jj++,k++ )
			 {	rdata[k]=img->data[ii][jj];	}
		 }

		fits_stat_raw_sky_skysigma(rdata,tn,stddev,&val,&sig);

		free(rdata);
		y=0.5*(double)(imax+imin-1);
		x=0.5*(double)(jmax+jmin-1);
		pdsigma[i][j].x=pdsky[i][j].x=x/(double)sx;
		pdsigma[i][j].y=pdsky[i][j].y=y/(double)sy;
		pdsky  [i][j].value=val;
		pdsigma[i][j].value=sig;
		pdsky  [i][j].weight=1.0;
		pdsigma[i][j].weight=1.0;
	  }
  }
 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int fits_stat_sky_like_isis(fitsimage *img,double stddev,double *rsky,double *rskysigma)
{
 int	i,j,k,sx,sy;
 double	*rawdata,sky,skysigma;

 sx=img->sx,sy=img->sy;
  
 rawdata=(double *)malloc(sizeof(double)*sx*sy);
 for ( i=0,k=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	rawdata[k]=img->data[i][j];
		k++;
	 }
  }
 sky=median(rawdata,sx*sy);
 for ( i=0 ; i<sx*sy ; i++ )
  {	rawdata[i]=fabs(rawdata[i]-sky);		}
 skysigma=median(rawdata,sx*sy);

 free(rawdata);

 if ( rsky      != NULL )	*rsky     =sky;
 if ( rskysigma != NULL )	*rskysigma=skysigma;

 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int fits_stat_sky_biquad(fitsimage *img,double stddev,double *rsky,double *rskysigma)
{
 double	**bqc,*rawdata,w,sum,sm2,avg,sct,med,rlevel;
 char	**mask;
 int	i,j,sx,sy,n,k,niter;

 if ( img==NULL )	return(1);
 if ( img->data==NULL )	return(1);
 sx=img->sx,sy=img->sy;
 if ( sx<=0 || sy<=0 )	return(1);

 rawdata=(double *)tensor_alloc_1d(double,sx*sy);
 bqc=(double **)tensor_alloc_2d(double,2*sx+1,2*sy+1);
 mask=(char **)tensor_alloc_2d(char,sx,sy);

 biquad_coeff(img->data,sx,sy,bqc,NULL);
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	w=biquad_scatter(bqc,j,i);
		bqc[2*i+1][2*j+1]=w;		/* strange, i know... */
	 }
  }
	
 for ( i=0 ; i<sy ; i++ )
  {	memset(mask[i],0,sx);		}

 niter=3;rlevel=3.0;

 for ( k=0 ; k<niter ; k++ )
  {	for ( i=0,n=0,sum=sm2=0.0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	if ( mask[i][j] )	continue;
			w=bqc[2*i+1][2*j+1];
			sum+=w;sm2+=w*w;rawdata[n]=w;
			n++;
		 }
	 }
	avg=sum/(double)n;
	sm2=sm2/(double)n-avg*avg;
	if ( sm2<=0.0 )	sct=0.0;
	else		sct=sqrt(sm2);
	med=median(rawdata,n);
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	if ( mask[i][j] )	continue;
			w=bqc[2*i+1][2*j+1];
			if ( fabs(w-med)>rlevel*sct )
				mask[i][j]=-1;
		 }
	 }
	if ( k==0 )	fits_mask_expand_false(mask,sx,sy,1,-1,-1,1);
  }				

 for ( i=0,n=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( mask[i][j] )	continue;
		w=bqc[2*i+1][2*j+1];
		if ( med-sct < w && w < med+sct )
			rawdata[n]=img->data[i][j],n++;
	 }
  }

 *rskysigma=med*1.8498;
 *rsky=median(rawdata,n);

 fits_mask_free(mask);
 tensor_free(bqc);
 tensor_free(rawdata);

 return(0);
}

/*****************************************************************************/
