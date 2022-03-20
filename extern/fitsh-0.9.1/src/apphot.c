/*****************************************************************************/
/* apphot.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Low-level aperture-photometry routines.				     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <fits/fits.h>

#include "fi.h"

#include "math/spline/spline.h"
#include "math/spline/biquad.h"
#include "math/spline/biquad-isc.h"
#include "math/intersec/intersec.h"
#include "math/intersec/intersec-cri.h"
#include "statistics.h"
#include "tensor.h"

#include "apphot.h"

/*****************************************************************************/

static int fracpixel_compare(const void *p1,const void *p2)
{
 if ( ((fracpixel *)p1)->flux < ((fracpixel *)p2)->flux )	return(-1);
 else								return(1);
}
int fracpixel_order(fracpixel *fp,int n)
{
 qsort((void *)fp,n,sizeof(fracpixel),fracpixel_compare);
 return(0);
}

#define	FRACPIXEL_MEDIAN_SPLINE

double fracpixel_median(fracpixel *fp,int n)
{
 int	i,j;
 double a,area,halfval,median;
 double	*x,*y,*y2;

 x=y=y2=NULL;
 fracpixel_order(fp,n);
 for ( i=0,area=0.0 ; i<n ; i++ )	area+=fp[i].area;

#ifdef	FRACPIXEL_MEDIAN_SPLINE
 x=(double *)malloc(sizeof(double)*(3*n));
 y=x+n,y2=y+n;
 halfval=a=0.0;
 for ( i=0,j=0 ; i<n ; i++ )
  {	double	wa;
	wa=fp[i].area;
	x[i]=a+0.5*wa;
	y[i]=fp[i].flux;
	a+=wa;
  }
 spline_coeff(x,y,n,NULL,NULL,y2);
 median=spline_inter(x,y,y2,n,0.5*area);
 free(x);
#else
 a=0.0;j=-1;
 halfval=0.5*area;
 for ( i=0 ; i<nbgx && j<0 ; i++ )
  {	double	wa;
	wa=fp[i].area;
	if ( a <= halfval && halfval <= a+wa )	j=i;
	a+=wa;
  }
 median=fp[i].flux;
#endif

 return(median);
}
int fracpixel_stat(fracpixel *fp,int n,double *rs,double *rmean,double *rsigma)
{
 double	s,sf,sff;
 int	i;
 s=sf=sff=0.0;
 for ( i=0 ; i<n ; i++ )
  {	double	wa,wf;
	wa=fp[i].area;
	wf=fp[i].flux;
	s  +=wa,
	sf +=wa*wf,
	sff+=wa*wf*wf;
  }
 sf/=s,sff/=s;
 if ( rs != NULL )	*rs    =s;
 if ( rmean != NULL )	*rmean =sf;
 if ( rsigma != NULL )	*rsigma=sqrt(sff-sf*sf);
 return(0);
}

int fracpixel_reject(fracpixel *fp,int n,double lower,double higher)
{
 int	i;
 for ( i=0 ; i<n ; )
  {	if ( fp[i].flux < lower || higher < fp[i].flux )
	 {	if ( i<n-1 )	memmove(fp+i,fp+i+1,sizeof(fracpixel)*(n-1-i));
		n--;
	 }
	else	i++;
  }
 return(n);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* #define	DEBUG_APPHOT_BACK */

int aperture_photometry_back_int(double **data,char **mask,int sx,int sy,double cx,double cy,
	apphotpar *par,double *rbgarea,double *rbgflux,double *rbgavg,double *rbgsigma)
{
 int		i,j,ix,iy,n;
 int		hsan;
 double 	ra,da;

 double		arann,flann,prevarann; 
 double		wr,ril2,rol2;
 double		bgavg,bgmean,bgsigma;
 dcircle	annin,annout;

 fracpixel	*bgxs;
 int		nbgx;

 int		n_iterations;
 double		reject_lower,reject_upper;
 
 if ( data==NULL )	return(1);

 ra  =par->ra,
 da  =par->da;

 hsan=(int)(ra+da+3.0);
 ix=(int)floor(cx);
 iy=(int)floor(cy);

 annin.x=annout.x=cx,
 annin.y=annout.y=cy;
    annin.radius=ra;
   annout.radius=ra+da;

 n=(2*hsan+1)*(2*hsan+1);
 bgxs=(fracpixel *)malloc(sizeof(fracpixel)*n);
 nbgx=0;
 wr=sqrt(0.5);
 ril2=(ra-wr)*(ra-wr),
 rol2=(ra+da+wr)*(ra+da+wr);
 flann=arann=0.0;
 for ( i=iy-hsan ; i<=iy+hsan ; i++ )
  { if ( i<0 || i>=sy )	continue;
    for ( j=ix-hsan ; j<=ix+hsan ; j++ )
     {	drectangle 	pixel;
	double		a,dx,dy,d2,p;

	if ( j<0 || j>=sx )	continue;

	if ( mask != NULL && mask[i][j] )	continue;

	pixel.x=(double)j,
	pixel.y=(double)i;

	p=data[i][j];

	dx=(double)j-cx,
	dy=(double)i-cy;
	d2=dx*dx+dy*dy;
	if ( d2 <= ril2 || rol2 <= d2 )
		continue;

	pixel.sx=pixel.sy=1.0;

	a=area_intersec_of_rect_circ(&pixel,&annout)-
	  area_intersec_of_rect_circ(&pixel,&annin );
	if ( a<=0.0 )	continue;

	bgxs[nbgx].area=a,
	bgxs[nbgx].flux=p;
	flann += a * p,
	arann += a;
	nbgx++;
     }
  }

 if ( nbgx==0 )
  {	free(bgxs);return(1);		}

/* Read iterative sigma rejection parameters from (apphotpar *)par: */
 n_iterations=par->bgm.rejniter;
 reject_lower=par->bgm.rejlower;
 reject_upper=par->bgm.rejupper;
/* Discard and reset illegal values of the sigma rejection parameters: */
 if ( n_iterations<1 )		n_iterations=1;
 if ( reject_lower<=0.0 )	n_iterations=1;
 if ( reject_upper<=0.0 )	n_iterations=1;

 prevarann=bgavg=0.0;

 #ifdef DEBUG_APPHOT_BACK
 fprintf(stderr,"[%g",arann);
 #endif
 bgsigma=0.0;
 for ( i=0 ; i<n_iterations ; i++ )
  {	
	if ( nbgx==0 )
	 {	free(bgxs);return(1);		}

	fracpixel_stat(bgxs,nbgx,&arann,&bgmean,&bgsigma);
	if ( i>0 && prevarann==arann )	break;
	prevarann=arann;
	
	bgavg=fracpixel_median(bgxs,nbgx);
	#ifdef DEBUG_APPHOT_BACK
	fprintf(stderr,"(%g,%g,%g)",arann,bgavg,bgsigma);
	#endif
	if ( i<n_iterations-1 )
	 {	double	rl,rh;
		rl=bgavg - bgsigma*reject_lower;
		rh=bgavg + bgsigma*reject_upper;
		nbgx=fracpixel_reject(bgxs,nbgx,rl,rh);
	 }
  }
 #ifdef DEBUG_APPHOT_BACK
 fprintf(stderr,"]\n");
 #endif

 *rbgarea=arann;
 *rbgflux=flann;
 *rbgavg=bgavg;
 *rbgsigma=bgsigma;

 free(bgxs); 

 return(0);
}

int get_histogram_peak(double *list,int n,double width,double *ret)
{
 double	fmin,fmax,cv,sum;
 int	kmin,kmax,lmin,i;

 if ( list==NULL || n<=0 )	return(0);
 
 cv=fmin=list[0];
 fmax=list[n-1];
 kmin=kmax=lmin=0;
 
 for ( i=0 ; i<n ; i++ )
  {	if ( list[i]>=cv+width )
	 {	if ( i-lmin > kmax-kmin )
		 {	kmin=lmin;
			kmax=i;
		 }
		lmin=i;
		cv+=width;
	 }
  }
 if ( n-lmin > kmax-kmin )
  {	kmin=lmin;
	kmax=n;
  }

 sum=0.0;
 for ( i=kmin ; i<kmax ; i++ )
  {	sum+=list[i];		}
 sum/=(double)(kmax-kmin);

 if ( ret != NULL )	*ret=sum;

 return(0); 
}

int aperture_photometry_back(double **data,char **mask,int sx,int sy,
	double cx,double cy,apphotpar *par,
	double *rbgarea,double *rbgflux,double *rbgavg,double *rbgsigma,
	int *ratot,int *rabad,char **ringmask)
{
 int		i,j,ix,iy,hsize,fsize,nbg,atot,abad;
 double 	ra,da,*bgs,*bds,dx,dy,dd,ril2,rol2;
 double		med,sig,std,sum,mean,mode,min,max,avg,flx,p,w;

 int		n_iterations,type,fbg,lbg,n;
 double		reject_lower,reject_upper;

 if ( data==NULL )	return(1);

 ra  =par->ra,
 da  =par->da;

 hsize=(int)(ra+da+2.0);
 fsize=2*hsize+1;
 ix=(int)floor(cx);
 iy=(int)floor(cy);

 bgs=(double *)malloc(sizeof(double)*fsize*fsize);
 nbg=0;
 ril2=ra*ra,
 rol2=(ra+da)*(ra+da);
 atot=abad=0;

 for ( i=iy-hsize ; i<=iy+hsize ; i++ )
  { 
    for ( j=ix-hsize ; j<=ix+hsize ; j++ )
     {	
	dx=(double)j+0.5-cx,
	dy=(double)i+0.5-cy;
	dd=dx*dx+dy*dy;
	if ( dd < ril2 || rol2 < dd )		continue;
	atot++;

	if ( i<0 || i>=sy || j<0 || j>=sx )	{ abad++;continue;	}
	if ( mask != NULL && mask[i][j]   )	{ abad++;continue;	}
	if ( ringmask!= NULL && ringmask[i][j] ){ abad++;continue;	}

	p=data[i][j];
	bgs[nbg]=p;
	nbg++;
     }
  }

 if ( ratot != NULL )	*ratot=atot;
 if ( rabad != NULL )	*rabad=abad;

 if ( nbg<=2 )
  {	free(bgs);return(1);		}

/* Read iterative sigma rejection parameters from (apphotpar *)par: */
 type=par->bgm.type;
 n_iterations=par->bgm.rejniter;
 reject_lower=par->bgm.rejlower;
 reject_upper=par->bgm.rejupper;
/* Discard and reset illegal values of the sigma rejection parameters: */
 if ( n_iterations<1 )		n_iterations=1;
 if ( reject_lower<=0.0 )	n_iterations=1;
 if ( reject_upper<=0.0 )	n_iterations=1;

 bds=(double *)malloc(sizeof(double)*fsize*fsize);
 fbg=0;lbg=nbg-1;
 median(bgs,nbg);

/*
 avg=0.0;
 get_histogram_peak(bgs,nbg,10.0,&avg);
 flx=0.0;
 for ( i=0 ; i<nbg ; i++ )
  {	flx+=bgs[i];		}
 sig=0.0;
*/

 flx=sig=std=avg=0.0;
 while ( n_iterations>0 )
  {	n=lbg-fbg+1;
	for ( i=fbg,sum=0.0 ; i<=lbg ; i++ )
	 {	sum+=bgs[i];			}
	min=bgs[fbg];
	max=bgs[lbg];
	mean=sum/(double)n;
	med=0.5*(bgs[fbg+n/2]+bgs[fbg+(n-1)/2]);
	flx=sum;
	switch ( type )
	 {   case BGTYPE_MEAN:
		avg=mean;
		break;
	     case BGTYPE_MEDIAN:
		avg=med;
		break;
	     case BGTYPE_MODE:
		mode=3.0*med-2.0*mean;
		if ( mode<min )	mode=min;
		if ( mode>max )	mode=max;
		avg=mode;
		break;
	     default:
		avg=mean;
		break;
	 }

	std=0.0;
	for ( i=0 ; i<n ; i++ )
	 {	w=avg-bgs[i+fbg];
	 	bds[i]=fabs(w);
		std+=w*w;
	 }
	if ( n>0 )	std/=(double)n;	
	else		std=0.0;
	sig=median(bds,n);

	n_iterations--;
	if ( n_iterations>0 && n>2 )
	 {	i=0;
		while ( fbg<lbg && bgs[fbg]<avg-reject_lower*sig ) fbg++,i++;
		while ( fbg<lbg && bgs[lbg]>avg+reject_upper*sig ) lbg--,i++;
		if ( i==0 )	n_iterations=0;
	 }
  }

 *rbgarea=(double)nbg;
 *rbgflux=flx;
 *rbgavg=avg;
 *rbgsigma=sqrt(std);

 free(bds);
 free(bgs);

 return(0);
}

double weighted_intersection(drectangle *pixel,dcircle *aperture,
	double **subpixeldata,int gx,int gy)
{
 int		i,j;
 double		a;
 drectangle	p;

 a=0.0;

 p.sx=pixel->sx/(double)gx;
 p.sy=pixel->sy/(double)gy;

 for ( i=0 ; i<gy ; i++ )
  {	p.y=pixel->y+(double)i*pixel->sy/(double)gy;
	for ( j=0 ; j<gx ; j++ )
	 {	p.x=pixel->x+(double)j*pixel->sx/(double)gx;
	 	a+=area_intersec_of_rect_circ(&p,aperture)*subpixeldata[i][j];
	 }
  }
 
 return(a);
}

typedef struct
 {	int	hsap;
	double	dx,dy,r0;
	double	**arr;
 } apphotcache;

#define		AC_TOLERANCE		(1e-7)

int aperture_photometry_flux(double **data,char **mask,int sx,int sy,
	double cx,double cy,double r0,double *rarea,double *rflux,
	apphot_out *out,double background,
	int *rrtot,int *rrbad,int *rrign,
	int maskignore,double **subpixeldata,int subg,char **xmask)
{
 int			i,j,ix,iy;
 int			hsap;
 int			cflag,flag,rtot,rbad,rign;

 double			flux,area,fcoeff[6]; 
 dcircle		aperture;

 static apphotcache	ac = { 0, 0.0, 0.0, 0.0, NULL };

 if ( data==NULL )	return(1);
 if ( rflux==NULL )	return(0);

 hsap=(int)(r0+2.0);
 ix=(int)floor(cx);
 iy=(int)floor(cy);

 aperture.x=cx,
 aperture.y=cy;
 aperture.radius=r0;
 
 area=flux=0.0;

 flag=0,rtot=rbad=rign=0;

 if ( subg<=0 )	subpixeldata=NULL;
 else if ( subpixeldata==NULL )	subg=0;

 if ( out != NULL )
  {	out->fw=0.0;
	out->fwx=out->fwy=0.0;
	out->fwxx=out->fwxy=out->fwyy=0.0;
  }

 if ( subpixeldata==NULL )
  {	dcircle		dc;
	drectangle	dr;
	double		dx,dy;

	dx=cx-(double)ix;
	dy=cy-(double)iy;
	if ( ! ( ac.r0==r0 && ac.hsap==hsap && 
	fabs(dx-ac.dx)<AC_TOLERANCE && fabs(dy-ac.dy)<AC_TOLERANCE ) )
	 {	ac.r0=r0;
		ac.dx=dx;
		ac.dy=dy;
		ac.hsap=hsap;
		if ( ac.arr != NULL )
		 {	tensor_free(ac.arr);
			ac.arr=NULL;
		 }
		ac.arr=tensor_alloc_2d(double,2*hsap+1,2*hsap+1);
		dr.sx=dr.sy=1.0;	
		dc.x=dx;
		dc.y=dy;
		dc.radius=r0;
		for ( i=-hsap ; i<=hsap ; i++ )
		 {	for ( j=-hsap ; j<=hsap ; j++ )
			 {	dr.x=(double)j;
				dr.y=(double)i;
				ac.arr[i+hsap][j+hsap]=area_intersec_of_rect_circ(&dr,&dc);
			 }
		 }
	 }

  }

 for ( i=iy-hsap ; i<=iy+hsap ; i++ )
  { for ( j=ix-hsap ; j<=ix+hsap ; j++ )
     {	drectangle 	pixel;
	double		a;

	pixel.x=(double)j,
	pixel.y=(double)i;
	pixel.sx=pixel.sy=1.0;
	if ( subpixeldata==NULL && ac.arr==NULL )
		a=area_intersec_of_rect_circ(&pixel,&aperture);
	else if ( subpixeldata==NULL )
		a=ac.arr[i-iy+hsap][j-ix+hsap];
	else
		a=weighted_intersection(&pixel,&aperture,subpixeldata,subg,subg);
	cflag=0;

	if ( out != NULL )
		intersec_cri_integrate_monoms(pixel.x-cx,pixel.y-cy,1.0,1.0,r0,fcoeff,2);

	if ( a>0.0 )
	 {	if ( i<0 || i>=sy || j<0 || j>=sx )
			cflag |= MASK_OUTER,a=0.0;
		else
		 {	if ( mask != NULL && mask[i][j] )
				cflag |= mask[i][j];
			if ( xmask != NULL && xmask[i][j] )
				cflag |= xmask[i][j];
		 }
	 }
	else	continue;

	if ( (! (cflag & MASK_OUTER)) && (cflag&(~maskignore))==0 )
	 {	flux += a * data[i][j];
		area += a;
		if ( out != NULL )
		 {	double	d;
			d=data[i][j]-background;
			out->fw   += fcoeff[0]*d;
			out->fwx  += fcoeff[1]*d;
			out->fwy  += fcoeff[2]*d;
			out->fwxx += fcoeff[3]*d;
			out->fwxy += fcoeff[4]*d;
			out->fwyy += fcoeff[5]*d;
		 }
	 }

	rtot++;
	if ( cflag )
	 {	rbad++;
		if ( (cflag&(~maskignore))==0 )	rign++;
	 }

	flag |= cflag;
     }
  }

 if ( rrtot != NULL )	*rrtot=rtot;
 if ( rrbad != NULL )	*rrbad=rbad;
 if ( rrign != NULL )	*rrign=rign;

 *rarea=area;
 *rflux=flux;

 return(flag);
}

int aperture_photometry(fitsimage *img,char **mask,
	double cx,double cy,apphotpar *par,double *rflux,double *rfluxerr)
{
 double		bgarea,bgflux,bgavg,bgsigma;
 double		area,flux,fluxerr;
 int		r;

 if ( img==NULL || img->data==NULL )	return(1);
 if ( img->sx<=0 || img->sy<=0 )	return(1);

 r=aperture_photometry_back(img->data,mask,img->sx,img->sy,cx,cy,par,&bgarea,&bgflux,&bgavg,&bgsigma,NULL,NULL,NULL);
 if ( r )	return(1);
 r=aperture_photometry_flux(img->data,mask,img->sx,img->sy,cx,cy,par->r0,&area,&flux,NULL,0.0,NULL,NULL,NULL,0,NULL,0,NULL);
 if ( r )	return(1);

 flux -= bgavg * area;
 *rflux =flux;

 if ( rfluxerr != NULL )
  {	fluxerr  = sqrt ( flux + area * (bgsigma*bgsigma)*(1.0+1.0/bgarea) );
	*rfluxerr=fluxerr/sqrt(par->gain);
  }

 return(0);
}

/*****************************************************************************/

int aperture_photometry_flux_biquad(double **c,char **mask,int sx,int sy,
	double cx,double cy,double r0,double *rarea,double *rflux,
	int *rrtot,int *rrbad,int *rrign,int maskignore,char **xmask)
{
 int		i,j,ix,iy,rtot,rbad,flag,rign,cflag;
 int		hsap;

 double		flux,area; 
 dcircle	aperture;

 if ( c==NULL )		return(1);
 if ( rflux==NULL )	return(0);

 hsap=(int)(r0+2.0);
 ix=(int)floor(cx);
 iy=(int)floor(cy);

 aperture.x=cx,
 aperture.y=cy;
 aperture.radius=r0;
 
 area=flux=0.0;
 rtot=rbad=rign=0;flag=0;

 for ( i=iy-hsap ; i<=iy+hsap ; i++ )
  { for ( j=ix-hsap ; j<=ix+hsap ; j++ )
     {	drectangle 	pixel;
	double		a,f;

	pixel.x=(double)j,
	pixel.y=(double)i;
	pixel.sx=pixel.sy=1.0;
	a=area_intersec_of_rect_circ(&pixel,&aperture);
	cflag=0;
	if ( a>0.0 )
	 {	if ( i<0 || i>=sy || j<0 || j>=sx )
			cflag |= MASK_OUTER,a=0.0;
		if (  mask != NULL && ! (cflag&MASK_OUTER) &&  mask[i][j] )
			cflag |=  mask[i][j];
		if ( xmask != NULL && ! (cflag&MASK_OUTER) && xmask[i][j] )
			cflag |= xmask[i][j];
	 }
	else	continue;

	if ( ! (cflag & MASK_OUTER) && (cflag&(~maskignore))==0 )
	 {	area += a;
		f=biquad_isc_int_pixel_circle(c,j,i,cx,cy,r0);
		flux += f;
	 }

	rtot++;
	if ( cflag )
	 {	rbad++;
		if ( (cflag&(~maskignore))==0 )	rign++;
	 }

	flag |= cflag;
     }
  }
 if ( rrtot != NULL )	*rrtot=rtot;
 if ( rrbad != NULL )	*rrbad=rbad;
 if ( rrign != NULL )	*rrign=rign;

 *rarea=area;
 *rflux=flux;

 return(flag);
}

/*****************************************************************************/
          
