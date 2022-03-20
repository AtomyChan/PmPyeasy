/*****************************************************************************/
/* fiphot.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line tool for performing photometry on FITS images.		     */
/*****************************************************************************/
#define	FI_PHOT_VERSION		"0.9z7"
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>

#include <fits/fits.h>

#include "longhelp.h"
#include "fi.h"

#include "fitsmask.h"
#include "io/iof.h"
#include "io/scanarg.h"
#include "io/tokenize.h"
#include "math/spline/biquad.h"
#include "math/fit/lmfit.h"
#include "math/poly.h"
#include "math/point.h"
#include "statistics.h"
#include "magnitude.h"

#include "basis.h"
#include "tensor.h"
#include "common.h"
#include "apphot.h"
#include "kernel.h"
#include "weight.h"

#include "fiphot.h"

#ifdef	HAVE_NO_CC_EXTENSION
#define	__extension__
#endif

/*****************************************************************************/

int	is_verbose,is_comment;
char	*progbasename;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int fprint_error(char *expr,...)
{
 va_list	ap;
 fprintf(stderr,"%s: error: ",progbasename);
 va_start(ap,expr);
 vfprintf(stderr,expr,ap);
 va_end(ap);
 fprintf(stderr,"\n");
 return(0);
}

int fprint_warning(char *expr,...)
{
 va_list	ap;
 fprintf(stderr,"%s: warning: ",progbasename);
 va_start(ap,expr);
 vfprintf(stderr,expr,ap);
 va_end(ap);
 fprintf(stderr,"\n");
 return(0);
}

/*****************************************************************************/

int get_kernel_block_coords(int id,int *rbx,int *rby)
{
 int	bx,by,fs,bl,hs;

 if ( id<=0 )	bx=by=0;
 else
  {	hs=1,fs=3;bl=8;id--;
	while ( 1 )
	 {	if ( id<bl )
		 {	if ( id<fs )		bx=-hs+id,by=-hs;
			else if ( id>=bl-fs )	bx=-hs+id-bl+fs,by=hs;
			else
			 {	id-=fs;
				by=-hs+1+id/2;
				if ( id%2==0 )	bx=-hs;
				else		bx=+hs;
			 }
			break;
		 }
		else
		 {	id-=bl,hs++,
			fs+=2,bl+=8;
		 }
	 };
  }
 *rbx=bx,
 *rby=by;
 return(0);
}
int get_kernel_block_id(int bx,int by)
{
 int	hs,b;
 if ( bx==0 && by==0 )	return(0);
 hs=0;
 if ( +bx>hs )	hs=+bx;
 if ( -bx>hs )	hs=-bx;
 if ( +by>hs )	hs=+by;
 if ( -by>hs )	hs=-by;
 b=(2*hs-1)*(2*hs-1);
      if ( by==-hs )	return(b+hs+bx);
 else if ( by==+hs )	return(b+7*hs-1+bx);
 else if ( bx==-hs )	return(b+2*hs+1+2*(hs+by-1)+0);
 else if ( bx==+hs )	return(b+2*hs+1+2*(hs+by-1)+1);
 else			return(-1);
}

/*****************************************************************************/

double optimal_aperture_xi_approx(double m)
{
 double	x,xi,b,d,c,q,w;
 x=log(m);
 b=2.316,d=-1.914,q=3.0;
 w=fabs(x-d);
 c=2.512862417;
 xi=sqrt(pow(b*b+pow(w,q),1.0/q)-(x-d)+c);
 return(xi);
}
double optimal_aperture_xi(double m,int n)
{
 double	xi,k,t;
 int	i;
 xi=optimal_aperture_xi_approx(m);
 for ( i=0 ; i<n ; i++ )
  {	t=2.0*m*(1.0+xi*xi)+1.0;
	k=0.5*(t-sqrt(t*t-8.0*m));
	xi=sqrt(-2.0*log(k));
  }
 return(xi);
}

double optimal_aperture(double g,double s,double bg,double bgar,double flux)
{
 double	m,xi,r;
 m=g*(s*s*bg*bg)*(1.0+1.0/bgar)*M_PI/flux;
 xi=optimal_aperture_xi(m,1);
 r=s*xi;
 return(r);
}

/*****************************************************************************/

int read_subpixel_file(FILE *fr,double ***rsubpixeldata,int *rsubg)
{
 double **subpixeldata,x,y,v;
 int	subg,n,i,j,ix,iy;
 char	buff[256],*cmd[4];
 point	*points;
 int	npoint;
 points=NULL;npoint=0;
 while ( ! feof(fr) )
  {	if ( fgets(buff,255,fr)==NULL )	break;
	remove_newlines_and_comments(buff);
	n=tokenize_spaces(buff,cmd,3);
	if ( n<3 )	continue;
	ix=iy=-1;
	sscanf(cmd[0],"%lg",&x);ix=(int)floor(x);
	sscanf(cmd[1],"%lg",&y);iy=(int)floor(y);
	if ( ! isfinite(x) || ! isfinite(y) )	continue;
	if ( ix<0 || iy<0 )			continue;
	sscanf(cmd[2],"%lg",&v);
	if ( ! isfinite(v) )			continue;
	points=(point *)realloc(points,sizeof(point)*(npoint+1));
	points[npoint].x=(double)ix,
	points[npoint].y=(double)iy;
	points[npoint].value=v;
	npoint++;
  };
 subg=0;
 for ( i=0 ; i<npoint ; i++ )
  {	ix=(int)floor(points[i].x),
	iy=(int)floor(points[i].y);
	if ( ix>subg )	subg=ix;
	if ( iy>subg )	subg=iy;
  }
 subg++;
 subpixeldata=tensor_alloc_2d(double,subg,subg);
 for ( i=0 ; i<subg ; i++ )
  { for ( j=0 ; j<subg ; j++ )
     {	subpixeldata[i][j]=1.0;		}
  }
 for ( i=0 ; i<npoint ; i++ )
  {	ix=(int)floor(points[i].x),
	iy=(int)floor(points[i].y);
	subpixeldata[iy][ix]=points[i].value;
  }
 free(points);

 *rsubpixeldata=subpixeldata;
 *rsubg=subg;
	
 return(0);
}
int normalize_subpixeldata(double **d,int g)
{
 int	i,j;
 double	s,n;
 for ( i=0,s=0.0,n=0.0 ; i<g ; i++ )
  {	for ( j=0 ; j<g ; j++,n=n+1.0 )
	 {	s+=d[i][j];		}
  }
 if ( s<=0.0 )	return(1);
 for ( i=0 ; i<g ; i++ )
  {	for ( j=0 ; j<g ; j++ )
	 {	d[i][j]=d[i][j]*n/s;	}
  }

 return(0);
}


/*****************************************************************************/

int calculate_magnitudes(photstar *ps,int np,magflux *mf0)
{
 int		i,j;
 photflux	*pf;
 magflux	*mf,mf1;

 for ( i=0 ; i<np ; i++ )
  {	for ( j=0 ; j<ps[i].n ; j++ )
	 {	pf=&ps[i].fluxes[j];
		if ( ps[i].use_ref && ps[i].rfflux != NULL )
		 {	mf1.magnitude=ps[i].ref_mag;
			mf1.intensity=ps[i].rfflux[j].flux;
			mf=&mf1;
		 }
		else
			mf=mf0;
		flux_to_mag_magerr(pf->flux,pf->fluxerr,mf,&pf->mag,&pf->magerr);
	 }
  }

 return(0); 
}

/*****************************************************************************/

int add_to_data_weight(double **data,char **mask,int sx,int sy,
	int hsize,int grid,weight *wg,double mul)
{
 int	i,j,fsize,ix,iy;

 if ( data==NULL || wg==NULL )	return(1);

 fsize=grid*(2*hsize+1);

 for ( i=0 ; i<fsize ; i++ )
  {	iy=(i/grid)-hsize+wg->iy;
	if ( iy<0 || iy>=sy )	continue;
	for ( j=0 ; j<fsize ; j++ )
	 {	ix=(j/grid)-hsize+wg->ix;
		if ( ix<0 || ix>=sx )			continue;
		if ( mask != NULL && mask[iy][ix] )	continue;
		data[iy][ix] += mul * wg->iarr[i][j];
	 }
  }

 return(0);
}
int add_to_image_weights(fitsimage *img,char **mask,weightlist *wl,double mul)
{
 int	i;

 if ( img==NULL || img->data==NULL )	return(1);
 if ( img->sx<=0 || img->sy<=0 )	return(1);
 if ( wl==NULL || wl->weights==NULL )	return(0);

 for ( i=0 ; i<wl->nweight ; i++ )
  {	add_to_data_weight(img->data,mask,img->sx,img->sy,
		wl->hsize,wl->grid,&wl->weights[i],mul);
  }

 return(0);
}

/*****************************************************************************/

/*****
typedef struct
{	fitsimage	*img;
	char		**mask;
	weightlist	*wl;
	int		weightusage;
	int		is_calc_opt_apert;
	apgeom		*inaps;
	int		numap;
	int		bkhsize;
	double		gain,sigma;
	xphotpar	*xpp;
	double		**bqc;
} phot_parallel_param;
*****/

int ringmask_subtract(char **ringmask,int sx,int sy,double x,double y,double r)
{
 int	i,j,i0,i1,j0,j1;
 double	r2,dx,dy;

 r2=r*r;
 i0=(int)(y-r-1.0); if ( i0<0  )	i0=0;
 i1=(int)(y+r+1.0); if ( i1>sy )	i1=sy;
 j0=(int)(x-r-1.0); if ( j0<0  )	j0=0;
 j1=(int)(x+r+1.0); if ( j1>sx )	j1=sx;
 for ( i=i0 ; i<i1 ; i++ )
  {	dy=(double)(i+0.5)-y;
	for ( j=j0 ; j<j1 ; j++ )
	 {	dx=(double)(j+0.5)-x;
		if ( dx*dx+dy*dy<=r2 && ringmask[i][j] )
			ringmask[i][j]--;
	 }
  }
 return(0);
}

int ringmask_coadd(char **ringmask,int sx,int sy,double x,double y,double r)
{
 int	i,j,i0,i1,j0,j1;
 double	r2,dx,dy;

 r2=r*r;
 i0=(int)(y-r-1.0); if ( i0<0  )	i0=0;
 i1=(int)(y+r+1.0); if ( i1>sy )	i1=sy;
 j0=(int)(x-r-1.0); if ( j0<0  )	j0=0;
 j1=(int)(x+r+1.0); if ( j1>sx )	j1=sx;
 for ( i=i0 ; i<i1 ; i++ )
  {	dy=(double)(i+0.5)-y;
	for ( j=j0 ; j<j1 ; j++ )
	 {	dx=(double)(j+0.5)-x;
		if ( dx*dx+dy*dy<=r2 && ringmask[i][j]<2 )
			ringmask[i][j]++;
	 }
  }
 return(0);
}

/* 
 dradius:
	- positive: use this value as disjoint area radius
	- zero: use the aperture radius as disjoint area radius
	- negative: use the inner radius of the annulus as disjoint area radius
*/

char **	ringmask_create(int sx,int sy,photstar *ps,int np,apgeom *inaps,int a,double dradius)
{
 char **ringmask;
 int	i,n;
 apgeom	*aps;
 double	r;

 ringmask=(char **)tensor_alloc_2d(char,sx,sy);
 for ( i=0 ; i<sy ; i++ )
  {	memset(ringmask[i],0,sx);		}

 for ( n=0 ; n<np ; n++ )
  {	
	if ( dradius>0.0 )
		r=dradius;
	else 
	 {	if ( inaps != NULL )	aps=inaps;
		else			aps=ps[n].inaps;
		if ( aps==NULL )	continue; /* however, unexpected. */
		if ( dradius<0.0 )	r=aps[a].ra;
		else			r=aps[a].r0;
	 }

 	ringmask_coadd(ringmask,sx,sy,ps[n].x,ps[n].y,r);
  }

 return(ringmask);
}

int do_photometry(fitsimage *img,char **mask,
	photstar *ps,int np,weightlist *wl,int weightusage,int is_calc_opt_apert,
	apgeom *inaps,int numap,
	spatialgain *sg,double sigma,xphotpar *xpp)
{
 int		i,j,k,jp,sx,sy,r;
 apphotpar	ap;
 apgeom		*aps;
 photflux	*pf,*wpf;

 double		bgarea,bgflux,bgmedian,bgsigma,gain;
 double		area,flux,fluxerr,cx,cy,dx,dy,cfd2;
 double		**bqc;
 weight		*ww;

 double		**aphdata,**aphwsum,grid2,igrid2;
 char		**aphmask;
 char		***ringmasks;
 int		hsize,fsize,grid,wfsize;

 if ( img==NULL || img->data==NULL )	return(1);
 sx=img->sx,sy=img->sy;
 if ( sx<=0 || sy<=0 )			return(1);

 if ( xpp->use_biquad )
  {	bqc=tensor_alloc_2d(double,sx*2+1,sy*2+1);
	if ( bqc==NULL )	return(-1);
	biquad_coeff(img->data,sx,sy,bqc,NULL);
	weightusage &= ~(USE_WEIGHT_SUBTRACTED|USE_WEIGHT_WEIGHTED);
  }
 else
	bqc=NULL;

 if ( wl==NULL || ! weightusage )	wl=NULL,weightusage=0;

 memcpy(&ap.bgm,&xpp->bgm,sizeof(bgmode));

 if ( wl != NULL )
  {	grid=wl->grid;
	wfsize=grid*(2*wl->hsize+1);
  }
 else
  {	grid=0;
	wfsize=0;
  }

 if ( weightusage & USE_WEIGHT_SUBTRACTED )
  { 	add_to_image_weights(img,mask,wl,-1.0);	
	cfd2=xpp->wconfdist*xpp->wconfdist;
  }
 else	
	cfd2=0.0;

 if ( weightusage & USE_WEIGHT_WEIGHTED )
  {	hsize=0;
	if ( inaps != NULL )
	 {	aps=inaps;
		for ( j=0 ; j<numap ; j++ )
		 {	k=(int)(aps[j].r0+1.0);
			if ( k>hsize )	hsize=k;
		 }
	 }
	else
	 {	for ( i=0 ; i<np ; i++ )
		 {	aps=ps[i].inaps;
			for ( j=0 ; j<numap && aps != NULL ; j++ )
			 {	k=(int)(aps[j].r0+1.0);
				if ( k>hsize )	hsize=k;
			 }
		 }
	 }

	fsize=grid*(2*hsize+1);
	aphdata=(double **)tensor_alloc_2d(double,fsize,fsize);
	aphmask=(char **)tensor_alloc_2d(char,fsize,fsize);
	aphwsum=(double **)tensor_alloc_2d(double,2*wl->hsize+1,2*wl->hsize+1);
  }
 else
  {	hsize=fsize=0;
	grid=1;
	aphdata=NULL;
	aphwsum=NULL;
	aphmask=NULL;
  }

 if ( xpp->is_disjoint_rings )
  {	ringmasks=(char ***)malloc(sizeof(char **)*numap);
	for ( j=0 ; j<numap ; j++ )
	 {	ringmasks[j]=ringmask_create(sx,sy,ps,np,inaps,j,-1.0);		}
  }
 else if ( xpp->is_disjoint_apertures )
  {	ringmasks=(char ***)malloc(sizeof(char **)*numap);
	for ( j=0 ; j<numap ; j++ )
	 {	ringmasks[j]=ringmask_create(sx,sy,ps,np,inaps,j,0.0);		}
  }
 else if ( xpp->disjoint_radius > 0.0 )
  {	ringmasks=(char ***)malloc(sizeof(char **)*numap);
	for ( j=0 ; j<numap ; j++ )
	 {	ringmasks[j]=ringmask_create(sx,sy,ps,np,inaps,j,xpp->disjoint_radius);		}
  }
 else
	ringmasks=NULL;

 grid2=(double)(grid*grid);
 igrid2=1.0/grid2;

 for ( i=0 ; i<np ; i++ )
  {	
	if ( inaps != NULL )	aps=inaps;
	else			aps=ps[i].inaps;

	if ( aps==NULL )	continue;	/* however, unexpected. */

	pf=(photflux *)malloc(sizeof(photflux)*numap);
	ps[i].fluxes=pf;
	ps[i].n=numap;

	cx=ps[i].x,
	cy=ps[i].y;

	gain=eval_2d_poly(cx,cy,sg->order,sg->coeff,0.5*(double)sx,0.5*(double)sy,0.5*(double)sx);
	if ( 0<sg->vmin && gain<sg->vmin )	gain=sg->vmin;
	ap.gain=gain;

	if ( weightusage & USE_WEIGHT_SUBTRACTED )
	 {	ww=weight_get_closest(wl,cx,cy);
		dx=cx-ww->x,dy=cy-ww->y;
		if ( cfd2<=0.0 || dx*dx+dy*dy<cfd2 )
			add_to_data_weight(img->data,mask,sx,sy,wl->hsize,wl->grid,ww,+1.0);
		else	
			ww=NULL;
	 }
	else
		ww=NULL;

/*
	if ( weightusage & USE_WEIGHT_WEIGHTED )
	 {	int	ix0,iy0,i,j;
		int	ix,iy,wx,wy;
		weight	*ww;

		ww=weight_get_closest(wl,cx,cy);

		ix0=((int)floor(cx))-hsize;
		iy0=((int)floor(cy))-hsize;

		for ( i=0 ; i<wfsize/grid ; i++ )
		 {	for ( j=0 ; j<wfsize/grid ; j++ )
			 {	aphwsum[i][j]=0.0;		}
		 }
		for ( i=0 ; i<wfsize ; i++ )
		 {	for ( j=0 ; j<wfsize ; j++ )
			 {	aphwsum[i/grid][j/grid]+=ww->iarr[i][j];	}
		 }
		
		for ( i=0 ; i<fsize ; i++ )
		 {	iy=iy0+i/grid;
			for ( j=0 ; j<fsize ; j++ )
			 {	ix=ix0+j/grid;
				if ( ix<0 || iy<0 || ix>=sx || iy>=sy )
				 {	aphmask[i][j] = MASK_OUTER;
					aphdata[i][j] = 0.0;
					continue;
				 }
				else if ( mask != NULL )
				 {	aphmask[i][j] = mask[iy][ix];	}
				else
				 {	aphmask[i][j] = 0;		}

				aphdata[i][j]=img->data[iy][ix]*igrid2;

				wx=j+grid*(-hsize+wl->hsize);
				wy=i+grid*(-hsize+wl->hsize);
				if ( wx>=0 && wy>=0 && wx<wfsize && wy<wfsize )
				 {	aphdata[i][j]+=ww->iarr[wy][wx];
					aphdata[i][j]-=aphwsum[wy/grid][wx/grid]*igrid2;
				 }
			
			 }
		 }
	 }
*/

	for ( j=0 ; j<numap ; j++ )
	 {	int		rtot,rbad,rign,atot,abad;
		apphot_out	out;
		double		background;

		ap.r0=aps[j].r0,
		ap.ra=aps[j].ra,
		ap.da=aps[j].da;

		wpf=&pf[j];		
		wpf->ag.r0=ap.r0;
		wpf->ag.ra=ap.ra;
		wpf->ag.da=ap.da;

		for ( k=0,jp=-1 ; k<j && jp<0 ; k++ )
		 {	if ( pf[k].ag.ra==ap.ra && pf[k].ag.da==ap.da )
			 {	jp=k;break;		}
		 }

		if ( jp>=0 )
		 {	bgarea=pf[jp].bgarea,
			bgflux=pf[jp].bgflux,
			bgmedian=pf[jp].bgmedian,
			bgsigma=pf[jp].bgsigma;
			atot=pf[jp].atot,
			abad=pf[jp].abad;
			r=0;
		 }
		else if ( (! xpp->use_sky) && ringmasks != NULL )
		 {	r=aperture_photometry_back(img->data,mask,sx,sy,cx,cy,
				&ap,&bgarea,&bgflux,&bgmedian,&bgsigma,
				&atot,&abad,ringmasks[j]);
		 }
		else if ( (! xpp->use_sky) )
		 {	r=aperture_photometry_back(img->data,mask,sx,sy,cx,cy,
				&ap,&bgarea,&bgflux,&bgmedian,&bgsigma,
				&atot,&abad,NULL);
		 }
		else
		 {	bgarea=M_PI*(2.0*ap.ra*ap.da+ap.da*ap.da);
			bgflux=xpp->sky*bgarea;
			bgmedian=xpp->sky;
			bgsigma=0.0;
			atot=(int)bgarea;
			abad=0;
			r=0;
		 }

		background=bgmedian;

		if ( r )
		 {	wpf=&pf[j];
			wpf->flag=MASK_NOBACKGROUND;
			wpf->rtot=0;
			wpf->rbad=0;
			wpf->rign=0;
			wpf->atot=0;
			wpf->abad=0;
			continue;
		 }

		wpf=&pf[j];
		wpf->bgmedian=bgmedian,
		wpf->bgsigma=bgsigma;
		wpf->bgflux=bgflux;
		wpf->bgarea=bgarea;
		wpf->atot=atot;
		wpf->abad=abad;

		/*
		if ( xpp->is_disjoint_apertures )
			ringmask_subtract(ringmasks[j],sx,sy,cx,cy,ap.ra);
		*/

		wpf=&pf[j];
		wpf->flag=0;

		if ( bqc != NULL )
		 {	/*
			if ( xpp->is_disjoint_apertures || xpp->disjoint_radius>0.0 )
				r=aperture_photometry_flux_biquad(bqc,mask,sx,sy,cx,cy,ap.r0,&area,&flux,&rtot,&rbad,&rign,xpp->maskignore,ringmasks[j]);
			else
			*/
			r=aperture_photometry_flux_biquad(bqc,mask,sx,sy,cx,cy,ap.r0,&area,&flux,&rtot,&rbad,&rign,xpp->maskignore,NULL);
		 }
		else if ( ! (weightusage & USE_WEIGHT_WEIGHTED) )
		 {	/*
			if ( xpp->is_disjoint_apertures || xpp->disjoint_radius>0.0 )
				r=aperture_photometry_flux(img->data,mask,sx,sy,cx,cy,ap.r0,&area,&flux,&rtot,&rbad,&rign,xpp->maskignore,xpp->subpixeldata,xpp->subg,ringmasks[j]);
			else
			*/
			r=aperture_photometry_flux(img->data,mask,sx,sy,cx,cy,ap.r0,&area,&flux,&out,background,&rtot,&rbad,&rign,xpp->maskignore,xpp->subpixeldata,xpp->subg,NULL);
		 }
		else
		 {	double	r0;
			double	ccx,ccy,cx0,cy0,bg;
			int	ix0,iy0,i,j;
			int	ix,iy,wx,wy;
			weight	*ww;

			cx0=floor(cx);
			cy0=floor(cy);
			ccx=(double)grid*((cx-cx0)+(double)hsize);
			ccy=(double)grid*((cy-cy0)+(double)hsize);
			r0 =(double)grid*ap.r0;

			ww=weight_get_closest(wl,cx,cy);

			ix0=(int)cx0-hsize;
			iy0=(int)cy0-hsize;

			for ( i=0 ; i<wfsize/grid ; i++ )
			 {	for ( j=0 ; j<wfsize/grid ; j++ )
				 {	aphwsum[i][j]=0.0;		}
			 }
			for ( i=0 ; i<wfsize ; i++ )
			 {	for ( j=0 ; j<wfsize ; j++ )
				 {	aphwsum[i/grid][j/grid]+=ww->iarr[i][j];	}
			 }

			bg=bgmedian*igrid2;
		
			for ( i=0 ; i<fsize ; i++ )
			 {	iy=iy0+i/grid;
				for ( j=0 ; j<fsize ; j++ )
				 {	ix=ix0+j/grid;
					if ( ix<0 || iy<0 || ix>=sx || iy>=sy )
					 {	aphmask[i][j] = MASK_OUTER;
						aphdata[i][j] = 0.0;
						continue;
					 }
					else if ( mask != NULL )
					 {	aphmask[i][j] = mask[iy][ix];	}
					else
					 {	aphmask[i][j] = 0;		}

					aphdata[i][j]=img->data[iy][ix]*igrid2-bg;

					wx=j+grid*(-hsize+wl->hsize);
					wy=i+grid*(-hsize+wl->hsize);
					if ( wx>=0 && wy>=0 && wx<wfsize && wy<wfsize )
					 {	aphdata[i][j] *= grid2*ww->iarr[wy][wx]/aphwsum[wy/grid][wx/grid];	}

					aphdata[i][j]+=bg;
			
				 }
			 }

			r=aperture_photometry_flux(aphdata,aphmask,fsize,fsize,ccx,ccy,r0,&area,&flux,&out,background,&rtot,&rbad,&rign,xpp->maskignore,xpp->subpixeldata,xpp->subg,NULL);
		 }
		wpf->rtot=rtot;
		wpf->rbad=rbad;
		wpf->rign=rign;
		wpf->flag |= r;	

		flux -= bgmedian*area*igrid2;
		wpf->flux=flux;
		if ( flux<=0.0 )
		 { 	wpf->fluxerr=0.0;
		 }
		else
		 {	fluxerr=sqrt((0.0<gain?flux/gain:0.0)+area*(bgsigma*bgsigma)*(1.0+1.0/bgarea)),
			wpf->fluxerr=fluxerr;
		 }

		if ( 0.0 < out.fw && 0 < flux )
		 {	double	mx,my,mxx,mxy,myy,sxx,syy,sxy,ac0,w;
			mx=out.fwx/out.fw;
			my=out.fwy/out.fw;
			wpf->cntr_x=mx+cx;
			wpf->cntr_y=my+cy;
			mxx=out.fwxx/out.fw;
			mxy=out.fwxy/out.fw;
			myy=out.fwyy/out.fw;
			sxx=mxx-mx*mx;
			sxy=mxy-mx*my;
			syy=myy-my*my;
			ac0=area*bgsigma/flux;
			wpf->cntr_x_err=sqrt((0.0<gain?(sxx)/(flux*gain):0.0)+(ac0*ac0)/(4.0*M_PI));
			wpf->cntr_y_err=sqrt((0.0<gain?(syy)/(flux*gain):0.0)+(ac0*ac0)/(4.0*M_PI));
			w=sxx*syy-sxy*sxy;
			if ( w <= 0.0 )	
			 {	wpf->cntr_width=0.0;
				wpf->cntr_w_err=-1.0;
			 }
			else	
			 {	double	sig2;
				double	a,b,c,ws,wd,wk;
				a=(sxx+syy)/2;
				b=(sxx-syy)/2;
				c=sxy;
				ws=sqrt(0.5*(a+sqrt(a*a-b*b-c*c)));
				wd=b/(2*ws);
				wk=c/(2*ws);
				wpf->cntr_width=ws;
				wpf->cntr_w_d=wd;
				wpf->cntr_w_k=wk;
				sig2=sqrt(w);
				wpf->cntr_w_err=sqrt((0.0<gain?(2*sig2*sig2/(flux*gain)):0.0)+ac0*ac0*area/(8.0*M_PI));
			 }
		 }
		else
		 {	wpf->cntr_x=0.0;
			wpf->cntr_y=0.0;
			wpf->cntr_x_err=-1.0;
			wpf->cntr_y_err=-1.0;
			wpf->cntr_width=-1.0;
			wpf->cntr_w_err=-1.0;
		 }


		/*
		if ( xpp->is_disjoint_apertures )
			ringmask_coadd(ringmasks[j],sx,sy,cx,cy,ap.ra);
		*/

	 } /* for ( j=0 ; j<numap ; j++ ) */

	if ( ww != NULL )
		add_to_data_weight(img->data,mask,sx,sy,wl->hsize,wl->grid,ww,-1.0);

	if ( is_calc_opt_apert && numap==1 && ( ! pf[0].flag ) )
	 {	double	r0opt;

		if ( 0.0<gain )
			r0opt=optimal_aperture(gain,sigma,pf[0].bgsigma,pf[0].bgarea,pf[0].flux);
		else
			r0opt=0.0;

		ps[i].optimal.r0=r0opt;
		ps[i].optimal.ra=pf[0].ag.ra;
		ps[i].optimal.da=pf[0].ag.da;
	 }
	else
	 {	ps[i].optimal.r0=0.0;
		ps[i].optimal.ra=0.0;
		ps[i].optimal.da=0.0;
	 }
  }

 if ( ringmasks != NULL )
  {	for ( j=numap-1 ; j>=0 ; j-- )
	 {	if ( ringmasks[j] != NULL )
			tensor_free(ringmasks[j]);	
	 }
	free(ringmasks);
	ringmasks=NULL;
  }

 if ( aphmask != NULL )	tensor_free(aphmask);
 if ( aphdata != NULL )	tensor_free(aphdata);
 if ( bqc != NULL )	tensor_free(bqc);

 if ( weightusage & USE_WEIGHT_SUBTRACTED )
 	add_to_image_weights(img,mask,wl,+1.0);

 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int do_subtracted_photometry(fitsimage *img,char **mask,photstar *ps,int np,
	apgeom *inaps,int numap,
	spatialgain *sg,xphotpar *xpp,kernellist *klist)
{
 int		i,j,k,jp,bx,by,sx,sy,r,atot,abad,rtot,rbad,rign,flag;
 apphotpar	ap;
 apgeom		*aps;
 int		n;
 photflux	*pf,*wpf;
 double		bgarea,bgflux,bgmedian,bgsigma,gain;
 double		area,flux,fluxerr,cx,cy,kx,ky,flux0,sflux,kflux;
 double		**bqc;
 photflux	**fluxarr;
 double		**kernarr;
 int		hsize,fsize;

 if ( img==NULL || img->data==NULL )	return(1);
 sx=img->sx,sy=img->sy;
 if ( sx<=0 || sy<=0 )			return(1);

 if ( xpp->use_biquad )
  {	bqc=tensor_alloc_2d(double,sx*2+1,sy*2+1);
	if ( bqc==NULL )	return(-1);
	biquad_coeff(img->data,sx,sy,bqc,NULL);
  }
 else
	bqc=NULL;

 memcpy(&ap.bgm,&xpp->bgm,sizeof(bgmode));

 if ( klist != NULL )
  {	hsize=0;
	for ( k=0 ; k<klist->nkernel ; k++ )
	 {	if ( klist->kernels[k].hsize>hsize )
			hsize=klist->kernels[k].hsize;
	 }
  }
 else
	hsize=0;

 fsize=2*hsize+1;
 fluxarr=(photflux **)tensor_alloc_2d(photflux,fsize,fsize);
 kernarr=(double   **)tensor_alloc_2d(double  ,fsize,fsize);

 for ( i=0 ; i<np ; i++ )
  {	if ( inaps != NULL )	aps=inaps,n=numap;
	else			aps=ps[i].inaps,n=ps[i].n;

	ps[i].n=n;
	pf=(photflux *)malloc(sizeof(photflux)*n);
	ps[i].fluxes=pf;

	cx=ps[i].x ,cy=ps[i].y;
	kx=ps[i].x ,ky=ps[i].y;

	gain=eval_2d_poly(cx,cy,sg->order,sg->coeff,0.5*(double)sx,0.5*(double)sy,0.5*(double)sx);
	if ( sg->vmin>0 && gain<sg->vmin )	gain=sg->vmin;
	ap.gain=gain;

	for ( j=0 ; j<n ; j++ )
	 {	
		wpf=&pf[j];

		wpf->flag=0;
		wpf->atot=0;
		wpf->abad=0;
		wpf->rtot=0;	
		wpf->rbad=0;
		wpf->rign=0;

		ap.r0=aps[j].r0,
		ap.ra=aps[j].ra,
		ap.da=aps[j].da;
		for ( k=0,jp=-1 ; k<j && jp<0 ; k++ )
		 {	if ( pf[k].ag.ra==ap.ra && pf[k].ag.da==ap.da )
			 {	jp=k;break;		}
		 }
		if ( jp>=0 )
		 {	bgarea  =pf[jp].bgarea,
			bgflux  =pf[jp].bgflux,
			bgmedian=pf[jp].bgmedian,
			bgsigma =pf[jp].bgsigma;
			atot    =pf[jp].atot;
			abad    =pf[jp].abad;
			r=0;
		 }
		else
			r=aperture_photometry_back(img->data,mask,sx,sy,cx,cy,&ap,&bgarea,&bgflux,&bgmedian,&bgsigma,&atot,&abad,NULL);


		if ( r )
		 {	wpf->flag=MASK_NOBACKGROUND;
			continue;
		 }

		flag=0;

		rtot=0;	
		rbad=0;
		rign=0;

		for ( by=-hsize ; by<=hsize ; by++ )
		 {	for ( bx=-hsize ; bx<=hsize ; bx++ )
			 {	double	a,f;
				int	rt,rb,ri;
				if ( xpp->use_biquad )
					r=aperture_photometry_flux_biquad(bqc,mask,sx,sy,cx+bx,cy+by,ap.r0,&a,&f,&rt,&rb,&ri,xpp->maskignore,NULL);
				else
					r=aperture_photometry_flux(img->data,mask,sx,sy,cx+bx,cy+by,ap.r0,&a,&f,NULL,0.0,&rt,&rb,&ri,xpp->maskignore,NULL,0,NULL);
				if ( rt>rtot )	rtot=rt;
				if ( rb>rbad )	rbad=rb;
				if ( ri>rign )	rign=ri;

				fluxarr[hsize+by][hsize+bx].flux=f-bgmedian*a;
				fluxarr[hsize+by][hsize+bx].flag=r;
				fluxarr[hsize+by][hsize+bx].bgarea=a;

				flag |= r;
			 }
		 }

		if ( flag )
		 {	wpf->flag |= flag;
			wpf->flux=0.0;
			wpf->fluxerr=0.0;
			continue;
		 }

		if ( klist == NULL )
		 {	flux=fluxarr[0][0].flux;
			kernarr[0][0]=1.0;
			kflux=1.0;
		 }

		else
		 {	int	c;
			kernel	*k;
			double	kcoeff;

			for ( by=0 ; by<fsize ; by++ )
			 {	for ( bx=0 ; bx<fsize ; bx++ )
				 {	kernarr[by][bx]=0.0;		}
			 }

			for ( c=0,k=klist->kernels ; c<klist->nkernel ; c++,k++ )
			 {	if ( k->type==KERNEL_BACKGROUND )
					continue;
				kcoeff=eval_2d_poly(kx,ky,k->order,k->coeff,klist->ox,klist->oy,klist->scale);

				if ( k->type==KERNEL_IDENTITY )
					kernarr[hsize][hsize]+=kcoeff;
				else
				 {	for ( by=-k->hsize ; by<=+k->hsize ; by++ )
					 {  for ( bx=-k->hsize ; bx<=+k->hsize ; bx++ )
					     {	kernarr[hsize+by][hsize+bx]+=kcoeff*k->image[k->hsize+by][k->hsize+bx];	}
					 }
				 }				
			 }

			kflux=0.0;
			for ( by=0 ; by<fsize ; by++ )
			 {	for ( bx=0 ; bx<fsize ; bx++ )
				 {	kflux += kernarr[by][bx];	}
			 }
		 }

		sflux=0.0;
		area=0.0;
		for ( by=0 ; by<fsize ; by++ )
		 {	for ( bx=0 ; bx<fsize ; bx++ )
			 {	sflux+=kernarr[by][bx]*fluxarr[by][bx].flux;
				area +=kernarr[by][bx]*fluxarr[by][bx].bgarea;
			 }
		 }
		sflux/=(kflux*kflux);
		area /=kflux;

		wpf->atot=atot;
		wpf->abad=abad;
		wpf->rtot=rtot;	
		wpf->rbad=rbad;
		wpf->rign=rign;

		wpf->bgmedian=bgmedian,
		wpf->bgsigma=bgsigma;
		wpf->bgflux=bgflux;
		wpf->bgarea=bgarea;

		wpf->ag.r0=ap.r0;
		wpf->ag.ra=ap.ra;
		wpf->ag.da=ap.da;

		if ( ps[i].rfflux[j].flag )
		 {	wpf->flag |= ps[i].rfflux[j].flag;
			wpf->flux=0.0;
			wpf->fluxerr=0.0;
			continue;
		 }

		flux0=ps[i].rfflux[j].flux;

		flux = flux0 + sflux;

		if ( flux<=0.0 || flux0<=0.0 )
		 {	wpf->flux=0.0,
			wpf->fluxerr=0.0;
		 }
		else
		 {	fluxerr=sqrt((0.0<gain?flux/gain:0.0)+area*(bgsigma*bgsigma)*(1.0+1.0/bgarea));
			wpf->flux=flux,
			wpf->fluxerr=fluxerr;
		 }

		wpf->flag=0;
	 }

	ps[i].n=n;
  }

 if ( bqc != NULL )	tensor_free(bqc);

 return(0);
}

/*****************************************************************************/

int do_magnitude_fit(photstar *stars,int nstar,
magfitparams *mfp,int sx,int sy,magfitstat *mfs)
{
 int		a,nvar,i,j,n,maxorder,k,l,v,iiter;
 double		**amatrix,*bvector,*fvars,*monoms,yvar;
 double		ox,oy,scale,dmag,cmul,w;
 photstar	*s;
 photflux	*fi,*fr;
 double		*dmags,*wghts,*smags;
 
 nvar=0;
 maxorder=0;
 for ( j=0 ; j<=mfp->norder ; j++ )
  {	nvar+=(mfp->orders[j]+1)*(mfp->orders[j]+2)/2;
	if ( maxorder<mfp->orders[j] )
		maxorder=mfp->orders[j];
  }

 /* n: total (maximal) number of apertures */
 n=0;
 for ( i=0 ; i<nstar ; i++ )
  {	if ( n < stars[i].n )	n=stars[i].n;		}

 amatrix=matrix_alloc(nvar);
 bvector=vector_alloc(nvar);
 fvars  =vector_alloc(nvar);
 monoms =vector_alloc((maxorder+1)*(maxorder+2)/2);

 ox=(double)sx/2.0;
 oy=(double)sy/2.0;
 scale=ox;

 dmags=(double *)malloc(sizeof(double)*nstar);
 smags=(double *)malloc(sizeof(double)*nstar);
 wghts=(double *)malloc(sizeof(double)*nstar);

 if ( mfs != NULL )
  {	mfs->nstar=nstar;
	mfs->naperture=n;
	mfs->ninit=(int *)malloc(sizeof(int)*n);
	mfs->nrejs=(int *)malloc(sizeof(int)*n);
	for ( a=0 ; a<n ; a++ )
	 {	mfs->ninit[a]=0;
		mfs->nrejs[a]=0;
	 }
  }

 for ( a=0 ; a<n ; a++ )
  {
	for ( i=0 ; i<nstar ; i++ )
	 {	wghts[i]=-1.0;
		dmags[i]=0.0;

		s=&stars[i];
		if ( s->n < a || s->fluxes==NULL || s->rfflux==NULL )
			continue;
		fi=&s->fluxes[a];
		fr=&s->rfflux[a];

		if ( fi->flag || fr->flag || fi->flux<=0.0 || fr->flux<=0.0 )
			continue;

		dmag=log(fi->flux/fr->flux);
		dmags[i]=dmag;
		w=sqrt(fr->flux);
		wghts[i]=w;

		mfs->ninit[a]++;
	 }

	for ( iiter=0 ; iiter<=mfp->niter ; iiter++ )

	 {	for ( k=0 ; k<nvar ; k++ )
		 {	for ( l=0 ; l<nvar ; l++ )
			 {	amatrix[k][l]=0.0;		}
			bvector[k]=0.0;
		 }

		for ( i=0 ; i<nstar ; i++ )
		 {	
			if ( wghts[i] <= 0.0 )
				continue;

			s=&stars[i];
			fi=&s->fluxes[a];
			fr=&s->rfflux[a];

			dmag=dmags[i];

			v=0;
			cmul=1.0;

			w=wghts[i];	/* weight */
	
			for ( j=0 ; j<=mfp->norder ; j++ )
			 {	eval_2d_monoms(s->x,s->y,mfp->orders[j],fvars+v,ox,oy,scale);
				l=(mfp->orders[j]+1)*(mfp->orders[j]+2)/2;
				for ( k=0 ; k<l ; k++ )
				 {	fvars[v] *= cmul;
					v++;
				 }
				cmul=cmul*s->ref_col;
			 }
			yvar=dmag;
	
			for ( k=0 ; k<nvar ; k++ )
			 {	for ( l=0 ; l<nvar ; l++ )
				 {	amatrix[k][l] += w*fvars[k]*fvars[l];	}
				bvector[k] += w*fvars[k]*yvar;
			 }
		 }

		solve_gauss(amatrix,bvector,nvar);

		if ( iiter < mfp->niter )
		 {	double	s0,s1,s2,d;

			for ( i=0 ; i<nstar ; i++ )
			 {	s=&stars[i];
				if ( wghts[i] <= 0.0 )
					continue;
				fi=&s->fluxes[a];
				fr=&s->rfflux[a];
	
				v=0;
				cmul=1.0;
				for ( j=0 ; j<=mfp->norder ; j++ )
				 {	eval_2d_monoms(s->x,s->y,mfp->orders[j],fvars+v,ox,oy,scale);
					l=(mfp->orders[j]+1)*(mfp->orders[j]+2)/2;
					for ( k=0 ; k<l ; k++ )
					 {	fvars[v] *= cmul;
						v++;
					 }
					cmul=cmul*s->ref_col;
				 }
	
				dmag=0.0;
				for ( k=0 ; k<nvar ; k++ )
				 {	dmag += fvars[k]*bvector[k];		}
				smags[i]=dmag;
			 }

			s0=s1=s2=0.0;			
			for ( i=0 ; i<nstar ; i++ )
			 {	if ( wghts[i] <= 0.0 )
					continue;
				d=dmags[i]-smags[i];		
				w=1.0;
				s0+=w;
				s1+=w*d;
				s2+=w*d*d;
			 }
			s1/=s0;
			s2/=s0;
			s2=s2-s1*s1;
			if ( s2<=0.0 )	s2=0.0;
			else		s2=sqrt(s2);

			for ( i=0 ; i<nstar ; i++ )
		 	 {	if ( wghts[i] <= 0.0 )
					continue;
				d=dmags[i]-smags[i];
				if ( fabs(d) >= s2*mfp->sigma )
				 {	wghts[i]=-1.0;
					mfs->nrejs[a]++;
				 }
			 }
		 }
	 }

	for ( i=0 ; i<nstar ; i++ )
	 {	s=&stars[i];
		if ( s->n < a || s->fluxes==NULL || s->rfflux==NULL )
			continue;
		fi=&s->fluxes[a];
		fr=&s->rfflux[a];

		if ( fi->flux<=0.0 || fr->flux<=0.0 )
			continue;

		v=0;
		cmul=1.0;
		for ( j=0 ; j<=mfp->norder ; j++ )
		 {	eval_2d_monoms(s->x,s->y,mfp->orders[j],fvars+v,ox,oy,scale);
			l=(mfp->orders[j]+1)*(mfp->orders[j]+2)/2;
			for ( k=0 ; k<l ; k++ )
			 {	fvars[v] *= cmul;
				v++;
			 }
			cmul=cmul*s->ref_col;
		 }

		dmag=dmags[i];
		for ( k=0 ; k<nvar ; k++ )
		 {	dmag -= fvars[k]*bvector[k];		}

		fi->flux=fr->flux*exp(dmag);
	
	 }

  }

 free(wghts);
 free(smags);
 free(dmags);

 vector_free(monoms);
 vector_free(fvars);
 vector_free(bvector);
 matrix_free(amatrix);

 return(0);
}

/*****************************************************************************/

int fprint_fiphot_usage(FILE *fw)
{
 fprintf(fw,
"Usage:\tfiphot [-h|--help] [-V|-W|--verbose] [-C|--comment]\n"
"\t[[-i|--input] <input>|-] [-o|--output <out.phot>|-]\n"
"\t-L|--input-list <input-list>|- [--output-list <output-list>|-]\n"
"\t{<image>[<F>]|-r <reference>[<F>] -s <sub1> [<sub2>...] [--frame <F>]}\n");
 fprintf(fw,
"General parameters:\n"
"\t[-a|--aperture[s] r01[:ra1[:da1]][,r02[:ra2[:da2]][,...]] ]\n"
"\t[--sky-fit {mean|mode|median,iterations=n,lower=l,upper=u,sigma=s}|\n"
"\t           |force=<sky>]\n"
"\t[-j|--disjoint-rings | --disjoint-apertures | -x|--disjoint-radius <r>]\n"
"\t[-z|--zoom <z>] [-k|--spline] [--aperture-mask-ignore <masks>]\n"
"\t[-g|--gain <gain>[,<polynomial-coeffs>]|auto] [--mag-flux <mag>,<flux>]\n"
"\t[-u <subpixelfile> [--normalize]]\n");
 fprintf(fw,
"Input and output format parameters:\n"
"\t[--col-xy <>,<>] [--col-ap <>[:<>[:<>]]...] [--col-id <>]\n"
"\t[--col-mag <> [--col-magerr <>] [--col-color <>]]\n"
"\t[-F|--format [XYIS-],[MmBbFfAaEeSsWwDK-]] [--nan-string <INDEF-string>]\n"
"\t[--serial <serial-string>]\n");
 fprintf(fw,
"Post-photometry magnitude transformation parameters (a.k.a. 'magfit'):\n"
"\t[--magfit orders=<c0>[:<c1>[:<c2>[:<c3>]]],niter=<n>,sigma=<s>]\n");
 fprintf(fw,
"Parameters of optimal aperture determination:\n"
"\t[-d|--skysigma <skysigma>] [-f|--fwhm <FWHM>]\n");
 fprintf(fw,
"Weight stamps:\n"
"\t[--input-weight <weight.fits>]\n"
"\t[--weight-usage [weighting],[subtraction],confidence=<confdist>]\n");
 return(0);
}

longhelp_entry fiphot_long_help[]=
{
 LONGHELP_OPTIONS,

 { "General options:", NULL },
 { "-h, --help",
	"Gives general summary about the command line options." },
 { "--long-help",
	"Gives a detailed list of command line options." },
 { "--version",
	"Gives some version information about the program." },
 { "-i, --input <image file>",
	"Name of the input FITS image file." },

 { "Simple aperture photometry:", NULL },
 { "-L, --input-list <input coordinate list file>",
	"Name of the input coordinate list file." },
 { "--col-xy <colx>,<coly>",
	"Column indices for centroid coordinates. The coordinates read from "
	"this file follows the native coordinating scheme (which is not "
	"the same as e.g. in IRAF), namely the lower-left corner of the "
	"lower-left pixel has the coordinate of (0,0) while the center of the "
	"lower-left pixel has the coordinate of (0.5,0.5). Programs like IRAF "
	"use the coordinate (1,1) for the center of the lower-left pixel. "
	"This coordinating princible can be overwritten by the --coord=* "
	"option." },
 { "--col-ap <A1>,<A2>,...", 
	"Column indices for various apertures. In each such column, there "
	"must be three colon-separated number, for the radius of the aperture, "
	"inner radius for the background annulus and the thickness of the "
	"annulus, respectively. This option is not mandatory, all of the "
	"objects can be measured with the same set of apertures, see also "
	"-a|--apertures for more details. " },
 { "--col-id <identifier column>",
	"Column index for object identifier. " },
 { "--col-mag, --col-magnitude <magnitude column>",
	"Column index for reference magnitude. " },
 { "--col-col, --col-color <color column>",
	"Column index for photometric color." },
 { "--col-err, --col-error <magnitude error column>",
	"Column index for magnitude uncertainty." },
 { "--coord=native, iraf, zero",
	"Use the native or IRAF-like coordinates (see --col-xy for further "
	"details), or the coordinate convention where the center of the "
	"lower-left pixel has a coordinate of (0,0)" },
 { "-z, --zoom <zoom level>",
	"Mutiply both the input centroid coordinates and aperture/annulus"
	"radii by the given integer factor. " },

 { "--serial <serial>",
	"Serial identifier for the whole photometry procedure. Can be any "
	"arbitrary string and used only in the formatted output (see "
	"-F|--format for more details)." },
 { "-F, --format, --format-output <objecttags>,<photometrytags>",
	"List of output format tags. The formatted (user-friendly) output "
	"photometry contains a few columns containing the data related to "
	"the object, which followed by a the per-aperture photometry "
	"results. See ``Format tags'' for the list of format tags used here." },
 { "--nan <nan-string>",
	"String which is used to denote bad photometry. By default, objects "
	"on which photometry cannot be performed (due to various reasons, e.g. "
	"the object is off from the image, the background level cannot be "
	"determined or there are bad pixels in the aperture itself), "
	"are marked by a simple dash ('-') in the output file. " },

 { "-M, --input-mask <image file>",
	"Input mask file to co-add to the mask of the input image. Useful for "
	"marking pixels to be ignored from the photometry process beyond the "
	"ones which are previously marked in the input image." },
 { "-a, --aperture, --apertures <list of apertures>",
	__extension__
	"List of apertures to be involved in the photometry. Each aperture "
	"is defined by three numbers: the radius of the aperture, and the "
	"inner radius and ``thickness'' of the annulus used for sky background "
	"estimation. The aperture specifications mustbe spearated by "
	"commas while these three numbers must be separated by colons. "
	"E.g. to perform aperture photometry on a series of apertures with "
	"a radius of 1.5, 2.0 and 2.5 pixels, where all of the annuli "
	"have an inner and outer radius of 6.5 and 12 pixels (i.e. the "
	"thickness is 5.5 pixels), one should write "
	"1.5:6.5:5.5,2.0:6.5:5.5,2.5:6.5:5.5 as an argument for this option." },
 { "-g, --gain, --gain-poly <gain polynomial>",
	__extension__
	"The polynomial describing the gain level throughout the image. "
	"Altough during the image readout process, the electron <-> ADU "
	"conversion ratio has a fixed value, during the calibration process "
	"when the vignetting is strong, the ADU levels may substantially "
	"change. The comma-separated numbers in the gain polynomial should "
	"denote the coefficients for the monomials 1, x, y, 1/2x^2, xy, ... "
	"(the standard order for 2 dimensional polynomial coefficients), "
	"where x and y are the normalized coordinates (i.e. zero at the "
	"center of the image, x = +/- 1 at the right/left edge of the image" 
	"and y is scaled appropriately keeping the aspect ratio). Note that"
	"the number of coefficients should be 1, 3, 6, 10, ... and so on, "
	"for zeroth, first, second, third... order variations, respectively." 
	"Specifying zero or negative gain will imply an ``infinitely large'' "
	"gain, thus data are treated as being affected only by instrumental "
	"noise and lack intrinsic photon noise. " },
 { "--gain-vmin <minimal gain>",
	"The minimal value for the gain level. If the polynomial describing "
	"the spatial gain variations is evaluated on the regular image domain "
	"and yielded a smaller value than this given number, this "
	"will be used as gain level. In certain optical systems, the vignetting "
	"can be well described by second-order polynomial coefficients except "
	"at the very corners of the image. In such a situation this minimal gain "
	"is quite useful. " },

 { "--mag-flux <mag>,<flux>",
	"Magnitude - flux conversion level. The specified magnitude will "
	"be equivalent to the specified flux level. " },
 { "--sky-fit <sky fitting parameters>",
	"This argument is followed by a set of parameters forthe "
	"sky (i.e. background level) fitting algorithm. See ``Sky fitting "
	"parameters'' below for more details. " },
 { "--aperture-mask-ignore <list of masks>",
	"This switch is followed by a space separated list of standard "
	"masks which should be ignored if such a pixel is marked so in the "
	"aperutre. In practice, it might be useful to put saturated objects "
	"into the set of ``good'' stars. " },

 { "-o, --output <output photometry file>",
	"Name of the output file containing the results of the aperture "
	"photometry. The format and content of this file can be arbitrarily "
	"set, see -F|--format." },
 { "--output-raw-photometry <output raw photometry>",
	__extension__
	"Name of the output file containing the all of the low-level photometric "
	"information in a fixed format. From this file, one can derive all of the "
	"quantities which are written to the ``normal'' output photometry "
	"file. The main purpose of this file is to be an input for image "
	"subtraction based photometry, i.e. the photometric information for "
	"the reference image is supposed to be stored in this format and "
	"the successive calls of `fiphot` on the subtracted residual images "
	"read this information in order to derive the final photometric "
	"information. See the subsection ``Photometry on subtracted/convolved "
	"images'' about more details on the image subtraction "
	"based photometry." },

 { "Note that the literal \"auto\" argument can also be used after the "
   "-g|--gain switch. In this case, `fiphot` tries to figure out the "
   "gain polynomial from the GAINPOLY and GAIN keywords (in this order) "
   "as well as the minimal value for the gain from the GAINVMIN keyword.",
   NULL }, { "", NULL },

 { "Format tags for generic object data:", NULL },
 { "I",	"identifier" },
 { "S", "serial identifier (set by --serial)" },
 { "X", "X coordinate of the centroid" },
 { "Y", "Y coordinate of the centroid" },
 { "-", "empty column" },
 
 { "Format tags for photometric data:", NULL },
 { "M", "magnitude" },
 { "m", "uncertainty in the magnitude" },
 { "B", "background level" },
 { "b", "background scatter" },
 { "F", "flux" },
 { "f", "flux uncertainty" },
 { "A", "same as ``F''" },
 { "a", "same as ``f''" },
 { "E", "flux in electrons" },
 { "e", "flux uncertainty in electrons" },
 { "X", "X coordinate of the fitted centroid" },
 { "x", "X coordinate uncertainty" },
 { "Y", "Y coordinate of the fitted centroid" },
 { "y", "Y coordinate uncertainty" },
 { "W", "Statistical profile size (S), in pixels" },
 { "D", "Statistical profile deviation parameter (D), in pixels" },
 { "K", "Statistical profile deviation parameter (K), in pixels" },
 { "w", "Uncertainty of statistical profile size" },
 { "-", "empty column" },
 
 { "Fine tuning of aperture photometry:", NULL },
 { "-j, --disjoint-annuli, --disjoint-rings",
	"During the bacground determination on the aperture annuli, omit "
	"the pixels which belongs to the annuli of other centriods. On very "
	"dense fields this might make the aperture photometry impossible since "
	"for some or many of the target centroids, the bacground level cannot "
	"be derived due to the lack of sufficient number of background pixels. " },
 { "-p, --disjoint-apertures",
	"During the bacground determination on the aperture annuli, omit "
	"the pixels which belongs to the apertures of other centriods. " },
 { "-x, --disjoint-radius <radius>",
	"During the bacground determination on the aperture annuli, omit "
	"the pixels which are closer to the other centroids than the "
	"specified radius." },
 { "-k, --spline",
	"Use a flux-conserving biquadratic spline interpolation surface "
	"for photometry purposes. This surface yields some sort of weighting "
	"within each pixel due to the properties of the spline interpolation." },
 { "-m, --magfit orders=<c0>[:<c1>[:<c2>[:<c3>]]],iterations=<d>,sigma=<s>",
	__extension__
	"Perform a magnitude transformation after the photometry. Currently "
	"impleneted only in the case of image subtraction photometry and "
	"always use the reference photometry (see --input-raw-photometry) "
	"as a reference for magnitude transformation too. This command line "
	"option must be followed by the parameters of the magnitude fit, "
	"namely the list the of spatial polynomial orders for the "
	"subsequent color orders (colon-separated list), and the number of "
	"(optional) outlier rejection iterations and its limit in standard "
	"deviation (sigma) units." },

 { "Photometry on subtracted/convolved images:", NULL },
 { "-P, --input-raw-photometry <input reference raw photometry>",
	"Name of the file containing the coordinate lists and the raw "
	"photometric information for the reference image." },
 { "-K, --input-kernel <input file with kernel solution>",
	"The kernel solution which resulted during the creation of the "
	"convolved and/or subtracted image. This information is also required "
	"for the proper self-consistent aperture photometry on subtracted images. "
	"Omitting this file will result an assumption for identical "
	"convolution transformation, which only appropriate if the subtracted "
	"image is created by a literal arithmetic subtraction (and not "
	"convolution based subtraction)." },

 { "Optimal aperture determination:", NULL },
 { "--output-list <output coordinate list file>",
	"The name of the output centroid list file, in which the radius "
	"of the optimal aperture is also stored. The optimal aperture is "
	"derived from the object flux itself, the gain level, the "
	"background noise and the FWHM of the objects. "
	"The background noise can be specified using the -d|--sky-noise "
	"argument (see there). The FWHM is given by the option -f|--fwhm." },
 { "-d, --skynoise <noise>",
	"Sky (bacground level) noise level in ADUs." },
 { "-f, --fwhm <FWHM>",
	"Full width at half magnitude (FWHM) for the stellar objects." },

 { "Sky fitting parameters:", NULL },
 { "mean",
	"Use the mean value of the pixels in the annulus as a sky value." },
 { "median",
	"Use the median value of the pixels in the annulus as a sky value." },
 { "mode",
	"Use the mode of the pixels in the annulus as a sky value." },
 { "iterations=<iterations>",
	"Do the specified number of iterations in order to reject outlier pixels." },
 { "lower, upper, sigma=<sigma level>",
	"Lower, upper and generic (symmetric) outlier level in the units of "
	"standard deviations." },
 { "force=<level>",
	"Use the forced constant level for sky background, with zero nominal "
	"scatter." },

 { NULL, NULL }

};

int fprint_fiphot_long_help(FILE *fw)
{
 fprintf(fw,
"Usage:\tfiphot [options] [<input>] [-o|--output <output>]\n"
"This program performs aperture photometry on normal or subtracted/convolved "
"images.\n");
 fprintf(fw,"\n");

 longhelp_fprint(fw,fiphot_long_help,0,-1);

 fprintf(fw,"\n");
 fprintf(fw,"Report bugs to <%s>\n",FI_MAINT_EMAIL);
 
 return(0);
}

/*****************************************************************************/

int parse_gainpoly_string(char *gainstr,spatialgain *sg)
{
 char	*gs,**tokens;
 int	i,p,o,n;

 gs=strdup(gainstr);
 tokens=tokenize_char_dyn(gs,',');
 for ( n=0 ; tokens[n] != NULL ; )	n++;
 sg->coeff=NULL;
 for ( o=0,p=0 ; p+o+1 <= n ; )
  {	sg->coeff=(double *)realloc(sg->coeff,sizeof(double)*(p+o+1));
	for ( i=0 ; i<=o ; i++,p++ )
	 {	sscanf(tokens[p],"%lg",&sg->coeff[p]);		}
	o++;
  }
 sg->order=o-1;

 free(tokens);
 free(gs);

 return(0);
}

int get_gain_info_from_header(spatialgain *sg,fitsheaderset *header)
{
 fitsheader	*fh;
 double		gain,vmin;

 if ( (fh=fits_headerset_get_header(header,"GAINPOLY",0)) != NULL && fh->vtype==FITS_VSTR )
  {	parse_gainpoly_string(fh->vstr,sg);	  }
 else if ( ! fits_headerset_get_as_double(header,"GAIN",&gain,!0) )
  {	sg->order=0;
 	sg->coeff=(double *)malloc(sizeof(double));
	sg->coeff[0]=gain;
  }

 if ( ! fits_headerset_get_as_double(header,"GAINVMIN",&vmin,!0) )
  {	if ( vmin>0.0 )
		sg->vmin=vmin;
	else
		sg->vmin=0.0;
  }
 else
	sg->vmin=0.0;

 return(0);
}

/*****************************************************************************/

int main(int argc,char *argv[])
{
 fits		*img,*refimg,*subimg;

 char		**mask,**submask;
 FILE		*fr;
 char 		*imgname,*inlistname,*inweightname,
		*outlistname,*outname,*outphotname,*inphotname;
 char		*refname,*subname,*subkernelfile,**inmasklist,
		*appar,*apcolpar,*skyfitpar,*wusagepar,*magfitpar;
 char		*oformat,*ofxy,*ofph,*nanstring,*serialstring,*maskignorestr,
		*gainstr;
 int		frameno;
 kernellist	klist;
 double		skysigma,sigma,fwhm;
 magflux	mf0;
 int		i,is_help,basistype;
 xphotpar	xpp;
 photstar	*ps;
 int	  	np;
 colread	col;
 apgeom		*inaps;
 int		ninap,numap,zoom;
 weightlist	wl_data,*wl;
 int		weightusage;
 char		*subpixfile;
 int		is_normalize;
 magfitparams	mfp;
 magfitstat	mfs;
 spatialgain	sg;
 int		is_auto_gain=0;

 progbasename=strrchr(argv[0],'/');
 if ( progbasename != NULL )	progbasename++;
 else				progbasename=argv[0];

 imgname=inlistname=outlistname=outname=NULL;
 subname=subkernelfile=inweightname=NULL;
 refname=NULL;inmasklist=NULL;
 outphotname=inphotname=NULL;
 subpixfile=NULL;is_normalize=0;
 wusagepar=oformat=appar=apcolpar=NULL,
 serialstring=NULL,nanstring="-";
 magfitpar=NULL;

 is_verbose=is_comment=is_help=0;frameno=0;
 gainstr=NULL;
 fwhm=SIG_FWHM;

 skyfitpar=NULL;
 xpp.bgm.type=BGTYPE_MEDIAN;
 xpp.bgm.rejniter=2,
 xpp.bgm.rejlower=xpp.bgm.rejupper=2.0;
 xpp.maskignore=MASK_SATURATED|MASK_HOT;
 xpp.use_biquad=0;
 xpp.wconfdist=0.0;
 xpp.use_sky=0;
 xpp.sky=0.0;
 xpp.is_disjoint_rings=0;
 xpp.is_disjoint_apertures=0;
 xpp.disjoint_radius=-1.0;

 maskignorestr=NULL;

 mf0.magnitude=10.0,
 mf0.intensity=10000.0;

 col.colx=1,col.coly=2; 
 col.colap=NULL;
 col.colid=col.colmag=col.colcol=-1;
 basistype=0;
 
 sg.vmin=0.0;

 zoom=1;

 i=scanarg(argc,argv,SCANARG_ALLOW_FLAGS,
	"--version:%NS-1f%q",&is_help,
        "--version-short|--short-version:%NS-2f%q",&is_help,
	"-h|--help:%f%q",&is_help,
	"--long-help|--help-long:%SN2f%q",&is_help,
	"-i|--input:%s",&imgname,
	"--frame:%d",&frameno,
	"-o|--output:%s",&outname,
	"-L|--input-list:%s",&inlistname,
	"--output-list:%s",&outlistname,
	"--output-raw-photometry:%s",&outphotname,
	"-R|--input-raw-photometry:%s",&inphotname,
	"-M|--input-mask:%t",&inmasklist,
	"-W|--input-weight:%s",&inweightname,
	"-w|--weight-usage:%s",&wusagepar,
	"-r|--reference:%s",&refname,
	"-s|--input-subtracted:%s",&subname,
	"-K|--input-kernel:%s",&subkernelfile,
	"-d|--skysigma|--sky-noise|--skynoise:%g",&skysigma,
	"-g|--gain|--gain-poly|--gainpoly:%s",&gainstr,
	"--gain-vmin|--gainvmin:%g",&sg.vmin,
	"-f|--fwhm:%g",&fwhm,
	"--mag-flux:%g,%g",&mf0.magnitude,&mf0.intensity,
	"--sky-fit:%s",&skyfitpar,
	"--aperture-mask-ignore:%s",&maskignorestr,
	"-a|--aperture|--apertures:%s",&appar,
	"-j|--disjoint-rings|--disjoint-annuli:%f",&xpp.is_disjoint_rings,
	"-p|--disjoint-apertures:%f",&xpp.is_disjoint_apertures,
	"-x|--disjoint-radius:%g",&xpp.disjoint_radius,
	"-k|--spline:%f",&xpp.use_biquad,
	"-m|--magfit:%s",&magfitpar,
	"--col-xy:%Pd,%Pd",&col.colx,&col.coly,
	"--col-ap|--col-aperture|--col-apertures:%s",&apcolpar,
	"--col-id|--col-identifier:%Pd",&col.colid,
	"--col-mag|--col-magnitude:%Pd",&col.colmag,
	"--col-col|--col-color:%Pd",&col.colcol,
	"--col-err|--col-error|--col-merr|--col-magerr:%Pd",&col.colerr,
        "--coord=native:%SN0f",&basistype,
        "--coord=iraf:%SN1f",&basistype,
        "--coord=zero:%SN-1f",&basistype,
	"-z|--zoom:%Pd",&zoom,
	"-F|--format|--format-output:%s",&oformat,
	"--nan|--nan-format|--nan-string:%s",&nanstring,
	"--serial:%s",&serialstring,
	"-u|--input-subpixel:%s",&subpixfile,
	"--normalize:%f",&is_normalize,
	"--comment:%i",&is_comment,"(C):%i",&is_comment,
	"--verbose:%i",&is_verbose,"(V):%i",&is_verbose,
	"-:%w",&imgname,
	"-*|+*:%e",
	"*:%w",&imgname,
	NULL);

 if ( i )		
  {	fprint_error("invalid command line argument near '%s'",argv[i]);
	return(1);
  }
 else if ( is_help<0 )
  {	fprint_generic_version(stdout,argv[0],"fiphot",FI_PHOT_VERSION,is_help);
	return(0);
  }
 else if ( is_help>1 )
  {	fprint_fiphot_long_help(stdout);
	return(0);
  }
 else if ( is_help )
  {	fprint_fiphot_usage(stdout);
	return(0);
  }

 col.colx--,col.coly--;
 if ( col.colid>0 )	col.colid--;
 if ( col.colmag>0 )	col.colmag--;
 if ( col.colcol>0 )	col.colcol--;

 if ( gainstr != NULL && strcmp(gainstr,"auto")==0 )
  {	sg.order=-1;
	sg.coeff=NULL;
	sg.vmin=0.0;
	is_auto_gain=1;
  }
 else if ( gainstr != NULL )
  {	parse_gainpoly_string(gainstr,&sg);
	is_auto_gain=0;
  }
 else
  {	sg.order=0;
	sg.coeff=(double *)malloc(sizeof(double)*1);
	sg.coeff[0]=1.0;
	sg.vmin=0.0;
	is_auto_gain=0;
  }

 sigma=fwhm/SIG_FWHM;
 col.colap=create_col_ap_data(apcolpar);

 if ( subkernelfile != NULL )
  {	fr=fopenread(subkernelfile);
	if ( fr==NULL )	
	 {	fprint_error("unable to open kernel file '%s'",subkernelfile);
		return(1);
	 }
	else if ( kernel_info_read(fr,&klist) )
	 {	fprint_error("unable to parse kernel file contents");
		return(1);
	 }
	kernel_init_images(&klist);
	fcloseread(fr);
  }

 if ( magfitpar != NULL )
  {	char	**cmd;
	int	i,k;
	double	d;

	cmd=tokenize_char_dyn(magfitpar,',');

	mfp.norder=-1;
	mfp.niter=0;
	mfp.sigma=3.0;
	for ( i=0 ; cmd != NULL && cmd[i] != NULL ; i++ )
	 {	if ( (k=sscanf(cmd[i],"orders=%d:%d:%d:%d",
		&mfp.orders[0],&mfp.orders[1],
		&mfp.orders[2],&mfp.orders[3])) >= 1 )
			mfp.norder=k-1;
		else if ( sscanf(cmd[i],"niter=%d",&k)==1 )	
			mfp.niter=k;
		else if ( sscanf(cmd[i],"sigma=%lg",&d)==1 && d>0.0 )
			mfp.sigma=d;
		else
		 {	fprint_error("invalid magnitude fit parameter near '%s'.\n",cmd[i]);
			return(1);
		 }
	 }

  }
 else
  {	mfp.norder=-1;
	mfp.niter=0;
	mfp.sigma=0.0;
  }
 
 if ( inphotname != NULL && refname != NULL )
  {	fprint_error("invalid combination of command line options (both input raw photometry file and reference image are specified)");
	return(1);
  }
 if ( inphotname != NULL && inlistname != NULL )
  {	fprint_error("invalid combination of command line options (both input raw photometry file and input list file are specified)");
	return(1);
  }

 if ( skyfitpar != NULL )
  {	int	is_sigma;
	double	sigma=2.0;
	is_sigma=0;
	i=scanpar(skyfitpar,SCANPAR_DEFAULT,
		"mean:"  SNf(BGTYPE_MEAN)  ,&xpp.bgm.type,
		"median:"SNf(BGTYPE_MEDIAN),&xpp.bgm.type,
		"mode:"  SNf(BGTYPE_MODE)  ,&xpp.bgm.type,
		"iterations:%d",&xpp.bgm.rejniter,
		"lower:%g",&xpp.bgm.rejlower,
		"upper:%g",&xpp.bgm.rejupper,
		"sigma:%g%f",&sigma,&is_sigma,
		"force:%g%f",&xpp.sky,&xpp.use_sky,
		NULL);
	if ( is_sigma )	xpp.bgm.rejlower=xpp.bgm.rejupper=sigma;
  }
		

 if ( oformat==NULL )
  {	ofxy="IXY";
	ofph="MmBb";
  }
 else
  {	ofxy=output_format_xy(oformat);
	ofph=output_format_ph(oformat);
	if ( ofxy==NULL || ofph==NULL )	
	 {	fprint_error("invalid output format specification in '%s'",oformat);
		return(1);
	 }
  }

 if ( maskignorestr != NULL )
  {	xpp.maskignore=parse_mask_flags(maskignorestr);
	if ( xpp.maskignore<0 )
	 {	fprint_error("invalid mask specification in '%s'",maskignorestr);
		return(1);
	 }
  }

 /* silently ignore unexpected '--zoom' values: */
 if ( zoom<1 )	zoom=1;	

 i=create_input_ap_data(appar,&inaps,&ninap,zoom);
 if ( i )	
  {	fprint_error("invalid aperture specification in '%s'",appar);
	return(1);
  }

 if ( inlistname==NULL && inphotname==NULL )
  {	fprint_error("no coordinate list has been specified");
	return(1);
  }
 if ( inlistname != NULL )
  {	fr=fopenread(inlistname);
	if ( fr==NULL )	
	 {	fprint_error("unable to open input list file '%s'",inlistname);
		return(1);
	 }
	read_input_star_list(fr,&ps,&np,&col,zoom);
	fcloseread(fr);
	i=check_apertures(ps,np);
	if ( i<0 )
	 {	fprint_error("invalid aperture information in input list file");
		return(1);
	 }
	else if ( ninap==0 && i==0 )
	 {	fprint_error("no aperture information has been found in input list file");
		return(1);
	 }
	else if ( ninap> 0 && i> 0 )
	 {	fprint_error("ambiguous aperture specifications");
		return(1);
	 }
	else if ( ninap>0 )
		numap=ninap;
	else
		numap=i,inaps=NULL;
  }
 else	numap=0;

 if ( outlistname != NULL && numap != 1 )
  {	fprint_error("only one aperture should be specified if optimal aperture estimation is requested");
	return(1);
  }
	
 if ( subpixfile != NULL )
  {	FILE	*fr;
	int	i,j;
	fr=fopenread(subpixfile);
	if ( fr==NULL )	
	 {	fprint_error("unable to open subpixel data file");
		return(1);
	 }
	read_subpixel_file(fr,&xpp.subpixeldata,&xpp.subg);
	fcloseread(fr);
	if ( is_normalize )
		normalize_subpixeldata(xpp.subpixeldata,xpp.subg);
	for ( i=0 ; i<xpp.subg ; i++ )
	 {	for ( j=0 ; j<xpp.subg ; j++ )
		 {	xpp.subpixeldata[i][j]=1.0/xpp.subpixeldata[i][j];	}
	 }
  }
 else
  {	xpp.subpixeldata=NULL;
	xpp.subg=0;
  }

 if ( inweightname != NULL )
  {	fits    *img;	
	if ( (fr=fopenread(inweightname))==NULL )
	 {	fprint_error("unable to open input weight file '%s'.",inweightname);
		return(1);
	 }
	img=fits_read(fr);
	if ( img==NULL )
	 {	fprint_error("unable to interpret input as FITS data.");
		return(1);
	 }
	wl=&wl_data;
	if ( weight_parse_fits(img,wl) )
	 {	fprint_error("unable to interpret input as weight data.");
		return(1);
	 }
	fits_free(img);
	weight_sort(wl);
  }
 else
	wl=NULL;

 weightusage=0;
 if ( wusagepar != NULL )
  {	i=scanpar(wusagepar,SCANPAR_DEFAULT,
		"subtraction:%0f",&weightusage,
		"weighting:%1f",&weightusage,
		"confidence:%g",&xpp.wconfdist,
		NULL);
	if ( i )	
	 {	fprint_error("invalid weight usage specification in '%s'",wusagepar);
		return(1);
	 }
  }
	

/******************************************************************************
   One input image is given: simple aperture photometry for one image and/or
   do an estimation for the optimal aperture for each star.
******************************************************************************/
 if ( imgname != NULL )
  {	int	is_calc_opt_apert;
	char	*basename;

	basename=fits_basename(imgname,&frameno);
	if ( (fr=fopenread(basename))==NULL )
	 {	fprint_error("unable to open input file '%s'.",basename);
		return(1);
	 }
	img=fits_read_frame_to_image(fr,frameno);
	fclose(fr);
	if ( img==NULL )
	 {	fprint_error("unable to interpret input as FITS data.");
		return(1);
	 }
	if ( img->i.dim != 2 )
	 {	fprint_error("image dimension differs from 2.");
		return(1);
	 }

	fits_rescale(img);

	if ( outlistname != NULL )	is_calc_opt_apert=1;
	else				is_calc_opt_apert=0;

	mask=fits_mask_read_from_header(&img->header,img->i.sx,img->i.sy,NULL);
	if ( inmasklist != NULL )
	 {	if ( join_masks_from_files(mask,img->i.sx,img->i.sy,inmasklist) )
		 {	fprint_error("unable to read (one of the) mask file(s)");
			return(1);
		 }		
	 }
	fits_mask_mark_nans(&img->i,mask,MASK_NAN);

	if ( is_auto_gain )
		get_gain_info_from_header(&sg,&img->header);
	do_photometry(&img->i,mask,ps,np,wl,weightusage,is_calc_opt_apert,inaps,numap,&sg,sigma,&xpp);

	fits_free(img);

	if ( outlistname != NULL )
	 {	FILE	*fw;
		fw=fopenwrite(outlistname);
		if ( fw==NULL )		
		 {	fprint_error("unable to create output list file '%s'",outlistname);
			return(1);
		 }
		write_output_star_list(fw,ps,np);
		fclosewrite(fw);
	 }

	if ( outname != NULL )
	 {	FILE	*fw;
		calculate_magnitudes(ps,np,&mf0);
		fw=fopenwrite(outname);
		if ( fw==NULL )	
		 {	fprint_error("unable to create output photometry file '%s'",outname);
			return(1);
		 }
		if ( is_comment )
		 {	fprintf(fw,"# Created by fiphot v%s (fi: %s)\n",FI_PHOT_VERSION,FI_VERSION);
			fprintf(fw,"# Invoked command:");
			for ( i=0 ; i<argc ; i++ )
			 {	if ( is_any_nasty_char(argv[i]) )
					fprintf(fw," \"%s\"",argv[i]);
				else
					fprintf(fw," %s",argv[i]);
			 }
			fprintf(fw,"\n");
		 }
		write_photometry(fw,ps,np,ofxy,ofph,&sg,basistype,nanstring,serialstring,img->i.sx,img->i.sy);
		fclosewrite(fw);
	 }
	if ( outphotname != NULL )
	 {	FILE	*fw;
		fw=fopenwrite(outphotname);
		if ( fw==NULL )	
		 {	fprint_error("unable to create output raw photometry file '%s'",outphotname);
			return(1);
		 }
		write_raw_photometry(fw,ps,np,basistype,nanstring);
		fclosewrite(fw);
	 }

  }
/******************************************************************************
   One reference image and one subtracted image are given: 
   simple aperture photometry on the reference image and after it, add the 
   all aperture photometry results of the subtracted images to the result 
   of the photometry of  the reference image.
******************************************************************************/

 else if ( ( refname != NULL || inphotname != NULL ) && subname != NULL )
  {	kernellist	*kl;
	char		*basename;
	int		i,sx,sy;

	refimg=NULL;
	mask=NULL;
	if ( refname != NULL )
	 {	
		basename=fits_basename(refname,&frameno);
		if ( (fr=fopenread(basename))==NULL )
		 {	fprint_error("unable to open reference image file '%s'.",basename);
			return(1);
		 }
		fclose(fr);
		if ( refimg==NULL )
		 {	fprint_error("unable to interpret reference image input as FITS data.");
			return(1);
		 }
		if ( refimg->i.dim != 2 )
		 {	fprint_error("reference image dimension differs from 2.");
			return(1);
		 }

		mask=fits_mask_read_from_header(&refimg->header,refimg->i.sx,refimg->i.sy,NULL);
		if ( inmasklist != NULL )
		 {	if ( join_masks_from_files(mask,refimg->i.sx,refimg->i.sy,inmasklist) )
			 {	fprint_error("unable to read (one of the) mask file(s)");
				return(1);
			 }
		 }

		logmsg(is_verbose,"Performing reference photometry ... ");
		if ( is_auto_gain )
			get_gain_info_from_header(&sg,&refimg->header);
		do_photometry(&refimg->i,mask,ps,np,wl,weightusage,0,inaps,numap,&sg,sigma,&xpp);
		logmsg(is_verbose,"done.\n");
		for ( i=0 ; i<np ; i++ )
		 {	ps[i].rfflux=ps[i].fluxes;
			ps[i].fluxes=NULL;
		 }
	 }
	else if ( inphotname != NULL )
	 {	fr=fopenread(inphotname);
		if ( fr==NULL )			
		 {	fprint_error("unable to open input raw photomery file '%s'",inphotname);
			return(1);
		 }
		read_raw_photometry(fr,&ps,&np);
		fcloseread(fr);
		numap=0;
		for ( i=0 ; i<np ; i++ )
		 {	if ( ps[i].n>numap )	numap=ps[i].n;		}
	 }

	basename=fits_basename(subname,&frameno);
	if ( (fr=fopenread(basename))==NULL )
	 {	fprint_error("unable to open subtracted image file '%s'.",basename);
		return(1);
	 }
	subimg=fits_read_frame_to_image(fr,frameno);
	fclose(fr);
	if ( subimg==NULL )
	 {	fprint_error("unable to interpret subtracted image input as FITS data.");
		return(1);
	 }
	if ( subimg->i.dim != 2 )
	 {	fprint_error("subtracted image dimension differs from 2.");
		return(1);
	 }

	sx=subimg->i.sx;
	sy=subimg->i.sy;

	if ( refimg != NULL && mask != NULL )
 	 {	submask=fits_mask_duplicate(mask,refimg->i.sx,refimg->i.sy);
		fits_mask_mask_from_header(submask,&subimg->header,refimg->i.sx,refimg->i.sy,NULL);
	 }
	else	submask=fits_mask_read_from_header(&subimg->header,subimg->i.sx,subimg->i.sy,NULL);

	if ( subkernelfile != NULL )
		kl=&klist;
	else
		kl=NULL;

	if ( is_auto_gain )
		get_gain_info_from_header(&sg,&subimg->header);

	logmsg(is_verbose,"Photometry on subtracted image '%s' ... ",subname);
	do_subtracted_photometry(&subimg->i,submask,ps,np,inaps,ninap,&sg,&xpp,kl);
	logmsg(is_verbose,"done.\n");

	fits_mask_free(submask);
	if ( refimg != NULL )	fits_free(refimg);
	if ( mask != NULL )	fits_mask_free(mask);

	if ( mfp.norder>=0 )
		do_magnitude_fit(ps,np,&mfp,sx,sy,&mfs);

	if ( outname != NULL )
	 {	FILE	*fw;
		int	a;
		calculate_magnitudes(ps,np,&mf0);
		fw=fopenwrite(outname);
		if ( fw==NULL )	
		 {	fprint_error("unable to create output photometry file '%s'",outname);
			return(1);
		 }
		if ( is_comment )
		 {	fprintf(fw,"# Created by fiphot v%s (fi: %s)\n",FI_PHOT_VERSION,FI_VERSION);
			fprintf(fw,"# Invoked command:");
			for ( i=0 ; i<argc ; i++ )
			 {	if ( is_any_nasty_char(argv[i]) )
					fprintf(fw," \"%s\"",argv[i]);
				else
					fprintf(fw," %s",argv[i]);
			 }
			fprintf(fw,"\n");
		 }
		write_photometry(fw,ps,np,ofxy,ofph,&sg,basistype,nanstring,serialstring,subimg->i.sx,subimg->i.sy);
		if ( mfp.norder>=0 && is_comment )
		 {	fprintf(fw,"# Magnitude transformation statistics - %d stars:\n",mfs.nstar);	
			for ( a=0 ; a<mfs.naperture ; a++ )
			 {	fprintf(fw,"#\taperture#%d: %5d %5d\n",a+1,mfs.ninit[a],mfs.nrejs[a]);	}
		 }
		fclosewrite(fw);
	 }

  }

 return(0);
}

