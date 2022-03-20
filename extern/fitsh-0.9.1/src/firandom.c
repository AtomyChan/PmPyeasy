/*****************************************************************************/
/* firandom.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line utility to create artifical images.			     */
/*****************************************************************************/
#define	FI_RANDOM_VERSION	"1.3b0"
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <stdarg.h>

#include <fits/fits.h>

#include "longhelp.h"
#include "fi.h"

#include "io/iof.h"
#include "io/scanarg.h"
#include "io/tokenize.h"
#include "math/point.h"

#include "basis.h"
#include "tensor.h"
#include "common.h"
#include "stars.h"
#include "psf.h"
#include "psf-io.h"
#include "weight.h"
#include "magnitude.h"
#include "history.h"
#include "fitsmask.h"
#include "maskdraw.h"

#include "firandom.h"

#ifdef  HAVE_NO_CC_EXTENSION 
#define __extension__ 
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

/* gain: fotons/adu, usually > 1... */
int draw_star_montecarlo_gauss(fitsimage *img,stargenparam *sgp,double x0,double y0,double dflux,double is,double id,double ik)
{
 int	sx,sy,ix,iy,subx,suby,flux;
 double	x,y,fadu;

 if ( img==NULL )	return(-1);
 if ( img->data==NULL )	return(-1);
 sx=img->sx,sy=img->sy;

 if ( sgp->is_intinelect )	flux=(int)dflux;
 else				flux=(int)(dflux*sgp->gain+0.5);

 if ( flux<=0 )		return(1);
 fadu=1.0/sgp->gain;

 if ( sgp->subpixeldata==NULL || sgp->subg==0 )
  {	for ( ; flux>0 ; flux-- )
	 {	get_gaussian_2d(x0,y0,is,id,ik,&x,&y);
		ix=(int)floor(x),
		iy=(int)floor(y);
		if ( ix<0 || iy<0 || ix>=sx || iy>=sy )	continue;
		img->data[iy][ix]+=fadu;
	 }
  }
 else
  {	for ( ; flux>0 ; flux-- )
	 {	get_gaussian_2d(x0,y0,is,id,ik,&x,&y);
		ix=(int)floor(x),
		iy=(int)floor(y);
		if ( ix<0 || iy<0 || ix>=sx || iy>=sy )	continue;
		subx=(int)((x-(double)ix+2)*(double)sgp->subg)%sgp->subg;
		suby=(int)((y-(double)iy+2)*(double)sgp->subg)%sgp->subg;
		img->data[iy][ix]+=sgp->subpixeldata[suby][subx]*fadu;
 	 }
  }
 return(0);
}

int draw_star_from_array(fitsimage *img,int ix0,int iy0,double **iarr,int hsize,int sg,double dflux,double **subpix,double gain)
{
 int	i,j,k,l,fsize,is_subgrid,sx,sy,si,sj,is_photnoise;
 double	sum,w,sw;

 fsize=(2*hsize+1)*sg;

 if ( img==NULL || img->data==NULL )	return(1);
 sx=img->sx,sy=img->sy;
 if ( sx<=0 || sy<=0 )			return(1);

 sum=0.0;
 for ( i=0 ; i<fsize ; i++ )
  {	for ( j=0 ; j<fsize ; j++ )
 	 {	sum+=iarr[i][j];		}
  }
 w=dflux/sum;
 for ( i=0 ; i<fsize ; i++ )
  {	for ( j=0 ; j<fsize ; j++ )
	 {	iarr[i][j]*=w;		}
  }

 if ( sg>1 && subpix!=NULL )	is_subgrid=1;
 else				is_subgrid=0;

 if ( gain > 0.0 )	is_photnoise=1;
 else			is_photnoise=0;

 for ( i=iy0-hsize,si=0 ; i<=iy0+hsize ; i++,si+=sg )
  {	if ( i<0 || i>=sy )	continue;
	for ( j=ix0-hsize,sj=0 ; j<=ix0+hsize ; j++,sj+=sg )
	 {	if ( j<0 || j>=sx )	continue;
		if ( ! is_subgrid )	w=iarr[si][sj];
		else
		 {	w=0.0;
			for ( k=0 ; k<sg ; k++ )
			 {	for ( l=0 ; l<sg ; l++ )
				 {	sw=subpix[k][l];
					w+=sw*iarr[si+k][sj+l];
				 }
			 }
		 }
		if ( is_photnoise && w>0.0 )
			w=get_gaussian(w,sqrt(w/gain));
		if ( w<=0.0 )	continue;

		img->data[i][j]+=w;
	 }
  }

 return(0);
}

int draw_star_integral_gauss(fitsimage *img,stargenparam *sgp,double x0,double y0,double dflux,double is,double id,double ik)
{
 int		hsize,grid,ix0,iy0,fsize;
 double		vgain,dgrid;
 static double	**iarr=NULL;
 static int	afsize=0;

 if ( img==NULL || img->data==NULL )	return(-1);

 if ( sgp->is_intinelect )	dflux/=sgp->gain;
 if ( sgp->is_photnoise )	vgain=sgp->gain;
 else				vgain=0.0;
 if ( dflux*sgp->gain<=0.0 )	return(0);

 hsize=(int)(is*sqrt(2.0*log(dflux*sgp->gain*sgp->nsuppress/(is*is)))+1.0);
 if ( hsize>1023 )	hsize=1023;	/* extraordinary... */

 if ( sgp->subpixeldata==NULL || sgp->subg<=0 )	grid=1;
 else						grid=sgp->subg;

 dgrid=(double)grid;
 fsize=(2*hsize+1)*grid;

 if ( fsize>afsize )
  {	afsize=fsize;
	if ( iarr != NULL )	tensor_free(iarr);
	iarr=(double **)tensor_alloc_2d(double,afsize,afsize);
  }

 ix0=(int)floor(x0),x0-=(double)ix0,
 iy0=(int)floor(y0),y0-=(double)iy0;

 star_draw_gauss(iarr,fsize,fsize,(x0+hsize)*dgrid,(y0+hsize)*dgrid,
		 is*dgrid,id*dgrid,ik*dgrid);

 draw_star_from_array(img,ix0,iy0,iarr,hsize,grid,dflux,sgp->subpixeldata,vgain);

 return(0);
}

int draw_star_integral_deviated(fitsimage *img,stargenparam *sgp,double x0,double y0,double dflux,double gs,int order,double *mom)
{
 int		hsize,grid,ix0,iy0,fsize;
 double		is,vgain,dgrid;
 static double	**iarr=NULL;
 static int	afsize=0;
 double		cgs,gmom[MAX_DEVIATION_COEFF],*cmom;

 if ( img==NULL || img->data==NULL )	return(-1);

 if ( sgp->is_intinelect )	dflux/=sgp->gain;
 if ( sgp->is_photnoise )	vgain=sgp->gain;
 else				vgain=0.0;
 if ( dflux*sgp->gain<=0.0 )	return(0);
 if ( gs<=0.0 )			return(0);

 ix0=(int)floor(x0),x0-=(double)ix0,
 iy0=(int)floor(y0),y0-=(double)iy0;

 is=1/sqrt(gs);
 hsize=(int)(sqrt(2.0*log(gs*dflux*sgp->gain*sgp->nsuppress)/gs)+1.0);

 if ( sgp->subpixeldata==NULL || sgp->subg<=0 )	grid=1;
 else						grid=sgp->subg;

 fsize=(2*hsize+1)*grid;
 dgrid=(double)grid;

 if ( fsize>afsize )
  {	afsize=fsize;
	if ( iarr  != NULL )	tensor_free(iarr);
	iarr =(double **)tensor_alloc_2d(double,afsize,afsize);
  }

 if ( grid>1 )
  {	double	igrid,w;
	int	o,l,j;

	igrid=1.0/dgrid;
	w=igrid*igrid;
	cgs=w*gs;
	for ( o=2,l=0 ; o<=order ; o++ )
	 {	for ( j=0 ; j<=o ; j++,l++ )
		 {	gmom[l]=mom[l]*w;
			w=w*igrid;
		 }
	 }
	cmom=gmom;
  }
 else
  {	cgs =gs;
	cmom=mom;
  }

 star_draw_deviated(iarr,fsize,fsize,(x0+hsize)*dgrid,(y0+hsize)*dgrid,cgs,order,cmom);

 draw_star_from_array(img,ix0,iy0,iarr,hsize,grid,dflux,sgp->subpixeldata,vgain);
 
 return(0);
}

int draw_star_integral_psf(fitsimage *img,stargenparam *sgp,double x0,double y0,
	double px,double py,
	double is,double id,double ik,double il,double dflux)
{
 static double	**iarr=NULL;
 static int	afsize=0;
 int		hsize,fsize,nvar,bx,by;
 int		ix0,iy0,grid;
 double		vgain,dgrid;
 double		det,nrm;
 psf		*p;

 if ( img==NULL || img->data==NULL )	return(-1);

 if ( sgp->is_intinelect )	dflux/=sgp->gain;
 if ( sgp->is_photnoise )	vgain=sgp->gain;
 else				vgain=0.0;
 if ( dflux*sgp->gain<=0.0 )	return(0);

 det=is*is+il*il-id*id-ik*ik;
 if ( is<=0.0 && det<=0.0 )	nrm=1.0;
 else				nrm=sqrt(is*is+il*il+id*id+ik*ik);

 ix0=(int)floor(x0),x0-=(double)ix0,
 iy0=(int)floor(y0),y0-=(double)iy0;

 if ( sgp->subpixeldata==NULL || sgp->subg<=0 )	grid=1;
 else						grid=sgp->subg;

 p=sgp->tpd;
 if ( p==NULL )	return(1);

 nvar=(p->order+1)*(p->order+2)/2;
 bx=p->grid*(2*p->hsize+1);
 by=p->grid*(2*p->hsize+1);

 hsize=(int)(nrm*(double)(p->hsize+2));
 fsize=grid*(2*hsize+1);
 dgrid=(double)grid;

 if ( fsize>afsize )
  {	if ( iarr  != NULL )	tensor_free(iarr);
	iarr =(double **)tensor_alloc_2d(double,fsize,fsize);
	afsize=fsize;
  }

 star_draw_psf(iarr,fsize,fsize,(x0+hsize)*dgrid,(y0+hsize)*dgrid,p,px,py,
	is*dgrid,id*dgrid,ik*dgrid,il*dgrid);

 draw_star_from_array(img,ix0,iy0,iarr,hsize,grid,dflux,sgp->subpixeldata,vgain);
 
 return(0); 
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int draw_starlist(fitsimage *img,stargenparam *sgp,star *stars,int nstar,int zoom)
{
 int	o,j,l,n,order;
 double	s,d,k,dflux,zf,iz,gs,mom[MAX_DEVIATION_COEFF],*wmom,wz;
 star	*ws;

 if ( zoom<1 )	zf=1.0,zoom=1;
 else		zf=(double)zoom;
 iz=1.0/(double)zoom;

 for ( n=0 ; n<nstar ; n++ )
  {	ws=&stars[n];

	dflux=ws->flux;
	if ( dflux<=0.0 )	continue;

	if ( ws->shape.model==SHAPE_GAUSS || ws->shape.model==SHAPE_ELLIPTIC )
	 {	s=ws->gsig;
		if ( ws->shape.model==SHAPE_ELLIPTIC )
			d=ws->gdel,k=ws->gkap;
		else	
			d=k=0;

		if ( ! sgp->method )
		 {	if ( sgp->is_photnoise )
				dflux=get_gaussian(dflux,sqrt(dflux/sgp->gain));
			draw_star_montecarlo_gauss(img,sgp,ws->location.gcx*zf,ws->location.gcy*zf,dflux,s*zf,d*zf,k*zf);
		 }
		else
			draw_star_integral_gauss(img,sgp,ws->location.gcx*zf,ws->location.gcy*zf,dflux,s*zf,d*zf,k*zf);
	 }
	else if ( ws->shape.model==SHAPE_DEVIATED )
	 {	if ( zoom>1 )
	 	 {	wz=iz*iz;
			gs=ws->shape.gs*wz;
			order=ws->shape.order;
			for ( o=2,l=0 ; o<=order ; o++ )
			 {	for ( j=0 ; j<=o ; j++,l++ )
				 {	mom[l]=ws->shape.mom[l]*wz;	}
				wz=wz*iz;
			 }
			wmom=mom;
		 }
		else
		 {	gs=ws->shape.gs;
			order=ws->shape.order;
			wmom=ws->shape.mom;
		 }
		draw_star_integral_deviated(img,sgp,ws->location.gcx*zf,ws->location.gcy*zf,dflux,gs,order,wmom);
	 }
	else if ( ws->shape.model==SHAPE_PSF )
	 {	draw_star_integral_psf(img,sgp,ws->location.gcx*zf,ws->location.gcy*zf,
		ws->location.gcx,ws->location.gcy,
		ws->shape.gs*zf,ws->shape.gd*zf,ws->shape.gk*zf,ws->shape.gl*zf,dflux);
	 }
  }
 return(0);		
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int quantize_image(fitsimage *img)
{
 int	i,j,sx,sy;

 if ( img==NULL || img->data==NULL )	return(1);
 sx=img->sx,sy=img->sy;
 if ( sx<=0 || sy<=0 )			return(1);
 
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	img->data[i][j]=floor(img->data[i][j]);		}
  }

 return(0);
}

int divide_image(fitsimage *img,double d)
{
 int	i,j,sx,sy;

 if ( img==NULL || img->data==NULL )	return(1);
 sx=img->sx,sy=img->sy;
 if ( sx<=0 || sy<=0 )			return(1);
 if ( d<=0.0 )				return(1);
 d=1.0/d;
 
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	img->data[i][j]=d*img->data[i][j];		}
  }

 return(0);
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

#define		MAX_COL	32

int read_input_list(FILE *fr,star **rstars,int *rnstar,inlistparam *ilp)
{
 star	*stars,*ws;
 int	nstar,i,n,maxcol,model,order;
 char	*rbuff,*cmd_static[MAX_COL],**cmd,**cmd_dynamic;
 double	cx,cy,cf,m1,m2,m3,s,d,k,gs,gd,gk,gl,mom[MAX_DEVIATION_COEFF],bshx,bshy;

 maxcol=0;
 if ( ilp->colx>maxcol )	maxcol=ilp->colx;
 if ( ilp->coly>maxcol )	maxcol=ilp->coly;
 if ( ilp->colmag>maxcol )	maxcol=ilp->colmag;
 if ( ilp->colflux>maxcol )	maxcol=ilp->colflux;
 switch ( ilp->listtype )
  {  case LISTTYPE_FEP: case LISTTYPE_SIG: case LISTTYPE_SDK:	/* 3 sh.cols */
	if ( ilp->colsh1>maxcol )	maxcol=ilp->colsh1;
	if ( ilp->colsh2>maxcol )	maxcol=ilp->colsh2;
	if ( ilp->colsh3>maxcol )	maxcol=ilp->colsh3;
	break;
     case LISTTYPE_SMM:						/* 2 sh.cols */
	if ( ilp->colsh1>maxcol )	maxcol=ilp->colsh1;
	if ( ilp->colsh2>maxcol )	maxcol=ilp->colsh2;
	break;
     case LISTTYPE_PDS:						/* 4 sh.cols */
	if ( ilp->colsh1>maxcol )	maxcol=ilp->colsh1;
	if ( ilp->colsh2>maxcol )	maxcol=ilp->colsh2;
	if ( ilp->colsh3>maxcol )	maxcol=ilp->colsh3;
	if ( ilp->colsh4>maxcol )	maxcol=ilp->colsh4;
	break;
     case LISTTYPE_PSF:						/* no sh.cols*/
	break;
     default:
	break;
  }
 maxcol++;

 if ( maxcol>MAX_COL )	cmd=cmd_dynamic=(char **)malloc(sizeof(char *)*(maxcol+1));
 else			cmd=cmd_static,cmd_dynamic=NULL;

 bshx=bshy=0.0;
 stars=NULL;nstar=0;

 while ( ! feof(fr) )
  {	rbuff=freadline(fr);

	if ( rbuff==NULL )
		break;
	else if ( basis_read_comment(rbuff,&bshx,&bshy) )
		continue;

	remove_newlines_and_comments(rbuff);
	n=tokenize_spaces(rbuff,cmd,maxcol);
	if ( n<maxcol )	continue;

	cx=cy=cf=0.0;
	i =sscanf(cmd[ilp->colx  ],"%lg",&cx);
	i+=sscanf(cmd[ilp->coly  ],"%lg",&cy);
	if ( ilp->colflux>=0 )
	 {	i+=sscanf(cmd[ilp->colflux],"%lg",&cf);		}
	else if ( ilp->colmag>=0 )
	 {	i+=sscanf(cmd[ilp->colmag],"%lg",&cf);
		cf=mag_to_flux(cf,&ilp->mf0);
	 }
	if ( i<3 || ! ( isfinite(cx) && isfinite(cy) && isfinite(cf) )	) continue;

	switch ( ilp->listtype )	
	 {   case LISTTYPE_FEP: case LISTTYPE_SIG: case LISTTYPE_SDK:
		gs=gd=gk=gl=0.0;order=0;	/*** fep, sdk or SDK ****/
		i =sscanf(cmd[ilp->colsh1],"%lg",&m1);
		i+=sscanf(cmd[ilp->colsh2],"%lg",&m2);
		i+=sscanf(cmd[ilp->colsh3],"%lg",&m3);
		if ( i<3 || ! ( isfinite(m1)  && isfinite(m2)  && isfinite(m3) ) ) continue;
		switch ( ilp->listtype )
		 {	case LISTTYPE_FEP:  fep_to_sdk(m1,m2,m3,&s,&d,&k);break;
			case LISTTYPE_SIG:  s=m1;d=m2;k=m3;break;
			case LISTTYPE_SDK: isdk_to_sdk(m1,m2,m3,&s,&d,&k);break;
			default: s=m1;d=m2;k=m3;break;
		 }
		model=SHAPE_ELLIPTIC;
		break;
	     case LISTTYPE_SMM:
		s=d=k=0.0;			/*** Smom ***************/
		i=sscanf(cmd[ilp->colsh1],"%lg",&gs);
		if ( i<1 || ! ( isfinite(gs) || gs<=0.0 ) ) continue;
		i=sscanf(cmd[ilp->colsh2],"%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg",
			&mom[0],&mom[1],&mom[2],&mom[3],&mom[4],&mom[5],&mom[6],
			&mom[7],&mom[8],&mom[9],&mom[10],&mom[11]);
		if ( i==12 )		order=4;
		else if ( i>=7 )	order=3;
		else if ( i>=3 )	order=2;
		else			continue;
		model=SHAPE_DEVIATED;
		gd=gk=gl=0.0;
		break;
	    case LISTTYPE_PSF:
		s=d=k=0.0; 
		gs=1.0;gd=gk=gl=0.0;
		model=SHAPE_PSF;
		order=0;
		break;
	    case LISTTYPE_PDS:
		s=d=k=0.0;
		i =sscanf(cmd[ilp->colsh1],"%lg",&gs);
		i+=sscanf(cmd[ilp->colsh2],"%lg",&gd);
		i+=sscanf(cmd[ilp->colsh3],"%lg",&gk);
		i+=sscanf(cmd[ilp->colsh4],"%lg",&gl);
		if ( i<4 || ! ( isfinite(gs)  && isfinite(gd)  && isfinite(gk) && isfinite(gl) ) ) continue;
		model=SHAPE_PSF;
		order=0;
		break;
	    default:
		model=0;
		order=0;
		break;
	 }
	if ( ! model )	continue;

	stars=(star *)realloc(stars,sizeof(star)*(nstar+1));

	ws=&stars[nstar];
	ws->location.gcx=cx,
	ws->location.gcy=cy,
	ws->flux=cf;
	ws->gsig=s;
	ws->gdel=d;
	ws->gkap=k;
	ws->shape.model=model;
	ws->shape.order=order;
	ws->shape.gs=gs;
	ws->shape.gd=gd;
	ws->shape.gk=gk;
	ws->shape.gl=gl;
	for ( i=0 ; i<(order+1)*(order+2)/2-3 ; i++ )
	 {	ws->shape.mom[i]=mom[i];		}
	for ( ; i<MAX_DEVIATION_COEFF ; i++ )
	 {	ws->shape.mom[i]=0.0;			}

	nstar++;
  };

 *rstars=stars;
 *rnstar=nstar;

 if ( cmd_dynamic != NULL )	free(cmd_dynamic);

 return(0);
}

int write_output_list(FILE *fw,star *stars,int nstar,int listtype,magflux *mf,int basistype)
{
 int	i;
 double	s,d,k,m1,m2,m3,bshx,bshy;

 if ( is_comment )
  {	fprintf(fw,"# Star list, dumped by firandom.\n");
	if ( mf==NULL )
		fprintf(fw,"#     X           Y         FLUX     ");
	else
		fprintf(fw,"#     X           Y         Magnitude");
	if ( listtype==LISTTYPE_FEP ) 
		fprintf(fw,"   FWHM       Ellipt.     PosAng. (--fep)");
	else if ( listtype==LISTTYPE_SIG )
		fprintf(fw,"   sigma      delta       kappa   (--sdk)");
	else if ( listtype==LISTTYPE_SDK )
		fprintf(fw,"   S          D           K       (--SDK)");
	fprintf(fw," IdNum\n");
  }
 basis_set_shift_coords(basistype,&bshx,&bshy);
 basis_write_comment(fw,bshx,bshy);

 for ( i=0 ; i<nstar ; i++ )
  {	fprintf(fw,"%11.5f %11.5f ",stars[i].location.gcx+bshx,stars[i].location.gcy+bshy);
	if ( mf==NULL )	fprintf(fw,"%11.3f ",stars[i].flux);
	else		fprintf(fw,"%11.4f ",flux_to_mag(stars[i].flux,mf));
	s=stars[i].gsig,
	d=stars[i].gdel,
	k=stars[i].gkap;
	switch ( listtype )
	 {	case 0 : sdk_to_fep(s,d,k,&m1,&m2,&m3);break;
		case 2 : sdk_to_isdk(s,d,k,&m1,&m2,&m3);break;
		default: m1=s,m2=d,m3=k;break;
	 }
	fprintf(fw,"%11g %11g %11g %11d\n",m1,m2,m3,i+1);
  }
 return(0);
}

/*****************************************************************************/

int fprint_firandom_usage(FILE *fw)
{
 fprintf(fw,
"Usage:\tfirandom [-o <out.fits>] [-h|--help] [-C|--comment] [-V|--verbose]\n"
"Common modelling parameters:\n"
"\t[-g|--gain <gain>]  [-m \"<sky-as-funct-of-x,y>\"] [-d <extra-noise>]\n"
"\t[--photon-noise|--no-photon-noise] [-U <subpixelfile> [--normalize]]\n"
"Output image size and data format parameters:\n"
"\t[-s|--size <sx>,<sy>] [-b|--bitpix <bitpix>]\n"
"\t[-D|--data bitpix=<bitpix>,bscale=<scale>,bzero=<zero>|<C-type>]\n");
 fprintf(fw,
"Input specification:\n"
"\t[-L <in-starlist>|-l \"<star-parameters>\"] [--output-list <out-starlist>]\n"
"\t[--{fep|sdk|SDK|Smom|psf|psfdst}{|-input|-output}] [-ip <input-psf>]\n"
"\t[--col-xy <>,<>] [--col-{flux|mag} <>] [--col-shape <>[,<>,<>[,<>]]]\n");
 fprintf(fw,
"Fine-tuning of modelling:\n"
"\t[--seed {<seed>|auto}] [--seed-{noise|spatial|photon} {<seed>|auto}]\n"
"\t[--mag-flux <mag>,<flux>]\n"
"\t[--[no-]monte-carlo|--[no-]integral [--noise-suppression <level>]]\n"
"\t[--[no-]quantize] [--[no-]adus|--[no-]electrons]\n");
 fprintf(fw,
"Handling masks:\n"
"\t[-q|--mask-block <maskflags>:<x0>,<y0>[:<x1>,<y1>]] [--mask-block ...]\n"
"\t[-om|--output-mask <outputmaskfile>\n");
 fprintf(fw,
"Weight output generation:\n"
"\t[--output-weight <weight.fits> [-w|--weight hsize=<hsize>,zoom=<z>]]\n");
 fprintf(fw,
"Default parameters:\n"
"\t-s 256,256 -m 500 -d 10 -b 16 --mag-flux 10,10000 --adus --quantize\n"
"\t--gain 1.0 --no-photon-noise --monte-carlo --noise-supression 1.0\n"
"\t--seed 0 --col-xy 1,2 --col-flux 3 --col-shape 4,5,6\n");
 fprintf(fw,
"Format of <star-parameters>: \n"
"\t[f=<FWHM>],[e=<ellip>],[p=<position-angle (degrees)>],\n"
"\t[s=<sigma>],[d=<delta>],[k=<kappa>],[S=<S>],[D=<D>],[K=<K>],\n"
"\t<N>*[{x|X}=<x>,{y|Y}=<y>,{m=<mag>|i=<flux>},[{f,e,p|s,d,k|S,D,K}=...]],\n"
"\t[f,e,p=...],[s,d,k=...],[S,D,K=...],<Num>*[...],...\n");

 return(0);
}

longhelp_entry firandom_long_help[]=
{
 LONGHELP_OPTIONS,

 { "General options:", NULL },
 { "-h, --help",
	"Gives general summary about the command line options." },
 { "--long-help",
	"Gives a detailed list of command line options." },
 { "--version",
	"Gives some version information about the program." },

 { "Creating artifical object lists:", NULL },
 { "-l, --list <list>",
	"Specifications for object list. The ``list'' parameter should be a "
	"set of comma separated tags, which can either be a value declaration "
	"or a repeat count followed by an expression between square brackets "
	"giving specifications for individual objects to be added to the list, "
	"also in the form of value declaration. The value declaration "
	"has the sintax of <variable>=<value>, where the variables can be the "
	"following: " },
 { "f",	"full width at half magnitude (FWHM) of the stellar "
	"profiles to be created" },
 { "e", "ellipticity of the stellar profiles" },
 { "p", "position angle of the stellar profiles" },
 { "s", "sigma parameter for the stellar profile (FWHM is roughly 2.35 * sigma)" },
 { "d", "delta (plus-shaped deviance) parameter for the stellar profile" },
 { "k", "kappa (cross-shaped deviance) parameter for the stellar profile" },
 { "S", "Gaussian momenum (a.k.a. profile sharpness parameter) for the stellar profile (S=1/sigma^2)" },
 { "D", "plus-shaped momentum for the stellar profile" },
 { "K", "cross-shaped momentum for the stellar profile" },
 { "x", "normalized X coordinate of the profile centroid (using the standard normalization)" },
 { "X", "absolute X coordinate of the profile centroid" },
 { "y", "normalized y coordinate of the profile centroid (using the standard normalization)" },
 { "Y", "absolute Y coordinate of the profile centroid" },
 { "m", "magnitude of the stellar object" },
 { "i", "flux of the stellar object" },

 { __extension__ 
   "One can use only one of the three equivalent set of profile shape "
   "parameters (i.e. f/e/p, s/d/k or S/D/K). See some more detailed documentation "
   "about these parameters. In the expressions which are in the square brackets, "
   "one can use arbitrary arithmetic expressions, using the standard basic "
   "arithmetical operators, elementary functions and the functions r(lo,hi) "
   "and g(mean,sigma) which results an uniformly distrbuted random number "
   "between ``lo'' and ``hi'' and a Gaussian random number with the specified "
   "``mean'' and ``sigma'' (standard deviation), respectively. "
   "In the expressions for magnitude or intensity, one can use the previously "
   "defined values for the centroid coordinates too. The variable ``n'' is "
   "increased between 0 and the repeat count during the evaluation of the "
   "square bracket expressions. ", NULL }, { "", NULL },

 { "--output-list <list file>",
	"Name of the list file where the object list created by the subsequent "
	"--list options are saved. " },
 { "--fep, --fep-output",
	"Save the shape parameters as FWHM, ellipticity and position angle to "
	"the output list file." },
 { "--sdk, --sdk-output",
	"Save the shape parameters as sigma, delta and kappa to "
	"the output list file." },
 { "--SDK, --SDK-output",
	"Save the shape parameters as Gaussian momenta to "
	"the output list file." },

 { "Creating artifical images:", NULL },
 { "-s, --size <sx>,<sy>",
	"Size of the image to be created." },
 { "-b, --bitpix <bitpix>",
        "Standard FITS output bitpix value." },
 { "-D, --data <spec>",
        "Output pixel data format specification." },
 { "-m, --sky <sky>",
	"Sky (background level) for the image. This can be either a constant "
	"or an arbitrary function of the x, y, X and Y coordinates (see above) "
	"for a backgroud with shows systematic variations. One can use the "
	"previously discussed r(lo,hi) and g(mean,sigma) functions, in order "
	"to add some sort of noise to the background. " },
 { "-d, --sky-noise <noise>",
	"Additional Gaussian noise, equivalent to the term ``+g(<noise>,0)'' "
	"added after the background level expression. " },
 { "--photon-noise, --no-photon-noise",
	"Emulate or disable the effect of photon noise on the individual "
	"stellar objects." },
 { "-l, --list <list>",
	"Specifications for object list (see above)." },

 { "-L, --input-list <list file>",
	"Name of the input list file from which the coordinates, shape "
	"parameters and intensities are read for the individual objects. " },
 { "--col-xy <colx>,<coly>",
	"Column indices for X and Y (absolute) centroid coordinates. "},
 { "--col-flux <flux column>",
	"Column index for flux (intensity)." },
 { "--col-shape <profile width>,[<profile shape 1>,<profile shape 2>]",
	"Column indices for stellar profile parameters. Either 1 or 3 columns "
	"should be specified following this command line switch. One shape parameter "
	"is interpreted as a profile size parameter where the 2 additional "
	"(optional) shape parameters describe the deviation from the symmetric "
	"profile. See also options --fep, --sdk or --SDK for more details. " },
 { "--fep, --fep-input",
	"Interpret the shape parameters read from the input list file as "
	"FWHM, ellipticity and position angle." },
 { "--sdk, --sdk-input",
	"Interpret the shape parameters read from the input list file as "
	"the sigma, delta and kappa parameters." },
 { "--SDK, --SDK-input",
	"Interpret the shape parameters read from the input list file as "
	"the Gaussian momenta parameters." },

 { "-S, --input-sky <sky list file>",
	"Name of the input file containing the sky level. This file "
	"should contain at least three columns: the two pixel coordinates "
	"and the sky vaule. See also --col-pixel and --col-value. " },
 { "--col-pixel <colx>,<coly>",
	"Column indices for X and Y (absolute) pixel coordinates. "},
 { "--col-value<sky value column>",
	"Column index for sky value (intensity)." },

 { "--mag-flux <mag>,<flux>",
	"Magnitude - flux conversion level. The specified magnitude will "
	"be equivalent to the specified flux level." },
 { "--integral, --no-monte-carlo",
	"Draw the stellar profiles to the image using exact integration. " },
 { "--monte-carlo, --no-integral",
	"Draw the stellar profiles to the image using a Monte-Carlo way. "
	"Note that using this Monte-Carlo method without additional photon "
	"noise emulation would result assymetric stellar profiles even when "
	"the profile would be symmetric. Use this option only when "
	"the --photon-noise option is also used, therefore the profiles "
	"are strained with photon noise either. " },
 { "--noise-suppression <level>",
	__extension__ 
	"If the profiles are drawn using exact integration, the profiles would be "
	"infintely large since an analytical Gaussian profile is positive "
	"on the whole image domain. In order to limit the integration boundaries, "
	"this level limits the size of the integration domain, by the following way. "
	"The expected level of the objects's own photon noise at the edges of the "
	"integration domain is smaller by this factor at least than the flux "
	"level. Higher suppression level results larger integration domain. "
	"In the case of additional photon noise, the default value of 1.0 "
	"is satisfactory. For images with no photon noise, this level should "
	"be increased appropriately. " },
 { "--quantize, --no-quantize",
	"Quantize the output images to integers or not. Note that altering "
	"this option yields somehow the same as when the bitpix value "
	"is altered. " },
 { "--adus, --no-electrons",
	"Use the input fluxes as ADUs instead of electrons (default). " },
 { "--electrons, --no-adus",
	"Use the input fluxes as electrons insead of ADUs. " },
 { "-g, --gain <gain>",
	"Electron/ADU ratio (gain)." }, 

 { "Random seeds:", NULL}, 

 { "--seed <seed>|auto",
	"Generic random seed for `firandom`. A literal ``auto'' argument yields "
	"a random seed derived from random sources available on the "
	"architecture (/dev/urandom, current time)." },
 { "--seed-noise <seed>|auto",
	"Specific random seed for creating background noise." },
 { "--seed-spatial <seed>|auto",
	"Specific random seed for creating random spatial coordinates, i.e. "
	"the random seed for functions in the --list arguments. " },
 { "--seed-photon <seed>|auto",
	"Specific random seed for photon noise." },
 
 { "Command line argument combinations:", NULL },

 { "--list <list> --output-list <list file>",
	"This combination creates only a list file based on the --list arguments." },
 { "--input-list <list file> --output-list <list file>",
	"This combination just filters and copies the relevant contents "
	"from the input list to the output list. The shape parameters might "
	"be converted, for example --SDK-input --fep-output would convert "
	"Gaussian momenta to FWHM, ellipticity and position angle. " },
 { "--list <list> --output <output image> [--output-list <list file>]",
	"This combination creates an artifical list of sources and then creates "
	"an artifical image with this newly created set of objects. By default, "
	"the list itself (incl. the centroid coordinates, shape parameters and "
	"intensities) is not saved unless an output list file is given." },
 { "--input-list <list file> --list <list>",
	"This combination is invalid, the centroid list must either be read "
	"from a file or created by the program invocation but lists cannot "
	"be merged this way. In such case, save the object list to a separete "
	"file and merge the files using standard tools. " },
	
 { NULL, NULL }

};


int fprint_firandom_long_help(FILE *fw)
{
 fprintf(fw,
"Usage:\tfirandom [options] [-o|--output <output>]\n"
"The main purpose of this program is to generate artifical object lists and/or \n"
"artifical (astronomical) images.\n");
 fprintf(fw,"\n");

 longhelp_fprint(fw,firandom_long_help,0,-1);

 fprintf(fw,"\n");
 fprintf(fw,"Report bugs to <%s>\n",FI_MAINT_EMAIL);

 return(0);
}

/*****************************************************************************/

int get_seed_arg(char *arg)
{
 int	r;
 if ( strcmp(arg,"auto")==0 )	return(-1);
 if ( sscanf(arg,"%d",&r)<1 )	return(0);
 if ( ! isfinite(r) )		return(0);
 if ( r>0 )	return(r);
 else		return(-r);
}

char	*rndfiles[]= { "/dev/random", "/dev/urandom", NULL };

int get_random_seed(void)
{
 FILE	*fr;
 int	seed,i;

 for ( i=0,fr=NULL ; rndfiles[i] != NULL && fr==NULL ; i++ )
  {	fr=fopen(rndfiles[i],"rb");
	if ( fr != NULL )	break;
  }

 if ( fr != NULL )
  {	(void)fread(&seed,1,sizeof(int),fr);
	fclose(fr);
  }
 else
  {	time_t	ctm;
	time(&ctm);
	seed=ctm;
  }

 return(seed);
}
int set_seed(int *rseed,char *strseed)
{
 if ( strseed==NULL )
  {	*rseed=0;
	return(0);
  }
 else if ( strcmp(strseed,"auto")==0 )
  {	*rseed=get_random_seed();
	if ( *rseed<0 )	*rseed=-(*rseed);
	return(0);
  }
 else
  {	if ( sscanf(strseed,"%d",rseed)==1 )	return(0);
	else					return(1);
  }
}


int main(int argc,char *argv[])
{
 int		i,sx,sy,aseed,nseed,pseed,sseed,is_help;
 double		stddev;
 char		*subpixelfile,*psffile,*outfile,*ilist,*olist,*oweight,*ifsky,
		*stararg,*bgarg,*weightarg,*outmaskname,**maskblocklist;
 char		*saseed,*snseed,*spseed,*ssseed,*fdpstring,*shapepar;
 int		subpixelfile_normalize,is_exact,zoom,wgthsize,wgthzoom;

 int		ilisttype, /* listtype=0: --fep : x,y,flux,FWHM,ELL,PA 	     */
		olisttype; /* listtype=1: --sdk : x,y,flux,sigma,delta,kappa */
			   /* listtype=2: --SDK : x,y,flux,S/D/K momenta     */
			   /* listtype=3: --Smom: x,y,flux,S, M20/M11/M02/...*/
			   /* listtype=4: --psf	: x,y,flux + external PSF!   */
			   /* listtype=5: --psfdst: x,y,flux + PSF + distort */

 int		basistype,is_out_mag;
 star		*stars;
 int		nstar;
 double		ox,oy,scale;
 fitsdataparam	fdp;
 stargenparam	sgp;
 starlistparam	slp;
 inlistparam	ilp;
 psf		tpd;
 char		**mask;
 int		col_pix_x,col_pix_y,col_value;

 progbasename=strrchr(argv[0],'/');
 if ( progbasename != NULL )	progbasename++;
 else				progbasename=argv[0];

 is_verbose=is_comment=is_help=0;
 sx=sy=256;
 fdp.bitpix=16;fdp.is_scale=0;fdp.bscale=1.0;fdp.bzero=0.0;fdp.nquantizebit=0;
 fdpstring=NULL;

 stddev=0.0;sgp.gain=1.0;sgp.nsuppress=100.0;
 sgp.subg=0;sgp.subpixeldata=NULL;
 sgp.is_photnoise=0; sgp.method=0;
 sgp.is_intinelect=sgp.dontquantize=0;
 is_exact=0;
 sgp.tpd=NULL;

 psffile=NULL;subpixelfile=NULL;subpixelfile_normalize=0;
 outfile=outmaskname=NULL;
 ifsky=ilist=olist=oweight=stararg=NULL;weightarg=NULL;
 basistype=ilisttype=olisttype=0; is_out_mag=0;

 maskblocklist=NULL;
 
 slp.mf0.intensity=10000.0,
 slp.mf0.magnitude=   10.0;
 aseed=sseed=nseed=pseed=0;bgarg=NULL;
 saseed=ssseed=snseed=spseed=NULL;

 ilp.listtype=ilisttype,
 ilp.colx=1,ilp.coly=2,
 ilp.colflux=3;
 ilp.colmag=0;
 col_pix_x=1;
 col_pix_y=2;
 col_value=3;
 shapepar=NULL;

 zoom=1;

 i=scanarg(argc,argv,SCANARG_ALLOW_FLAGS,
	"--version:%NS-1f%q",&is_help,
        "--version-short|--short-version:%NS-2f%q",&is_help,
	"-h|--help:%f%q",&is_help,
	"--long-help|--help-long:%SN2f%q",&is_help,
	"-o|--output:%s",&outfile,
	"--output-mask:%s",&outmaskname,
	"-L|--input-list:%s",&ilist,
	"-S|--input-sky:%s",&ifsky,
	"-P|--input-psf:%s",&psffile,
	"-l|--list:%Cr",&stararg,
	"--output-list:%s",&olist,
	"--output-weight:%s",&oweight,
	"-w|--weight:%s",&weightarg,

	"--fep:"    SNf(LISTTYPE_FEP)SNf(LISTTYPE_FEP),&ilisttype,&olisttype,
	"--fep-input:"	  SNf(LISTTYPE_FEP),&ilisttype,
	"--fep-output:"	  SNf(LISTTYPE_FEP),&olisttype,
	"--sdk:"    SNf(LISTTYPE_SIG)SNf(LISTTYPE_SIG),&ilisttype,&olisttype,
	"--sdk-input:"	  SNf(LISTTYPE_SIG),&ilisttype,
	"--sdk-output:"	  SNf(LISTTYPE_SIG),&olisttype,
	"--SDK:"    SNf(LISTTYPE_SDK)SNf(LISTTYPE_SDK),&ilisttype,&olisttype,
	"--SDK-input:"	  SNf(LISTTYPE_SDK),&ilisttype,
	"--SDK-output:"	  SNf(LISTTYPE_SDK),&olisttype,
	"--Smom:"   SNf(LISTTYPE_SMM)SNf(LISTTYPE_SMM),&ilisttype,&olisttype,
	"--Smom-input:"	  SNf(LISTTYPE_SMM),&ilisttype,
	"--Smom-output:"  SNf(LISTTYPE_SMM),&olisttype,
	"--psf:"    SNf(LISTTYPE_PSF)SNf(LISTTYPE_PSF),&ilisttype,&olisttype,
	"--psf-input:"	  SNf(LISTTYPE_PSF),&ilisttype,
	"--psf-output:"	  SNf(LISTTYPE_PSF),&olisttype,
	"--psfdst:" SNf(LISTTYPE_PDS)SNf(LISTTYPE_PDS),&ilisttype,&olisttype,
	"--psfdst-input:" SNf(LISTTYPE_PDS),&ilisttype,
	"--psfdst-output:"SNf(LISTTYPE_PDS),&olisttype,

	"--coord=native:%SN0f",&basistype,
	"--coord=iraf:%SN1f",&basistype,
	"--coord=zero:%SN-1f",&basistype,

	"--col-xy:%Pd,%Pd",&ilp.colx,&ilp.coly,
	"--col-flux|--col-int:%Pd%SN0f",&ilp.colflux,&ilp.colmag,
	"--col-mag:%Pd%SN0f",&ilp.colmag,&ilp.colflux,
	"--col-shape:%s",&shapepar,
	"--col-pixel:%Pd,%Pd",&col_pix_x,&col_pix_y,
	"--col-value:%Pd",&col_value,
	"--write-flux:%SN0f",&is_out_mag,
	"--write-mag|--write-magnitude:%SN1f",&is_out_mag,

	"-q|--mask-block:%Dt",&maskblocklist,

	"-s|--size:%d,%d",&sx,&sy,
	"-b|--bitpix:%d",&fdp.bitpix,
	"-D|--data:%s",&fdpstring,
	"-d|--sky-noise:%g",&stddev,
	"-m|--sky|--background:%s",&bgarg,
	"-g|--gain:%g",&sgp.gain,
	"-z|--zoom:%d",&zoom,
	"--no-photon-noise:%SN0f",&sgp.is_photnoise,
	"--photon-noise:%SN1f",&sgp.is_photnoise,
	"--monte-carlo|--no-integral:%SN0f",&sgp.method,
	"--integral|--no-monte-carlo:%SN1f",&sgp.method,
	"--exact:%f",&is_exact,
	"--noise-suppression:%g",&sgp.nsuppress,
	"--adus|--no-electrons:%SN0f",&sgp.is_intinelect,
	"--electrons|--no-adus:%SN1f",&sgp.is_intinelect,
	"--no-quantize:%SN1f",&sgp.dontquantize,
	"--quantize:%SN0f",&sgp.dontquantize,
	"--mag-flux:%g,%g",&slp.mf0.magnitude,&slp.mf0.intensity,
	"-U|--input-subpixel-structure:%s",&subpixelfile,
	"--normalize:%f",&subpixelfile_normalize,
	"--seed:%s",&saseed,
	"--seed-noise:%s",&snseed,
	"--seed-spatial:%s",&ssseed,
	"--seed-photon:%s",&spseed,
	"--comment:%f",&is_comment,"(C):%f",&is_comment,
	"--verbose:%i",&is_verbose,"(V):%f",&is_verbose, 
	"*:%e",
	NULL);

 if ( i )	
  {	fprint_error("invalid command line argument near '%s'",argv[i]);
	return(1);
  }
 else if ( is_help<0 )
  {	fprint_generic_version(stdout,argv[0],"firandom",FI_RANDOM_VERSION,is_help);
	return(0);
  }
 else if ( is_help>1 )
  {	fprint_firandom_long_help(stdout);
	return(0);
  }
 else if ( is_help )
  {	fprint_firandom_usage(stdout);
	return(0);
  }

 ilp.listtype=ilisttype;
 ilp.colx--,ilp.coly--;
 if ( ilp.colmag>0 && ilp.colflux>0 )
  {	fprint_error("invalid combination of command line arguments");
	return(1);
  }
 else if ( ilp.colflux>0 )	ilp.colflux--,ilp.colmag=-1;
 else if ( ilp.colmag>0  )	ilp.colmag--,ilp.colflux=-1;
 else
  {	fprint_error("invalid combination of command line arguments");
	return(1);
  }
 if ( ilp.colx<0 || ilp.coly<0 )
  {	fprint_error("invalid column specifications for coordinates");
	return(1);
  }
 col_pix_x--,col_pix_y--;
 col_value--;

 if ( shapepar==NULL )
  {	ilp.colsh1=3,ilp.colsh2=4,ilp.colsh3=5,ilp.colsh4=6;		}
 else
  {	int	j=0;
	switch ( ilisttype )
	 {   case LISTTYPE_FEP: case LISTTYPE_SIG: case LISTTYPE_SDK:
		i=sscanf(shapepar,"%d,%d,%d",&ilp.colsh1,&ilp.colsh2,&ilp.colsh3);
		if ( i<3 ) 	j=1;
		else
		 {	ilp.colsh1--,ilp.colsh2--,ilp.colsh3--;
			if ( ilp.colsh1<0 || ilp.colsh2<0 || ilp.colsh3<0 )
		 	 {	fprint_error("invalid column specifications for shape parameters");
				return(1);
			 }
		 }
		break;
	     case LISTTYPE_SMM:
	 	i=sscanf(shapepar,"%d,%d",&ilp.colsh1,&ilp.colsh2);
		if ( i<2 ) 	j=1;
		else
		 {	ilp.colsh1--,ilp.colsh2--;
			if ( ilp.colsh1<0 || ilp.colsh2<0 )
		 	 {	fprint_error("invalid column specifications for shape parameters");
				return(1);
			 }
		 }
		break;
	     case LISTTYPE_PDS:
	 	i=sscanf(shapepar,"%d,%d,%d,%d",&ilp.colsh1,&ilp.colsh2,&ilp.colsh3,&ilp.colsh4);
		if ( i<4 ) 	j=1;
		else
		 {	ilp.colsh1--,ilp.colsh2--,ilp.colsh3--,ilp.colsh4--;
			if ( ilp.colsh1<0 || ilp.colsh2<0 || ilp.colsh3<0 || ilp.colsh4<0 )	
			 {	fprint_error("invalid column specifications for shape parameters");
				return(1);
			 }
		 }
		break;	
	 }
	if ( j )
	 {	fprint_error("invalid column specifications for shape parameters in '%s'",shapepar);
		return(1);
	 }
  }

 if ( parse_fits_data_param(fdpstring,&fdp) )
  {	fprint_error("invalid pixel data format");
	return(1);
  }

 wgthsize=3,wgthzoom=1;
 if ( weightarg != NULL )
  {	i=scanpar(weightarg,SCANPAR_DEFAULT,
		"hsize:%d",&wgthsize,
		"zoom:%d",&wgthzoom,
		NULL);
	if ( i || wgthsize<0 || wgthzoom<=0 )	
	 {	fprint_error("invalid weight specification in '%s",weightarg);
		return(1);
	 }
  }

 if ( psffile != NULL )
  {	FILE	*fr;
	if ( (fr=fopenread(psffile))==NULL )
	 {	fprint_error("unable to open input PSF file '%s'.",psffile);
		return(1);
	 }
	psf_read_fits(fr,&tpd);
	fclose(fr);
	sgp.tpd=&tpd;
  }
 else	sgp.tpd=NULL;

 if ( is_exact )
  {	sgp.dontquantize=1;
	sgp.method=1;
	sgp.is_photnoise=0;
	fdp.bitpix=-32;
	fdp.nquantizebit=0;
	bgarg=NULL;
	stddev=0.0;
  }

 if ( saseed != NULL )	snseed=spseed=ssseed=saseed;

 i=0;
 i |= set_seed(&nseed,snseed);
 i |= set_seed(&pseed,spseed);
 i |= set_seed(&sseed,ssseed);
 if ( i )	
  {	fprint_error("invalid random seed value(s)");
	return(1);
  }

 ox=0.5*(double)sx,
 oy=0.5*(double)sy,
 scale=0.5*(double)sx;

 if ( subpixelfile != NULL )
  {	FILE	*fr;
	fr=fopenread(subpixelfile);
	if ( fr==NULL )		
	 {	fprint_error("unable to open subpixel data file");
		return(1);
	 }
	read_subpixel_file(fr,&sgp.subpixeldata,&sgp.subg);
	fcloseread(fr);
	if ( subpixelfile_normalize )
		normalize_subpixeldata(sgp.subpixeldata,sgp.subg);
  }

 if ( ilist != NULL && stararg != NULL )
  {	fprint_error("ambiguous input source specifications");
	return(1);
  }

 stars=NULL;nstar=0;
 if ( ilist != NULL )
  {	FILE	*fr;
	fr=fopenread(ilist);
	if ( fr==NULL )			
	 {	fprint_error("unable to open input list file");
		return(1);
	 }
	ilp.mf0.magnitude=slp.mf0.magnitude;
	ilp.mf0.intensity=slp.mf0.intensity;
	read_input_list(fr,&stars,&nstar,&ilp);
	fcloseread(fr);
  }
 if ( stararg != NULL )
  {	int	r;
	slp.ox   =ox,
	slp.oy   =oy,
	slp.scale=scale;
	slp.sx   =sx,
	slp.sy   =sy;
	r=replace_limiters(stararg);
	if ( replace_limiters(stararg) || create_input_list(stararg,&slp,&stars,&nstar,basistype,sseed) )
	 {	fprint_error("unable to parse input source argument");
		return(1);
	 }
  }

 mask=fits_mask_create_empty(sx,sy);
 for ( i=0 ; maskblocklist != NULL && maskblocklist[i] != NULL ; i++ )
  {	if ( maskdraw_parse_and_draw(mask,sx,sy,maskblocklist[i])>0 )
	 {	fprint_warning("invalid mask block specification '%s', skipped.\n",maskblocklist[i]);	}
  }

 if ( outfile != NULL )
  {	FILE	*fw;
	fits	*img;
	char	buff[64];

	img=fits_create();
	fits_set_standard(img,NULL);
	fits_alloc_image(img,sx*zoom,sy*zoom);
	img->i.bit=fdp.bitpix;
	img->i.curr.bscale=1.0;
	img->i.curr.bzero=0.0;
	if ( fdp.is_scale )	img->i.read.bscale=fdp.bscale,img->i.read.bzero=fdp.bzero;
	else if ( img->i.bit<0 )img->i.read.bscale=1.0,img->i.read.bzero=0.0;
	else			img->i.read.bscale=1.0,img->i.read.bzero=(1<<(img->i.bit-1));
	fits_set_image_params(img);
	sprintf(buff,"Created by firandom, version: %s",FI_RANDOM_VERSION);
	fits_set_origin(img,buff,NULL);
	fits_history_export_command_line(img,"firandom",FI_RANDOM_VERSION,argc,argv);
	fits_header_export_command_line(img,"FIRANDOM",NULL,NULL,argc,argv);
	fits_set_header_double(img,"GAIN",FITS_SH_FIRST,sgp.gain,"photon/ADU value modelled");

	srand48((long int)(nseed));
	if ( create_background(&img->i,bgarg,stddev,ox,oy,scale,zoom) ) 
	 {	fprint_error("symbolic or syntax error in bacground expression");
		return(1);
	 }

	if ( ifsky != NULL )
	 {	FILE	*fr;
		fr=fopenread(ifsky);
		if ( fr==NULL )
		 {	fprint_error("unable to open input sky file");
			return(1);
		 }
		while ( ! feof(fr) )
		 {	char	*line,**cmd;
			int	n,x,y;
			double	w;
			if ( (line=freadline(fr))==NULL )
				break;
			remove_newlines_and_comments(line);
			cmd=tokenize_spaces_dyn(line);
			for ( n=0 ; cmd != NULL && cmd[n] != NULL ; )	n++;
			if ( n<=0 )
			 {	if ( cmd != NULL )	free(cmd);
				free(line);
				continue;
			 }
			if ( n<=col_pix_x || n<=col_pix_y || n<=col_value )
			 {	free(cmd);
				free(line);
				continue;
			 }
			if ( sscanf(cmd[col_pix_x],"%d",&x)<1 || 
			     sscanf(cmd[col_pix_y],"%d",&y)<1 ||
			     sscanf(cmd[col_value],"%lg",&w)<1 )
			 {	free(cmd);
				free(line);
				continue;
			 }
			/* x--,y--; */ /* LR=(0,0)! */ /* see also fiinfo.c */
			if ( 0<=x && x<sx && 0<=y && y<sy )
			 {	img->i.data[y][x]=w;
				/* fprintf(stderr,"n=%d x=%d y=%d w=%g sx=%d sy=%d\n",n,x,y,w,sx,sy); */
			 }
			free(cmd);
			free(line);
		 }
		fcloseread(fr);
	 }

	if ( sgp.is_intinelect )	divide_image(&img->i,sgp.gain);

	if ( stars != NULL && nstar>0 )
	 {	srand48((long int)(pseed));
		draw_starlist(&img->i,&sgp,stars,nstar,zoom);
	 }

	if ( ! sgp.dontquantize )
		quantize_image(&img->i);
	else if ( fdp.nquantizebit>0 )
		fits_image_quantize(&img->i,fdp.nquantizebit);

	fits_backscale(img,img->i.read.bscale,img->i.read.bzero);
	mark_integerlimited_pixels(&img->i,mask,img->i.bit,1,MASK_OVERSATURATED,MASK_OVERSATURATED);

	if ( mask != NULL && (!fits_mask_is_clean(mask,sx,sy)) )
		fits_mask_export_as_header(&img->header,1,mask,sx,sy,NULL);

	fw=fopenwrite(outfile);
	if ( fw==NULL )			
	 {	fprint_error("unable to create output file '%s'",outfile);
		return(1);
	 }
 	fits_write(fw,img); 
	fclosewrite(fw);

	fits_free(img);
  }

 if ( outmaskname != NULL )
  {	FILE	*fw;
	fits	*img;
	char	buff[64];

	img=fits_create();
	fits_set_standard(img,NULL);
	img->i.bit=fdp.bitpix;
	img->i.curr.bscale=1.0;
	img->i.curr.bzero=0.0;
	if ( fdp.is_scale )	img->i.read.bscale=fdp.bscale,img->i.read.bzero=fdp.bzero;
	else if ( img->i.bit<0 )img->i.read.bscale=1.0,img->i.read.bzero=0.0;
	else			img->i.read.bscale=1.0,img->i.read.bzero=(1<<(img->i.bit-1));
	fits_set_image_params(img);
	sprintf(buff,"Created by firandom, version: %s",FI_RANDOM_VERSION);
	fits_set_origin(img,buff,NULL);
	fits_history_export_command_line(img,"firandom",FI_RANDOM_VERSION,argc,argv);
	fits_header_export_command_line(img,"FIRANDOM",NULL,NULL,argc,argv);
	fits_set_header_double(img,"GAIN",FITS_SH_FIRST,sgp.gain,"photon/ADU value modelled");

	if ( mask != NULL )
		fits_mask_export_as_header(&img->header,1,mask,sx,sy,NULL);

	fw=fopenwrite(outmaskname);
	if ( fw==NULL )		
	 {	fprint_error("unable to create output mask file");
		return(1);
	 }
 	fits_write(fw,img); 
	fclosewrite(fw);

	fits_free(img);
  }

 if ( mask != NULL )
	fits_mask_free(mask);

 if ( oweight != NULL )
  {	FILE		*fw;
	fits		*img;
	char		buff[64];
	weightlist	wl;

	weight_draw(&wl,stars,nstar,wgthsize,wgthzoom,sgp.tpd);
	img=weight_fits_create(&wl);

	sprintf(buff,"Created by firandom, version: %s",FI_RANDOM_VERSION);
	fits_set_origin(img,buff,NULL);
	fits_header_export_command_line(img,"FIRANDOM",NULL,NULL,argc,argv);

	fw=fopenwrite(oweight);
	if ( fw==NULL )			
	 {	fprint_error("unable to create output weight file '%s'",oweight);
		return(1);
	 }
 	fits_write(fw,img); 
	fclosewrite(fw);

	fits_free(img);
  }

 if ( olist != NULL )
  {	FILE	*fw;
	fw=fopenwrite(olist);
	if ( fw==NULL )			
	 {	fprint_error("unable to create output list file '%s'",olist);
		return(1);
	 }
	if ( is_out_mag )
		write_output_list(fw,stars,nstar,olisttype,&slp.mf0,basistype);
	else
		write_output_list(fw,stars,nstar,olisttype,NULL,basistype);
	fclosewrite(fw);
  }

 return(0);
}
