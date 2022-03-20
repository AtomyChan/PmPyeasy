/*****************************************************************************/
/* star-psf.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Functions related to star PSF fitting to a FITS image. 		     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2006, Pal, A. (apal@szofi.elte.hu), part of the 'FI' package	     */
/*****************************************************************************/

#define	__STAR_PSF_C	1

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <fits/fits.h>

#include "fitsmask.h"
#include "math/fit/lmfit.h"
#include "math/poly.h"
#include "math/spline/biquad.h"
#include "math/spline/biquad-isc.h"
#include "link/linkblock.h"
#include "tensor.h"

#include "fitsh.h"
#include "common.h"
#include "stars.h"

typedef struct
 {	int	hsize;
	int	grid;
	double	**bqc;
	double	**bqcdx,
		**bqcdy;
 } dcpsf;

void psf_2d_exact_xy(void *xpnt,double *a,double *yy,double *dyda,void *p)
{
 double	x,y,x0,y0,x1,y1,x2,y2,px1,py1,px2,py2,	
	grid,offs,maxv,i0,ix,iy,is_normal;
 dcpsf  *dp=(dcpsf *)p;

 x=((ipoint *)xpnt)->x,x0=a[2],x1=x-x0,x2=x1+1.0;
 y=((ipoint *)xpnt)->y,y0=a[3],y1=y-y0,y2=y1+1.0;

 grid=(double)dp->grid;
 offs=0.5+(double)dp->hsize;
 maxv=grid*2.0*offs;

 px1=grid*(x1+offs),py1=grid*(y1+offs);
 px2=grid*(x2+offs),py2=grid*(y2+offs);

 if ( px1<0.0 )	px1=0.0;
 if ( px1>maxv)	px1=maxv;
 if ( px2<0.0 )	px2=0.0;
 if ( px2>maxv)	px2=maxv;
 if ( py1<0.0 )	py1=0.0;
 if ( py1>maxv)	py1=maxv;
 if ( py2<0.0 )	py2=0.0;
 if ( py2>maxv)	py2=maxv;
 is_normal=(px1<px2 && py1<py2);

 if ( is_normal )	i0=biquad_isc_int_rectangle(dp->bqc,px1,py1,px2,py2);
 else			i0=0.0;

 *yy=a[0]*i0+a[1];

 if ( dyda != NULL )
  {	dyda[0]=i0;
	dyda[1]=1.0;
	if ( is_normal )
	 {	ix=biquad_isc_int_rectangle(dp->bqcdx,px1,py1,px2,py2),
		iy=biquad_isc_int_rectangle(dp->bqcdy,px1,py1,px2,py2);
	 }
	else
	 {	ix=0.0,
		iy=0.0;
	 }
	dyda[2]=-ix*a[0]*grid;
	dyda[3]=-iy*a[0]*grid;
  }
}

int fit_psf(int nipoint,double *yvals,void **fitpnt,candidate *wc,
	starlocation *loc,psffit *pfp,psf *tpd)
{
 void	(*psf_abxy)(void *,double *,double *,double *,void *);
 dcpsf	dp;
 int	grid,hsize,order,bx,by,nvar;
 int	i,j,k;
 double	***btarr,**coeff,*cpoly,afp[4],lam,w,sum,ox,oy,scale;

 psf_abxy=psf_2d_exact_xy;

 grid =dp.grid =tpd->grid;
 hsize=dp.hsize=tpd->hsize;
 order=tpd->order;
 nvar=(order+1)*(order+2)/2;
 ox=tpd->ox,
 oy=tpd->oy,
 scale=tpd->scale;

 bx=grid*(2*hsize+1),
 by=grid*(2*hsize+1);

 coeff=(double  **)tensor_alloc_2d(double,bx,by);
 btarr=(double ***)tensor_alloc_3d(double,2*bx+1,2*by+1,3);
 dp.bqc  =btarr[0];
 dp.bqcdx=btarr[1];
 dp.bqcdy=btarr[2];

 cpoly=(double *)tensor_alloc_1d(double,nvar);
 sum=0.0;
 for ( i=0 ; i<by ; i++ )
  {	for ( j=0 ; j<bx ; j++ )
	 {	for ( k=0 ; k<nvar ; k++ )
		 {	cpoly[k]=tpd->coeff[k][i][j];	}
		w=eval_2d_poly(wc->cx,wc->cy,order,cpoly,ox,oy,scale);
		coeff[i][j]=w;
		sum+=w;
	 }
  }
 tensor_free(cpoly);
 w=1.0/sum;
 for ( i=0 ; i<by ; i++ )
  {	for ( j=0 ; j<bx ; j++ )
	 {	coeff[i][j] *= w; 	}
  }
		
 biquad_coeff(coeff,bx,by,dp.bqc,NULL);
 biquad_diff_x(dp.bqc,bx,by,dp.bqcdx,NULL);
 biquad_diff_y(dp.bqc,bx,by,dp.bqcdy,NULL);

 afp[0]=wc->amp*sum;
 afp[1]=wc->bg;
 afp[2]=wc->cx;
 afp[3]=wc->cy;
 lam=0.001;

 for ( i=0 ; i<pfp->iterations ; i++ )
	lam=nlm_fit_base(fitpnt,yvals,afp,NULL,psf_abxy,4,nipoint,&dp,lam,10.0);

/* lin_fit(fitpnt,yvals,afp,NULL,psf_abxy,2,nipoint,&dp,NULL); */

 /*!!*/
 for ( i=0 ; i<nipoint ; i++ )
  {	double	val;
	psf_abxy(fitpnt[i],afp,&val,NULL,&dp);
	yvals[i]=val-afp[1];
  }
/* fprintf(stderr,"bg: %g -> %g\n",wc->bg,afp[1]); */

 loc->gamp=afp[0];
 loc->gbg=afp[1];
 loc->gcx=afp[2];
 loc->gcy=afp[3];

 tensor_free(btarr);
 tensor_free(coeff);

 return(0);
}


int fit_star_psf_native(fitsimage *img,char **mask,candidate *cands,int ncand,star **rstars,int *rnstar,psffit *pfp,psf *tpd)
{
 star		*stars,*ws;
 candidate	*wc;
 int		nstar,n,ix,iy,bsize;
 int		i;
 void		**fitpnt;
 double		*yvals;
 starlocation	loc;
 ipoint		*ipoints;

 if ( rstars != NULL )  *rstars=NULL;
 if ( rnstar != NULL )  *rnstar=0;
 
 if ( ncand == 0 )      return(0);
 else if ( ncand < 0 )  return(1);

 bsize=0;
 for ( n=0 ; n<ncand ; n++ )
  {     if ( cands[n].nipoint>bsize )
                bsize=cands[n].nipoint;
  }

 fitpnt=(void  **)malloc(sizeof(void *)*bsize);
 yvals =(double *)malloc(sizeof(double)*bsize);

 stars=NULL;
 nstar=0;

 for ( n=0 ; n<ncand ; n++ )
  {	wc=&cands[n];

        if ( wc->marked )                               continue;
        if ( wc->nipoint==0 || wc->ipoints==NULL )      continue;

        ipoints=wc->ipoints;
        for ( i=0 ; i<wc->nipoint ; i++ )
         {	ix=wc->ipoints[i].x,
		iy=wc->ipoints[i].y;
		yvals[i]=img->data[iy][ix];
		fitpnt[i]=(void *)(&ipoints[i]);
         }

	i=fit_psf(wc->nipoint,yvals,fitpnt,wc,&loc,pfp,tpd);
	if ( i )	continue;	/* fit failed */

	/*!!*/
        for ( i=0 ; i<wc->nipoint ; i++ )
         {	ix=wc->ipoints[i].x,
		iy=wc->ipoints[i].y;
		img->data[iy][ix]-=yvals[i];
         }


	stars=(star *)realloc(stars,sizeof(star)*(nstar+1));
	ws=&stars[nstar];
	nstar++;

	memcpy(&ws->location,&loc,sizeof(starlocation));
	ws->shape.model=SHAPE_PSF;
	ws->shape.order=0;
	ws->shape.gs=1.0;
	ws->shape.gd=ws->shape.gk=ws->shape.gl=0.0;
	for ( i=0 ; i<MAX_DEVIATION_COEFF ; i++ )
	 {	ws->shape.mom[i]=0.0;			}

	ws->gsig=ws->gdel=ws->gkap=0.0;
	ws->gfwhm=ws->gellip=ws->gpa=0.0;

	ws->flux=ws->location.gamp;	/* flux of psf is unity */

	ws->marked=0;
	ws->cand=wc;
  }

/*
 {	FILE	*fw;
	fw=fopen("psffit-temp.fits","wb");
	fits_write(fw,img);
	fclose(fw);
 }
*/

 if ( rstars != NULL )	*rstars=stars;
 if ( rnstar != NULL )	*rnstar=nstar;
 
 return(0);
}

int fit_star_psf(fitsimage *img,char **mask,candidate *cands,int ncand,star **rstars,int *rnstar,psffit *pfp,psf *tpd)
{
 star		*stars;
 candidate	*wc;
 int		nstar,n,b,ix,iy;
 linkblock	*lbc,*wl;
 int		hsize;

 hsize=tpd->hsize;

/* fprintf(stderr,"[1]\n"); */

 lbc=(linkblock *)malloc(sizeof(linkblock)*ncand);
 for ( n=0 ; n<ncand ; n++ )
  {	wc=&cands[n];
	ix=(int)floor(wc->cx);
	iy=(int)floor(wc->cy);
	wl=&lbc[n];
	wl->x1=ix-hsize,wl->x2=ix+hsize+1;
	wl->y1=iy-hsize,wl->y2=iy+hsize+1;
	wl->flag=0;
  }
/* fprintf(stderr,"[2]\n"); */
 linkblock_connect(lbc,ncand);
/* fprintf(stderr,"[3]\n"); */

 b=0;

 stars=NULL;
 nstar=0;

 for ( n=0 ; n<ncand ; n++ )
  {	wl=&lbc[n];
	if ( wl->flag )	continue;
	while ( wl->prev != NULL )	wl=wl->prev;
	for ( ; wl != NULL ; wl=wl->next )
	 {	fprintf(stdout,"%d %d %d %d\n",wl->x1,wl->y1,wl->x2-wl->x1,0);
		fprintf(stdout,"%d %d %d %d\n",wl->x1,wl->y2,wl->x2-wl->x1,0);
		fprintf(stdout,"%d %d %d %d\n",wl->x1,wl->y1,0,wl->y2-wl->y1);
		fprintf(stdout,"%d %d %d %d\n",wl->x2,wl->y1,0,wl->y2-wl->y1);
		wl->flag=1;			}
	b++;
	fprintf(stdout,"\n");
	stars=(star *)realloc(stars,sizeof(star)*(nstar+1));
	nstar++;
  }
/* fprintf(stderr,"b=%d n=%d\n",b,n); */

 if ( rstars != NULL )	*rstars=stars;
 if ( rnstar != NULL )	*rnstar=nstar;
 
 return(0);
}

