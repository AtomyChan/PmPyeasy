/*****************************************************************************/
/* psf.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Functions related to PSF determination. The currently implemented 	     */
/* algorithms are the following:					     */
/*  - native (theoretically wrong, however, asymptotically works fine for    */
/*    wide stars, when FWHM is definitely larger than the pixel size);       */
/*  - integral (the method is somehow degenerated, the reason is not known;  */
/*    however, if an additional constraint (see: `kappa`) is applied, it     */
/*    seems to work, but the returned PSF is narrower than it should be);    */
/*  - circle (it works fine, but somehow it is also degenerated... still     */
/*    it is not known why, but the effect is negligable).		     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2004, 2006; Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>
#include <time.h>

#include <fits/fits.h>

#include "fitsh.h"

#include "fitsmask.h"
#include "math/spline/biquad.h"
#include "math/spline/biquad-isc.h"
#include "math/fit/lmfit.h"
#include "math/poly.h"
#include "math/intersec/intersec-cri.h"

#include "tensor.h"

#include "psf.h"
#include "psf-base.h"
#include "psf-determine.h"

/*****************************************************************************/

static double interpolate_subpixel(double **arr,int ix,int iy,double x0,double y0,double x1,double y1)
{
 double	ret;

 ret=arr[iy][ix]*(x1-x0)*(y1-y0);

 return(ret);
}

static double interpolate(double **arr,double x0,double y0,double x1,double y1)
{
 double	w,ret;
 int	sig,ix0,ix1,iy0,iy1,x,y;
 sig=1;
 if ( x0>x1 )	sig=-sig,w=x0,x0=x1,x1=w;
 if ( y0>y1 )	sig=-sig,w=y0,y0=y1,y1=w;
 ix0=(int)floor(x0),iy0=(int)floor(y0);
 ix1=(int)floor(x1),iy1=(int)floor(y1);
 x0-=(double)ix0,y0-=(double)iy0;
 x1-=(double)ix1,y1-=(double)iy1;
 if ( ix0==ix1 )
  {	if ( iy0==iy1 )
		ret=interpolate_subpixel(arr,ix0,iy0,x0,y0,x1,y1);
	else
	 {	ret=interpolate_subpixel(arr,ix0,iy0,x0,y0,x1,1.0);
		for ( y=iy0+1 ; y<iy1 ; y++ )
		 	ret+=interpolate_subpixel(arr,ix0,y,x0,0.0,x1,1.0);	
		if ( y1>0.0 )
			ret+=interpolate_subpixel(arr,ix0,iy1,x0,0.0,x1,y1);
	 }
  }
 else
  {	if ( iy0==iy1 )
	 {	ret=interpolate_subpixel(arr,ix0,iy0,x0,y0,1.0,y1);
		for ( x=ix0+1 ; x<ix1 ; x++ )
			ret+=interpolate_subpixel(arr,x,iy0,0.0,y0,1.0,y1);	
		if ( x1>0.0 )
			ret+=interpolate_subpixel(arr,ix1,iy0,0.0,y0,x1,y1);
	 }
	else
	 {	ret=interpolate_subpixel(arr,ix0,iy0,x0,y0,1.0,1.0);
		for ( x=ix0+1 ; x<ix1 ; x++ )
			ret+=interpolate_subpixel(arr,x,iy0,0.0,y0,1.0,1.0);
		if ( x1>0.0 )
			ret+=interpolate_subpixel(arr,ix1,iy0,0.0,y0,x1,1.0);
		for ( y=iy0+1 ; y<iy1 ; y++ )
		 {	ret+=interpolate_subpixel(arr,ix0,y,x0,0.0,1.0,1.0);
			for ( x=ix0+1 ; x<ix1 ; x++ )
				ret+=arr[y][x];
			if ( x1>0.0 )
				ret+=interpolate_subpixel(arr,ix1,y,0.0,0.0,x1,1.0);
		 }
		if ( y1>0.0 )
		 {	ret+=interpolate_subpixel(arr,ix0,iy1,x0,0.0,1.0,y1);
			for ( x=ix0+1 ; x<ix1 ; x++ )
				ret+=interpolate_subpixel(arr,x,iy1,0.0,0.0,1.0,y1);
			if ( x1>0.0 )
				ret+=interpolate_subpixel(arr,ix1,iy1,0.0,0.0,x1,y1);
		 }
	 }
  }
 if ( sig>0 )	return(ret);
 else		return(-ret);
}

int psf_determine_native(fitsimage *img,char **mask,
	psfcandidate *cands,int ncand,int is_subtracted,
	int hsize,int grid,int order,psf *out,int use_biquad)
{
 int		i,j,k,l,n,nvar,sx,sy,t;
 int		fsize,wwd;
 double		***psfstack,*pmonom,**amatrix,*bvector;
 double		ox,oy,scale,igrid;
 double		**bqc;
 psfcandidate	*wc;
 ipoint		*wi;
 double		dx,dy,x1,y1,x2,y2,z,gbg,gamp,weight;
 int		ix1,ix2,iy1,iy2;
 
 if ( img==NULL )	return(1);
 if ( img->data==NULL )	return(1);
 sx=img->sx,sy=img->sy;

 if ( use_biquad )
  {	bqc=(double **)tensor_alloc_2d(double,2*sx+1,2*sy+1);
	biquad_coeff(img->data,sx,sy,bqc,NULL);
  }
 else
	bqc=NULL;

 fsize=2*hsize+1;
 wwd=fsize*grid;

 nvar=(order+1)*(order+2)/2;

 psfstack=(double ***)tensor_alloc_3d(double,wwd,wwd,nvar);
 
 pmonom=vector_alloc(nvar);
 bvector=vector_alloc(nvar);
 amatrix=matrix_alloc(nvar);

 ox=(double)sx*0.5,
 oy=(double)sy*0.5,
 scale=(double)sx;

 igrid=1.0/(double)grid;

 for ( i=0 ; i<wwd ; i++ )
  { for ( j=0 ; j<wwd ; j++ )
     {
	dx=igrid*(double)j;
	dy=igrid*(double)i;

	for ( k=0 ; k<nvar ; k++ )
	 {	for ( l=0 ; l<nvar ; l++ )
	 	 {	amatrix[k][l]=0.0;	}
		bvector[k]=0.0;
	 }

	t=0;	
	for ( n=0 ; n<ncand ; n++ )
	 {
		wc=&cands[n];
		if ( wc->nipoint <= 0 || wc->ipoints==NULL )	continue;

		x1=wc->cx+dx-((double)hsize+0.5),x2=x1+igrid;
		y1=wc->cy+dy-((double)hsize+0.5),y2=y1+igrid;
		ix1=(int)x1,iy1=(int)y1;
		ix2=(int)x2,iy2=(int)y2;
		if ( ix1<0 || ix2>=sx || iy1<0 || iy2>=sy )	continue;

		k=0;
		for ( l=0,wi=wc->ipoints ; l<wc->nipoint && k<4 ; l++,wi++ )
		 {	if ( wi->x==ix1 && wi->y==iy1 )	k++;
			if ( wi->x==ix1 && wi->y==iy2 )	k++;
			if ( wi->x==ix2 && wi->y==iy1 )	k++;
			if ( wi->x==ix2 && wi->y==iy2 )	k++;
		 }
		if ( k<4 )	continue;

		gbg =wc->bg;
		gamp=wc->amp;

		eval_2d_monoms(wc->cx,wc->cy,order,pmonom,ox,oy,scale);

		if ( bqc != NULL )	z=biquad_isc_int_rectangle(bqc,x1,y1,x2,y2)*(double)(grid*grid);
		else			z=interpolate(img->data,x1,y1,x2,y2)*(double)(grid*grid);

		z=(z-gbg)/gamp;
		weight=gamp;

		for ( k=0 ; k<nvar ; k++ )
		 {	for ( l=0 ; l<nvar ; l++ )
			 {	amatrix[k][l]+=pmonom[k]*pmonom[l]*weight;	}
			bvector[k]+=pmonom[k]*z*weight;
		 }
		t++;
	 }

	if ( t<nvar )
	 {	for ( k=0 ; k<nvar ; k++ )
		 {	psfstack[k][i][j]=0.0;		}
	 }
	else
	 {	solve_gauss(amatrix,bvector,nvar);
		for ( k=0 ; k<nvar ; k++ )
		 {	psfstack[k][i][j]=bvector[k];		}
	 }
		
     }
  }

 matrix_free(amatrix);
 vector_free(bvector);
 vector_free(pmonom);

 if ( bqc != NULL )	tensor_free(bqc);

 out->coeff=psfstack;
 out->hsize=hsize;
 out->grid=grid;
 out->order=order;
 out->ox=ox,out->oy=oy;
 out->scale=scale;

 return(0);
}

/*****************************************************************************/

typedef struct
 {	int	indx;
	double	area;
 } idlist;

int psf_determine_integral(fitsimage *img,char **mask,
	psfcandidate *cands,int ncand,int is_subtracted,
	int hsize,int grid,int order,psf *out,double kappa)
{
 int		i,j,k,l,n,m,p,q,nvar,nfv,t,sx,sy,ttot;
 int		fsize,wwd;
 double		***psfstack,*pmonom,**amatrix,*bvector;
 double		ox,oy,scale,dgrid,poffs,kfactor,sumweight,flux,tflux,tamp;
 psfcandidate	*wc;
 ipoint		*wi;
 double		x1,y1,x2,y2,wx,wy,z,gbg,gamp,weight,area,w,wk,wl;
 int		ix1,ix2,iy1,iy2;
 idlist		*ils;
 int		nil;

 /*order=0;*/	/* okay, it is fixed for now. */

 fprintf(stderr,"[1] clock = %12.3f\n",(double)clock()/(double)(CLOCKS_PER_SEC));

 if ( img==NULL )	return(1);
 if ( img->data==NULL )	return(1);
 sx=img->sx,sy=img->sy;

 fsize=2*hsize+1;
 wwd=fsize*grid;

 nvar=(order+1)*(order+2)/2;

 psfstack=(double ***)tensor_alloc_3d(double,wwd,wwd,nvar);
 nfv=wwd*wwd*nvar;	/* lot of variables to fit... */
 
 pmonom =vector_alloc(nvar);
 bvector=vector_alloc(nfv);
 amatrix=matrix_alloc(nfv);

 ils=(idlist *)tensor_alloc_1d(idlist,wwd*wwd);

 ox=(double)sx*0.5,
 oy=(double)sy*0.5,
 scale=(double)sx;

 for ( i=0 ; i<nfv ; i++ )
  {	for ( j=0 ; j<nfv ; j++ )
	 {	amatrix[i][j]=0.0;		}
	bvector[i]=0.0;
  }

 for ( n=0 ; n<ncand ; n++ )
  {	wc=&cands[n];
	if ( wc->yvals==NULL )
	 {	is_subtracted=0;
		break;
	 }
  }

 if ( is_subtracted )
  {	for ( n=0 ; n<ncand ; n++ )
	 {	wc=&cands[n];
		for ( m=0 ; m<wc->nipoint ; m++ )
		 {	wi=&wc->ipoints[m];
			if ( mask != NULL && mask[wi->y][wi->x] )	continue;
			if ( mask != NULL && mask[wi->y][wi->x] )	continue;
			img->data[wi->y][wi->x] -= wc->yvals[m];
		 }
	 }
  }

 dgrid=(double)grid;
 poffs=(double)hsize+0.5;

 ttot=0; sumweight=0.0;
 tflux=tamp=0.0;
 for ( n=0 ; n<ncand ; n++ )
  {	wc=&cands[n];
	if ( wc->nipoint<=0 || wc->ipoints==NULL )	continue;
	gbg =wc->bg,

	/*gamp=wc->amp;*/ /* obsoleted */

	flux=0.0;
	for ( m=0 ; m<wc->nipoint ; m++ )
	 {	wi=&wc->ipoints[m];
		if ( mask != NULL && mask[wi->y][wi->x] )	continue;
		flux+=img->data[wi->y][wi->x]-gbg;
		if ( is_subtracted )
		 {	flux+=wc->yvals[m];		}
	 }

	gamp=flux;
	if ( gamp<=0.0 )	continue;

	tflux+=flux;
	tamp +=wc->amp;

	weight=gamp;

	eval_2d_monoms(wc->cx,wc->cy,order,pmonom,ox,oy,scale);

	for ( m=0 ; m<wc->nipoint ; m++ )
	 {	wi=&wc->ipoints[m];
		if ( mask != NULL && mask[wi->y][wi->x] )	continue;

		z=img->data[wi->y][wi->x];
		if ( is_subtracted )
		 {	z+=wc->yvals[m];	}
		
		z=(z-gbg)/gamp;
		/*fprintf(stdout,"%d %d %g\n",wi->x,wi->y,z);*/

		x1=dgrid*(((double)wi->x)-(wc->cx)+poffs),x2=x1+dgrid;
		y1=dgrid*(((double)wi->y)-(wc->cy)+poffs),y2=y1+dgrid;
		if ( x1<0.0 )	x1=0.0;
		if ( x1>wwd )	x1=wwd;
		if ( x2<0.0 )	x2=0.0;
		if ( x2>wwd )	x2=wwd;
		if ( y1<0.0 )	y1=0.0;
		if ( y1>wwd )	y1=wwd;
		if ( y2<0.0 )	y2=0.0;
		if ( y2>wwd )	y2=wwd;
		if ( x1==x2 || y1==y2 )	continue;
		
		ix1=(int)x1,iy1=(int)y1;
		ix2=(int)x2,iy2=(int)y2;
		nil=0;area=0.0;
		for ( i=iy1 ; i<=iy2 ; i++ )
		 {	if ( iy1==iy2 )		wy=y2-y1;
			else if ( i==iy1 )	wy=1.0-(y1-(double)iy1);
			else if ( i==iy2 )	wy=(y2-(double)iy2);
			else			wy=1.0;
			if ( wy<=0.0 )		continue;
			for ( j=ix1 ; j<=ix2 ; j++ )
			 {	if ( ix1==ix2 )		wx=x2-x1;
				else if ( j==ix1 )	wx=1.0-(x1-(double)ix1);
				else if ( j==ix2 )	wx=(x2-(double)ix2);
				else			wx=1.0;
				if ( wx<=0.0 )		continue;
				t=i*wwd+j;
				w=wx*wy;
				area+=w;
				ils[nil].indx=t;
				ils[nil].area=w;
				nil++;
			 }
		 }

		if ( nvar>1 )
		 {	t=wwd*wwd;
			for ( i=0 ; i<nil ; i++ )
			 {  k =ils[i].indx,
			    wk=ils[i].area;
			    for ( j=0 ; j<nil ; j++ )
			     {	l =ils[j].indx,
				wl=ils[j].area;
				for ( p=0 ; p<nvar ; p++ )
				 {   for ( q=0 ; q<nvar ; q++ )
					amatrix[t*p+k][t*q+l]+=pmonom[p]*pmonom[q]*wk*wl*weight;
				 }	
			     }
			    for ( p=0 ; p<nvar ; p++ )
				bvector[t*p+k]+=pmonom[p]*wk*z*area*weight;
			 }	
		 }
		else
		 {	for ( i=0 ; i<nil ; i++ )
			 {  k =ils[i].indx,
			    wk=ils[i].area;
			    for ( j=0 ; j<nil ; j++ )
			     {	l =ils[j].indx,
				wl=ils[j].area;
				amatrix[k][l]+=wk*wl*weight;
			     }
			    bvector[k]+=wk*z*area*weight;
			 }
		 }

		ttot++;
		sumweight+=weight;
		/*fprintf(stderr,"nil=%d, area = %g\n",nil,area);*/
 	 }
  }

 fprintf(stderr,"[2] clock = %12.3f\n",(double)clock()/(double)(CLOCKS_PER_SEC));

 kfactor=kappa*sumweight;
 if ( kfactor>0.0 )
  {	int	t0,t1,t2;
	for ( k=0 ; k<nvar ; k++ )
	 {	t0=k*wwd*wwd;
		for ( i=0 ; i<wwd ; i++ )
		 {	for ( j=0 ; j<wwd-1 ; j++ )
			 {	t1=t0+(i*wwd+j),
				t2=t1+1;
				amatrix[t1][t1]+=kfactor;
				amatrix[t2][t1]-=kfactor;
				amatrix[t1][t2]-=kfactor;
				amatrix[t2][t2]+=kfactor;
			 }
		 }
		for ( j=0 ; j<wwd ; j++ )
		 {	for ( i=0 ; i<wwd-1 ; i++ )
			 {	t1=t0+(i*wwd+j),
				t2=t1+wwd;
				amatrix[t1][t1]+=kfactor;
				amatrix[t2][t1]-=kfactor;
				amatrix[t1][t2]-=kfactor;
				amatrix[t2][t2]+=kfactor;
			 }
		 }
	 }
  }

 fprintf(stderr,"[3] clock = %12.3f\n",(double)clock()/(double)(CLOCKS_PER_SEC));

/* fprintf(stderr,"[solve - 1 - ttot=%d, n=%d]\n",ttot,nfv); */
 i=solve_gauss(amatrix,bvector,nfv);
/* fprintf(stderr,"[solve - 2 ret=%d]\n",i); */

 fprintf(stderr,"[4] clock = %12.3f\n",(double)clock()/(double)(CLOCKS_PER_SEC));

 t=0;
 wx=tflux/tamp;
 for ( k=0 ; k<nvar ; k++ )
  {	for ( i=0 ; i<wwd ; i++ )
	 {	l=wwd*(k*wwd+i);
		for ( j=0 ; j<wwd ; j++,t++ )
		 {	psfstack[k][i][j]=bvector[l+j]*wx;		}
	 }
  }

 tensor_free(ils);

 matrix_free(amatrix);
 vector_free(bvector);
 vector_free(pmonom);

 if ( is_subtracted )
  {	for ( n=0 ; n<ncand ; n++ )
	 {	wc=&cands[n];
		for ( m=0 ; m<wc->nipoint ; m++ )
		 {	wi=&wc->ipoints[m];
			if ( mask != NULL && mask[wi->y][wi->x] )	continue;
			if ( mask != NULL && mask[wi->y][wi->x] )	continue;
			img->data[wi->y][wi->x] += wc->yvals[m];
		 }
	 }
  }

 out->coeff=psfstack;
 out->hsize=hsize;
 out->grid=grid;
 out->order=order;
 out->ox=ox,out->oy=oy;
 out->scale=scale;

 fprintf(stderr,"[5] clock = %12.3f\n",(double)clock()/(double)(CLOCKS_PER_SEC));

 return(0);
}

/*****************************************************************************/

typedef struct
 {	double	rmax;
	int	order;

	double	weight;

	int	nvar,indx,icount;
 } ring;

static int poly_to_fourier(int order,double *poly,double *four)
{
 int	o,k,f,s,b;

 if ( order<0 )	return(1);

 four[0]=poly[0];
 poly++,four++;

 for ( o=1 ; o<=order ; o++ )
  {	four[0]=0.0;
	four[1]=0.0;
	for ( k=0,f=0,s=1,b=1 ; k<=o ; k++ )
	 {	four[f]+=s*b*poly[k];
		f=1-f;
		if ( f==0 )	s=-s;
		b=b*(o-k)/(k+1);
	 }
	four+=2;
	poly+=o+1;
  }
 
 return(0);
}

int psf_determine_circle(fitsimage *img,char **mask,
	psfcandidate *cands,int ncand,int is_subtracted,
	int hsize,int grid,int order,psf *out,double circwd,int circorder)
{
 int		c,i,j,k,l,n,m,p,q,nvar,nfv,t,sx,sy,circnvar,circpvar;
 int		fsize,wwd,nil;
 double		***psfstack,*pmonom,**amatrix,*bvector,**crints,
		*crpoly,*crfour;
 double		ox,oy,scale,dgrid,igrid,poffs,flux,tflux,tamp;
 psfcandidate	*wc;
 ipoint		*wi;
 double		x1,y1,x2,y2,z,gbg,gamp,weight,area,wk,wl,wx,cwg,cr;
 idlist		*ils;
 ring		*rings,*rw;
 int		nring,nrfvar,mxcorder;

 order=0;		/* okay, it is fixed for now.		*/
 /*circorder=0;*/	/* like this, which won't be trivial...	*/

 if ( img==NULL || img->data==NULL )	return(1);
 sx=img->sx,sy=img->sy;

 fsize=2*hsize+1;
 wwd=fsize*grid;

 nring=(int)(((double)hsize+0.5)/circwd);/* number of circle rings    	     */

 rings=(ring *)malloc(sizeof(ring)*nring);
 for ( i=0 ; i<nring ; i++ )
  {	rings[i].rmax=circwd*(double)(i+1);
	/*rings[i].weight=(double)((i+1)*(i+1));*/
	rings[i].weight=1.0;
	rings[i].order=circorder;
  }

 nrfvar=0;mxcorder=0;
 for ( i=0 ; i<nring ; i++ )
  {	rings[i].nvar=1+2*rings[i].order;
	rings[i].indx=nrfvar;
	nrfvar+=rings[i].nvar;		/* total number of Fourier variables */
	if ( rings[i].order>mxcorder )	mxcorder=rings[i].order;
  }

 nvar=(order+1)*(order+2)/2;		/* spatial variation of PSF	     */
 nfv=nvar*nrfvar;			/* total number of variables to fit  */

 circpvar=(mxcorder+1)*(mxcorder+2)/2;
 circnvar=1+2*mxcorder;
 
 pmonom =vector_alloc(nvar);
 bvector=vector_alloc(nfv);
 amatrix=matrix_alloc(nfv);
 crints =(double **)tensor_alloc_2d(double,nring,circnvar);
 crpoly =(double *)tensor_alloc_1d(double,circpvar);
 crfour =(double *)tensor_alloc_1d(double,circnvar);

 ils=(idlist *)tensor_alloc_1d(idlist,nrfvar);

 ox=(double)sx*0.5,
 oy=(double)sy*0.5,
 scale=(double)sx;

 tflux=tamp=0.0;

 if ( is_subtracted )
  {	for ( n=0 ; n<ncand ; n++ )
	 {	wc=&cands[n];
		for ( m=0 ; m<wc->nipoint ; m++ )
		 {	wi=&wc->ipoints[m];
			if ( mask != NULL && mask[wi->y][wi->x] )	continue;
			if ( mask != NULL && mask[wi->y][wi->x] )	continue;
			img->data[wi->y][wi->x] -= wc->yvals[m];
		 }
	 }
  }


 for ( i=0 ; i<nfv ; i++ )
  {	for ( j=0 ; j<nfv ; j++ )
	 {	amatrix[i][j]=0.0;		}
	bvector[i]=0.0;
  }

 dgrid=(double)grid;
 igrid=1.0/dgrid;
 poffs=(double)hsize+0.5;

 for ( n=0 ; n<ncand ; n++ )
  {	wc=&cands[n];
	if ( wc->nipoint<=0 || wc->ipoints==NULL )	continue;
	gbg =wc->bg;

	/*gamp=wc->amp;*/ /* obsoleted */

	flux=0.0;
	for ( m=0 ; m<wc->nipoint ; m++ )
	 {	wi=&wc->ipoints[m];
		if ( mask != NULL && mask[wi->y][wi->x] )	continue;
		flux+=img->data[wi->y][wi->x]-gbg;
		if ( is_subtracted )
		 {	flux+=wc->yvals[m];		}
	 }

	gamp=flux;
	if ( gamp<=0.0 )	continue;

	tflux+=flux;
	tamp +=wc->amp;

	weight=gamp;

	eval_2d_monoms(wc->cx,wc->cy,order,pmonom,ox,oy,scale);

	for ( m=0 ; m<wc->nipoint ; m++ )
	 {	wi=&wc->ipoints[m];
		if ( mask != NULL && mask[wi->y][wi->x] )	continue;

		z=img->data[wi->y][wi->x];
		if ( is_subtracted )
		 {	z+=wc->yvals[m];	}

		z=(z-gbg)/gamp;

		x1=(((double)wi->x)-(wc->cx)+poffs),x2=x1+1.0;
		y1=(((double)wi->y)-(wc->cy)+poffs),y2=y1+1.0;

		if ( x1<0.0 )	x1=0.0;
		if ( x1>fsize )	x1=fsize;
		if ( x2<0.0 )	x2=0.0;
		if ( x2>fsize )	x2=fsize;
		if ( y1<0.0 )	y1=0.0;
		if ( y1>fsize )	y1=fsize;
		if ( y2<0.0 )	y2=0.0;
		if ( y2>fsize )	y2=fsize;

		if ( x1==x2 || y1==y2 )	continue;

		x1-=poffs,x2-=poffs;
		y1-=poffs,y2-=poffs;

		for ( c=0 ; c<nring ; c++ )
		 {	cr=rings[c].rmax;
			intersec_cri_integrate_monoms(x1,y1,x2-x1,y2-y1,cr,crpoly,circorder);
			poly_to_fourier(circorder,crpoly,crfour);
			for ( p=0 ; p<circnvar ; p++ )
			 {	crints[p][c]=crfour[p];		}
		 }
		for ( p=0 ; p<circnvar ; p++ )
		 {	for ( c=nring-2 ; c>=0 ; c-- )
			 {	crints[p][c+1]-=crints[p][c];		}
		 }

		nil=0;wx=0.0;
		for ( c=0 ; c<nring ; c++ )
		 {	rw=&rings[c];
			for ( p=0 ; p<rw->nvar ; p++ )
			 {	if ( crints[p][c] == 0.0 )	continue;
				ils[nil].indx=rw->indx+p;
				ils[nil].area=crints[p][c];
				nil++;
				if ( rw->weight>wx )	wx=rw->weight;
			 }
		 }

		cwg=weight*wx;
		
		area=(x2-x1)*(y2-y1);

		if ( nvar>1 )
		 {	for ( i=0 ; i<nil ; i++ )
			 {	k =ils[i].indx,
				wk=ils[i].area;
				for ( j=0 ; j<nil ; j++ )
				 {	l =ils[j].indx,
					wl=ils[j].area;
					for ( p=0 ; p<nvar ; p++ )
					 {  for ( q=0 ; q<nvar ; q++ )
					  	amatrix[nrfvar*p+k][nrfvar*q+l]+=pmonom[p]*pmonom[q]*wk*wl*cwg;
					 }
				 }
				for ( p=0 ; p<nvar ; p++ )
					bvector[nrfvar*p+k]+=pmonom[p]*wk*z*area*cwg;
			 }
			 
	 	 }
		else
		 {	for ( i=0 ; i<nil ; i++ )
			 {	k =ils[i].indx,
				wk=ils[i].area;
				for ( j=0 ; j<nil ; j++ )
				 {	l =ils[j].indx,
					wl=ils[j].area;
					amatrix[k][l]+=wk*wl*cwg;
				 }
				bvector[k]+=wk*z*area*cwg;
			 }
		 }
 	 }
  }

 i=solve_gauss(amatrix,bvector,nfv);

 psfstack=(double ***)tensor_alloc_3d(double,wwd,wwd,nvar);

 t=0;
 wx=dgrid*dgrid*tflux/tamp;	/* murphy-factor ;] */

 for ( k=0 ; k<nvar ; k++ )
  {	for ( i=0 ; i<wwd ; i++ )
	 {	for ( j=0 ; j<wwd ; j++,t++ )
		 {	x1=(igrid*(double)j-poffs),x2=x1+igrid;
			y1=(igrid*(double)i-poffs),y2=y1+igrid;

			for ( c=0 ; c<nring ; c++ )
			 {	cr=rings[c].rmax;
				intersec_cri_integrate_monoms(x1,y1,x2-x1,y2-y1,cr,crpoly,circorder);
				poly_to_fourier(circorder,crpoly,crfour);
				for ( p=0 ; p<circnvar ; p++ )
				 {	crints[p][c]=crfour[p];		}
			 }
			for ( p=0 ; p<circnvar ; p++ )
			 {	for ( c=nring-2 ; c>=0 ; c-- )
				 {	crints[p][c+1]-=crints[p][c];		}
			 }

			area=0.0;
			for ( c=0 ; c<nring ; c++ )
			 {	rw=&rings[c];
				for ( p=0 ; p<rw->nvar ; p++ )
				 {	area+=crints[p][c]*bvector[k*nrfvar+rw->indx+p];		}
			 }

			psfstack[k][i][j]=area*wx;
		 }
	 }
  }

 tensor_free(ils);

 tensor_free(crfour);
 tensor_free(crpoly);
 tensor_free(crints);
 matrix_free(amatrix);
 vector_free(bvector);
 vector_free(pmonom);

 free(rings);

 if ( is_subtracted )
  {	for ( n=0 ; n<ncand ; n++ )
	 {	wc=&cands[n];
		for ( m=0 ; m<wc->nipoint ; m++ )
		 {	wi=&wc->ipoints[m];
			if ( mask != NULL && mask[wi->y][wi->x] )	continue;
			if ( mask != NULL && mask[wi->y][wi->x] )	continue;
			img->data[wi->y][wi->x] += wc->yvals[m];
		 }
	 }
  }


 out->coeff=psfstack;
 out->hsize=hsize;
 out->grid=grid;
 out->order=order;
 out->ox=ox,out->oy=oy;
 out->scale=scale;

 return(0);
}

/*****************************************************************************/

int psf_determine(fitsimage *img,char **mask,
	psfcandidate *cands,int ncand,int is_subtracted,
	psfdetermine *pd,psf *p)
{
 int	hsize,grid,order,ret;
 
 if ( pd==NULL )	return(-1);
 hsize=pd->hsize;
 grid =pd->grid;
 order=pd->order;

 if ( hsize<0 || grid<=0 || order<0 )	return(-1);

 switch ( pd->type )
  {  case PSF_DET_NATIVE:
	ret=psf_determine_native(img,mask,cands,ncand,is_subtracted,
		hsize,grid,order,p,
		pd->param.native.use_biquad);
	break;
     case PSF_DET_INTEGRAL:
	ret=psf_determine_integral(img,mask,cands,ncand,is_subtracted,
		hsize,grid,order,p,
		pd->param.integral.kappa);
	break;
     case PSF_DET_CIRCLE:
	ret=psf_determine_circle(img,mask,cands,ncand,is_subtracted,
		hsize,grid,order,p,
		pd->param.circle.width,pd->param.circle.order);
	break;
     default:
	ret=-1;
  }

 if ( ret )	return(ret);

 if ( pd->is_symmetrize )	psf_symmetrize(p);

 return(0);
}

/*****************************************************************************/

int drawback_psf(ipoint *ipoints,int nipoint,double *yvals,
	double x0,double y0,double amp,psf *p,double mul)
{
 static double	**bqc=NULL,**coeff=NULL,*cpoly=NULL;
 static int	abx=0,aby=0,anvar=0;
 int		i,j,k,nvar,bx,by;
 double		x1,y1,x2,y2,px,py,afm,afo,sum;

 if ( ipoints==NULL || yvals==NULL )
	return(-1);
 if ( p==NULL )
	return(-1);

 nvar=(p->order+1)*(p->order+2)/2;
 bx=p->grid*(2*p->hsize+1);
 by=p->grid*(2*p->hsize+1);

 if ( bx>abx || by>aby )
  {	if ( bqc != NULL )	tensor_free(bqc); 
	if ( coeff != NULL )	tensor_free(coeff); 
	coeff=(double **)tensor_alloc_2d(double,bx,by);
	bqc  =(double **)tensor_alloc_2d(double,2*bx+1,2*by+1);
	abx=bx,aby=by;
  }
 if ( nvar>anvar )
  {	if ( cpoly != NULL )	tensor_free(cpoly); 
	cpoly=(double *)tensor_alloc_1d(double,nvar);
	anvar=nvar;
  }

 px=x0,
 py=y0;
 sum=0.0;
 for ( i=0 ; i<by ; i++ )
  {	for ( j=0 ; j<bx ; j++ )
	 {	for ( k=0 ; k<nvar ; k++ )
		 {	cpoly[k]=p->coeff[k][i][j];		}
		coeff[i][j]=eval_2d_poly(px,py,p->order,cpoly,p->ox,p->oy,p->scale);
		sum+=coeff[i][j];
	 }
  }
 amp=mul*amp/sum;
 for ( i=0 ; i<by ; i++ )
  {	for ( j=0 ; j<bx ; j++ )
	 {	coeff[i][j] *= amp;	}
  }

 biquad_coeff(coeff,bx,by,bqc,NULL);

 afm=(double)p->grid;
 afo=(double)p->grid*(0.5+(double)p->hsize);
 for ( i=0 ; i<nipoint ; i++ )
  {	x1=(double)(ipoints[i].x+0)-x0;
	x2=(double)(ipoints[i].x+1)-x0;
	y1=(double)(ipoints[i].y+0)-y0;
	y2=(double)(ipoints[i].y+1)-y0;

	x1=afm*x1+afo;
	x2=afm*x2+afo;
	y1=afm*y1+afo;
	y2=afm*y2+afo;
	if ( x1<0.0 )	x1=0.0;
	if ( x1>=bx )	x1=bx;
	if ( x2<0.0 )	x2=0.0;
	if ( x2>=bx )	x2=bx;
	if ( y1<0.0 )	y1=0.0;
	if ( y1>=by )	y1=by;
	if ( y2<0.0 )	y2=0.0;
	if ( y2>=by )	y2=by;
	if ( x1<x2 && y1<y2 )
		yvals[i] += biquad_isc_int_rectangle(bqc,x1,y1,x2,y2);
  }

 return(0); 
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int psf_bgamp_fit(fitsimage *img,char **mask,
	psfcandidate *cands,int ncand,int is_subtracted,psf *p)
{
 int		n,m;
 psfcandidate	*wc;
 ipoint		*wi;
 double		*yvals,*pvals,maa,mab,mbb,va,vb,det,x,y,a,b;
 int		nyval,nipoint;

 nyval=0;
 for ( n=0 ; n<ncand ; n++ )
  {	wc=&cands[n];
	if ( wc->yvals==NULL )
		is_subtracted=0;
	else if ( wc->nipoint>nyval )
		nyval=wc->nipoint;
  }

 yvals=(double *)malloc(sizeof(double)*nyval);
 pvals=(double *)malloc(sizeof(double)*nyval);

 if ( is_subtracted )
  {	for ( n=0 ; n<ncand ; n++ )
	 {	wc=&cands[n];
		for ( m=0 ; m<wc->nipoint ; m++ )
		 {	wi=&wc->ipoints[m];
			if ( mask != NULL && mask[wi->y][wi->x] )	continue;
			if ( mask != NULL && mask[wi->y][wi->x] )	continue;
			img->data[wi->y][wi->x] -= wc->yvals[m];
		 }
	 }
  }

 for ( n=0 ; n<ncand ; n++ )
  {	wc=&cands[n];
	nipoint=wc->nipoint;

	if ( is_subtracted )
	 {	for ( m=0 ; m<nipoint ; m++ )
		 {	wi=&wc->ipoints[m];
			if ( mask != NULL && mask[wi->y][wi->x] )	continue;
			yvals[m] = img->data[wi->y][wi->x]+wc -> yvals[m];
		 }
	 }
	else
	 {	for ( m=0 ; m<nipoint ; m++ )
		 {	wi=&wc->ipoints[m];
			if ( mask != NULL && mask[wi->y][wi->x] )	continue;
			yvals[m] = img->data[wi->y][wi->x];
		 }
	 }

	for ( m=0 ; m<nipoint ; m++ )	pvals[m]=0.0;
	drawback_psf(wc->ipoints,nipoint,pvals,wc->cx,wc->cy,+1.0,p,+1.0);

	maa=mab=mbb=0.0;
	va=vb=0.0;
	for ( m=0 ; m<nipoint ; m++ )
	 {	x=pvals[m];
		y=yvals[m];
		maa+=x*x;
		mab+=x;
		mbb+=1.0;
		va+=x*y;
		vb+=y;
	 }
	det=maa*mbb-mab*mab;
	a=(va*mbb-vb*mab)/det;
	b=(maa*vb-mab*va)/det;
	
	wc->bg =b;
	wc->amp=a;
  }
 
 if ( is_subtracted )
  {	for ( n=0 ; n<ncand ; n++ )
	 {	wc=&cands[n];
		for ( m=0 ; m<wc->nipoint ; m++ )
		 {	wi=&wc->ipoints[m];
			if ( mask != NULL && mask[wi->y][wi->x] )	continue;
			if ( mask != NULL && mask[wi->y][wi->x] )	continue;
			img->data[wi->y][wi->x] += wc->yvals[m];
		 }
	 }
  }

 free(pvals);
 free(yvals);

 return(0);
}

/*****************************************************************************/

