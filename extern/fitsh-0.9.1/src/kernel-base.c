/*****************************************************************************/
/* kernel-base.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Base functions for manipulating convolution kernels.			     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include <fits/fits.h>

#include "fi.h"

#include "math/fit/lmfit.h"
#include "math/poly.h"
#include "io/tokenize.h"
#include "io/iof.h"
#include "tensor.h"

#include "kernel.h"

int	kernel_verbose=2;
#define	KERNEL_LOG	1
                                            
/*****************************************************************************/

int kernel_set_verbosity(int vlevel)
{
 if ( vlevel>0 )	kernel_verbose=vlevel;
 else			kernel_verbose=0;
 return(0);
}

/*****************************************************************************/

/* hermite():
   Evaluates H_n(x). (Currently, it's obsolente...)			     */
double hermite(int n,double x)
{
 double	hp,hpp,h;
 int	k;
 switch ( n )
  {	case 0 :	return(1.0);
	case 1 :	return(2.0*x);
	case 2 :	return(4.0*x*x-2.0);
	case 3 :	return(x*(8.0*x*x-12.0));
	default :
		if ( n<0 )	return(0.0);
		h=hp=2.0*x,hpp=1.0;
		for ( k=1 ; k<n ; k++ )
		 {	h=2.0*(x*hp-(double)k*hpp);
			hpp=hp,hp=h;
		 }	
		return(h);
  }	
}

/* monom():
   Evaluates x^n.							     */
double monom(int n,double x)
{
 double w,t;
 w=1.0;t=x;
 while ( n>0 )
  {	if ( n&1 )	w=w*t;
	t=t*t,n=n/2;
  };
 return(w);
}

/* eval_gaussian():
   Evaluates exp(-1/2*(u*u+v*v)/(...))*u^(...)*v^(...), depending on the 
   parameters (sigma, basis bx and by) of the Gaussian kernel 'k'.	     */
double eval_gaussian(kernel *k,double u,double v)
{
 double	w;
 w=exp(-(u*u+v*v)/(2.0*k->sigma*k->sigma))*monom(k->bx,u)*monom(k->by,v); 
 return(w);
}

/* kernel_image_norm():
   Norms kernel image 'k->image' so sum of the kernel pixels will be 'sum'.  */
int kernel_image_norm(kernel *k,double sum)
{
 int	i,j;
 double	w;

 if ( k->image==NULL )	return(1);
	
 w=0.0;
 for ( i=0 ; i<=2*k->hsize ; i++ )
  {	for ( j=0 ; j<=2*k->hsize ; j++ )
 	 {	w=w+k->image[i][j];		}
  }
 w=sum/w;
 for ( i=0 ; i<=2*k->hsize ; i++ )
  {	for ( j=0 ; j<=2*k->hsize ; j++ )
	 {	k->image[i][j]=k->image[i][j]*w;		 }
  } 

 return(0);
}

/* kernel_image_calc_gaussian():
   Calculates 'k->image' for a given 'k->hsize' using the parametes of a 
   Gaussian kernel. Uses a grid with a size of 'subg' times 'subg' for
   the numerical integration on each pixel.				     */
int kernel_image_calc_gaussian(kernel *k)
{
 int	i,j,hsize,bx,by,subg;
 hsize=k->hsize;
 bx=k->bx,
 by=k->by;
 k->image=matrix_alloc(2*hsize+1);
 subg=10;
 for ( i=0 ; i<=2*hsize ; i++ )
  { for ( j=0 ; j<=2*hsize ; j++ )
     {	double	w,dx,dy;
	int	m,n;
	w=0.0;
	for ( m=0 ; m<subg ; m++ )
	 {  dy=(double)(1+2*m-subg)/(double)(2*subg);
	    for ( n=0 ; n<subg ; n++ )
	     {	dx=(double)(1+2*n-subg)/(double)(2*subg);
		w+=eval_gaussian(k,(double)(j-hsize)+dx,(double)(i-hsize)+dy);
	     }
	 }
	k->image[i][j]=w/(double)(subg*subg);
     }
  }
 if ( bx%2==0 && by%2==0 )
  {	kernel_image_norm(k,1.0);		}
 k->offset=0.0;
 return(0);
}
int kernel_image_calc_linear(kernel *k)
{
 int	i,j,mx,px,py;

 px=k->bx,py=k->by;
 mx=abs(px);
 if ( abs(py)>mx ) mx=abs(py);
 k->hsize=mx;
 k->image=matrix_alloc(2*mx+1);
 for ( i=0 ; i<=2*mx ; i++ )
  {	for ( j=0 ; j<=2*mx ; j++ )
	 {	k->image[i][j]=0.0;	}
  }
 k->image[mx+ 0][mx+ 0]=-0.5;
 k->image[mx+py][mx+px]=+0.5;
 k->offset=0.0;
 return(0);
}

int kernel_image_subtract(kernel *k,kernel *subk)
{
 int	i,j,hsk,hss,mh;
 
 mh=hsk=k->hsize;
 hss=subk->hsize;
 if ( hss<mh )	mh=hss;

 for ( i=-mh ; i<=mh ; i++ )
  {	for ( j=-mh ; j<=mh ; j++ )
	 {	k->image[i+hsk][j+hsk]-=subk->image[i+hss][j+hss];	}
  }

 return(0);
}

/*****************************************************************************/

double convolve_point(fitsimage *img,kernel *k,int x,int y)
{
 int	i,j,hsize;
 double	w;
 w=0.0;
 if ( k->type==KERNEL_UNKNOWN || k->type==KERNEL_GAUSSIAN )
  {	hsize=k->hsize;
	for ( i=-hsize ; i<=hsize ; i++ )
	 {	for ( j=-hsize ; j<=hsize ; j++ )
		 {	w+=img->data[y-i][x-j]*k->image[hsize+i][hsize+j]; }
	 }
	w+=k->offset;
  }
 else if ( k->type==KERNEL_BACKGROUND )
  {	w=1.0;					}
 else if ( k->type==KERNEL_IDENTITY )
  {	w=img->data[y][x];			}
 else if ( k->type==KERNEL_DDELTA )
  {	int	bx,by;
	bx=k->bx,by=k->by;
	w=0.5*(img->data[y-by][x-bx]-img->data[y][x]);
  }
 else	fprintf(stderr,"???\n");
 return(w);
}

int convolve_image(fitsimage *img,char **mask,kernel *k,fitsimage *out)
{
 int	i,j,sx,sy;

 if ( out->sx != img->sx || out->sy != img->sy )	return(1);
 sx=img->sx,sy=img->sy;

 for ( i=0 ; i<sy ; i++ )
  { for ( j=0 ; j<sx ; j++ )
     {	if ( mask[i][j] )	continue;
	out->data[i][j]=convolve_point(img,k,j,i);
     }
  }
		
 return(0);
}

/*****************************************************************************/

typedef struct _stamp_struct
 {	int	x,y;
	int	sx,sy;
	double	xc,yc;
 } stamp;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
					
static int fit_kernel_poly_coefficients_subblock(fitsimage *refimg,fitsimage *img,
	char **mask,double **weight,
	stamp *stamps,int nstamp,kernellist *klist,kernellist *xlist)
{
 int	b,i,j,k,l,nvar,nkernel,nvp,fkernel;
 int	sx,sy,mxorder;
 double	**amatrix,*bvector;
 double	**nsamatrix,*nsbvector,*nsconv;
 int	*koffset,*foffset;
 
 double	*monoms;
 kernel	*kernels,*wk;

 double	ox,oy,scale;

 sx=refimg->sx,sy=refimg->sy;
 if ( img->sx != sx || img->sy != sy )	return(1);

 ox=0.5*(double)(sx);
 oy=0.5*(double)(sy);
 scale=0.5*(double)sx;

 for ( i=0 ; i<klist->nkernel ; i++ )
  {	klist->kernels[i].target=0;		}
 if ( xlist != NULL )
  {	for ( i=0 ; xlist->kernels != NULL && i<xlist->nkernel ; i++ )
	 {	kernel	*x,*k;
		int	nk;
		x=&xlist->kernels[i];
		if ( x->type == KERNEL_BACKGROUND || x->type == KERNEL_IDENTITY )
			continue;
		nk=klist->nkernel;
		klist->kernels=(kernel *)realloc(klist->kernels,sizeof(kernel)*(nk+1));
		k=&klist->kernels[nk];
		memcpy(k,x,sizeof(kernel));
		k->target=1;
		k->flag=1;
		klist->nkernel++;
	 }
  }

 kernels=klist->kernels;
 nkernel=klist->nkernel;

 koffset=(int *)malloc(sizeof(int)*nkernel);
 foffset=(int *)malloc(sizeof(int)*nkernel);
 nvar=0;mxorder=0;fkernel=0;

 for ( i=0 ; i<nkernel ; i++ )	/* Let's count how many kernels do have to   */
  {	int	order;		/* be fitted. Foffset[k] refers to the kth   */
				/* kernel to be fitted, totally there are    */
				/* fkernel such kernels. The 'flag' indica-  */
				/* tes the status of a kernel. If it is      */
				/* negative, the kernel should be ignored.   */
				/* If it is zero, the kernel should be taken */
				/* into account during the fit with the appr-*/
				/* opriate order and known coefficients. If  */
				/* it is positive, the coefficients should   */
				/* be fitted up to the given order 'order'.  */
				
	koffset[i]=-1;		
	order=kernels[i].order;
	if ( order>mxorder )	mxorder=order;

	if ( kernels[i].flag<=0 )	continue;
	koffset[i]=nvar;
	nvar += (order+1)*(order+2)/2;
	foffset[fkernel]=i;
	fkernel++;
  }
 if ( fkernel==0 )		/* There's nothing to be fitted. */
  {	free(koffset);
	free(foffset);
	return(0);
  }
 nvp=(mxorder+1)*(mxorder+2)/2;

#if KERNEL_LOG != 0
 logmsg(kernel_verbose>=2," -> allocating desired memory "
 "[nkernel=%d,fkernel=%d,nvar=%d]...",nkernel,fkernel,nvar);
#endif
	amatrix=matrix_alloc(nvar);
	bvector=vector_alloc(nvar);
	monoms=vector_alloc(nvp);
	nsamatrix=matrix_alloc(fkernel);
	nsbvector=vector_alloc(fkernel);
	nsconv=vector_alloc(nkernel);
#if KERNEL_LOG != 0
 logmsg(kernel_verbose>=2,"done.\n");
 logmsg(kernel_verbose>=2," -> calculating least-squares matrix...");
#endif

 for ( i=0 ; i<nvar ; i++ )			/* Initialization of the LS  */
  {	for ( j=0 ; j<nvar ; j++ )		/* matrix and vector...	     */
	 {	amatrix[i][j]=0.0;	}
	bvector[i]=0.0;
  }

 for ( b=0 ; b<nstamp ; b++ )
  {	int	imin,imax,jmin,jmax,npix;
	double	xc,yc,w;
	
	jmin=stamps[b].x,
	imin=stamps[b].y,
	jmax=jmin+stamps[b].sx,
	imax=imin+stamps[b].sy;
	xc=stamps[b].xc,
	yc=stamps[b].yc;

	/* Calculate monomials for the spatial distibution. */
	eval_2d_monoms(xc,yc,mxorder,monoms,ox,oy,scale);

	for ( k=0 ; k<fkernel ; k++ )
 	 {	for ( l=0 ; l<fkernel ; l++ )
		 {	nsamatrix[k][l]=0.0;		}
		nsbvector[k]=0.0;
	 }

	for ( i=imin,npix=0 ; i<imax ; i++ )
	 {  for ( j=jmin ; j<jmax ; j++ )
	     {
		int	fk,fl;
		double	d;

		if ( mask[i][j] )	/* Ignore nasty pixels (bad pixels or*/
			continue;	/* pixels close to the boundaries    */
					/* where the convolution cannot be   */
					/* performed).			     */

		for ( k=0 ; k<nkernel ; k++ )	/* Build the LS-basis...     */
		 {  if ( kernels[k].flag>=0 )
		     {	if ( ! kernels[k].target )
			    nsconv[k]=convolve_point(refimg,&kernels[k],j,i);
			else
			    nsconv[k]=-convolve_point(img,&kernels[k],j,i);
		     }
		    else
			nsconv[k]=0.0;
		 }

		if ( weight != NULL ) 	w=weight[i][j];
		else			w=1.0;

		if ( w<=0.0 )	continue;

		d=img->data[i][j];	/* The kernels with flag==0 should   */
					/* be taken into account here...     */
		for ( k=0,wk=kernels ; k<nkernel ; k++,wk++ )
		 {	int	kn;
			if ( wk->flag != 0 )	continue;
			kn=(wk->order+1)*(wk->order+2)/2;
			for ( l=0 ; l<kn ; l++ )
				d-=nsconv[k]*wk->coeff[l]*monoms[l];
		 }

					/* Build the LS matrix, using the    */
					/* elements of the LS-basis nsconv[].*/

		for ( k=0 ; k<fkernel ; k++ )
		 {  fk=foffset[k];
		    for ( l=k ; l<fkernel ; l++ )
		     {	fl=foffset[l];
			nsamatrix[k][l]+=nsconv[fk]*nsconv[fl]*w;
		     }
		    nsbvector[k]+=d*nsconv[fk]*w;
		 }

		npix++;

	     }
	 }

	if ( npix <= 0 )	continue;
	
	for ( k=0 ; k<fkernel ; k++ )
	 {	for ( l=0 ; l<k ; l++ )
		 {	nsamatrix[k][l]=nsamatrix[l][k];		}
	 }


	for ( k=0 ; k<fkernel ; k++ )		/* Build the space-varying   */
	 {  int	korder,koff,knv,kk,fk;		/* LS-matrix, using constant */
	    fk=foffset[k];			/* LS-matrices in the given  */
	    korder=kernels[fk].order,		/* blocks. So, the size of   */
	    knv=(korder+1)*(korder+2)/2,	/* the blocks stamps[] should*/
            koff=koffset[fk];			/* be small enough to yield a*/
	    for ( l=0 ; l<fkernel ; l++ )	/* successful fit.	     */
	     {	int	lorder,loff,lnv,ll,fl;
		fl=foffset[l];
		lorder=kernels[fl].order,
	        lnv=(lorder+1)*(lorder+2)/2,
                loff=koffset[fl];
		for ( kk=0 ; kk<knv ; kk++ )
		 {  for ( ll=0 ; ll<lnv ; ll++ )
		     {	amatrix[koff+kk][loff+ll]+=nsamatrix[k][l]*monoms[kk]*monoms[ll];	}
		 }
	     }
	    for ( kk=0 ; kk<knv ; kk++ )	/* Build the spac-varying    */
						/* LS-vector, like above...  */
	     {	bvector[koff+kk]+=nsbvector[k]*monoms[kk];	}
	 }

  }
#if KERNEL_LOG != 0
 logmsg(kernel_verbose>=2,"done.\n");
 logmsg(kernel_verbose>=2," -> solving the linear least-squares equation...");
#endif
 solve_gauss(amatrix,bvector,nvar);
#if KERNEL_LOG != 0
 logmsg(kernel_verbose>=2,"done.\n");
 logmsg(kernel_verbose>=2," -> saving fitted parameters to the kernel description array...");
#endif

 klist->type=1;
 klist->ox=ox,klist->oy=oy;
 klist->scale=scale;
 for ( i=0 ; i<fkernel ; i++ )
  {	int	korder,knv;
	int	fi;

	fi=foffset[i];
	wk=&kernels[fi];

	korder=wk->order;
	knv=(korder+1)*(korder+2)/2;

	if ( wk->coeff != NULL )	
		free(wk->coeff);
	wk->coeff=(double *)malloc(sizeof(double)*knv);

	k=koffset[fi];
	for ( j=0 ; j<knv ; j++ )
	 {	wk->coeff[j]=bvector[k+j];		}
  }

#if KERNEL_LOG != 0
 logmsg(kernel_verbose>=2,"done.\n");
#endif

 matrix_free(nsamatrix);
 vector_free(nsbvector);
 vector_free(nsconv);

 vector_free(monoms);
 vector_free(bvector);
 matrix_free(amatrix);

 free(foffset);
 free(koffset);

 if ( xlist != NULL )
  {	xlist->nkernel=0;
	for ( i=0 ; i<nkernel ; )
	 {	wk=&kernels[i];
		if ( wk->target )
		 {	memcpy(&xlist->kernels[xlist->nkernel],wk,sizeof(kernel));
			xlist->nkernel++;
			if ( i<nkernel-1 )
				memmove(kernels+i,kernels+i+1,sizeof(kernel)*(nkernel-i-1));
			nkernel--;
		 }
		else
			i++;
	 }
	xlist->ox=ox,xlist->oy=oy;
	xlist->scale=scale;
	xlist->type=1;
  }

 klist->kernels=kernels;
 klist->nkernel=nkernel;

 return(0);
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int fit_kernel_poly_coefficients_block(fitsimage *refimg,fitsimage *img,
	char **mask,double **weight,int nbx,int nby,kernellist *klist,kernellist *xlist)
{
 int	bi,bj,b,ret,sx,sy;
 int	nstamp;
 stamp	*stamps;

 nstamp=nbx*nby;
 stamps=(stamp *)malloc(sizeof(stamp)*nstamp);

 sx=img->sx,
 sy=img->sy;
 
 for ( bi=0 ; bi<nby ; bi++ )
  { for ( bj=0 ; bj<nbx ; bj++ )
     {	int	imin,imax,jmin,jmax;
	double	xc,yc;

	imin=bi*sy/nby,imax=(bi+1)*sy/nby;
	jmin=bj*sx/nbx,jmax=(bj+1)*sx/nbx;
	xc=(double)(jmax+jmin-1)*0.5;
	yc=(double)(imax+imin-1)*0.5;

	b=bi*nbx+bj;

	stamps[b].x=jmin,
	stamps[b].y=imin;
	stamps[b].sx=jmax-jmin,
	stamps[b].sy=imax-imin;
	
	stamps[b].xc=xc;	
	stamps[b].yc=yc;
     }
  }

 ret=fit_kernel_poly_coefficients_subblock(refimg,img,mask,weight,stamps,nstamp,klist,xlist);

 free(stamps);

 return(ret);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int fit_kernel_poly_coefficients_stamp(fitsimage *refimg,fitsimage *img,
	char **mask,double **weight,rectan *rects,int nrect,kernellist *klist,kernellist *xlist)
{
 int	b,ret;
 int	nstamp;
 stamp	*stamps;

 nstamp=nrect;
 stamps=(stamp *)malloc(sizeof(stamp)*nstamp);
 
 for ( b=0 ; b<nrect ; b++ )
  {	int	imin,imax,jmin,jmax;
	double	xc,yc;

	imin=rects[b].y1,imax=rects[b].y2;
	jmin=rects[b].x1,jmax=rects[b].x2;
	xc=(double)(jmax+jmin-1)*0.5;
	yc=(double)(imax+imin-1)*0.5;

	stamps[b].x=jmin,
	stamps[b].y=imin;
	stamps[b].sx=jmax-jmin,
	stamps[b].sy=imax-imin;
	
	stamps[b].xc=xc;	
	stamps[b].yc=yc;
  }

 ret=fit_kernel_poly_coefficients_subblock(refimg,img,mask,weight,stamps,nstamp,klist,xlist);

 free(stamps);

 return(ret);
}

/*****************************************************************************/

int convolve_with_kernel_set_poly(fitsimage *refimg,char **mask,kernellist *klist,fitsimage *outimg)
{
 int	i,j,sx,sy;
 double	ox,oy,scale;

 sx=refimg->sx,sy=refimg->sy;
 if ( outimg->sx != sx || outimg->sy != sy )	return(1);

 ox=klist->ox,
 oy=klist->oy;
 scale=klist->scale;


 for ( i=0 ; i<sy ; i++ )
  { 
#if KERNEL_LOG != 0
    logmsg(kernel_verbose>=2,"\rPerforming the linear transformation... [%4d/%4d] ",i+1,sy);
#endif
    for ( j=0 ; j<sx ; j++ )
     {	double	c,p,w;
	int	n;
	kernel	*k;

	if ( mask[i][j] )
	 {	outimg->data[i][j]=0.0;
		continue;
	 }

	c=0.0;
	for ( n=0,k=klist->kernels ; n<klist->nkernel ; n++,k++ )
	 {	w=eval_2d_poly((double)j,(double)i,k->order,k->coeff,ox,oy,scale);
		p=convolve_point(refimg,k,j,i);
		c += w*p;
 	 }
	
	outimg->data[i][j]=c;
     }
  }
#if KERNEL_LOG != 0
 logmsg(kernel_verbose>=2,"done.\n");
#endif

 return(0);
}


int convolve_with_kernel_set(fitsimage *refimg,char **mask,kernellist *klist,fitsimage *outimg)
{
 int	r;
 if ( klist->type ) r=convolve_with_kernel_set_poly(refimg,mask,klist,outimg);
 else		    r=-1;
 return(r);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int convolve_to_subtracted_poly(fitsimage *ref,fitsimage *img,char **mask,
	kernellist *klist,kernellist *xlist,fitsimage *outimg)
{
 int	i,j,sx,sy;
 double	ox,oy,scale;

 sx=ref->sx,sy=ref->sy;
 if ( outimg->sx != sx || outimg->sy != sy )	return(1);
 if ( img->sx != sx || img->sy != sy )		return(1);

 ox=klist->ox,
 oy=klist->oy;
 scale=klist->scale;

 for ( i=0 ; i<sy ; i++ )
  { 
#if KERNEL_LOG != 0
    logmsg(kernel_verbose>=2,"\rPerforming the linear transformation... [%4d/%4d] ",i+1,sy);
#endif
    for ( j=0 ; j<sx ; j++ )
     {	double	c,p,w;
	int	n;
	kernel	*k;

	if ( mask[i][j] )
	 {	outimg->data[i][j]=0.0;
		continue;
	 }

	c=img->data[i][j];
	if ( xlist != NULL )
	 {	for ( n=0,k=xlist->kernels ; n<xlist->nkernel ; n++,k++ )
		 {	w=eval_2d_poly((double)j,(double)i,k->order,k->coeff,ox,oy,scale);
			p=convolve_point(img,k,j,i);
			c += w*p;
	 	 }
	 }

	for ( n=0,k=klist->kernels ; n<klist->nkernel ; n++,k++ )
	 {	w=eval_2d_poly((double)j,(double)i,k->order,k->coeff,ox,oy,scale);
		p=convolve_point(ref,k,j,i);
		c -= w*p;
 	 }
	
	outimg->data[i][j]=c;
     }
  }
#if KERNEL_LOG != 0
 logmsg(kernel_verbose>=2,"done.\n");
#endif

 return(0);
}


int convolve_to_subtracted(fitsimage *ref,fitsimage *img,char **mask,
	kernellist *klist,kernellist *xlist,fitsimage *outimg)
{
 int	r;
 if ( klist->type ) r=convolve_to_subtracted_poly(ref,img,mask,klist,xlist,outimg);
 else		    r=-1;
 return(r);
}

/*****************************************************************************/

