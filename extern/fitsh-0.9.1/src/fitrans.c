/*****************************************************************************/
/* fitrans.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line utility for transforming FITS images (including both 	     */
/*   geometrical and spatial transformations). Currently implemented methods */
/*   are listed below.							     */
/* Geometrical transformations:						     */
/*  - zooming (inclding spline subpixel interpolation, w/ flux conservation) */
/*  - magnifiying (like zooming but no spline subpixel interpolation)        */
/*  - shrinking (w/ exact flux conservation)				     */
/*  - generic geometric transformation (incl. flux conservation)	     */
/*	* [-m] interpolate + steplike: very fast, no exact flux conservation */
/*	* [-l] integrate   + steplike: slower, exact flux conservation (EFC) */
/*	* [-c] interpolate + spline  : fast, use splines but no EFC	     */
/*	* [-k] integrate   + spline  : very slow but results the best w/ EFC */
/*  - flipping/mirroring both in X and Y directions			     */
/*  - trimming or expanding images (can be combined with the above methods)  */
/* Spatial transformations (these do not affect the geometry /size, etc./):  */
/*  - noise/signal level estimation 					     */
/*  - large scale smoothing of images					     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2004, 2007, 2008; Pal, A. (apal@szofi.elte.hu)			     */
/*****************************************************************************/
#define	FI_TRANS_VERSION	"0.9z6"
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include <fits/fits.h>

#include "longhelp.h"
#include "fi.h"

#include "fitsmask.h"
#include "io/iof.h"
#include "io/scanarg.h"
#include "io/tokenize.h"
#include "io/format.h"
#include "math/fit/lmfit.h"
#include "math/poly.h"
#include "fbase.h"
#include "math/spline/biquad.h"
#include "math/spline/biquad-isc.h"
#include "math/spline/bicubic.h"
#include "math/spline/spline.h"

#include "statistics.h"
#include "transform.h"
#include "tensor.h"
#include "common.h"
#include "weight.h"
#include "history.h"

#ifdef  HAVE_NO_CC_EXTENSION 
#define __extension__ 
#endif 

/*****************************************************************************/

#define		TRANS_INTERPOLATE		0x00
#define		TRANS_INTEGRATE			0x01
#define		TRANS_STEPLIKE			0x00
#define		TRANS_SPLINE			0x02

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

#define		SMOOTH_NONE			0
#define		SMOOTH_SPLINE			1
#define		SMOOTH_POLYNOMIAL		2

#define		SMOOTH_PREFILTER_NONE		0
#define		SMOOTH_PREFILTER_MEAN		1
#define		SMOOTH_PREFILTER_MEDIAN		2

typedef struct
 {	int	type;			/* see SMOOTH_*			     */
	int	xorder,yorder;		/* polynomial and/or spoine order    */
	int	filter;			/* see SMOOTH_PREFILTER_*	     */
	int	fxhsize,fyhsize;	/* filter blocx x/y half size	     */
	double	frejratio;		/* filter rejection level ratio	(<1) */
	int	niter;			/* number of iterations		     */
	double	lower,upper;		/* rejection level in sigma	     */
	int	is_mean_unity;		/* scale result to have unity mean   */ 
	int	is_detrend;		/* detrend image instead of smoothing*/
 } smooth;

/*****************************************************************************/

int min4(int a,int b,int c,int d)
{
 if ( b<a ) a=b;
 if ( c<a ) a=c;
 if ( d<a ) a=d;
 return(a);
}
int max4(int a,int b,int c,int d)
{
 if ( b>a ) a=b;
 if ( c>a ) a=c;
 if ( d>a ) a=d;
 return(a);
}

/*****************************************************************************/

int fprint_fitrans_usage(FILE *fw)
{
 fprintf(fw,
"Usage:\tfitrans [-h|--help] [-C|--comment] [-V|--verbose]\n"
"\t[[-i] <in>[<F>] [--frame <F>] [-M|--input-mask <mask>]] [-o <out>]\n"
"General (analytical) transformation:\n"
"\t[-T <transform-file>|-t <transform> [-m|-l|-c|-k] [--reverse]] \n"
"Layer extraction or slicing data cubeimages into disctint layers:\n"
"\t[-y|--layer <layer>]\n"
"\t[-x|--explode <basename> [-y|--first-layer <n>]]\n"
"Zooming/shrinking, shifting:\n"
"\t[-z|--zoom <factor>] [-g|--magnify <factor>]\n"
"\t[-r|--shrink <factor> [-d|--median] [--optimistic-masking]]\n"
"\t[-e|--shift <dx>,<dy>]\n");

 fprintf(fw,
"Noise estimation:\n"
"\t[-n|--noise]\n");

 fprintf(fw,
"Large scale smoothing:\n"
"\t[-a|--smooth {spline|polynomial},[xy]order=<order>,[unity],[detrend]\n"
"\t\t[mean|median],[xy]hsize=<halfsize>,[rejection=<ratio>],\n"
"\t\t[iterations=<iterations>],[lower|upper|sigma=<level>]]\n");

 fprintf(fw,
"Other optional parameters (override the default size/depth/offset values):\n"
"\t[-s|--size [-]<sx>,[-]<sy> [--flip-[x][y]]] [-f|--offset <ox>,<oy>]\n"
"\t[-b|--bitpix <bitpix>]\n"
"\t[-D|--data bitpix=<bitpix>,bscale=<scale>,bzero=<zero>|<C-type>]\n");
 
 return(0);
}

longhelp_entry fitrans_long_help[]=
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
 { "-o, --output <image file>",
	"Name of the output FITS image file." },
 { "-b, --bitpix <bitpix>",
	"Standard FITS output bitpix value." },
 { "-D, --data <spec>",
	"Output pixel data format specification." },

 { "Spatial transformations:", NULL },
 { "-T, --input-transformation <transformation file>",
	"Name of the file which contains the transformation description."
	"Such a file can be created e.g. by the programs `grtrans` or `grmatch`. "
	"This file contains basically the same set of <keyword> = <value> "
	"pairs as it is used after the -t|--transformation option (see there)." },
 { "-t, --transformation <transformation>",
	"Comma-separated list of parameters for the spatial transformation, "
	"see section ``Parameters for spatial transformations'' below." },
 { "-e, --shift <dx>,<dy>",
	"Imply a transformation that shifts the image by <dx>,<dy>."	},

 { "--reverse, --inverse",
	"Apply the inverse transformation to the image rather than the "
	"original one. " },
 { "-m",
	"Simple linear interpolation between pixels, with no exact flux "
	"conservation (just a multiplication by the Jacobian of "
	"the transformation)." },
 { "-l",
	"Linear interpolation between the pixels involving exact flux "
	"conservation by integrating on the image surface." },
 { "-c",
	"Bicubic spline interpolation between pixels, with no exact flux "
	"conservation (just a multiplication by the Jacobian of "
	"the transformation)." },
 { "-k",
	"Interpolation by integrationg the flux on a biquadratic interpolation "
	"surface, yielding exact flux conservation. " },

 { "-s, --size <sx>,<sy>",
	"The size of the output image if it should differ from the original "
	"image size." },
 { "-f, --offset <x>,<y>",
	"Zero-point coordinate of the output image in the input image. " },

 { "Parameters for spatial transformations:", NULL },
 { "type=<type>",
	"Type of the transformation. In the actual implementation, the "
	"only supported type for a transformation is ``polynomial''. " },
 { "order=<order>",
	"Polynomial order for the transformation." },
 { "dxfit=<coefficients>",
	"Comma-separated list of the polynomial coefficients for the "
	"X coordinate. The number of coefficients must be 1, 3, 6, ... for " 
	"the orders 0, 1, 2, ... respectively." },
 { "dyfit=<coefficients>",
	"Comma-separated list of the polynomial coefficients for the "
	"Y coordinate." },

 { "Other simple spatial geometric transformations:", NULL },
 { "-z, --zoom <factor>",
	"Zoom the image by the given (integer) factor, involving a biquadratic "
	"subpixel-level interpolation and therefore exact flux conservation. " },
 { "-r, --shrink <factor>",
	"Shrink the image by the given (integer) factor. " },
 { "-d, --median",
	"Use a median-based averaging during the shrinking operation. " },
 { "--optimistic-masking",
	"Imply some optimism during the shrinking operation: masked pixels "
	"are ignored during the averaging process and the final mask will "
	"be computed in a complement manner. " },
 { "-g, --magnify <factor>",
	"Same as zooming the image but there is no subpixel-level interpolation." },

 { "Large-scale image smoothing:", NULL },
 { "-a, --smooth <parameters>",
	"Perform a smoothing on the image. The parameters of the smoothing "
	"are the following:" },
 { "spline",
	"Do a spline interpolation smoothing" },
 { "polynomial",
	"Do a polynomial interpolation smoothing" },
 { "[xy]order=<order>",
	"Spatial order of the smoothing function. The order in the X and Y "
	"coordinates can be set independently, by setting ``xorder=...'' or "
	"``yorder=...''." },
 { "unity",
	"Scale the resulting smoothed image to have a mean of 1." },
 { "detrend",
	"The resulting image will be the original image divided by the "
	"best fit smoothed surface." },
 { "[xy]hsize=<halfsize>",
	"Do a box filtering  with the given halfsize." },
 { "mean",
	"Use the mean value of the pixels for the box filtering." },
 { "median",
	"Use the median value of the pixels for the box filtering." },
 { "iterations=<iterations>",
	"Number of iterations to reject outlier pixels from the box." },
 { "lower, upper, sigma=<sigma>",
	"Lower, upper or symmetric rejection level in the units of "
	"standard deviation." },
 
 { "Noise estimation:", NULL },
 { "-n, --noise",
	"Derive an image which reflects the ``noise level'' of the image. " },

 { "Slicing or exploding data cube images:", NULL },
 { "-y, --layer <layer>",
	"Layer (z-axis index) of the desired image slice." },
 { "-x, --explode <basename>",
	"Explode the input image into individual planar (two dimensional) "
	"FITS image. The basename must contain at least one printf-like "
	"tag of %d, %i, %o, %x or %X that is replaced by the appropriate "
	"layer number index. " },
 { "-y, --first-layer <n>",
	"Use the specified value for the first layer index. The subsequent "
	"layer indices are incremented normally. By default, the index of "
	"the first data cube layer is 0. " },
 
 { NULL, NULL }
};


int fprint_fitrans_long_help(FILE *fw)
{
 fprintf(fw,
"Usage:\tfitrans [transformation and options] <input> [-o|--output <output>]\n"
"The main purpose of this program is to perform specific or generic geometric\n"
"transformations on the input image.\n");
 fprintf(fw,"\n");

 longhelp_fprint(fw,fitrans_long_help,0,-1);

 fprintf(fw,"\n");
 fprintf(fw,"Report bugs to <%s>\n",FI_MAINT_EMAIL);

 return(0);

}

/*****************************************************************************/

/* prefilter_image():

  Do a moving block filtering on the data matrix 'data' (with the respecive
  size of 'sx' by 'sy') optionally taking in to account the mask 'mask' marking
  the bad and/or unexpected points. The resulted block filtered image is 
  then stored in 'fltd' (with the same size) and the appropriate mask is 
  stored in 'fmsk' (if it is not NULL). The following smooth parameters
  are used from the object 'sp':
	sp->filter	(see SMOOTH_PREFILTER_* definitions)
	sp->f[xy]hsize	(block halfsizes, both 0 means do nothing)
	sp->frejratio	(rejection ratio if SMOOTH_PREFILTER_MEAN is used)
*/

typedef struct
 {	int	x,y;
	double	val;
 } ppixel;

static int ppixel_compare(const void *v1,const void *v2)
{
 ppixel	*p1=(ppixel *)v1;
 ppixel	*p2=(ppixel *)v2;
 if ( p1->val < p2->val )
	return(-1);
 else
	return(1);
}
int ppixel_sort(ppixel *pps,int npp)
{
 if ( pps != NULL && npp>0 )
	qsort(pps,npp,sizeof(ppixel),ppixel_compare);

 return(0);
}

int prefilter_image(double **data,int sx,int sy,double **fltd,
	smooth *sp,char **mask,char **fmsk)
{
 int	hx,hy;
 int	i,j,k,l,k0,k1;
 ppixel	*pps,*ppt,*ppw,*ppn,**apr;
 double	sumpp;
 int	npp,tpp,nnp;
 int	*ran;

 if ( data==NULL || fltd==NULL )
	return(-1);

 if ( ! sp->filter || (sp->fxhsize <= 0 && sp->fyhsize <= 0) )
  {	for ( i=0 ; i<sy ;  i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	fltd[i][j]=data[i][j];
			if ( fmsk != NULL && mask != NULL )
				fmsk[i][j]=mask[i][j];
			else if ( fmsk != NULL )
				fmsk[i][j]=0;
		 }
	 }
	return(0);
  }

 hx=sp->fxhsize;
 hy=sp->fyhsize;

 pps=(ppixel *)malloc(sizeof(ppixel)*(2*hx+1)*(2*hy+1));
 ppt=(ppixel *)malloc(sizeof(ppixel)*(2*hx+1)*(2*hy+1));
 ppn=(ppixel *)malloc(sizeof(ppixel)*(2*hy+1));
 apr=(ppixel **)tensor_alloc_2d(ppixel,2*hy+1,sx);
 ran=(int *)malloc(sizeof(int)*sx);

 for ( j=0 ; j<sx ; j++ )
  {	l=0;
	for ( k=0 ; k<=hy && k<sy ; k++ )
	 {	if ( mask != NULL && mask[k][j] )
			continue;
		apr[j][l].x=j;
		apr[j][l].y=k;
		apr[j][l].val=data[k][j];
		l++;
	 }
	ppixel_sort(apr[j],l);
	ran[j]=l;
  }

 for ( i=0 ; i<sy ; i++ )
  {	npp=0;
	sumpp=0.0;

	k0=-hy;
	if ( i+k0<0 )	k0=-i;
	k1=+hy;
	if ( i+k1>=sy )	k1=sy-1-i;

	for ( k=k0 ; k<=k1 ; k++ )
	 {	for ( l=0 ; l<=hx && l<sx ; l++ )
		 {	if ( mask != NULL && mask[i+k][l] )
				continue;
			pps[npp].x=l;
			pps[npp].y=i+k;
			sumpp+=(pps[npp].val=data[i+k][l]);
			npp++;
		 }
	 }
	ppixel_sort(pps,npp);

	for ( j=0 ; j<sx ; j++ )
	 {	
		if ( npp>0 )
		 {	if ( sp->filter==SMOOTH_PREFILTER_MEDIAN )
			 {	fltd[i][j]=0.5*(pps[(npp-1)/2].val+pps[(npp)/2].val);	}
			else if ( sp->filter==SMOOTH_PREFILTER_MEAN && sp->frejratio<=0.0 )
			 {	fltd[i][j]=sumpp/npp;					}
			else if ( sp->filter==SMOOTH_PREFILTER_MEAN )
			 {	int	k,l0,l1;
				l0=(int)((double)npp*(sp->frejratio/2.0));
				l1=npp-l0;
				fltd[i][j]=0.0;
				for ( k=l0 ; k<l1 ; k++ )
				 {	fltd[i][j]+=pps[k].val;		}
				if ( l1>l0 )
					fltd[i][j]/=(double)(l1-l0);
				else
				 {	fltd[i][j]=0.0;
					fmsk[i][j]=MASK_FAULT;
				 }
			 }
			else /* unimplemented mode: */
			 {	fltd[i][j]=0.0;
				fmsk[i][j]=MASK_FAULT;
			 }
		 }
		else
		 {	fltd[i][j]=0.0;
			fmsk[i][j]=MASK_FAULT;	/* TBD: combination of masks */
		 }

		if ( ! (j<sx-1) )
			continue;

		/*
		nnp=0;
		if ( j+hx+1<sx )
		 {	for ( k=k0 ; k<=k1 ; k++ )
			 {	if ( mask != NULL && mask[i+k][j+hx+1] )
					continue;
				ppn[nnp].x=j+hx+1;
				ppn[nnp].y=i+k;
				ppn[nnp].val=data[i+k][j+hx+1];
				nnp++;
			 }
			ppixel_sort(ppn,nnp);
		 }
		*/

		if ( j+hx+1<sx )
		 {	ppw=apr[j+hx+1];
			nnp=ran[j+hx+1];
		 }
		else
		 {	ppw=NULL;
			nnp=0;
		 }

		tpp=0;
		sumpp=0.0;
		for ( k=0,l=0 ; k<npp ; k++ )
		 {	if ( pps[k].x > j-hx )
			 {	while ( l<nnp && ppw[l].val<pps[k].val )
				 {	ppt[tpp]=ppw[l];
					sumpp+=ppw[l].val;
					tpp++;
					l++;
				 }
				ppt[tpp]=pps[k];
				sumpp+=pps[k].val;
				tpp++;
			 }
		 }
		while ( l<nnp )
		 {	ppt[tpp]=ppw[l];
			sumpp+=ppw[l].val;
			tpp++;
			l++;
		 }

		ppw=pps,pps=ppt,ppt=ppw;
		npp=tpp;
	 }


	if ( ! (i<sy-1) )
		continue;

	for ( j=0 ; j<sx ; j++ )
	 {	l=0;

		if ( i+hy+1<sy && ( mask==NULL || (! mask[i+hy+1][j]) ) )
		 {	ppn[0].x=j;
			ppn[0].y=i+hy+1;
			ppn[0].val=data[i+hy+1][j];
			nnp=1;
		 }
		else
			nnp=0;

		tpp=0;

		for ( k=0,l=0 ; k<ran[j] ; k++ )
		 {	if ( apr[j][k].y > i-hy )
			 {	while ( l<nnp && ppn[l].val<apr[j][k].val )
				 {	ppt[tpp]=ppn[l];
					tpp++;
					l++;
				 }
				ppt[tpp]=apr[j][k];
				tpp++;
			 }
		 }
		while ( l<nnp )
		 {	ppt[tpp]=ppn[l];
			tpp++;
			l++;
		 }

		memcpy(apr[j],ppt,sizeof(ppixel)*tpp);
		ran[j]=tpp;
	 }
  }

 free(ran);
 tensor_free(apr);
 free(ppn);	
 free(ppt);
 free(pps);

 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* smooth_image():
  This call fits a smooth function to the data matrix 'data' (with the
  respiective size of 'sx' by 'sy') taking into account the optional mask 
  'mask' for marking bad or unexpected points. The resulted fitted image
  is stored in 'fltd' (with the same size) and the appropriate
  mask is stored in 'fmsk' (if it is not NULL). The following smooth parameters
  are used from the object 'sp':
	sp->type	(see SMOOTH_* definitions)
	sp->[xy]order	(x and y spatial variation orders)
	sp->niter	(sigma rejection iterations)
	sp->lower,upper	(rejection level for outlyers)
*/

int smooth_image(double **data,int sx,int sy,double **fltd,
	smooth *sp,char **mask,char **omsk)
{
 int	xorder,yorder;
 double	**fxbase,**fybase;
 double	*fvars,*bvector,**amatrix,w,f;
 int	nvar,nxvar,nyvar;
 char	**tmsk;
 int	i,j,iiter,k,l;
 
 if ( data==NULL || fltd==NULL )
	return(-1);

 xorder=sp->xorder;
 yorder=sp->yorder;

 fxbase=(double **)tensor_alloc_2d(double,sx,xorder+1);
 fybase=(double **)tensor_alloc_2d(double,sy,yorder+1);

 if ( sp->type==SMOOTH_SPLINE ) 
  {	fbase_spline(fxbase,xorder,sx);
	fbase_spline(fybase,yorder,sy);
  }
 else if ( sp->type==SMOOTH_POLYNOMIAL )
  {	fbase_polynomial(fxbase,xorder,sx);
	fbase_polynomial(fybase,yorder,sy);
  }
 else
  {	xorder=0;
	yorder=0;
	fbase_polynomial(fxbase,0,sx);
	fbase_polynomial(fybase,0,sy);
  }

 nxvar=xorder+1;
 nyvar=yorder+1;
 nvar=nxvar*nyvar;

 tmsk=(char **)tensor_alloc_2d(char,sx,sy);
 for ( i=0 ; i<sy ; i++ )
  {	if ( mask != NULL )
		memcpy(tmsk[i],mask[i],sx);
	else
		memset(tmsk[i],0,sx);
  }

 fvars  =vector_alloc(nvar);
 bvector=vector_alloc(nvar);
 amatrix=matrix_alloc(nvar);

 for ( iiter=0 ; iiter<=(sp->niter>0?sp->niter:0) ; iiter++ )
  {	for ( k=0 ; k<nvar ; k++ )
	 {	for ( l=0 ; l<nvar ; l++ )
		 {	amatrix[k][l]=0.0;	 }
		bvector[k]=0.0;
	 }
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	if ( tmsk[i][j] )
				continue;
			for ( k=0 ; k<nyvar ; k++ )
			 {	for ( l=0 ; l<nxvar ; l++ )
				 {	fvars[k*nxvar+l]=fybase[k][i]*fxbase[l][j];	}
			 }
			w=1.0;
			f=data[i][j];
			for ( k=0 ; k<nvar ; k++ )
			 {	for ( l=0 ; l<nvar ; l++ )
				 {	amatrix[k][l] += w*fvars[k]*fvars[l];	}
				bvector[k] += w*f*fvars[k];
			 }
		 }
	 }

	solve_gauss(amatrix,bvector,nvar);

	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	for ( k=0 ; k<nyvar ; k++ )
			 {	for ( l=0 ; l<nxvar ; l++ )
				 {	fvars[k*nxvar+l]=fybase[k][i]*fxbase[l][j];	}
			 }
			f=0.0;
			for ( k=0 ; k<nvar ; k++ )
				f+=bvector[k]*fvars[k];
			fltd[i][j]=f;
		 }
	 }

	if ( iiter < sp->niter )
	 {	double	s2,sig;
		s2=0.0;
		for ( i=0 ; i<sy ; i++ )
		 {	for ( j=0 ; j<sx ; j++ )
			 {	f=data[i][j]-fltd[i][j];
				s2+=f*f;
			 }
		 }
		s2/=(double)(sx*sy);
		if ( s2>0.0 )	sig=sqrt(s2);
		else		break;

		for ( i=0 ; i<sy ; i++ )
		 {	for ( j=0 ; j<sx ; j++ )
			 {	f=data[i][j]-fltd[i][j];
				if ( f < -sig*sp->lower || +sig*sp->upper < f )
					tmsk[i][j]=1;
			 }
		 }
	 }

  }

 for ( i=0 ; i<sy ; i++ )
  {	memset(omsk[i],0,sx);		  }

 matrix_free(amatrix);
 vector_free(bvector);
 vector_free(fvars);

 tensor_free(tmsk);

 free(fybase);
 free(fxbase);

 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int combined_smooth_image(fitsimage *img,char **mask,fitsimage *out,char **omsk,smooth *sp)
{
 int	sx,sy;
 double	**data;
 char	**fmsk;

 if ( img==NULL  || out==NULL  )	return(1);
 if ( mask==NULL || omsk==NULL )	return(1);
 sx=img->sx,
 sy=img->sy;
 if ( sx<=0 || sy<=0 )			return(1);
 if ( out->sx != sx || out->sy != sy )	return(1);

 if ( sp->filter )
  {	data=(double **)tensor_alloc_2d(double,sx,sy);
	fmsk=(char **)tensor_alloc_2d(char,sx,sy);
	prefilter_image(img->data,sx,sy,data,sp,mask,fmsk);
  }
 else
  {	data=img->data;
	fmsk=mask;
  }

 if ( sp->type )
	smooth_image(data,sx,sy,out->data,sp,fmsk,omsk);
 else
  {	int	i,j;
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	out->data[i][j]=data[i][j];
			omsk[i][j]=fmsk[i][j];
		 }
	 }
  }

 if ( sp->is_mean_unity )
  {	double	s,n;
	int	i,j;

	s=n=0.0;
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	if ( omsk[i][j] )
				continue;
			n+=1.0;
			s+=out->data[i][j];
		 }
	 }

	if ( s>0.0 )
	 {	s=n/s;
		for ( i=0 ; i<sy ; i++ )
		 {	for ( j=0 ; j<sx ; j++ )
			 {	out->data[i][j] *= s;		}
		 }
	 }
  }
 else if ( sp->is_detrend )
  {	double	s,n,d;
	int	i,j;

	s=n=0.0;
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	if ( omsk[i][j] )
				continue;
			n+=1.0;
			s+=out->data[i][j];
		 }
	 }
	s=s/n;
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	d=out->data[i][j];
			if ( d>0.0 )
				out->data[i][j]=img->data[i][j]*s/d;
			else
				out->data[i][j]=0.0;
		 }
	 }
  }

 if ( sp->filter )
  {	tensor_free(fmsk);
	tensor_free(data);
  }

 return(0);
}

/*****************************************************************************/

int interpolate_image_linear(fitsimage *img,char **mask,double x,double y,double *ret)
{
 int	ix,iy,flag;
 double	u,v;

 ix=(int)floor(x),iy=(int)floor(y);
 u=x-(double)ix,v=y-(double)iy;
 if ( x<0 || y<0 || ix>=img->sx-2 || iy>=img->sy-2 )	return(MASK_OUTER);
 if ( mask != NULL )
	flag=mask[iy][ix]|mask[iy+1][ix+1]|mask[iy][ix+1]|mask[iy+1][ix];	
 else
	flag=0;

 *ret=  img->data[iy+0][ix+0]*(1-u)*(1-v)+
	img->data[iy+0][ix+1]*( +u)*(1-v)+
	img->data[iy+1][ix+0]*(1-u)*( +v)+
	img->data[iy+1][ix+1]*( +u)*( +v);

 return(flag);
}

int interpolate_image_bicubic(double **c,int sx,int sy,char **mask,double fx,double fy,double *ret)
{
 int	ix,iy,flag;

 ix=(int)floor(fx),iy=(int)floor(fy);
 if ( fx<0 || fy<0 || ix>=sx-2 || iy>=sy-2 )	return(MASK_OUTER);
 if ( mask != NULL )
	flag=mask[iy][ix]|mask[iy+1][ix+1]|mask[iy][ix+1]|mask[iy+1][ix];	
 else	
	flag=0;

 *ret=bicubic_inter(c,fx,fy);

 return(flag);
}

int apply_transformation_polynomial_interpolate
	(fitsimage *img,char **mask,fitsimage *outimg,char **outmask,
	int ofx,int ofy,transformation *tf,int is_spline,int is_invert)
{
 int	i,j,sx,sy,nsx,nsy,flag,order;
 double	x,y,nx,ny,w,bshx,bshy,ox,oy,scale;
 double	*jxx,*jxy,*jyx,*jyy,cjxx,cjxy,cjyx,cjyy,cj;
 double	**splinecoeff;

 sx =img->sx   ,sy =img->sy;
 nsx=outimg->sx,nsy=outimg->sy;
 ox=tf->ox,oy=tf->oy,scale=tf->scale;
 order=tf->order;

 if ( is_spline )
  {	splinecoeff=(double **)tensor_alloc_2d(double,2*sx,2*sy);
	if ( splinecoeff==NULL )	return(-1);
	logmsg(is_verbose,"Calculating spline coefficients... ");
	bicubic_coeff(img->data,sx,sy,splinecoeff,mask);
	logmsg(is_verbose,"done.\n");
  }
 else	
	splinecoeff=NULL;
 
 transformation_get_jacobi(tf,&jxx,&jxy,&jyx,&jyy);

 logmsg(is_verbose,"Transforming... "); 
 bshx=tf->bshx,
 bshy=tf->bshy;

 for ( i=0 ; i<nsy ; i++ )
  { for ( j=0 ; j<nsx ; j++ )
     {	x=(double)(j-ofx)+bshx,
	y=(double)(i-ofy)+bshy;
	
	if ( ! is_invert )
		transformation_eval_normal_2d(x,y,tf,&nx,&ny);
	else
		transformation_eval_invert_2d(x,y,tf,&nx,&ny,jxx,jxy,jyx,jyy);

	nx-=bshx,ny-=bshy;
	cjxx=eval_2d_poly(x,y,order-1,jxx,ox,oy,scale);
	cjxy=eval_2d_poly(x,y,order-1,jxy,ox,oy,scale);
	cjyx=eval_2d_poly(x,y,order-1,jyx,ox,oy,scale);
	cjyy=eval_2d_poly(x,y,order-1,jyy,ox,oy,scale);

	if ( ! is_invert )	cj=cjxx*cjyy-cjxy*cjyx;
	else			cj=1.0/(cjxx*cjyy-cjxy*cjyx);

	if ( is_spline )
		flag=interpolate_image_bicubic(splinecoeff,sx,sy,mask,nx,ny,&w);
  	else
		flag=interpolate_image_linear(img,mask,nx,ny,&w);

	if ( ! flag )	outmask[i][j]=0,outimg->data[i][j]=w*cj;
	else		outmask[i][j]=flag,outimg->data[i][j]=0.0;
     }
    logmsg(is_verbose>=2,"\rTransforming... (%d/%d)... ",i+1,nsy);
  }
 logmsg(is_verbose,"done.\n");

 if ( splinecoeff != NULL )
	tensor_free(splinecoeff);
 if ( jyy != NULL )	free(jyy);
 if ( jyx != NULL )	free(jyx);
 if ( jxy != NULL )	free(jxy);
 if ( jxx != NULL )	free(jxx);

 return(0);
}

int apply_transformation_polynomial_integrate
	(fitsimage *img,char **mask,fitsimage *outimg,char **outmask,int ofx,int ofy,
	transformation *tf,int is_spline,int is_invert)
{
 int	i,j,sx,sy,nsx,nsy,flag,order;
 double	x,y,nx,ny,dx1,dy1,dx2,dy2,dx3,dy3,dx4,dy4,w,bshx,bshy,ox,oy,scale;
 double	**bqc,**dxyline,*wd;
 double	*jxx,*jxy,*jyx,*jyy;
 int	idx1,idy1,idx2,idy2,idx3,idy3,idx4,idy4,idxl,idxh,idyl,idyh,ix,iy;

 sx =img->sx   ,sy =img->sy;
 nsx=outimg->sx,nsy=outimg->sy;
 ox=tf->ox,oy=tf->oy,scale=tf->scale;
 order=tf->order;

 if ( is_spline )
  {	bqc=(double **)tensor_alloc_2d(double,2*sx+1,2*sy+1);
	if ( bqc==NULL )	return(-1);
	logmsg(is_verbose,"Calculating biquadratic spline coefficients... ");
	biquad_coeff(img->data,sx,sy,bqc,mask);
	logmsg(is_verbose,"done.\n");
  }
 else
	bqc=NULL;
 
 if ( is_invert )
	transformation_get_jacobi(tf,&jxx,&jxy,&jyx,&jyy);
 else
	jxx=jxy=jyx=jyy=NULL;

 logmsg(is_verbose,"Transforming... "); 
 bshx=tf->bshx,
 bshy=tf->bshy;

 dxyline=(double **)tensor_alloc_2d(double,nsx+1,4);

 y=(double)(0-ofy);
 for ( j=0 ; j<=nsx ; j++ )
  {	x=(double)(j-ofx);
	if ( ! is_invert )
		transformation_eval_normal_2d(x+bshx,y+bshy,tf,&nx,&ny);
	else
		transformation_eval_invert_2d(x+bshx,y+bshy,tf,&nx,&ny,jxx,jxy,jyx,jyy);
	dxyline[0][j]=nx-bshx,
	dxyline[1][j]=ny-bshy;
  }
 for ( i=0 ; i<nsy ; i++ )
  { 
    y=(double)(i+1-ofy);
    for ( j=0 ; j<=nsx ; j++ )
     {	x=(double)(j-ofx);
	if ( ! is_invert )
		transformation_eval_normal_2d(x+bshx,y+bshy,tf,&nx,&ny);
	else
		transformation_eval_invert_2d(x+bshx,y+bshy,tf,&nx,&ny,jxx,jxy,jyx,jyy);
	dxyline[2][j]=nx-bshx,
	dxyline[3][j]=ny-bshy;
     }

    for ( j=0 ; j<nsx ; j++ )
     {	dx1=dxyline[0][j],dx2=dxyline[0][j+1],
	dy1=dxyline[1][j],dy2=dxyline[1][j+1];
	dx3=dxyline[2][j],dx4=dxyline[2][j+1],
	dy3=dxyline[3][j],dy4=dxyline[3][j+1];

	idx1=(int)dx1,idx2=(int)dx2,idx3=(int)dx3,idx4=(int)dx4;
	idy1=(int)dy1,idy2=(int)dy2,idy3=(int)dy3,idy4=(int)dy4;

	flag=0;
   	     if ( dx1<0.0 || dy1<0.0 || dx2<0.0 || dy2<0.0 ) flag|=MASK_OUTER;
	else if ( dx3<0.0 || dy3<0.0 || dx4<0.0 || dy4<0.0 ) flag|=MASK_OUTER;	
	else if ( dx1>=sx || dy1>=sy || dx2>=sx || dy2>=sy ) flag|=MASK_OUTER;
	else if ( dx3>=sx || dy3>=sy || dx4>=sx || dy4>=sy ) flag|=MASK_OUTER;
	else if ( mask != NULL )
	 {	idxl=min4(idx1,idx2,idx3,idx4);
		if ( idxl<0   )	flag|=MASK_OUTER,idxl=0;
		idxh=max4(idx1,idx2,idx3,idx4);
		if ( idxh>=sx )	flag|=MASK_OUTER,idxh=sx-1;
	 	idyl=min4(idy1,idy2,idy3,idy4);
		if ( idyl<0   )	flag|=MASK_OUTER,idyl=0;
		idyh=max4(idy1,idy2,idy3,idy4);
		if ( idyh>=sy )	flag|=MASK_OUTER,idyh=sy-1;
		for ( iy=idyl ; iy<=idyh ; iy++ )
		 {	for ( ix=idxl ; ix<=idxh ; ix++ )
			 {	flag |= mask[iy][ix];		}
		 }
	 }
	if ( ! flag )
	 {	outmask[i][j]=0;
		if ( is_spline && bqc != NULL )
			w=biquad_isc_int_triangle(bqc,1,dx1,dy1,dx2,dy2,dx3,dy3,sx,sy)+
			  biquad_isc_int_triangle(bqc,1,dx2,dy2,dx4,dy4,dx3,dy3,sx,sy);
		else
			w=biquad_isc_int_triangle(img->data,0,dx1,dy1,dx2,dy2,dx3,dy3,sx,sy)+
			  biquad_isc_int_triangle(img->data,0,dx2,dy2,dx4,dy4,dx3,dy3,sx,sy);
		outimg->data[i][j]=w;
	 }
	else	
	 {	outmask[i][j]=flag;
		outimg->data[i][j]=0.0;
	 }
     }
    wd=dxyline[0],dxyline[0]=dxyline[2],dxyline[2]=wd;
    wd=dxyline[1],dxyline[1]=dxyline[3],dxyline[3]=wd;
    logmsg(is_verbose>=2,"\rTransforming... (%d/%d)... ",i+1,nsy);
  }
 logmsg(is_verbose,"done.\n");

 if ( bqc != NULL )	tensor_free(bqc);
 if ( dxyline != NULL )	tensor_free(dxyline);

 if ( jyy != NULL )	free(jyy);
 if ( jyx != NULL )	free(jyx);
 if ( jxy != NULL )	free(jxy);
 if ( jxx != NULL )	free(jxx);

 return(0);
}

int apply_transformation
	(fitsimage *img,char **mask,fitsimage *outimg,char **outmask,
	int ofx,int ofy,transformation *tf,int int_method,int is_invert)
{
 int	r;

 if ( tf->type == 1 )
  {	if ( int_method & TRANS_INTEGRATE )
		r=apply_transformation_polynomial_integrate(img,mask,outimg,outmask,ofx,ofy,tf,(int_method&TRANS_SPLINE),is_invert);
	else
		r=apply_transformation_polynomial_interpolate(img,mask,outimg,outmask,ofx,ofy,tf,(int_method&TRANS_SPLINE),is_invert);
  }
 else
	r=-1;

 return(r);

}

/*****************************************************************************/

int zoom_image
	(fitsimage *img,char **mask,fitsimage *outimg,char **outmask,
	int ofx,int ofy,int scale)
{
 int	sx,sy,nsx,nsy;
 int	i,j,k,l,m,ii,jj;
 double	**wret,**bqc;

 if ( img==NULL || outimg==NULL )	return(1);
 if ( mask==NULL || outmask==NULL )	return(1);
 sx=img->sx,
 sy=img->sy;
 nsx=outimg->sx;
 nsy=outimg->sy;
 if ( sx==0 || sy==0 || nsx==0 || nsy==0 )	return(0);

 wret=(double **)tensor_alloc_2d(double,scale,scale);
 bqc =(double **)tensor_alloc_2d(double,2*sx+1,2*sy+1);
 logmsg(is_verbose,"Calculating biquadratic coefficients... ");
 biquad_coeff(img->data,sx,sy,bqc,mask);
 logmsg(is_verbose,"done.\n");

 logmsg(is_verbose==1,"Zooming... ");
 for ( i=0 ; i<nsy/scale ; i++ )
  {	ii=i+ofy;
	for ( j=0 ; j<nsx/scale ; j++ )
	 {	jj=j+ofx;
		m=0;
		if ( ii<0 || jj<0 || ii>sy-1 || jj>sx-1 )
		 {	for ( k=0 ; k<scale ; k++ )
			 {  for ( l=0 ; l<scale ; l++ )
			     {	wret[k][l]=0.0;			}
			 }
			m|=MASK_OUTER;
		 }
		else if ( mask[ii][jj] )
		 {	for ( k=0 ; k<scale ; k++ )
			 {  for ( l=0 ; l<scale ; l++ )
			     {	wret[k][l]=0.0;			}
			 }
			m|=mask[ii][jj];
		 }
		else 
	 	 {	biquad_isc_int_block_subpixels(bqc,jj,ii,scale,scale,wret);
			m=0;
		 }

		for ( k=0 ; k<scale ; k++ )
		 {  for ( l=0 ; l<scale ; l++ )
		     {	outimg->data[i*scale+k][j*scale+l]=wret[k][l];
			outmask[i*scale+k][j*scale+l]=m;
		     }
		 }
	 }
	logmsg(is_verbose>=2,"\rZooming... (%d/%d)... ",i+1,nsy/scale);
  }
 logmsg(is_verbose,"done.\n");

 tensor_free(bqc);
 tensor_free(wret);

 return(0); 
}

int zoom_raw_image
	(fitsimage *img,char **mask,fitsimage *outimg,char **outmask,
	int ofx,int ofy,int scale)
{
 int	sx,sy,nsx,nsy;
 int	i,j,k,l,m,ii,jj;
 double	wr,sfactor;

 if ( img==NULL || outimg==NULL )	return(1);
 if ( mask==NULL || outmask==NULL )	return(1);
 sx=img->sx,
 sy=img->sy;
 nsx=outimg->sx;
 nsy=outimg->sy;
 if ( sx==0 || sy==0 || nsx==0 || nsy==0 )	return(0);

 sfactor=1.0/(double)(scale*scale);

 for ( i=0 ; i<nsy/scale ; i++ )
  {	ii=i+ofy;
	for ( j=0 ; j<nsx/scale ; j++ )
	 {	jj=j+ofx;
		m=0;
		if ( ii<0 || jj<0 || ii>sy-1 || jj>sx-1 )
		 {	wr=0.0;
			m|=MASK_OUTER;
		 }
		else if ( mask[ii][jj] )
		 {	wr=0.0;
			m|=mask[ii][jj];
		 }
		else 
	 	 {	wr=img->data[ii][jj]*sfactor;
			m=0;
		 }

		for ( k=0 ; k<scale ; k++ )
		 {  for ( l=0 ; l<scale ; l++ )
		     {	outimg->data[i*scale+k][j*scale+l]=wr;
			outmask[i*scale+k][j*scale+l]=m;
		     }
		 }
	 }
  }

 return(0); 
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int shrink_image
	(fitsimage *img,char **mask,fitsimage *outimg,char **outmask,
	int ofx,int ofy,int scale,int do_median,int do_avg_mode)
{
 int	sx,sy,nsx,nsy;
 int	i,j,k,l,ii,jj,m,t;
 double	w,*medarr;

 if ( img==NULL || outimg==NULL )	return(1);
 if ( mask==NULL || outmask==NULL )	return(1);
 sx=img->sx,
 sy=img->sy;
 nsx=outimg->sx;
 nsy=outimg->sy;
 if ( sx==0 || sy==0 || nsx==0 || nsy==0 )	return(0);

 if ( do_median )
	medarr=(double *)malloc(sizeof(double)*(scale*scale));
 else
	medarr=NULL;

 for ( i=0 ; i<nsy ; i++ )
  {	ii=i*scale+ofy;
	for ( j=0 ; j<nsx ; j++ )
	 {	jj=j*scale+ofx;
		if ( ii<0 || jj<0 || ii>sy-scale || jj>sx-scale )
		 {	outimg->data[i][j]=0.0;
			outmask[i][j]=MASK_OUTER;
		 }
		else if ( medarr != NULL )
		 {	m=0;
			if ( do_avg_mode==0 )
			 {	for ( k=0,t=0 ; k<scale ; k++ )
				 {  for ( l=0 ; l<scale ; l++ )
				     {	m|=mask[ii+k][jj+l];
					medarr[t]=img->data[ii+k][jj+l];
					t++;
				     }
				 }
				outmask[i][j]=m;
				outimg->data[i][j]=median(medarr,t)*(double)(scale*scale);
			 }
			else if ( do_avg_mode==1 )
			 {	for ( k=0,t=0 ; k<scale ; k++ )
				 {  for ( l=0 ; l<scale ; l++ )
				     {	m|=~mask[ii+k][jj+l];
					if ( ! mask[ii+k][jj+l] )
					 {	medarr[t]=img->data[ii+k][jj+l];
						t++;
					 }
				     }
				 }
				outmask[i][j]=(~m)&MASK_ALL;
				outimg->data[i][j]=median(medarr,t)*(double)(scale*scale);
			 }
			else
			 {	outmask[i][j]=MASK_ALL;
				outimg->data[i][j]=0.0;
			 }
		 }
		else
		 {	w=0.0;m=0;
			for ( k=0 ; k<scale ; k++ )
			 {  for ( l=0 ; l<scale ; l++ )
			     {	m|=mask[ii+k][jj+l];
				w+=img->data[ii+k][jj+l];
			     }
			 }
			outmask[i][j]=m;
			outimg->data[i][j]=w;
		 }
	 }
  }

 if ( medarr != NULL )
	free(medarr);

 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int scatter_image
	(fitsimage *img,char **mask,fitsimage *outimg,char **outmask)
{
 int	sx,sy,nsx,nsy;
 int	i,j;
 double	**bqc;

 if ( img==NULL || outimg==NULL )	return(1);
 if ( mask==NULL || outmask==NULL )	return(1);
 sx=img->sx,
 sy=img->sy;
 nsx=outimg->sx;
 nsy=outimg->sy;
 if ( sx==0 || sy==0 || nsx==0 || nsy==0 )	return(0);
 if ( sx != nsx || sy != nsy )			return(1);

 bqc =(double **)tensor_alloc_2d(double,2*sx+1,2*sy+1);
 logmsg(is_verbose,"Calculating biquadratic coefficients... ");
 biquad_coeff(img->data,sx,sy,bqc,mask);
 logmsg(is_verbose,"done.\n");

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( ! mask[i][j] )
			outmask[i][j]=0,outimg->data[i][j]=biquad_scatter(bqc,j,i);
		else
			outmask[i][j]=mask[i][j],outimg->data[i][j]=0.0;
	 }
  }

 tensor_free(bqc);

 return(0); 
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int trim_image
	(fitsimage *img,char **mask,fitsimage *outimg,char **outmask,
	int ofx,int ofy)
{
 int	sx,sy,nsx,nsy;
 int	i,j,ii,jj,m;

 if ( img==NULL || outimg==NULL )	return(1);
 if ( mask==NULL || outmask==NULL )	return(1);
 sx=img->sx,
 sy=img->sy;
 nsx=outimg->sx;
 nsy=outimg->sy;
 if ( sx==0 || sy==0 || nsx==0 || nsy==0 )	return(0);

 for ( i=0 ; i<nsy ; i++ )
  {	ii=i+ofy;
	for ( j=0 ; j<nsx ; j++ )
	 {	jj=j+ofx;
		m=0;
		if ( ii<0 || jj<0 || ii>=sy || jj>=sx )
		 {	m|=MASK_OUTER,
			outimg->data[i][j]=0.0;
		 }
		else
		 {	m|=mask[ii][jj];
			outimg->data[i][j]=img->data[ii][jj];
		 }
		outmask[i][j]=m;
	 }
  }

 return(0);
}

/*****************************************************************************/

int main(int argc,char *argv[])
{
 FILE		*fw,*fr;
 int		i,is_help,is_flip_x,is_flip_y,is_invert,is_shift;
 fits		*img,*outimg;
 int		ofx,ofy,nsx,nsy,frameno,layer,do_median,do_avg_mode;
 char		*outimgfile,*inimgfile,*transfile,*transparam,*basename,
		**inmasklist,*inweightfile;
 int		zratio,gratio,sratio,bqscatt,int_method;
 char		**mask,**outmask,*fdpstring;
 char		*smoothparam;
 char		*explode_basename;
 transformation	trf_data,*trf=&trf_data;
 fitsdataparam	fdp;
 weightlist	wl_data,*wl;
 smooth		sp;
 double		shift_dx,shift_dy;

 progbasename=strrchr(argv[0],'/');
 if ( progbasename != NULL )	progbasename++;
 else				progbasename=argv[0];

 outimgfile=inimgfile=transfile=inweightfile=transparam=NULL;inmasklist=NULL;
 is_comment=is_verbose=is_help=0;zratio=sratio=gratio=bqscatt=0;
 ofx=ofy=is_flip_x=is_flip_y=nsx=nsy=0;int_method=is_invert=0;
 do_median=do_avg_mode=0;
 frameno=0;
 smoothparam=NULL;
 explode_basename=NULL;

 is_shift=0;
 shift_dx=shift_dy=0.0;

 fdp.bitpix=0;fdp.is_scale=0;fdp.bscale=1.0;fdp.bzero=0.0;fdpstring=NULL;

 sp.type=SMOOTH_NONE;
 sp.filter=SMOOTH_PREFILTER_NONE;

 layer=-1;

 i=scanarg(argc,argv,SCANARG_ALLOW_FLAGS,
	"--version:%NS-1f%q",&is_help,
        "--version-short|--short-version:%NS-2f%q",&is_help,
        "--long-help|--help-long:%SN2f%q",&is_help,
	"-h|--help:%f%q",&is_help,
	"-i|--input:%s",&inimgfile,
	"--frame:%d",&frameno,
	"-o|--output:%s",&outimgfile,
	"-t|--transformation:%s",&transparam,
	"-T|--input-transformation:%s",&transfile,
	"-W|--input-weight:%s",&inweightfile,
        "--reverse|--inverse:%f",&is_invert,
	"-m:%SN0f",&int_method,		/* interpolate + steplike	*/ 
	"-l:%SN1f",&int_method,		/* integrate   + steplike	*/
	"-c:%SN2f",&int_method, 	/* interpolate + spline		*/
	"-k:%SN3f",&int_method,		/* integrate   + spline		*/
	"-M|--input-mask:%t",&inmasklist,
	"-z|--zoom:%d",&zratio,
	"-g|--magnify:%d",&gratio,
	"-r|--shrink:%d",&sratio, 
	"-d|--median:%f",&do_median,
	"--optimistic-masking:%SN1f",&do_avg_mode,
	"-n|--noise:%f",&bqscatt,
	"-a|--smooth:%s",&smoothparam,
	"-f|--offset:%d,%d",&ofx,&ofy,
	"-s|--size:%d,%d",&nsx,&nsy,
	"-e|--shift:%f%g,%g",&is_shift,&shift_dx,&shift_dy,
	"-y|--layer|--first-layer:%d",&layer,
	"-x|--explode:%s",&explode_basename,
	"--flip-x:%f",&is_flip_x,
	"--flip-y:%f",&is_flip_y,
	"--flip-xy|--flip-yx:%f%f",&is_flip_x,&is_flip_y,
	"-b|--bitpix:%d",&fdp.bitpix, 
	"-D|--data:%s",&fdpstring,
	"--comment:%f",&is_comment,"(C):%f",&is_comment,
	"--verbose:%i",&is_verbose,"(V):%i",&is_verbose,
	"-:%w",&inimgfile,
	"-*|+*:%e",
	"*:%w",&inimgfile,
	NULL);

 if ( i )		
  {	fprint_error("invalid command line argument near '%s'",argv[i]);
	return(1);
  }
 else if ( is_help<0 )
  {	fprint_generic_version(stdout,argv[0],"fitrans",FI_TRANS_VERSION,is_help);
	return(0);
  }
 else if ( is_help>1 )
  {	fprint_fitrans_long_help(stdout);
	return(0);
  }
 else if ( is_help )
  {	fprint_fitrans_usage(stdout);
	return(0);
  }

 is_invert=!is_invert;	/* negate to be self-consistent with grtrans!	*/
			/* [inserted: 2006.05.09]			*/

 if ( parse_fits_data_param(fdpstring,&fdp) )
  {	fprint_error("invalid pixel data format");
	return(1);
  }

 if ( smoothparam != NULL )
  {	int	order,hsize;
	double	sigma;

	order=-1;
	sigma=-1.0;
	hsize=-1;

	sp.type=SMOOTH_NONE;
	sp.xorder=sp.yorder=0;
	sp.filter=SMOOTH_PREFILTER_NONE;
	sp.fxhsize=sp.fyhsize=0;
	sp.niter=0;
	sp.lower=3.0;
	sp.upper=3.0;
	sp.is_mean_unity=0;
	sp.is_detrend=0;

	i=scanpar(smoothparam,SCANPAR_DEFAULT,
		"spline:"	SNf(SMOOTH_SPLINE), &sp.type,
		"polynomial:" 	SNf(SMOOTH_POLYNOMIAL), &sp.type,
		"order:%d",	&order,
		"xorder:%d",	&sp.xorder,
		"yorder:%d",	&sp.yorder,
		"mean:"		SNf(SMOOTH_PREFILTER_MEAN), &sp.filter,
		"median:"	SNf(SMOOTH_PREFILTER_MEDIAN), &sp.filter,
		"hsize:%d",	&hsize,
		"xhsize:%d",	&sp.fxhsize,
		"yhsize:%d",	&sp.fyhsize,
		"rejection:%g",	&sp.frejratio,
		"iterations:%d",&sp.niter,
		"lower:%g",	&sp.lower,
		"upper:%g",	&sp.upper,
		"sigma:%g",	&sigma,
		"unity:%f",	&sp.is_mean_unity,
		"detrend:%f",	&sp.is_detrend,
		NULL);
	if ( i )
	 {	fprint_error("invalid smooth mode parameter in '%s'.",smoothparam);
		return(1);
	 }
	if ( order>=0 )		sp.xorder=sp.yorder=order;
	if ( sigma>=0.0 )	sp.lower =sp.upper =sigma;
	if ( hsize>=0 )		sp.fxhsize=sp.fyhsize=hsize;
  }

 if ( inweightfile != NULL )
  {	fits	*img;
	if ( (fr=fopenread(inweightfile))==NULL )
	 {	fprint_error("unable to open input weight file '%s'.",inweightfile);
		return(1);
	 }
	img=fits_read(fr);
	if ( img==NULL )
	 {	fprint_error("unable to interpret input weight file as FITS data.");
		return(1);
	 }
	wl=&wl_data;
	if ( weight_parse_fits(img,wl) )
	 {	fprint_error("unable to interpret input weight file properly.");
		return(1);
	 }
	fits_free(img);
  }
 else
	wl=NULL;

 if ( explode_basename != NULL )
  {	int		i,l,nlayer,sx,sy;
	fitsimage	*fi;

	if ( ! format_check_if_formatted(explode_basename,"dxXio") ) 
	 {	fprint_error("invalid explode basename template '%s' (perhaps '%%d' is missing)",explode_basename);
		return(1);
	 }
	basename=fits_basename(inimgfile,&frameno);
	fr=fopenread(basename);
	if ( fr==NULL )	
	 {	fprint_error("unable to access/read input file");
		return(1);
	 }
	if ( (img=fits_seek_frame_to_image(fr,0)) == NULL )
	 {	fprint_error("unable to parse input file as a FITS image");
		return(1);
	 }
	if ( layer<0 )	layer=0;
	
	fi=&img->i;

	nlayer=1;
	for ( l=2 ; l<fi->dim ; l++ )
	 {	nlayer *= fi->naxis[l];		}

	sx=fi->naxis[0];
	if ( fi->dim>1 )
		sy=fi->naxis[1];
	else
		sy=1;

	for ( l=0 ; l<nlayer ; l++ )
	 {	char	*filename;
		char	buff[256];
		fits	*img;
		FILE	*fw;

		filename=format_replace(explode_basename,0,                    
                        'd',FORMAT_INTEGER,l+layer,
                        'i',FORMAT_INTEGER,l+layer,
                        'x',FORMAT_HEX_INTEGER,l+layer,
                        'X',FORMAT_HXC_INTEGER,l+layer,
                        'o',FORMAT_OCT_INTEGER,l+layer,
			0);

		img=fits_create();
		fits_set_standard(img,NULL);
		fits_alloc_image(img,sx,sy);
		img->i.bit=fi->bit;
		img->i.curr.bscale=fi->curr.bscale;
		img->i.curr.bzero=fi->curr.bzero;
		img->i.read.bscale=fi->read.bscale;
		img->i.read.bzero=fi->read.bzero;
		fits_set_image_params(img);
		sprintf(buff,"Created by fitrans, version: %s",FI_TRANS_VERSION);
		fits_set_origin(img,buff,NULL);
		fits_history_export_command_line(img,"fitrans",FI_TRANS_VERSION,argc,argv);

		for ( i=0 ; i<sy ; i++ )
		 {	fits_read_image_line(fr,sx,fi->bit,img->i.data[i]);	}

		fw=fopenwrite(filename);
		if ( fw==NULL )			
		 {	fprint_error("unable to create output file '%s'",filename);
			return(1);
		 }
	 	fits_write(fw,img); 
		fclosewrite(fw);

		fits_free(img);

		free(filename);
	 }

	fcloseread(fr);
	return(0);
  }

 basename=fits_basename(inimgfile,&frameno);
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
 if ( img->i.dim != 2 && layer<0 )
  {	fprint_error("image dimension differs from 2 and no layer has been specified.");
	return(1); 
  }
 else if ( img->i.dim == 3 && layer>=0 && layer<img->i.naxis[2] )
  {	fitsimage	fi;
	int		i,sx,sy;
	double		***dcube;
	sx=img->i.naxis[0];
	sy=img->i.naxis[1];
	fits_image_alloc(&fi,sx,sy);
	fi.bit=img->i.bit;
	fi.curr=img->i.curr;
	fi.read=img->i.read;
	dcube=(double ***)img->i.vdata;
	for ( i=0 ; i<sy ; i++ )
	 {	memcpy(fi.data[i],dcube[layer][i],sizeof(double)*sx);	}
	fits_image_free(&img->i);
	img->i=fi;
	fits_image_set_params(&img->header,&img->i);
  }
 fits_rescale(img);
 
 logmsg(is_verbose,"[%d,%d]\n",img->i.sx,img->i.sy);

/* Create the mask for the valid pixels: */
 mask=fits_mask_read_from_header(&img->header,img->i.sx,img->i.sy,NULL);
 if ( inmasklist != NULL )
  {	if ( join_masks_from_files(mask,img->i.sx,img->i.sy,inmasklist) )
	 {	fprint_error("unable to read one of the input mask files");
		return(1);
	 }
  }
 fits_mask_mark_nans(&img->i,mask,MASK_NAN);

 if ( nsx<0 )	is_flip_x=!is_flip_x,nsx=-nsx;
 if ( nsy<0 )	is_flip_y=!is_flip_y,nsy=-nsy;

 if ( nsx==0 && nsy==0 )
  {	if ( zratio>1 )
		nsx=img->i.sx*zratio,
		nsy=img->i.sy*zratio;
	else if ( gratio>1 )
		nsx=img->i.sx*gratio,
		nsy=img->i.sy*gratio;
	else if ( sratio>1 )
		nsx=(img->i.sx+sratio-1)/sratio,
		nsy=(img->i.sy+sratio-1)/sratio;
	else
		nsx=img->i.sx,
		nsy=img->i.sy;
  }

/* Create new FITS image for the output: */
 outimg=fits_duplicate_empty(img);
 if ( nsx != img->i.sx || nsy != img->i.sy )
  {	fits_free_image(outimg);
	fits_alloc_image(outimg,nsx,nsy);
	fits_reset_image(outimg);
	fits_set_header_integer(outimg,"NAXIS1",FITS_SH_FIRST,nsx,NULL);
	fits_set_header_integer(outimg,"NAXIS2",FITS_SH_FIRST,nsy,NULL);
  }
 outmask=fits_mask_create_empty(nsx,nsy);

/* Read transformation file, if it has been specified in the command line: */
 if ( transfile != NULL )
  {	FILE	*ft;
	ft=fopenread(transfile);
	if ( ft==NULL )	
	 {	fprint_error("unable to open input transformation file '%s'",transfile);
		return(1);
	 }
	i=transformation_read_data(ft,trf);
	if ( i )
	 {	fprint_error("unable to parse the contents of input transformation file '%s'",transfile);
		return(1);
	 }
	fcloseread(ft);
  }
 else if ( transparam != NULL )
  {	i=transformation_parse_params(transparam,trf);
	if ( i )		
	 {	fprint_error("unable to parse transformation string '%s'",transparam);
		return(1);
	 }
  }
 else if ( is_shift )
  {	trf->type=TRANS_POLYNOMIAL;
	trf->order=1;
	trf->nval=2;
	trf->ox=trf->oy=0.0;
	trf->scale=1.0;
	trf->bshx=trf->bshy=0.0;
 	trf->vfits=(double **)tensor_alloc_2d(double,3,2);
	trf->vfits[0][0]=shift_dx;
	trf->vfits[0][1]=1.0;
	trf->vfits[0][2]=0.0;
	trf->vfits[1][0]=shift_dy;
	trf->vfits[1][1]=0.0;
	trf->vfits[1][2]=1.0;
  }
 else
	trf=NULL;

 if ( trf != NULL && trf->nval != 2 )
  {	fprint_error("input transformation is not an 2D -> 2D transformation");
	return(1);
  }


/* Apply the required transformation: */
 if ( trf != NULL )
	apply_transformation(&img->i,mask,&outimg->i,outmask,ofx,ofy,trf,int_method,is_invert);
 else if ( zratio>1 )
	zoom_image(&img->i,mask,&outimg->i,outmask,ofx,ofy,zratio);
 else if ( gratio>1 )
	zoom_raw_image(&img->i,mask,&outimg->i,outmask,ofx,ofy,gratio);
 else if ( sratio>1 )
	shrink_image(&img->i,mask,&outimg->i,outmask,ofx,ofy,sratio,do_median,do_avg_mode);
 else if ( bqscatt )
	scatter_image(&img->i,mask,&outimg->i,outmask);
 else if ( sp.type || sp.filter )
	combined_smooth_image(&img->i,mask,&outimg->i,outmask,&sp);
 else 
	trim_image(&img->i,mask,&outimg->i,outmask,ofx,ofy);

 fits_mask_free(mask);
 if ( fdp.bitpix )
  {	outimg->i.bit=fdp.bitpix;		}
 if ( fdp.is_scale )	
  {	outimg->i.read.bscale=fdp.bscale,
	outimg->i.read.bzero=fdp.bzero;
  }
 else if ( outimg->i.bit<0 )
  {	outimg->i.read.bscale=1.0,
	outimg->i.read.bzero=0.0;
  }
 fits_set_image_params(outimg);
 fits_backscale(outimg,outimg->i.read.bscale,outimg->i.read.bzero);

 mark_integerlimited_pixels(&outimg->i,outmask,outimg->i.bit,1,MASK_OVERSATURATED,MASK_OVERSATURATED);
 
/* Write the output image: */
 if ( is_flip_y )
  {	double	*pd;
	char	*pc;
	int	i;
	for ( i=0 ; i<nsy/2 ; i++ )
	 {	pd=outimg->i.data[i],
		outimg->i.data[i]=outimg->i.data[nsy-i-1],
		outimg->i.data[nsy-i-1]=pd;
		pc=outmask[i],
		outmask[i]=outmask[nsy-i-1],
		outmask[nsy-i-1]=pc;
	 }
  }
 if ( is_flip_x )
  {	double	wd;
	int	i,j,wc;
	for ( i=0 ; i<nsy ; i++ )
	 {	for ( j=0 ; j<nsx/2 ; j++ )
		 {	wd=outimg->i.data[i][j],
			outimg->i.data[i][j]=outimg->i.data[i][nsx-1-j];
			outimg->i.data[i][nsx-1-j]=wd;
			wc=outmask[i][j],
			outmask[i][j]=outmask[i][nsx-1-j],
			outmask[i][nsx-1-j]=wc;
		 }
	 }
  }

 fits_history_export_command_line(outimg,"fitrans",FI_TRANS_VERSION,argc,argv);

 fits_mask_export_as_header(&outimg->header,1,outmask,nsx,nsy,NULL);
 if ( outimgfile==NULL )	fw=stdout;
 else				fw=fopenwrite(outimgfile);
 if ( fw==NULL )		
  {	fprint_error("unable to create output image file '%s'",outimgfile);
	return(1);
  }

 fits_write(fw,outimg);
 fclosewrite(fw);

 return(0);
}

/*****************************************************************************/

