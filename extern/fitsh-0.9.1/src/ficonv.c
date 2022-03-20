/*****************************************************************************/
/* ficonv.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line user interface for (FITS) image convolution and subtraction. */
/*****************************************************************************/
#define	FI_CONV_VERSION		"0.4"
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include <fits/fits.h>

#include "longhelp.h"
#include "fi.h"

#include "math/spline/spline.h"
#include "math/spline/biquad.h"
#include "math/spline/bicubic.h"
#include "statistics.h"
#include "math/fit/lmfit.h"
#include "math/poly.h"
#include "fitsmask.h"
#include "io/iof.h"
#include "io/scanarg.h"
#include "io/tokenize.h"

#include "tensor.h"
#include "common.h"
#include "kernel.h"
#include "history.h"

/*****************************************************************************/

int	is_comment,	/* write comments to some of the output files or not */
	is_verbose;	/* verbosity level				     */
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

#define		FM_NORMAL		0x00
#define		FM_MASKED		0x01
#define		FM_WEIGHTED		0x02
#define		FM_BGITERATIVE		0x04

/*****************************************************************************/

int affine_img_transformation_spline_fit(fitsimage *ref,fitsimage *img,char **mask,int bx,int by,double *coeff)
{
 int	i,j,k,l,v,sx,sy,nvar;
 double	**xbase,**ybase;
 double	**xnspc,**ynspc;
 double	*tmp,x,y;
 double	**amatrix,*bvector,*fvars,f,w;

 if ( ref==NULL || img==NULL || coeff==NULL )
	return(-1);
 else if ( (sx=ref->sx) != img->sx || (sy=ref->sy) != img->sy )
	return(-1);
 else if ( bx<=0 || by<=0 )
	return(-1);

 nvar=2*(bx+1)*(by+1);

 xbase=tensor_alloc_2d(double,sx,bx+1);
 ybase=tensor_alloc_2d(double,sy,by+1);
 xnspc=tensor_alloc_2d(double,bx+1,bx+1);
 ynspc=tensor_alloc_2d(double,by+1,by+1);
 tmp  =tensor_alloc_1d(double,(bx>by?bx:by)+1);

 for ( i=0 ; i<=by ; i++ )
  {	for ( k=0 ; k<=by ; k++ )
	 {	tmp[k]=0.0;		}
	tmp[i]=1.0;
	natspline_coeff(tmp,by+1,ynspc[i]);
	for ( k=0 ; k<sy ; k++ )
	 {	y=(double)by*((double)k+0.5)/(double)sy;
		ybase[i][k]=natspline_inter(tmp,ynspc[i],by+1,y);
	 }
  }

 for ( j=0 ; j<=bx ; j++ )
  {	for ( l=0 ; l<=bx ; l++ )
	 {	tmp[l]=0.0;		}
	tmp[j]=1.0;
	natspline_coeff(tmp,bx+1,xnspc[j]);
	for ( l=0 ; l<sx ; l++ )
	 {	x=(double)bx*((double)l+0.5)/(double)sx;
		xbase[j][l]=natspline_inter(tmp,xnspc[j],bx+1,x);
	 }
  }

 amatrix=matrix_alloc(nvar);
 bvector=vector_alloc(nvar);
 fvars  =vector_alloc(nvar);

 for ( i=0 ; i<nvar ; i++ )
  {	for ( j=0 ; j<nvar ; j++ )
	 {	amatrix[i][j]=0.0;		}
	bvector[i]=0.0;
  }

 for ( k=0 ; k<sy ; k++ )
  {	for ( l=0 ; l<sx ; l++ )
	 {	if ( mask != NULL && mask[k][l] )
			continue;
		v=0;
		w=ref->data[k][l];
		for ( i=0 ; i<=by ; i++ )
		 {	for ( j=0 ; j<=bx ; j++ )
			 {	f=xbase[j][l]*ybase[i][k];
				fvars[v+0]=f*w;
				fvars[v+1]=f;
				v+=2;
			 }
		 }
		for ( i=0 ; i<nvar ; i++ )
		 {	for ( j=0 ; j<nvar ; j++ )
			 {	amatrix[i][j] += fvars[i]*fvars[j];	}
			bvector[i] += fvars[i]*img->data[k][l];
		 }
	 }
  }

 solve_gauss(amatrix,bvector,nvar);

 for ( i=0 ; i<nvar ; i++ )
	coeff[i]=bvector[i];

 vector_free(fvars);
 vector_free(bvector);
 matrix_free(amatrix);

 tensor_free(tmp);
 tensor_free(ynspc);
 tensor_free(xnspc);
 tensor_free(ybase);
 tensor_free(xbase);

 return(0);
}

int affine_img_transformation_spline_eval(fitsimage *ref,char **mask,int bx,int by,double *coeff)
{
 double	**bcscl,**bcoff,**tbscl,**tboff;
 int	i,j,k,l,v,sx,sy;
 double	scl,off,x,y;

 bcscl=(double **)tensor_alloc_2d(double,2*(bx+1),2*(by+1));
 bcoff=(double **)tensor_alloc_2d(double,2*(bx+1),2*(by+1));
 tbscl=(double **)tensor_alloc_2d(double,bx+1,by+1);
 tboff=(double **)tensor_alloc_2d(double,bx+1,by+1);

 v=0;
 for ( i=0 ; i<=by ; i++ )
  {	for ( j=0 ; j<=bx ; j++ )
	 {	tbscl[i][j]=coeff[v+0];
		tboff[i][j]=coeff[v+1];
		v+=2;
	 }
  }

 bicubic_coeff(tbscl,bx+1,by+1,bcscl,NULL);
 bicubic_coeff(tboff,bx+1,by+1,bcoff,NULL);

 tensor_free(tboff);
 tensor_free(tbscl);

 sx=ref->sx;
 sy=ref->sy;

 for ( k=0 ; k<sy ; k++ )
  {	y=(double)by*((double)k+0.5)/(double)sy;

	for ( l=0 ; l<sx ; l++ )
	 {	if ( mask != NULL && mask[k][l] )
			continue;

		x=(double)bx*((double)l+0.5)/(double)sx;

		scl=bicubic_inter(bcscl,x,y);
		off=bicubic_inter(bcoff,x,y);
	
		ref->data[k][l]=ref->data[k][l]*scl+off;
	 }
  }
 
 tensor_free(bcoff);
 tensor_free(bcscl);

 return(0);
}

/*****************************************************************************/

char ** create_level_mask(fitsimage *img,char **inmask,double mlev1,double klev1,double mlev2,double klev2)
{
 char	**mask;
 double	**bqc,*ps,sc,th1,th2;
 int	i,j,sx,sy,np;

 sx=img->sx,
 sy=img->sy;
 bqc=(double **)tensor_alloc_2d(double,2*sx+1,2*sy+1);
 biquad_coeff(img->data,sx,sy,bqc,inmask);
 ps=tensor_alloc_1d(double,sx*sy);
 mask=fits_mask_duplicate(inmask,sx,sy);
 np=0;
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( mask[i][j] )	continue;
		sc=biquad_scatter(bqc,j,i);
		bqc[2*i+1][2*j+1]=ps[np]=sc;
		np++;
	 }
  }
 median(ps,np);

 if ( klev1>0.0 ) i=(int)((double)(np-1)*mlev1),th1=ps[i]*klev1;
 else		  th1=0.0;
 if ( klev2>0.0 ) i=(int)((double)(np-1)*mlev2),th2=ps[i]*klev2;
 else		  th2=0.0;

 tensor_free(ps);

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( mask[i][j] )	continue;
		sc=bqc[2*i+1][2*j+1];
		if ( sc<th1 && th1>0.0 )	mask[i][j]=-1;
		else if ( sc>th2 && th2>0.0 )	mask[i][j]=+1;
		else				mask[i][j]= 0;
	 }
  }	
 tensor_free(bqc);
 return(mask);
}

char ** create_foreground_mask(fitsimage *img,char **inmask,double mlev,double klev)
{
 char	**ret;
 ret=create_level_mask(img,inmask,0.0,0.0,mlev,klev);
 return(ret); 
}
char ** create_background_mask(fitsimage *img,char **inmask,double mlev,double klev)
{
 char	**ret;
 int	i,j;
 ret=create_level_mask(img,inmask,mlev,klev,0.0,0.0);
 for ( i=0 ; i<img->sy ; i++ )
  {	for ( j=0 ; j<img->sx ; j++ )
	 {	if ( ret[i][j]<0 )	ret[i][j]=1;	}
  }
 return(ret); 
}

double **create_image_weight(fitsimage *img,char **inmask)
{
 double	**bqc,**wght;
 int	i,j,sx,sy;
 sx=img->sx,sy=img->sy;
 bqc=(double **)tensor_alloc_2d(double,2*sx+1,2*sy+1);
 biquad_coeff(img->data,sx,sy,bqc,inmask);
 wght=(double **)tensor_alloc_2d(double,sx,sy);
 if ( inmask==NULL )
  {	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	wght[i][j]=biquad_scatter(bqc,j,i);	 }
	 }
  }
 else
  {	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	if ( inmask[i][j] )
				wght[i][j]=0.0;
			else
				wght[i][j]=biquad_scatter(bqc,j,i);
		 }	
	 }
  }
 tensor_free(bqc);
 return(wght);
}

/*****************************************************************************/

int mark_stamps(fitsimage *img,rectan *rcts,int nrct)
{
 int	i;
 rectan	*r; 
 double	color;

 color=10000.0;
 for ( i=0 ; i<nrct ; i++ )
  {	r=&rcts[i];
	fits_image_draw_line(img,r->x1  ,r->y1  ,r->x2-1,r->y1  ,color,0x3333);
	fits_image_draw_line(img,r->x2-1,r->y1  ,r->x2-1,r->y2-1,color,0x3333);
	fits_image_draw_line(img,r->x2-1,r->y2-1,r->x1  ,r->y2-1,color,0x3333);
	fits_image_draw_line(img,r->x1  ,r->y2-1,r->x1  ,r->y1  ,color,0x3333);
  }
 return(0);
}

/*****************************************************************************/

int make_subtracted_image(fitsimage *cnv,char **mask,fitsimage *img,char **mask_img,fitsimage *out,char **mask_out,kernellist *xlist)
{
 fitsimage	xkcimg;
 int		i,j,sx,sy;
 char		**imask,**aimask;

 if ( cnv==NULL || cnv->data==NULL )	return(1);
 if ( img==NULL || img->data==NULL )	return(1);
 if ( out==NULL || out->data==NULL )	return(1);
 sx=cnv->sx,sy=cnv->sy;
 
 if ( xlist != NULL && xlist->kernels != 0 && xlist->nkernel>0 )
  {	kernel	*k;
	int	ihsize;

	ihsize=0;
	for ( i=0,k=xlist->kernels ; i<xlist->nkernel ; i++,k++ )
	 {	if ( k->hsize > ihsize ) ihsize=k->hsize;	}
	imask=aimask=fits_mask_expand_false(mask_img,sx,sy,ihsize,-1,-1,1);
	fits_image_duplicate(&xkcimg,img,1);
	convolve_with_kernel_set(img,imask,xlist,&xkcimg);
  }
 else	
  {	xkcimg.data=NULL;
	aimask=NULL;
	imask=mask_img;
  }
		
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( ! mask[i][j] && ! imask[i][j] )
		 {	out->data[i][j]=img->data[i][j]-cnv->data[i][j];
			if ( xkcimg.data != NULL )
				out->data[i][j]+=xkcimg.data[i][j];
			mask_out[i][j]=0;
		 }
		else
		 {	out->data[i][j]=0.0;
			mask_out[i][j]=-1;
		 }
	 }
  }
 if ( xkcimg.data != NULL )	fits_image_free(&xkcimg);
 if ( aimask != NULL )		fits_mask_free(aimask);

 return(0);
}

/*****************************************************************************/

int create_weights(fitsimage *img,char **mask,char **levmask,double **wcc,int sign)
{
 int	i,j,sx,sy;
 sx=img->sx,sy=img->sy;
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( mask[i][j] )			wcc[i][j]=0.0;
		else if ( sign>0 && levmask[i][j]>0 )	wcc[i][j]=1.0;
		else if ( sign<0 && levmask[i][j]<0 )	wcc[i][j]=1.0;
		else if ( sign==0 && levmask[i][j] )	wcc[i][j]=1.0;
		else					wcc[i][j]=0.0;
	 }
  }
 return(0);
}

int fit_kernels_native(fitsimage *ref,fitsimage *img,char **mask_conv,
	kernellist *klist,kernellist *xlist,int method,int bdc)
{
 int	i,j,l,sx,sy;
 kernel	*bgkernel,*k;

 sx=ref->sx,sy=ref->sy;

 bgkernel=NULL;
 for ( i=0,k=klist->kernels ; i<klist->nkernel ; i++,k++ )
  {	if ( k->type==KERNEL_BACKGROUND )	bgkernel=k;
	k->flag=1;
  }
 if ( xlist != NULL )
  {	for ( i=0,k=xlist->kernels ; i<xlist->nkernel ; i++,k++ )
	 {	k->flag=1;	}
  }

 if ( method & FM_WEIGHTED )	method|=FM_MASKED;	/* -w implies -m ! */

 if ( bgkernel != NULL && ( method & FM_MASKED ) )
  {	double		**wcc,w;
	kernel		*k;
	char		**fgmask;
	
	wcc=tensor_alloc_2d(double,sx,sy);

	fgmask=create_foreground_mask(ref,mask_conv,0.5,3.0);

	for ( i=0 ; i<klist->nkernel ; i++ )
	 {	k=&klist->kernels[i];k->flag=1;		}

	l=0;
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	if ( mask_conv[i][j] )
				wcc[i][j]=0.0;
			else if ( fgmask[i][j] )
		 	 {	if ( method & FM_WEIGHTED )
					w=ref->data[i][j];
				else	
					w=1.0;					
				wcc[i][j]=w;
				l++;
			 }
			else	wcc[i][j]=0.0;
		 }
	  }		
	logmsg(is_verbose>=1,"Foreground pixels: %d/%d\n",l,sx*sy);
	fit_kernel_poly_coefficients_block(ref,img,mask_conv,wcc,bdc,bdc,klist,xlist);
	fits_mask_free(fgmask);
	tensor_free(wcc);
  }

 else
 	fit_kernel_poly_coefficients_block(ref,img,mask_conv,NULL,bdc,bdc,klist,xlist);

 return(0);
}


int fit_kernels(fitsimage *ref,char **mask_ref,fitsimage *img,char **mask_img,
	char **inmask,kernellist *klist,kernellist *xlist,int method,int bdc,
	int niter,double rlevel,double gain)
{
 int	hsize,i,sx,sy;
 char	**mask,**mask_conv;
 kernel	*k;

 if ( ref==NULL || img==NULL )		return(1);
 if ( ref->data==NULL || img->data==NULL )	return(1);
 sx=ref->sx,sy=ref->sy;
 if ( sx != img->sx || sy != img->sy )		return(1);

 if ( xlist==NULL || xlist->nkernel<=0 || xlist->kernels==NULL )
	xlist=NULL;

 logmsg(is_verbose>=1,"Fitting kernel coefficients ...\n");

 hsize=0;
 for ( i=0,k=klist->kernels ; i<klist->nkernel ; i++,k++ )
  {	if ( k->hsize > hsize )		hsize=k->hsize;		}
 if ( xlist != NULL )
  {	for ( i=0,k=xlist->kernels ; k != NULL && i<xlist->nkernel ; i++,k++ )
	 {	if ( k->hsize > hsize )	hsize=k->hsize;		}
  }

 mask=fits_mask_create_empty(sx,sy);
 fits_mask_and(mask,sx,sy,mask_img);
 fits_mask_and(mask,sx,sy,mask_ref);
 if ( inmask != NULL )
 	fits_mask_and(mask,sx,sy,inmask);
 mask_conv=fits_mask_expand_false(mask,sx,sy,hsize,-1,-1,1);

 if ( niter<=0 )	niter=1;
 else			niter++;

 while ( niter>0 )
  {	fit_kernels_native(ref,img,mask_conv,klist,xlist,method,bdc);
	niter--;
	if ( niter>0 )
	 {	fitsimage	sub;
		int		i,j,rej,tot;
		char		**ms;
		double		s,s2,w,n,sig,nz;

		ms=mask_conv;

		fits_image_duplicate(&sub,ref,1);
		convolve_to_subtracted(ref,img,ms,klist,xlist,&sub);
		n=s=s2=0.0;
		for ( i=0 ; i<sy ; i++ )
		 {	for ( j=0 ; j<sx ; j++ )
			 {	if ( ! ms[i][j] )	continue;
				w=sub.data[i][j];
				s+=w,s2+=w*w,n+=1.0;
			 }
		 }
		s/=n,s2/=n;
		sig=s2-s*s;
		fprintf(stderr,"s=%g,sig=%g\n",s,sqrt(sig));
		rej=tot=0;
		for ( i=0 ; i<sy ; i++ )
		 {	for ( j=0 ;j<sx ;  j++ )
			 {	if ( ! ms[i][j] )	continue;
				tot++;
				nz=sqrt(sig+ref->data[i][j]/gain);
				if ( fabs(sub.data[i][j])>rlevel*nz )
					ms[i][j]=0,rej++;
			 }
		 }
		logmsg(is_verbose>=1,"Rejected: %d from %d\n",rej,tot);
		fits_image_free(&sub);
	 }
  };

 tensor_free(mask_conv);
 fits_mask_free(mask);

 return(0);
}

/*****************************************************************************/

int stamp_parse_argument(char *stamparg,char **fitmask,int sx,int sy)
{
 char	*tmpstamparg,**cmd;
 int	k,ix,iy,is,i,j;
 double	x,y,s;

 tmpstamparg=strdup(stamparg);

 cmd=tokenize_char_dyn(tmpstamparg,':');
 for ( k=0 ; cmd != NULL && cmd[k] != NULL ; k++ )
  {	if ( sscanf(cmd[k],"%lg,%lg,%lg",&x,&y,&s)<3 )
		continue;
	if ( s<=0.0 )
		continue;
	ix=(int)x,
	iy=(int)y;
	is=(int)s;
	for ( i=iy-is ; i<=iy+is ; i++ )
	 {	if ( i<0 || i>=sy )	continue;
		for ( j=ix-is ; j<=ix+is ; j++ )
		 {	if ( j<0 || j>=sx )	continue;
			fitmask[i][j] |= 0x80;
		 }
	 }
  }
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( fitmask[i][j] & 0x80 )
			fitmask[i][j] &= 0x7f;
		else
			fitmask[i][j] = (fitmask[i][j] & 0x7f) | MASK_OUTER;
	 }
  }

 if ( cmd != NULL )	free(cmd);
 free(tmpstamparg);

 return(0);
}
int stamp_read_file(FILE *fr,char **fitmask,int sx,int sy)
{
 int	ix,iy,is,i,j;
 double	x,y,s;
 char	*line;

 while ( ! feof(fr) )
  {	line=freadline(fr);
	if ( line==NULL )
		break;
	remove_newlines_and_comments(line);
	if ( sscanf(line,"%lg %lg %lg",&x,&y,&s)<3 )
	 {	free(line);
		continue;
	 }
	free(line);
	if ( s<=0.0 )
		continue;
	ix=(int)x,
	iy=(int)y;
	is=(int)s;
	for ( i=iy-is ; i<=iy+is ; i++ )
	 {	if ( i<0 || i>=sy )	continue;
		for ( j=ix-is ; j<=ix+is ; j++ )
		 {	if ( j<0 || j>=sx )	continue;
			fitmask[i][j] |= 0x80;
		 }
	 }
  }
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( fitmask[i][j] & 0x80 )
			fitmask[i][j] &= 0x7f;
		else
			fitmask[i][j] = (fitmask[i][j] & 0x7f) | MASK_OUTER;
	 }
  }

 return(0);
}


/*****************************************************************************/

int fprint_ficonv_usage(FILE *fw)
{
 fprintf(fw,
"Usage:\tficonv [-C|--comment] [-V|--verbose] [-h|--help]\n"
"\t[-r <reference.fits>] [-i <input-to-fit.fits>\n"
"\t[-m|--masked [-w|--weighted]] [-b|--background-iterative]\n"
"\t[-n|--iterations <niter>] [-s|--rejection-level <rejlevel>]]\n"
"\t{--input-kernel <input-kernels>|-k \"<input-kernel-set>\"}\n"
"\t[--output-kernel <output-kernels>]\n");
 fprintf(fw,
"\t[--input-extra-kernel <pdcnv-kernels>|-x \"<pdcnv-kernel-set>\"]"
"\t[--output-extra-kernel <output-pdcnv>]\n"
"\t[-o|--output <output-convolved.fits>]\n"
"\t[--output-subtracted <output-subtracted.fits>]\n"
"\t[-a <add-to-convolved.fits>] [-d|--divide <blocks>]\n"
"\t[-t|--stamp[s] <x1>,<y1>,<hs1>[:<x2>,<y2>,<hs2>...]]\n"
"\t[--input-stamps <stamplist.dat>]\n");

 return(0);
}

longhelp_entry ficonv_long_help[]=
{
 LONGHELP_OPTIONS,
 { "General options:", NULL },
 { "-h, --help",
	"Gives general summary about the command line options." },
 { "--long-help",
	"Gives a detailed list of command line options." },
 { "--version",
	"Gives some version information about the program." },
 { "-r, --input-reference <image file>",
	"Name of the reference FITS image file." },
 { "-i, --input <image file>",
	"Name of the input FITS image file (required only for kernel fitting)." },
 { "-k, --kernel <kernel set>",
	"List of kernel bases used for fitting convolution kernel." 
	"See also ``Kernel specifications'' below for the format of this "
	"<kernel set> argument." },
 { "--input-kernel <file>",
	"Name of the file containing kernel bases. The kernel bases in this "
	"file should have no associated coefficients if convolution fitting "
	"is done, otherwise the kernel basis file must contain the "
	"convolution coefficients" },
 { "--output-kernel <file>",
	"Name of the file where the coefficients for the kernel bases are "
	"stored after convolution kernel fitting" },
 { "-o, -oc, --output, --output-convolved <image file>", 
	"Name of the output file which is the reference image convolved with "
	"the kernel solution (which can either be a previously fitted and now "
	"read from a file or the result of the current fit)" },
 { "--output-subtracted <image file>",
	"The difference between the input image and the convolved reference "
	"image." },
 { "-a, --add-to <image file>", 
	"This optionally specified file is added to the convolved image" },
 { "-M, --input-mask <fits>",
	"Input mask file to co-add to output image mask." },
 { "-n, --iterations <iterations>",
	"Use an iterative fit with the rejection of the outlier pixels."
	"The maximum number of iterations should be specified with this "
	"command line argument" },
 { "-s, --rejection-level <sigma>",
	"Rejection level in standard deviation units." },

 { "Kernel specifications (each separated with a semi-colon, ``;''):", NULL },
 { "i/<spatial order>",
	"identity kernel (a.k.a. ``flux term'') with the specified order "
	"of polynomial spatial variations" },
 { "b/<spatial order>",
	"constant offset kernel (a.k.a. ``background term'') with the "
	"specified order of polynomial spatial variations" },
 { "d=<size>/<spatial order>",
	"discrete kernel with the half-size of <size> and the specified "
	"order of polynomial spatial variations" },
 { "g=<size>,<sigma>,<order>/<spatial order>", 
	"Gaussian kernel with the half-size of <size>, standard deviation "
	"of <sigma> and Hermite basis order of <order>, with the specified "
	"order of polynomial spatial variations" },

 { NULL, NULL }

};

int fprint_ficonv_long_help(FILE *fw)
{
 fprintf(fw,
"Usage:\tficonv [options] [-r <reference>] [-i <input>] [-o <output>]\n"
"The purpose of this program is twofold. First, using a reference and input\n"
"FITS image, the program tries to figure out the best-fit convolution kernel\n"
"which transforms the reference image to the input image. The convolution\n"
"kernel solution is saved as a separate, human-readable file. Second, using\n"
"an existing such kernel solution file, `ficonv` convolves the reference\n"
"image and optionally co-adds an additional image to the result.\n\n");
 fprintf(fw,
"The program figures out the desired mode (i.e. whether to fit a convolution\n"
"kernel or use an existing kernel solution to convolve an image) from the\n"
"presence or lack of the command line options.\n\n" );

 longhelp_fprint(fw,ficonv_long_help,0,-1);

 fprintf(fw,"\n");
 fprintf(fw,"Report bugs to <%s>\n",FI_MAINT_EMAIL);

 return(0);
}

/*****************************************************************************/

int main(int argc,char *argv[])
{
 int		i,j,sx,sy,bdc;
 int		is_help,method,niter;
 double		gain,rejectlevel;
 int		psx,psy;

 char		*imgname,*refname,*ikname,*okname,*ixname,*oxname,
 		*oconvname,*osubname,*addname,**inmasklist,*kernelarg,
		*stamparg,*stampfile,*splinestampfile,*xkarg;
 kernellist	klist_data,*klist=&klist_data,
		xlist_data,*xlist=&xlist_data;

 fits		*img,*refimg,*outimg;
 char		**mask_img,**mask_ref,**inmask,**fitmask,**spfmask;

 progbasename=strrchr(argv[0],'/');
 if ( progbasename != NULL )	progbasename++;
 else				progbasename=argv[0];

 is_comment=is_verbose=is_help=0;
 addname=imgname=refname=ikname=okname=ixname=oxname=oconvname=osubname=NULL;
 stamparg=kernelarg=xkarg=NULL; inmasklist=NULL;
 stampfile=splinestampfile=NULL;
 method=0; gain=0.0; niter=0;rejectlevel=3.0;
 bdc=32;
 psx=psy=-1;

 i=scanarg(argc,argv,SCANARG_ALLOW_FLAGS,
	"--version:%NS-1f%q",&is_help,
        "--version-short|--short-version:%NS-2f%q",&is_help,
	"-h|--help:%f%q",&is_help,
	"--long-help|--help-long:%SN2f%q",&is_help,
	"-o|--output|--output-convolved:%s",&oconvname,
	"-g|--gain:%g",&gain,
	"--output-subtracted:%s",&osubname,
	"-i|--input|-f|--fit:%s",&imgname,
	"-M|--input-mask:%t",&inmasklist,
	"-r|--reference:%s",&refname,
	"-a|--add-to:%s",&addname,
	"-t|--stamp|--stamps:%s",&stamparg,
	"--input-stamp|--input-stamps:%s",&stampfile,
	"-p|--pre-spline:%d,%d",&psx,&psy,
	"--input-spline-stamp|--input-spline-stamps:%s",&splinestampfile,
	"--input-kernel-list:%s",&ikname,
	"--output-kernel-list:%s",&okname,
	"-k|--kernel:%s",&kernelarg,
	"--input-extra-kernel-list:%s",&ixname,
	"--output-extra-kernel-list:%s",&oxname,
	"-x|--extra-kernel:%s",&xkarg,
	"-m|--masked:%0f",&method,
	"-w|--weighted:%1f",&method,
	"-b|--background-iterative:%2f",&method,
	"-n|--iterations:%d",&niter,
	"-s|--rejection-level:%g",&rejectlevel,
	"-d|--divide:%d",&bdc,
	"--comment:%f",&is_comment,"(C):%f",&is_comment,
	"--verbose:%i",&is_verbose,"(V):%i",&is_verbose,
	"*:%e",
	NULL);

 if ( i )
  {	fprint_error("invalid command line argument near '%s'",argv[i]);
	return(1);
  }
 else if ( is_help<0 )
  {	fprint_generic_version(stdout,argv[0],"ficonv",FI_CONV_VERSION,is_help);
	return(0);
  }
 else if ( is_help>1 )
  {	fprint_ficonv_long_help(stdout);
	return(0);
  }
 else if ( is_help )
  {	fprint_ficonv_usage(stdout);
	return(0);
  }

 kernel_set_verbosity(is_verbose);

 img=refimg=NULL;sx=sy=0;
 
 if ( imgname != NULL )
  {	FILE	*fr;
	fr=fopenread(imgname);
	if ( fr==NULL )	
	 {	fprint_error("input image: unable to open file '%s'",imgname);
		return(1);
	 }
	img=fits_read(fr);
	fcloseread(fr);
	if ( img==NULL )
	 {	fprint_error("input image: unable to read image contents");
		return(1);
	 }
	if ( img->i.dim != 2 )
	 {	fprint_error("input image: mage is not a 2D one");
		return(1);
	 }
	if ( img->i.sx <= 0 || img->i.sy <= 0 )
	 {	fprint_error("input image: invalid/unexpected size");
		return(1);
	 }
	if ( fits_rescale(img) )
	 {	fprint_error("input image: unable to rescale image");
		return(1);
	 }
	sx=img->i.sx,
	sy=img->i.sy;
  }

 if ( refname != NULL )
  {	FILE	*fr;
 	fr=fopenread(refname);
	if ( fr==NULL )				
	 {	fprint_error("reference image: unable to open file '%s'",imgname);
		return(1);
	 }
	refimg=fits_read(fr);
	fcloseread(fr);
	if ( refimg==NULL )			
	 {	fprint_error("reference image: unable to read image contents");
		return(1);
	 }
	if ( refimg->i.dim != 2 )		
	 {	fprint_error("reference image: mage is not a 2D one");
		return(1);
	 }
	if ( refimg->i.sx <=0 || refimg->i.sy <= 0 )	
	 {	fprint_error("reference image: invalid/unexpected size");
		return(1);
	 }
	if ( fits_rescale(refimg) )		
	 {	fprint_error("reference image: unable to rescale image");
		return(1);
	 }
	if ( img != NULL && ( refimg->i.sx != sx || refimg->i.sy != sy ) )
	 {	fprint_error("reference and input image sizes differ");
		return(1);
	 }
	sx=refimg->i.sx,
	sy=refimg->i.sy;
  }

 if ( refimg != NULL && gain<=0.0 )
	fits_get_gain(refimg,&gain);
 if ( gain<=0.0 )
	gain=1.0;
  
 if ( inmasklist != NULL )
  {	inmask=fits_mask_create_empty(sx,sy);
	if ( join_masks_from_files(inmask,sx,sy,inmasklist) )
	 {	fprint_error("unable to read one of the input mask files");
		return(1);
	 }
  }
 else	inmask=NULL;

 if ( stamparg != NULL )
  {	fitmask=fits_mask_create_empty(sx,sy);
	stamp_parse_argument(stamparg,fitmask,sx,sy);
  }
 else if ( stampfile != NULL )
  {	FILE	*fr;
	fitmask=fits_mask_create_empty(sx,sy);
	fr=fopenread(stampfile);
	if ( fr==NULL )
	 {	fprint_error("stamp file: unable to open file '%s'", stampfile);
		return(1);
	 }
	stamp_read_file(fr,fitmask,sx,sy);
	fcloseread(fr);
  }
 else
	fitmask=NULL;

 if ( splinestampfile != NULL && psx>0 && psy>0 )
  {	FILE	*fr;
	spfmask=fits_mask_create_empty(sx,sy);
	fr=fopenread(splinestampfile);
	if ( fr==NULL )
	 {	fprint_error("spline stamp file: unable to open file '%s'",splinestampfile);
		return(1);
	 }
	stamp_read_file(fr,spfmask,sx,sy);
	fcloseread(fr);
  }
 else
	spfmask=NULL;

 if ( sx>0 && sy>0 )
	logmsg(is_verbose>=1,"[%d,%d]\n",sx,sy);

 if ( kernelarg != NULL && ikname != NULL )
  {	fprint_error("ambiguous kernel definitions");
	return(1);
  }
 if ( xkarg != NULL && ixname != NULL )
  {	fprint_error("ambiguous reference kernel definitions");
	return(1);
  }
 if ( kernelarg == NULL && ikname == NULL )
  {	fprint_error("convolution information is not known/specified");
	return(1);
  }
 if ( xkarg == NULL && ixname == NULL )
	oxname=NULL;	/* ignore -ox if -ix or -x is not specified */
 if ( addname != NULL && imgname != NULL )
  {	fprint_error("invalid combination of input image specifications");
	return(1);
  }

 klist->nkernel=0,
 klist->kernels=NULL;
 if ( kernelarg != NULL )
  {	int	r;
	r=create_kernels_from_kernelarg(kernelarg,klist);
	klist->type=0;
	if ( r )
	 {	fprint_error("unable to parse kernel specification list");
		return(1);
	 }
  }
 if ( ikname != NULL )
  {	FILE	*fr;
	int	r;
	fr=fopenread(ikname);
	if ( fr==NULL )	
	 {	fprint_error("unable to open input kernel list file '%s'",ikname);
		return(1);
	 }
	r=kernel_info_read(fr,klist);
	fcloseread(fr);
	if ( r )
	 {	fprint_error("unable to read or parse contents of input kernel list file '%s'",ikname);
		return(1);
	 }
  }
 for ( i=0 ; i<klist->nkernel ; i++ )
  {	klist->kernels[i].target=0;		}

 xlist->nkernel=0;
 xlist->kernels=NULL;
 if ( xkarg != NULL )
  {	int	r;
	r=create_kernels_from_kernelarg(xkarg,xlist);
	xlist->type=0;
	if ( r )
	 {	fprint_error("unable to parse reference kernel specification list");
		return(1);
	 }
  }
 if ( ixname != NULL )
  {	FILE	*fr;
	int	r;
	fr=fopenread(ixname);
	if ( fr==NULL )
	 {	fprint_error("unable to open reference kernel list file '%s'",ixname);
		return(1);
	 }
	r=kernel_info_read(fr,xlist);
	fcloseread(fr);
	if ( r )
	 {	fprint_error("unable to read or parse contents of reference kernel list file '%s'",ixname);
		return(1);
	 }
  }
 for ( i=0 ; xlist->kernels != NULL && i<xlist->nkernel ; i++ )
  {	xlist->kernels[i].target=1;		}

 logmsg(is_verbose>=1,"Number of kernels: %d+%d.\n",klist->nkernel,xlist->nkernel);
 
 kernel_init_images(klist);
 kernel_init_images(xlist);

 if ( img != NULL )	mask_img=fits_mask_read_from_header(&img->header,sx,sy,NULL);
 else			mask_img=NULL;
 if ( refimg != NULL )	mask_ref=fits_mask_read_from_header(&refimg->header,sx,sy,NULL);
 else			mask_ref=NULL;
	
 /****************************************************************************/
 if ( img != NULL && refimg != NULL )
  {	
	if ( klist->type != 0 )	
	 {	fprint_error("kernel specification contains fitted amplitudes");
		return(1);
	 }

	if ( psx>0 && psy>0 )
	 {	double	*coeff;

		if ( spfmask==NULL )
			spfmask=fits_mask_create_empty(sx,sy);
		if ( inmask != NULL )
			fits_mask_and(spfmask,sx,sy,inmask);
		fits_mask_and(spfmask,sx,sy,mask_ref);
		fits_mask_and(spfmask,sx,sy,mask_img);

		coeff=(double *)malloc(sizeof(double)*2*(psx+1)*(psy+1));

		affine_img_transformation_spline_fit(&img->i,&refimg->i,spfmask,psx,psy,coeff);
		affine_img_transformation_spline_eval(&img->i,mask_img,psx,psy,coeff);

		free(coeff);
	 }

	if ( fitmask==NULL )
		fitmask=fits_mask_create_empty(sx,sy);
	if ( inmask != NULL )
		fits_mask_and(fitmask,sx,sy,inmask);

	fit_kernels(&refimg->i,mask_ref,&img->i,mask_img,fitmask,klist,xlist,method,bdc,niter,rejectlevel,gain);
  }
 if ( fitmask != NULL )
  {	fits_mask_free(fitmask);
	fitmask=NULL;
  }
 if ( spfmask != NULL )
  {	fits_mask_free(spfmask);
	spfmask=NULL;
  }

 if ( refimg != NULL && ( oconvname != NULL || osubname != NULL ) )
  {	
	int	hsize;
	char	**mask,**mask_conv;

	logmsg(is_verbose>=1,"Convolving ...\n");
	if ( ! klist->type )
	 {	fprint_error("kernel specification lacks amplitudes");
		return(1);
	 }

	hsize=0;
	for ( i=0 ; i<klist->nkernel ; i++ )
	 {	if ( klist->kernels[i].type==KERNEL_BACKGROUND || klist->kernels[i].type==KERNEL_IDENTITY )
			continue;
		if ( klist->kernels[i].hsize > hsize )
			hsize=klist->kernels[i].hsize;
	 }
	
	mask=fits_mask_create_empty(sx,sy);
	fits_mask_and(mask,sx,sy,mask_ref);
	if ( inmask != NULL )	fits_mask_and(mask,sx,sy,inmask);
	mask_conv=fits_mask_expand_false(mask,sx,sy,hsize,-1,-1,1);

	outimg=fits_duplicate_empty(refimg);
	convolve_with_kernel_set(&refimg->i,mask_conv,klist,&outimg->i);

	if ( addname != NULL )
	 {	FILE	*fr;
		fits	*addimg;
		char	**addmask;
		while ( 1 )
		 {	fr=fopenread(addname);
			if ( fr == NULL )
				break;
			addimg=fits_read(fr);
			fcloseread(fr);
			if ( addimg==NULL )
				break;
			if ( addimg->i.data==NULL || addimg->i.sx != sx || addimg->i.sy != sy || addimg->i.dim != 2 )
				break;
			if ( fits_rescale(addimg) )
				break;
			addmask=fits_mask_read_from_header(&addimg->header,sx,sy,NULL);
			for ( i=0 ; i<sy ; i++ )
			 {  for ( j=0 ; j<sx ; j++ )
			     {	if ( ! addmask[i][j] && ! mask_conv[i][j] )
					outimg->i.data[i][j]+=addimg->i.data[i][j],
					mask[i][j]=0;
				else
					mask[i][j]=-1;
			     }
			 }
			fits_mask_free(addmask);
			break;
		 };
	 }
	else
	 {	for ( i=0 ; i<sy ; i++ )
		 {	for ( j=0 ; j<sx ; j++ )
			 {	mask[i][j]=mask_conv[i][j];		}
		 }
	 }

        fits_history_export_command_line(img,"ficonv",FI_CONV_VERSION,argc,argv);

	if ( oconvname != NULL )
	 {	FILE	*fw;
		fw=fopenwrite(oconvname);
		if ( fw==NULL )	
		 {	fprint_error("unable to create output image '%s'",oconvname);
			return(1);
		 }
		mark_integerlimited_pixels(&outimg->i,mask,outimg->i.bit,1,MASK_OVERSATURATED,MASK_OVERSATURATED);
		fits_mask_export_as_header(&outimg->header,1,mask,sx,sy,NULL);
		fits_write(fw,outimg);
		fclosewrite(fw);
	 }
					
	if ( osubname != NULL && img != NULL )
	 {	FILE	*fw;

		make_subtracted_image(&outimg->i,mask,&img->i,mask_img,&outimg->i,mask,xlist);

		fw=fopenwrite(osubname);
		if ( fw==NULL )	
		 {	fprint_error("unable to create output subtracted image '%s'",osubname);
			return(1);
		 }
		fits_copy_full_header(outimg,img);
		fits_mask_export_as_header(&outimg->header,1,mask,sx,sy,NULL);
		fits_write(fw,outimg);
		fclosewrite(fw);
	 }
	fits_mask_free(mask_conv);
	fits_mask_free(mask);
  }

 if ( okname != NULL )
  {	FILE	*fw;
	fw=fopenwrite(okname);
	if ( fw==NULL )	
	 {	fprint_error("unable to create output kernel list file '%s'",okname);
		return(1);
	 }
	kernel_info_write(fw,klist);
	fclosewrite(fw);
  }
 if ( oxname != NULL )
  {	FILE	*fw;
	fw=fopenwrite(oxname);
	if ( fw==NULL )	
	 {	fprint_error("unable to create output reference kernel list file '%s'",oxname);
		return(1);
	 }
	kernel_info_write(fw,xlist);
	fclosewrite(fw);
  }

 if ( mask_ref != NULL )	fits_mask_free(mask_ref);
 if ( mask_img != NULL )	fits_mask_free(mask_img);
 if ( inmask != NULL )		fits_mask_free(inmask);

 return(0);
}

