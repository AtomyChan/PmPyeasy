/*****************************************************************************/
/* fiarith.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line user interface for basic arithmetics on FITS images.	     */
/*****************************************************************************/
#define	FI_ARITH_VERSION	"1.0pre3"
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include <psn/psn.h>
#include <fits/fits.h>

#include "fi.h"

#include "math/spline/biquad.h"
#include "math/dft/pbfft.h"
#include "fitsmask.h"
#include "psn/psn-general.h"
#include "io/iof.h"
#include "io/scanarg.h"

#include "history.h"
#include "longhelp.h"

#include "imgtrans.h"
#include "common.h"
#include "tensor.h"

/*****************************************************************************/

typedef struct
 {	char	filename[240];
	char	symname[12];
	fits	*img;
	int	is_read;
 } operand;

typedef struct
 {	int	sx,sy;		/* size of the (target) image(s) */
	operand *ops;		/* set of image operands	 */
	int	nop;		/* number of image operands	 */
	psn	**udfs;		/* set of user-defined functions */
	int	nudf;		/* number of user-def. functions */
	fits	*hdrsave;	/* output will inherit this hdr! */
	char	**mask;
 } evaldata;

/*****************************************************************************/

#define		IO_ADD		1	/* + */
#define		IO_SUB		2	/* - */
#define		IO_MUL		3	/* * */
#define		IO_DIV		4	/* / */
#define		IO_CHS		6	/* change sign, unary variant of '-' */

#define		IF_MIN		32
#define		IF_MAX		33
#define		IF_MEA		34

#define		IF_ABS		36
#define		IF_SQRT		37
#define		IF_SQ		38
#define		IF_NORM		39
#define		IF_SIGN		40
#define		IF_THETA	41

#define		IF_LAP		56
#define		IF_SCT		57
#define		IF_SMT		58
#define		IF_WGH		59
#define		IF_CORR		60

#define		UDF_OFFSET	64

psnsym		psn_subimg_var[] = {	{ T_VAR, 0, "x", 0 },
					{ T_VAR, 1, "y", 0 },
					{ T_VAR, 2, "X", 0 },
					{ T_VAR, 3, "Y", 0 },
					{ T_VAR, 4, "a", 0 },
					{ T_VAR, 5, "b", 0 },
					{ T_VAR, 6, "c", 0 },
					{ T_VAR, 7, "d", 0 },
					{ T_VAR, 8, "e", 0 },
					{ T_VAR, 9, "f", 0 },
					{ T_VAR,10, "g", 0 },
					{ T_VAR,11, "h", 0 },
					{ 0,0,NULL,0 }
				   };

psnsym		psn_img_op[]	= {	{ T_OP, IO_ADD, "+", TO_INFIX  },
					{ T_OP, 0     , "+", TO_PREFIX },
					{ T_OP, IO_SUB, "-", TO_INFIX  },
					{ T_OP, IO_CHS, "-", TO_PREFIX },
					{ T_OP, IO_MUL, "*", TO_INFIX  },
					{ T_OP, IO_DIV, "/", TO_INFIX  },

					{0,0,NULL,0}
				  };

psnsym		psn_img_fn[]	= {	
					{ T_FN, IF_MIN, "min"    , -1  },
					{ T_FN, IF_MAX, "max"    , -1  },
					{ T_FN, IF_MEA, "mean"   , -1  },
					{ T_FN, IF_ABS, "abs"	 ,  1  },
					{ T_FN, IF_SIGN,"sign"	 ,  1  },
					{ T_FN, IF_THETA,"theta" ,  1  },
					{ T_FN, IF_NORM, "norm"	 ,  1  },
					{ T_FN, IF_SQRT,"sqrt"	 ,  1  },
					{ T_FN, IF_SQ,  "sq"	 ,  1  },
					{ T_FN, IF_LAP, "laplace",  1  },
					{ T_FN, IF_SCT, "scatter",  1  },
					{ T_FN, IF_SMT, "smooth" ,  3  },
					{ T_FN, IF_WGH, "weight" ,  3  },
					{ T_FN, IF_CORR,"corr"	 ,  2  },
					{ 0,0,NULL,0 } 
				 };


psnprop		psn_img_prop[] = {	{ IO_ADD,  2, 20, ASSOC_LEFT },
					{ IO_SUB,  2, 20, ASSOC_LEFT },
					{ IO_MUL,  2, 21, ASSOC_LEFT },
					{ IO_DIV,  2, 21, ASSOC_LEFT },
					{ IO_CHS,  1, 22 },
					{ IF_MIN, -1, 0 },
					{ IF_MAX, -1, 0 },
					{ IF_MEA, -1, 0 },	
					{ IF_ABS,  1, 0 },
					{ IF_SIGN, 1, 0 },
					{ IF_THETA,1, 0 },
					{ IF_NORM, 1, 0 },
					{ IF_SQRT, 1, 0 },
					{ IF_SQ,   1, 0 },
					{ IF_LAP,  1, 0 },
					{ IF_SCT,  1, 0 },
					{ IF_SMT,  3, 0 },
					{ IF_WGH,  3, 0 },
					{ IF_CORR, 2, 0 },
					{ 0, 0, 0,0  } };

/*****************************************************************************/

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

typedef struct
 {	fits	*img;
	double	val;
 } imgstack;

/*****************************************************************************/

fits *newimage_empty(evaldata *ed)
{
 fits	*img;

 img=fits_create();
 fits_alloc_image(img,ed->sx,ed->sy);

 fits_set_header_boolean(img,"SIMPLE",FITS_SH_FIRST,1,"Fits standard");
 fits_set_header_integer(img,"BITPIX",FITS_SH_FIRST,-32,NULL);
 fits_set_header_integer(img,"NAXIS",FITS_SH_FIRST,2,NULL);
 fits_set_header_integer(img,"NAXIS1",FITS_SH_FIRST,ed->sx,NULL);
 fits_set_header_integer(img,"NAXIS2",FITS_SH_FIRST,ed->sy,NULL);

 img->i.bit=-32;
 img->i.curr.bscale=1.0;
 img->i.curr.bzero=0.0;
 img->i.read.bscale=1.0;
 img->i.read.bzero=0.0;

 return(img);
}
fits *newimage(evaldata *ed,double d)
{
 fits	*img;
 int	i,j;

 img=newimage_empty(ed);

 for ( i=0 ; i<img->i.sy ; i++ )
  {	for ( j=0 ; j<img->i.sx ; j++ )
	 {	img->i.data[i][j]=d;		}
  }	

 return(img);
}

/*****************************************************************************/

int scatter(fits *img)
{ 
 double	**bqc;
 int	i,j,sx,sy;

 sx=img->i.sx,
 sy=img->i.sy;

 bqc=(double **)tensor_alloc_2d(double,2*sx+1,2*sy+1);
 biquad_coeff(img->i.data,sx,sy,bqc,NULL);
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	img->i.data[i][j]=biquad_scatter(bqc,j,i);		}
  }
 tensor_free(bqc);
 return(0);
}

/*****************************************************************************/

int evaluate_min(imgstack *stack,int narg,int sx,int sy)
{
 int	i,j,k;
 double	w,c;
 fits	*img;
 if ( narg <= 0 )
  {	stack[0].img=NULL;stack[0].val=0.0;return(0);	  }
 for ( k=0,img=NULL ; k<narg && img==NULL ; k++ )
  {	img=stack[k].img;	}
 if ( img==NULL )
  {	w=stack[0].val;
	for ( k=1 ; k<narg ; k++ )
	 {	if ( stack[k].val<w )	w=stack[k].val;		}
	stack[0].val=w;
	return(0);
  }
 for ( i=0,w=0.0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	for ( k=0,w=0.0 ; k<narg ; k++ )
		 {	if ( stack[k].img != NULL )
				c=stack[k].img->i.data[i][j];
			else
				c=stack[k].val;

			if ( k==0 )	w=c;
			else if ( c<w )	w=c;
		 }
		img->i.data[i][j]=w;
	 } 
  }
 for ( k=0 ; k<narg ; k++ )
  {	if ( stack[k].img != NULL && stack[k].img != img )
		fits_free(stack[k].img);
  }
 stack[0].img=img;stack[0].val=0.0;
 return(0);
}

int evaluate_max(imgstack *stack,int narg,int sx,int sy)
{
 int	i,j,k;
 double	w,c;
 fits	*img;
 if ( narg <= 0 )
  {	stack[0].img=NULL;stack[0].val=0.0;return(0);	  }
 for ( k=0,img=NULL ; k<narg && img==NULL ; k++ )
  {	img=stack[k].img;	}
 if ( img==NULL )
  {	w=stack[0].val;
	for ( k=1 ; k<narg ; k++ )
	 {	if ( stack[k].val>w )	w=stack[k].val;		}
	stack[0].val=w;
	return(0);
  }
 for ( i=0,w=0.0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	for ( k=0,w=0.0 ; k<narg ; k++ )
		 {	if ( stack[k].img != NULL )
				c=stack[k].img->i.data[i][j];
			else
				c=stack[k].val;
			if ( k==0 )	w=c;
			else if ( c>w )	w=c;
		 }
		img->i.data[i][j]=w;
	 } 
  }
 for ( k=0 ; k<narg ; k++ )
 {	if ( stack[k].img != NULL && stack[k].img != img )
		fits_free(stack[k].img);
  }
 stack[0].img=img;stack[0].val=0.0;
 return(0);
}

int evaluate_mean(imgstack *stack,int narg,int sx,int sy)
{
 int	i,j,k;
 double	w,c,div;
 fits	*img;
 if ( narg <= 0 )
  {	stack[0].img=NULL;stack[0].val=0.0;return(0);	  }
 div=1.0/(double)narg;
 for ( k=0,img=NULL ; k<narg && img==NULL ; k++ )
  {	img=stack[k].img;	}
 if ( img==NULL )
  {	w=0.0;
	for ( k=0 ; k<narg ; k++ )
	 {	w+=stack[k].val;		}
	stack[0].val=w*div;
	return(0);
  }
 for ( i=0,w=0.0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	for ( k=0,w=0.0 ; k<narg ; k++ )
		 {	if ( stack[k].img != NULL )
				c=stack[k].img->i.data[i][j];
			else
				c=stack[k].val;
			w+=c;
		 }
		img->i.data[i][j]=w*div;
	 } 
  }
 for ( k=0 ; k<narg ; k++ )
  {	if ( stack[k].img != NULL && stack[k].img != img )
		fits_free(stack[k].img);
  }
 stack[0].img=img;stack[0].val=0.0;
 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int evaluate_smooth(evaldata *ed,imgstack *stack,double sigma,double dhsize)
{
 int	hsize,fsize,i,j,sx,sy,k,l;
 double	**kernel,sum,w,km,x0,x1,y0,y1;
 char	**nmask;
 fits	*img,*nimg;
 if ( stack->img==NULL )	return(0);	/* nothing to be done */
 hsize=(int)dhsize;
 if ( hsize<=0 )		return(0);

 fsize=2*hsize+1;
 kernel=tensor_alloc_2d(double,fsize,fsize);

 km=1.0/(sqrt(2.0)*sigma);
 for ( i=-hsize ; i<=hsize ; i++ )
  {	y0=km*((double)i-0.5);
	y1=km*((double)i+0.5);
	for ( j=-hsize ; j<=hsize ; j++ )
	 {	x0=km*((double)j-0.5);
		x1=km*((double)j+0.5);
		kernel[i+hsize][j+hsize]=(erf(x1)-erf(x0))*(erf(y1)-erf(y0));
	 }
  }

 sum=0.0;
 for ( i=0 ; i<fsize ; i++ )
  {	for ( j=0 ; j<fsize ; j++ )
	 {	sum+=kernel[i][j];		}
  }
 for ( i=0 ; i<fsize ; i++ )
  {	for ( j=0 ; j<fsize ; j++ )
	 {	kernel[i][j]/=sum;		}
  }

 img=stack->img;
 sx=ed->sx,
 sy=ed->sy;
 nmask=fits_mask_expand_false(ed->mask,sx,sy,hsize,-1,-1,1);
 nimg=fits_duplicate(img);
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( nmask[i][j] )
		 {	if ( ! ed->mask[i][j] )	nimg->i.data[i][j]=0.0;
			continue;
		 }
		w=0.0;
		for ( k=-hsize ; k<=hsize ; k++ )
		 {	for ( l=-hsize ; l<=hsize ; l++ )
			 {	w+=img->i.data[i-k][j-l]*
					kernel[hsize+k][hsize+l];
			 }
		 }
		nimg->i.data[i][j]=w;
	 }
  }
 fits_mask_free(ed->mask);
 ed->mask=nmask;
 fits_free(img);
 stack->img=nimg;

 tensor_free(kernel);
 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int evaluate_weight(evaldata *ed,imgstack *simg,imgstack *swgh,double dbsize)
{
 int	bsize,sx,sy,i,j,bx,by,k,l;
 fits	*img,*wgh;
 double	sm,ii,iw,sw;

 /* nothing to be done if any of the images (target image and/or weight)   */
 /* are constant, so (imgstack *)->img == NULL (no associated FITS image): */
 if ( simg->img==NULL || swgh->img==NULL )	return(0);

 sx=ed->sx,
 sy=ed->sy;
 img=simg->img;
 wgh=swgh->img;
 
 bsize=(int)dbsize;

 /* if block size is less than or equal to 1, do nothing. */
 if ( bsize<=1 )	return(0);

 for ( i=0 ; i<sy ; i+=bsize )
  {	by=bsize;
	if ( i+by>sy )	by=sy-i;
	for ( j=0 ; j<sx ; j+=bsize )
	 {	bx=bsize;
		if ( j+bx>sx )	bx=sx-j;
		sm=sw=0.0;
		for ( k=0 ; k<by ; k++ )
		 {	for ( l=0 ; l<bx ; l++ )
			 {	ii=img->i.data[i+k][j+l];
				iw=wgh->i.data[i+k][j+l];
				sm+=(ii-iw);
				sw+=iw;
			 }
		 }
		if ( sw<=0.0 )	continue;
		sm=sm/(double)(bx*by);
		for ( k=0 ; k<by ; k++ )
		 {	for ( l=0 ; l<bx ; l++ )
			 {	iw=wgh->i.data[i+k][j+l];
				img->i.data[i+k][j+l] = iw+sm;
			 }
		 }
	 }
  }
	
 return(0); 
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int evaluate_correlation(evaldata *ed,imgstack *s1,imgstack *s2)
{
 int		sx,sy,i,j;
 double		val;

 if ( s1->img==NULL && s2->img==NULL )
  {	s1->val = s1->val * s2->val;
	return(0);
  }
 else if ( s1->img==NULL || s2->img==NULL )
  {	if ( s1->img != NULL )
	 {	s2->img=s1->img;
		s1->val=s2->val;
		s1->img=NULL;
	 }
	sx=ed->sx,
	sy=ed->sy;
	val=0.0;
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	val+=s2->img->i.data[i][j];		}
	 }
	s1->val = s1->val * val;
	return(0);
  }
 else
  {	complex		**c1,**c2,*c,w;
	fitsimage	*i1,*i2;
	double		norm;

	i1=&s1->img->i;
	i2=&s2->img->i;

	sx=ed->sx;
	sy=ed->sy;

	c1=tensor_alloc_2d(complex,sx,sy);
	c2=tensor_alloc_2d(complex,sx,sy);
	c =tensor_alloc_1d(complex,sy);

	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
	 	 {	c1[i][j].re=i1->data[i][j];
			c1[i][j].im=0.0;
			c2[i][j].re=i2->data[i][j];
			c2[i][j].im=0.0;
		 }
		pbfft_conv(c1[i],sx,0);
		pbfft_conv(c2[i],sx,0);
	 }

	for ( j=0 ; j<sx ; j++ )
	 {	for ( i=0 ; i<sy ; i++ )
		 {	c[i].re=c1[i][j].re;
			c[i].im=c1[i][j].im;
	  	 }
		pbfft_conv(c,sy,0);
		for ( i=0 ; i<sy ; i++ )
		 {	c1[i][j].re=c[i].re;
			c1[i][j].im=c[i].im;
	  	 }

		for ( i=0 ; i<sy ; i++ )
		 {	c[i].re=c2[i][j].re;
			c[i].im=c2[i][j].im;
	  	 }
		pbfft_conv(c,sy,0);
		for ( i=0 ; i<sy ; i++ )
		 {	c2[i][j].re=c[i].re;
			c2[i][j].im=c[i].im;
	  	 }
	 }
	
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	w.re=+c1[i][j].re*c2[i][j].re+c1[i][j].im*c2[i][j].im;
			w.im=+c1[i][j].re*c2[i][j].im-c1[i][j].im*c2[i][j].re;
			c1[i][j].re=w.re;
			c1[i][j].im=w.im;
		 }
	 }
	/*
	c1[0][0].re=0.0;
	c1[0][0].im=0.0;
	*/

	for ( j=0 ; j<sx ; j++ )
	 {	for ( i=0 ; i<sy ; i++ )
		 {	c[i].re=c1[i][j].re;
			c[i].im=c1[i][j].im;
	  	 }
		pbfft_conv(c,sy,1);
		for ( i=0 ; i<sy ; i++ )
		 {	c1[i][j].re=c[i].re;
			c1[i][j].im=c[i].im;
	  	 }
	 }

	for ( i=0 ; i<sy ; i++ )
	 {	pbfft_conv(c1[i],sx,1);		}

	norm=1.0/((double)sx*(double)sy);

	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	w.re=c1[i][j].re;
			w.im=c1[i][j].im;
			i1->data[i][j]=norm*w.re;
		 }
	 }
	
	tensor_free(c);
	tensor_free(c2);
	tensor_free(c1);

	return(0);
  }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

fits *evaluate(evaldata *ed,psn *pseq,int hid)
{
 FILE	  *fr;
 imgstack stack[32]; 
 fits	  *img;
 psnterm  *seq;
 int	  i,j,k,sp,sx,sy,udfid,narg,iac;
 psn	  *udf;
 double	  uvars[16],res,crc,w;

 sx=ed->sx,
 sy=ed->sy;

 sp=0;ed->hdrsave=NULL;
 for ( seq=pseq->terms ; seq->type ; seq++ )
  { narg=seq->minor;
    switch( seq->type )
    {  
     case T_CONST:
	stack[sp].img=NULL;
	stack[sp].val=pseq->cons[seq->major];
	sp++;
  	break;
     case T_SCONST:
	stack[sp].img=NULL;
	stack[sp].val=(double)(seq->major);
	sp++;
	break;
     case T_VAR:
	k=seq->major;
	if ( ed->ops[k].is_read )
	 {	img=ed->ops[k].img;		}
	else
	 {	fr=fopen(ed->ops[k].filename,"rb");
		img=fits_read(fr);
		fits_rescale(img);
		fclose(fr);
	 }
	if ( k==hid )
	 {	ed->hdrsave=fits_create();
		fits_copy_full_header(ed->hdrsave,img);
	 }
	stack[sp].img=img;
	stack[sp].val=0.0;
	sp++;
	break;
     case T_OP: case T_FN:
	switch ( seq->major )
	 {   case IO_ADD:
		if ( stack[sp-2].img!=NULL && stack[sp-1].img!=NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	stack[sp-2].img->i.data[i][j]+=
					stack[sp-1].img->i.data[i][j];
				 }
			 }
			fits_free(stack[sp-1].img);
		 }
		else if ( stack[sp-2].img==NULL && stack[sp-1].img!=NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	stack[sp-1].img->i.data[i][j]+=stack[sp-2].val;	}
			 }
			stack[sp-2].img=stack[sp-1].img;
		 }
		else if ( stack[sp-2].img!=NULL && stack[sp-1].img==NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	stack[sp-2].img->i.data[i][j]+=stack[sp-1].val;	}
			 }
		 }
		else
		 {	stack[sp-2].val += stack[sp-1].val;	}
		sp--;
		break;
	     case IO_SUB:
		if ( stack[sp-2].img!=NULL && stack[sp-1].img!=NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	stack[sp-2].img->i.data[i][j]-=
					stack[sp-1].img->i.data[i][j];
				 }
			 }
			fits_free(stack[sp-1].img);
		 }
		else if ( stack[sp-2].img==NULL && stack[sp-1].img!=NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	stack[sp-1].img->i.data[i][j]=stack[sp-2].val-stack[sp-1].img->i.data[i][j];	}
			 }
			stack[sp-2].img=stack[sp-1].img;
		 }
		else if ( stack[sp-2].img!=NULL && stack[sp-1].img==NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	stack[sp-2].img->i.data[i][j]-=stack[sp-1].val;	}
			 }
		 }
		else
		 {	stack[sp-2].val -= stack[sp-1].val;	}
		sp--;
		break;
	     case IO_MUL:
		if ( stack[sp-2].img!=NULL && stack[sp-1].img!=NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	stack[sp-2].img->i.data[i][j]*=
					stack[sp-1].img->i.data[i][j];
				 }
			 }
			fits_free(stack[sp-1].img);
		 }
		else if ( stack[sp-2].img==NULL && stack[sp-1].img!=NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	stack[sp-1].img->i.data[i][j]*=stack[sp-2].val;	}
			 }
			stack[sp-2].img=stack[sp-1].img;
		 }
		else if ( stack[sp-2].img!=NULL && stack[sp-1].img==NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	stack[sp-2].img->i.data[i][j]*=stack[sp-1].val;	}
			 }
		 }
		else
		 {	stack[sp-2].val *= stack[sp-1].val;	}
		sp--;
		break;
	     case IO_DIV:
		if ( stack[sp-2].img!=NULL && stack[sp-1].img!=NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	w=stack[sp-1].img->i.data[i][j];
					if ( w != 0.0 )
						stack[sp-2].img->i.data[i][j]/=w;
					else
						stack[sp-2].img->i.data[i][j]=0.0;
				 }
			 }
			fits_free(stack[sp-1].img);
		 }
		else if ( stack[sp-2].img==NULL && stack[sp-1].img!=NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	w=stack[sp-1].img->i.data[i][j];
					if ( w != 0.0 )
						stack[sp-1].img->i.data[i][j]=stack[sp-2].val/w;
					else
						stack[sp-1].img->i.data[i][j]=0.0;
				 }
			 }
			stack[sp-2].img=stack[sp-1].img;
		 }
		else if ( stack[sp-2].img!=NULL && stack[sp-1].img==NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	stack[sp-2].img->i.data[i][j]/=stack[sp-1].val;	}
			 }
		 }
		else
		 {	stack[sp-2].val /= stack[sp-1].val;	}
		sp--;
		break;
	     case IO_CHS:
		if ( stack[sp-1].img != NULL )
	 	 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
					stack[sp-1].img->i.data[i][j]=-stack[sp-1].img->i.data[i][j];
			 }
		 }
		else
		 {	stack[sp-1].val=-stack[sp-1].val;	}
		break;
	     case IF_MIN:
		evaluate_min(&stack[sp-narg],narg,sx,sy);
		sp+=1-narg;
		break;
	     case IF_MAX:
		evaluate_max(&stack[sp-narg],narg,sx,sy);
		sp+=1-narg;
		break;
	     case IF_MEA:
		evaluate_mean(&stack[sp-narg],narg,sx,sy);
		break;
	     case IF_ABS:
		if ( stack[sp-1].img != NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
					stack[sp-1].img->i.data[i][j]=fabs(stack[sp-1].img->i.data[i][j]);
			 }
		 }
		else
			stack[sp-1].val=fabs(stack[sp-1].val);
		break;
	     case IF_SIGN:
		if ( stack[sp-1].img != NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	double	d;
					d=stack[sp-1].img->i.data[i][j];
					stack[sp-1].img->i.data[i][j]=(d<0?-1.0:d>0?+1.0:0.0);
				 }
			 }
		 }
		else
		 {	double	d;
			d=stack[sp-1].val;
			stack[sp-1].val=(d<0?-1.0:d>0?+1.0:0.0);
		 }
		break;
	     case IF_THETA:
		if ( stack[sp-1].img != NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	double	d;
					d=stack[sp-1].img->i.data[i][j];
					stack[sp-1].img->i.data[i][j]=(d<0?-1.0:+1.0);
				 }
			 }
		 }
		else
		 {	double	d;
			d=stack[sp-1].val;
			stack[sp-1].val=(d<0?-1.0:+1.0);
		 }
		break;
	     case IF_NORM:
		if ( stack[sp-1].img != NULL )
		 {	double	norm,n;
			norm=0.0,n=0.0;
			for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	norm += abs(stack[sp-1].img->i.data[i][j]);
					n    += 1.0;
				 }
			 }
			fits_free(stack[sp-1].img);
			stack[sp-1].img=NULL;
			stack[sp-1].val=norm/n;
		 }
		else
			stack[sp-1].val=fabs(stack[sp-1].val);
		break;
	     case IF_SQRT:
		if ( stack[sp-1].img != NULL )
		 {	double	d;
			for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	d=stack[sp-1].img->i.data[i][j];
					if ( d>0.0 )	d=sqrt(d);
					else		d=0.0;
					stack[sp-1].img->i.data[i][j]=d;
				 }
			 }
		 }
		else
		 {	double	d;
			d=stack[sp-1].val;
			if ( d>0.0 )	d=sqrt(d);
			else		d=0.0;
			stack[sp-1].val=d;
		 }
		break;
	     case IF_SQ:
		if ( stack[sp-1].img != NULL )
		 {	double	d;
			for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	d=stack[sp-1].img->i.data[i][j];
					stack[sp-1].img->i.data[i][j]=d*d;
				 }
			 }
		 }
		else
		 {	double	d;
			d=stack[sp-1].val;
			stack[sp-1].val=d*d;
		 }
		break;
	     case IF_LAP:
		if ( stack[sp-1].img != NULL )
			cyclic_laplace_of_image(&stack[sp-1].img->i);
		else
			stack[sp-1].val=0.0;
		break;
	     case IF_SCT:
		if ( stack[sp-1].img != NULL )
			scatter(stack[sp-1].img);
		else
			stack[sp-1].val=0.0;
		break;
	     case IF_SMT:
		evaluate_smooth(ed,&stack[sp-3],stack[sp-2].val,stack[sp-1].val);
		if ( stack[sp-2].img != NULL )	fits_free(stack[sp-2].img);
		if ( stack[sp-1].img != NULL )	fits_free(stack[sp-1].img);
		sp+=1-3;
		break;
	     case IF_WGH:
		evaluate_weight(ed,&stack[sp-3],&stack[sp-2],stack[sp-1].val);
		if ( stack[sp-2].img != NULL )	fits_free(stack[sp-2].img);
		if ( stack[sp-1].img != NULL )	fits_free(stack[sp-1].img);
		sp+=1-3;
		break;
	     case IF_CORR:
		evaluate_correlation(ed,&stack[sp-2],&stack[sp-1]);
		if ( stack[sp-1].img != NULL )	fits_free(stack[sp-1].img);
		sp+=1-2;
		break;
	     default :
		if ( seq->major < UDF_OFFSET )	break;
		udfid=seq->major-UDF_OFFSET;			
		udf=ed->udfs[udfid];
		narg=seq->minor;
		img=newimage(ed,0.0);
		iac=-1;crc=0.0;
		for ( i=0 ; i<sy ; i++ )
		 {	for ( j=0 ; j<sx ; j++ )
			 {	uvars[0]=(double)(2*j-sx)/(double)sx;
				uvars[1]=(double)(2*i-sy)/(double)sx;
				uvars[2]=(double)j;
				uvars[3]=(double)i;
				for ( k=0 ; k<narg ; k++ )
				 {	if ( stack[sp-narg+k].img != NULL )
						uvars[4+k]=stack[sp-narg+k].img->i.data[i][j];
					else
						uvars[4+k]=stack[sp-narg+k].val;
				 }
				k=psn_double_calc(udf,psn_general_funct,&res,uvars);
				if ( k )	res=0.0;
				if ( iac<0 )	crc=res,iac=0;
				else if ( ! iac && crc != res )	iac=1;
				img->i.data[i][j]=res;
			 }
		 }
		for ( k=0 ; k<narg ; k++ )
		 {	if ( stack[sp-narg+k].img != NULL )
				fits_free(stack[sp-narg+k].img);
		 }
		sp-=narg;
		if ( ! iac )
		 {	fits_free(img);
			stack[sp].img=NULL;
			stack[sp].val=crc;
		 }
		else
		 {	stack[sp].img=img;
			stack[sp].val=0.0;
		 }
		sp++;
	 }
	break;
     }
  }
 if ( sp != 1 )				return(NULL);
 else if ( stack[0].img != NULL )	return(stack[0].img);
 else					return(newimage(ed,stack[0].val));
}

/*****************************************************************************/

int fprint_fiarith_usage(FILE *fw)
{
 fprintf(fw,
"Usage:\tfiarith [-h|--help|--long-help] [--version]\n"
"\t\"<expression>\" [-o|--output <out.fits>] \n"
"\t[-b|--bitpix <output-bitpix>] [-e|--header <inherit-header-from.fits>]\n"
"\t[-D|--data bitpix=<bitpix>,bscale=<scale>,bzero=<zero>|<C-type>]\n"
"\t[-M|--input-mask <mask> [-im <...>]] [--output-mask <mask>]\n"
"\t[-s|--size <sx,sy>] [-a|--apply-mask]\n");
 fprintf(fw,
"Available image operators:\n"
"\t+, -, *, /.\n");
 fprintf(fw,
"Available image functions:\n"
"\tmin(...), max(...), mean(...), abs(.), sqrt(.), sq(.),\n"
"\tsign(.), theta(.),\n"
"\tlaplace(.), scatter(.), corr(.,.), \n"
"\tsmooth(.,<sigma>,<size>), weight(.,.,<blocksize>).\n");
 fprintf(fw,
"Available general operators and functions: see ./src/psn/psn-general.[ch]\n");
 return(0);
}

longhelp_entry fiarith_long_help[]=
{
 LONGHELP_OPTIONS, 
 { "General options:", NULL },
 { "-h, --help",
	"Gives general summary about the command line options." },
 { "--long-help",
	"Gives a detailed list of command line options." },
 { "--version",
	"Gives some version information about the program." },
 { "-o, --output <fits>",
	"The name of the output file (omitting or specifing '-' yields the "
	"output to be written to stdout)." },
 { "-s, --size <sx>,<sy>",
	"Horizontal and vertical size of the output image (meaningless "
	"if <expression> contains any existing FITS image)." },
 { "-M, --input-mask <fits>",
	"Input mask file to co-add to output image." },
 { "--output-mask <fits>",
	"Name of output mask file to create." },
 { "-D, --data <spec>",
	"Output pixel data format specification." },
 { "-b, --bitpix <bitpix>",
	"Standard FITS output bitpix value." },
 { "-a, --apply-mask",
	"Resets pixel values to zero if masked." },
 { "-e, --header <fits>",
	"Inherit output header from the specified FITS file." },

 { "Available arithmetic operators:", NULL },
 { "+",
	"Addition." },
 { "-",
	"Subtraction." },
 { "*",
	"Multiplication." },
 { "/",
	"Division. Dividing by zero yields zero automatically." },

 { "Available arithmetic functions:", NULL },
 { "sin(.), cos(.), tan(.)",
	"Trigonometric functions." },
 { "asin(.), acos(.), atan(.), atan2(.,.), arg(.,.)",
	"Inverse trigonometric functions." },
 { "abs(.), sq(.), sqrt(.), exp(.), log(.)",
	"Absolute value, square, Square root, exponential and natural logarithm functions." },
 
 { "Available image functions:", NULL },
 { "min(....)",
	"The per pixel minimal value of pixel intensities of the images "
	"listed (separated by comma) in the argument of the function." },
 { "max(....)",
	"The per pixel maximal value of pixel intensities of the images "
	"listed (separated by comma) in the argument of the function." },
 { "mean(...)",
	"The per pixel mean value of pixel intensities of the images "
	"listed (separated by comma) in the argument of the function." },
 { "norm(.)", 
	"The mean of the pixel values on the image." },
 { "sign(.)",
	"The per pixel sign function." },
 { "theta(.)",
	"The per pixel Heaviside step function." },
 { "laplace(.)",
	"The Laplace transform of the image." },
 { "scatter(.)",
	"Noise level estimation for the image." },
 { "corr(.,.)",
	"Correlation between two images." },
 { "smooth(<image>,<sigma>,<size>)",
	"The smoothed variant of the first argument (threated as an image), "
	"which is yielded by the convolution with a Gaussian profile with a "
	"standard deviation of <sigma>, evaluated on a kernel stamp "
	"with the size of <size>." },

 { "Examples:",NULL },
 { "fiarith \"'img1.fits'-'orig.fits'\" -o diff.fits",
	"Calculates the per pixel difference between images ``img1.fits'' "
	"and ``orig.fits'' and stores the result in ``diff.fits''." },
 { "fiarith -s 512,512 \"137.42\" -o constant.fits",
	"Creates a new image with the size of 512 by 512 pixels with the constant "
	"value of 137.42." },
 { "fiarith \"'-'*26\" -o -",
	"Multiplies the pixel intensities of the image read from the standard "
	"input by 26 and writes the result to the standard output. The ``-o -'' can "
	"either be omitted, since by default the output is written to stdout." },
 { "fiarith -s 512,512 \"['input.fits'](a/(1-0.3*(X^2-Y^2)))\" -o flattened.fits",
	"This invocation evaluates the ``a/(1-0.3*(X^2+Y^2))'' expression "
	"on each pixel of the image ``input.fits'', where ``a'' referes to the "
	"pixel value itself and X and Y are the normalized spatial coordinates. "
	"The normalization is always done as X=(2*x-SX)/(SX) and Y=(2*y-SY)/SX, " 
	"i.e. the X and Y coordinates are zero at the center of the images, X has "
	"a unity absolute value at the left and right edges and Y is scaled "
	"properly to X." },
 { "fiarith \"('raw.fits'-'bias.fits'-0.5*'dark.fits')*(norm('flat.fits')/'flat.fits')\" -o calibrated.fits",
	"This operation performs a simple bias, dark and flat calibration on the "
	"input raw image ``raw.fits'', expecting ``bias.fits'', ``dark.fits'' and "
	"``flat.fits'' to be the master bias, dark and flat calibration frames, "
	"respectively; and assuming an exposure time for the raw image which is "
	"the half of the exposure time of the master dark frame. "
	"Due to the ``norm(...)'' function, the mean intensity level of the image "
	"is kept to be the same as it is in the raw image." },

 { NULL,NULL }
};

int fprint_fiarith_long_help(FILE *fw)
{
 fprintf(fw,"Usage: fiarith [options] \"<expression>\" [-o|--output <output.fits>]\n");
 fprintf(fw,
"The purpose of this program is to evaluate <expression> which is intented to\n"
"contain arithmetic operations and functions on FITS images (of the same size).\n"
"See the section ``Examples'' below for examples on the syntax of <expression>\n\n");

 longhelp_fprint(fw,fiarith_long_help,0,-1);

 fprintf(fw,"\n");
 fprintf(fw,"Report bugs to <%s>\n",FI_MAINT_EMAIL);

 return(0);
}

int main(int argc,char *argv[])
{
 fits		*img;
 FILE		*fr,*fw;
 char		*headerfile,*outfile,**inmasklist,*outmaskname,*fdpstring;
 int		i,j,headeropid,nsx,nsy,apply_mask;
 char		*expr,name[256],buff[256];
 psnsym		*mysyms[5],*myfunctsyms[4],*csym,*udfsym;
 psn		*w,*pexpr;
 evaldata	ed_data,*ed=&ed_data;
 int		is_help;
 fitsdataparam  fdp;

 progbasename=strrchr(argv[0],'/');
 if ( progbasename != NULL )	progbasename++;
 else				progbasename=argv[0];

 outmaskname=headerfile=outfile=NULL;expr=NULL;
 inmasklist=NULL;nsx=nsy=0;apply_mask=0;is_help=0;

 fdpstring=NULL;
 fdp.bitpix=-32;fdp.is_scale=0;fdp.bscale=1.0;fdp.nquantizebit=0;

 i=scanarg(argc,argv,SCANARG_WILDCARD_DASH, 
	"--version:%NS-1f%q",&is_help,
        "--version-short|--short-version:%NS-2f%q",&is_help,
	"-h|--help|--short-help|--help-short:%f%q",&is_help,
	"--long-help|--help-long:%SN2f%q",&is_help,
	"-o|--output:%s",&outfile,
	"-e|--header:%s",&headerfile,
	"-b|--bitpix:%d",&fdp.bitpix,
	"-D|--data:%s",&fdpstring,
	"-s|--size:%d,%d",&nsx,&nsy,
	"-M|--input-mask:%t",&inmasklist,
	"--output-mask:%s",&outmaskname,
	"-a|--apply-mask:%f",&apply_mask,
	"-[a-zA-Z]|--*:%e", 
	"*:%a",&expr,
	NULL);

 if ( i )	
  {	fprint_error("invalid command line argument near '%s'",argv[i]);
	return(1);
  }
 if ( is_help<0 )
  {	fprint_generic_version(stdout,argv[0],"fiarith",FI_ARITH_VERSION,is_help);
	return(0);
  }
 else if ( is_help>1 )
  {	fprint_fiarith_long_help(stdout);
	return(0);
  }
 else if ( is_help )
  {	fprint_fiarith_usage(stdout);
	return(0);
  }

 if ( expr==NULL || strlen(expr)<=0 )
  {	fprint_error("no expression has been defined");
	return(1);
  }

 if ( parse_fits_data_param(fdpstring,&fdp) )
  {	fprint_error("invalid pixel data format");
	return(1);
  }

 ed->nop =0,ed->ops =NULL;
 ed->nudf=0;ed->udfs=NULL;
 headeropid=-1;ed->mask=NULL;
 udfsym=NULL;

 myfunctsyms[0]=psn_subimg_var;
 myfunctsyms[1]=psn_general_op;
 myfunctsyms[2]=psn_general_fn;
 myfunctsyms[3]=NULL;

 for ( i=0 ; i<strlen(expr) ; )
  {	if ( expr[i]=='\'' )
	 {	int	j,k,is_read,imgid;

		name[0]=0;j=-1;k=0;
		j=i;i++;
		for ( k=0 ; i<strlen(expr) ; i++ )
		 {	if ( expr[i]=='\'' )	break;
			name[k]=expr[i];k++;
		 }
		if ( expr[i]=='\'' )
			i++;
		else
		 {	fprint_error("unmatched input image quotation mark in arthmetic expression");
			return(1);
		 }

		name[k]=0;

		if ( strcmp(name,"-")==0 )	is_read=1;
		else				is_read=0;

		fr=fopenread(name);
		if ( fr==NULL )
		 {	fprint_error("file '%s' can't be opened for reading",name);
			return(1);
		 }
		if ( ! is_read )
		 {	img=fits_create();
			fits_read_header(fr,img);
			if ( img==NULL )
			 {	fprint_error("file '%s' has invalid format",name);
				return(1);
			 }
			k=fits_get_image_params(img);
			if ( k || img->i.dim != 2 )
			 {	fprint_error("file '%s' has invalid header",name);
				return(1);
			 }
	 	 }
		else
		 {	img=fits_read(fr);
			if ( img==NULL || img->i.dim != 2 )
			 {	fprint_error("unable to read data from file '%s'",name);
				return(1);
			 }
			fits_rescale(img);
		 }
		fclose(fr);

		if ( ed->nop==0 )
		 {	ed->ops=(operand *)malloc(sizeof(operand));
			ed->sx=img->i.sx,
			ed->sy=img->i.sy;
			ed->mask=fits_mask_read_from_header(&img->header,ed->sx,ed->sy,NULL);
			fits_mask_mark_nans(&img->i,ed->mask,MASK_NAN);
		 }
		else
		 {	ed->ops=(operand *)realloc(ed->ops,sizeof(operand)*(ed->nop+1));
			if ( ed->sx != img->i.sx || ed->sy != img->i.sy )
			 {	fprintf(stderr,"File '%s' has different size (not %dx%d), exiting.\n",name,ed->sx,ed->sy);
				exit(1);
			 }
			fits_mask_mask_from_header(ed->mask,&img->header,ed->sx,ed->sy,NULL);
		 }
		imgid=ed->nop+1;

		sprintf(buff,"img%d",imgid);
		strcpy(ed->ops[ed->nop].symname,buff);
		strcpy(ed->ops[ed->nop].filename,name);
		ed->ops[ed->nop].img=img;
		ed->ops[ed->nop].is_read=is_read;
		memmove(expr+j,expr+i,strlen(expr+i)+1);
		memmove(expr+j+strlen(buff),expr+j,strlen(expr+j)+1);
		memcpy(expr+j,buff,strlen(buff));
		i=j+strlen(buff);
		ed->nop++;
		if ( headerfile != NULL )
		 {	if ( strcmp(name,headerfile)==0 )
				headeropid=ed->nop-1;
		 }
	 }
	else if ( expr[i]=='[' )
	 {	int	j,k,l,m,n,level;

		j=i;i++;level=1;n=1;
		for ( k=0 ; i<strlen(expr) ; i++ )
		 {	if ( expr[i]=='[' )			level++;
			else if ( expr[i]==']' )		level--;
			if ( level==0 )				break;
			else if ( expr[i]==',' && level==1 )	n++;
			k++;
		 }
		if ( expr[i] != ']' )
		 {	fprint_error("arithmetic expression: unmatched '['");
			return(1);
		 }
		expr[j]='(',expr[i]=')';
		i++;
		while ( expr[i]==32 || expr[i]==9 || expr[i]==10 )	i++;
		if ( expr[i] != '(' )
		 {	fprint_error("arithmetic expression: expected '(' is missing");
			return(1);
		 }
		l=i;i++;level=1;
		for ( m=0 ; i<strlen(expr) ; i++ )
		 {	if ( expr[i]=='(' )		level++;
			else if ( expr[i]==')' )	level--;
			if ( level==0 )			break;
			m++;
		 }
		if ( expr[i] != ')' )
		 {	fprint_error("arithmetic expression: expected ')' is missing");
			return(1);
		 }
		i++;
		m+=2;	/* include ( and ) */
		k=expr[i],expr[i]=0;
		w=psn_conv_string(expr+l,myfunctsyms);
		expr[i]=k;
		memmove(expr+l,expr+l+m,strlen(expr+l+m)+1);
		i-=m;
		if ( w==NULL )	
		 {	fprint_error("subexpression: symbolic error");
			return(1);
		 }
		psn_init(w,psn_general_prop);
		pexpr=psn_conv(w,psn_general_prop);
		psn_free(w);
		if ( pexpr==NULL || psn_test(pexpr) )
		 {	fprint_error("subexpression: parse error");
			return(1);
		 }
		ed->udfs=(psn **)realloc(ed->udfs,sizeof(psn *)*(ed->nudf+1));
		ed->udfs[ed->nudf]=pexpr;
		sprintf(name,"udf%d",ed->nudf+1);
		m=strlen(name);
		l=strlen(expr);
		expr=(char *)realloc(expr,l+m+1);
		memmove(expr+j+m,expr+j,l-j+1);
		memcpy(expr+j,name,m);
		i=j;
		k=ed->nudf;
		udfsym=(psnsym *)realloc(udfsym,sizeof(psnsym)*(k+2));
		udfsym[k].type=T_FN;
		udfsym[k].major=UDF_OFFSET+k;
		udfsym[k].name=(char *)malloc(m+1);
		strcpy(udfsym[k].name,name);
		udfsym[k].minor=n;
		ed->nudf++;
	 }
	else
	 {	i++;
		continue;
	 }
  }

 if ( udfsym != NULL )
  {	udfsym[ed->nudf].type=0;
	udfsym[ed->nudf].major=0;
	udfsym[ed->nudf].name=NULL;
	udfsym[ed->nudf].minor=0;
  }

 if ( ed->nop==0 )
  {	ed->sx=nsx,ed->sy=nsy;
	if ( ed->sx<=0 || ed->sy<=0 )
	 {	fprint_error("unexpected or invalid image size");
		return(1);
	 }
  }

 if ( ed->mask==NULL )
  {	ed->mask=fits_mask_create_empty(ed->sx,ed->sy);		}
 if ( inmasklist != NULL )
  {	if ( join_masks_from_files(ed->mask,ed->sx,ed->sy,inmasklist) )
	 {	fprint_error("unable to read one of the input masks");
		return(1);
	 }
  }

 if ( headerfile != NULL && headeropid<0 )
  {	fprintf(stderr,"Warning: file specified after -e not found in the expression, the\n");
	fprintf(stderr,"         header of the output may not contain the expected information!\n");
  }

 csym=(psnsym *)malloc(sizeof(psnsym)*(ed->nop+1));
 for ( i=0 ; i<ed->nop ; i++ )
  {	csym[i].type=T_VAR;
	csym[i].major=i;
	csym[i].name=ed->ops[i].symname;
  }
 csym[ed->nop].type=0;
 csym[ed->nop].major=0;
 csym[ed->nop].name=NULL;
 csym[ed->nop].minor=0;

 mysyms[0]=psn_img_op;
 mysyms[1]=psn_img_fn;
 mysyms[2]=csym;
 if ( udfsym != NULL )
  {	mysyms[3]=udfsym;
	mysyms[4]=NULL;
  }
 else
	mysyms[3]=NULL;

 w=psn_conv_string(expr,mysyms);
 if ( w==NULL )
  {	fprint_error("arithmetic expression: symbolic error");
	return(1);
  }
 psn_init(w,psn_img_prop);
 pexpr=psn_conv(w,psn_img_prop);
 psn_free(w);
 if ( pexpr==NULL || psn_test(pexpr) )
  {	fprint_error("arithmetic expression: parse error");
	return(1);
  }

 img=evaluate(ed,pexpr,headeropid);

 psn_free(pexpr);

 img->i.bit=fdp.bitpix;
 img->i.curr.bscale=1.0;
 img->i.curr.bzero=0.0;
 if ( fdp.is_scale )		img->i.read.bscale=fdp.bscale,img->i.read.bzero=fdp.bzero;
 else if ( img->i.bit<0 )	img->i.read.bscale=1.0,img->i.read.bzero=0.0;
 else				img->i.read.bscale=1.0,img->i.read.bzero=(1<<(img->i.bit-1));

 if ( ed->hdrsave != NULL )
	fits_copy_full_header(img,ed->hdrsave);

 fits_set_image_params(img);
/* mask=fits_mask_create_floyd(img,ed->sx,ed->sy,1,100); */ /* for testing */

 fits_history_export_command_line(img,"fiarith",FI_ARITH_VERSION,argc,argv);

 fits_backscale(img,img->i.read.bscale,img->i.read.bzero);
 mark_integerlimited_pixels(&img->i,ed->mask,img->i.bit,1,MASK_OVERSATURATED,MASK_OVERSATURATED);

 fits_mask_export_as_header(&img->header,1,ed->mask,ed->sx,ed->sy,NULL);
 for ( i=0 ; i<ed->sy && apply_mask ; i++ )
  {	for ( j=0 ; j<ed->sx ; j++ )
	 {	if ( ed->mask[i][j] )	img->i.data[i][j]=0.0;	}
  }
 
 if ( fdp.nquantizebit>0 )
	fits_image_quantize(&img->i,fdp.nquantizebit);

 if ( outfile==NULL )	fw=stdout;
 else			fw=fopenwrite(outfile);
 if ( fw==NULL )
  {	fprint_error("unable to create output image");
	return(1);
  }
 fits_write(fw,img);
 fclosewrite(fw);

 if ( outmaskname != NULL )
  {	fw=fopenwrite(outmaskname);
	if ( fw==NULL )
	 {	fprint_error("unable to create output mask file");
		return(1);
	 }
	fits_write_header(fw,img);
	fclosewrite(fw);
  }

 fits_free(img);

 return(0);
}

/*****************************************************************************/
                 
