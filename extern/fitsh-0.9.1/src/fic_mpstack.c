/*****************************************************************************/
/* fic_mpstack.c							     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* A contributed utility to the `fitsh` package aided for faint minor planet */
/* search projects.							     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2004, 2007, 2008, 2010; Pal, A. (apal@szofi.net)			     */
/*****************************************************************************/
#define	FI_C_MPSTACK_VERSION	"0.1"
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
#include "fbase.h"
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

int fprint_fic_mpstack_usage(FILE *fw)
{
 fprintf(fw,
"Usage:\tfic_mpstack [-h|--help] [-C|--comment] [-V|--verbose]\n"
"\t-f|--offset <ox>,<oy> -s|--size <sx>,<sy>\n"
"\t-t|--shift <minx>:<stepx>:<maxx>,<miny>:<stepy>:<maxy>\n"
"\t[-b|--border <borderwidth>] [-o|--output <output.fits>]\n"
"\t<input.fits> ...\n");
 
 return(0);
}

longhelp_entry fic_mpstack_long_help[]=
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
 { "-it, --input-transformation <transformation file>",
	"Name of the file which contains the transformation description."
	"Such a file can be created e.g. by the programs `grtrans` or `grmatch`. "
	"This file contains basically the same set of <keyword> = <value> "
	"pairs as it is used after the -t|--transformation option (see there)." },
 { "-t, --transformation <transformation>",
	"Comma-separated list of parameters for the spatial transformation, "
	"see section ``Parameters for spatial transformations'' below." },

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


int fprint_fic_mpstack_long_help(FILE *fw)
{
 fprintf(fw,
"Usage:\tfic_mpstack [stacking options] <inputs> [-o|--output <output>]\n"
"The main purpose of this program is to perform a special kind of \n"
"4 dimensional transformation and image stacking that aids faint minor planet \n"
"searches.\n");
 fprintf(fw,"\n");

 longhelp_fprint(fw,fic_mpstack_long_help,0,-1);

 fprintf(fw,"\n");
 fprintf(fw,"Report bugs to <%s>\n",FI_MAINT_EMAIL);

 return(0);

}

/*****************************************************************************/

typedef struct
 {	double	**image;
	char	**mask;
	double	**coeff;
 } stack;

int fic_mpstack_create(double **out,stack *images,int nimage,
	int nsx,int nsy,int ssx,int ssy,
	int rx,int minx,int stepx,int ry,int miny,int stepy,
	int bwidth,double bgvalue)
{
 double	*arr,m;
 int	i,j,k,l,ox,oy,n,w;
 double	dx,dy,ax,ay,ix,iy;
 stack	*s;

 arr=tensor_alloc_1d(arr,nimage);

 for ( i=0 ; i<ry ; i++ )
  {	oy=i*(nsy+2*bwidth);
	dy=miny+i*stepy;
	for ( j=0 ; j<rx ; j++ )
	 {	ox=j*(nsx+2*bwidth);
		dx=minx+j*stepx;
		for ( k=0 ; k<nsy+2*bwidth ; k++ )
		 {	for ( l=0 ; l<nsx+2*bwidth ; l++ )
			 {	out[oy+k][ox+l]=bgvalue;		}
		 }
		for ( k=0 ; k<nsy ; k++ )
		 {	for ( l=0 ; l<nsx ; l++ )
			 {	n=0;
				for ( w=0 ; w<nimage ; w++ )
				 {	s=&images[w];
					ax=dx*((double)w/(double)(nimage-1)-0.5);
					ay=dy*((double)w/(double)(nimage-1)-0.5);
					ix=l+ax-minx;
					iy=k+ay-miny;
					if ( ! ( 0<=ix && ix<=ssx-1 && 0<=iy && iy<=ssy-1 ) )
						continue;
					arr[n]=bicubic_inter(s->coeff,ix,iy);
					n++;
				 }
				if ( n>0 )
					m=median(arr,n);
				else
					m=0.0;
				out[oy+bwidth+k][ox+bwidth+l]=m;
			 }
		 }
	 }
  }

 tensor_free(arr);

 return(0);
}

int main(int argc,char *argv[])
{
 FILE		*fw;
 fits		*outimg;
 int		i,is_help;
 char		**inputfiles;
 int		ofx,ofy,nsx,nsy,tsx,tsy;
 int		minx,stepx,maxx,
		miny,stepy,maxy;
 int		rx,ry,ssx,ssy,bwidth,asx,asy;
 char		*outimgfile;
 double		bgvalue;

 stack		*images;
 int		nimage;

 progbasename=strrchr(argv[0],'/');
 if ( progbasename != NULL )	progbasename++;
 else				progbasename=argv[0];

 inputfiles=NULL;
 outimgfile=NULL;
 is_comment=is_verbose=is_help=0;
 ofx=ofy=nsx=nsy=0;

 bwidth=0;
 minx=stepx=maxx=0;
 miny=stepy=maxy=0;

 bgvalue=0;

 i=scanarg(argc,argv,SCANARG_ALLOW_FLAGS,
	"--version:%NS-1f%q",&is_help,
        "--version-short|--short-version:%NS-2f%q",&is_help,
        "--long-help|--help-long:%SN2f%q",&is_help,
	"-h|--help:%f%q",&is_help,
	"-f|--offset:%d,%d",&ofx,&ofy,
	"-s|--size:%d,%d",&nsx,&nsy,
	"-t|--shift:%d:%d:%d,%d:%d:%d",&minx,&stepx,&maxx,&miny,&stepy,&maxy,
	"-b|--border:%d",&bwidth,
	"-o|--output:%s",&outimgfile,
	"-w|--background:%g",&bgvalue,
	"-*|+*:%e",
	"*:%l",&inputfiles,
	NULL);

 if ( i )		
  {	fprint_error("invalid command line argument near '%s'",argv[i]);
	return(1);
  }
 else if ( is_help<0 )
  {	fprint_generic_version(stdout,argv[0],"fic_mpstack",FI_C_MPSTACK_VERSION,is_help);
	return(0);
  }
 else if ( is_help>1 )
  {	fprint_fic_mpstack_long_help(stdout);
	return(0);
  }
 else if ( is_help )
  {	fprint_fic_mpstack_usage(stdout);
	return(0);
  }

 if ( nsx<=0 || nsy<=0 )
  {	fprint_error("undefined or unexpected size values");
	return(1);
  }

 if ( stepx<=0 || stepy<=0 )
  {	fprint_error("undefined or unexpected grid step values");
	return(1);
  }
 if ( maxx<minx )	i=minx,minx=maxx,maxx=i;
 if ( maxy<miny )	i=miny,miny=maxy,maxy=i;
 rx=(maxx-minx)/stepx+1;
 ry=(maxy-miny)/stepy+1;

 asx=rx*(2*bwidth+nsx);
 asy=ry*(2*bwidth+nsy); 

 ssx=nsx-minx+maxx;
 ssy=nsy-miny+maxy;

 for ( i=0 ; inputfiles != NULL && inputfiles[i] != NULL ; )	i++;
 nimage=i;
 if ( nimage<2 )
  {	fprint_error("too few images: we need at least two ones");
	return(1);
  }

 images=(stack *)malloc(sizeof(stack)*nimage);

 tsx=tsy=0;
 for ( i=0 ; i<nimage ; i++ )
  {	FILE	*fr;
	fits	*img;
	char	**mask;
	stack	*s;
	int	k,l,ix,iy;

	if ( (fr=fopenread(inputfiles[i]))==NULL )
	 {	fprint_error("unable to open input file '%s'",inputfiles[i]);
		return(1);
	 }	
	img=fits_read(fr);
	fcloseread(fr);
	if ( img==NULL )
	 {	fprint_error("unable to interpret input from '%s' as FITS data",inputfiles[i]);
		return(1); 
	 }
	else if ( img->i.dim != 2 )
	 {	fprint_error("image dimension differs from 2");
		return(1); 
	 }
	if ( i==0 )
	 {	tsx=img->i.sx;
		tsy=img->i.sy;
	 }
	else if ( tsx != img->i.sx || tsy != img->i.sy )
	 {	fprint_error("FITS data '%s' has different size",inputfiles[i]);
		return(1);
	 }
	fits_rescale(img);
	mask=fits_mask_read_from_header(&img->header,img->i.sx,img->i.sy,NULL);

	s=&images[i];

	s->image=tensor_alloc_2d(double,ssx,ssy);
	s->mask=fits_mask_create_empty(ssx,ssy);
	s->coeff=tensor_alloc_2d(double,2*ssx,2*ssy);

	for ( k=0 ; k<ssy ; k++ )
	 {	iy=ofy+miny+k;
		for ( l=0 ; l<ssx ; l++ )
		 {	ix=ofx+minx+l;
			if ( ! ( 0<=ix && ix<tsx && 0<=iy && iy<tsy ) )
			 {	s->mask[k][l]=MASK_OUTER;
				s->image[k][l]=0.0;
			 }
			else
			 {	s->mask[k][l]=mask[iy][ix];
				s->image[k][l]=img->i.data[iy][ix];
			 }
		 }
	 }	

	fits_mask_free(mask);
	fits_free(img);

  }

 for ( i=0 ; i<nimage ; i++ )
  {	stack	*s;
	s=&images[i];
	bicubic_coeff(s->image,ssx,ssy,s->coeff,s->mask);
  }

 outimg=fits_create();
 fits_set_standard(outimg,NULL);
 outimg->i.bit=-32;
 outimg->i.curr.bscale=1.0;
 outimg->i.curr.bzero=0.0;
 fits_alloc_image(outimg,asx,asy);

 fic_mpstack_create(outimg->i.data,images,nimage,nsx,nsy,ssx,ssy,rx,minx,stepx,ry,miny,stepy,bwidth,bgvalue);
 
 outimg->i.read.bscale=1.0,
 outimg->i.read.bzero=0.0;
 fits_set_image_params(outimg);
 fits_backscale(outimg,outimg->i.read.bscale,outimg->i.read.bzero);
 fits_history_export_command_line(outimg,"fic_mpstack",FI_C_MPSTACK_VERSION,argc,argv);
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

