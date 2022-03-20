/*****************************************************************************/
/* fiinfo.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line user interface to get some statistics from FITS data.	     */
/*****************************************************************************/
#define	FI_INFO_VERSION		"0.9d1"
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>
#include <stdint.h>
#include <inttypes.h>

#include <fits/fits.h>

#include "longhelp.h"
#include "fi.h"

#include "fitsmask.h"
#include "math/spline/biquad.h"
#include "math/spline/spline.h"
#include "math/fit/lmfit.h"
#include "statistics.h"
#include "math/poly.h"
#include "math/polyfit.h"
#include "io/iof.h"
#include "io/scanarg.h"
#include "io/tokenize.h"
#include "link/linkpoint.h"
#include "tensor.h"
#include "common.h"

#include "fiinfo.h"

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

int fprint_fiinfo_usage(FILE *fw)
{
 fprintf(fw,
"Usage:\tfiinfo [<in>[<F>] [--ignore-mask]] [--frame <F>]] \n"
"\t[-o|--output <out.dat>] [-h|--help] [--[long-]summary]\n");
 fprintf(fw,
"\t[-s|--statistics [mean|median],[iterations=<n>],[lower|upper|sigma=<s>]]\n"
"\t[-d min,max,mean,stddev,sky,skysigma,bgdump,bgfit,bgsky,bgskysigma] [-n]\n"
"\t[-b|--box <num-of-blocks-for-all-bg...>] [-a|--order <order-for-bgfit>]\n"
"\t[-od|--output-dump <pixels.dump> [-m|--dump--mask]]\n");
 fprintf(fw,
"Generating PGM/PPM images:\n"
"\t[-op|--output-pnm <img.pnm>]\n"
"\t[--p[pg]m {[defaults],[reverse],[contrast=<c>],[brightness=<b>]\n"
"\t  [palette=<palette-spec>],[linear|log|squared|sqrt|histequ],\n"
"\t  [minmax|percentage=<%%>|min=<l>,max=<h>|\n"
"\t   zscale|zmax|zmin,[zcontrast=<zc>]]}]]\n");

 return(0);
}

longhelp_entry fiinfo_long_help[]=
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
 { "--summary, --long-summary",
	"Give a summary about the content structure of the FITS file. Namely, "
	"the main characeristics and dimensions of the primary image, followed "
	"by the list of optional extensions and their main properties." },
 { "-s, --statistics <list of statistics>",
	"Calculate basic statistics for the image. See ``Statistics options'' "
	"below for available statistics methods. The statistics itself "
	"containts four numbers: the total number of pixels involved in the "
	"calculations, the final number of pixels used for the statistics "
	"(which might be smaller than the previous one if the outliers are "
	"rejected), an average value and a scatter." },
 { "-d, --data <list of derived image data>",
	"Calculate some other more quantities related to astronomical images. "
	"This option should be followed by a comma-separated list of quantities. "
	"See ``Image characteristics'' below for more details about these." },
 { "-b, --box <number of blocks>",
	"This option specifies the number of blocks, which is used to divide "
	"the input image. Some quantities (see ``Image characteristics'') "
	"can be derived on a per block basis either." },
 { "-a, --order <order>",
	"Order of polynomial spatial variations in some derived image "
	"characterization quantities (see also ``Image characteristics'')." },
 { "-n, --newline",
	"In the output, each quantity should be written in separate lines. "
	"By default, the output is a single line, containing the desired "
	"quantities or statistics." },
 { "--ignore-mask",
	"Completely ignore the mask associated to the input image." },
 { "-od, --output-dump <file>",
	"Name of an output file in which a raw image dump is written. "
	"Each line of this file contains 3 or 4 columns: X, Y coordinates "
	"and flux, optionally followed by the associated mask flag "
	"(see also ``-m|--dump-mask'')." },
 { "-m, --dump-mask",
	"The raw image dump specified by ``-od|--output-dump'' should contain "
	"the masking information beyond the coordinates and intensities." },
 { "-op, --output-pnm, --output-ppm, --output-pgm",
	"Name of an output file in which the image is stored in a variant of "
	"PNM format. These images are intented to be a kind of ``human "
	"visible'' images, appropriately scaled for normal displays. "
	"These images are stored in PNM format, which is an easily parseable "
	"(thus raw, uncompressed) format, supported by many graphic programs "
	"(and by the NETPBM package). Such an image conversation always "
	"results data loss. "
	"See also options ``--pgm'' or ``--ppm'' for further details." },
 { "--pgm <PGM specific conversion options>",
	"This command line argument is followed by a comma-separated list "
	"of options, which specifies the scaling and other properties of the "
	"output image. The resulted image will be a grey-scale (PGM) image, "
	"even if a color palette is requested. See ``PNM specifications'' below." },
 { "--ppm <PPM specific conversion options>",
	"This command line argument is followed by a comma-separated list "
	"of options, which specifies the scaling and other properties of the "
	"output image. The resulted image will be a tru-color (PPM) image, "
	"even if a greyscale colormap is requested. See also "
	"``PNM specifications'' below for more details." },
 
 { "Statistics options: ", NULL },
 { "mean",
	"The mean value of the pixel intensities." },
 { "median",
	"The median value of the pixel intensities." },
 { "iterations",
	"Reject the outlier pixels before doing any statistics." },
 { "lower=<sigma>, upper=<sigma>, sigma=<sigma>",
	"Lower, upper or common rejection level, in the units of "
	"standard deviation (which is derived around the mean or median value, "
	"depending on the request of the user)." },

 { "Image characteristics:" ,NULL },
 { "min, max",
	"Minimal and maximal pixel intensities on the image." },
 { "mean",
	"Mean intensity level." },
 { "stddev",
	"Standard deviation." },
 { "sky",
	"Sky background level." },
 { "skysigma",
	"Sky background scatter." },

 { "PNM specifications:", NULL },
 { "linear",
	"Use a linear intensity scaling." },
 { "log",
	"Use a logarithmic intensity scaling." },
 { "squared",
	"Use a squared intensity scaling." },
 { "sqrt",
	"Use a square root intensity scaling." },
 { "histequ",
	"Use a histogram equalized intensity scaling." },
 { "minmax",
	"Use the minimal and maximal pixel intensities for scaling boundaries." },
 { "percentage=<%>",
	"Use the minimal and maximal values of the innermost "
	"specified percent of the pixel intensities." },
 { "min=<min>, max=<max>",
	"Use the specified minimal and maximal values for scaling boundaries." },
 { "zscale",
	"Use the ``zscale'' algorithm to determine scaling boundaries." },
 { "zmax, zmin",
	"Use the ``zmax'' or ``zmin'' algorithm to determine scaling boundaries." },
 { "zcontrast=<zcontrast>",
	"Use the specified contrast value to determine the scaling boundaries "
	"in the case of ``zscale'', ``zmax'' or ``zmin'' methods. The default "
	"value is 0.25." },
 { "reverse",
	"Use an inverted color map." },
 { "contrast=<C>, brightness=<B>",
	"Use the specified values for adjusting the final contrast and "
	"brightness. The default values are 1 and 0.5, respectively, "
	"according to the standard image contrast and "
	"brightness level definitions." },
 { "palette=<color1>:<color2>:<color3>:...",
	"Specify an alternate color map. Each color should be a hexadecimal "
	"representation of a given color, i.e. it should be in one of the "
	"forms of G, GG, RGB or RRGGBB, denoting 4 bit grey, 8 bit grey, "
	"3x4 bit truecolor or 3x8 bit truecolor representation, respectively. "
	"The color map gradient will be continuous if the colors are separated by "
	"colons. Jumps in the gradient can be defined by separating the "
	"successive colors by a slash, ``/''." },
 { "Note that the syntax followed by the ``--pgm'' or ``--ppm'' command line "
   "arguments is exactly the same for both options. However, color images "
   "will be converted to greyscale if ``--pgm'' is specified, and vice versa, "
   "``--ppm'' always yields a PPM format, even if the color gradient is "
   "merely a grayscale one. ", NULL },
 
 { NULL,NULL }
};
int fprint_fiinfo_long_help(FILE *fw)
{
 fprintf(fw,
"Usage:\tfiinfo [options] [-i <input>] <outputs>\n"
"The main purpose of the `fiinfo` program is to give some information about "
"the FITS files (primarily FITS images, but output dump is supported for "
"tables and binary tables also).\n");
 fprintf(fw,"\n");

 longhelp_fprint(fw,fiinfo_long_help,0,-1);

 fprintf(fw,"\n");
 fprintf(fw,"Report bugs to <%s>\n",FI_MAINT_EMAIL);

 return(0);
}
/*
void fprint_help(FILE *fw)
{
 fprintf(fw,
"<in.fits>: name of input file (optional; by default, data read from stdin).\n"
"<out.dat>: name of output file (optional; by default, data written to stdout).\n"
"-h|--help: prints info about the usage of the program.\n"
"-n       : prints a newline character after each data (except 'bdump').\n"
"-d <...> : list of statistics to be calculated (and written to the output):\n");
 fprintf(fw,
" * min, max, mean, stddev: the minimal, maximal and the mean values of the\n"
"   pixels in the image, and the standard deviaton of the pixels;\n"
" * sky, skysigma: naive approximation for the background and the scatter of\n"
"   the background (reliable only for sparse fields);\n");
 fprintf(fw,
" * bdump, bfit, bsky, bskysigma: divides the image into\n"
"   B times B blocks (B can be defined with -b|--box, default value is 16),\n"
"   estimates the desired parameters for each block and averages them.\n"
"   The 'bgfit' method fits a polynomial to the sky values (up to the order\n"
"   defined by -a|--order), the 'bgdump' method dumps the estimated sky\n"
"   and skysigma values for each block.\n");
}
*/

/*****************************************************************************/

int format_check(char *f)
{
 while ( *f=='-' || *f=='+' )	f++;
 while ( isdigit(*f) )		f++;
 if ( *f=='+' )			f++;
 if ( *f=='.' )			f++;
 if ( *f=='+' )			f++;
 while ( isdigit(*f) )		f++;
 if ( *f=='e' || *f=='E' || *f=='f' || *f=='F' || *f=='g' || *f=='G' )
	return(0);
 else
	return(1);
}

/*****************************************************************************/

int strmask(char *mbuff,int k)
{
 strcpy(mbuff,"--------");
 if ( k & MASK_FAULT )		mbuff[0]='f';
 if ( k & MASK_HOT )		mbuff[1]='h';
 if ( k & MASK_COSMIC )		mbuff[2]='c';
 if ( k & MASK_OUTER )		mbuff[3]='o';
 if ( k & MASK_OVERSATURATED )	mbuff[4]='s';
 if ( k & MASK_LEAKED )		mbuff[5]='l';
 if ( k & MASK_SATURATED )	mbuff[6]='S';
 if ( k & MASK_INTERPOLATED )	mbuff[7]='i';
 return(0);
}

#define		DAT_MIN		 0
#define		DAT_MAX		 1
#define		DAT_MEAN	 2
#define		DAT_STDDEV	 3
#define		DAT_MEDIAN	11
#define		DAT_SKY		 4
#define		DAT_SKYSIGMA	 5
#define		DAT_BGDUMP	 6
#define		DAT_BGFIT	 8
#define		DAT_BGSKY	 9
#define		DAT_BGSKYSIGMA	10

typedef struct
 {	char	*name;
	int	code;
 } dataname;

static dataname	datanamelist[]=
 {	{ "min",	DAT_MIN		},
	{ "max",	DAT_MAX		},
	{ "mean",	DAT_MEAN	},
	{ "stddev",	DAT_STDDEV	},
	{ "median",	DAT_MEDIAN	},
	{ "sky",	DAT_SKY		},
	{ "skysigma",	DAT_SKYSIGMA	},
	{ "bgdump",	DAT_BGDUMP	},
	{ "bgfit",	DAT_BGFIT	},
	{ "bgsky",	DAT_BGSKY	},
	{ "bgskysigma",	DAT_BGSKYSIGMA	},
	{ NULL,		-1		}
 };


int parse_dump_data(char *dumpstr,int **rdats,int *rndat)
{
 int		*dats;
 int		ndat,l,j;
 char		buff[16];
 dataname	*dnames,*wdd;

 dats=NULL;
 ndat=0;

 dnames=datanamelist;

 while ( *dumpstr )
  {	l=0;
	while ( *dumpstr && *dumpstr != ',' && l<15 )
	 {	buff[l]=*dumpstr,
		l++,dumpstr++;
	 };
	buff[l]=0;
	if ( *dumpstr==',' )	dumpstr++;

	for ( wdd=dnames,j=-1 ; wdd->name != NULL && j<0 ; wdd++ )
	 {	if ( strcmp(buff,wdd->name)==0 )
			j=wdd->code;
	 }
	if ( j<0 )
	 {	if ( dats != NULL )	free(dats);
		return(1);
	 }

	dats=(int *)realloc(dats,sizeof(int)*(ndat+1));
	dats[ndat]=j;
	ndat++;
  };

 if ( rdats != NULL )	*rdats=dats;
 if ( rndat != NULL )	*rndat=ndat;

 return(0);
}

/*****************************************************************************/

typedef struct
 {	int	*dats;
	int	ndat;
	int	nbox,norder;
 } infoimage;

int fits_stat_v1_image(fitsimage *fi,char **mask,FILE *fw,int nwl,infoimage *ii)
{
 int		i,j;
 int		is_med,is_mord,is_bg,is_bgdump;
 double		min,max,mean,stddev,sky,skysigma,*pfit,med;
 double		bgsky,bgskysigma;

 is_med=is_mord=is_bg=is_bgdump=0;
 for ( i=0 ; i<ii->ndat && ii->dats != NULL ; i++ )
  {	j=ii->dats[i];
	if ( j==11 )				is_med=1;
	else if ( j==4 || j==5 )		is_mord=1;
	else if ( j==6 )			is_bg=is_bgdump=1;
	else if ( j==8 || j==9 || j==10 )	is_bg=1;
  }

/* Basic statistics: min, max, mean and standard deviation: */
 fits_stat_basic(fi,mask,&min,&max,&mean,&stddev);

/* Statistics which require the pixel data to be ordereded (eg. median): */
 if ( is_mord )
	fits_stat_sky_biquad(fi,stddev,&sky,&skysigma);

 if ( is_med )	med=fits_stat_median(fi);
 else		med=0.0;

/* Statistics which require the image to be divided into boxes: */
 if ( is_bg )
  {	point	**pdsky,**pdsigma;
	double	*rawdata;
	int	i,j,k,nbx,nby;

	nbx=nby=ii->nbox;

	pdsky  =(point **)tensor_alloc_2d(point,nbx,nby);
	pdsigma=(point **)tensor_alloc_2d(point,nbx,nby);

	fits_stat_background(fi,nbx,nby,pdsky,pdsigma,stddev);

	if ( is_bgdump )
	 {	for ( i=0 ; i<nby ; i++ )
		 {  for ( j=0 ; j<nbx ; j++ )
			fprintf(fw,"%11g %11g %11g %11g\n",
			pdsky[i][j].x,pdsky[i][j].y,
			pdsky[i][j].value,pdsigma[i][j].value);
		 }
	 }

	pfit =(double *)malloc(sizeof(double)*((ii->norder+1)*(ii->norder+2)/2));
	fit_2d_poly(&pdsky[0][0],nbx*nby,ii->norder,pfit,0.5,0.5,1.0);

	rawdata=(double *)malloc(sizeof(double)*(nbx*nby));
	for ( i=0,k=0 ; i<nby ; i++ )
	 {  for ( j=0 ; j<nbx ; j++,k++ )
		rawdata[k]=pdsky[i][j].value;
	 }
	bgsky=median(rawdata,nbx*nby);
	for ( i=0,k=0 ; i<nby ; i++ )
	 {  for ( j=0 ; j<nbx ; j++,k++ )
		rawdata[k]=pdsigma[i][j].value;
	 }
	bgskysigma=median(rawdata,nbx*nby);
	free(rawdata);

	tensor_free(pdsigma);
	tensor_free(pdsky);
  }
 else
  {	pfit=NULL;bgsky=bgskysigma=0.0;			}


/* Dump output data to the stream fw: */
 for ( i=0 ; i<ii->ndat ; i++ )
  {	switch ( ii->dats[i] )
	 {   case DAT_MIN:
	 	fprintf(fw,"%11g ",min);
		break;
	     case DAT_MAX:
		fprintf(fw,"%11g ",max);
		break;
	     case DAT_MEAN:
		fprintf(fw,"%11g ",mean);
		break;
	     case DAT_STDDEV:
		fprintf(fw,"%11g ",stddev);
		break;
	     case DAT_SKY:
		fprintf(fw,"%11g ",sky);
		break;
	     case DAT_SKYSIGMA:
		fprintf(fw,"%11g ",skysigma);
		break;
	     case DAT_BGFIT:
		for ( j=0 ; j<(ii->norder+1)*(ii->norder+2)/2 ; j++ )
		 {	fprintf(fw,"%11g ",pfit[j]);		}
		break;
	     case DAT_BGSKY:
		fprintf(fw,"%11g ",bgsky);
		break;
	     case DAT_BGSKYSIGMA:
		fprintf(fw,"%11g ",bgskysigma);
		break;
	     case DAT_MEDIAN:
		fprintf(fw,"%11g ",med);
		break;
	 }
	
	if ( nwl )	fprintf(fw,"\n");
  }
 if ( ! nwl && ii->ndat>0 )	fprintf(fw,"\n");

 return(0);
}

/*****************************************************************************/

typedef struct
 {	double	floydratio;
	int	use_median;
	int	niter;
	double	lower,upper;
 } imgstatparam;

static int data_add_point(double **rdata,int *rn,double d)
{
 double	*data;
 int	n;

 data=*rdata;
 n=*rn;

 if ( ! (n%256) )
	data=(double *)realloc(data,sizeof(double)*(n+256));

 data[n]=d;
 n++;

 *rdata=data;
 *rn=n;

 return(0);
}

int fits_stat_v2_image(fitsimage *img,char **inmask,FILE *fw,imgstatparam *isp,int is_comment)
{
 int	sx,sy;
 char	**mask;
 int	i,j;
 double	*data,center,s2,d,sig,d0,d1,sum;
 int	n,iiter,i0,i1,k;

 sx=img->sx;
 sy=img->sy;

 if ( isp->floydratio>0.0 && isp->floydratio<1.0 )
  {	int	a,b;
	a=isp->floydratio*32768.0;
	b=32768;

	mask=fits_mask_create_floyd(sx,sy,a,b,1);
	if ( inmask != NULL )
	 {	for ( i=0 ; i<sy ; i++ )
		 {	for ( j=0 ; j<sx ; j++ )
			 {	mask[i][j] |= inmask[i][j];		}
		 }
	 }
  }
 else
	mask=inmask;

 data=NULL; 
 n=0; 
 sum=0.0;
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( mask != NULL && mask[i][j] )
			continue;
		data_add_point(&data,&n,img->data[i][j]);
		sum+=img->data[i][j];
	 }
  }

 if ( isp->use_median )
	center=median(data,n);
 else
  {	center=sum/(double)n;
	median(data,n);
  }

 i0=0;
 i1=n;
 sig=0.0;

 for ( iiter=0 ; iiter <= (isp->niter>0?isp->niter:0) ; iiter++ )
  {	s2=0.0;
	for ( i=i0 ; i<i1 ; i++ )
	 {	d=data[i]-center;
		s2+=d*d;
	 }
	sig=sqrt(s2/(double)(i1-i0));

	if ( iiter<isp->niter )
	 {	k=0;
		d0=center-sig*isp->lower;
		while ( i0<n && data[i0]<d0 )	sum-=data[i0],i0++,k++;
		d1=center+sig*isp->upper;
		while ( 0<i1 && d1<data[i1-1] )	i1--,sum-=data[i1],k++;
		if ( isp->use_median )
			center=0.5*(data[(i1+i0)/2]+data[(i1+i0-1)/2]);
		else
			center=sum/(double)(i1-i0);

		if ( k <= 0 )
			break;
	 }
  }

 if ( is_comment )
	fprintf(fw,"# Num    used    average   scatter \n");

 fprintf(fw," %8d %8d %12.7g %12.7g\n",n,i1-i0,center,sig);

 if ( mask != NULL && mask != inmask )
 	fits_mask_free(mask);
	
 return(0);
}

/*****************************************************************************/

int fitsttable_dump(FILE *fw,fitsttable *ft)
{
 unsigned char	*line;
 int		i,j,p,len;
 char		*buff;

 if ( fw==NULL )	return(1);
 if ( ft==NULL || ft->data==NULL || ft->tfields==NULL )	return(1);

 buff=(char *)malloc(ft->rowsize+1);
 if ( buff==NULL )	return(-1);

 for ( i=0 ; i<ft->nrow ; i++ )
  {	line=ft->data[i];
	for ( j=0 ; j<ft->ntfield ; j++ )
	 {	p=ft->tfields[j].colindex;
		if ( j<ft->ntfield-1 )	len=ft->tfields[j+1].colindex-p;
		else			len=ft->rowsize-p;
		strncpy(buff,(char *)line+p,len);
		buff[len]=0;
		if ( j )	fprintf(fw," %s",buff);
		else		fprintf(fw,"%s",buff);
	 }
	fprintf(fw,"\n");
  }

 free(buff);

 return(0);
}

int fitsbtable_dump(FILE *fw,fitsbtable *fb)
{
 unsigned char	*line,*ptr,m;
 int		i,j,k,offset,repeat,basesize,form;
 char		*buff;

 if ( fw==NULL )	return(1);
 if ( fb==NULL || fb->data==NULL || fb->bfields==NULL )	return(1);

 buff=(char *)malloc(fb->rowsize+1);
 if ( buff==NULL )	return(-1);

 for ( i=0 ; i<fb->nrow ; i++ )
  {	line=fb->data[i];
	offset=0;
	for ( j=0 ; j<fb->nbfield ; j++ )
	 {	repeat=fb->bfields[j].repeat;
		form=fb->bfields[j].form;
		basesize=fits_bintable_form_basesize(form);
		if ( basesize<0 )	break;	/* unexpected TFORM... */
		else if ( basesize>0 && form != FTF_STRING )
		 {	for ( k=0 ; k<repeat ; k++ )
			 {	ptr=line+offset;
				switch ( form )
				 {   case FTF_LOGICAL:
					if ( *ptr )	fprintf(fw,"T");
					else		fprintf(fw,"F");
					break;
				     case FTF_BYTE:
					fprintf(fw,"%3d",*ptr);
					break;
				     case FTF_SHORT:
					fprintf(fw,"%6d",*(short *)ptr);
					break;
				     case FTF_LONG:
					if ( sizeof(long)==4 )
						fprintf(fw,"%11ld",*(long *)ptr);
					else	
						fprintf(fw,"%11d",*(int *)ptr);
					break;
				     case FTF_LONGLONG:
					fprintf(fw,"%21lld",*(long long *)ptr);
					break;
				     case FTF_FLOAT:
					fprintf(fw,"%12.6g",*(float *)ptr);
					break;
				     case FTF_DOUBLE:
					fprintf(fw,"%15.13g",*(double *)ptr);
					break;
				     case FTF_FLOATCOMPLEX:
					fprintf(fw,"%12.6g,%12.6g",((float *)ptr)[0],((float *)ptr)[1]);
					break;
				     case FTF_DOUBLECOMPLEX:
					fprintf(fw,"%15.13g,%15.13g",((double *)ptr)[0],((double *)ptr)[1]);
					break;
				 }
				fprintf(fw," ");	
				offset+=basesize;
			 }
		 }	
		else if ( basesize>0 )
		 {	memcpy(buff,line+offset,repeat);
			buff[repeat]=0;
			for ( k=0 ; k<repeat ; k++ )
			 {	if ( buff[k]==0 )	buff[k]=32;	}
			fprintf(fw,"%s ",buff);
			offset+=repeat*basesize;
		 }
		else 
		 {	for ( k=0 ; k<repeat ; k++ )
			 {	m=128>>(k%8);
				if ( line[offset+(k/8)] & m )	fprintf(fw,"1");
				else				fprintf(fw,"0");
			 }
			offset+=(repeat+7)/8;
		 }
	 }
	fprintf(fw,"\n");
  }

 free(buff);

 return(0);
}

/*****************************************************************************/

int fprint_fitsimage_summary(FILE *fw,fitsimage *fi)
{
 int	i;

 fprintf(fw,"array:dim=%d:",fi->dim);
 for ( i=0 ; i<fi->dim ; i++ )
  {	if ( i>0 )	fprintf(fw,"x");
	fprintf(fw,"%d",fi->naxis[i]);
  }
 fprintf(fw," bitpix=%d (%s)\n",fi->bit,fits_image_bitpix_cname(fi->bit));

 return(0);
}

static int cut_spaces(char *in,char *out)
{
 while ( isspace(*in) )			in++;
 while ( *in && ! isspace(*in) )	*out=*in,out++,in++;
 *out=0;
 return(0);
}

int fprint_fitsttable_summary(FILE *fw,fitsttable *ft,int nwl)
{
 int	i;
 char	buff[32],*format,*p;

 fprintf(fw,"rows=%d rowsize=%d ",ft->nrow,ft->rowsize);
 fprintf(fw,"fields:");
 for ( i=0 ; i<ft->ntfield ; i++ )
  {	if ( i>0 )	fprintf(fw,",");
	format=ft->tfields[i].format;
	while ( *format && isspace(*format) )	format++;
	strncpy(buff,format,31);buff[31]=0;
	p=strchr(buff,32);
	if ( p != NULL )	*p=0;
	fprintf(fw,"%s",buff);
  }
 fprintf(fw,"\n");
 return(0);
}
int fprint_fitsbtable_summary(FILE *fw,fitsbtable *fb,int nwl)
{
 int	i;
 char	type[32],unit[32],null[32],*cname;

 fprintf(fw,"rows=%d rowsize=%d ",fb->nrow,fb->rowsize);
 if ( ! nwl )
  {	fprintf(fw,"fields:");
	for ( i=0 ; i<fb->nbfield ; i++ )
	 {	if ( i>0 )	fprintf(fw,",");
		fprintf(fw,"%dx%c",fb->bfields[i].repeat,fb->bfields[i].form);
	 }
	fprintf(fw,"\n");
  }
 else
  {	fprintf(fw,"nfield=%d\n",fb->nbfield);
	for ( i=0 ; i<fb->nbfield ; i++ )
	 {	fprintf(fw,"\t\t\t\t%2dx",fb->bfields[i].repeat);
		cname=fits_bintable_form_cname(fb->bfields[i].form);
		cut_spaces(fb->bfields[i].type,type);
		cut_spaces(fb->bfields[i].unit,unit);
		cut_spaces(fb->bfields[i].null,null);
		fprintf(fw," (%s)",cname);

		if ( type[0] )
		 {	fprintf(fw,"%s",type);				}
		else
		 {	fprintf(fw,"'%c'",fb->bfields[i].form);		}

		if ( unit[0] )
		 {	fprintf(fw,"\t[%s]",unit);			}
		fprintf(fw,"\n");
	 }
  }
	 
 return(0);
}

int fprint_fits_summary(FILE *fw,fits *img,int nwl)
{
 int	is_primary,i,fcnt;

 if ( img->i.vdata != NULL )	is_primary=1;
 else				is_primary=0;

 if ( ! is_primary )
  {	fprintf(fw,"Primary      : NONE\n");
	fcnt=0;
  }
 else
  {	fprintf(fw,"Primary   [1]: IMAGE   : ");
	fprint_fitsimage_summary(fw,&img->i);
	fcnt=1;
  }
 for ( i=0 ; i<img->nxtn ; i++ )
  {	fitsheader	*fh;
	fh=fits_headerset_get_header(&img->xtns[i].header,"EXTNAME",0);
	if ( fh != NULL && fh->vtype==FITS_VSTR )
		fprintf(fw,"Extension ['%s']: ",fh->vstr);
	else
		fprintf(fw,"Extension [%d]: ",fcnt+1);
	switch ( img->xtns[i].type )
	 {   case FITS_EXT_IMAGE:
		fprintf(fw,"IMAGE   : ");
		fprint_fitsimage_summary(fw,&img->xtns[i].x.i);
		break;
	     case FITS_EXT_TABLE:
		fprintf(fw,"TABLE   : ");
		fprint_fitsttable_summary(fw,&img->xtns[i].x.t,nwl);
		break;
	     case FITS_EXT_BINTABLE:
		fprintf(fw,"BINTABLE: ");
		fprint_fitsbtable_summary(fw,&img->xtns[i].x.b,nwl);
		break;
	     default:
		fprintf(fw,"UNKNOWN\n");
		break;
	 }
	fcnt++;
  }
 return(0);
}
int main(int argc,char *argv[])
{
 fits		*img;

 FILE		*fw,*fr;
 int		i,j,sx,sy,is_help,frameno,is_primary,is_summary;
 char		*file_in,*file_out,*file_inbase,
		*file_outpnm,*file_outdump,
		*statstr,*dumpstr,*dumpformat,*pnmmstr,*file_outbg;
 int		ignore_mask,is_ppm,dump_mask;
 int		nwl;
 char		**mask;
 pnmparam	pp;
 infoimage	ii;
 imgstatparam	isp;

 progbasename=strrchr(argv[0],'/');
 if ( progbasename != NULL )	progbasename++;
 else				progbasename=argv[0];

 is_help=is_verbose=is_comment=0;

 file_outbg=file_outpnm=file_in=file_out=file_outdump=NULL;
 
 statstr=NULL;
 ii.nbox=16,ii.norder=2; 
 ii.dats=NULL,ii.ndat=0;
 
 dumpstr=NULL;nwl=0; 
 pnmmstr=NULL;ignore_mask=0;frameno=-1,is_summary=0;
 dump_mask=0;

 dumpformat="%12g";

 i=scanarg(argc,argv,SCANARG_ALLOW_FLAGS,
 	"--version:%NS-1f%q",&is_help,
        "--version-short|--short-version:%NS-2f%q",&is_help,
	"-h|--help:%f%q",&is_help,
	"--long-help|--help-long:%SN2f%q",&is_help,
	"-i|--input:%s",&file_in,
	"--frame:%d",&frameno,
	"--summary:%f",&is_summary,
	"--long-summary:%f%f",&is_summary,&nwl,
	"-o|--output:%s",&file_out,
/* I */	"-od|--output-dump:%s",&file_outdump,
	"-F|--format:%s",&dumpformat,
/* I */	"-ob|--output-background:%s",&file_outbg,
	"-m|--dump-mask:%f",&dump_mask,
/* I */	"-op|--output-pnm|--output-ppm|--output-pgm:%s",&file_outpnm,
	"--pgm:%SN0f%s",&is_ppm,&pnmmstr,
	"--ppm:%SN1f%s",&is_ppm,&pnmmstr,
	"-n|--newline|-l|--long:%f",&nwl,
	"--ignore-mask:%f",&ignore_mask,
	"-b|--box:%d",&ii.nbox,
	"-a|--order:%d",&ii.norder,
/* I */	"-d|--data:%s",&dumpstr,
	"-s|--statistics:%s",&statstr,
	"--comment:%i",&is_comment,"(C):%i",&is_comment,
	"--verbose:%i",&is_verbose,"(V):%i",&is_verbose,
	"-:%w",&file_in,
	"-*|+*:%e",
	"*:%w",&file_in,
	NULL);

 if ( i )		
  {	fprint_error("invalid command line argument near '%s'",argv[i]);
	return(1);
  }
 else if ( is_help<0 )
  {	fprint_generic_version(stdout,argv[0],"fiinfo",FI_INFO_VERSION,is_help);
	return(0);
  }
 else if ( is_help>1 )
  {	fprint_fiinfo_long_help(stdout);
	return(0);
  }
 else if ( is_help )
  {	fprint_fiinfo_usage(stdout);
	return(0);
  }

 if ( *dumpformat=='%' )	dumpformat++;
 if ( format_check(dumpformat) || strlen(dumpformat)>=30 )
  {	fprint_error("unexpected format token '%s'.\n",dumpformat);
	return(1);
  }

 if ( dumpstr != NULL )
  {	i=parse_dump_data(dumpstr,&ii.dats,&ii.ndat);
	if ( i )
	 {	fprint_error("invalid data format '%s'",dumpstr);
		return(1);
	 }
  }
 if ( statstr != NULL )
  {	double	sigma;

	isp.floydratio=0.0;
	isp.niter=0;
	isp.lower=3.0;
	isp.upper=3.0;
	isp.use_median=0;
	sigma=0.0;

	i=scanpar(statstr,SCANPAR_DEFAULT,
		"ratio:%g",&isp.floydratio,
		"iterations:%d",&isp.niter,
		"sigma:%g",&sigma,
		"lower:%g",&isp.lower,
		"upper:%g",&isp.upper,
		"median:%SN1f",&isp.use_median,
		"mean:%SN0f",&isp.use_median,
		NULL);

	if ( i )	
  	 {	fprint_error("invalid statistical parameter in '%s'.\n",statstr);
		return(1);
	 }

	if ( sigma>0.0 )
	 {	isp.lower=sigma;
		isp.upper=sigma;
	 }
   }

 pp.is_invert=0;pp.minmaxmethod=0,pp.scalemethod=0;
 pp.mmin_set=pp.mmax_set=0;
 pp.zcontrast=0.25; pp.percentage=90.0;
 pp.contrast=1.0; pp.brightness=0.5;
 pp.palette=NULL;pp.ncol=0;
 pp.is_color=is_ppm;
 pp.is_flip=pp.is_mirror=0;
 if ( pnmmstr != NULL )
  {	char	*palstr;
	int	dummy;

	palstr=NULL;
	dummy=0;

	j=scanpar(pnmmstr,SCANPAR_DEFAULT,
		"default|defaults:%f",&dummy,
		"linear:"	SNf(SCALE_LINEAR) ,&pp.scalemethod,
		"histequ:"	SNf(SCALE_HISTEQU),&pp.scalemethod,
		"log:"		SNf(SCALE_LOG)    ,&pp.scalemethod,
		"sqrt:"		SNf(SCALE_SQRT)   ,&pp.scalemethod,
		"squared:"	SNf(SCALE_SQUARED),&pp.scalemethod,
		"minmax:"	SNf(MM_MINMAX)    ,&pp.minmaxmethod,
		"zscale:"	SNf(MM_ZSCALE)    ,&pp.minmaxmethod,
		"zmin:"		SNf(MM_ZMIN)      ,&pp.minmaxmethod,
		"zmax:"		SNf(MM_ZMAX)      ,&pp.minmaxmethod,
		"percentage:"	SNf(MM_PERCENTAGE)"%g",&pp.minmaxmethod,&pp.percentage,
		"min:"		SNf(MM_MANUAL)"%f%g",&pp.minmaxmethod,&pp.mmin_set,&pp.manmin,
		"max:"		SNf(MM_MANUAL)"%f%g",&pp.minmaxmethod,&pp.mmax_set,&pp.manmax,
		"flip:%f",&pp.is_flip,
		"flop|mirror:%f",&pp.is_mirror,
		"palette:%s",&palstr,
		"zcontrast:%g",&pp.zcontrast,
		"reverse:%f",&pp.is_invert,
		"contrast:%g",&pp.contrast,
		"brightness:%g",&pp.brightness,
		NULL);

	if ( j )
	 {	fprint_error("invalid PNM specifications in '%s'",pnmmstr);
		return(1);
	 }

	if ( pp.percentage<0 )		pp.percentage=1.0;
	else if ( pp.percentage>100.0 )	pp.percentage=100.0;
	if ( pp.contrast<0.0 )		pp.contrast  =0.0;
	if ( pp.brightness<0.0 )	pp.brightness=0.0;
	else if ( pp.brightness>1.0 )	pp.brightness=1.0;

	if ( palstr != NULL && parse_palette(palstr,&pp.palette,&pp.ncol) )
	 {	fprint_error("invalid palette format");		}

  }

/* Read input image from stdin or 'file_in': */
 file_inbase=fits_basename(file_in,&frameno);
 fr=fopenread(file_inbase);
 if ( fr==NULL )
  {	fprint_error("unable to open input file '%s'.",file_inbase);
	return(1);
  }

 if ( frameno<0 || is_summary )	
	img=fits_read(fr);
 else
	img=fits_read_frame_as_extension(fr,frameno);

 fcloseread(fr);

 if ( img==NULL )
  {	fprint_error("unable to read input file as FITS data.");
	return(1);
  }

 if ( img->i.vdata != NULL )	is_primary=1;
 else				is_primary=0;
 if ( is_summary )
  {	fprint_fits_summary(stdout,img,nwl);
	return(0);
  }
 else if ( img->nxtn>=2 || (img->nxtn>=1 && is_primary) )
  {	fprintf(stderr,"Warning: FITS data contains multiple arrays or extensions, summary reported.\n");
	fprint_fits_summary(stdout,img,nwl);
	return(0);
  }

 else if ( is_primary || (img->xtns != NULL && img->xtns[0].type==FITS_EXT_IMAGE ) )
  {	fitsimage	*fi;

	if ( is_primary )
	 {	fi=&img->i;
		fits_image_get_params(&img->header,fi);
	 }
	else
	 {	fi=&img->xtns[0].x.i;
		fits_image_get_params(&img->xtns[0].header,fi);
	 }
	
	fits_image_rescale(fi);

	if ( ignore_mask )	mask=NULL;

	else if ( is_primary )
	 {	mask=fits_mask_read_from_header(&img->header,fi->sx,fi->sy,NULL);
		fits_mask_mark_nans(fi,mask,MASK_NAN);
	 }
	else
	 {	mask=fits_mask_read_from_header(&img->xtns[0].header,fi->sx,fi->sy,NULL);
		fits_mask_mark_nans(fi,mask,MASK_NAN);
	 }
	
	if ( file_out==NULL )	fw=stdout;
	else			fw=fopenwrite(file_out);
	if ( fw==NULL )	
	 {	fprint_error("unable to create output file '%s'",file_out);
		return(1);
	 }

	if ( ii.ndat>0 )
		fits_stat_v1_image(fi,mask,fw,nwl,&ii);
	else if ( statstr != NULL )
		fits_stat_v2_image(fi,mask,fw,&isp,is_comment);

	sx=fi->sx,
	sy=fi->sy;

	if ( file_outdump==NULL && dump_mask )
	 {	int	maskcount[128];
		int	i,j,k;
		char	mbuff[16];
	
		for ( k=0 ; k<128 ; k++ )
		 {	maskcount[k]=0;		}
		if ( mask==NULL )	maskcount[0]=sx*sy;
		else
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	maskcount[mask[i][j]&0x7F]++;		}
			 }
		 }
		for ( k=0 ; k<128 ; k++ )
		 {	if ( maskcount[k]<=0 )	continue;
			strmask(mbuff,k);
			fprintf(fw,"%s %9d\n",mbuff,maskcount[k]);
		 }
	 }
	
	fclosewrite(fw);

	if ( file_outdump != NULL )
	 {	char	mbuff[32];
		int	k;
	
		fw=fopenwrite(file_outdump);
		if ( fw==NULL )		
		 {	fprint_error("unable to create output file '%s'",file_outdump);
			return(1);
		 }
		for ( i=0 ; i<sy ; i++ )
		 {	for ( j=0 ; j<sx ; j++ )
			 {	fprintf(fw,"%5d %5d ",j,i);	/* LL=(0,0)! */
				sprintf(mbuff,"%%%s",dumpformat);
				fprintf(fw,mbuff,fi->data[i][j]);
				if ( dump_mask )
				 {	if ( mask != NULL )	k=mask[i][j];
					else			k=0;
					strmask(mbuff,k);
					fprintf(fw," %s\n",mbuff);
				 }
				else	fprintf(fw,"\n");
			 }
		 }
		fclosewrite(fw);
	 }

	if ( file_outpnm != NULL )
	 {	fw=fopenwrite(file_outpnm);
		if ( fw==NULL )		
		 {	fprint_error("unable to create output PNM image file '%s'",file_outpnm);
			return(1);
		 }
		fits_image_rescale(fi);
		fitsimage_dump_pnm(fi,mask,fw,&pp);
		fclosewrite(fw);
	 }
	if ( file_outbg != NULL )
	 {	fits	*bgimg;
		bgimg=fits_duplicate(img);
		create_link_background(fi,mask,&bgimg->i,0,0);
		fw=fopenwrite(file_outbg);
		if ( fw==NULL )		
		 {	fprint_error("unable to create output background file '%s'",file_outbg);
			return(1);
		 }
		fits_write(fw,bgimg);
		fclosewrite(fw);
		fits_free(bgimg);
	 }

	/* Release the mask 'mask': */
	if ( mask != NULL )	fits_mask_free(mask);
  }
 else if ( img->xtns != NULL && img->xtns[0].type==FITS_EXT_TABLE )
  {	if ( file_out != NULL )	fw=fopenwrite(file_out);
	else			fw=stdout;
	if ( fw==NULL )		
	 {	fprint_error("unable to create output file '%s'",file_out);
		return(1);
	 }
	fitsttable_dump(fw,&img->xtns[0].x.t);
	fclosewrite(fw);
  }
 else if ( img->xtns != NULL && img->xtns[0].type==FITS_EXT_BINTABLE )
  {	if ( file_out != NULL )	fw=fopenwrite(file_out);
	else			fw=stdout;
	if ( fw==NULL )		
	 {	fprint_error("unable to create output file '%s'",file_out);
		return(1);
	 }
	fitsbtable_dump(fw,&img->xtns[0].x.b);
	fclosewrite(fw);
  }

/* Release the FITS data 'img': */
 fits_free(img);

 return(0);
}
