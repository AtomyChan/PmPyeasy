/*****************************************************************************/
/* ficalib.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line tool for calibrating astronomical images.		     */
/*****************************************************************************/
#define	FI_CALIB_VERSION	"0.9"
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <fits/fits.h>

#include "fi.h"

#include "fitsmask.h"
#include "statistics.h"
#include "io/iof.h"
#include "io/scanarg.h"
#include "combine.h"
#include "tensor.h"
#include "common.h"
#include "history.h"
#include "str.h"
#include "math/spline/spline.h"
#include "math/fit/lmfit.h"
#include "fbase.h"
#include "math/splinefit.h"
#include "math/polyfit.h"
#include "io/tokenize.h"
#include "longhelp.h"

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

int fprint_ficalib_usage(FILE *fw)
{
 fprintf(fw,
"Usage:\tficalib [-C|--comment] [-V|--verbose] [-h|--help|--long-help]\n");
 fprintf(fw,
"Main calibration master images and input/output specification:\n"
"\t[-B|--input-master-bias <master-bias>]\n"
"\t[-D|--input-master-dark <master-dark>]\n"
"\t[-F|--input-master-flat <master-flat>]\n"
"\t[-i|--input <input> [<input2>...] \n"
"\t-r|--rewrite-output-name '[<from>|]<to>'\n"
"\t-o|--output <output> [<output2>...] [--[no-]clobber]\n");
 fprintf(fw,
"\t[-b|--bitpix <B>] [-M|--input-mask <mask>]\n"
"\t[--data bitpix=<bitpix>,bscale=<scale>,bzero=<zero>|<C-type>]\n");
 fprintf(fw,
"Overscan correction:\n"
"\t[--overscan area=<x1>:<y1>:<x2>:<y2>,{spline|polynomial},order=<o>,\n"
"\t\titerations=<n>,sigma=<s>]\n");
 fprintf(fw,
"Saturation levels and leaking:\n"
"\t[-s|--saturation <saturation-level>]\n"
"\t[--leak-{left-right|lower-upper|any}|--lr|--lu|--an]\n");
 fprintf(fw,
"Trimming:\n"
"\t[--image <x1>:<y1>:<x2>:<y2> [--trim|--no-trim]]\n");
 fprintf(fw,
"More fine-tunings: scaling, exposure time corrections, ...:\n"
"\t[--post-scale <mean> | --post-multiply <factor>]\n"
"\t[-g|--gain <originalgain>] [--gain-order <spatial-order>]\n"
"\t[--[no-]exptime-correction]\n");
 fprintf(fw,
"Mosaic images:\n"
"\t[--mosaic size=(<xsize>,<ysize>), \n"
"\t          ( name=<EXTNAME1>,offset=(<X1>,<Y1>),image=(<trim1>),\n"
"\t            [overscan=(<overscan1>)[,overscan=(...)]] ),\n"
"\t          ( name=<EXTNAME2>,offset=(<X2>,<Y2>),image=(<trim2>),\n"
"\t            [overscan=(<overscan2>)[,overscan=(...)]] ),\n"
"\t          ( [...] ) , ... ]\n");
 fprintf(fw,
"Other corrections:\n"
"\t[--horizontal-stripe-removal]\n");
 fprintf(fw,
"FITS header keywords (for automatic processing):\n"
"\t[--header exptime=<H>,date=<H>,time=<H>,jd=<H>,hjd=<H>,jy=<H>,hjy=<H>\n"
"\t\tgain=<H>,gainpoly=<H>,gainvmin=<H>,extname=<H>]\n"
"\t[--no-history]\n");
 fprintf(fw,
"Default headers: exptime=EXPTIME,date=DATE-OBS,time=TIME-OBS,\n"
"\t\tjd=JD,hjd=HJD,jy=JY,hjy=HJY,\n"
"\t\tgain=GAIN,gainpoly=GAINPOLY,gainvmin=GAINVMIN,extname=EXTNAME\n");
		
 return(0);
}

longhelp_entry ficalib_long_help[]=
{
 LONGHELP_OPTIONS,
 { "General options:", NULL },
 { "-h, --help",
	"Gives general summary about the command line options." },
 { "--long-help",
	"Gives a detailed list of command line options." },
 { "--version",
	"Gives some version information about the program." },
 { "-i, --input <inputs>",
	"List of input images (up to the next option)." },
 { "-o, --output <outputs>",
	"List of output images (up to the next option)." },
 { "-r, --rewrite-output-name '[<from>|]<to>'",
	"Rewrite the output file names according to the rule specified by "
	"the argument of this switch. The <to> part must contain a single "
	"asterix (*) that is replaced by the original file name or, "
	"if the <from> part is specified and it can be matched by the "
	"filename, the appropriate part of the filename is written to "
	"the output filename. For instance, '*.fits|*.calib.fits' will "
	"replace the `fits` extension to `calib.fits` if the original "
	"filename has a fits extension, otherwise, concatenates the "
	"`calib.fits` as an extension to the filename. Note that the <from> "
	"and <to> part must be separated by a vertical bar (the character "
	"used in UNIX shells for creating pipes), therefore either "
	"this character must be escaped or put into quotation marks "
	"(according to the current shell). "
	"If the -o|--output list of filenames are not, but this "
	"argument is specified, the list of original output files will be "
	"the same as the list of input filenames. " },

 { "Generic calibration options:", NULL },
 { "--overscan <overscan>",
	"Specifications for overscan corrections (see also ``Overscan "
	"correction parameters'' below)." },
 { "--mosaic size=(<sx>,<sy>),(<mosaic section 1>),...",
	"Mosaic image topology specfications (see also ``Mosaic specification "
	"parameters'' below)." },
 { "-M, --input-mask <fits>",
	"Input mask file to co-add to output image." },
 { "-B, --input-master-bias <bias>",
	"Master bias image (if any)." },
 { "-D, --input-master-dark <dark>",
	"Master dark image (if any)." },
 { "-F, --input-master-flat <flat>",
	"Master flat field image (if any)." },
 { "-s, --saturation <level>",
	"Saturation level." },
 { "--leak-left-right, --leak-lower-upper, --leak-any",
	"Pixel/readout ``blooming'' direction." },
 { "-g, --gain <gain>",
	"Nominal gain value." },
 { "--gain-order <order>",
	"Spatial polynomial order to describe gain variations." },
 { "--exptime-correction",
	"Do exposure time correction during dark subtraction." },
 { "--no-exptime-correction",
	"Don't do any exposure time correction for darks." },
 { "--no-clobber",
	"Don't overwrite existing calibrated images." },

 { "Output postprocessing and pixel depts:", NULL },
 { "--image <x1>:<y1>:<x2>:<y2>",
	"Final image area." },
 { "--trim",
	"Trim output to final image area." },
 { "--no-trim",
	"Don't trim output to final image area (default)." },
 { "--post-scale <mean>",
	"Scale output to have the specified mean value." },
 { "--post-multiply <factor>",
	"Multiply output image by the given factor." },
 { "--data <spec>",
	"Output pixel data format specification." },
 { "-b, --bitpix <bitpix>",
	"Standard FITS output bitpix value." },

 { "Overscan correction parameters:", NULL },
 { "area=<x1>:<y1>:<x2>:<y2>",
	"Overscan area." },
 { "spline",
	"Use cubic splines for overscan modelling." },
 { "polynomial",
	"Use polynomials for overscan modelling." },
 { "order=<order>",
	"Order for overscan variation (zero: constant)." },
 { "iterations=<n>",
	"Number of iterations for rejecting overscan outliers." },
 { "sigma=<s>",
	"Overscan outlier limit in units of stddev." },

 { "Other corrections:", NULL },
 { "--horizontal-stripe-removal", 
	"Remove horizontal stripes." },

 { "Mosaic specification parameters:", NULL },
 { "size=(<xsize>,<ysize>)",
	"Vertical and horizontal dimensions of the final image." },
 { "name=<extension name>",
	"FITS extenison name of the given mosaic section." },
 { "offset=(<x0>,<y0>)",
	"Offset of the given mosaic section in the final image." },
 { "image=(<x1>:<y1>:<x2>:<y2>)",
	"Trim area of the given mosaic section." },
 { "overscan=(<overscan>)",
	"Overscan parameters for the given mosaic section, in the same "
	"syntax and format as used after ``--overscan''." },

 { NULL,NULL },

};

int fprint_ficalib_long_help(FILE *fw)
{
 fprintf(fw,"Usage: ficalib [options] [-i|--input <inputs>] [-o|--output <outputs>]\n");
 fprintf(fw,
"The purpose of this program is to do various calibration steps on the input\n"
"images, including mosaic frame joining, generic overscan corrections, \n"
"bias and (exp. time scaled) dark subtraction and flat field corrections.\n\n");

 longhelp_fprint(fw,ficalib_long_help,0,-1);

 fprintf(fw,"\n");
 fprintf(fw,"Report bugs to <%s>\n",FI_MAINT_EMAIL);
		
 return(0);
}

/*****************************************************************************/

#define		OVERSCAN_NONE		0
#define		OVERSCAN_SPLINE		1
#define		OVERSCAN_POLY		2

typedef struct
 {	int	x0,y0;
	int	sx,sy; 
 } trim;

typedef struct
 {	trim	area;
	int	type,order;
	int	niter;
	double	sigma;
 } overscan;

typedef struct
 {	trim		image;
	overscan	*overscans;
	int		noverscan;
 } geometry;

typedef struct
 {	char		name[16];
	geometry	g;
	int		x0,y0;
 } section;

typedef struct
 {	int		sx,sy;
	section		*sections;
	int		nsection;
 } mosaic;

typedef struct
 {	int	sx,sy;
	char	**mask;
	char	**inmasklist;

	int	is_trim;
	int	is_exptime_correction;
	int	is_calc_gainpoly;
	int	is_combined_dark;
	double	saturation;
	double	postscale,postmultiply;
	int	reset_method;
	int	is_remove_horizontal_stripes;

	int	argc;
	char	**argv;
	char	*inputbias;
	char	*inputdark;
	char	*inputflat;

	double	gain;
	int	gainorder;
	double	*flatpoly;
	double	flatresd;
	double	flatvmin;
 } calibdata;

typedef struct
 {	char	*exptime;
	char	*date,*time;	
	char	*gain,*gainpoly,*gainvmin;
	char	*jd,*hjd,*jy,*hjy;
	char	*extname;
 } headername;

/*****************************************************************************/

void srand48_auto()
{
 int		fd;
 long int	seed;

 fd=open("/dev/urandom",O_RDONLY);
 if ( fd>=0 )
  {	(void)read(fd,&seed,sizeof(seed));
	close(fd);
  }
 else
	seed=time(NULL)*65536+getpid();

 srand48(seed);
}

int mkstemp_default_perm(char *templ)
{ 
 int	fd,c,max_attempts;
 size_t	len;
 char	*tst,*p;

 if ( templ==NULL || (len=strlen(templ))<6 )
	return(-1);
 tst=templ+len-6;
 for ( p=tst ; *p ; p++ )
  {	if ( *p != 'X' )
		return(-1);
  }

 max_attempts=1024*1024*1024;
 while ( max_attempts>0 )
  {	for ( p=tst ; *p ; p++ )
	 {	c=lrand48()%(62);
		*p=(c<10?c+'0':c<36?c-10+'a':c-36+'A');
	 }
	fd=open(templ,O_CREAT|O_RDWR|O_EXCL,S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH|S_IWOTH);
	if ( fd>=0 )
		return(fd);
	max_attempts--;
  };
 return(-1);
}

/*****************************************************************************/

int overscan_model(int sy,double *spy,double *spw,double *spm,overscan *o,fitsheaderset *header)
{
 double	f,sig,sw,**fbase,w;
 int	i,j,n,k,order,nvar,used,avail;
 double	**amatrix,*bvector,*fvars;

 switch ( o->type )
  {  case OVERSCAN_SPLINE:
	order=o->order;
 	fbase=(double **)tensor_alloc_2d(double,sy,order+1);
	fbase_spline(fbase,order,sy);
	break;
     case OVERSCAN_POLY:
	order=o->order;
 	fbase=(double **)tensor_alloc_2d(double,sy,order+1);
	fbase_polynomial(fbase,order,sy);
	break;
     default:
	order=0;
 	fbase=(double **)tensor_alloc_2d(double,sy,order+1);
	fbase_polynomial(fbase,order,sy);
	break;
  }

 nvar=order+1;
 amatrix=matrix_alloc(nvar);
 bvector=vector_alloc(nvar);
 fvars  =vector_alloc(nvar);

 avail=0;
 for ( k=0 ; k<sy ; k++ )
  {	if ( spw[k] > 0.0 )
		avail++;
  }
 used=avail;
 
 for ( n=0 ; n<=(o->niter<0?0:o->niter) ; n++ )
  {	
	for ( i=0 ; i<nvar ; i++ )
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	amatrix[i][j]=0.0;	}
		bvector[i]=0.0;
	 }

	for ( k=0 ; k<sy ; k++ )
	 {	w=spw[k];
		for ( i=0 ; i<nvar ; i++ )
		 {	fvars[i]=fbase[i][k];	}
		f=spy[k];
		for ( i=0 ; i<nvar ; i++ )
		 {	for ( j=0 ; j<nvar ; j++ )
			 {	amatrix[i][j]+=w*fvars[i]*fvars[j];	}
			bvector[i]+=w*f*fvars[i];
		 }
	 }

	solve_gauss(amatrix,bvector,nvar);

	for ( k=0 ; k<sy ; k++ )
	 {	f=0.0;
		for ( i=0 ; i<nvar ; i++ )
		 {	f+=bvector[i]*fbase[i][k];	}
		spm[k]=f;
	 }

	if ( n < o->niter )
	 {	sig=0.0;
		sw=0.0;
		for ( k=0 ; k<sy ; k++ )
		 {	f=spy[k]-spm[k];
			w=spw[k];
			sig+=w*f*f;
			sw +=w;
		 }
		sig/=sw;
		if ( sig>0.0 )	sig=sqrt(sig);
		else		sig=0.0;
		for ( k=0 ; k<sy ; k++ )
		 {	f=spy[k]-spm[k];
			if ( fabs(f) >= o->sigma * sig )
			 {	spw[k]=0.0;
				used--;
			 }
		 }
	 }
  }

 sig=0.0;
 sw=0.0;
 for ( k=0 ; k<sy ; k++ )
  {	f=spy[k]-spm[k];
	w=spw[k];
	sig+=w*f*f;
	sw +=w;
  }
 sig/=sw;
 if ( sig>0.0 )	sig=sqrt(sig);
 else		sig=0.0;

 /* just export the overscan information */
 if ( header != NULL )
  {	char	buff[80],*type;
	int	l;
	switch ( o->type )
	 {   case OVERSCAN_SPLINE:
		type="spline";
		break;
	     case OVERSCAN_POLY:
		type="polynomial";
		break;
	     default:
		type="const";
		break;
	 }
	l=snprintf(buff,80,"ficalib: %s o=%d c=",type,o->order);
	for ( i=0 ; i<nvar ; i++ )
	 {	if ( i )
			l+=snprintf(buff+l,(80-l>0?80-l:0),",");
		l+=snprintf(buff+l,(80-l>0?80-l:0),"%.1f",bvector[i]);
	 }
	l+=snprintf(buff+l,(80-l>0?80-l:0)," s=%.1f n=%d/%d/%d",sig,sy,sy-avail,sy-used);

	fits_headerset_set_string(header,"OVERSCAN",FITS_SH_ADD,buff,NULL);
  }

 vector_free(fvars);
 vector_free(bvector);
 matrix_free(amatrix);

 tensor_free(fbase);

 return(0);
}

int overscan_do_vertical(fitsimage *img,char **mask,int x0,int y0,int sx,int sy,int ox,int wx,overscan *o,compar *cp,fitsheaderset *header)
{
 int	i,j,lp,np,x,y;
 double	*line,*spw,*spy,*spm;

 spw =(double *)malloc(sizeof(double)*sy);
 spy =(double *)malloc(sizeof(double)*sy);
 line=(double *)malloc(sizeof(double)*wx);

 np=0;
 for ( i=0 ; i<sy ; i++ )
  {	lp=0;
	for ( j=0 ; j<wx ; j++ )
	 {	y=y0+i;
		x=ox+j;
		if ( x<0 || x>=img->sx || y<0 || y>=img->sy )
			continue;
		else if ( mask != NULL && mask[y][x] ) 
			continue;
		line[lp]=img->data[y][x];
		lp++;
	 }
	if ( lp>0 )
	 {	spw[i]=1.0;
		spy[i]=combine_points(line,lp,cp);
		np++;
	 }
	else
	 {	spw[i]=0.0;
		spy[i]=0.0;
	 }
  }
 
 if ( np < o->order+1 )	/* overscan failed due to insufficient num of points */
  {	free(line);
	free(spy);
	free(spw);
	return(-1);	
  }

 spm  =(double *)malloc(sizeof(double)*sy);

 overscan_model(sy,spy,spw,spm,o,header);

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	img->data[y0+i][x0+j] -= spm[i];	}
  }

 free(spm);

 free(line);
 free(spy);
 free(spw);
 
 return(0);
}

int overscan_do_horizontal(fitsimage *img,char **mask,int x0,int y0,int sx,int sy,int oy,int wy,overscan *o,compar *cp,fitsheaderset *header)
{
 int	i,j,lp,np,x,y;
 double	*line,*spw,*spy,*spm;

 spw =(double *)malloc(sizeof(double)*sx);
 spy =(double *)malloc(sizeof(double)*sx);
 line=(double *)malloc(sizeof(double)*wy);

 np=0;
 for ( i=0 ; i<sx ; i++ )
  {	lp=0;
	for ( j=0 ; j<wy ; j++ )
	 {	x=x0+i;
		y=oy+j;
		if ( x<0 || x>=img->sx || y<0 || y>=img->sy )
			continue;
		else if ( mask != NULL && mask[y][x] ) 
			continue;
		line[lp]=img->data[y][x];
	 }
	if ( lp>0 )
	 {	spw[i]=1.0;
		spy[i]=combine_points(line,lp,cp);
		np++;
	 }
	else
	 {	spw[i]=0.0;
		spy[i]=0.0;
	 }
  }

 if ( np < o->order+1 )	/* overscan failed due to insufficient num of points */
  {	free(line);
	free(spy);
	free(spw);
	return(-1);	
  }

 spm  =(double *)malloc(sizeof(double)*sx);

 overscan_model(sx,spy,spw,spm,o,header);

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	img->data[y0+i][x0+j] -= spm[j];	}
  }

 free(spm);

 free(line);
 free(spy);
 free(spw);
 
 return(0);
}

int overscan_check(trim *area,trim *o)
{
 if ( o->y0==area->y0 && o->sy==area->sy && 
 abs(2*(o->x0-area->x0)+(o->sx-area->sx))>=o->sx+area->sx ) 
	return(+1);	/* vertical overscan	*/
 else if ( o->x0==area->x0 && o->sx==area->sx && 
 abs(2*(o->y0-area->y0)+(o->sy-area->sy))>=o->sy+area->sy )
	return(-1);	/* horizontal overscan	*/
 else
	return(0);	/* invalid overscan	*/
}

int overscan_correction(fitsimage *img,char **mask,trim *area,overscan *scans,int nscan,compar *cp,fitsheaderset *header)
{
 int		ret,r,i;
 overscan	*o;

 ret=0;
 for ( i=0 ; i<nscan ; i++ )
  {	o=&scans[i];
	if ( o->area.y0==area->y0 && o->area.sy==area->sy && 
	abs(2*(o->area.x0-area->x0)+(o->area.sx-area->sx))>=o->area.sx+area->sx )
	 {	r=overscan_do_vertical(img,mask,area->x0,area->y0,area->sx,area->sy,o->area.x0,o->area.sx,o,cp,header);
		if ( ! r )	ret++;
	 }
	else if ( o->area.x0==area->x0 && o->area.sx==area->sx && 
	abs(2*(o->area.y0-area->y0)+(o->area.sy-area->sy))>=o->area.sy+area->sy )
	 {	r=overscan_do_horizontal(img,mask,area->x0,area->y0,area->sx,area->sy,o->area.y0,o->area.sy,o,cp,header);
		if ( ! r )	ret++;
	 }

  }

 if ( header != NULL )
  {	int	l;
	char	buff[80];
	l=snprintf(buff,80,"ficalib: total_overscans=%d/%d",ret,nscan);
	fits_headerset_set_string(header,"OVERSCNS",FITS_SH_ADD,buff,NULL);
  }
 
 return(ret);
}

/*****************************************************************************/

/*
fits *combine_more_images(char **files,int nfile,int frameno,
	compar *cp,calibdata *cd,presubdata *presubs,int npresub)
{
 FILE	*fr;
 fits	*img,*ret;
 comimg	*inputs;
 int	i,j,sx,sy;
 char	**outmask,*basename;

 inputs=(comimg *)malloc(sizeof(comimg)*nfile);
 for ( i=0 ; i<nfile ; i++ )
  {	basename=fits_basename(files[i],&frameno);
	if ( (fr=fopenread(basename))==NULL )
	 {	fprint_error("unable to open input file '%s'.",basename);
		return(NULL);
	 }
	img=fits_seek_frame_to_image(fr,frameno);
	if ( img==NULL )
	 {	fprint_error("unable to interpret input as FITS data.");
		return(NULL);
	 }
	if ( img->i.dim != 2 )
	 {	fprint_error("image dimension differs from 2.");
		return(NULL);
	 }
	inputs[i].img=img;
	inputs[i].fr=fr;
  }

 if ( cd->sx<=0 || cd->sy<=0 )	cd->sx=sx=inputs[0].img->i.sx,
				cd->sy=sy=inputs[0].img->i.sy;
 else				sx=cd->sx,sy=cd->sy;

 for ( i=0,j=0 ; i<nfile && ! j ; i++ )
  {     if ( inputs[i].img->i.sx != sx || inputs[i].img->i.sy != sy )
		j=1;
  }
 if ( j )
  {	fprint_error("size of input images differ");
	exit(1);
  }

 if ( cd->mask==NULL )
  {	cd->mask=fits_mask_create_empty(sx,sy);
	if ( cd->inmasklist != NULL )
	 {     if ( join_masks_from_files(cd->mask,sx,sy,cd->inmasklist) )
		 {	fprint_error("unable to read one of the input mask files");
			exit(1);
		 }	
	 }
  }

 outmask=fits_mask_create_empty(sx,sy);
 ret=fits_create();
 fits_copy_full_header(ret,inputs[0].img);
 fits_alloc_image(ret,sx,sy);
 fits_reset_image(ret);
 ret->i.curr.bscale=1.0;
 ret->i.curr.bzero=0.0;
 ret->i.bit=-32;

 combine_images_from_files(inputs,nfile,ret,cp,cd->mask,outmask,presubs,npresub,0);

 combine_cleanup(inputs,nfile);
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	cd->mask[i][j] |= outmask[i][j];	}
  }
 fits_mask_free(outmask);
 
 free(inputs);

 return(ret);
}
*/

/*****************************************************************************/

fitsextension *calibrate_mosaic_get_image_by_extname(fits *img,char *name,headername *hn)
{
 fitsextension	*fx;
 int		i,l;
 fitsheader	*hd;

 for ( i=0 ; i<img->nxtn ; i++ )
  {	if ( img->xtns[i].type != FITS_EXT_IMAGE )
		continue;
	hd=fits_headerset_get_header(&img->xtns[i].header,hn->extname,0);
	if ( hd==NULL || hd->vtype != FITS_VSTR )
		continue;
	l=strlen(name);
	if ( l>0 && memcmp(name,hd->vstr,l)==0 && 0<=hd->vstr[l] && hd->vstr[l]<=32 )
	 {	fx=&img->xtns[i];
		return(fx);
	 }
  }
 return(NULL);		
}

int calibrate_mosaic_check_sections(fits *img,mosaic *m,headername *hn)
{
 int		i,nvalid,bitpix;
 fitsextension	*fx;

 if ( img==NULL || img->xtns==NULL || img->nxtn<=0 || m==NULL || m->sections==NULL )
	return(0);

 bitpix=0;
 nvalid=0;
 for ( i=0 ; i<m->nsection ; i++ )
  {	fx=calibrate_mosaic_get_image_by_extname(img,m->sections[i].name,hn);
	if ( fx != NULL && fx->x.i.dim==2 && fx->x.i.sx>0 && fx->x.i.sy>0 )
	 {	nvalid++;
		if ( ! bitpix )	
			bitpix=fx->x.i.bit;
		else if ( bitpix>0 && fx->x.i.bit<0 )
			bitpix=fx->x.i.bit;
		else if ( bitpix>0 && fx->x.i.bit>bitpix )
			bitpix=fx->x.i.bit;
		else if ( bitpix<0 && fx->x.i.bit<bitpix )
			bitpix=fx->x.i.bit;
	 }
  }

 if ( nvalid==m->nsection )
	return(bitpix);
 else
	return(0);
	
}

double calibrate_get_flat_mean(fitsimage *flat,char **mask)
{
 int	i,j,sx,sy;
 double	s0,s1;

 sx=flat->sx;
 sy=flat->sy;

 s0=s1=0.0;
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( mask != NULL && mask[i][j] )
			continue;
		s1+=flat->data[i][j];
		s0+=1.0;
	 }
  }

 if ( s0 <= 0.0 )
	return(-1.0);
 else
	return(s1/s0);
}

int calibrate_image(char *input,char *output,
	geometry *g,mosaic *m,calibdata *cd,compar *cpover,fitsdataparam *fdp,headername *hn,
	fitsimage *bias,fitsimage *dark,fitsimage *flat,
	double flatmean,double darktime,int no_history)
{
 FILE		*fr,*fw;
 fits		*img;
 char		*tmpoutput;
 int		sx,sy,i,j,bitpix;
 int		had_bias,had_dark,had_flat;
 double		dark_exp_scale,pscale;
 char		**mask;
 double		cexptime;
 double		gain;
 int		is_trim,is_mark_saturated,is_do_overscan;

 fr=fopenread(input);
 if ( fr==NULL )
  {	fprint_error("unable to open file '%s' for calibration, file skipped",input);
	return(-1);
  }
 img=fits_read(fr);
 fcloseread(fr);
 if ( img==NULL )
  {	fprint_error("unable to parse file '%s' as a FITS image, skipped",input);
	return(-1);
  }

 /* we have a mosaic image: */
 if ( (bitpix=calibrate_mosaic_check_sections(img,m,hn)) != 0 )
  {	fits		*out;
	int		i,k,l,x,y,ix,iy;
	fitsextension	*fx;
	fitsimage	*fi;
	char		**tmask;
	trim		*t;

	out=fits_create();
	/* fits_set_standard(out,NULL); */
	fits_headerset_copy(&out->header,&img->header);
	fits_alloc_image(out,m->sx,m->sy);
	out->i.bit=bitpix;
	out->i.curr.bscale=1.0;
	out->i.curr.bzero=0.0;
	if ( out->i.bit<0 )	out->i.read.bscale=1.0,out->i.read.bzero=0.0;
	else			out->i.read.bscale=1.0,out->i.read.bzero=(1<<(out->i.bit-1));
	fits_set_image_params(out);
	/* fits_headerset_merge(&out->header,&img->header,1); */ /* TBD */

	mask=fits_mask_create_empty(m->sx,m->sy);
	for ( k=0 ; k<m->sy ; k++ )
	 {	memset(mask[k],MASK_OUTER,m->sx);		}

	for ( i=0 ; i<m->nsection ; i++ )
	 {	fx=calibrate_mosaic_get_image_by_extname(img,m->sections[i].name,hn);
		if ( fx==NULL )		/* should not happen */
			continue;

		fi=&fx->x.i;

		tmask=fits_mask_read_from_header(&fx->header,fi->sx,fi->sy,NULL);
		fits_image_rescale(fi);

		t=&m->sections[i].g.image;

		if ( cd->saturation > 0.0 )
			mark_saturated_pixels(fi,tmask,NULL,cd->saturation,cd->reset_method);

		 /* do the real overscan correction(s) on the image: */
		if ( m->sections[i].g.noverscan>0 )
			overscan_correction(fi,tmask,t,m->sections[i].g.overscans,m->sections[i].g.noverscan,cpover,&fx->header);

		/* fprintf(stderr,"[3] m=(%d,%d) fi=(%d,%d) size=(%d,%d) offset=(%d,%d)\n",m->sx,m->sy,fi->sx,fi->sy,t->sx,t->sy,t->x0,t->y0); */

		for ( k=0 ; k<t->sy ; k++ )
		 {	y=k+t->y0;
			if ( ! ( 0<=y && y<fi->sy ) )
				continue;
			iy=m->sections[i].y0+k;
			if ( ! ( 0<=iy && iy<m->sy ) )
				continue;
			for ( l=0 ; l<t->sx ; l++ )
			 {	x=l+t->x0;
				if ( ! ( 0<=x && x<fi->sx ) )
					continue;
				ix=m->sections[i].x0+l;
				if ( ! ( 0<=ix && ix<m->sx ) )
					continue;
				out->i.data[iy][ix] = fi->data[y][x];
				mask[iy][ix] = tmask[y][x];
			 }
		 }
	
		fits_mask_free(tmask);
	 }

	is_trim=0;		/* handling mosaics was far more complicated ;) */
	is_mark_saturated=0;	/* already done */
	is_do_overscan=0;	/* already done */

	fits_free(img);		/* drop the original image... */
	img=out;		/* ... and use the freshly prepared one */
  }

 /* we have a simple image, hopefully: */
 else if ( img->i.dim != 2 )
  {	fprint_error("file '%s' is not a 2 dimensional FITS image, skipped",input);
	fits_free(img);
	return(-1);
  }
 else
  {
	fits_rescale(img);

	/* extract additional information from the header... */
	mask=fits_mask_read_from_header(&img->header,img->i.sx,img->i.sy,NULL);

	is_trim=cd->is_trim;
	is_mark_saturated=1;
	is_do_overscan=1;
  }

 if ( cd->gain > 0.0 )
	gain=cd->gain;
 else if ( fits_headerset_get_as_double(&img->header,hn->gain,&gain,0) )
	gain=1.0;

 if ( fits_headerset_get_as_double(&img->header,hn->exptime,&cexptime,0) )
 	cexptime=-1.0;

 if ( (! cd->is_exptime_correction) && cexptime>0.0 && darktime>0.0 &&
 (cd->is_combined_dark || bias==NULL) && cexptime != darktime )
 	fprint_warning("image '%s' cannot be properly calibrated due to the lack of bias and/or the explicit usage of --no-exptime-correction",input);
 	
 /* ... done: beyond this point we do not need the header img->header */

 /*   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   */

 /* start of the real calibration (up to end, this part would be preserved   */
 /* in the cases when ficalib would handle FITS with multiple extensions)    */

 /* mark saturated pixels: */
 if ( cd->saturation > 0.0 && is_mark_saturated )
	mark_saturated_pixels(&img->i,mask,NULL,cd->saturation,cd->reset_method);

 /* do the real overscan correction(s) on the image: */
 if ( g->noverscan>0 && is_do_overscan )
 	overscan_correction(&img->i,mask,&g->image,g->overscans,g->noverscan,cpover,&img->header);

 /* trim after the overscan correction (optional): */
 if ( is_trim )
  {	/* this is the only strange part: if we have multiple frames/ext's   */
 	/* in the future, these lines should operate on the extension array  */
	/* instead on the primary image of the general FITS container:	     */

	sx=g->image.sx;
	sy=g->image.sy;

	fits_mask_trim(&mask,img->i.sx,img->i.sy,g->image.x0,g->image.y0,sx,sy,MASK_OUTER);
	fits_image_trim(&img->i,g->image.x0,g->image.y0,sx,sy,0.0);
   }
 else
  {	sx=img->i.sx;
	sy=img->i.sy;
  }

 if ( join_masks_from_files(mask,sx,sy,cd->inmasklist) )
  {	fprint_error("unable to read one of the input mask files");
	exit(1); 
  }

 /* adjoin masks: cd->mask comes from the master frames */
 if ( cd->mask != NULL )
  {	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
			mask[i][j] |= cd->mask[i][j];
	 }
  }

 /* do the bias subtraction (really simple): */
 if ( bias != NULL && bias->sx==sx && bias->sy==sy )
  {	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	img->i.data[i][j] -= bias->data[i][j];		}
	 }
	had_bias=1;
  }
 else
	had_bias=0;

 /* do the dark frame subtraction (optionally scaled by the exposure times): */
 if ( dark != NULL && dark->sx==sx && dark->sy==sy )
  {	double	dfactor;
	if ( ! cd->is_exptime_correction )
		dfactor=1.0;
	else if ( cexptime>0.0 && darktime>0.0 )
		dfactor=cexptime/darktime;
	else
		dfactor=1.0;
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	img->i.data[i][j] -= dfactor*dark->data[i][j];	}
	 }
	had_dark=1;
	dark_exp_scale=dfactor;
  }
 else
  {	had_dark=0;
	dark_exp_scale=0.0;
  }

 /* do the flat correction, mark invalid pixels as 'outer' */
 if ( flat != NULL && flat->sx==sx && flat->sy==sy && flatmean>0.0 )
  {	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	if ( flat->data[i][j] > 0.0 )
				img->i.data[i][j] *= (flatmean/flat->data[i][j]);
			else
			 {	img->i.data[i][j] = 0.0;
				mask[i][j] |= MASK_OUTER;  /* or something else? */
			 }
		 }
	 }
	had_flat=1;
  }
 else
	had_flat=0;

 if ( cd->postmultiply != 0.0 )
  {	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	img->i.data[i][j] *= cd->postmultiply;		}
	 }
	pscale=cd->postmultiply;
  }
 else if ( cd->postscale != 0.0 )
  {	double	s0,s1,m;
	s0=s1=0.0;
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	if ( mask[i][j] )
				continue;
			s0+=1.0;
			s1+=img->i.data[i][j];
		 }
	 }
	m=cd->postscale*s0/s1;
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	img->i.data[i][j] *= m;		}
	 }
	pscale=m;
  }
 else
	pscale=0.0;

 if ( cd->is_remove_horizontal_stripes )
  {	double	*arr,avg;
	arr=(double *)malloc(sizeof(double)*sx);
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	arr[j]=img->i.data[i][j];	 }
		avg=median(arr,sx);
		for ( j=0 ; j<sx ; j++ )
		 {	img->i.data[i][j] -= avg;	 }
	 }
	free(arr);
	arr=NULL;
  }
	
 /* the end of the real calibration per standalone image */

 /*   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   */

 if ( output==NULL || strcmp(output,"-")==0 || strcmp(output,"/dev/stdout")==0 )
	tmpoutput=NULL;
 else
  {	tmpoutput=NULL;
	strappendf(&tmpoutput,"%s.XXXXXX",output);
  }

 if ( tmpoutput==NULL )
	fw=stdout;
 else
  {	int	fh;
	fh=mkstemp_default_perm(tmpoutput);
	fw=fdopen(fh,"wb");
  }

 if ( fw==NULL )
  {	fprint_error("unable to create file '%s'",output);
	if ( tmpoutput != NULL )	free(tmpoutput);
	fits_mask_free(mask);
	fits_free(img);
	return(-1);
  }

 if ( is_trim )
  {	char		buff[80];
	fitsheaderset	*header;
	double		crpix1,crpix2;

	header=&img->header;
	
	sprintf(buff,"[%d:%d,%d:%d]",
		g->image.x0+1,g->image.x0+g->image.sx,
		g->image.y0+1,g->image.y0+g->image.sy);
	fits_headerset_set_integer(header,"WCSDIM",FITS_SH_ADD,1,NULL);
	fits_headerset_set_double(header,"LTM1_1",FITS_SH_ADD,1.0,NULL);
	fits_headerset_set_double(header,"LTM2_2",FITS_SH_ADD,1.0,NULL);
	fits_headerset_set_string(header,"WAT0_001",FITS_SH_ADD,"system=physical",NULL);
	fits_headerset_set_string(header,"WAT1_001",FITS_SH_ADD,"wtype=linear",NULL);
	fits_headerset_set_string(header,"WAT2_001",FITS_SH_ADD,"wtype=linear",NULL);
	fits_headerset_set_string(header,"CCDSEC",FITS_SH_ADD,buff,NULL);
	fits_headerset_set_double(header,"LTV1",FITS_SH_ADD,-(double)g->image.x0,NULL);
	fits_headerset_set_double(header,"LTV2",FITS_SH_ADD,-(double)g->image.y0,NULL);
	fits_headerset_delete_all(header,"DATASEC");
	fits_headerset_delete_all(header,"BIASSEC");
	fits_headerset_delete_all(header,"TRIMSEC");

	if ( (! fits_headerset_get_as_double(header,"CRPIX1",&crpix1,0)) &&
	     (! fits_headerset_get_as_double(header,"CRPIX1",&crpix2,0)) )
	 {	crpix1-=g->image.x0;
		crpix2-=g->image.y0;
		fits_headerset_set_double(header,"CRPIX1",FITS_SH_FIRST,crpix1,NULL);
		fits_headerset_set_double(header,"CRPIX2",FITS_SH_FIRST,crpix2,NULL);
	 }
  }

 if ( had_bias || cd->is_combined_dark )
  {	fitsheaderset	*header;
	char		buff[80];
	header=&img->header;
	snprintf(buff,80,"ficalib: zero: bias=\"%s\"",cd->inputbias);
	fits_headerset_set_string(header,"ZEROCOR",FITS_SH_ADD,buff,NULL);
  }

 if ( had_dark )
  {	fitsheaderset	*header;
	char		buff[80];
	header=&img->header;
	snprintf(buff,80,"ficalib: dark: expscale=%.3f dark=\"%s\"",dark_exp_scale,cd->inputdark);
	fits_headerset_set_string(header,"DARKCOR",FITS_SH_ADD,buff,NULL);
  }

 if ( had_flat )
  {	fitsheaderset	*header;
	char		buff[80],bcmm[80];
	header=&img->header;
	snprintf(buff,80,"ficalib: flat: mean=%.1f flat=\"%s\"",flatmean,cd->inputflat);
	fits_headerset_set_string(header,"FLATCOR",FITS_SH_ADD,buff,NULL);
	if ( cd->flatpoly != NULL )
	 {	int	i,n,nvar;
		nvar=(cd->gainorder+1)*(cd->gainorder+2)/2;
		n=0;
		for ( i=0 ; i<nvar ; i++ )
		 {	n+=snprintf(buff+n,80,(i>0?",%.2f":"%.2f"),
			gain*cd->flatpoly[i]/flatmean);
		 }
		sprintf(bcmm,"[%.2f]",gain*cd->flatresd/flatmean);
		fits_headerset_set_string(header,hn->gainpoly,FITS_SH_ADD,buff,bcmm);
		fits_headerset_set_double(header,hn->gainvmin,FITS_SH_ADD,gain*cd->flatvmin/flatmean,NULL);
	 }

  }

 if ( pscale != 0.0 )
  {	fitsheaderset	*header;
	char		buff[80];
	header=&img->header;
	snprintf(buff,80,"ficalib: post-multiply: factor=%.4f",pscale);
	fits_headerset_set_string(header,"PSTSCALE",FITS_SH_ADD,buff,NULL);
	fits_headerset_set_double(header,"MULTIPLD",FITS_SH_ADD,pscale,"ficalib: post-multiply");
  }

 if ( ! no_history )
	fits_history_export_command_line(img,"ficalib",FI_CALIB_VERSION,cd->argc,cd->argv);

 if ( fdp->bitpix )
  {     img->i.bit=fdp->bitpix;               }
 if ( fdp->is_scale )
  {     img->i.read.bscale=fdp->bscale,
        img->i.read.bzero=fdp->bzero;
  }
 else if ( img->i.bit<0 )
  {     img->i.read.bscale=1.0,
        img->i.read.bzero=0.0;
  }

 fits_set_image_params(img);

 fits_backscale(img,img->i.read.bscale,img->i.read.bzero); 
 mark_integerlimited_pixels(&img->i,mask,img->i.bit,1,MASK_OVERSATURATED,MASK_OVERSATURATED);
 fits_mask_export_as_header(&img->header,1,mask,sx,sy,NULL);

 fits_write(fw,img);

 fits_mask_free(mask);
 fits_free(img);

 fclosewrite(fw);
 if ( tmpoutput != NULL )
  {	int	r;
	if ( (r=link(tmpoutput,output)) && errno==EEXIST )
	 {	unlink(output);
		r=link(tmpoutput,output);
	 }
	if ( ! r )
		unlink(tmpoutput);
	free(tmpoutput);
	tmpoutput=NULL;
  }

 return(0);
}


/*****************************************************************************/

int calibrate_fit_flatpoly(double **data,int sx,int sy,char **inmask,
int order,double *flatpoly,double *flatresd,double *flatvmin)
{
 int	i,j,k,l,nvar;
 double	x,y,s2,s0,sig,z,dz;
 double	**amatrix,*bvector,*fvars,vmin;
 int	niter,iiter;
 double	sigma;
 char	**mask;

 if ( data==NULL || sx<=0 || sy<=0 || flatpoly==NULL || order<0 )
	return(-1);

 nvar=(order+1)*(order+2)/2;

 amatrix=matrix_alloc(nvar);
 bvector=vector_alloc(nvar);
 fvars  =vector_alloc(nvar);

 niter=2;
 sigma=2.0;

 mask=fits_mask_create_empty(sx,sy);

 if ( inmask != NULL )
  {	for ( i=0 ; i<sy ; i++ )
	 {	memcpy(mask[i],inmask[i],sx);		}
  }

 for ( iiter=0 ; iiter<=niter ; iiter++ )
  {	for ( k=0 ; k<nvar ; k++ )
	 {	for ( l=0 ; l<nvar ; l++ )
		 {	amatrix[k][l]=0.0;		}
		bvector[k]=0.0;
	 }

	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	if ( mask[i][j] )
				continue;
			x=(double)j+0.5;
			y=(double)i+0.5;
			z=data[i][j];

			eval_2d_monoms(x,y,order,fvars,0.5*(double)sx,0.5*(double)sy,0.5*(double)sx);

			for ( k=0 ; k<nvar ; k++ )
			 {	for ( l=0 ; l<nvar ; l++ )
				 {	amatrix[k][l] += fvars[k]*fvars[l];	}
				bvector[k] += fvars[k]*z;
			 }
		 }
	 }
	
	solve_gauss(amatrix,bvector,nvar);
	for ( k=0 ; k<nvar ; k++ )
	 {	flatpoly[k]=bvector[k];		}
	
	s2=s0=0.0;
	vmin=0.0;
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	if ( mask[i][j] )
				continue;
			x=(double)j+0.5;
			y=(double)i+0.5;
	
			z=eval_2d_poly(x,y,order,flatpoly,0.5*(double)sx,0.5*(double)sy,0.5*(double)sx);
			dz=data[i][j]-z;
			s2+=dz*dz;
			if ( s0<=0.0 || data[i][j]<=vmin )
				vmin=data[i][j];
			s0+=1.0;
		 }
	 }
	if ( s0>0 && s2>0 )
		sig=sqrt(s2/s0);
	else
		sig=0.0;

	if ( flatresd != NULL )
		*flatresd=sig;
	if ( flatvmin != NULL )
		*flatvmin=vmin;

	if ( iiter >= niter )
		continue;

	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	if ( mask[i][j] )
				continue;
			x=(double)j+0.5;
			y=(double)i+0.5;
	
			z=eval_2d_poly(x,y,order,flatpoly,0.5*(double)sx,0.5*(double)sy,0.5*(double)sx);
			dz=data[i][j]-z;
			if ( ! ( -sig*sigma <= dz && dz <= +sig*sigma ) )
				mask[i][j] |= 0x80;
		 }
	 }

  }

 fits_mask_free(mask);

 vector_free(fvars);
 vector_free(bvector);
 matrix_free(amatrix);

 return(0);
}

static int count_list(char **list)
{
 int	i;
 if ( list==NULL )	return(0);
 for ( i=0 ; list[i] != NULL ; ) i++;
 return(i);
}

static int match_filename(char *in,char *pattern,char **ret)
{
 char	*p;

 while ( (*in)==(*pattern) && *pattern )
	in++,pattern++;

 if ( *pattern=='*' )
	pattern++;
 else if ( *pattern || *in )
	return(0);
 else
  {	if ( ret != NULL ) *ret=NULL;
	return(1);
  }

 if ( ! (*pattern) )
  {	if ( ret != NULL )
		*ret=strdup(in);
	return(1);
  }

 p=in;
 while ( *p )
  {	if ( strcmp(p,pattern)==0 )
	 {	if ( ret != NULL )
		 {	*ret=(char *)malloc((p-in)+1);
			memcpy(*ret,in,(p-in));
			(*ret)[p-in]=0;
		 }
		return(1);
	 }
	p++;
  };

 return(0);
}

char *rewrite_filename(char *in,char *from,char *to)
{
 char	*intoken;

 if ( strchr(from,'*')==NULL || strchr(from,'*') != strrchr(from,'*') )
	return(NULL);
 else if ( strchr(to,'*') != NULL && strchr(to,'*') != strrchr(to,'*') )
	return(NULL);

 if ( ! ( match_filename(in,from,&intoken) && intoken != NULL ) )
	intoken=in;

 if ( strchr(to,'*')==NULL )
  {	char	*ret;
	ret=NULL;
	strappendf(&ret,"%s%s",intoken,to);
	if ( intoken != in ) free(intoken);
	return(ret);
  }
 else
  {	char	*p1,*p2,*ret;
	p1=strdup(to);
	p2=strchr(p1,'*');
	*p2=0,p2++;
	ret=NULL;
	strappendf(&ret,"%s%s%s",p1,intoken,p2);
	free(p1);
	if ( intoken != in ) free(intoken);
	return(ret);
  }

}

/*****************************************************************************/

int ficalib_parse_imagetrim_data(char *imgarea,trim *image)
{
 int	t;
	
 t=sscanf(imgarea,"%d:%d:%d:%d",&image->x0,&image->y0,&image->sx,&image->sy);
 if ( t<4 )
	return(1);

 if ( image->sx<image->x0 )
  {	t=image->sx,image->sx=image->x0,image->x0=t;	}
 if ( image->sy<image->y0 )
  {	t=image->sy,image->sy=image->y0,image->y0=t;	}
 image->sx=image->sx-image->x0+1;
 image->sy=image->sy-image->y0+1;

 return(0);
}

int ficalib_parse_overscan_data(char *overscanspec,overscan *o)
{
 int	i,t,x1,y1,x2,y2;
 int	niter,order,type;
 double	sigma;
 char	*area;

 x1=y1=x2=y2=0;
 niter=0;
 sigma=3.0;
 order=2;
 type=OVERSCAN_SPLINE;
 area=NULL;

 i=scanpar(overscanspec,SCANPAR_DEFAULT,
	"area:%s",&area,
	"iterations:%d",&niter,
	"sigma:%g",&sigma,
	"order:%d",&order,
	"polynomial:" SNf(OVERSCAN_POLY) , &type,
	"spline:" SNf(OVERSCAN_SPLINE)   , &type,
	NULL);

 if ( area != NULL )
  {	if ( sscanf(area,"%d:%d:%d:%d",&x1,&y1,&x2,&y2)<4 )
	i=1;
  }	
 if ( i )
	return(1);

 if ( x2<x1 )	t=x1,x1=x2,x2=t;
 if ( y2<y1 )	t=y1,y1=y2,y2=t;

 o->area.x0=x1;
 o->area.sx=x2-x1+1;
 o->area.y0=y1;
 o->area.sy=y2-y1+1;
 o->type=type;
 o->order=order;
 o->niter=niter;
 o->sigma=sigma;

 return(0);
}

int ficalib_parse_mosaic_data(char *mstr,mosaic *m)
{
 char		*wms,**mtokens,*wss,**stokens,*cmd[4];
 int		i,j,k,l,n;
 section	*s;

 wms=strdup(mstr);
 remove_spaces_and_comments(wms);
 for ( i=0 ; wms[i] ; i++ )
  {	if ( wms[i]=='(' )	wms[i]='[';
	else if ( wms[i]==')' )	wms[i]=']';
  }

 m->nsection=0;
 m->sections=NULL;

 for ( i=0,l=0 ; wms[i] ; i++ )
  {	if ( wms[i]=='[' )	l++;
	else if ( wms[i]==']' )	l--;
	else if ( wms[i]==',' || wms[i]==';' )
	 {	if ( l>0 )	wms[i]=',';
		else		wms[i]=';';
	 }
  }
 if ( l != 0 )
  {	free(wms);
	return(1);
  }
 mtokens=tokenize_char_dyn(wms,';'); 
 if ( mtokens==NULL )
  {	free(wms);
	return(1);
  }
 else if ( mtokens[0]==NULL )
  {	free(mtokens);
	free(wms);
	return(1);
  }

 for ( n=0 ; mtokens[n] != NULL ; n++ )
  {	wss=mtokens[n];
	l=strlen(wss);
	if ( l>2 && (wss[0]=='[') && (wss[l-1]==']') )
	 {	int		x,y;
		overscan	o;

		wss[l-1]=0;
		wss++;

		m->sections=(section *)realloc(m->sections,sizeof(section)*(m->nsection+1));
		s=&m->sections[m->nsection];
		m->nsection++;
		s->x0=0;
		s->y0=0;
		s->g.image.x0=0;
		s->g.image.y0=0;
		s->g.image.sx=0;
		s->g.image.sy=0;
		s->g.overscans=NULL;
		s->g.noverscan=0;

		for ( i=0,l=0 ; wss[i] ; i++ )
		 {	if ( wss[i]=='[' )	l++;
			else if ( wss[i]==']' )	l--;
			else if ( wss[i]==',' || wss[i]==';' )
		 	 {	if ( l>0 )	wss[i]=',';
				else		wss[i]=';';
			 }
		 }
		if ( l != 0 )
		 {	free(mtokens);free(wms);
			return(1);
		 }
		stokens=tokenize_char_dyn(wss,';');
		if ( stokens==NULL )
		 {	free(mtokens);free(wms);
			return(1);
		 }
		else if ( stokens[0]==NULL )
		 {	free(stokens);free(mtokens);free(wms);
			return(1);
		 }
		for ( j=0 ; stokens[j] != NULL ; j++ )
		 {	k=tokenize_char(stokens[j],cmd,'=',2);
			if ( k != 2 )
			 {	free(stokens);free(mtokens);free(wms);
				return(1);
			 }
			l=strlen(cmd[1]);
			if ( strcmp(cmd[0],"name")==0 )
			 {	strncpy(s->name,cmd[1],16);
				s->name[15]=0;
			 }
			else if ( strcmp(cmd[0],"offset")==0 && sscanf(cmd[1],"[%d,%d]",&x,&y)==2 )
			 {	s->x0=x;
				s->y0=y;
			 }
			else if ( strcmp(cmd[0],"image")==0 && cmd[1][0]=='[' && cmd[1][l-1]==']' )
			 {	cmd[1][l-1]=0;
				k=ficalib_parse_imagetrim_data(cmd[1]+1,&s->g.image);
				if ( k )
				 {	free(stokens);free(mtokens);free(wms);
					return(1);
				 }
			 }
			else if ( strcmp(cmd[0],"overscan")==0 && cmd[1][0]=='[' && cmd[1][l-1]==']' )
			 {	cmd[1][l-1]=0;
				
				k=ficalib_parse_overscan_data(cmd[1]+1,&o);
				if ( k )
				 {	free(stokens);free(mtokens);free(wms);
					return(1);
				 }
				s->g.overscans=(overscan *)realloc(s->g.overscans,sizeof(overscan)*(s->g.noverscan+1));
				s->g.overscans[s->g.noverscan]=o;
				s->g.noverscan++;
			 }
			else
			 {	free(stokens);free(mtokens);free(wms);
				return(1);
			 }
	 	 }

		free(stokens);
	 }
	else
	 {	int	x,y;
		k=tokenize_char(mtokens[n],cmd,'=',2);
		if ( k != 2 )
		 {	free(mtokens);free(wms);
			return(1);
		 }
		else if ( strcmp(cmd[0],"size")==0 && sscanf(cmd[1],"[%d,%d]",&x,&y)==2 )
		 {	m->sx=x;
			m->sy=y;
		 }
		else
		 {	free(mtokens);free(wms);
			return(1);
		 }
	 }
		
  }


 free(mtokens);
 free(wms);

 return(0);
}

/*****************************************************************************/

int main(int argc,char *argv[])
{
 char		**inputlist,**outputlist,**inmasklist,*rewriterule;
 fits		*bias,*dark,*flat;
 double		flatmean,darktime;
 FILE		*fr;
 int		ninimg;
 int		i,n,is_help,apply_mask,no_rewrite,no_history;
 char		*fdpstring,*overmodestr,*imgarea,*gainstr,
		*hdrdefstr,**overscanlist,*mosaicstr;
 compar		cpover;
 fitsdataparam	fdp;
 calibdata	cd;
 headername	hn;
 geometry	g;
 mosaic		m;

 progbasename=strrchr(argv[0],'/');
 if ( progbasename != NULL )	progbasename++;
 else				progbasename=argv[0];

 is_help=is_verbose=is_comment=no_rewrite=no_history=0;
 inputlist=outputlist=NULL;
 rewriterule=NULL;
 cd.inputbias=cd.inputdark=cd.inputflat=NULL;
 inmasklist=NULL;

 apply_mask=0;
 fdp.bitpix=-32;fdp.is_scale=0;fdp.bscale=1.0;fdp.bzero=0.0;
 fdpstring=NULL,hdrdefstr=NULL,overscanlist=NULL;

 cpover.mode=COM_MODE_MED;
 cpover.ignore_flag=0;
 cpover.logicalmethod=0;
 cpover.niter=0;
 cpover.lower=cpover.upper=3.0;
 overmodestr=NULL;

 hn.exptime="EXPTIME";
 hn.date="DATE-OBS";
 hn.time="TIME-OBS";
 hn.jd="JD";
 hn.hjd="HJD";
 hn.jy="JY";
 hn.hjy="HJY";
 hn.gain="GAIN";
 hn.gainpoly="GAINPOLY";
 hn.gainvmin="GAINVMIN";
 hn.extname="EXTNAME";

 imgarea=NULL;
 g.image.x0=0;
 g.image.y0=0;
 g.image.sx=0;
 g.image.sy=0;
 g.overscans=NULL;
 g.noverscan=0;

 cd.argc=argc;
 cd.argv=argv;

 cd.is_trim=0;
 cd.is_exptime_correction=1;
 cd.saturation=0;
 cd.reset_method=0;
 cd.postscale=0.0;
 cd.postmultiply=0.0;
 cd.is_remove_horizontal_stripes=0;

 cd.gain=1.0;
 cd.gainorder=0;
 gainstr=NULL;

 m.sections=NULL;
 m.nsection=0;
 m.sx=m.sy=0;

 mosaicstr=NULL;

 i=scanarg(argc,argv,SCANARG_ALLOW_FLAGS,
	"--version:%NS-1f%q",&is_help,
        "--version-short|--short-version:%NS-2f%q",&is_help,
	"-h|--help|--short-help|--help-short:%f%q",&is_help,
	"--long-help|--help-long:%SN2f%q",&is_help,
	"-b|--bitpix:%d",&fdp.bitpix,
	"--data:%s",&fdpstring,
	"-B|--input-master-bias:%s",&cd.inputbias,
	"-D|--input-master-dark:%s",&cd.inputdark,
	"-F|--input-master-flat:%s",&cd.inputflat,
	"--overscan-combination:%s",&overmodestr,
	"--overscan:%Dt",&overscanlist,
	"--image:%s",&imgarea,
	"--trim:%SN1f",&cd.is_trim,
	"--no-trim:%SN0f",&cd.is_trim,
	"--mosaic:%Cr",&mosaicstr,
	"--exptime-correction:%SN1f",&cd.is_exptime_correction,
	"--no-exptime-correction:%SN0f",&cd.is_exptime_correction,
	"-c|--scale|--post-scale:%g",&cd.postscale,
	"-m|--multiply|--post-multiply:%g",&cd.postmultiply,
        "--leak-lower-upper|--lu:%N1f",&cd.reset_method,
        "--leak-left-right|--lr:%N2f",&cd.reset_method,
        "--leak-any|--an:%N3f",&cd.reset_method,
	"--horizontal-stripe-removal:%f",&cd.is_remove_horizontal_stripes,
        "-s|--saturation:%g",&cd.saturation,
	"-i|--input:%Mt",&inputlist,
	"-o|--output:%Mt",&outputlist,
	"-r|--rewrite-output-name:%s",&rewriterule,
	"-M|--input-mask:%t",&inmasklist,
	"-a|--apply-mask:%f",&apply_mask,
	"--no-rewrite|--no-clobber:%f",&no_rewrite,
	"--clobber:%SN0f",&no_rewrite,
	"--no-history:%f",&no_history,
	"--history:%SN0f",&no_history,
	"-g|--gain:%s",&gainstr,
	"--gain-order:%d",&cd.gainorder,
	"--header:%s",&hdrdefstr,
	"--comment:%f",&is_comment,"(C):%f",&is_comment,
	"--verbose:%i",&is_verbose,"(V):%i",&is_verbose,
	"-*|+*:%e",
	"*:%l",&inputlist,
	NULL);

 if ( i ) 
  {	fprint_error("invalid command line argument near '%s'",argv[i]);
	return(1);
  }
 else if ( is_help<0 )
  {	fprint_generic_version(stdout,argv[0],"ficalib",FI_CALIB_VERSION,is_help);
	return(0);
  }
 else if ( is_help>1 )
  {	fprint_ficalib_long_help(stdout);
	return(0);
  }
 else if ( is_help )
  {	fprint_ficalib_usage(stdout);
	return(0);
  }

 if ( gainstr != NULL )
  {	/* automatic gain: */
	if ( strcmp(gainstr,"auto")==0 )
		cd.gain=0.0;
	/* note: by default (w/o any -g|--gain) the gain is 1 even if there  */
	/* is a GAIN keyword-value pair in the header...!		     */
	else if ( sscanf(gainstr,"%lg",&cd.gain)<1 )
		cd.gain=1.0;				
  }

 if ( imgarea != NULL && ficalib_parse_imagetrim_data(imgarea,&g.image) )
  {	fprint_error("unable to parse image trim data '%s'",imgarea);
	return(1);
  }

 if ( cd.postmultiply != 0.0 && cd.postscale != 0.0 )
  {	fprint_error("post-multiply and post-scale operations are conflicting\n");
	return(1);
  }

 if ( parse_fits_data_param(fdpstring,&fdp) )
  {	fprint_error("invalid pixel data format");
	return(1);
  }

 if ( combine_parse_mode(overmodestr,&cpover) )
  {	fprint_error("invalid overscan mode parameters in '%s'",overmodestr);
	return(1);
  }

 if ( hdrdefstr != NULL )
  {	i=scanpar(hdrdefstr,SCANPAR_DEFAULT,
		"exptime:%m",	&hn.exptime,
		"date:%m",	&hn.date,
		"time:%m",	&hn.time,
		"jd:%m",	&hn.jd,
		"hjd:%m",	&hn.hjd,
		"jy:%m",	&hn.jy,
		"hjy:%m",	&hn.hjy,
		"gain:%m",	&hn.gain,
		"gainpoly:%m",	&hn.gainpoly,
		"gainvmin:%m",	&hn.gainvmin,
		NULL);
	if ( i )
	 {	fprint_error("invalid header keyword assignment definition in '%s'",hdrdefstr);
		return(1);
	 }
  }

 if ( overscanlist != NULL )
  {	char		**os;
	overscan	o;

	for ( os=overscanlist ; *os != NULL ; os++ )
	 {	i=ficalib_parse_overscan_data(*os,&o);
		if ( i )
		 {	fprint_error("invalid overscan mode specification '%s'",*os);
			return(1);
		 }
		g.overscans=(overscan *)malloc(sizeof(overscan)*(g.noverscan+1));
		g.overscans[g.noverscan]=o;
		g.noverscan++;
	 }
  }

 if ( mosaicstr != NULL && ficalib_parse_mosaic_data(mosaicstr,&m) )
  {	fprint_error("unable to parse mosaic data '%s'",mosaicstr);
	return(1);
  }

 if ( ( g.image.sx <= 0 || g.image.sy <= 0 ) && g.noverscan>0 )
  {	fprint_warning("no trim section has been specified, overscans are ignored");
	if ( g.overscans != NULL )	free(g.overscans);
	g.noverscan=0;
	g.overscans=NULL;
  }
 else if ( g.image.sx>0 && g.image.sy>0 && g.noverscan>0 )
  {	trim	*o;
	for ( i=0 ; i<g.noverscan ; i++ )
	 {	o=&g.overscans[i].area;
		if ( ! overscan_check(&g.image,o) )
		 {	fprint_warning("overscan geometry %d:%d:%d:%d cannot be adjusted to the trim area, ignoring overscan correction request\n",
				o->x0,o->y0,o->x0+o->sx-1,o->y0+o->sy-1);
		 }
	 }
  }
	
 ninimg=count_list(inputlist);
 if ( outputlist==NULL )
	outputlist=inputlist;

 i=count_list(outputlist);
 if ( i != ninimg )
  {	fprint_error("number of input and output images differ");
	return(1);
  }
 if ( rewriterule != NULL )
  {	char	**newlist,*from,*to;
	if ( strchr(rewriterule,'|') == NULL )
	 {	from="*";
		to=rewriterule;
	 }
	else if ( strchr(rewriterule,'|') != strrchr(rewriterule,'|') )
	 {	fprint_error("invalid re-write pattern '%s'.",rewriterule);
		return(1);
	 }
	else
	 {	char	*p;
		int	flen;
		p=strchr(rewriterule,'|');
		flen=p-rewriterule;
		if ( flen<=0 )
		 {	from="*";
			to=p+1;
		 }
		else
		 {	from=(char *)malloc(flen+1);
			memcpy(from,rewriterule,flen);
			from[flen]=0;
			to=p+1;
		 }
	 }
	newlist=(char **)malloc(sizeof(char *)*(ninimg+1));
	for ( i=0 ; i<ninimg ; i++ )
	 {	newlist[i]=rewrite_filename(outputlist[i],from,to);
		if ( newlist[i]==NULL )
		 {	fprint_error("unexpected rewrite-rule '%s' for the filename '%s'.\n",rewriterule,outputlist[i]);
			return(1);
		 }
		/* fprintf(stderr,"%s -> %s\n",outputlist[i],newlist[i]); */
	 }
	outputlist=newlist;
	/* return(0); */
  }

 cd.sx=cd.sy=0;
 cd.mask=NULL; 

 cd.inmasklist=inmasklist;

 /* reading bias image if any */
 if ( cd.inputbias != NULL )	
  {	if ( (fr=fopenread(cd.inputbias))==NULL )
	 {	fprint_error("unable to open input bias image file '%s'",cd.inputbias);
		return(1);
	 }
	bias=fits_read(fr);
	fcloseread(fr);
	if ( bias==NULL )
	 {	fprint_error("unable to parse input bias file as FITS data.");
		return(1);
	 }
	if ( bias->i.dim != 2 )	
	 {	fprint_error("bias image dimension differs from 2.");
		return(1);
	 }
	fits_rescale(bias);
	cd.sx=bias->i.sx,
	cd.sy=bias->i.sy;
	cd.mask=fits_mask_read_from_header(&bias->header,cd.sx,cd.sy,NULL);
  }
 else	bias=NULL;

 /* reading dark image if any */
 if ( cd.inputdark != NULL )		/* read from file (-iD ...)  */
  {	if ( (fr=fopenread(cd.inputdark))==NULL )
	 {	fprint_error("unable to open input dark image file '%s'",cd.inputdark);
		return(1);
	 }
	dark=fits_read(fr);
	if ( dark==NULL )
	 {	fprint_error("unable to parse input dark file as FITS data.");
		return(1);
	 }
	if ( dark->i.dim != 2 )	
	 {	fprint_error("dark image dimension differs from 2.");
		return(1);
	 }
	if ( cd.sx>0 && cd.sy>0 && ( dark->i.sx != cd.sx || dark->i.sy != cd.sy ) )
	 {	fprint_error("dark image size differs from the expected image size.");
		return(1);
	 }
	fits_rescale(dark);
	if ( cd.mask==NULL )
 	 {	cd.sx=dark->i.sx,
		cd.sy=dark->i.sy;
		cd.mask=fits_mask_read_from_header(&dark->header,cd.sx,cd.sy,NULL);
	 }	
	else
		fits_mask_mask_from_header(cd.mask,&dark->header,cd.sx,cd.sy,NULL);
  }
 else	dark=NULL;

 if ( (! cd.is_exptime_correction) && dark != NULL && bias != NULL )
  {	int	i,j;
	for ( i=0 ; i<cd.sy ; i++ )
	 {	for ( j=0 ; j<cd.sx ; j++ )
		 {	dark->i.data[i][j] += bias->i.data[i][j];	}
	 }
	fits_free(bias);	/* we do it because of this step, mainly */
	bias=NULL;
	cd.is_combined_dark=1;
  }
 else
	cd.is_combined_dark=0;

 /* reading flat image if any */
 if ( cd.inputflat != NULL )		/* read from file (-iF ...)  */
  {	if ( (fr=fopenread(cd.inputflat))==NULL )
	 {	fprint_error("unable to open input flat image file '%s'",cd.inputflat);
		return(1);
	 }
	flat=fits_read(fr);
	if ( flat==NULL )
	 {	fprint_error("unable to parse input flat file as FITS data.");
		return(1);
	 }
	if ( flat->i.dim != 2 )	
	 {	fprint_error("flat image dimension differs from 2.");
		return(1);
	 }
	if ( cd.sx>0 && cd.sy>0 && ( flat->i.sx != cd.sx || flat->i.sy != cd.sy ) )
	 {	fprint_error("flat image size differs from the expected image size.");
		return(1);
	 }
	fits_rescale(flat);
	if ( cd.mask==NULL )
 	 {	cd.sx=flat->i.sx,
		cd.sy=flat->i.sy;
		cd.mask=fits_mask_read_from_header(&flat->header,cd.sx,cd.sy,NULL);
	 }	
	else
		fits_mask_mask_from_header(cd.mask,&flat->header,cd.sx,cd.sy,NULL);
  }
 else	flat=NULL;

 if ( flat != NULL )
  {	flatmean=calibrate_get_flat_mean(&flat->i,cd.mask);
	cd.flatpoly=(double *)malloc(sizeof(double)*(cd.gainorder+1)*(cd.gainorder+2)/2);
	calibrate_fit_flatpoly(flat->i.data,cd.sx,cd.sy,cd.mask,cd.gainorder,cd.flatpoly,&cd.flatresd,&cd.flatvmin);
  }
 else
  {	flatmean=0.0;
	cd.flatpoly=NULL;
  }

 if ( dark != NULL )
  {	if ( fits_headerset_get_as_double(&dark->header,hn.exptime,&darktime,0) )
		darktime=-1.0;
  }
 else
	darktime=-1.0;	

 srand48_auto();

 /* do the calibration for all input images: */
 for ( n=0 ; n<ninimg ; n++ )
  {	fitsimage	*ibias,*idark,*iflat;
	struct stat	st;

	/* just to have fitsimage's instead of fits's (will be better for    */
	/* further developments when ficalib can handle mutiple extensions)  */
	if ( bias != NULL )	ibias=&bias->i;
	else			ibias=NULL;
	if ( dark != NULL )	idark=&dark->i;
	else			idark=NULL;
	if ( flat != NULL )	iflat=&flat->i;
	else			iflat=NULL;

	/* check whether the output file exists (and it is a regular file): */
	if ( no_rewrite && (!stat(outputlist[n],&st)) && S_ISREG(st.st_mode) )
	 {	if ( is_verbose )
		 {	fprint_warning("output file '%s' exists, re-calibration is skipped",outputlist[n]);	}
		continue;
	 }

	/* the main calibration stuff per file: */
	calibrate_image(inputlist[n],outputlist[n],&g,&m,&cd,&cpover,
	&fdp,&hn,ibias,idark,iflat,flatmean,darktime,no_history);

  } /* for ( n=0 ; n<ninimg ; n++ ) */

 return(0);
}

/*****************************************************************************/
                        
