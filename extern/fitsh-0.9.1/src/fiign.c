/*****************************************************************************/
/* fiign.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line user interface for ignorig false pixels from a (FITS) image. */
/*****************************************************************************/
#define	FI_IGN_VERSION		"0.9e"
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
#include "math/poly.h"
#include "statistics.h"
#include "math/spline/biquad.h"
#include "math/spline/biquad-isc.h"
#include "io/iof.h"
#include "io/scanarg.h"
#include "io/tokenize.h"
#include "maskdraw.h"

#include "tensor.h"
#include "common.h"
#include "history.h"

/*****************************************************************************/

int	is_verbose,is_comment;

/*****************************************************************************/

typedef struct
 {	int	match;
	int	value;
	int	reset;
	int	set;
 } maskconvert;

typedef struct
 {	int	x0,y0;
	int	sx,sy;
	int	mask;
 } maskblock;

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

int ignore_cosmics(fitsimage *img,char **mask,double th_low,double th_high,int is_repl,double skysigma)
{
 int	i,j,k,l,hsize,fsize,ii,jj;
 int	sx,sy;
 double	arr[25],med,s,s2,sig,w;
 int	is_low,is_high;

 if ( img==NULL )	return(1);
 if ( img->data==NULL )	return(0);
 sx=img->sx,sy=img->sy;

 is_low=is_high=0;
 if ( th_low  > 0.0 )	is_low =1;
 if ( th_high > 0.0 )	is_high=1;
 
 hsize=1;
 fsize=2*hsize+1;

 for ( i=hsize ; i<sy-hsize ; i++ )
  { for ( j=hsize ; j<sx-hsize ; j++ )
     {	l=0;
	for ( k=0,ii=jj=-hsize ; k<fsize*fsize ; k++ )
	 {	if ( ! (mask[i+ii][j+jj]&MASK_FAULT) && ! ( ii == 0 && jj == 0 ) )
		 {	w=img->data[i+ii][j+jj];
			arr[l]=w;l++;
		 }
		if ( jj==hsize )	ii++,jj=-hsize;
		else			jj++;
	 }
	if ( l<fsize*fsize-1 )	continue;

	med=median(arr,l);
	s=s2=0.0;
	for ( k=0 ; k<l ; k++ )	s+=arr[k],s2+=arr[k]*arr[k];
	s/=l,s2/=l;sig=sqrt(s2-s*s);

	if ( skysigma>0.0 && sig>2*skysigma )	continue;

	if ( skysigma>0.0 )	sig=skysigma;

	if ( is_low  && img->data[i][j] < s - th_low *sig )
	 {	mask[i][j] |= MASK_COSMIC;
		if ( is_repl )
		 {	img->data[i][j]=s,
			mask[i][j] |= MASK_INTERPOLATED;
		 }
	 }
	if ( is_high && img->data[i][j] > s + th_high*sig )
	 {	mask[i][j] |= MASK_COSMIC;
		if ( is_repl )
		 {	img->data[i][j]=s,
			mask[i][j] |= MASK_INTERPOLATED;
		 }
	 }
     }
  }
 return(0);
}
/*
int ignore_cosmics_biquad(fits *img,char **mask)
{
 int	sx,sy,i,j,k,l,subg,n,mi,mj;
 double	**bqc,**subc,mmin,nbmax,nbmin,pval,w,s,s2;

 if ( mask==NULL )	return(1);
 if ( img==NULL )	return(1);
 if ( img->data==NULL )	return(1);
 sx=img->sx,sy=img->sy;
 if ( sx<=0 || sy<=0 )	return(1);
 
 subg=4;
 bqc=tensor_alloc_2d(double,2*sx+1,2*sy+1);
 if ( bqc==NULL )	return(-1);
 subc=tensor_alloc_2d(double,subg,subg);
 if ( subc==NULL )	return(-1);

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( mask[i][j] )	mask[i][j]=1;	}
  }
 biquad_coeff(img->data,sx,sy,bqc,mask);
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( ! mask[i][j] )	continue;
		biquad_isc_int_block_subpixels(bqc,j,i,subg,subg,subc);
		mmin=subc[0][0];
		for ( k=0 ; k<subg ; k++ )
		 {	for ( l=0 ; l<subg ; l++ )
			 {	if ( subc[k][l]<mmin )	mmin=subc[k][l];	}
		 }
		mmin=mmin*(double)(subg*subg);
		pval=img->data[i][j];
		n=0;s=s2=nbmax=nbmin=0.0;mi=mj=0;
		if ( j>0 && mask[i][j-1] )
		 {	w=img->data[i][j-1];
			if ( n==0 || w>nbmax )	nbmax=w,mi=i,mj=j-1;
			if ( n==0 || w<nbmin )	nbmin=w;
			s+=w,s2+=w*w,n++;
		 }
		if ( j<sx-1 && mask[i][j+1] )
		 {	w=img->data[i][j+1];
			if ( n==0 || w>nbmax )	nbmax=w,mi=i,mj=j+1;
			if ( n==0 || w<nbmin )	nbmin=w;
			s+=w,s2+=w*w,n++;
		 }
		if ( i>0 && mask[i-1][j] )
		 {	w=img->data[i-1][j];
			if ( n==0 || w>nbmax )	nbmax=w,mi=i-1,mj=j;
			if ( n==0 || w<nbmin )	nbmin=w;
			s+=w,s2+=w*w,n++;
		 }
		if ( i<sy-1 && mask[i+1][j] )
		 {	w=img->data[i+1][j];
			if ( n==0 || w>nbmax )	nbmax=w,mi=i+1,mj=j;
			if ( n==0 || w<nbmin )	nbmin=w;
			s+=w,s2+=w*w,n++;
		 }
		fprintf(stdout,"%12g %12g %12g %12g\n",nbmin,nbmax,mmin,pval);
		if ( n>=2 && mmin<nbmin && mmin<pval && nbmax>pval )
		 {	s-=nbmax,s2-=nbmax*nbmax;
			s/=(double)(n-1),s2/=(double)(n-1);
			s2=sqrt(s2-s*s);
			if ( 3.0*(nbmin-mmin)<(nbmax-pval) && 100.0*s2<nbmax-nbmin )
			 {	mask[mi][mj]=0;
				fprintf(stdout,"Ign: %d,%d\n",mj+1,mi+1);
			 }
		 }
	 }
  }

 tensor_free(subc);
 tensor_free(bqc);
 return(0);
}
*/

/*****************************************************************************/

int fprint_fiign_usage(FILE *fw)
{
 fprintf(fw,
"Usage:\tfiign [[-i] <in>[<F>] [--frame <F>]] [-o <out>]\n"
"\t[-h|--help] [-V|--verbose] [-C|--comment]\n");
 fprintf(fw,
"General options:\n"
"\t[--ignore-mask | -M|--input-mask <mask>]\n"
"\t[-a|--apply-mask [-m|--mask-value <value>]] [--output-mask <mask>]\n");
 fprintf(fw,
"Mask conversion:\n"
"\t[--convert <match>:<value>:<reset>:<set> [--convert <>:<>:<>:<>] ...]\n");
 fprintf(fw,
"Marking saturated pixels:\n"
"\t[-s <saturation-level>|-S <saturation_image.fits>]\n"
"\t[--leak-{left-right|lower-upper|any}|--lr|--lu|--an]\n");
 fprintf(fw,
"General ignorance:\n"
"\t[-n|--ignore-nonpositive] [-g|--ignore-negative] [-z|--ignore-zero]\n");
 fprintf(fw,
"Cosmics removal and/or interpolation:\n"
"\t[-c|--ignore-cosmics [-r|--replace-cosmics]]\n");

 return(0);
}

longhelp_entry fiign_long_help[]=
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
	"Name of the output FITS image file (can be the same as the input "
	"image file)." },

 { "Generic pixel masking options:", NULL },
 { "-n, --ignore-nonpositive",
	"Mask pixels with non-positive values." },
 { "-g, --ignore-negative",
	"Mask pixels with negative values." },
 { "-z, --ignore-zero",
	"Mask pixels with a value of zero." },
 { "--ignore-mask",
	"Completely ignore mask associated to the input image." },
 { "-M, --input-mask <image file>",
	"Input mask file to co-add to output image mask." },
 { "-a, --apply-mask",
	"Apply the mask to the image, i.e. set the pixel values with non-null "
	"mask to be zero (by default, or any other value specified by "
	"``-m|--mask-value'')." },
 { "-m, --mask-value <value>",
	"Override the default pixel value (zero) during explicit marking of "
	"masked pixels (see also ``-a|--apply-mask'')." },

 { "Marking saturated pixels:", NULL },
 { "-s <saturation-level>",
	"Saturation level." },
 { "-S <image file>",
	"Image containing saturation level on a per pixel basis." },
 { "--leak-left-right, --leak-lower-upper, --leak-any",
	"Readout direction, i.e. orientation of ``blooming'' stripes." },

 { "Mask conversion:", NULL },
 { "--convert <match>:<value>:<reset>:<set>",
	"Convert masks: from a mask which matches to the <match>:<value> pair, "
	"i.e. the masks with the type of <match> have a value of <value> "
	"the masks specified by <reset> are cleared and the masks specified by "
	"<set> are set. The <match>, <value>, <reset> and <set> tags are "
	"comma-separated list of mask names (see below)." },

 { "Mask names:", NULL },
 { "none",
	"no mask at all" },
 { "clear",
	"same as ``none''" },
 { "fault",
	"mask for faulty pixels" },
 { "hot",
	"mask for hot pixels" },
 { "cosmic",
	"mask for marking cosmic pixels" },
 { "outer",
	"pixels originating from out of image areas" },
 { "oversaturated",
	"oversaturated pixels" },
 { "bloomed",
	"``bloomed'' pixels (i.e. not oversaturated but neighbouring pixel(s) may be so)" },
 { "saturated",
	"oversaturated or bloomed pixels" },
 { "interpolated",
	"pixels havin an interpolated value (e.g. hot or cosmic pixels are "
	"replaced by the average value of the surrounding pixels)." },
 { NULL, NULL }
};


int fprint_fiign_long_help(FILE *fw)
{
 fprintf(fw,
"Usage:\tfiign [options and operations] [<input>] [-o|--output <output>]\n"
"In the context of FITS image data processing, ``masks'' are per pixel\n"
"associated meta-data, representing the state of the given pixel. In\n"
"general, pixels considered to be somehow ``bad'' are marked with these masks\n"
"in order to exclude or use only with caution during further processing.\n");
 fprintf(fw,
"These masks can either mark the initial state of the given pixels (e.g.\n"
"pixels can be marked as hot or bad pixels, which describes the detector itself\n"
"and not the individual scientific or calibration frames), or masks can be\n"
"added during the subsequent steps of the processing (e.g. saturated pixels,\n"
"``outer'' pixels).\n\n");
 fprintf(fw,
"The purpose of the `fiign` program is to give a low-level\n"
"access to these masks. Altough the operations on the images automatically\n"
"yields the respective operations on the masks (e.g. if an image is transformed\n"
"or trimmed, the associated mask will also be transformed or trimmed\n"
"with the same geometry), with this program the masks can be manipulated\n"
"arbitrarily.\n\n");

 longhelp_fprint(fw,fiign_long_help,0,-1);

 fprintf(fw,"\n");
 fprintf(fw,"Report bugs to <%s>\n",FI_MAINT_EMAIL);
 
 return(0);
}

/*****************************************************************************/


int main(int argc,char *argv[])
{
 fits		*img,*satimg;
 char		**mask;

 FILE		*fw,*ft,*fr;
 int		i,j,sx,sy,is_help,is_no_mask;
 double		saturation,skysigma,threshold,thlow,thhigh,maskvalue;
 char 		*satimgname,*inimgname,*outimgname,*outmaskname,**inmasklist,
		*ignmasknames,**convertlist,*basename,**maskblocklist;
 int		ign_cosmics,rep_cosmics,ign_num,frameno,
		reset_method,apply_mask;
 maskconvert	*mcls;
 int		nmcl;

 progbasename=strrchr(argv[0],'/');
 if ( progbasename != NULL )	progbasename++;
 else				progbasename=argv[0];

 saturation=0.0;threshold=4.0;skysigma=0.0;
 thlow=10.0;thhigh=50.0;
 is_no_mask=is_verbose=is_comment=is_help=0;
 reset_method=0;apply_mask=0;ign_num=0;ign_cosmics=rep_cosmics=0;
 inimgname=outimgname=outmaskname=satimgname=NULL;frameno=0;
 ignmasknames=NULL;inmasklist=convertlist=NULL;maskvalue=0.0;
 maskblocklist=NULL;
 mcls=NULL;nmcl=0;

 i=scanarg(argc,argv,SCANARG_ALLOW_FLAGS,
	"--version:%NS-1f%q",&is_help,
        "--version-short|--short-version:%NS-2f%q",&is_help,
	"-h|--help:%f%q",&is_help,
	"--long-help|--help-long:%SN2f%q",&is_help,
	"-i|--input:%s",&inimgname,
	"-q|--mask-block:%Dt",&maskblocklist,
	"--frame:%d",&frameno,
	"-o|--output:%s",&outimgname,
	"--output-mask:%s",&outmaskname,
	"-M|--input-mask:%t",&inmasklist,
	"-d:%g",&skysigma,
	"-t|--threshold:%g",&threshold,
	"--threshold-low:%g",&thlow,
	"--threshold-high:%g",&thhigh,
	"-c|--ignore-cosmics:%f",&ign_cosmics,
	"-r|--replace-cosmics:%f%f",&ign_cosmics,&rep_cosmics,
	"-z|--ignore-zero:%N1f",&ign_num,
	"-g|--ignore-negative:%N2f",&ign_num,
	"-n|--ignore-nonpositive:%N3f",&ign_num,
	"-a|--apply-mask:%f",&apply_mask,
	"-m|--mask-value:%g",&maskvalue,
	"--no-mask|--ignore-mask:%f",&is_no_mask,
	"--mask-ignore:%s",&ignmasknames,
	"--convert:%t",&convertlist,
	"--lu:%N1f",&reset_method,
	"--lr:%N2f",&reset_method,
	"--an:%N3f",&reset_method,
	"-s|--saturation:%g",&saturation,
	"-S|--saturation-image:%s",&satimgname,
	"--comment:%i",&is_comment,"(C):%i",&is_comment,
	"--verbose:%i",&is_verbose,"(V):%i",&is_verbose,
	"-:%w",&inimgname,
	"-*|+*:%e",
	"*:%w",&inimgname,
	NULL);

 if ( i )		
  {	fprint_error("invalid command line argument near '%s'",argv[i]);
	return(1);
  }
 if ( is_help<0 )
  {	fprint_generic_version(stdout,argv[0],"fiign",FI_IGN_VERSION,is_help);
	return(0);
  }
 else if ( is_help>1 )
  {	fprint_fiign_long_help(stdout);
	return(0);
  }
 else if ( is_help )
  {	fprint_fiign_usage(stdout);
	return(0);
  }

 if ( convertlist != NULL )
  {	char		*cc,*cmd[5];
	int		i,k,n,match,value,reset,set,is_any;
	maskconvert	*wm;

	for ( i=0 ; convertlist[i] != NULL ; i++ )
	 {	cc=strdup(convertlist[i]);
		/*fprintf(stderr,"[%s]\n",convertlist[i]);*/
		n=tokenize_char(cc,cmd,':',4);
		match=value=reset=set=-1,is_any=0;
		if ( n==4 )
		 {	match=parse_mask_flags(cmd[0]);
			if ( strcmp(cmd[1],"any")==0 )	is_any=1,value=0;
			else	is_any=0,value=parse_mask_flags(cmd[1]);
			reset=parse_mask_flags(cmd[2]);
			set  =parse_mask_flags(cmd[3]);
		 }
		free(cc);
		if ( match<=0 || value<0 || reset<0 || set<0 || ( ! is_any && value==0 ) )
		 {	fprintf(stderr,"Warning: invalid conversion specification '%s', skipped.\n",convertlist[i]);
			continue;	
		 }
		if ( is_any )
		 {	for ( k=1 ; k & MASK_ALL ; k=k<<1 )
			 {	if ( ! ( k&match ) )	continue;
				mcls=(maskconvert *)realloc(mcls,sizeof(maskconvert)*(nmcl+1));
				wm=&mcls[nmcl];
				wm->match=k;
				wm->value=k;
				wm->reset=reset;
				wm->set=set;
				nmcl++;
			 }
		 }			
		else
		 {	mcls=(maskconvert *)realloc(mcls,sizeof(maskconvert)*(nmcl+1));
			wm=&mcls[nmcl];
			wm->match=match;
			wm->value=value;
			wm->reset=reset;
			wm->set=set;
			nmcl++;
		 }
	 }
  }

 basename=fits_basename(inimgname,&frameno);
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
 if ( img->i.dim != 2 )
  {	fprint_error("image dimension differs from 2.");
	return(1); 
  }
 fits_rescale(img);

 sx=img->i.sx,
 sy=img->i.sy;
 if ( is_no_mask )
	mask=fits_mask_create_empty(sx,sy);
 else
	mask=fits_mask_read_from_header(&img->header,sx,sy,NULL);

 if ( inmasklist != NULL )
  {	if ( join_masks_from_files(mask,sx,sy,inmasklist) )
	 {	fprint_error("unable to read one of the input mask files");
		return(1);
	 }
  }
 fits_mask_mark_nans(&img->i,mask,MASK_NAN);

 if ( ignmasknames != NULL )
  {	int	flag;
	flag=parse_mask_flags(ignmasknames);
	if ( flag<0 )	
	 {	fprint_error("invalid mask specification '%s'",ignmasknames);
		return(1);
	 }
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	mask[i][j] &= ~flag;		}
	 }
  }



 if ( satimgname != NULL )
  {	ft=fopenread(satimgname);
	if ( ft==NULL )	
	 {	fprint_error("unable to open saturation level image '%s'",satimgname);
		return(1);
	 }
	satimg=fits_read(ft);
	if ( satimg==NULL )
	 {	fprint_error("unable to read staturation image file '%s'",satimgname);
		return(1);
 	 }
	if ( satimg->i.dim != 2 )
	 {	fprint_error("saturation image is not 2D");
		return(1);
	 }
	fcloseread(ft);
	if ( sx != satimg->i.sx || sy != satimg->i.sy )
	 {	fprint_error("input and saturation image size differ");
		return(1);
	 }
  }
 else	satimg=NULL;

 if ( satimg != NULL )
	mark_saturated_pixels(&img->i,mask,&satimg->i,1.0,reset_method);
 else if ( saturation > 0.0 )
	mark_saturated_pixels(&img->i,mask,NULL,saturation,reset_method);

 for ( i=0 ; maskblocklist != NULL && maskblocklist[i] != NULL ; i++ )
  {     if ( maskdraw_parse_and_draw(mask,sx,sy,maskblocklist[i])>0 )
         {      fprint_warning("invalid mask block specification '%s', skipped.\n",maskblocklist[i]);   }
  }

 if ( ign_num )
  {	int	cm;
	for ( i=0 ; i<sy ; i++ )
	 { for ( j=0 ; j<sx ; j++ )
	    {	if ( img->i.data[i][j]==0.0 )		cm=1;
		else if ( img->i.data[i][j]< 0.0 )	cm=2;
		else					cm=0;
		if ( cm & ign_num )	mask[i][j] |= MASK_FAULT;
	    }
	 }
  }

 if ( ign_cosmics )
	ignore_cosmics(&img->i,mask,thlow,thhigh,rep_cosmics,skysigma);

 for ( i=0 ; i<nmcl ; i++ )
  {	maskconvert	*wm;
	int		m;
	wm=&mcls[i];
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	m=mask[i][j];
			if ( (m & wm->match)==(wm->value & wm->match) )
				mask[i][j]=(m&(~wm->reset))|wm->set;
		 }
	 }
  }
			

 if ( apply_mask )
  {	for ( i=0 ; i<sy ; i++ )
	 { for ( j=0 ; j<sx ; j++ )
	    {	if ( mask[i][j] )	img->i.data[i][j]=maskvalue;	}
	 }
  }

 fits_history_export_command_line(img,"fiign",FI_IGN_VERSION,argc,argv);

 fits_set_image_params(img);
 fits_backscale(img,img->i.read.bscale,img->i.read.bzero);
 mark_integerlimited_pixels(&img->i,mask,img->i.bit,1,MASK_OVERSATURATED,MASK_OVERSATURATED);
 fits_mask_export_as_header(&img->header,1,mask,sx,sy,NULL);

 if ( outimgname==NULL )	fw=stdout;
 else				fw=fopenwrite(outimgname);
 if ( fw==NULL )		
  {	fprint_error("unable to create output image file '%s'",outimgname);
	return(1);
  }
 fits_write(fw,img);
 fclosewrite(fw);

 if ( outmaskname != NULL )
  {	fw=fopenwrite(outmaskname);
	if ( fw==NULL )		
	 {	fprint_error("unable to create output mask file '%s'",outmaskname);
		return(1);
	 }
	fits_write_header(fw,img);
	fclosewrite(fw);
  }

 fits_free(img);
 if ( satimg != NULL )	fits_free(satimg);
 fits_mask_free(mask);

 return(0);
}

