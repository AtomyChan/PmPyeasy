/*****************************************************************************/
/* ficombine.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line user interface for combining more FITS images.		     */
/*****************************************************************************/
#define	FITSH_FICOMBINE_VERSION	"0.9"
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include <fits/fits.h>

#include "longhelp.h"
#include "fitsh.h"

#include "fitsmask.h"
#include "statistics.h"
#include "io/iof.h"
#include "io/scanarg.h"

#include "tensor.h"
#include "common.h"
#include "combine.h"

#include "history.h"


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

int fprint_ficombine_usage(FILE *fw)
{
 fprintf(fw,
"Usage:\tficombine [-h|--help|--long-help|--wiki-help] [--version[-short]]\n"
"\t<input> ... [-o|--output <output.fits>]\n"
"\t[-C|--comment] [-V|--verbose] [--[no-]history] [-M|--input-mask <mask>]\n"
"\t[-b|--bitpix <bitpix>] [--max-memory <max-memory>]\n"
"\t[--data bitpix=<bitpix>,bscale=<scale>,bzero=<zero>|<C-type>]\n");
 fprintf(fw,
"\t[-m|--mode mean|median|rejection|sum|squaresum|scatter,or|and,\n"
"\t ignorenegative,iterations=<n>,[lower|upper|sigma]=<limit>,min,max,\n"
"\t truncated|winsorized,lowest=<cnt>,highest=<cnt>]\n"
"\t[-n|--ignore-negative] [--logical-or|--logical-and]\n");
 return(0);
}

longhelp_entry ficombine_long_help[]=
{
 LONGHELP_OPTIONS,
 { "General options:", NULL },
 { "-h, --help",
	"Gives general summary about the command line options." },
 { "--long-help, --help-long",
	"Gives a detailed list of command line options." },
 { "--wiki-help, --help-wiki, --mediawiki-help, --help-mediawiki",
	"Gives a detailed list of command line options in Mediawiki format." },
 { "--version, --version-short, --short-version",
	"Gives some version information about the program." },
 { "-o, --output <fits>",
	"The name of the output file (omitting or specifing '-' yields "
	"the output to be written to stdout)." },
 { "--history, --no-history",
	"Incude or exclude the command line of the invoking `ficombine` call. "
	"By default, the command line is added into the FITS headers (FI_HSTRY) "
	"but if there are enormous amount of input files, one can use "
	"--no-history to omit the respective lengthy keywords." },
 { "--data <spec>",
	"Output pixel data format specification." },
 { "-b, --bitpix <bitpix>",
	"Standard FITS output bitpix value." },
 { "-M, --input-mask <fits>",
	"Input mask file to co-add to output image." },

 { "Generic combination specification:", NULL },
 { "-m, --mode <mode>",
	"Use the specified mode to combine images." },
 { "-n, --ignore-negative",
	"Ignore (i.e. mask) pixels with negative values." },
 { "--logical-or",
	"Use logical \"or\" combination between masks." },
 { "--logical-and",
	"Use logical \"and\" combination between masks." }, 

 { "Combination modes (a comma-separated list of these should follow -m|--mode):", NULL },
 { "mean",
	"The mean value of the pixels." },
 { "median",
	"The median value of the pixels." },
 { "min, minimum",
	"The minimum value of the pixels." },
 { "max, maximum",
	"The maximum value of the pixels." },
 { "rejmed, rejection",
	"Median value after rejecting outliers." },
 { "rejmean",
	"Mean value after rejecting outliers." },
 { "iterations=<n>",
	"Number of iterations to reject outliers." },
 { "lower=<l>, upper=<u>, sigma=<s>",
	"Outlier limit (lower, upper, both) in standard deviation. Note that " 
	"setting sigma to some value is equivalent to setting 'lower' and 'upper' "
	"to the same value simultaneously." },
 { "truncated",
	"Truncated mean." },
 { "winsorized",
	"Winsorized mean." },
 { "lowest=<cnt>, highest=<cnt>, discard=<cnt>",
	"Reject the 'lowest' and 'highest' count of points when computing "
	"truncated or winsorized mean. Setting 'discard' to some value "
	"is equivalent to setting 'lowest' and 'highest' to this value "
	"simultaneously. " },
 { "sum",
	"Sum of the pixel values." },
 { "squaresum",
	"Sum for the squarers of the pixel values." },
 { "scatter, stddev",
	"Pixel scatter (standard deviation)." },
 { "or",
	"Use logical \"or\" combination between masks." },
 { "and",
	"Use logical \"and\" combination between masks." },
 { "ignorenegative",
	"Ignore (i.e. mask) pixels with negative values." },
 { "ignorezero",
	"Ignore (i.e. mask) pixels with a zero value." },
 { "ignorenegative",
	"Ignore (i.e. mask) pixels with positive values." },

 { NULL,NULL }
};

int fprint_ficombine_long_help(FILE *fw,int is_wiki)
{
 char	*synopsis=
	"ficombine [options] <input images> [-o|--output <output>]";
 char	*description=
	"The purpose of this program is to combine the input images (with the same "
	"sizes) to a single image. This combination refers to a kind of \"averaging\" "
	"for the images, however, other modes are also available.";

 fprint_generic_long_help(fw,is_wiki,ficombine_long_help,synopsis,description);
 
 return(0);
}



/*****************************************************************************/

int main(int argc,char *argv[])
{
 char			**inputfiles,**inmasklist,*outname,*outmaskname,
			*fdpstring,*imagmodestr,*maxmemstr;

 char			**mask,**outmask;

 int			ninput,sx,sy;
 fits			*outimg;
 FILE			*fw;
 comimg			*inputs;
 int			i,j,is_help,is_add_history,apply_mask;
 combine_parameters	cp_data,*cp=&cp_data;
 fitsdataparam		fdp;
 size_t			maxmem;

 progbasename=strrchr(argv[0],'/');
 if ( progbasename != NULL )	progbasename++;
 else				progbasename=argv[0];

 is_help=is_verbose=is_comment=0;
 inmasklist=inputfiles=NULL;
 outmaskname=outname=NULL;
 imagmodestr=NULL;
 is_add_history=1;

 combine_parameters_reset(cp);
 cp->mode=COM_MODE_MEAN;
 cp->ignore_flag=0;
 cp->logicalmethod=0; 

 apply_mask=0;
 fdp.bitpix=0;fdp.is_scale=0;fdp.bscale=1.0;fdp.bzero=0.0;fdpstring=NULL;
 maxmemstr=NULL;

 i=scanarg(argc,argv,SCANARG_ALLOW_FLAGS,
	"--version:%NS-1f%q",&is_help,
        "--version-short|--short-version:%NS-2f%q",&is_help,
	"-h|--help|--short-help|--help-short:%f%q",&is_help,
	"--long-help|--help-long:%SN2f%q",&is_help,
	"--mediawiki-help|--help-mediawiki|--wiki-help|--help-wiki:%SN3f%q",&is_help,
	"--history:%SN1f",&is_add_history,
	"--no-history:%SN0f",&is_add_history,
	"-o|--output:%s",&outname,
	"-b|--bitpix:%d",&fdp.bitpix,
	"--data:%s",&fdpstring,
	"--max-memory:%s",&maxmemstr,
	"-m|--mode:%s",&imagmodestr,
	"-M|--input-mask:%Dt",&inmasklist,
	"--output-mask:%s",&outmaskname,
	"-a|--apply-mask:%f",&apply_mask,
	"-n|--ignore-negative:%0f",&cp->ignore_flag,
	"--logical-or:%SN0f",&cp->logicalmethod,
	"--logical-and:%SN1f",&cp->logicalmethod,
	"--comment:%f",&is_comment,"(C):%f",&is_comment,
	"--verbose:%i",&is_verbose,"(V):%i",&is_verbose,
	"-*|+*:%e",
	"*:%l",&inputfiles,
	NULL);

 if ( i )
  {	fprint_error("invalid command line argument near '%s'",argv[i]);
	return(1); 
  }
 else if ( is_help<0 )
  {	fprint_generic_version(stdout,argv[0],"ficombine",FITSH_FICOMBINE_VERSION,is_help);
	return(0);
  }
 else if ( 1<is_help )
  {	fprint_ficombine_long_help(stdout,2<is_help);
	return(0);
  }
 else if ( is_help )
  {	fprint_ficombine_usage(stdout);
	return(0);
  }

 if ( combine_parse_mode(imagmodestr,cp) )
  {	fprint_error("invalid combination mode '%s'",imagmodestr);
	return(1);
  }

 if ( parse_fits_data_param(fdpstring,&fdp) )
  {	fprint_error("invalid pixel data format");
	return(1);
  }

 if ( maxmemstr != NULL )
  {	maxmem=parse_max_memory_string(maxmemstr);
	if ( maxmem<=1 )
	 {	fprint_error("invalid maximum memory specification '%s'.\n",maxmemstr);
		return(1);
	 }
  }
 else
	maxmem=0;

 if ( inputfiles==NULL )	ninput=0;
 else	
  {	for ( ninput=0 ; inputfiles[ninput] != NULL ; ninput++ ) ;	}
 if ( ninput==0 )	return(0);	/* no input files: exit successfully */

 inputs=(comimg *)malloc(sizeof(comimg)*ninput);
 
 for ( i=0 ; i<ninput ; i++ )
  {	FILE	*fr;
	fits	*img;
	
	if ( (fr=fopenread(inputfiles[i]))==NULL )
	 {	fprint_error("unable to open input file '%s'.",inputfiles[i]);
		return(1);
	 }
	img=fits_seek_frame_to_image(fr,0);
	if ( img==NULL )
	 {	fprint_error("unable to parse input as FITS data.");
		return(1);
	 }
	else if ( img->i.dim<2 || img->i.dim>3 )
	 {	fprint_error("image dimension must be 2 or 3.");
		return(1);
	 }
	else if ( img->i.dim==3 )
	 {	fits_alloc_image_gen(img,img->i.dim,img->i.naxis);
		fits_read_image(fr,img);
		fclose(fr);
		fits_rescale(img);
		inputs[i].img=img;
		inputs[i].fr=NULL;
	 }
	else
	 {	inputs[i].img=img;
		inputs[i].fr=fr;
		img->i.vdata=NULL;	
	 }
  }

 sx=inputs[0].img->i.sx,
 sy=inputs[0].img->i.sy;
 for ( i=1,j=0 ; i<ninput && ! j ; i++ )
  {	if ( inputs[i].img->i.sx != sx || inputs[i].img->i.sy != sy )
		j=1;
  }
 if ( j )
  {	fprint_error("size of input images differ");
	return(1);
  }

 mask=fits_mask_create_empty(sx,sy);
 if ( inmasklist != NULL )
  {	if ( join_masks_from_files(mask,sx,sy,inmasklist) )
	 {	fprint_error("unable to read one of the input mask files");
		return(1);
	 }
  }
 outmask=fits_mask_create_empty(sx,sy);

 outimg=fits_create();
 fits_copy_full_header(outimg,inputs[0].img);
 fits_alloc_image(outimg,sx,sy);
 fits_reset_image(outimg);
 outimg->i.curr.bscale=1.0;
 outimg->i.curr.bzero=0.0;
 outimg->i.bit=-32;

 combine_images_from_files(inputs,ninput,outimg,cp,mask,outmask,NULL,0,maxmem);

 if ( fdp.bitpix )              outimg->i.bit=fdp.bitpix;
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

 fits_set_header_integer(outimg,"NCOMBINE",FITS_SH_FIRST,ninput,"number of images used by ficombine");
 fits_set_header_string (outimg,"MCOMBINE",FITS_SH_FIRST,(imagmodestr?imagmodestr:"unspecified"),"combination mode used by ficombine");

 fits_mask_export_as_header(&outimg->header,1,outmask,sx,sy,NULL);
 if ( apply_mask )
  {	for ( i=0 ; i<sy ; i++ )
	 { for ( j=0 ; j<sx ; j++ )
	    {	if ( mask[i][j] )	outimg->i.data[i][j]=0.0;	}
	 }
  }
 
 if ( is_add_history )
	fits_history_export_command_line(outimg,"ficombine",FITSH_FICOMBINE_VERSION,argc,argv);

 if ( outname==NULL )	fw=stdout;
 else			fw=fopenwrite(outname);
 if ( fw==NULL )
  {	fprint_error("unable to create output image");
	return(1);
  }
 fits_write(fw,outimg);
 fclosewrite(fw);
 if ( outmaskname != NULL )
  {	fw=fopenwrite(outmaskname);
	if ( fw==NULL )	
	 {	fprint_error("unable to create output mask file");
		return(1);
	 }
	fits_write_header(fw,outimg);
	fclosewrite(fw);
  }

 fits_free(outimg);
 combine_cleanup(inputs,ninput);

 return(0);
}

/*****************************************************************************/

