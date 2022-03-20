/*****************************************************************************/
/* fic_mpstack.c							     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* A contributed utility to the `fitsh` package aided for faint minor planet */
/* search projects.							     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2004, 2007, 2008, 2010, 2015; Pal, A. (apal@szofi.net)		     */
/*****************************************************************************/
#define	FITSH_FIC_MPSTACK_VERSION	"0.1"
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/time.h>

#include <fits/fits.h>

#include "longhelp.h"
#include "fitsh.h"

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

typedef struct
 {	fits	*img;
	char	**mask;
	double	time;
 } stack;

#define			MODE_MEAN		1
#define			MODE_MEDIAN		2

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
"\t[-i|--input-list <input-list-file>] [-m|--mode {mean|median}]\n"
"\t[-c|--input-coefficients <coeff-file>]\n"
"\t-f|--offset <ox>,<oy> -s|--size <sx>,<sy>\n");
 
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
 
 { NULL, NULL }
};


int fprint_fic_mpstack_long_help(FILE *fw)
{
 fprintf(fw,
"Usage:\tfic_mpstack --input-list <list> --input-coefficients <coefflist>\n"
"The main purpose of this program is to perform a special kind of \n"
"4 dimensional transformation and image stacking that aids faint minor planet \n"
"searches.\n");
 fprintf(fw,"\n");

 longhelp_fprint(fw,fic_mpstack_long_help,0,-1);

 fprintf(fw,"\n");
 fprintf(fw,"Report bugs to <%s>\n",FITSH_MAINT_EMAIL);

 return(0);

}

/*****************************************************************************/


int main(int argc,char *argv[])
{
 FILE		*fr,*fw;

 int		i,is_help;
 char		*inlistfile,*incoefffile,*modestr;
 int		ofx,ofy,nsx,nsy,tsx,tsy;
 int		mode;

 stack		*images;
 int		nimage;

 progbasename=strrchr(argv[0],'/');
 if ( progbasename != NULL )	progbasename++;
 else				progbasename=argv[0];

 inlistfile=incoefffile=NULL;
 is_comment=is_verbose=is_help=0;
 ofx=ofy=nsx=nsy=0;

/*
  input list format: 
	[1] filename
	[2] time offset (usually in days)

*/

/*
  input coefficient format
	[1] Order
	[2] x'
	[3] x''
	[4] y'
	[5] y''
	[6] output filename
*/

 modestr=NULL;

 i=scanarg(argc,argv,SCANARG_ALLOW_FLAGS,
	"--version:%NS-1f%q",&is_help,
        "--version-short|--short-version:%NS-2f%q",&is_help,
        "--long-help|--help-long:%SN2f%q",&is_help,
	"-V|--verbose:%f",&is_verbose,
	"-h|--help:%f%q",&is_help,
	"-i|--input-list:%s",&inlistfile,
	"-c|--input-coefficients:%s",&incoefffile,
	"-f|--offset:%d,%d",&ofx,&ofy,
	"-s|--size:%d,%d",&nsx,&nsy,
	"-m|--mode:%s",&modestr,
	"*:%e",
	NULL);

 if ( i )		
  {	fprint_error("invalid command line argument near '%s'",argv[i]);
	return(1);
  }
 else if ( is_help<0 )
  {	fprint_generic_version(stdout,argv[0],"fic_mpstack",FITSH_FIC_MPSTACK_VERSION,is_help);
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

 if ( modestr ) 
  {	if ( strcmp(modestr,"mean")==0 )
		mode=MODE_MEAN;
	else if ( strcmp(modestr,"median")==0 )
		mode=MODE_MEDIAN;
	else
	  {	fprint_error("invalid mode '%s'",modestr);
		return(1);
	 }
  }
 else
	mode=MODE_MEDIAN;

 if ( inlistfile==NULL || strcmp(inlistfile,"-")==0 )
	fr=stdin;
 else if ( (fr=fopen(inlistfile,"rb"))==NULL )
  {	fprint_error("unable to open input list file '%s'",inlistfile);
	return(1);
  }

 nimage=0;
 images=NULL;

 tsx=tsy=0;

 while ( ! feof(fr) )
  {	char	buff[256],*cmd[16];
	int	n,i;
	FILE	*ff;

	if ( fgets(buff,256,fr)==NULL )
		break;
	buff[255]=0;
	for ( i=0 ; buff[i] ; i++ )
	 {	if ( buff[i]=='#' )	
		 {	buff[i]=0;
			break;
		 }
	 }
	n=tokenize_spaces(buff,cmd,16);
	if ( n<2 )	continue;
	
	ff=fopen(cmd[0],"rb");
	if ( ff==NULL )
	 {	fprint_error("unable to open input file '%s'",cmd[0]);
		return(1);
	 }

	images=(stack *)realloc(images,sizeof(stack)*(nimage+1));
	if ( (images[nimage].img=fits_read(ff))==NULL )
	 {	fprint_error("unable to interpret data from file '%s' as FITS data",cmd[0]);
		return(1);
	 }
	else if ( images[nimage].img->i.dim != 2 )
	 {	fprint_error("image dimension from file '%s' differs from 2",cmd[0]);
		return(1); 
	 }

	if ( nimage==0 )
	 {	tsx=images[nimage].img->i.sx;
		tsy=images[nimage].img->i.sy;
	 }
	else if ( tsx != images[nimage].img->i.sx || tsy != images[nimage].img->i.sy )
	 {	fprint_error("FITS data '%s' has different size",cmd[0]);
		return(1);
	 }

	fits_rescale(images[nimage].img);
	
	if ( sscanf(cmd[1],"%lg",&images[nimage].time)<1 )
	 {	fprint_error("unable to interpret '%s' as time stamp",buff[1]);
		return(1);
	 }

	images[nimage].mask=fits_mask_read_from_header(&images[nimage].img->header,tsx,tsy,NULL);

	fclose(ff);

	nimage++;
	
  }

 fclose(fr);

 if ( ! ( 0<tsx && 0<tsy ) )
 	return(0);
 if ( nsx<=0 || nsy<=0 )
	nsx=tsx,nsy=tsy;

 if ( incoefffile==NULL || strcmp(incoefffile,"-")==0 )
	fr=stdin;
 else if ( (fr=fopen(incoefffile,"rb"))==NULL )
  {	fprint_error("unable to open input coefficient file '%s'",incoefffile);
	return(1);
  }

 while ( ! feof(fr) )
  {	char	buff[4096],*cmd[16],*outfile;
	double	xpoly[16],ypoly[16],*arr,*cxy;
	int	i,j,k,n,order,*idx,*idy;
	fits	*out;
	char	**mask;
	struct	timeval	tv0,tv1;

	if ( fgets(buff,4096,fr)==NULL )
		break;
	buff[4095]=0;
	for ( i=0 ; buff[i] ; i++ )
	 {	if ( buff[i]=='#' )	
		 {	buff[i]=0;
			break;
		 }
	 }
	n=tokenize_spaces(buff,cmd,16);
	if ( n<2 )
		continue;
	sscanf(cmd[0],"%d",&order);
	if ( ! ( 1<=order && order<=4 )	)
		continue;
	if ( 2+order*2 != n )
		continue;
	j=0;
	for ( i=0 ; i<order ; i++ )
	 {	j+=sscanf(cmd[1+0*order+i],"%lg",&xpoly[i]);
		j+=sscanf(cmd[1+1*order+i],"%lg",&ypoly[i]);
	 }
	if ( j<2*order )
		continue;
	outfile=cmd[1+2*order];

	gettimeofday(&tv0,NULL);

	arr=(double *)malloc(sizeof(double)*nimage);
	idx=(int *)malloc(sizeof(int)*nimage);
	idy=(int *)malloc(sizeof(int)*nimage);
	cxy=(double *)malloc(sizeof(double)*nimage*4);

	out=fits_create();
	fits_set_standard(out,NULL);
	out->i.bit=-32;
	out->i.curr.bscale=1.0;
	out->i.curr.bzero=0.0;
	fits_alloc_image(out,nsx,nsy);
	mask=fits_mask_create_empty(nsx,nsy);

	for ( k=0 ; k<nimage ; k++ )
	 {	stack	*s;
		int	o;
		double	t,dx,dy,adx,ady;

		s=&images[k];
		dx=dy=0;
		t=s->time;
		for ( o=0 ; o<order ; o++ )
		 {	dx += xpoly[o]*t;
			dy += ypoly[o]*t;
			t  = t*s->time;
		 }
		idx[k]=(int)floor(dx);
		idy[k]=(int)floor(dy);
		adx=dx-floor(dx);
		ady=dy-floor(dy);
		cxy[4*k+0]=(1-adx)*(1-ady);
		cxy[4*k+1]=( adx )*(1-ady);
		cxy[4*k+2]=(1-adx)*( ady );
		cxy[4*k+3]=( adx )*( ady );
	 }
		
	for ( i=0 ; i<nsy ; i++ )
	 {	for ( j=0 ; j<nsx ; j++ )
		 {	int	np;
			double	sp;

			np=0;
			sp=0.0;
			for ( k=0 ; k<nimage ; k++ )
			 {
				stack	*s;
				double	v;
				int	ix,iy;

				s=&images[k];

				ix=j-ofx+idx[k];
				if ( ! ( 0<=ix && ix<tsx-1 ) )
					continue;
				iy=i-ofy+idy[k];
				if ( ! ( 0<=iy && iy<tsy-1 ) )
					continue;
				if ( s->mask[iy][ix] || s->mask[iy][ix+1] || s->mask[iy+1][ix] || s->mask[iy+1][ix+1] )
					continue;
				
				v=s->img->i.data[iy  ][ix  ]*cxy[4*k+0]+
				  s->img->i.data[iy  ][ix+1]*cxy[4*k+1]+
				  s->img->i.data[iy+1][ix  ]*cxy[4*k+2]+
				  s->img->i.data[iy+1][ix+1]*cxy[4*k+3];

				sp+=v;
				arr[np]=v;
				np++;
			 }
			if ( 0<np )
				out->i.data[i][j]=(mode==MODE_MEAN?sp/np:median(arr,np));
			else
			 {	out->i.data[i][j]=0.0;
				mask[i][j]=MASK_OUTER;
			 }
		 }
	 }
 
	out->i.read.bscale=1.0,
	out->i.read.bzero=0.0;
	fits_set_image_params(out);
	fits_backscale(out,out->i.read.bscale,out->i.read.bzero);
	fits_history_export_command_line(out,"fic_mpstack",FITSH_FIC_MPSTACK_VERSION,argc,argv);
	fits_mask_export_as_header(&out->header,1,mask,nsx,nsy,NULL);
	if ( outfile==NULL )	fw=stdout;
	else			fw=fopenwrite(outfile);
	if ( fw==NULL )		
	 {	fprint_error("unable to create output image file '%s'",outfile);
		return(1);
	 }
	fits_write(fw,out);
	fclosewrite(fw);
	
	fits_mask_free(mask);
	fits_free(out);

	free(cxy);
	free(idy);
	free(idx);
	free(arr);

	gettimeofday(&tv1,NULL);
	tv1.tv_usec -= tv0.tv_usec;
	tv1.tv_sec  -= tv0.tv_sec;
	if ( tv1.tv_usec < 0 )
	 {	tv1.tv_usec += 1000000;
		tv1.tv_sec  -= 1;
	 }
	if ( is_verbose )
		fprintf(stderr,"%s: done: %d.%.6d\n",outfile,(int)tv1.tv_sec,(int)tv1.tv_usec);
  }


 return(0);
}

/*****************************************************************************/

