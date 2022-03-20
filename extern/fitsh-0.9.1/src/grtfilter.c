/*****************************************************************************/
/* grtfilter.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* An implementation for the parameter decorrelation and trend filtering     */
/* algorithm.								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2008; Pal, A. (apal@szofi.elte.hu); based on 2005, MNRAS, 356, 557    */
/*****************************************************************************/
#define	FI_GRTFILTER_VERSION	"0.1"
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include "fi.h"

#include "io/iof.h"
#include "io/scanarg.h"
#include "io/tokenize.h"

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

typedef struct
 {	char		*key;
	char		*line;
	double		mag;
	double		err;
 } lcpoint;

typedef struct
 {	lcpoint		*lcpoints;
	int		nlcpoint;
	double		mean;
 } lightcurve;

typedef struct
 {	lightcurve	*templatecurves;
	int		ntemplatecurve;
	double		**fmatrix;
	double		**imatrix;
 } tfatemplate;

typedef struct
 {	int		flt_niter;
	double		flt_sigma;
	int		ilc_niter;
	double		ilc_sigma;
 } filterparam;

typedef struct
 {	int		o_xy;
	int		o_fxfy;
	int		o_sdk;
	int		o_hr;
	int		o_am;
	int		minfreedom;
 } epdmodel;

typedef struct
 {	int		x,y;
	int		s,d,k;
	int		hr,am;
 } epdcolumn;

/*****************************************************************************/

int fprint_grtfilter_usage(FILE *fw)
{
 fprintf(fw,
"Usage:\tgrtfilter [-h|--help] [-V|--verbose]\n"
"\t[-t [@]<templates>] {-l <list> | [-i <input> ...]  [-o <output> ...]}\n"
"\t[--write-epd]\n");
 fprintf(fw,
"Common options:\n"
"\t[--col-[input-]key <>[,<>,...]]\n"
"\t[--col-xy <>,<>] [--col-hram <>,<>] [--col-sdk <>,<>,<>]\n"
"\t[--col-{mag|err|model} <>[,<>[,...]]] [--model {mean|median|user}]\n");
 fprintf(fw,
"External parameter decorrelation tuning:\n"
"\t[--epd [xy=<order>,fxfy=<order>,ha=<order>,am=<order>,sdk=<order>]\n"
"\t       [freedom=<minimum-degrees-of-freedom>]]\n"
"\t[--distinct-blocks <pattern1>[,<pattern2>,[...]]]\n");
 fprintf(fw,
"Filtering:\n"
"\t[--{template|input}-filter [iterations=<n>,sigma=<sigma]]\n");
 return(0);
}

int fprint_grtfilter_long_help(FILE *fw)
{
 return(0);
}

/*****************************************************************************/

int extract_more_columns(char *colstr,int **rcols)
{
 int	ncol,*cols,c; 
 char	*dup,**cmd;

 dup=strdup(colstr);
 cmd=tokenize_char_dyn(dup,','); 
 cols=NULL;
 for ( ncol=0 ; cmd[ncol] != NULL ; ncol++ )
  {	cols=(int *)realloc(cols,sizeof(int)*(ncol+1));
	if ( sscanf(cmd[ncol],"%d",&c)<1 )
	 {	free(cmd);free(dup);
		return(-1);
	 }
	cols[ncol]=c;
  }
 free(cmd);
 free(dup);

 if ( rcols != NULL )	*rcols=cols;

 return(ncol);
}

int main(int argc,char *argv[])
{
 int		i,is_help;
 char		*fltstr,*ilcstr;
 char		**templatelist,**inputlist,**outputlist;
 char		*listfile;
 char		*epdmodelstr,*colmagstr,*colerrstr,*colmodstr;
 epdcolumn	ecol;
 epdmodel	emod;

 int		nmag;	/* number of magnitudes to apply epd/tfa */ 
 int		*cmags,*cerrs,*cmods;	
 int		ctkey,cikey;

 filterparam	fp;

 progbasename=strrchr(argv[0],'/');
 if ( progbasename != NULL )	progbasename++;
 else				progbasename=argv[0];

 fltstr=NULL;
 ilcstr=NULL;

 templatelist=NULL;
 inputlist=NULL;
 outputlist=NULL;
 listfile=NULL;

 epdmodelstr=NULL;
 ecol.x=ecol.y=-1;
 ecol.s=ecol.d=ecol.k=-1;
 ecol.hr=ecol.am=-1; 

 ctkey=1;
 cikey=1;
 colmagstr="2";
 colerrstr=NULL;
 colmodstr=NULL;
 
 is_help=0;

 i=scanarg(argc,argv,SCANARG_ALLOW_FLAGS,
	"--version:%NS-1f%q",&is_help,
        "--version-short|--short-version:%NS-2f%q",&is_help,
	"-h|--help|--short-help|--help-short:%f%q",&is_help,
	"--long-help|--help-long:%SN2f%q",&is_help,
	"-i|--input|--inputs:%Mt",&inputlist,
	"-t|--template|--templates:%Mt",&templatelist,
	"-o|--output|--outputs:%Mt",&outputlist,
	"-l|--list-file|--input-list-file:%s",&listfile,
	"--col-xy:%Pd,%Pd",&ecol.x,&ecol.y,
	"--col-hram:%Pd,%Pd",&ecol.hr,&ecol.am,
	"--col-hr:%Pd",&ecol.hr,
	"--col-sdk:%Pd,%Pd,%Pd",&ecol.s,&ecol.d,&ecol.k,
	"--col-mag|--col-magnitude|--col-magnitudes:%s",&colmagstr,
	"--col-err|--col-magnitude-error:%s",&colerrstr,
	"--col-model:%s",&colmodstr,
	"--col-key|--col-template-key:%Pd",&ctkey,
	"--col-input-key:%Pd",&cikey,
	"--epd|--epd-model:%s",&epdmodelstr,
	
	"--template-filter:%s",&fltstr,
	"--input-filter:%s",&ilcstr,

	NULL);

 if ( i )
  {	fprint_error("invalid command line argument near '%s'",argv[i]);
	return(1);
  }
 else if ( is_help<0 )
  {	fprint_generic_version(stdout,argv[0],"grtfilter",FI_GRTFILTER_VERSION,is_help);
	return(0);
  } 
 else if ( is_help>1 )
  {	fprint_grtfilter_long_help(stdout);
	return(0);
  }
 else if ( is_help )
  {	fprint_grtfilter_usage(stdout);
	return(0);
  }

 if ( fltstr != NULL )
  {	i=scanpar(fltstr,SCANPAR_DEFAULT,
		"iterations|niter:%d" ,&fp.flt_niter,
		"sigma|limit|level:%g",&fp.flt_sigma,
		NULL);
	if ( i )
	 {	fprint_error("invalid parameter in '%s'.\n",fltstr);
		return(1);
	 }
  }
 else
  {	fp.flt_niter=0;
	fp.flt_sigma=3.0;
  }
 if ( ilcstr != NULL )
  {	i=scanpar(ilcstr,SCANPAR_DEFAULT,
		"iterations|niter:%d" ,&fp.ilc_niter,
		"sigma|limit|level:%g",&fp.ilc_sigma,
		NULL);
	if ( i )
	 {	fprint_error("invalid parameter in '%s'.\n",ilcstr);
		return(1);
	 }
  }
 else
  {	fp.ilc_niter=0;
	fp.ilc_sigma=3.0;
  }

 emod.o_xy=0;
 emod.o_fxfy=0;
 emod.o_hr=0;
 emod.o_am=0;
 emod.o_sdk=0;
 emod.minfreedom=-1;

 if ( epdmodelstr != NULL )
  {	i=scanpar(epdmodelstr,SCANPAR_DEFAULT,
		"xy:%d",&emod.o_xy,
		"fxfy:%d",&emod.o_fxfy,
		"hr:%d",&emod.o_hr,
		"am:%d",&emod.o_am,
		"sdk:%d",&emod.o_sdk,
		"freedom:%d",&emod.minfreedom,
		NULL);
	if ( i )
	 {	fprint_error("invalid parameter in '%s'.\n",ilcstr);
		return(1);
	 }
  }

 nmag=extract_more_columns(colmagstr,&cmags);
 if ( nmag<=0 )
  {	fprint_error("unable to parse magnitude column declaration argument '%s'.\n",colmagstr);
	return(1);
  }

 if ( colerrstr != NULL )
  {	i=extract_more_columns(colerrstr,&cerrs);
	if ( i<=0 )
	 {	fprint_error("unable to parse magnitude error column declaration argument '%s'.\n",colmagstr);
		return(1);
	 }
	else if ( i != nmag )
	 {	fprint_error("number of magnitude and error columns differ.\n");
		return(1);
	 }
  }
 else
	cerrs=NULL;

 if ( colmodstr != NULL )
  {	i=extract_more_columns(colmodstr,&cmods);
	if ( i<=0 )
	 {	fprint_error("unable to parse magnitude model column declaration argument '%s'.\n",colmagstr);
		return(1);
	 }
	else if ( i != nmag )
	 {	fprint_error("number of magnitude and model columns differ.\n");
		return(1);
	 }
  }
 else
	cmods=NULL;



 return(0);
}
