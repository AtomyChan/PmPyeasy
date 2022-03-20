/*****************************************************************************/
/* lfit.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Symbolic fitting & arithmetic evaluating utility. Supported (fit) methods:*/
/*  - NONE: simple evaluation of a set of functions			     */
/*  - CLLS: classical linear least squares method			     */
/*  - NLLM: nonlinear Levenberg-Marquardt least squares methods 	     */
/*  - MCMC: Markov chain Monte-Carlo method 				     */
/*  - MCHI: mapping of values of chi2 on a grid 			     */
/*  - EMCE: error estimation based on synthetic Monte-Carlo re-fitting	     */
/*  - DHSX: downhill simplex minimalization 				     */
/*  - XMMC: extended Markov chain Monte-Carlo method 			     */
/*  - LMND: Levenberg-Marquardt with numerical approximation for derivatives */
/*  - FIMA: Fisher Information Matrix Analysis (not really a fit method)     */
/* The set of built-in functions are definded in the `lfit-builtin.[ch]`     */
/* modules. Additional functions can easily be registered by using	     */
/* the function lfit_register_function() (see _FI_SOURCE part for example).  */
/* Further extensions are also available runtime, using dynamically loaded   */
/* external shared libraries (see -d|--dynamic option or `linear.c` for ex.) */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 1996, 2002, 2004-2005, 2006, 2007-2008, 2009; 			     */
/* Pal, Andras (apal@szofi.net)						     */
/*****************************************************************************/

#define	LFIT_VERSION		"7.9"
#define	LFIT_LASTCHANGE		"2009.05.25"

/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>

#include <time.h>
#include <sys/time.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include <fcntl.h>

#include "longhelp.h"

#include <lfit/lfit.h>

#include "lfit-info.h"
#include "lfit-builtin.h"

#ifndef	LFIT_ENABLE_DYNAMIC_EXTENSIONS
#define	LFIT_ENABLE_DYNAMIC_EXTENSIONS		1
#endif

#ifdef	LFIT_ENABLE_DYNAMIC_EXTENSIONS
#include <dlfcn.h>
#endif

#if defined	_CHDET_SOURCE
	#include "psn.h"
	#include "str.h"
	#include "iof.h"
	#include "tokenize.h"
	#include "scanarg.h"
	#include "lmfit.h"
	#include "pbfft.h"
	#include "downhill.h"
#elif defined	_FI_SOURCE
	#include "str.h"
	#include <psn/psn.h>
	#include "io/iof.h"
	#include "io/tokenize.h"
	#include "io/scanarg.h"
	#include "math/fit/lmfit.h"
	#include "math/dft/pbfft.h"
	#include "math/fit/downhill.h"
#else
	#include <psn/psn.h>
	#include "str.h"
	#include "iof.h"
	#include "tokenize.h"
	#include "scanarg.h"
	#include "lmfit.h"
	#include "pbfft.h"
	#include "downhill.h"
#endif

/*****************************************************************************/

#define		LFIT_DEFAULT_FORMAT		"%12g"
#define		LFIT_DEFAULT_CORR_FORMAT	"%6.3f"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define		RESIDUAL_BIASED		0
#define		RESIDUAL_UNBIASED	1
#define		RESIDUAL_CHI2		2
#define		RESIDUAL_BIASED_CHI2	2
#define		RESIDUAL_UNBIASED_CHI2	3

/*****************************************************************************/

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841968
#endif

/*****************************************************************************/

typedef struct
 {	char	*name;
	psn	*pmacro;
	int	major;
	int	argnum;
 } psnmacro;

typedef struct
 {	lfitfunction	*lff;
	int		fmajor;
	int		nmajor;
	double		*diff;
 } lffreg;

typedef struct
 {	psnsym		* pl_sym;
	int		nsym;
	psnprop		* pl_prop;
	int		nprop;
	psnfunct	* pl_funct;
	int		nfunct;
	psndiff		* pl_diff;
	int		ndiff;
	psnsymeval	* pl_symeval;
	int		nsymeval;
	psnmacro	* pl_macro;
	char		**symbols;
	int		nsymbol;
	short		**pl_simp;
	lffreg		*lffregs;
	int		nlffreg;
 } lfitpsnglobal;

lfitpsnglobal	lpg = 
 {	NULL,
	0,
	NULL,
	0,
	NULL,
	0,
	NULL,
	0,
	NULL,
	0,
	NULL,
	NULL,
	0,
	psn_lfit_simp,
	NULL,
	0
 };

/*****************************************************************************/

int error_vpush(char **rstr,char *prefix,char *fmt,va_list ip)
{
 int		n,l,p,size;
 char		*str,*nstr;
 /*  va_list	ap; */

 if ( rstr==NULL )	
	return(1);

 str=NULL;
 l=0;
 size=128;
 str=realloc(str,l+size);
 if ( str==NULL )      return(-1);
 while ( 1 )
  {     /* va_copy(ap,ip); */
	n=vsnprintf((str)+l,size,fmt,ip);
	/* va_end(ap); */
        if ( n>-1 && n<size )
                break;
        else if ( n>-1 )
                size=n+1;       
        else
                size=size*2;    
        if ( (nstr=realloc(str,l+size))==NULL )
         {	free(str);
		return(-1);
	 }
	else
		str=nstr;
  };

 if ( *rstr==NULL )
	*rstr=str;
 else
  {	l=strlen(*rstr);
	if ( prefix==NULL )	p=0; 
	else			p=strlen(prefix);
	*rstr=realloc((*rstr),l+p+n+1);
	memmove((*rstr)+n+p,(*rstr),l+1);
	memcpy((*rstr),str,n);
	if ( p>0 )	memcpy((*rstr)+n,prefix,p);
	free(str);
  }
	
 return(0);    
}

int error_push(char **rstr,char *prefix,char *fmt,...)
{
 va_list	ip;
 int		r;

 va_start(ip,fmt);
 r=error_vpush(rstr,prefix,fmt,ip);
 va_end(ip);

 return(r);
}

int error_free(char **rstr)
{
 if ( rstr==NULL )
	return(-1);
 else if ( *rstr==NULL )
	return(1);
 else
  {	free(*rstr);
	return(0);
  }
}

int lfit_error_push(char **rstr,char *fmt,...)
{
 va_list	ip;
 int		r;

 va_start(ip,fmt);
 r=error_vpush(rstr,": ",fmt,ip);
 va_end(ip);

 return(r);
}

int lfit_error_fprint_error(char **rstr)
{
 lfit_error_push(rstr,"error");
 lfit_error_push(rstr,"lfit");
 fprintf(stderr,"%s.\n",*rstr);
 return(0);
}

int lfit_error_fprint_warning(char **rstr)
{
 lfit_error_push(rstr,"error");
 lfit_error_push(rstr,"lfit");
 fprintf(stderr,"%s.\n",*rstr);
 return(0);
}

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

char *errmsgs[]=
 {	"",	/* dummy */
/* 1 */	"invalid command line argument (see --help for help)",
/* 2 */	"unable to open input file",
/* 3 */	"unable to create output file",
/* 4 */	"definition of fit variables missing, use -v 'vars' to define",
/* 5 */	"column definitions missing, use -c 'cols' to define",
/* 6 */	"fit function is missing, use -f 'func' to define one",
/* 7 */	"invalid function (symbolic error)",
/* 8 */	"invalid function (parse error)",
/* 9 */	"fatal - inconsistency in the psn library!",
/*10 */	"too few lines for fitting",
/*11 */	"singular matrix",
/*12 */	"invalid dependent expression (symbolic error)",
/*13 */	"invalid dependent expression (parse error)",
/*14 */	"fit function is non-linear in one of the fitting parameters",
/*15 */	"dependent value missing or ambigous, use -y 'expr' or -f '...=expr' to define",
/*16 */	"non-differentable operator in the fitting function",
/*17 */ "unable to create output list file for fitted lines",
/*18 */ "unable to create output list file for rejected lines",
/*19 */ "invalid, singular or non-linear constraint",
/*20 */ "constraints can be used only with linear fits",
/*21 */ "too many constraints",
/*22 */ "invalid column specification (note that each column can only be defined once)",
/*23 */ "none of the fit or dependent expressions contain variables to be fitted",
/*24 */ "both fit and dependent expression contain variables to be fitted",
/*25 */ "unable to create output list file",
/*26 */ "unable to create variable list file",
/*27 */ "unable to create MCMC output file",
/*28 */	"invalid combination of input fit data",
/*29 */	"both error and weight have been defined, use only one of them",
/*30 */ "input file name is missing, use -i<key> 'filename' to define",
/*31 */ "some of the data blocks are used for fitting, some of them are not - be consequent",
/*32 */ "invalid chi-weight format or value",
/*33 */ "invalid variable format string",
/*34 */	"invalid string for variable differences",
	NULL
 };

/*
int effective_strlen(char *s)
{
 int	r;
 for ( r=0 ; *s != 0 ; s++ )
  {	if ( *s == '\t' )
		r=(r+8)&(~7);
	else
		r++;
  }
 return(r);
}
*/

int fprint_lfit_function_description(FILE *fw,int width,int w,int pad,char *desc)
{
 int	i;

 for ( i=w ; i<pad ; i++ )
  {	fprintf(fw," ");		}

 if ( width<=0 )
	fprintf(fw,"%s\n",desc);
 else
  {	char	*dd,**cmd;
	int	n,w,p,l;
	dd=strdup(desc);
	cmd=tokenize_spaces_dyn(dd);
	p=pad;w=0;
	for ( n=0 ; cmd != NULL && cmd[n] != NULL ; n++ )
	 {	l=strlen(cmd[n]);
		if ( w>0 )	l++;
		if ( l>width-p )
		 {	fprintf(fw,"\n");
			for ( i=0 ; i<pad ; i++ )
			 {	fprintf(fw," ");	}
			if ( w>0 )	l--;
			w=0;
			p=pad;
		 }
		if ( w>0 )
			fprintf(fw," %s",cmd[n]);
		else
			fprintf(fw,"%s",cmd[n]);
		p+=l;
		w++;
	 };
	fprintf(fw,"\n");
	if ( cmd != NULL )	free(cmd);
	free(dd);
  }

 return(0);
}

int fprint_lfit_function_info(FILE *fw,int width,psnlfit *pl)
{
 int	w;

 if ( pl->info==NULL )
	return(1);
 else if ( pl->type==T_OP )
  {	if ( pl->minor==TO_INFIX )
	 {	fprintf(fw,"(.)%s(.)",pl->name);
		w=3+strlen(pl->name)+3;
	 }
	else if ( pl->minor==TO_PREFIX )
	 {	fprintf(fw,"%s(.)",pl->name);	
		w=strlen(pl->name)+3;
	 }
	else if ( pl->minor==TO_SUFFIX )
	 {	fprintf(fw,"(.)%s",pl->name);	
		w=3+strlen(pl->name);
	 }
	else
		return(1);
	if ( pl->info != NULL && pl->info->description != NULL )
		fprint_lfit_function_description(fw,width,w,16,pl->info->description);
	else
		fprintf(fw,"\n");
	return(0);
  }
 else if ( pl->type==T_FN )
  {	int	i;
	fprintf(fw,"%s(",pl->name);
	w=strlen(pl->name)+1;
	for ( i=0 ; i<pl->argnum ; i++ )
	 {	if ( i>0 )
		 {	fprintf(fw,",");
			w++;
		 }
		fprintf(fw,".");
		w++;
	 }
	fprintf(fw,")");
	w++;
	if ( pl->info != NULL && pl->info->description != NULL )
		fprint_lfit_function_description(fw,width,w,16,pl->info->description);
	else
		fprintf(fw,"\n");
	return(0);
  }

 else
	return(1);
}

int fprint_lfit_function_sublist(FILE *fw,int width,psnlfit *pl)
{
 for ( ; pl->name != NULL ; pl++ )
  { 	fprint_lfit_function_info(fw,width,pl);		}
 return(0);
}

int fprint_lfit_function_list(FILE *fw,...)
{
 va_list	ap;
 int		width;
 psnlfit	*pl;

#ifdef	TIOCGWINSZ
 if ( isatty(fileno(fw)) )
  {	struct	winsize	ws;
	if ( ! ioctl(fileno(fw),TIOCGWINSZ,&ws) )
		width=ws.ws_col-1;
	else
		width=0;
  }	
 else
#endif
	width=0;

 fprintf(fw,"# List of built-in operators and functions\n");
 fprintf(fw,"# Name\t\tDescription\n");
 fprintf(fw,"# [1]\t\t[2]\n");

 va_start(ap,fw);
 while ( (pl=va_arg(ap,psnlfit *)) != NULL )
  {	fprint_lfit_function_sublist(fw,width,pl);		}
 va_end(ap);

 return(0);
}

void exit_with_usage(int t)
{
 if ( t )	fprint_error("%s",errmsgs[t]);
 else		fprint_lfit_usage(stderr);
 exit(t);
}

/*****************************************************************************/

int lfit_register_symbol_exists(lfitpsnglobal *lpg,char *name)
{
 char	**s;
 int	n;
 for ( s=lpg->symbols,n=lpg->nsymbol ; n>0 && s != NULL && (*s) != NULL ; s++,n-- )
  {	if ( strcmp((*s),name)==0 )
		return(1);
  }
 return(0);
}
int lfit_register_symbol_add(lfitpsnglobal *lpg,char *name)
{
 lpg->symbols=(char **)realloc(lpg->symbols,sizeof(char *)*(lpg->nsymbol+1));
 lpg->symbols[lpg->nsymbol]=name;
 lpg->nsymbol++;
 return(0);
}

int lfit_register_check_code(lfitpsnglobal *lpg,int major)
{
 psnprop	*p;
 int		n;
 for ( p=lpg->pl_prop,n=lpg->nprop ; n>0 && p != NULL && p->major>=0 ; p++,n-- )
  {	if ( p->major==major )
		return(1);
  }
 return(0);
}

int lfit_register_function(char *name,int major,int argnum,
	int   (*funct)(double *top_of_the_stack_pointer),
	short * diffrule,
	char  * symevalstring)
{
 if ( lfit_register_check_code(&lpg,major) )
  {	fprint_error("internal: registering function '%s': major code %d aready in use",name,major);
	return(1);
  }

 if ( major>=0 && name != NULL )
  {	lpg.pl_sym=(psnsym *)realloc(lpg.pl_sym,sizeof(psnsym)*(lpg.nsym+2));
	lpg.pl_sym[lpg.nsym].type=T_FN;
	lpg.pl_sym[lpg.nsym].major=major;
	lpg.pl_sym[lpg.nsym].name=name;
	lpg.pl_sym[lpg.nsym].minor=argnum;
	memset(&lpg.pl_sym[lpg.nsym+1],0,sizeof(psnsym));
	lpg.nsym++;
	lfit_register_symbol_add(&lpg,name);
  }

 if ( funct != NULL )
  {	lpg.pl_prop=(psnprop *)realloc(lpg.pl_prop,sizeof(psnprop)*(lpg.nprop+2));
	lpg.pl_prop[lpg.nprop].major=major;
	lpg.pl_prop[lpg.nprop].argnum=argnum;
	lpg.pl_prop[lpg.nprop].precedency=0;
	lpg.pl_prop[lpg.nprop].associativity=0;
	memset(&lpg.pl_prop[lpg.nprop+1],0,sizeof(psnprop));
	lpg.nprop++;

	lpg.pl_funct=(psnfunct *)realloc(lpg.pl_funct,sizeof(psnfunct)*(lpg.nfunct+2));
	lpg.pl_funct[lpg.nfunct].major=major;
	lpg.pl_funct[lpg.nfunct].funct=funct;
	memset(&lpg.pl_funct[lpg.nfunct+1],0,sizeof(psnfunct));
	lpg.nfunct++;
  }

 if ( diffrule != NULL )
  {	lpg.pl_diff=(psndiff *)realloc(lpg.pl_diff,sizeof(psndiff)*(lpg.ndiff+2));
	lpg.pl_diff[lpg.ndiff].major=major;
	lpg.pl_diff[lpg.ndiff].oplist=diffrule;
	memset(&lpg.pl_diff[lpg.ndiff+1],0,sizeof(psndiff));
	lpg.ndiff++;
  }

 if ( symevalstring != NULL )
  {	lpg.pl_symeval=(psnsymeval *)realloc(lpg.pl_symeval,sizeof(psnsymeval)*(lpg.nsymeval+2));
	lpg.pl_symeval[lpg.nsymeval].major=major;
	lpg.pl_symeval[lpg.nsymeval].string=symevalstring;
	lpg.pl_symeval[lpg.nsymeval].strength=argnum;
	lpg.pl_symeval[lpg.nsymeval].affixation=0;
	memset(&lpg.pl_symeval[lpg.nsymeval+1],0,sizeof(psnsymeval));
	lpg.nsymeval++;
  }

 return(0);
}

int lfit_register_operator(char *name,int major,int minor,int argnum,
	int   precedency,int associativity,
	int   (*funct)(double *top_of_the_stack_pointer),
	short * diffrule,
	char  * symevalstring,int strength,int affixation)
{
 if ( lfit_register_check_code(&lpg,major) )
  {	fprint_error("internal: registering operator  '%s': major code %d aready in use",name,major);
	return(1);
  }

 if ( major>=0 && name != NULL )
  {	lpg.pl_sym=(psnsym *)realloc(lpg.pl_sym,sizeof(psnsym)*(lpg.nsym+2));
	lpg.pl_sym[lpg.nsym].name=name;
	lpg.pl_sym[lpg.nsym].type=T_OP;
	lpg.pl_sym[lpg.nsym].major=major;
	lpg.pl_sym[lpg.nsym].minor=minor;
	memset(&lpg.pl_sym[lpg.nsym+1],0,sizeof(psnsym));
	lpg.nsym++;
	lfit_register_symbol_add(&lpg,name);
  }

 if ( funct != NULL )
  {	lpg.pl_prop=(psnprop *)realloc(lpg.pl_prop,sizeof(psnprop)*(lpg.nprop+2));
	lpg.pl_prop[lpg.nprop].major=major;
	lpg.pl_prop[lpg.nprop].argnum=argnum;
	lpg.pl_prop[lpg.nprop].precedency=precedency;
	lpg.pl_prop[lpg.nprop].associativity=associativity;
	memset(&lpg.pl_prop[lpg.nprop+1],0,sizeof(psnprop));
	lpg.nprop++;

	lpg.pl_funct=(psnfunct *)realloc(lpg.pl_funct,sizeof(psnfunct)*(lpg.nfunct+2));
	lpg.pl_funct[lpg.nfunct].major=major;
	lpg.pl_funct[lpg.nfunct].funct=funct;
	memset(&lpg.pl_funct[lpg.nfunct+1],0,sizeof(psnfunct));
	lpg.nfunct++;
  }

 if ( diffrule != NULL )
  {	lpg.pl_diff=(psndiff *)realloc(lpg.pl_diff,sizeof(psndiff)*(lpg.ndiff+2));
	lpg.pl_diff[lpg.ndiff].major=major;
	lpg.pl_diff[lpg.ndiff].oplist=diffrule;
	memset(&lpg.pl_diff[lpg.ndiff+1],0,sizeof(psndiff));
	lpg.ndiff++;
  }

 if ( symevalstring != NULL )
  {	lpg.pl_symeval=(psnsymeval *)realloc(lpg.pl_symeval,sizeof(psnsymeval)*(lpg.nsymeval+2));
 	lpg.pl_symeval[lpg.nsymeval].major=major;
 	lpg.pl_symeval[lpg.nsymeval].string=symevalstring;
 	lpg.pl_symeval[lpg.nsymeval].strength=strength;
	lpg.pl_symeval[lpg.nsymeval].affixation=affixation;
	memset(&lpg.pl_symeval[lpg.nsymeval+1],0,sizeof(psnsymeval));
  }

 return(0);
}

int lfit_register_internal(psnlfit *pl)
{
 int	ret;

 if ( pl->type==T_OP )
  {	ret=lfit_register_operator(pl->name,pl->major,pl->minor,pl->argnum,
		pl->precedency,pl->associativity,
		pl->funct,
		pl->diff,
		pl->string,pl->strength,pl->affixation);
  }
 else if ( pl->type==T_FN )
  {	ret=lfit_register_function(pl->name,pl->major,pl->argnum,
		pl->funct,
		pl->diff,
		pl->string);
  }
 else
	ret=2;

 return(ret);
}
		
/*****************************************************************************/

int	is_verbose;

/*****************************************************************************/

typedef struct
 {	FILE	*f_write;	/* [-o ]  default output		*/
	FILE	*f_lused;	/* [-of] used lines			*/
	FILE	*f_lrejd;	/* [-or] rejected lines			*/
	FILE	*f_lall;	/* [-oa] all lines			*/
	FILE	*f_expr;	/* [-ox] expression			*/
	FILE	*f_vval;	/* [-ov] variable=value args.		*/
 } fitout;

#define		FIT_METHOD_NONE	0	/* no fitting, just evaluation	     */
#define		FIT_METHOD_CLLS	1	/* classical linear least squares    */
#define		FIT_METHOD_NLLM	2	/* nonlinear Levenberg-Marquardt     */
#define		FIT_METHOD_MCMC	3	/* Markov Chain Monte-Carlo	     */
#define		FIT_METHOD_MCHI	4	/* mapping of chi^2 space on a grid  */
#define		FIT_METHOD_EMCE	5	/* Monte-Carlo error estimation	     */
#define		FIT_METHOD_DHSX	6	/* downhill simplex	 	     */
#define		FIT_METHOD_XMMC	7	/* extended Markov Chain Monte-Carlo */
#define		FIT_METHOD_LMND	8	/* L-M with numerical derivatives    */
#define		FIT_METHOD_FIMA 9	/* Fisher information matrix analysis*/

typedef struct
 {	int	fit_method;	/* see FIT_METHOD_* definitions above	     */
	int	niter;		/* number of sigmareject iterations	     */
	double	sigma;		/* rejection level in the units of std dev   */
	int	mc_iterations;	/* EMCE/MCMC/XMMC iterations		     */

	struct
	 {	struct
		 {	int	use_fisher_sx;	/* use inital Fisher simplex */
		 } dhsx;
		struct
		 {	int	sub_method;	/* fit method of the substeps*/
			int	skip_initial_fit; /* skip initial EMCE fit   */
		 } emce;
		struct
		 {	double	lambda,		/* initial \lambda 	     */
				lambda_mpy;	/* \lambda multiplier...     */
			int	max_iter;	/* max LM iterations	     */
			int	numeric_derivs;	/* LMND vs. NLLM	     */
		 } nllm;
		struct
		 {	int	do_montecarlo;	/* do a MC dataset	     */
			int	write_original;	/* write original values     */
			int	write_uncert;	/* write uncertainties	     */
			int	write_correl;	/* write correlation matrix  */
		 } fima;
		struct
		 {	int	use_gibbs;	/* do use Gibbs sampler */
			int	count_accepted;	/* count accepted transitions*/
		 } mcmc;
		struct
		 {	int	skip_initial_fit; /* skip initial (dhsx) fit */
			int	is_adaptive;	/* use adaptive XMMC         */
			int	count_accepted;	/* count accepted transitions*/
			int	window;		/* autocorr. length window   */
			int	niter;		/* additional iterations     */
		 } xmmc;
	 } tune;	/* this is not an union since some complex methods */
			/* require more parameter sets (e.g. EMCE+NLLM)    */
	
 } fitparam;

typedef struct
 {	int		nc;		/* total number of linear constraints*/
	double		**cmatrix;	/* [nc]x[nvar] matrix of lin. constrs*/
	double		*cvector;	/* [nvar] vector of constraint values*/
	double		**invlmatrix;	/* inverse of the [nc+nvar]x[nc+nvar]*/
					/* matrix (Lagrange-multiplicators)  */
	double		**cproject;	/* constraint projection matrix	     */
 } fitconstraint;

typedef struct
 {	psn		*funct;		/* function describing the constraint*/
	int		ctype;		/* type				     */
 } dcfunct;

typedef struct
 {	int		ndcfunct;	/* number of domain constraints	     */
	dcfunct		*dcfuncts;
 } domconstraint;

#define		VAR_IS_CONSTANT	0x01	/* is it only a constant? 	     */
#define		VAR_MIN		0x02	/* MCMC minimal value is set & used  */
#define		VAR_MAX		0x04	/* MCMC maximal value is set & used  */
#define		VAR_STEP	0x08	/* chi2 mapping grid resolution	     */
#define		VAR_ERR		0x10	/* has an initial error be defined?  */

typedef struct
 {	char	*name;		/* name of the variable to be fitted	     */
	char	format[16];	/* (*)printf-like output format		     */
	double	init;		/* initial value (used by NLLM, DHSX, *MC)   */
	double	ierr;		/* initial error (used by MCMC, DHSX)	     */
	double	imin,imax;	/* allowed range: minimal and maximal values */
	double	istp;		/* grid step of the chi2 mapping	     */
	double	diff;		/* difference used by LMND (or EMCE+LMND)    */
	int	flags;		/* some flags (see VAR_*)		     */
	int	is_linear;	/* all of the functions linear in this var...*/
	int	is_separated;	/* ..., and if so, separate this variable    */
 } variable;

typedef struct
 {	char	*name;		/* name of the derived variable		     */
	char	format[16];	/* (*)printf-like output format		     */
	char	*expr;		/* expression for definition function	     */
	psn	*funct;		/* definition function			     */
	psn	**diff;		/* partial derivatives with respect to 	     */
				/* the (fitted) variables		     */
 } dvariable;

typedef struct
 {	char	*name;		/* the name of the given column		     */
 } column;

typedef struct
 {	char	format[16];	/* *printf like output format		     */
 } dumpexpr;

typedef struct
 {	int	dbidx;		/* data block index respecting to this row   */
	char	*line;		/* orig. input line (as read from the file)  */
	double	*x;		/* pointer to the numeric column data	     */
 } fitinputrow;	

typedef struct
 {	char	*key;		/* <key> of the given data block	     */

	char	*colarg,	/* -c<key> ...				     */
		*fncarg,	/* -f<key> ...				     */
		*deparg,	/* -y<key> ...				     */
		*errarg,	/* -e<key> ... or -w<key> ...		     */
		*inparg;	/* -i<key> ...				     */
	int	errtype;	/* -e (error, 0) or -w (weight, !0) ?	     */
	
	int	ncol;		/* number of columns (derived from colarg)   */
	column	*cols;		/* column specifications (i.e. names)	     */

	psnsym	*colsym;	/* column symbols for the given data block   */

	psn	*funct,		/* PSN sequence of the function to be fitted */
		*functorig,	/* original (affine but nonlinear) version   */
		**diff;		/* partial derivatives of the function	     */

	lfitfunction
		*lfunct;
	int	*map_var;
	int	*map_idv;

	psn	*dep,		/* PSN sequence of dependent expression (-y) */
		*err;		/* PSN sequence of the error/weight expr.    */
	int	*fchain;	/* argument chain for 'funct'		     */

	int	is_linear;	/* is 'funct' linear in the fit variables?   */

	int	offset;		/* offset of the related data in the arraies */
	int	size;		/* number of related data rows...	     */

	double	dbscatter;	/* scatter (RMS) respecting to the datablock */
	double	dbredchi2;	/* reduced chi2 value resp to the datablk    */
	double	xnoise;		/* extra (red?) noise added by the EMCE	     */
 } datablock;

typedef struct
 {	char		corrfm[16];	/* correlation format (printf-like)  */
	psnsym		*varsym;	/* variable symbol information	     */
	datablock	*datablocks;	/* individual data blocks [DBs]	     */
	int		ndatablock;	/* number of data blocks 	     */
	int		is_linear;	/* are all of the functions linear?  */
	int		maxncol;	/* max number of columns of each DBs */
	fitparam	parameters;	/* general fit parameters	     */
	fitconstraint	fconstraint;	/* fit constraints		     */
	domconstraint	dconstraint;	/* domain constraints (for MC)	     */
 } lfitdata;

typedef struct
 {	int		nvar;
	double		*wvars;
	psnfunct	*functs;
	lfitdata	*lf;
 } fitfunctdata;

typedef struct
 {	char		*file;
	void		*handle;
 } lfitdynlib;

/*****************************************************************************/

char *key_seek_real_position(char *key,int chrskip)
{
 if ( key==NULL )
 	return(NULL);
 else if ( chrskip>=0 )
	return(key+chrskip);
 else if ( *key != '-' )
	return(key);
 else if ( *(key+1) != '-' )
	return(key+2);
 else
  {	key+=2;
	while ( *key && *key != '-' )	key++;
	if ( ! (*key) )
		return(NULL);
	key++;
	return(key);
  }
}

int key_add_list(char ***rkeys,int *rnkey,char **klist,int chrskip)
{
 char	**keys,**p,*key;
 int	i,nkey;

 keys=*rkeys;
 nkey=*rnkey;

 for ( p=klist ; klist != NULL && (*p) ; p++ )
  {	key=(*p)+chrskip;
	if ( ! (*key) )	
		continue;
	for ( i=0 ; i<nkey ; i++ )
	 {	if ( strcmp(keys[i],key)==0 )
			break;
	 }
	if ( i<nkey )
		continue;
	keys=(char **)realloc(keys,sizeof(char *)*(nkey+1));
	keys[nkey]=key;
	nkey++;
  }

 *rkeys=keys;
 *rnkey=nkey;

 return(0);
}

char *key_search_argument(char **keylist,int chrskip,char *key,char **arglist)
{
 char	*wkey;

 if ( keylist==NULL || arglist==NULL || key==NULL )
	return(NULL);
 while ( *keylist != NULL && *arglist != NULL )
  {	wkey=(*keylist)+chrskip;
	if ( strcmp(wkey,key)==0 )
		return(*arglist);
	else
	 {	keylist++;
		arglist++;
	 }
  }
 return(NULL);
}

datablock *datablock_search_by_key(lfitdata *lf,char *key)
{
 int		i;
 datablock	*db;
 for ( i=0,db=lf->datablocks ; i<lf->ndatablock ; i++,db++ )
  {	if ( strcmp(db->key,key)==0 )
		return(db);
  }
 return(NULL);
}

/*****************************************************************************/

/* constraint_initialize_lambda_matrix(): 
   This function initializes the force matrix 'invlmatrix' which is used 
   by constraint_force(). The matrix 'invlmatrix' should be pre-allocated 
   and must have a size of [nc+nvar] x [nc+nvar]			     */
int constraint_initialize_lambda_matrix(int nc,int nvar,
	double **cmatrix,double **invlmatrix)
{
 int	i,j;
 double	**lmatrix;

 lmatrix=invlmatrix;

 for ( i=0 ; i<nvar ; i++ )
  {	for ( j=0 ; j<nvar ; j++ )
	 {	lmatrix[i][j]=0.0;		}
	lmatrix[i][i]=1.0;
  }
 for ( i=0 ; i<nc ; i++ )
  {	for ( j=0 ; j<nvar ; j++ )
	 {	lmatrix[nvar+i][j]=cmatrix[i][j];
		lmatrix[j][nvar+i]=cmatrix[i][j];
	 }
  }
 for ( i=0 ; i<nc ; i++ )
  {	for ( j=0 ; j<nc ; j++ )
	 {	lmatrix[nvar+i][nvar+j]=0.0;	}
  }

 if ( invert_gauss(invlmatrix,nvar+nc) )
	return(1);
 else
	return(0);
}

int constraint_initialize_proj_matrix(int nc,int nvar,
	double **cmatrix,double **proj)
{
 int	i,j,k;
 double	*wp,w;

 if ( proj==NULL || nvar<=0 )
	return(-1);

 for ( i=0 ; i<nvar ; i++ )
  {	for ( j=0 ; j<nvar ; j++ )
	 {	proj[i][j]=0.0;		}
	proj[i][i]=1.0;
  }

 wp=vector_alloc(nvar);
 for ( k=0 ; k<nc && cmatrix != NULL ; k++ )
  {	w=0.0;
	for ( i=0 ; i<nvar ; i++ )
	 {	wp[i]=0.0;
		for ( j=0 ; j<nvar ; j++ )
		 {	wp[i]+=proj[i][j]*cmatrix[k][j];	}
		w+=wp[i]*wp[i];
	 }
	for ( i=0 ; i<nvar ; i++ )
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	proj[i][j]-=wp[i]*wp[j]/w;	}
	 }
  }
 vector_free(wp);

 return(0);
}

/* constraint_force():
   This function modifies the variable list 'vars' to the nearest value
   ('nearest' is meant as nearest in Eucledian distance) which also 
   satisfy the constraints described by 'invlmatrix' and 'cvector'. Note that
   the 'invlmatrix' contains the constraint information implicitly: it
   should be an [nvar+nc]x[nvar+nc] matrix, initialized by the function 
   constraint_initialize_lambda_matrix(). The value of 'nc' is the total number
   of constraints.							     */
int constraint_force(double *vars,int nvar,int nc,
	double **invlmatrix,double *cvector)
{
 double	*kvector;
 int	i,j;

 kvector=vector_alloc(nvar+nc);

 for ( i=0 ; i<nvar ; i++ )
  {	kvector[i]=vars[i];		}
 for ( i=0 ; i<nc ; i++ )
  {	kvector[nvar+i]=cvector[i];	}

 for ( i=0 ; i<nvar ; i++ )
  {	vars[i]=0.0;
	for ( j=0 ; j<nvar+nc ; j++ )
	 {	vars[i]+=invlmatrix[i][j]*kvector[j];		}
  }

 vector_free(kvector);

 return(0);
}



/*****************************************************************************/

int lfit_psn_chain_optimize(psn *pseq,int *chain)
{
 psnterm	*seq,*w;
 int		c,n;

 for ( seq=pseq->terms,c=0 ; seq->type ; seq++,c++ )
  {	if ( (n=chain[c])<0 )
		continue;
	w=&pseq->terms[n];
	if ( ! ( w->type==T_OP || w->type==T_FN ) )
		chain[c]=-1;
	else if ( ! ( w->major==O_MUL || w->major==O_DIV ) )
		chain[c]=-1;
	else if ( ! ( n>c+1) )
		chain[c]=-1;
  }

 return(0);
}

#define		PSNSTACKBLOCK		64

int lfit_psn_double_calc(psn *pseq,int *chain,psnfunct *flist,
	double *result,double *vars)
{
 double		stack_automatic[PSNSTACKBLOCK],*stack_dynamic,*stack;
 int		i,nopt,narg,sp,c,n;
 psnterm	*seq;
 int		stacklen;

 stack=stack_automatic;
 stack_dynamic=NULL;
 stacklen=PSNSTACKBLOCK;

 nopt=0;sp=0;
 for ( seq=pseq->terms,c=0 ; seq->type ; seq++,c++ )
  {	switch( seq->type )
	 {   case T_CONST:
		stack[sp]=pseq->cons[seq->major],sp++;
	  	break;
	     case T_VAR:
		stack[sp]=vars[seq->major],sp++;
		break;
	     case T_SCONST:
		stack[sp]=(double)(seq->major),sp++;
		break;
	     case T_STACKVAR:
		stack[sp]=stack[seq->major],sp++;
		if ( seq->major+1>nopt )	nopt=seq->major+1;
		break;
	     case T_OP: case T_FN:

		if ( seq->major >= LFIT_F_MAX_BUILTIN )
		 {	lffreg	*lr;
			double	ret;

			narg=seq->minor;
			i=seq->cache;

			if ( i<0 )
			 {	for ( i=0 ; i<lpg.nlffreg ; i++ )
				 {	if ( lpg.lffregs[i].fmajor <= seq->major && seq->major < lpg.lffregs[i].nmajor )
						break;
				 }
				if ( lpg.nlffreg <= i )
				 {	if ( stack_dynamic != NULL )
						free(stack_dynamic);
					return(psnerrno=PENOTFOUND);
				 }
				seq->cache=i;
			 }
			lr=&lpg.lffregs[i];
			i=seq->major-lr->fmajor;
			if ( i==0 )
			 {	double	x;
				if ( lr->lff->function(stack+sp-(lr->lff->nvar+lr->lff->nidv),stack+sp-(lr->lff->nidv),&x,NULL) )
				 {	if ( stack_dynamic != NULL )
						free(stack_dynamic);
					return(psnerrno=PENUMERICAL);
				 }
				ret=x;
				/* fprintf(stderr,"ret=%g\n",ret); */
			 }
			else if ( i>0 )
			 {	if ( lr->diff==NULL )
					lr->diff=(double *)malloc(sizeof(double)*lr->lff->nvar);

				if ( lr->lff->function(stack+sp-(lr->lff->nvar+lr->lff->nidv),stack+sp-(lr->lff->nidv),&ret,lr->diff) )
				 {	if ( stack_dynamic != NULL )
						free(stack_dynamic);
					return(psnerrno=PENUMERICAL);
				 }
				ret=lr->diff[i-1];
			 }
			else
				ret=0.0;

			sp-=narg;
			stack[sp]=ret;
			sp++;

			break;
		 }
		switch ( seq->major )
		 {   case O_ADD:
			sp--;stack[sp-1]+=stack[sp];
			break;
		     case O_SUB:
			sp--;stack[sp-1]-=stack[sp];
			break;
		     case O_MUL:
			sp--;stack[sp-1]*=stack[sp];
			break;
		     case O_DIV:
			if ( stack[sp-1]==0.0 )
			 {	if ( stack_dynamic != NULL )
					free(stack_dynamic);
				return(psnerrno=PENUMERICAL);
			 }
			sp--;stack[sp-1]/=stack[sp];
			break;
		     case O_POW:
			sp--;stack[sp-1]=pow(stack[sp-1],stack[sp]);
			break;
		     case O_CHS:
			stack[sp-1]=-stack[sp-1];
			break;
		     case O_RCP:
			if ( stack[sp-1]==0.0 )
			 {	if ( stack_dynamic != NULL )
					free(stack_dynamic);
				return(psnerrno=PENUMERICAL);
			 }
			stack[sp-1]=1.0/stack[sp-1];
			break;
		     case O_PSQ:
			stack[sp-1]*=stack[sp-1];
			break;
		     case F_SQR:
			if ( stack[sp-1]<0.0 )
			 {	if ( stack_dynamic != NULL )
					free(stack_dynamic);
				return(psnerrno=PENUMERICAL);
			 }
			stack[sp-1]=sqrt(stack[sp-1]);
			break;
		     default:
			narg=seq->minor;
			i=seq->cache;
			if ( i<0 )
			 {	i=psn_cache_search(seq->major,flist);
				if ( i<0 )
				 {	if ( stack_dynamic != NULL )
						free(stack_dynamic);
					return(psnerrno=PENOTFOUND);
				 }
				seq->cache=i;
			 }
			if ( flist[i].funct(stack+sp) )
			 {	if ( stack_dynamic != NULL )
					free(stack_dynamic);
				return(psnerrno=PENUMERICAL);
			 }
			sp+=1-narg;
			break;
		 }
		break;
	 }
	while ( chain != NULL && (n=chain[c])>=0 && stack[sp-1]==0.0 )
	 {	seq=&pseq->terms[n];
		for ( ; c<=n ; c++ )
		 {	if ( pseq->terms[c].type==T_STACKVAR && 
			pseq->terms[c].major+1>nopt )	
				nopt=pseq->terms[c].major+1;
		 }
		c=n;
	 }
	if ( sp>=stacklen )
	 {  if ( stacklen==PSNSTACKBLOCK )
	     {	int	i;
		stacklen+=PSNSTACKBLOCK;
		stack_dynamic=(double *)malloc(sizeof(double)*stacklen);
		for ( i=0 ; i<PSNSTACKBLOCK ; i++ )
		 {	stack_dynamic[i]=stack_automatic[i];	}
	     }
	    else
	     {	stacklen+=PSNSTACKBLOCK;
		stack_dynamic=(double *)realloc(stack_dynamic,
		sizeof(double)*stacklen);
	     }
	    stack=stack_dynamic;
	 }
  }

 for ( i=0 ; i<sp-nopt ; i++ )
  {	result[i]=stack[nopt+i];		}

 if ( stack_dynamic != NULL )
	free(stack_dynamic);

 return(psnerrno=0);
}

/*****************************************************************************/

int fwrite_vars(FILE *fw,double *bvector,variable *vars,int nvar)
{
 int	i;

 for ( i=0 ; i<nvar ; i++ )
  {	fprintf(fw,vars[i].format,bvector[i]);
	fprintf(fw," ");
  }

 return(0);
}
int fwrite_dvars(FILE *fw,double *bvector,dvariable *vars,int nvar)
{
 int	i;

 for ( i=0 ; i<nvar ; i++ )
  {	fprintf(fw,vars[i].format,bvector[i]);
	fprintf(fw," ");
  }

 return(0);
}

int fwrite_vars_chi2(FILE *fw,double *bvector,
	variable *vars,int nvar,double chi2)
{
 int	i;

 for ( i=0 ; i<nvar ; i++ )
  {	fprintf(fw,vars[i].format,bvector[i]);
	fprintf(fw," ");
  }

 if ( chi2>=0.0 )
	fprintf(fw,"%12.5f\n",chi2);
 else
	fprintf(fw,"\n");

 return(0);
}

int fprint_variable_name_list(FILE *fw,variable *vars,int nvar)
{
 int	i,w,l;
 char	*format;
 for ( i=0 ; i<nvar ; i++ )
  {	format=vars[i].format;
	if ( *format == '%' )	format++;
	if ( sscanf(format,"%d",&w)<1 )
		w=0;
	l=strlen(vars[i].name);
	for ( ; l<w ; l++ )
	 	fprintf(fw," ");
	fprintf(fw," %s",vars[i].name);
  }
 return(0);
}

int fprint_derived_variable_name_list(FILE *fw,dvariable *vars,int nvar)
{
 int	i,w,l;
 char	*format;
 for ( i=0 ; i<nvar ; i++ )
  {	format=vars[i].format;
	if ( *format == '%' )	format++;
	if ( sscanf(format,"%d",&w)<1 )
		w=0;
	l=strlen(vars[i].name);
	for ( ; l<w ; l++ )
	 	fprintf(fw," ");
	fprintf(fw," %s",vars[i].name);
  }
 return(0);
}

static int num_of_digits(int n)
{
 int	d;
 for ( d=0 ; n>0 ; n/=10 )
	d++;
 return(d);
}

int fprint_variable_name_index(FILE *fw,variable *vars,int nvar)
{
 int	i,w,l;
 char	*format;
 for ( i=0 ; i<nvar ; i++ )
  {	format=vars[i].format;
	if ( *format == '%' )	format++;
	if ( sscanf(format,"%d",&w)<1 )
		w=0;
	l=num_of_digits(i+1)+2;
	for ( ; l<w ; l++ )
	 	fprintf(fw," ");
	fprintf(fw," [%d]",i+1);
  }
 return(0);
}

int fprint_derived_variable_name_index(FILE *fw,dvariable *vars,int nvar,int offset)
{
 int	i,w,l;
 char	*format;
 for ( i=0 ; i<nvar ; i++ )
  {	format=vars[i].format;
	if ( *format == '%' )	format++;
	if ( sscanf(format,"%d",&w)<1 )
		w=0;
	l=num_of_digits(i+1+offset)+2;
	for ( ; l<w ; l++ )
	 	fprintf(fw," ");
	fprintf(fw," [%d]",i+1+offset);
  }
 return(0);
}

/*****************************************************************************/

/* fit_function(): 
   This is the main callback, the general right-hand side of the least 
   squares or MCMC fitting, to which the independent data should be fitted.  */
void fit_function(void *x,double *a,double *ry,double *rdy,void *param)
{
 int		i,ncol;
 datablock	*db;
 fitfunctdata	*ffd=(fitfunctdata  *)param;
 fitinputrow	*fir=(fitinputrow *)x;

 db=&ffd->lf->datablocks[fir->dbidx];

 if ( db->funct != NULL )
  {	for ( i=0 ; i<ffd->nvar ; i++ )
	 {	ffd->wvars[i]=a[i];			}
	ncol=db->ncol;
	for ( i=0 ; i<ncol ; i++ )
	 {	ffd->wvars[i+ffd->nvar]=fir->x[i];	}

	if ( ry != NULL )
	 {	lfit_psn_double_calc(db->funct,db->fchain,ffd->functs,ry,ffd->wvars);
		/* fprintf(stderr,"y=%g\n",(*ry)); */
	 }

	if ( rdy != NULL )
	 {	for ( i=0 ; i<ffd->nvar ; i++ )
		 {	lfit_psn_double_calc(db->diff[i],NULL,ffd->functs,&rdy[i],ffd->wvars); 
			/* fprintf(stderr,"dy[%d]=%g\n",i,rdy[i]); */
		 }
	 }
  }

 else if ( db->lfunct != NULL )
  {	for ( i=0 ; i<db->lfunct->nvar ; i++ )
	 {	ffd->wvars[i]=a[db->map_var[i]];			}
	for ( i=0 ; i<db->lfunct->nidv ; i++ )
	 {	ffd->wvars[i+ffd->nvar]=fir->x[db->map_idv[i]];		}
	db->lfunct->function(ffd->wvars,ffd->wvars+ffd->nvar,ry,rdy);
	if ( rdy != NULL )
	 {	for ( i=0 ; i<ffd->nvar ; i++ )
		 {	ffd->wvars[i]=0.0;			}
		for ( i=0 ; i<db->lfunct->nvar ; i++ )
		 {	ffd->wvars[db->map_var[i]]=rdy[i];	}
		memcpy(rdy,ffd->wvars,sizeof(double)*ffd->nvar);
	 }
  }
 
}

/* fit_function_optional_derivatives(): 
   Same as above but not all of the partial derivatives are evaluated. */
void fit_function_optional_derivatives(void *x,double *a,double *ry,double *rdy,void *param,
	int nidx,int *idxs)
{
 int		i,j,ncol;
 datablock	*db;
 fitfunctdata	*ffd=(fitfunctdata  *)param;
 fitinputrow	*fir=(fitinputrow *)x;

 if ( nidx<=0 || idxs==NULL || rdy==NULL )	/* imply no request */
  {	rdy=NULL;
	nidx=0;
	idxs=NULL;
  }

 db=&ffd->lf->datablocks[fir->dbidx];

 if ( db->funct != NULL )
  {	for ( i=0 ; i<ffd->nvar ; i++ )
	 {	ffd->wvars[i]=a[i];			}
	ncol=db->ncol;
	for ( i=0 ; i<ncol ; i++ )
	 {	ffd->wvars[i+ffd->nvar]=fir->x[i];	}

	if ( ry != NULL )
	 {	lfit_psn_double_calc(db->funct,db->fchain,ffd->functs,ry,ffd->wvars);	}

	for ( j=0 ; j<nidx ; j++ )
	 {	i=idxs[j];
		lfit_psn_double_calc(db->diff[i],NULL,ffd->functs,&rdy[i],ffd->wvars); 
	 }
  }

 /* well, yes, the lfit external API does not support specifications
    for optional derivatives, so this is the same as in fit_function(). */
 else if ( db->lfunct != NULL )
  {	for ( i=0 ; i<db->lfunct->nvar ; i++ )
	 {	ffd->wvars[i]=a[db->map_var[i]];			}
	for ( i=0 ; i<db->lfunct->nidv ; i++ )
	 {	ffd->wvars[i+ffd->nvar]=fir->x[db->map_idv[i]];		}
	db->lfunct->function(ffd->wvars,ffd->wvars+ffd->nvar,ry,rdy);
	if ( rdy != NULL )
	 {	for ( i=0 ; i<ffd->nvar ; i++ )
		 {	ffd->wvars[i]=0.0;			}
		for ( i=0 ; i<db->lfunct->nvar ; i++ )
		 {	ffd->wvars[db->map_var[i]]=rdy[i];	}
		memcpy(rdy,ffd->wvars,sizeof(double)*ffd->nvar);
	 }
  }
 
}

/*****************************************************************************/

int write_fit_results(fitout *off,int errdump,double *bvector,double *evector,
	int nvar,int nrow,int resdump,int is_dump_delta,
	void **fitpnt,double *fiteval,fitfunctdata *ffd,
	double *fitwght,double *fitdep,fitinputrow *fitinputrows,
	lfitdata *lf,variable *vars)
{
 int	i;

 if ( off->f_write != NULL )
  {	switch ( errdump )
	 {  case 1:
		fwrite_vars(off->f_write,bvector,vars,nvar);
		fprintf(off->f_write,"\n");	
		fwrite_vars(off->f_write,evector,vars,nvar);
		fprintf(off->f_write,"\n");
		break;
	     case 2:
		for ( i=0 ; i<nvar ; i++ )
		 {	fprintf(off->f_write,"%11g %11g ",bvector[i],evector[i]);		}
		fprintf(off->f_write,"\n");	
		break;
	     case 3:
		for ( i=0 ; i<nvar ; i++ )
		 {	fprintf(off->f_write,"%11g %11g\n",bvector[i],evector[i]);	}
		break;
	     default:
		fwrite_vars(off->f_write,bvector,vars,nvar);
		fprintf(off->f_write,"\n");
		break;
	 }
	if ( resdump>=0 )
	 {	double	d,sdd,sig,w;
		sdd=0.0;
		for ( i=0 ; i<nrow ; i++ )
		 {	fit_function(fitpnt[i],bvector,&fiteval[i],NULL,ffd);
			w=fitwght[i];
			d=fiteval[i]-fitdep[i];
			if ( resdump & RESIDUAL_CHI2 )
				sdd+=d*d*w;
			else
				sdd+=d*d;
		 }
		if ( resdump & RESIDUAL_UNBIASED )
		 {	if ( nrow>nvar )
				sdd/=(double)(nrow-nvar);
			else
				sdd=0.0;
		 }
		else
			sdd/=(double)nrow;

		if ( resdump & RESIDUAL_CHI2 )
			sig=sdd;
		else
			sig=sqrt(sdd);

		fprintf(off->f_write,"%11g\n",sig);
	 }
  }

 for ( i=0 ; i<nrow && ( off->f_lused != NULL || off->f_lrejd != NULL || off->f_lall != NULL ) ; i++ )
  {	double	delta;
	int	is_last_in_block;

	if ( i<nrow-1 && fitinputrows[i].dbidx != fitinputrows[i].dbidx )
		is_last_in_block=1;
	else
		is_last_in_block=0;

	if ( is_dump_delta )
	 {	fit_function(fitpnt[i],bvector,&fiteval[i],NULL,ffd);
		delta=(fitdep[i]-fiteval[i]);
	 }
	else
		delta=0.0;

	if ( off->f_lused != NULL && fitwght[i]>0.0 )
	 {	fprintf(off->f_lused,"%s",fitinputrows[i].line);
		if ( is_dump_delta<0 )	fprintf(off->f_lused,"  ## %12g",delta);
		if ( is_dump_delta>0 )	fprintf(off->f_lused,"     %12g",delta);
		fprintf(off->f_lused,"\n");
		if ( is_last_in_block )	fprintf(off->f_lused,"\n");
	 }
	if ( off->f_lrejd != NULL && fitwght[i]<=0.0 )
	 {	fprintf(off->f_lrejd,"%s",fitinputrows[i].line);
		if ( is_dump_delta<0 )	fprintf(off->f_lrejd,"  ## %12g",delta);
		if ( is_dump_delta>0 )	fprintf(off->f_lrejd,"     %12g",delta);
		fprintf(off->f_lrejd,"\n");
		if ( is_last_in_block )	fprintf(off->f_lrejd,"\n");
	 }
	if ( off->f_lall != NULL )
	 {	fprintf(off->f_lall,"%s",fitinputrows[i].line);
		if ( is_dump_delta<0 )	fprintf(off->f_lall,"  ## %12g",delta);
		if ( is_dump_delta>0 )	fprintf(off->f_lall,"     %12g",delta);
		fprintf(off->f_lall,"\n");
		if ( is_last_in_block )	fprintf(off->f_lall,"\n");
	 }

  }

 if ( off->f_expr != NULL ) 
  {	char		*functstr;
	int		fmajor,i,r,d;
	psnterm		*p;
	psn		*pf;
	datablock	*db;
	psnsym		*mysyms[8];

	mysyms[0]=lf->varsym,
	mysyms[1]=NULL,
	mysyms[2]=lpg.pl_sym,
	mysyms[3]=NULL;

	for ( d=0 ; d<lf->ndatablock ; d++ )
	 {	db=&lf->datablocks[d];
		mysyms[1]=db->colsym;
	
		pf=db->functorig;
		for ( i=0,fmajor=0 ; i<nvar ; i++ )
		 {	r=psn_cons_append(pf,bvector[i]);
			if ( i==0 )	fmajor=r;
		 }
		for ( i=0,p=pf->terms ; i<pf->nterm ; i++,p++ )
		 {	if ( p->type==T_VAR && p->major<nvar )
			 {	p->type=T_CONST;
				p->major=fmajor+p->major;
				p->minor=0;
				p->cache=0;
			 }
		 }
		functstr=psn_convert_symbolic(pf,lpg.pl_prop,lpg.pl_symeval,mysyms,NULL);
		if ( functstr != NULL )
		 {	fprintf(off->f_expr,"%s\n",functstr);
			free(functstr);
		 }
		else
			fprintf(off->f_expr,"-\n");
	 }
  }

 if ( off->f_vval != NULL ) 
  {	for ( i=0 ; i<nvar ; i++ )
	 {	fprintf(off->f_vval,"%s=%g",vars[i].name,bvector[i]);
		if ( i<nvar-1 )	fprintf(off->f_vval,",");
		else		fprintf(off->f_vval,"\n");
 	 }
  }

 return(0);
}

char **	read_fit_file(FILE *fr)
{
 char	**lines,*cline;
 int	i,ln;

 ln=0;

 lines=NULL;
 while ( ! feof(fr) )
  {	cline=freadline(fr);
	if ( cline==NULL )
		break;

	for ( i=0 ; cline[i] ; )
	 {	if ( cline[i]==13 || cline[i]==10 )
			cline[i]=0;
		else
			i++;
	 }
	lines=(char **)realloc(lines,sizeof(char *)*(ln+1));
	lines[ln]=cline;
	ln++;
  }
 lines=(char **)realloc(lines,sizeof(char *)*(ln+1));
 lines[ln]=NULL;
 ln++;
 
 return(lines);
}

int read_fit_data(FILE *fr,char **lines,int dbidx,column *cols,int ncol,
	int nvar,int *rnrow,int maxncol,
	double **rfitdata,double **rfitdep,double **rfitwght,fitinputrow **rfitinputrows,
	psn *pdep,psn *perr,int errtype)
{
 int		ln,nrow,cnc,i,j,is_cont;
 double		*fitdata,*fitdep,*fitwght;
 fitinputrow	*fitinputrows;
 char		*rbuff,**lvarstr,*cline;
 double		f,weight;
 double		*wvars,*lvars;

 if ( fr == NULL && lines == NULL )
	return(-1);
 else if ( fr != NULL && lines != NULL )
	return(-1);

 wvars=(double *)malloc(sizeof(double)*(nvar+maxncol));
 lvars=(double *)malloc(sizeof(double)*(nvar+maxncol));

 nrow=*rnrow;

 fitdata =*rfitdata;
 fitdep  =*rfitdep;
 fitwght =*rfitwght;
 fitinputrows=*rfitinputrows;

 rbuff=NULL;lvarstr=NULL;

 ln=0;
 while ( (fr != NULL && ! feof(fr)) || (lines != NULL && lines[ln] != NULL) )
  {	ln++;
	if ( rbuff   != NULL )	{ free(rbuff);  rbuff=NULL;   }
	if ( lvarstr != NULL )	{ free(lvarstr);lvarstr=NULL; }

	if ( fr != NULL )	rbuff=freadline(fr);
	else			rbuff=strdup(lines[ln-1]);

	if ( rbuff==NULL )	break;
	cline=strdup(rbuff);

	/* parsing, token extaction, ...: */
	for ( i=0 ; cline[i] ; )
	 {	if ( cline[i]==13 || cline[i]==10 )
			cline[i]=0;
		else
			i++;
	 }

	for ( i=0 ; i<strlen(rbuff) ; i++ )
	 {	if ( rbuff[i]==',' || rbuff[i]==';' )	rbuff[i]=32;
		if ( rbuff[i]=='#' )			rbuff[i]=0;
	 }
	
	lvarstr=tokenize_spaces_dyn(rbuff);
	if ( lvarstr==NULL || lvarstr[0]==NULL )
	 {	free(cline);
		/* if ( lvarstr != NULL )	free(lvarstr); */
		continue;
	 }
	for ( cnc=0 ; lvarstr[cnc] != NULL ; )	cnc++;

	if ( cnc<ncol )
	 {	if ( is_verbose>=0 )
			fprint_warning("missing column in line %d, skipped",ln);
		free(cline);
		continue;
	 }

	/* we check the used columns for valid (numerical) content: */
	is_cont=0;
	for ( i=0 ; i<ncol && (!is_cont) ; i++ )
	 {	if ( cols[i].name != NULL && cols[i].name[0] )
		 {	j=sscanf(lvarstr[i],"%lg",&lvars[i]);
			if ( j<1 || (! isfinite(lvars[i])) )
			 {	if ( is_verbose>=0 )
					fprint_warning("inappropriate field '%s' in line %d, skipped",lvarstr[i],ln);
				is_cont=1;
			 }
		 }
		else	lvars[i]=0.0;
	 }
	if ( is_cont )
	 {	free(cline);
		continue;
	 }

	for ( i=0 ; i<nvar ; i++ )
	 {	wvars[i]=0;		}
	for ( i=0 ; i<ncol ; i++ )
	 {	wvars[i+nvar]=lvars[i];	}

	/* evaluating the weights and errors: */
	if ( perr != NULL )
	 {	lfit_psn_double_calc(perr,NULL,lpg.pl_funct,&f,wvars);
		if ( errtype )		weight=f;		/* ...weight */
		else if ( f>0.0 )	weight=1.0/f;		/* sigma...  */
		else			weight=0.0;		/* unexpected*/
		/*
		if ( ! errtype )	weight=1.0/f;
		else			weight=f;
		*/
	 }
	else	weight=1.0;

	/*if ( weight == 0.0 )	continue;*/

	/* evaluating the current term: */
	lfit_psn_double_calc(pdep,NULL,lpg.pl_funct,&f,wvars);

	fitdata =(double *)realloc(fitdata ,(nrow+1)*sizeof(double)*(maxncol));
	fitdep  =(double *)realloc(fitdep  ,(nrow+1)*sizeof(double));
	fitwght =(double *)realloc(fitwght ,(nrow+1)*sizeof(double));
	fitinputrows=(fitinputrow *)realloc(fitinputrows,(nrow+1)*sizeof(fitinputrow));

	j=(nrow)*(maxncol);
	for ( i=0 ; i<ncol ; i++ )
	 {	fitdata[j+i]=lvars[i];	}
	for ( ; i<maxncol ; i++ )
	 {	fitdata[j+i]=0.0;	}

	fitinputrows[nrow].line =cline;
	fitinputrows[nrow].dbidx=dbidx;
	fitinputrows[nrow].x    =NULL;
	fitdep [nrow]=f;
	fitwght[nrow]=weight;

	nrow++;
  }

 if ( lvarstr != NULL )	free(lvarstr);
 if ( rbuff != NULL )	free(rbuff);

 free(lvars);
 free(wvars);

 *rnrow=nrow;
 *rfitdata =fitdata;
 *rfitdep  =fitdep;
 *rfitwght =fitwght;

 *rfitinputrows=fitinputrows;

 return(0);
}

int read_all_fit_data(lfitdata *lf,int nvar,
	double **rfitdata,double **rfitdep,double **rfitwght,fitinputrow **rfitinputrows,
	int *rnrow)
{
 int		d,nstdin;
 datablock	*db;
 FILE		*fr;
 char		**stdinlines,**lines;

 nstdin=0;
 for ( d=0 ; d<lf->ndatablock ; d++ )
  {	db=&lf->datablocks[d];

	if ( db->inparg==NULL || strcmp(db->inparg,"-")==0 )
		nstdin++;
  }
 if ( nstdin >= 2 )
 	stdinlines=read_fit_file(stdin);
 else
	stdinlines=NULL;

 for ( d=0 ; d<lf->ndatablock ; d++ )
  {	db=&lf->datablocks[d];

	db->offset=*rnrow;

	if ( db->inparg==NULL || strcmp(db->inparg,"-")==0 )
	 {	if ( stdinlines != NULL )
		 {	lines=stdinlines;
			fr=NULL;
		 }
		else
		 {	lines=NULL;
			fr=stdin;
		 }
	 }
	else
	 {	fr=fopen(db->inparg,"rb");
		if ( fr==NULL )
		 {	fprint_warning("unable to read file '%s' to read fit data",db->inparg);
			continue;
		 }
		lines=NULL;
	 }
	read_fit_data(fr,lines,d,db->cols,db->ncol,nvar,rnrow,lf->maxncol,
		rfitdata,rfitdep,rfitwght,rfitinputrows,
		db->dep,db->err,db->errtype);

	if ( fr != NULL && fileno(fr) != fileno(stdin) )
		fclose(fr);

	db->size=(*rnrow)-db->offset;
  }

 if ( stdinlines != NULL )
  {	char	**line;
	for ( line=stdinlines ; *line ; line++ )
	 {	free(*line);		}
	free(stdinlines);
	stdinlines=NULL;
  }

 return(0);
}

/* fit_linear_or_nonlinear():
   This function is called fit_linear_or_nonlinear() because the cruical
   part of the whole thing is not the fit itself (even linear or nonlinear) 
   but the 'auxiliary' stuff like reading/parsing the data, do the rejection, 
   take into account the effects of the contstraints, whatsoever... The fitting
   is done in a standalone library, see lmfit.[ch].			     */
int fit_linear_or_nonlinear(lfitdata *lf,fitout *off,
	variable *vars,int nvar,int errtype,
	int errdump,int resdump,int is_dump_delta,int is_weighted_sigma,int is_numeric_derivs)
{
 int		i,j,k,n,isig,nreject,nc;

 void		**fitpnt;
 double		*fitdata,*fitdep,*fiteval,*fitwght,lam,lam_mpy;
 double		*wvars;
 fitinputrow	*fitinputrows;
 int		nrow;
 double		**amatrix,*bvector,*evector,*xvector;
 fitfunctdata	ffd;
 constraint	ccstat,*cc;
 fitparam	*fp;

 fp=&lf->parameters;

 if ( lf->fconstraint.nc>0 )	/* we have some constraints */
  {	nc= ccstat.nc= lf->fconstraint.nc;
	ccstat.cmatrix=lf->fconstraint.cmatrix;
	ccstat.cvector=lf->fconstraint.cvector;
 	cc=&ccstat;
  }
 else			/* we don't have any constraint */
  {	cc=NULL;
	nc=0;
  }

 fitdata =NULL;
 fitdep  =NULL;
 fitwght =NULL;
 fitinputrows=NULL;
 nrow=0;

 read_all_fit_data(lf,nvar,&fitdata,&fitdep,&fitwght,&fitinputrows,&nrow);

 if ( nrow < nvar-nc )	
  {	fprint_error("too few lines for fitting");
	exit(1);
  }

 for ( i=0 ; i<nrow ; i++ )
  {	fitwght[i] = fitwght[i]*fitwght[i];		}

 amatrix=matrix_alloc(nvar);	/* least squares matrix */
 bvector=vector_alloc(nvar);	/* right-hand side of the last squares equ.  */
 evector=vector_alloc(nvar);	/* the error (inverse of the vairance) matrix*/
 xvector=vector_alloc(nvar);	/* differences for numerical derivatives     */

 for ( i=0 ; i<nvar ; i++ )
  {	bvector[i]=0.0;
	xvector[i]=vars[i].diff;
	for ( j=0 ; j<nvar ; j++ )	amatrix[i][j]=0.0;
  } 

 wvars=(double *)malloc(sizeof(double)*(nvar+lf->maxncol));

 fitpnt=(void **)malloc(nrow*sizeof(void *));
 for ( i=0 ; i<nrow ; i++ )
  {	fitinputrows[i].x=&fitdata[i*lf->maxncol];
	fitpnt[i]=(void *)(&fitinputrows[i]);
  }

 ffd.nvar=nvar;
 ffd.wvars=wvars;
 ffd.functs=lpg.pl_funct;
 ffd.lf=lf;

 fiteval=(double *)malloc(sizeof(double)*nrow);

 for ( isig=0 ; isig<=fp->niter ; isig++ )
  {	
	if ( lf->is_linear )
	 {	if ( errdump<=0 )
			i=lin_fit_con(fitpnt,fitdep,bvector,fitwght,fit_function,nvar,nrow,&ffd,cc,NULL);
		else
			i=lin_fit_con(fitpnt,fitdep,bvector,fitwght,fit_function,nvar,nrow,&ffd,cc,evector);
		if ( i )		
		 {	fprint_error("singular matrix");
			exit(1);
		 }
	 }
	else 
	 {	for ( i=0 ; i<nvar ; i++ )
		 {	bvector[i]=vars[i].init;		}
		lam    =fp->tune.nllm.lambda;
		lam_mpy=fp->tune.nllm.lambda_mpy;
		for ( n=0 ; n<fp->tune.nllm.max_iter ; n++ )
		 {	if ( is_verbose>=2 )
			 {	fwrite_vars(stderr,bvector,vars,nvar);
				fprintf(stderr,"(%11g)\n",lam);
			 }
			if ( ! is_numeric_derivs)
				lam=nlm_fit_base_con(fitpnt,fitdep,bvector,fitwght,fit_function,nvar,nrow,&ffd,cc,lam,lam_mpy);
			else
				lam=nlm_fit_nmdf_con(fitpnt,fitdep,bvector,fitwght,fit_function,nvar,nrow,&ffd,cc,lam,lam_mpy,xvector);
		 }
	 }

	if ( isig<fp->niter )
	 {	double	d,s,sdd,sig,maxdev,w;
		
		s=sdd=0.0;
		for ( i=0 ; i<nrow ; i++ )
		 {	fit_function(fitpnt[i],bvector,&fiteval[i],NULL,&ffd);
			d=fiteval[i]-fitdep[i];
			if ( is_weighted_sigma )
			 {	w=fitwght[i];
				s  +=w;
				sdd+=d*d*w;
			 }
			else
			 {	s+=1.0;
				sdd+=d*d;
			 }
		 }
		sdd/=s;
		sig =sqrt(sdd);

		maxdev=sig*fp->sigma;
		/*fprintf(stderr,"fp->sigma=%g\n",fp->sigma);*/
		for ( i=0,j=0,k=0 ; i<nrow ; i++ )
		 {	if ( fabs(fiteval[i]-fitdep[i]) > maxdev && fitwght[i]>0.0 )
			 {	fitwght[i]=0.0;
				j++;
			 }
			if ( fitwght[i]>0.0 )
				k++;
		 }
		nreject=j;
		if ( is_verbose>0 )
			fprintf(stderr,"Iteration#%3d: rejected/be used/all: %4d/%4d/%4d\n",isig+1,j,k,nrow);
	 }
	else
		nreject=0;

	if ( is_verbose>0 && fp->niter>0 && ( isig==fp->niter || nreject==0 ) )
	 {	for ( i=0,k=0 ; i<nrow ; i++ )
		 {	if ( fitwght[i]>0.0 )	k++;	}
		fprintf(stderr,"In the last  : used/all: %4d/%4d\n",k,nrow);
	 }
	
	if ( ! ( nreject>0 ) )	break;
  }

 if ( ! lf->is_linear )	errdump=0;

 write_fit_results(off,errdump,bvector,evector,nvar,nrow,resdump,is_dump_delta,
	fitpnt,fiteval,&ffd,fitwght,fitdep,fitinputrows,lf,vars);

 if ( fiteval != NULL )		free(fiteval);
 if ( fitpnt  != NULL )		free(fitpnt);

 vector_free(xvector);
 vector_free(evector);
 vector_free(bvector);
 matrix_free(amatrix);

 if ( fitinputrows != NULL )
  {	for ( i=0 ; i<nrow ; i++ )
	 {	free(fitinputrows[i].line);	}
	free(fitinputrows);
  }
 if ( fitwght  != NULL )	free(fitwght);
 if ( fitdep   != NULL )	free(fitdep);
 if ( fitdata  != NULL )	free(fitdata);

 return(0);
}

/*****************************************************************************/

#define ROTATE(a,i,j,k,l) 		\
 	do				\
	 {	g=a[i][j];		\
		h=a[k][l];		\
		a[i][j]=g-s*(h+g*tau);	\
		a[k][l]=h+s*(g-h*tau);	\
	 } while(0)			

/* jacobi_eigenvectors():
   This function determines the eigenvectors and eigenvalues of the 'n'x'n'
   symmetric matrix 'a'. The eigenvalues are stored in 'd' while the
   normalized eigenvectors are stored in the columns(!) of 'v' (thus
   the i'th eigenvector is v[0][i]...v[n-1][i]). The value 'nrot' is used 
   to return the total number of Jacobi rotations during the diagonalization.
   Note that the 'a' matrix is (partially, but) destroyed during this call.  */
int jacobi_eigenvectors(double **a, int n, double *d, double **v, int *nrot)
{
 int	j,iq,ip,i;
 double	tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

 b=(double *)malloc(sizeof(double)*n);
 z=(double *)malloc(sizeof(double)*n);

 for ( ip=0 ; ip<n ; ip++ )
  {	for ( iq=0 ; iq<n ; iq++ ) 
		v[ip][iq]=0.0;
	v[ip][ip]=1.0;
  }
 for ( ip=0 ; ip<n ; ip++ )  
  {	b[ip]=d[ip]=a[ip][ip];
	z[ip]=0.0;
  }

 *nrot=0;
 for ( i=1 ; i<=50 ; i++ )
  {	sm=0.0;
	for ( ip=0 ; ip < n-1 ; ip++ )
	 {	for ( iq=ip+1 ; iq<n ; iq++ )
			sm += fabs(a[ip][iq]);
	 }
	if ( sm <= 0.0 )
	 {	free(z);
		free(b);
		return(0);
	 }

	if ( i < 4 )
		tresh=0.2*sm/(n*n);
	else
		tresh=0.0;

 	for ( ip=0 ; ip<n-1 ; ip++ )
	 {	for ( iq=ip+1 ; iq<n ; iq++ )
		 {	g=100.0*fabs(a[ip][iq]);
			if ( i > 4 && (fabs(d[ip])+g) == fabs(d[ip]) && (fabs(d[iq])+g) == fabs(d[iq]) )
				a[ip][iq]=0.0;
			else if ( fabs(a[ip][iq]) > tresh )
			 {	h=d[iq]-d[ip];
				if ( fabs(h)+g == fabs(h) )
					t=(a[ip][iq])/h;
				else
				 {	theta=0.5*h/(a[ip][iq]);
					t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
					if ( theta < 0.0 ) t=-t;
				 }
				c=1.0/sqrt(1+t*t);
				s=t*c;
				tau=s/(1.0+c);
				h=t*a[ip][iq];
				z[ip] -= h;
				z[iq] += h;
				d[ip] -= h;
				d[iq] += h;
				a[ip][iq]=0.0;
				for ( j=0 ; j<ip-1 ; j++ )
				 {	ROTATE(a,j,ip,j,iq);	}
				for ( j=ip+1 ; j<=iq-1 ;j++ )
				 {	ROTATE(a,ip,j,j,iq);	}
				for ( j=iq+1 ; j<n ; j++ )
				 {	ROTATE(a,ip,j,iq,j);	}
				for ( j=0 ; j<n ; j++ )
				 {	ROTATE(v,j,ip,j,iq);	}
				++(*nrot);
			}
		}
	}
	for ( ip=0 ; ip<n ; ip++ )
	 {	b[ip] += z[ip];
		d[ip]=b[ip];
		z[ip]=0.0;
	}
  }
 return(1);
}

#undef ROTATE

/* matrix_root():
   This function returns the square root of the 'n'x'n' symmetrix matrix
   'matrix'. If 'rematrix' is not NULL, the function saves the set of 
   eigenvectors of the root matrix multiplied by the appropriate eigenvalues.*/
double	** matrix_root(double **matrix,int n,double ***rematrix)
{
 double	**root,**orig,**eigenvectors,*eigenvalues;
 int	i,j,k,nrot;

 root=matrix_alloc(n);
 eigenvectors=matrix_alloc(n);
 orig=matrix_alloc(n);
 eigenvalues =vector_alloc(n);

 for ( i=0 ; i<n ; i++ )
  {	for ( j=0 ; j<n ; j++ )
	 {	orig[i][j]=matrix[i][j];
		root[i][j]=0.0;
	 }
  }

 jacobi_eigenvectors(orig,n,eigenvalues,eigenvectors,&nrot);
/*
 for ( i=0 ; i<n ; i++ )
  {	fprintf(stderr,"lambda[%d]=%g\n",i,eigenvalues[i]);	}
 for ( i=0 ; i<n ; i++ )
  {	for ( j=0 ; j<n ; j++ )
		fprintf(stderr,"%12g ",eigenvectors[j][i]);
	fprintf(stderr,"\n");
  }
*/

 for ( i=0 ; i<n ; i++ )
  {	if ( eigenvalues[i]<0.0 )	eigenvalues[i]=0.0;
	else				eigenvalues[i]=sqrt(eigenvalues[i]);
  }
 
 for ( i=0 ; i<n ; i++ )
  {	for ( j=0 ; j<n ; j++ )
	 {	for ( k=0 ; k<n ; k++ )
		 {	root[i][j] += eigenvectors[i][k]*eigenvectors[j][k]*eigenvalues[k];	}
	 }
  }

 if ( rematrix != NULL )
  {	for ( i=0 ; i<n ; i++ )
	 {	for ( k=0 ; k<n ; k++ )
		 {	eigenvectors[i][k] *= eigenvalues[k];		}
	 }
  }
 
 vector_free(eigenvalues);
 matrix_free(orig);

 if ( rematrix==NULL )
	matrix_free(eigenvectors);
 else
	*rematrix = eigenvectors;

 return(root);
}

/* matrix_inverse_root():
   This function returns the inverse of the square root of the 'n'x'n' 
   symmetrix matrix 'matrix'. */
double	** matrix_inverse_root(double **matrix,int n,int mev)
{
 double	**iroot,**orig,**eigenvectors,*eigenvalues,**pmp;
 int	i,j,k,nrot;

 iroot=matrix_alloc(n);
 eigenvectors=matrix_alloc(n);
 orig=matrix_alloc(n);
 eigenvalues =vector_alloc(n);

 for ( i=0 ; i<n ; i++ )
  {	for ( j=0 ; j<n ; j++ )
	 {	orig[i][j]=matrix[i][j];
		iroot[i][j]=0.0;
	 }
  }

 jacobi_eigenvectors(orig,n,eigenvalues,eigenvectors,&nrot);

 pmp=matrix_alloc(n);
 for ( i=0 ; i<n ; i++ )
  {	for ( j=0 ; j<n ; j++ )
	 {	pmp[i][j]=eigenvectors[j][i];		}
  }

 for ( i=0 ; i<n-1 ; )
  {	if ( eigenvalues[i]<eigenvalues[i+1] )
	 {	double	w,*a;
		w=eigenvalues[i],
		eigenvalues[i]=eigenvalues[i+1],
		eigenvalues[i+1]=w;
		a=pmp[i],pmp[i]=pmp[i+1],pmp[i+1]=a;
		if ( i>0 )
			i--;
		else
			i++;
	 }
	else
		i++;
  }

 for ( i=0 ; i<n ; i++ )
  {	if ( eigenvalues[i]<=0.0 )	eigenvalues[i]=0.0;
	else				eigenvalues[i]=1.0/sqrt(eigenvalues[i]);
  }

 if ( mev<0 )	mev=n; 
 for ( i=0 ; i<n ; i++ )
  {	for ( j=0 ; j<n ; j++ )
	 {	for ( k=0 ; k<mev ; k++ )
		 {	iroot[i][j] += pmp[k][i]*pmp[k][j]*eigenvalues[k];	}
	 }
  }

 matrix_free(pmp);
 vector_free(eigenvalues);
 matrix_free(orig);
 matrix_free(eigenvectors);

 return(iroot);
}

/* matrix_inverse_normal():
   This function returns the inverse of the 'n'x'n' symmetrix matrix 'matrix'. */
double	** matrix_inverse_normal(double **matrix,int n,int mev)
{
 double	**iroot,**orig,**eigenvectors,*eigenvalues,**pmp;
 int	i,j,k,nrot;

 iroot=matrix_alloc(n);
 eigenvectors=matrix_alloc(n);
 orig=matrix_alloc(n);
 eigenvalues =vector_alloc(n);

 for ( i=0 ; i<n ; i++ )
  {	for ( j=0 ; j<n ; j++ )
	 {	orig[i][j]=matrix[i][j];
		iroot[i][j]=0.0;
	 }
  }

 jacobi_eigenvectors(orig,n,eigenvalues,eigenvectors,&nrot);

 pmp=matrix_alloc(n);
 for ( i=0 ; i<n ; i++ )
  {	for ( j=0 ; j<n ; j++ )
	 {	pmp[i][j]=eigenvectors[j][i];		}
  }

 for ( i=0 ; i<n-1 ; )
  {	if ( eigenvalues[i]<eigenvalues[i+1] )
	 {	double	w,*a;
		w=eigenvalues[i],
		eigenvalues[i]=eigenvalues[i+1],
		eigenvalues[i+1]=w;
		a=pmp[i],pmp[i]=pmp[i+1],pmp[i+1]=a;
		if ( i>0 )
			i--;
		else
			i++;
	 }
	else
		i++;
  }

 for ( i=0 ; i<n ; i++ )
  {	if ( eigenvalues[i]<=0.0 )	eigenvalues[i]=0.0;
	else				eigenvalues[i]=1.0/eigenvalues[i];
  }

 if ( mev<0 )	mev=n; 
 for ( i=0 ; i<n ; i++ )
  {	for ( j=0 ; j<n ; j++ )
	 {	for ( k=0 ; k<mev ; k++ )
		 {	iroot[i][j] += pmp[k][i]*pmp[k][j]*eigenvalues[k];	}
	 }
  }

 matrix_free(pmp);
 vector_free(eigenvalues);
 matrix_free(orig);
 matrix_free(eigenvectors);

 return(iroot);
}

/*****************************************************************************/

/* get_gaussian():
   This function returns a random value from a Gaussian (normal) distribution
   with the mean value of 'mean' and the standard deviation of 'stddev'.     */
double get_gaussian(double mean,double stddev)
{
 double u,v,r;
 u=sqrt(-2.0*log(drand48()));
 v=2*M_PI*drand48();
 r=u*cos(v);
 return(r*stddev+mean);
}

typedef struct
 {	int	n;
	double	sct;
	double	chi;
 } chistat;



double data_get_chisquare(void **fitpnt,double *fitdep,double *fitwght,double *fiteval,
	int nrow,double *bvector,fitfunctdata *ffd,int is_reduced)
{
 double		chi2,dy,f;
 int		i,n,ndatablock;
 chistat	*lchi;
 fitinputrow	*fir;

 ndatablock=ffd->lf->ndatablock;
 lchi=(chistat *)malloc(sizeof(chistat)*ndatablock);
 for ( i=0 ; i<ndatablock ; i++ )
  {	lchi[i].chi=0.0;
	lchi[i].n  =0;
  }

 if ( fiteval != NULL )
  {	for ( i=0 ; i<nrow ; i++ )
	 {	fir=fitpnt[i];
		fit_function(fitpnt[i],bvector,&fiteval[i],NULL,ffd);
		dy=(fiteval[i]-fitdep[i])*fitwght[i];
		lchi[fir->dbidx].chi+=dy*dy;
		lchi[fir->dbidx].n++;
	 } 
  }
 else
  {	for ( i=0 ; i<nrow ; i++ )
	 {	fir=fitpnt[i];
		fit_function(fitpnt[i],bvector,&f,NULL,ffd);
		dy=(f-fitdep[i])*fitwght[i];
		lchi[fir->dbidx].chi+=dy*dy;
		lchi[fir->dbidx].n++;
	 } 
  }

 chi2=0.0;

 n=0;
 for ( i=0 ; i<ndatablock ; i++ )
  {	chi2+=lchi[i].chi;
	n+=lchi[i].n;
  }

 if ( is_reduced )
	chi2/=(double)n;

 free(lchi);

 return(chi2);
}

double data_get_chisquare_separated_linear(void **fitpnt,double *fitdep,double *fitwght,double *fiteval,
	int nrow,double *ivector,fitfunctdata *ffd,int is_reduced,int nvar,int *idx,int nseparated,double *ovector)
{
 double		chi2,dy,y,f,w,z;
 int		i,k,l,n,ndatablock;
 chistat	*lchi;
 fitinputrow	*fir;
 double		**amatrix,*bvector,*fvars,*fdiff,*tvector;

 if ( nseparated <= 0 || idx==NULL || nrow<nseparated )
  {	if ( ovector != NULL )	/* just to be compatible... */
	 {	for ( i=0 ; i<nvar ; i++ )
		 {	ovector[i]=ivector[i];		}
	 }
	chi2=data_get_chisquare(fitpnt,fitdep,fitwght,fiteval,nrow,ivector,ffd,is_reduced);
	return(chi2);
  }

 tvector=vector_alloc(nvar);

 amatrix=matrix_alloc(nseparated);
 bvector=vector_alloc(nseparated);
 fvars  =vector_alloc(nseparated);
 fdiff  =vector_alloc(nvar);

 for ( i=0 ; i<nvar ; i++ )
  {	tvector[i]=ivector[i];		}

 for ( k=0 ; k<nseparated ; k++ )
  {	tvector[idx[k]]=0.0;
	for ( l=0 ; l<nseparated ; l++ )
	 {	amatrix[k][l]=0.0;		}
	bvector[k]=0.0;
  }

 for ( i=0 ; i<nrow ; i++ )
  {	fit_function_optional_derivatives(fitpnt[i],tvector,&f,fdiff,ffd,nseparated,idx);
	y=fitdep[i]-f;
	w=fitwght[i];
	for ( k=0 ; k<nseparated ; k++ )
	 {	fvars[k]=fdiff[idx[k]];		}
	for ( k=0 ; k<nseparated ; k++ )
	 {	z=w*fvars[k];
		for ( l=0 ; l<=k ; l++ )
		 {	amatrix[k][l] += z*fvars[l];		}
		bvector[k] += z*y;
	 }
  }
 for ( k=0 ; k<nseparated ; k++ )
  {	for ( l=k+1 ; l<nseparated ; l++ )
	 {	amatrix[k][l]=amatrix[l][k];		}
  }

 /* we silently ignore the failure of the matrix inversion...:*/
 if ( solve_gauss(amatrix,bvector,nseparated) )
  {	for ( i=0 ; i<nvar ; i++ )
	 {	tvector[i]=ivector[i];			}
  }
 else
  {	for ( k=0 ; k<nseparated ; k++ )
	 {	tvector[idx[k]]=bvector[k];		}
  }		

 vector_free(fdiff);
 vector_free(fvars);
 vector_free(bvector);
 matrix_free(amatrix);

 ndatablock=ffd->lf->ndatablock;
 lchi=(chistat *)malloc(sizeof(chistat)*ndatablock);
 for ( i=0 ; i<ndatablock ; i++ )
  {	lchi[i].chi=0.0;
	lchi[i].n  =0;
  }

 if ( fiteval != NULL )
  {	for ( i=0 ; i<nrow ; i++ )
	 {	fir=fitpnt[i];
		fit_function(fitpnt[i],tvector,&fiteval[i],NULL,ffd);
		dy=(fiteval[i]-fitdep[i])*fitwght[i];
		lchi[fir->dbidx].chi+=dy*dy;
		lchi[fir->dbidx].n++;
	 } 
  }
 else
  {	for ( i=0 ; i<nrow ; i++ )
	 {	fir=fitpnt[i];
		fit_function(fitpnt[i],tvector,&f,NULL,ffd);
		dy=(f-fitdep[i])*fitwght[i];
		lchi[fir->dbidx].chi+=dy*dy;
		lchi[fir->dbidx].n++;
	 } 
  }

 chi2=0.0;

 n=0;
 for ( i=0 ; i<ndatablock ; i++ )
  {	chi2+=lchi[i].chi;
	n+=lchi[i].n;
  }

 if ( is_reduced )
	chi2/=(double)n;

 free(lchi);

 /* save best linear fit 'tvector' to 'ovector' only if we need: */
 if ( ovector != NULL )
  {	for ( i=0 ; i<nvar ; i++ )
	 {	ovector[i]=tvector[i];		}
  }
 /* note that 'ovector' can be in fact the same as 'ivector'... */

 vector_free(tvector);

 return(chi2);
}

int data_get_scatters(void **fitpnt,double *fitdep,double *fitwght,double *fiteval,
	int nrow,double *bvector,fitfunctdata *ffd,int is_normalized)
{
 double		dy,de,f;
 int		i,ndatablock;
 chistat	*lchi;
 fitinputrow	*fir;

 ndatablock=ffd->lf->ndatablock;
 lchi=(chistat *)malloc(sizeof(chistat)*ndatablock);

 for ( i=0 ; i<ndatablock ; i++ )
  {	lchi[i].sct=0.0;
	lchi[i].chi=0.0;
	lchi[i].n  =0;
  }

 for ( i=0 ; i<nrow ; i++ )
  {	fir=fitpnt[i];
	if ( fiteval != NULL )
	 {	fit_function(fitpnt[i],bvector,&fiteval[i],NULL,ffd);
		dy=(fiteval[i]-fitdep[i]);
	 }
	else
	 {	fit_function(fitpnt[i],bvector,&f,NULL,ffd);
		dy=(f-fitdep[i]);
	 }
	/* fprintf(stderr,">> %12g\n",dy); */
	de=dy*fitwght[i];
	lchi[fir->dbidx].sct+=dy*dy;
	lchi[fir->dbidx].chi+=de*de;
	lchi[fir->dbidx].n++;
  } 

 for ( i=0 ; i<ndatablock ; i++ )
  {	if ( lchi[i].n > 0 )
	 {	ffd->lf->datablocks[i].dbscatter=sqrt(lchi[i].sct/lchi[i].n);
		ffd->lf->datablocks[i].dbredchi2=sqrt(lchi[i].chi/lchi[i].n);
	 }
	else
	 {	ffd->lf->datablocks[i].dbscatter=0.0;
		ffd->lf->datablocks[i].dbredchi2=0.0;
	 }
  }	

 free(lchi);

 return(0);
}

double **data_get_inversecovariance(void **fitpnt,double *fitdep,double *fitwght,
	int nrow,int nvar,double *bvector,fitfunctdata *ffd)
{
 double		y,*dy,**cmatrix,isig2;
 int		i,k,l,ndatablock;
 fitinputrow	*fir;

 cmatrix=matrix_alloc(nvar);
 dy=vector_alloc(nvar);

 for ( k=0 ; k<nvar ; k++ )
  {	for ( l=0 ; l<nvar ; l++ )
	 {	cmatrix[k][l]=0.0;		}
  }

 ndatablock=ffd->lf->ndatablock;

 for ( i=0 ; i<nrow ; i++ )
  {	fir=fitpnt[i];
	fit_function(fitpnt[i],bvector,&y,dy,ffd);
	isig2=fitwght[i]*fitwght[i];
	for ( k=0 ; k<nvar ; k++ )
	 {	for ( l=k ; l<nvar ; l++ )
		 {	cmatrix[k][l] += dy[k]*dy[l]*isig2;
		 }
	 }
  }

 for ( k=0 ; k<nvar ; k++ )
  {	for ( l=0 ; l<k ; l++ )
	 {	cmatrix[k][l]=cmatrix[l][k];		}
  }


 vector_free(dy);

 return(cmatrix);
}

double **data_get_rootcovariance(void **fitpnt,double *fitdep,double *fitwght,
	int nrow,int nvar,double *bvector,fitfunctdata *ffd,int nc,double **cproject)
{
 double		**icmatrix,**rootcv;

 icmatrix=data_get_inversecovariance(fitpnt,fitdep,fitwght,nrow,nvar,bvector,ffd);

 if ( cproject != NULL && 0<nc )
  {	double	**t1,**t2;
	int	i,j,k;
	t1=matrix_alloc(nvar);
	t2=matrix_alloc(nvar);
	for ( i=0 ; i<nvar ; i++ ) /* t1 = icmatrix * cproject */
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	t1[i][j]=0.0;
			for ( k=0 ; k<nvar ; k++ )
			 {	t1[i][j]+=icmatrix[i][k]*cproject[k][j];	}
		 }
	 }
	for ( i=0 ; i<nvar ; i++ ) /* t2 = cproject * (icmatrix * cproject) */
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	t2[i][j]=0.0;
			for ( k=0 ; k<nvar ; k++ )
			 {	t2[i][j]+=cproject[i][k]*t1[k][j];	}
		 }
	 }
	rootcv=matrix_inverse_root(t2,nvar,nvar-nc);
	matrix_free(t2);
	matrix_free(t1);
  }
 else
 	rootcv=matrix_inverse_root(icmatrix,nvar,-1);

 matrix_free(icmatrix);

 return(rootcv);
}

double **data_get_covariance(void **fitpnt,double *fitdep,double *fitwght,
	int nrow,int nvar,double *bvector,fitfunctdata *ffd)
{
 double		**icmatrix;

 icmatrix=data_get_inversecovariance(fitpnt,fitdep,fitwght,nrow,nvar,bvector,ffd);

 if ( ffd->lf->fconstraint.nc>0 )
  {	double	**t1,**t2;
	int	i,j,k;
	t1=matrix_alloc(nvar);
	t2=matrix_alloc(nvar);
	for ( i=0 ; i<nvar ; i++ ) /* t1 = icmatrix * cproject */
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	t1[i][j]=0.0;
			for ( k=0 ; k<nvar ; k++ )
			 {	t1[i][j]+=icmatrix[i][k]*ffd->lf->fconstraint.cproject[k][j];	}
		 }
	 }
	for ( i=0 ; i<nvar ; i++ ) /* t2 = cproject * (icmatrix * cproject) */
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	t2[i][j]=0.0;
			for ( k=0 ; k<nvar ; k++ )
			 {	t2[i][j]+=ffd->lf->fconstraint.cproject[i][k]*t1[k][j];	}
		 }
	 }
	matrix_free(icmatrix);
	icmatrix=matrix_inverse_normal(t2,nvar,nvar-ffd->lf->fconstraint.nc);
	matrix_free(t2);
	matrix_free(t1);
	return(icmatrix);
  }
 else
  {	if ( invert_gauss(icmatrix,nvar) )
	 {	matrix_free(icmatrix);	/* singular: fit is degenerated! */
		return(NULL);
	 }
	else
		return(icmatrix);
  }

}

complex * data_get_noise_dft(void **fitpnt,double *fitdep,double *fitwght,double *fiteval,
	int nrow,double *bvector,fitfunctdata *ffd)
{
 double		dy;
 int		i,d,ndatablock;
 datablock	*db;
 fitinputrow	*fir;
 complex	*res;

 ndatablock=ffd->lf->ndatablock;

 res=(complex *)malloc(sizeof(complex)*nrow);

 for ( i=0 ; i<nrow ; i++ )
  {	fir=fitpnt[i];
	fit_function(fitpnt[i],bvector,&fiteval[i],NULL,ffd);
	dy=(fiteval[i]-fitdep[i]);
	res[i].re=dy;
	res[i].im=0;
  } 

 for ( d=0 ; d<ndatablock ; d++ )
  {	db=&ffd->lf->datablocks[d];
	pbfft_conv(res+db->offset,db->size,0);
  }

 return(res);
}

int data_perturb_wtn(void **fitpnt,double *fiteval,double *fitvary,int nrow,fitfunctdata *ffd)
{
 int		i;
 fitinputrow	*fir;
 double		noise;
 datablock	*db;

 for ( i=0 ; i<nrow ; i++ )
  {	fir=fitpnt[i];
	db=&ffd->lf->datablocks[fir->dbidx];
	noise=db->dbscatter;
	if ( db->xnoise>0.0 )
		noise=sqrt(noise*noise+db->xnoise*db->xnoise);
	fitvary[i]=get_gaussian(fiteval[i],noise);
  }

 return(0);
}

int data_perturb_dft(void **fitpnt,double *fiteval,double *fitvary,int nrow,fitfunctdata *ffd,
	complex *resdft)
{
 int		d,i;
 fitinputrow	*fir;
 double		noise,pc,ps,phase,mul;
 datablock	*db;
 complex	*prtift;

 prtift=(complex *)malloc(sizeof(complex)*nrow);
 for ( i=0 ; i<nrow ; i++ )
  {	phase=2.0*M_PI*drand48();
	pc=cos(phase);
	ps=sin(phase);
	prtift[i].re=pc*resdft[i].re-ps*resdft[i].im;
	prtift[i].im=ps*resdft[i].re+pc*resdft[i].im;
  }

 for ( d=0 ; d<ffd->lf->ndatablock ; d++ )
  {	db=&ffd->lf->datablocks[d];
	pbfft_conv(prtift+db->offset,db->size,1);
  }

 mul=sqrt(2.0);

 for ( i=0 ; i<nrow ; i++ )
  {	fir=fitpnt[i];
	db=&ffd->lf->datablocks[fir->dbidx];
	if ( db->xnoise>0.0 )	noise=get_gaussian(0,db->xnoise);
	else			noise=0.0;
	fitvary[i]=fiteval[i]+mul*prtift[i].re/db->size+noise;
  }

 free(prtift);

 return(0);
}


int find_chi2_errors(void **fitpnt,double *fitdep,double *fitwght,
	int nrow,double *bvector,fitfunctdata *ffd,
	variable *vars,int nvar,
	double *errleft,double *erright,double cchi)
{
 int	i,sign,k;
 double	chi2;
 double	*dvector,delta,dmin,dmax;

 dvector=vector_alloc(nvar);
 for ( i=0 ; i<nvar ; i++ )
	dvector[i]=bvector[i];

 for ( i=0 ; i<nvar ; i++ )
  {	for ( sign=-1 ; sign<=1 ; sign+=2 )
	 {	delta=vars[i].ierr;
		if ( delta > 0.0 )
		 {	dmin=delta;
			for ( k=0 ; k<50 ; k++ )
			 {	dvector[i]=bvector[i]+sign*dmin;
				chi2=data_get_chisquare(fitpnt,fitdep,fitwght,NULL,nrow,dvector,ffd,1);
				if ( chi2<cchi )	break;
				dmin=dmin/2;
			 }
			if ( k==50 )	delta=0.0;
			dmax=dmin*2;
			for ( k=0 ; k<50 && delta>0.0 ; k++ )
			 {	dvector[i]=bvector[i]+sign*dmax;
				chi2=data_get_chisquare(fitpnt,fitdep,fitwght,NULL,nrow,dvector,ffd,1);
				if ( chi2>cchi )	break;
				dmax=dmax*2;
			 }
			if ( k==50 )	delta=0.0;
			for ( k=0 ; k<20 && delta>0.0 ; k++ )
			 {	delta=(dmin+dmax)/2.0;
				dvector[i]=bvector[i]+sign*delta;
				chi2=data_get_chisquare(fitpnt,fitdep,fitwght,NULL,nrow,dvector,ffd,1);
				if ( chi2>cchi )	dmax=delta;
				else			dmin=delta;
			 }
		 }
		if ( sign<0 )	errleft[i]=delta;
		else		erright[i]=delta;
	 } 
	dvector[i]=bvector[i];
  }

 vector_free(dvector);

 return(0);
}


/* fit_markov_chain_monte_carlo():
   This function is the core of the Markov-chain Monte Carlo (MCMC) method.  */
int fit_markov_chain_monte_carlo(lfitdata *lf,fitout *off,variable *vars,int nvar)
{
 int		i,n,nc,nrc,ccnt,is_not_accepted;

 void		**fitpnt;
 double		*fitdata,*fitdep,*fiteval,*fitwght;
 double		*wvars,chi2,nchi,mchi,vx,u,alpha,dalpha;
 fitinputrow	*fitinputrows;
 int		nrow,cgibbs;
 double		*bvector,*dvector,*evector,*mvector;
 double		*vclist; 
 fitfunctdata	ffd;
 FILE		*fw;
 fitparam	*fp;
 
 fp=&lf->parameters;

 nc=lf->fconstraint.nc;	/* number of constraints */
 nrc=nc;	/* number of "real" constraints, i.e. not explicit fixations */
 for ( i=0 ; i<nvar ; i++ )	
  {	if ( vars[i].flags & VAR_IS_CONSTANT )
		nrc--;
  }

 fitdata =NULL;
 fitdep  =NULL;
 fitwght =NULL;
 fitinputrows=NULL;
 nrow=0;

 bvector=vector_alloc(nvar);	/* set of unknown variables	*/
 dvector=vector_alloc(nvar);	/* perturbed variable vector	*/
 evector=vector_alloc(nvar);	/* error vector			*/
 mvector=vector_alloc(nvar);	/* vector of minimal values	*/

 read_all_fit_data(lf,nvar,&fitdata,&fitdep,&fitwght,&fitinputrows,&nrow);

 if ( nrow < nvar-nc )	
  {	fprint_error("too few lines for fitting");
	exit(1);
  }

 wvars=(double *)malloc(sizeof(double)*(nvar+lf->maxncol));

 fitpnt=(void **)malloc(nrow*sizeof(void *));
 for ( i=0 ; i<nrow ; i++ )
  {	fitinputrows[i].x=&fitdata[i*lf->maxncol];
	fitpnt[i]=(void *)(&fitinputrows[i]);
  }

 ffd.nvar=nvar;
 ffd.wvars=wvars;
 ffd.functs=lpg.pl_funct;
 ffd.lf=lf;

 fiteval=(double *)malloc(sizeof(double)*nrow);

 for ( i=0 ; i<nvar ; i++ )
  {	bvector[i]=vars[i].init;			}

 ccnt=1;

/*
   if we have constraints, the variables are pushed to the nearest point
   of the intersections of the hyperplanes defined by the constraints
   (i.e. the variables are forced to satisfy the constraints with the 
   minimal shift distance).
 */
 if ( nc > 0 )
  {	for ( i=0 ; i<nvar ; i++ )
	 {	if ( vars[i].flags & VAR_IS_CONSTANT )
	 	 {	bvector[i]=vars[i].init;		}
	 }
	if ( nrc>0 )
	 {	constraint_force(bvector,nvar,nc,
		lf->fconstraint.invlmatrix,lf->fconstraint.cvector);
	 }
  }

 chi2=data_get_chisquare(fitpnt,fitdep,fitwght,fiteval,nrow,bvector,&ffd,0);
 mchi=chi2;
 for ( i=0 ; i<nvar ; i++ )
  {	mvector[i]=bvector[i];		}
 /* fprintf(stderr,"chi2=%g fitwght[0]=%g\n",chi2,fitwght[0]); */

 vclist=NULL;

 fw=off->f_write;

 if ( fw != NULL )
	fwrite_vars_chi2(fw,bvector,vars,nvar,chi2);

 /* reset the Gibbs auxiliary counter: */
 cgibbs=0;

 /* we have fp->mc_iterations iterations in a single chain: */
 for ( n=0 ; (lf->parameters.tune.mcmc.count_accepted?ccnt:n)<fp->mc_iterations ; n++ )
  {	
	/* we assume accaptance, i.e. not out-of-boundary data */
	is_not_accepted=0;
	if ( ! fp->tune.mcmc.use_gibbs )
	 {	for ( i=0 ; i<nvar ; i++ )
		 {	if ( ! (vars[i].flags & VAR_IS_CONSTANT) )
			 {	/* vary all values by a Gaussian sampler: */
				dvector[i]=get_gaussian(bvector[i],vars[i].ierr);
				if ( (vars[i].flags & VAR_MIN) && dvector[i] < vars[i].imin )
				 {	is_not_accepted=1; /* o-o-b */
					dvector[i]=vars[i].imin;
				 }
				if ( (vars[i].flags & VAR_MAX) && dvector[i] > vars[i].imax )
				 {	is_not_accepted=1; /* o-o-b */
					dvector[i]=vars[i].imax;
				 }
			 }
			else
				dvector[i]=bvector[i];
		 }
	 }
	else 
	 {	/* copy everything... */
		for ( i=0 ; i<nvar ; i++ )
		 {	dvector[i]=bvector[i];		}

		while ( vars[cgibbs].flags & VAR_IS_CONSTANT )
			cgibbs=(cgibbs+1)%nvar;

		/* ... except the cgibbs'th variable which is varied now: */
		i=cgibbs;
		dvector[i]=get_gaussian(bvector[i],vars[i].ierr);
		if ( (vars[i].flags & VAR_MIN) && dvector[i] < vars[i].imin )
		 {	is_not_accepted=1;
			dvector[i]=vars[i].imin;
		 }
		if ( (vars[i].flags & VAR_MAX) && dvector[i] > vars[i].imax )
		 {	is_not_accepted=1;
			dvector[i]=vars[i].imax;
		 }

		/* increase Gibbs counter by one: */
		cgibbs=(cgibbs+1)%nvar;
	 }

	/* the next point is not accepted (i.e. it is an out of bound point) */
	if ( is_not_accepted )
		continue;

	/* again, constraints are forced to the variables: */	
	if ( nc > 0 )
	 {	for ( i=0 ; i<nvar ; i++ )
		 {	if ( vars[i].flags & VAR_IS_CONSTANT )
		 	 {	dvector[i]=vars[i].init;		}
		 }
		if ( nrc>0 )
		 {	constraint_force(dvector,nvar,nc,
			lf->fconstraint.invlmatrix,lf->fconstraint.cvector);
		 }
	 }

	/* we get a new chi square value */
	nchi=data_get_chisquare(fitpnt,fitdep,fitwght,fiteval,nrow,dvector,&ffd,0);

	/* just save it if it is minimal (useful as the fitted values) and   */
	/* to normalize the chisquare function (i.e. multiply the errors by  */
	/* some values if desired...)					     */
	if ( nchi<mchi )
	 {	mchi=nchi;
		for ( i=0 ; i<nvar ; i++ )
		 {	mvector[i]=dvector[i];		}
	 }

	vx=exp(-0.5*(nchi-chi2));
	/* vx=exp(-0.5*(nchi-mchi)); */
	alpha=(vx<1.0?vx:1.0);

	/* it is accepted, good: */
	if ( alpha>0.0 && (u=drand48()) <= alpha )
	 {	for ( i=0 ; i<nvar ; i++ )
		 {	bvector[i]=dvector[i];		}
		chi2=nchi;
		/* increase the acceptance counter: */
		ccnt++;
		if ( fw != NULL )
			fwrite_vars_chi2(fw,bvector,vars,nvar,chi2);
	 }
	/* not accepted, set is_not_accepted to true is a symbolic thing...  */
	else
		is_not_accepted=1;	
  }

 if ( off->f_write != NULL )
  {	fw=off->f_write;

	if ( n>0 )
	 {	alpha=(double)ccnt/(double)n;
		dalpha=sqrt((double)ccnt)/(double)n;
	 }
	else
	 {	alpha=0.0;
		dalpha=0.0;
	 }
	fprintf(fw,"#\n");
	fprintf(fw,"# Accepted transitions / total iterations: %d/%d\n",ccnt,n);
	fprintf(fw,"# Total acceptance ratio: %.5f +/- %.5f\n",alpha,dalpha);

	fprintf(fw,"#\n# chi^2 values:\n");
	mchi=data_get_chisquare(fitpnt,fitdep,fitwght,fiteval,nrow,mvector,&ffd,1);
	fprintf(fw,"#\tminimal: %g\n",mchi);
	fprintf(fw,"# Appropriate values for this chi^2:\n# ");
	fwrite_vars(fw,mvector,vars,nvar);
	fprintf(fw,"\n");

  }

 if ( fiteval != NULL )		free(fiteval);
 if ( fitpnt  != NULL )		free(fitpnt);
 if ( fitinputrows != NULL )
  {	for ( i=0 ; i<nrow ; i++ )
	 {	free(fitinputrows[i].line);	}
	free(fitinputrows);
  }
 if ( fitwght  != NULL )	free(fitwght);
 if ( fitdep   != NULL )	free(fitdep);
 if ( fitdata  != NULL )	free(fitdata);

 vector_free(mvector);
 vector_free(evector);
 vector_free(dvector);
 vector_free(bvector);

 return(0);
}

/*****************************************************************************/

typedef struct
 {	int	n,curr;
	double	imin,istp;
	int	varindex;
 } mapaxis;

typedef struct
 {	double	chi2;	/* chi2 value of the given grid point	*/
	int	next;	/* index of the maximal neighbour point	*/
	int	lmin;	/* local minimum in the grid		*/
 } gridpoint;

int map_chi2_downlink(mapaxis *axes,int ndim,gridpoint *griddata)
{
 int	i,j,k,t,m,igrid,ngrid,next;
 int	*axdiff;
 double	mchi;
 
 ngrid=1;
 for ( i=0 ; i<ndim ; i++ )
  {	axes[i].curr=0;
	ngrid *= axes[i].n;
  }

 axdiff=(int *)malloc(sizeof(int)*ndim);

 igrid=0;
 while ( 1 )
  {
	next=-1;
	for ( i=0 ; i<ndim ; i++ )
		axdiff[i]=-1;		
	mchi=0.0;
	while ( 1 )
	 {	t=0,m=1;
		for ( i=0 ; i<ndim ; i++ )
		 {	k=axes[i].curr+axdiff[i];
			if ( k<0 || k>=axes[i].n )
				break;
			t+=k*m;
			m*=axes[i].n;
		 }
		if ( i >= ndim )
		 {	if ( next<0 || griddata[t].chi2<mchi )
				next=t,mchi=griddata[t].chi2;
		 }
		for ( i=0 ; i<ndim ; i++ )
		 {	axdiff[i]++;
			if ( axdiff[i] <= 1 )
				break;
			else
				axdiff[i]=-1;
		 }
		if ( i>=ndim )
			break;
	 }
	griddata[igrid].next=next;

	for ( i=0 ; i<ndim ; i++ )
	 {	axes[i].curr++;
		if ( axes[i].curr < axes[i].n )
			break;
		else
		 	axes[i].curr=0;
	 }
	if ( i>=ndim )
		break;

	igrid++;
  };

 for ( i=0 ; i<ngrid ; i++ )
  {	for ( j=i ; griddata[j].next != j ; )
		j=griddata[j].next;
	griddata[i].lmin=j;
  }

 free(axdiff);
 
 return(0);
}

int map_chi2_partitions(gridpoint *griddata,int ngrid,int **rmins,int *rnmin)
{
 int	*mins,*gidx;
 int	i,nmin,lmin;

 mins=NULL;
 nmin=0;
 gidx=(int *)malloc(sizeof(int)*ngrid);
 for ( i=0 ; i<ngrid ; i++ )
  {	gidx[i]=-1;		}
 for ( i=0 ; i<ngrid ; i++ )
  {	lmin=griddata[i].lmin;
	if ( gidx[lmin] < 0 )
	 {	mins=(int *)realloc(mins,sizeof(int)*(nmin+1));
		mins[nmin]=lmin;
		nmin++;
		gidx[lmin]=lmin;
	 }
  }
 free(gidx);

 if ( rmins != NULL )	*rmins=mins;
 if ( rnmin != NULL )	*rnmin=nmin;

 return(0);
}

int map_chi2_fit_minimum(mapaxis *axes,int ndim,gridpoint *griddata,int g,
int bdist,double **mmatrix,double *mvector,double *rmscalar)
{
 int	i,j,n,m,t,s,k,ngrid,*axdiff,fdim;
 double	**amatrix,*bvector,**imatrix,*ivector,iscalar,*fvars,*xvars,*xoffs,y;

 fdim=1+ndim+ndim*(ndim+1)/2;
 amatrix=matrix_alloc(fdim);
 bvector=vector_alloc(fdim);
 fvars  =vector_alloc(fdim);
 xvars  =vector_alloc(ndim);
 xoffs  =vector_alloc(ndim);
 axdiff=(int *)malloc(sizeof(int)*ndim);
 imatrix=matrix_alloc(ndim);
 ivector=vector_alloc(ndim);

 ngrid=1;
 for ( i=0 ; i<ndim ; i++ )
  {	n=axes[i].n;
	axdiff[i]=-bdist;
	ngrid *= n;
  }

 for ( i=0 ; i<fdim ; i++ )
  {	for ( j=0 ; j<fdim ; j++ )
	 {	amatrix[i][j]=0.0;		}
	bvector[i]=0.0; 
  }

 s=g;
 for ( i=0 ; i<ndim ; i++ )
  {	k=s%axes[i].n;
	/* xoffs[i]=axes[i].imin+k*axes[i].istp; */
	xoffs[i]=(double)k;
	s/=axes[i].n;
  }

 while ( 1 )
  {	s=g;
	t=0,m=1;
	for ( i=0 ; i<ndim ; i++ )
	 {	n=axes[i].n;
		k=s%n+axdiff[i];
		/* xvars[i]=axes[i].imin+k*axes[i].istp-xoffs[i]; */
		xvars[i]=(double)k-xoffs[i];
		s/=n;
		t+=k*m;
		m*=n;
	 }
	y=griddata[t].chi2;

	fvars[0]=1.0;
	for ( i=0 ; i<ndim ; i++ )
		fvars[i+1]=xvars[i];
	k=ndim+1;
	for ( i=0 ; i<ndim ; i++ )
	 {	fvars[k]=0.5*xvars[i]*xvars[i];
		k++;
		for ( j=i+1 ; j<ndim ; j++ )
		 {	fvars[k]=xvars[i]*xvars[j];
			k++;
		 }
	 }

	for ( i=0 ; i<fdim ; i++ )
	 {	for ( j=0 ; j<fdim ; j++ )
		 {	amatrix[i][j]+=fvars[i]*fvars[j];	}
		bvector[i]+=fvars[i]*y;
	 }

	for ( i=0 ; i<ndim ; i++ )
	 {	axdiff[i]++;
		if ( axdiff[i] <= bdist )
			break;
		else
			axdiff[i]=-bdist;
	 }
	if ( i>=ndim )
		break;
  }

 solve_gauss(amatrix,bvector,fdim);

 iscalar=bvector[0];
 for ( i=0 ; i<ndim ; i++ )
	ivector[i]=bvector[i+1];
 k=ndim+1;
 for ( i=0 ; i<ndim ; i++ )
  {	imatrix[i][i]=bvector[k];
	k++;
	for ( j=i+1 ; j<ndim ; j++ )
	 {	imatrix[i][j]=imatrix[j][i]=bvector[k];
		k++;
	 }
  }

 for ( i=0 ; i<ndim ; i++ )
  {	for ( j=0 ; j<ndim ; j++ )
	 {	mmatrix[i][j]=imatrix[i][j];	}
  }

 solve_gauss(imatrix,ivector,ndim);
 for ( i=0 ; i<ndim ; i++ )
  {	/* mvector[i]=-ivector[i]+xoffs[i]; */
	mvector[i]=axes[i].imin+(xoffs[i]-ivector[i])*axes[i].istp;
  }

 vector_free(ivector);
 matrix_free(imatrix);
 free(axdiff);
 vector_free(xoffs);
 vector_free(xvars);
 vector_free(fvars);
 vector_free(bvector);
 matrix_free(amatrix);
 
 return(0);
}

int fit_map_chi2(lfitdata *lf,fitout *off,variable *vars,int nvar)
{
 int		i,j,n;
 void		**fitpnt;
 double		*fitdata,*fitdep,*fitwght;
 double		*wvars,chi2;
 fitinputrow	*fitinputrows;
 int		nrow,ndim,ngrid,igrid;
 mapaxis	*axes,*cm;
 gridpoint	*griddata;
 
 double		*bvector;
 fitfunctdata	ffd;
 FILE		*fw;
 fitparam	*fp;
 
 fp=&lf->parameters;

 fitdata =NULL;
 fitdep  =NULL;
 fitwght =NULL;
 fitinputrows=NULL;
 nrow=0;

 ndim=0;
 ngrid=0;
 axes=NULL;
 for ( i=0 ; i<nvar ; i++ )
  {	if ( ! (vars[i].flags & VAR_STEP) )
		continue;
	n=1+(int)((vars[i].imax+0.5*vars[i].istp-vars[i].imin)/vars[i].istp);
	if ( n <= 0 )
		continue;
	axes=(mapaxis *)realloc(axes,sizeof(mapaxis)*(ndim+1));
	cm=&axes[ndim];
	cm->n=n;
	cm->curr=0;
	cm->imin=vars[i].imin;
	cm->istp=vars[i].istp;
	cm->varindex=i;
	ndim++;
	if ( ngrid )	ngrid=ngrid*cm->n;
	else		ngrid=cm->n;
  }
 if ( ngrid>0 && ndim>0 )
	griddata=(gridpoint *)malloc(sizeof(gridpoint)*ngrid);
 else
	griddata=NULL;

 bvector=vector_alloc(nvar);	/* set of unknown variables	*/

 read_all_fit_data(lf,nvar,&fitdata,&fitdep,&fitwght,&fitinputrows,&nrow);

 wvars=(double *)malloc(sizeof(double)*(nvar+lf->maxncol));

 fitpnt=(void **)malloc(nrow*sizeof(void *));
 for ( i=0 ; i<nrow ; i++ )
  {	fitinputrows[i].x=&fitdata[i*lf->maxncol];
	fitpnt[i]=(void *)(&fitinputrows[i]);
  }

 ffd.nvar=nvar;
 ffd.wvars=wvars;
 ffd.functs=lpg.pl_funct;
 ffd.lf=lf;

 for ( i=0 ; i<nvar ; i++ )
  {	if ( vars[i].flags & VAR_STEP )
		bvector[i]=vars[i].imin;
	else
		bvector[i]=vars[i].init;
  }

 igrid=0;

 fw=off->f_write;

 while ( 1 )
  {	
	chi2=data_get_chisquare(fitpnt,fitdep,fitwght,NULL,nrow,bvector,&ffd,1);
	if ( fw != NULL )
	 {	fprintf(fw," ");
		fwrite_vars(fw,bvector,vars,nvar);
		fprintf(fw,"%12g\n",chi2);
	 }

	if ( griddata != NULL )
		griddata[igrid].chi2=chi2;	

	for ( i=0 ; i<ndim ; i++ )
	 {	axes[i].curr++;
		j=axes[i].varindex;
		bvector[j]=axes[i].imin+axes[i].istp*axes[i].curr;
		if ( axes[i].curr < axes[i].n )
			break;
		else
		 {	bvector[j]=axes[i].imin;
			axes[i].curr=0;
		 }
	 }
	if ( i >= ndim )
		break;

	igrid++;
  };

 if ( griddata != NULL && fw != NULL )
  {	int	k,s,v;
	int	*mins,nmin;
	int	bdist,b;
	double	**mmatrix,*mvector,mscalar;

	map_chi2_downlink(axes,ndim,griddata);
	mins=NULL;
	map_chi2_partitions(griddata,ngrid,&mins,&nmin);

	mmatrix=matrix_alloc(ndim);
	mvector=vector_alloc(ndim);

	if ( nmin>1 )	fprintf(fw,"# %d minima found:\n",nmin);
	else		fprintf(fw,"# %d minimum found:\n",nmin);

	for ( n=0 ; n<nmin ; n++ )
	 {	fprintf(fw,"#\n");
		k=mins[n];
		bdist=-1;
		for ( i=0 ; i<ndim ; i++ )
		 {	s=k%axes[i].n;
			b=axes[i].n-1-s;
			if ( s<b )	b=s;
			if ( bdist<0 || b<bdist )	bdist=b;
			v=axes[i].varindex;
			bvector[v]=axes[i].imin+s*axes[i].istp;
			k/=axes[i].n;
		 }
		fprintf(fw,"# - local minimum on the grid:\n");
		fprintf(fw,"#");
		fwrite_vars(fw,bvector,vars,nvar);
		fprintf(fw,"%12g\n",griddata[mins[n]].chi2);
		if ( ! ( bdist>0 ) )
		 {	fprintf(fw,"#\t(boundary minimum, highly possible to be false)\n");
			continue;
		 }
		fprintf(fw,"#\t(distinct minimum, distance from boundary: %d)\n",bdist);

		fprintf(fw,"# - interpolated minimum:\n");
		map_chi2_fit_minimum(axes,ndim,griddata,mins[n],3,
			mmatrix,mvector,&mscalar);
		fprintf(fw,"#");
		for ( i=0 ; i<ndim ; i++ )
		 {	v=axes[i].varindex;
			bvector[v]=mvector[i];
		 }
		fwrite_vars(fw,bvector,vars,nvar);
		fprintf(fw,"\n");

	 }

	vector_free(mvector);
	matrix_free(mmatrix);

	if ( mins != NULL )	free(mins);
  }

 if ( fitpnt  != NULL )		free(fitpnt);
 if ( fitinputrows != NULL )
  {	for ( i=0 ; i<nrow ; i++ )
	 {	free(fitinputrows[i].line);	}
 	free(fitinputrows);
  }
 if ( fitwght  != NULL )	free(fitwght);
 if ( fitdep   != NULL )	free(fitdep);
 if ( fitdata  != NULL )	free(fitdata);

 vector_free(bvector);

 if ( griddata != NULL )	free(griddata);
 if ( axes != NULL )		free(axes);

 return(0);
 
}

/*****************************************************************************/

double find_monte_carlo_minimum(void **fitpnt,double *fitdep,double *fitwght,
	int nrow,double *mvector,fitfunctdata *ffd,
	double *evector,int downlink_max_iter,int max_decrease_rejection,double decrease_error_ratio,
	lfitdata *lf,variable *vars,int nvar)
{
 double	*dvector;
 double	mchi,nchi;
 int	i,n,nc,nrc;
 int	decrease_error_counter,is_not_accepted;

 nc=lf->fconstraint.nc;	/* number of constraints */
 nrc=nc;	/* number of "real" constraints, i.e. not explicit fixations */
 for ( i=0 ; i<nvar ; i++ )	
  {	if ( vars[i].flags & VAR_IS_CONSTANT )
		nrc--;
  }

 dvector=vector_alloc(nvar);

 if ( nc > 0 )
  {	for ( i=0 ; i<nvar ; i++ )
	 {	if ( vars[i].flags & VAR_IS_CONSTANT )
	 	 {	mvector[i]=vars[i].init;		}
	 }
	if ( nrc>0 )
	 {	constraint_force(mvector,nvar,nc,
		lf->fconstraint.invlmatrix,lf->fconstraint.cvector);
	 }
  }

 mchi=data_get_chisquare(fitpnt,fitdep,fitwght,NULL,nrow,mvector,ffd,1);
 
 decrease_error_counter=0;

 for ( n=0 ; n<downlink_max_iter ; n++ )
  {	
	if ( decrease_error_counter>=max_decrease_rejection )
	 {	decrease_error_counter=0;
		for ( i=0 ; i<nvar ; i++ )
		 {	evector[i] *= decrease_error_ratio;	}
		/* fprintf(stderr,"decrease at n=%d (%g)\n",n,mchi);  */
	 }

	/* we assume accaptance, i.e. not out-of-boundary data */
	is_not_accepted=0;
	for ( i=0 ; i<nvar ; i++ )
	 {	if ( ! (vars[i].flags & VAR_IS_CONSTANT) )
		 {	/* vary all values by a Gaussian sampler: */
			dvector[i]=get_gaussian(mvector[i],evector[i]);
			if ( (vars[i].flags & VAR_MIN) && dvector[i] < vars[i].imin )
			 {	is_not_accepted=1; /* o-o-b */
				dvector[i]=vars[i].imin;
			 }
			if ( (vars[i].flags & VAR_MAX) && dvector[i] > vars[i].imax )
			 {	is_not_accepted=1; /* o-o-b */
				dvector[i]=vars[i].imax;
			 }
		 }
		else
			dvector[i]=mvector[i];
	 }

	/* the next point is not accepted (i.e. it is an out of bound point) */
	if ( is_not_accepted )
	 {	decrease_error_counter++;
		continue;
 	 }

	/* again, constraints are forced to the variables: */	
	if ( nc > 0 )
	 {	for ( i=0 ; i<nvar ; i++ )
		 {	if ( vars[i].flags & VAR_IS_CONSTANT )
		 	 {	dvector[i]=vars[i].init;		}
		 }
		if ( nrc>0 )
		 {	constraint_force(dvector,nvar,nc,
			lf->fconstraint.invlmatrix,lf->fconstraint.cvector);
		 }
	 }

	/* we get a new chi square value */
	nchi=data_get_chisquare(fitpnt,fitdep,fitwght,NULL,nrow,dvector,ffd,1);

	/* if we do only downlink minimalization, we implicitly not accept   */
	/* all of the cases when the value of \chi^2 increases:		     */
	if ( nchi<mchi )
	 {	mchi=nchi;
		for ( i=0 ; i<nvar ; i++ )
		 {	mvector[i]=dvector[i];		}
		decrease_error_counter=0;
	 }
	else
		decrease_error_counter++;
  }
 
 vector_free(dvector);

 return(mchi);
}

/*****************************************************************************/

int get_separate_index(variable *vars,int nvar,int **ridx)
{
 int	i,nseparate,*idx;

 idx=(int *)malloc(sizeof(int)*nvar);
 for ( i=0,nseparate=0 ; i<nvar ; i++ )
  {	if ( vars[i].is_separated )
	 {	idx[nseparate]=i;
		nseparate++;
	 }
  }
 if ( nseparate <= 0 )
  {	free(idx);
	idx=NULL;
  }

 if ( ridx != NULL )	*ridx=idx;
 
 return(nseparate);
}

/*****************************************************************************/

typedef struct
 {	void 		**fitpnt;
	double		*fitdep;
	double		*fitwght;
	int		nrow;
	fitfunctdata	*ffd;
	int		nvar;
	int		*idx;
	int		nseparated;
 } lfitdhsxdata;

double lfit_downhill_simplex_funct(void *param,double *bvector)
{
 lfitdhsxdata	*ldd=(lfitdhsxdata *)param;
 double		chi2;

 if ( ldd->idx != NULL && 0<ldd->nseparated )
	chi2=data_get_chisquare_separated_linear(ldd->fitpnt,ldd->fitdep,ldd->fitwght,NULL,ldd->nrow,bvector,ldd->ffd,1,ldd->nvar,ldd->idx,ldd->nseparated,NULL);
 else
	chi2=data_get_chisquare(ldd->fitpnt,ldd->fitdep,ldd->fitwght,NULL,ldd->nrow,bvector,ldd->ffd,1);

 return(chi2);
}

int get_error_matrix(double **ematrix,int nvar,double *evector,double **icmatrix,int nc,double **cproject)
{
 double	**pmp,**m1,**m2,w,*a,**rcmatrix;
 int	i,j,k,nrot;
 double	*eigenvalues,**eigenvectors;

 if ( ( evector==NULL && icmatrix==NULL ) || ( evector!=NULL && icmatrix!=NULL ) )
	return(-1);

 if ( nc>0 )
  {
	pmp=matrix_alloc(nvar);
	m1 =matrix_alloc(nvar);
	m2 =matrix_alloc(nvar);

	if ( icmatrix==NULL )
	 {	for ( i=0 ; i<nvar ; i++ )
		 {	for ( j=0 ; j<nvar ; j++ )
			 {	m2[i][j]=0.0;		}
			if ( evector[i] != 0.0 )
				m2[i][i]=1.0/(evector[i]*evector[i]);
			else
				m2[i][i]=1.0;
		 }
	 }
	else
	 {	for ( i=0 ; i<nvar ; i++ )
		 {	for ( j=0 ; j<nvar ; j++ )
			 {	m2[i][j]=icmatrix[i][j];		}
		 }
 	 }

	for ( i=0 ; i<nvar ; i++ )
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	m1[i][j]=0.0;
			for ( k=0 ; k<nvar ; k++ )
			 {	m1[i][j]+=cproject[i][k]*m2[k][j];	}
		 }
	 }
	for ( i=0 ; i<nvar ; i++ )
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	pmp[i][j]=0.0;
			for ( k=0 ; k<nvar ; k++ )
			 {	pmp[i][j]+=m1[i][k]*cproject[k][j];	}
		 }
	 }
	matrix_free(m2);
	matrix_free(m1);

	eigenvectors=matrix_alloc(nvar);
	eigenvalues =vector_alloc(nvar);

	jacobi_eigenvectors(pmp,nvar,eigenvalues,eigenvectors,&nrot);
	for ( i=0 ; i<nvar ; i++ )
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	pmp[i][j]=eigenvectors[j][i];		}
	 }

	for ( i=0 ; i<nvar-1 ; )
	 {	if ( eigenvalues[i]<eigenvalues[i+1] )
		 {	w=eigenvalues[i],
			eigenvalues[i]=eigenvalues[i+1],
			eigenvalues[i+1]=w;
			a=pmp[i],pmp[i]=pmp[i+1],pmp[i+1]=a;
			if ( i>0 )
				i--;
			else
				i++;
		 }
		else
			i++;
	 }

	for ( i=0 ; i<nvar-nc ; i++ )
	 {	if ( eigenvalues[i]>0.0 )
			eigenvalues[i]=1.0/sqrt(eigenvalues[i]);
		else
			eigenvalues[i]=0.0;

		for ( j=0 ; j<nvar ; j++ )
		 {	ematrix[i][j]=eigenvalues[i]*pmp[i][j];	
			/* fprintf(stderr,"%12g ",ematrix[i][j]); */
		 }
		/* fprintf(stderr,"\n"); */
	 }

	vector_free(eigenvalues);
	matrix_free(eigenvectors);

	matrix_free(pmp);
  }

 else
  {	if ( evector != NULL )
	 {	for ( i=0 ; i<nvar ; i++ )
		 {	for ( j=0 ; j<nvar ; j++ )
			 {	ematrix[i][j] = (i==j?evector[i]:0.0);		}
		 }
	 }
	else
	 {	rcmatrix=matrix_inverse_root(icmatrix,nvar,-1);
		for ( i=0 ; i<nvar ; i++ )
		 {	for ( j=0 ; j<nvar ; j++ )
			 {	ematrix[i][j] = rcmatrix[i][j];		}
		 }
		matrix_free(rcmatrix);
	 }
  }

 return(0);
}

double find_downhill_simplex_minimum(void **fitpnt,double *fitdep,double *fitwght,
	int nrow,double *mvector,fitfunctdata *ffd,
	double **ematrix,int nvar,int nsimplex,int *idx,int nseparated)
{
 double		**p,chi2;
 int		i,j,r;
 lfitdhsxdata	ldd;
 fitconstraint	*fc;

 p=matrix_alloc(nvar+1);

 fc=&ffd->lf->fconstraint;

 if ( fc->nc>0 )
  {	constraint_force(mvector,nvar,fc->nc,
		fc->invlmatrix,fc->cvector);
  }
 for ( j=0 ; j<nvar ; j++ )
  {	p[0][j]=mvector[j];	}

 for ( i=0 ; i<nsimplex; i++ )
  {	for ( j=0 ; j<nvar ; j++ )
	 {	p[i+1][j]=mvector[j] + ematrix[i][j];		}
  }

 ldd.fitpnt=fitpnt;
 ldd.fitdep=fitdep;
 ldd.fitwght=fitwght;
 ldd.nrow=nrow;
 ldd.ffd=ffd;

 ldd.nvar=nvar;
 if ( idx != NULL && 0<nseparated )
  {	ldd.idx=idx;
	ldd.nseparated=nseparated;
  }
 else
  {	ldd.idx=NULL;
	ldd.nseparated=0;
  }

 r=downhill_simplex(p,nvar,nsimplex,lfit_downhill_simplex_funct,&ldd,1e-10,1000);

/*
 if ( r<0 )
  {	fprintf(stderr,"downhill_simplex(): failed.\n");	}
*/

 for ( j=0 ; j<nvar ; j++ )
	mvector[j]=p[0][j];

 matrix_free(p);
 
 chi2=data_get_chisquare(fitpnt,fitdep,fitwght,NULL,nrow,mvector,ffd,1);

 return(chi2);
}

/*****************************************************************************/

/* fit_error_monte_carlo_estimation():
   This function is the core of the Error Monte-Carlo Estimation (EMCE).  */
int fit_error_monte_carlo_estimation(lfitdata *lf,fitout *off,variable *vars,int nvar)
{
 int		i,n,nc,nrc;

 void		**fitpnt;
 double		*fitdata,*fitdep,*fitvary,*fiteval,*fitwght;
 double		*wvars,mchi,nchi;
 fitinputrow	*fitinputrows;
 int		nrow,downlink_max_iter,max_decrease_rejection;
 double		decrease_error_ratio;
 double		*bvector,*dvector,*evector,*mvector,**ematrix,*xvector;
 fitfunctdata	ffd;
 FILE		*fw;
 fitparam	*fp;
 int		errdump,resdump,is_dump_delta,use_dft_rednoise=0;
 complex	*resdft;
 constraint	ccstat,*cc;

 errdump=0; 
 resdump=0;
 is_dump_delta=0;
 
 fp=&lf->parameters;

 nc=lf->fconstraint.nc;	/* number of constraints */
 nrc=nc;	/* number of "real" constraints, i.e. not explicit fixations */
 for ( i=0 ; i<nvar ; i++ )	
  {	if ( vars[i].flags & VAR_IS_CONSTANT )
		nrc--;
  }

 if ( lf->fconstraint.nc>0 )	/* we have some constraints */
  {	ccstat.nc=     lf->fconstraint.nc;
	ccstat.cmatrix=lf->fconstraint.cmatrix;
	ccstat.cvector=lf->fconstraint.cvector;
	cc=&ccstat;
  }
 else				/* we don't have any constraint */
	cc=NULL;

 fitdata =NULL;
 fitdep  =NULL;
 fitwght =NULL;
 fitinputrows=NULL;
 nrow=0;

 bvector=vector_alloc(nvar);	/* set of unknown variables		*/
 dvector=vector_alloc(nvar);	/* perturbed variable vector		*/
 evector=vector_alloc(nvar);	/* error vector				*/
 mvector=vector_alloc(nvar);	/* vector of minimal values		*/
 xvector=vector_alloc(nvar);	/* vector for numerical derivatives	*/

 read_all_fit_data(lf,nvar,&fitdata,&fitdep,&fitwght,&fitinputrows,&nrow);

 if ( nrow < nvar-nc )	
  {	fprint_error("too few lines for fitting");
	exit(1);
  }

 wvars=(double *)malloc(sizeof(double)*(nvar+lf->maxncol));

 fitvary=(double *)malloc(nrow*sizeof(double));
 fitpnt=(void **)malloc(nrow*sizeof(void *));
 for ( i=0 ; i<nrow ; i++ )
  {	fitinputrows[i].x=&fitdata[i*lf->maxncol];
	fitpnt[i]=(void *)(&fitinputrows[i]);
  }

 ffd.nvar=nvar;
 ffd.wvars=wvars;
 ffd.functs=lpg.pl_funct;
 ffd.lf=lf;

 fiteval=(double *)malloc(sizeof(double)*nrow);

 for ( i=0 ; i<nvar ; i++ )
  {	mvector[i]=vars[i].init;
	/* fprintf(stderr,"bvector[%d]=%g:%g\n",i,bvector[i],vars[i].ierr); */
  }

 downlink_max_iter=1000;
 max_decrease_rejection=10; 
 decrease_error_ratio=0.5;

 for ( i=0 ; i<nvar ; i++ )
  {	evector[i]=vars[i].ierr;

	xvector[i]=vars[i].diff;
  }

 if ( fp->tune.emce.sub_method==FIT_METHOD_DHSX )
  {	ematrix=matrix_alloc(nvar);
	get_error_matrix(ematrix,nvar,evector,NULL,lf->fconstraint.nc,lf->fconstraint.cproject);
  }
 else
	ematrix=NULL;

 if ( fp->tune.emce.skip_initial_fit )
	mchi=-1;
 else if ( fp->tune.emce.sub_method==FIT_METHOD_CLLS )
  {	lin_fit_con(fitpnt,fitdep,mvector,fitwght,
	fit_function,nvar,nrow,&ffd,cc,NULL);
	mchi=data_get_chisquare(fitpnt,fitdep,fitwght,NULL,nrow,mvector,&ffd,1);
  }
 else if ( fp->tune.emce.sub_method==FIT_METHOD_NLLM )
  {	double	lam;
	lam=fp->tune.nllm.lambda;
	for ( i=0 ; i<fp->tune.nllm.max_iter ; i++ )
	 {	lam=nlm_fit_base_con(fitpnt,fitdep,mvector,fitwght,
		fit_function,nvar,nrow,&ffd,cc,
		lam,fp->tune.nllm.lambda_mpy);
	 }
	mchi=data_get_chisquare(fitpnt,fitdep,fitwght,NULL,nrow,mvector,&ffd,1);
  }
 else if ( fp->tune.emce.sub_method==FIT_METHOD_LMND )
  {	double	lam;
	lam=fp->tune.nllm.lambda;
	for ( i=0 ; i<fp->tune.nllm.max_iter ; i++ )
	 {	lam=nlm_fit_nmdf_con(fitpnt,fitdep,mvector,fitwght,
		fit_function,nvar,nrow,&ffd,cc,
		lam,fp->tune.nllm.lambda_mpy,xvector);
	 }
	mchi=data_get_chisquare(fitpnt,fitdep,fitwght,NULL,nrow,mvector,&ffd,1);
  }
 else if ( fp->tune.emce.sub_method==FIT_METHOD_MCMC )
  {	mchi=find_monte_carlo_minimum(fitpnt,fitdep,fitwght,nrow,mvector,&ffd,
		evector,downlink_max_iter,max_decrease_rejection,decrease_error_ratio,
		lf,vars,nvar);
  }
 else /* FIT_METHOD_DHSX: default */
  {	mchi=find_downhill_simplex_minimum(fitpnt,fitdep,fitwght,nrow,
		mvector,&ffd,ematrix,nvar,nvar-nc,NULL,0);
  }

 data_get_scatters(fitpnt,fitdep,fitwght,fiteval,
	nrow,mvector,&ffd,0);

 if ( use_dft_rednoise )
	resdft=data_get_noise_dft(fitpnt,fitdep,fitwght,fiteval,
		nrow,mvector,&ffd);
 else
	resdft=NULL;
 
 fw=off->f_write;

 if ( fw != NULL )
  {	fprintf(fw,"# Initial result (nonperturbed, original data):\n");
	if ( mchi>=0.0 )
		fwrite_vars_chi2(fw,mvector,vars,nvar,mchi);
	else
	 {	fwrite_vars_chi2(fw,mvector,vars,nvar,0.0);
		fprintf(fw,"# (explicitly untouched initial input values)\n");
	 }
	fprintf(fw,"# Scatters (RMS, reduced \\chi^2): \n");
	for ( i=0 ; i<lf->ndatablock ; i++ )
	 {	fprintf(fw,"# [%s] %12g %12.8f\n",
			lf->datablocks[i].key,
			lf->datablocks[i].dbscatter,
			lf->datablocks[i].dbredchi2);
	 }
	if ( fp->mc_iterations>0 )
	 {	fprintf(fw,"# Additional fits:\n");		}
	fflush(fw);
  }

 for ( n=0 ; n<fp->mc_iterations ; n++ )
  {	
	if ( resdft != NULL )
		data_perturb_dft(fitpnt,fiteval,fitvary,nrow,&ffd,resdft);
	else
		data_perturb_wtn(fitpnt,fiteval,fitvary,nrow,&ffd);

	for ( i=0 ; i<nvar ; i++ )
	 {	evector[i]=vars[i].ierr;
		bvector[i]=mvector[i];
	 }

	if ( fp->tune.emce.sub_method==FIT_METHOD_CLLS )
	 {	lin_fit_con(fitpnt,fitvary,bvector,fitwght,fit_function,nvar,nrow,&ffd,cc,NULL);
		nchi=data_get_chisquare(fitpnt,fitvary,fitwght,NULL,nrow,bvector,&ffd,1);
	 }
	else if ( fp->tune.emce.sub_method==FIT_METHOD_NLLM )
	 {	double	lam;
		lam=fp->tune.nllm.lambda;
		for ( i=0 ; i<fp->tune.nllm.max_iter ; i++ )
		 {	lam=nlm_fit_base_con(fitpnt,fitvary,bvector,fitwght,
			fit_function,nvar,nrow,&ffd,cc,lam,fp->tune.nllm.lambda_mpy);
		 }
		nchi=data_get_chisquare(fitpnt,fitvary,fitwght,NULL,nrow,bvector,&ffd,1);
	 }
	else if ( fp->tune.emce.sub_method==FIT_METHOD_LMND )
	 {	double	lam;
		lam=fp->tune.nllm.lambda;
		for ( i=0 ; i<fp->tune.nllm.max_iter ; i++ )
		 {	lam=nlm_fit_nmdf_con(fitpnt,fitvary,bvector,fitwght,
			fit_function,nvar,nrow,&ffd,cc,
			lam,fp->tune.nllm.lambda_mpy,xvector);
		 }
		nchi=data_get_chisquare(fitpnt,fitvary,fitwght,NULL,nrow,bvector,&ffd,1);
	 }
	else if ( fp->tune.emce.sub_method==FIT_METHOD_MCMC )
	 {	nchi=find_monte_carlo_minimum(fitpnt,fitvary,fitwght,nrow,bvector,&ffd,
			evector,downlink_max_iter,max_decrease_rejection,decrease_error_ratio,
			lf,vars,nvar);
	 }	
	else /* FIT_METHOD_DHSX: default */
	 {	nchi=find_downhill_simplex_minimum(fitpnt,fitvary,fitwght,nrow,
			bvector,&ffd,ematrix,nvar,nvar-nc,NULL,0);
	 }
/*
	nchi=find_monte_carlo_minimum(fitpnt,fitvary,fitwght,nrow,bvector,&ffd,
		evector,downlink_max_iter,max_decrease_rejection,decrease_error_ratio,
		lf,vars,nvar);
*/

	if ( fw != NULL )
	 {	fwrite_vars_chi2(fw,bvector,vars,nvar,nchi);
		fflush(fw);
	 }
  }

 if ( resdft != NULL )		free(resdft);
	
 if ( fiteval != NULL )		free(fiteval);
 if ( fitpnt  != NULL )		free(fitpnt);
 if ( fitvary  != NULL )	free(fitvary);
 if ( fitinputrows != NULL )
  {	for ( i=0 ; i<nrow ; i++ )
	 {	free(fitinputrows[i].line);	}
	free(fitinputrows);
  }
 if ( fitwght  != NULL )	free(fitwght);
 if ( fitdep   != NULL )	free(fitdep);
 if ( fitdata  != NULL )	free(fitdata);

 vector_free(xvector);
 vector_free(mvector);
 if ( ematrix != NULL )		free(ematrix);
 vector_free(evector);
 vector_free(dvector);
 vector_free(bvector);

 return(0);
}

/*****************************************************************************/

int fprint_matrix(FILE *fw,double **matrix,int n,char *prefix)
{
 int	i,j;
 for ( i=0 ; i<n ; i++ )
  {	if ( prefix != NULL )
	 {	fprintf(fw,"%s",prefix);		}
	for ( j=0 ; j<n ; j++ )
	 {	fprintf(fw," %10.4g",matrix[i][j]);	}
  	fprintf(fw,"\n");
  }
 return(0);
}

static double xmmc_acceptances[]=
{
	1.00000000, /*  0 */
        0.70485636, /*  1 */
        0.55268644, /*  2 */
        0.45024928, /*  3 */
        0.37386438, /*  4 */
        0.31446319, /*  5 */
        0.26662749, /*  6 */
        0.22743129, /*  7 */
        0.19500379, /*  8 */
        0.16785605, /*  9 */
        0.14492830, /* 10 */
        0.12545244, /* 11 */
        0.10881852, /* 12 */
        0.09456381, /* 13 */
        0.08238798, /* 14 */
        0.07191076, /* 15 */
        0.06276243, /* 16 */
        0.05490247, /* 17 */
        0.04799478, /* 18 */
        0.04209247, /* 19 */
        0.03680494, /* 20 */
        0.03233213, /* 21 */
        0.02838227, /* 22 */
        0.02490090, /* 23 */
        0.02184710, /* 24 */
        0.01915704, /* 25 */
        0.01688550, /* 26 */
        0.01482915, /* 27 */
        0.01303621, /* 28 */
        0.01140750, /* 29 */
        0.00988681  /* 30 */
};

double xmmc_theoretical_acceptance(int nvar)
{
 if ( nvar<0 )
	return(-1.0);
 else if ( nvar<=30 )
	return(xmmc_acceptances[nvar]);
 else
	return(0.0);	/* TBD */
}	

double get_correlation_length(double *arr,int win)
{
 double	half,delta;
 int	i;

 if ( arr==NULL || win<1 )
	return(-1.0);
 else if ( arr[0]<=0.0 )
	return(0.0);

 half=arr[0]/2.0;
 for ( i=1 ; i<win ; i++ )
  {	if ( arr[i]<=half )
	 {	delta=(arr[i-1]-half)/(arr[i-1]-arr[i]);
		return((double)(i-1)+delta);
	 }
  }
 return(-1.0);
}

/* fi	t_extended_markov_chain_mc():
   This function is the core of the extended Markov-chain MC (XMMC) method.  */
int fit_extended_markov_chain_mc(lfitdata *lf,fitout *off,variable *vars,int nvar)
{
 int		i,j,n,nc,nrc,ccnt,is_not_accepted;

 void		**fitpnt;
 double		*fitdata,*fitdep,*fiteval,*fitwght;
 double		*wvars,chi2,nchi,mchi,vx,u,alpha,dalpha;
 fitinputrow	*fitinputrows;
 int		nrow,window,corr_c,iiter,niter;
 double		*bvector,*dvector,*mvector,*gvector,*tvector;
 double		**rcmatrix,**cmatrix;
 double		**dcmatrix,*dcvector,*ddvector,*dsvector;
 double		**corrlendd,**corr_ftmp,**corr_prev,**corr_carr,*corr_sum1,*corr_sum2;
 fitfunctdata	ffd;
 FILE		*fw;
 fitparam	*fp;
 int		nseparated,*idx;
 fitconstraint	ifc;

 fp=&lf->parameters;

 nc=lf->fconstraint.nc;	/* number of constraints */
 nrc=nc;	/* number of "real" constraints, i.e. not explicit fixations */
 nseparated=0;
 for ( i=0 ; i<nvar ; i++ )	
  {	if ( vars[i].flags & VAR_IS_CONSTANT )
		nrc--;
	if ( vars[i].is_separated )
		nseparated++;
  }
 
 fitdata =NULL;
 fitdep  =NULL;
 fitwght =NULL;
 fitinputrows=NULL;
 nrow=0;

 bvector=vector_alloc(nvar);	/* set of unknown variables	*/
 dvector=vector_alloc(nvar);	/* perturbed variable vector	*/
 mvector=vector_alloc(nvar);	/* vector of minimal values	*/
 gvector=vector_alloc(nvar);	/* gaussian random variables	*/
 tvector=vector_alloc(nvar);	/* temporary vector		*/

 dcmatrix=matrix_alloc(nvar);
 dcvector=vector_alloc(nvar);
 ddvector=vector_alloc(nvar);
 dsvector=vector_alloc(nvar);

 window=lf->parameters.tune.xmmc.window;
 corrlendd=matrix_alloc_gen(window,3*nvar);
 corr_ftmp=corrlendd+0*nvar;
 corr_prev=corrlendd+1*nvar;
 corr_carr=corrlendd+2*nvar;
 corr_sum1=vector_alloc(nvar);
 corr_sum2=vector_alloc(nvar);

 read_all_fit_data(lf,nvar,&fitdata,&fitdep,&fitwght,&fitinputrows,&nrow);

 if ( nrow < nvar-nc )	
  {	fprint_error("too few lines for fitting");
	exit(1);
  }

 wvars=(double *)malloc(sizeof(double)*(nvar+lf->maxncol));

 fitpnt=(void **)malloc(nrow*sizeof(void *));
 for ( i=0 ; i<nrow ; i++ )
  {	fitinputrows[i].x=&fitdata[i*lf->maxncol];
	fitpnt[i]=(void *)(&fitinputrows[i]);
  }

 ffd.nvar=nvar;
 ffd.wvars=wvars;
 ffd.functs=lpg.pl_funct;
 ffd.lf=lf;

 fiteval=(double *)malloc(sizeof(double)*nrow);

 for ( i=0 ; i<nvar ; i++ )
  {	bvector[i]=vars[i].init;			}

 fw=off->f_write;

 nseparated=get_separate_index(vars,nvar,&idx);

 /* now we create the internal fit constraint (ifc) matrix & projector: */
 ifc.nc=nc+nseparated;
 ifc.cmatrix=matrix_alloc(nvar);
 /* explicit constraints are inherited: */
 for ( i=0 ; i<nc ; i++ )
  {	for ( j=0 ; j<nvar ; j++ )
	 {	ifc.cmatrix[i][j]=lf->fconstraint.cmatrix[i][j];		}
  }
 /* linear terms imply other constraints which contribute to the projector: */
 for ( i=0 ; i<nseparated ; i++ )
  {	for ( j=0 ; j<nvar ; j++ )
 	 {	ifc.cmatrix[nc+i][j]=0.0;		}
	ifc.cmatrix[nc+i][idx[i]]=1.0;
  }
 ifc.cproject=matrix_alloc(nvar);
 constraint_initialize_proj_matrix(ifc.nc,nvar,ifc.cmatrix,ifc.cproject);

 if ( ! lf->parameters.tune.xmmc.skip_initial_fit )
  {	double	**ematrix;
	double	**icmatrix;

	ematrix=matrix_alloc(nvar);
	icmatrix=data_get_inversecovariance(fitpnt,fitdep,fitwght,nrow,nvar,bvector,&ffd);
	get_error_matrix(ematrix,nvar,NULL,icmatrix,ifc.nc,ifc.cproject);
	matrix_free(icmatrix);
	find_downhill_simplex_minimum(fitpnt,fitdep,fitwght,nrow,bvector,&ffd,ematrix,nvar,nvar-ifc.nc,idx,nseparated);
	matrix_free(ematrix);
  }

 if ( ! lf->parameters.tune.xmmc.is_adaptive )
	rcmatrix=data_get_rootcovariance(fitpnt,fitdep,fitwght,nrow,nvar,bvector,&ffd,ifc.nc,ifc.cproject);
 else
	rcmatrix=NULL;

/*
 fprint_matrix(stderr,lf->fconstraint.cproject,nvar,"# ");
 fprint_matrix(stderr,rcmatrix,nvar,"> ");
*/

/*
   if we have constraints, the variables are pushed to the nearest point
   of the intersections of the hyperplanes defined by the constraints
   (i.e. the variables are forced to satisfy the constraints with the 
   minimal shift distance).
 */

 if ( nc > 0 )
  {	for ( i=0 ; i<nvar ; i++ )
	 {	if ( vars[i].flags & VAR_IS_CONSTANT )
	 	 {	bvector[i]=vars[i].init;		}
	 }
	if ( nrc>0 )
	 {	constraint_force(bvector,nvar,nc,
		lf->fconstraint.invlmatrix,lf->fconstraint.cvector);
	 }
  }

 chi2=data_get_chisquare_separated_linear(fitpnt,fitdep,fitwght,fiteval,nrow,bvector,&ffd,0,nvar,idx,nseparated,bvector);
 mchi=chi2;
 for ( i=0 ; i<nvar ; i++ )
  {	 mvector[i]=bvector[i];
	ddvector[i]=bvector[i];
  }

 if ( ! lf->parameters.tune.xmmc.skip_initial_fit && fw != NULL )
  {	fprintf(fw,"# Downhill simplex best fit value:\n");
	fprintf(fw," ");
	fwrite_vars_chi2(fw,bvector,vars,nvar,chi2);
	fflush(fw);
  }

 niter=(lf->parameters.tune.xmmc.niter>0?lf->parameters.tune.xmmc.niter:0);

 for ( iiter=0 ; iiter <= niter ; iiter++ )
  {
	if ( fw != NULL && niter>0 )
	 {	fprintf(fw,"# XMMC: iteration %d:\n",iiter);
		fflush(fw);
	 }
	if ( fw != NULL )
	 {	fprintf(fw,"# XMMC values:\n");
		fflush(fw);
	 }

	for ( i=0 ; i<nvar ; i++ )
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	dcmatrix[i][j]=0.0;		}
		dsvector[i]=0.0;
	 }

	for ( i=0 ; i<nvar ; i++ )
	 {	for ( j=0 ; j<window ; j++ )
		 {	corr_ftmp[i][j]=0.0;
			corr_prev[i][j]=0.0;
			corr_carr[i][j]=0.0;
		 }
		corr_sum1[i]=0.0;
		corr_sum2[i]=0.0;
	 }
	corr_c=0;

	/* we have fp->mc_iterations iterations in a single chain: */
	for ( n=0,ccnt=0 ; (lf->parameters.tune.xmmc.count_accepted?ccnt:n)<fp->mc_iterations ; n++ )
	 {	
		/* we assume accaptance, i.e. not out-of-boundary data */
		is_not_accepted=0;
		for ( j=0 ; j<nvar ; j++ )
		 {	gvector[j]=get_gaussian(0.0,1.0);		}

		/* this is the main point of the XMMC method: */
		if ( lf->parameters.tune.xmmc.is_adaptive )
			rcmatrix=data_get_rootcovariance(fitpnt,fitdep,fitwght,nrow,nvar,bvector,&ffd,ifc.nc,ifc.cproject);

		for ( i=0 ; i<nvar ; i++ )
		 {	dvector[i]=bvector[i];
			for ( j=0 ; j<nvar ; j++ )
			 {	dvector[i]+=rcmatrix[i][j]*gvector[j];		}
		 }

		if ( lf->parameters.tune.xmmc.is_adaptive )
		 {	matrix_free(rcmatrix);
			rcmatrix=NULL;
		 }
	
		for ( i=0 ; i<nvar ; i++ )
		 {
			if ( (vars[i].flags & VAR_MIN) && dvector[i] < vars[i].imin )
			 {	is_not_accepted=1; /* o-o-b */
				dvector[i]=vars[i].imin;
			 }
			if ( (vars[i].flags & VAR_MAX) && dvector[i] > vars[i].imax )
			 {	is_not_accepted=1; /* o-o-b */
				dvector[i]=vars[i].imax;
			 }
		 }

		for ( i=0 ; i<lf->dconstraint.ndcfunct && (!is_not_accepted) ; i++ )
		 {	double	w;
			lfit_psn_double_calc(lf->dconstraint.dcfuncts[i].funct,NULL,lpg.pl_funct,&w,dvector);
			if ( w<=0.0 )
				is_not_accepted=1;
		 }
	
		/* the next point is not accepted (i.e. it is an out of bound point) */
		if ( is_not_accepted )
		 {	/* fprintf(stderr,"qqriq\n"); */
			continue;
		 }

		/* again, constraints are forced to the variables: */	
		if ( nc > 0 )
		 {	for ( i=0 ; i<nvar ; i++ )
			 {	if ( vars[i].flags & VAR_IS_CONSTANT )
			 	 {	dvector[i]=vars[i].init;		}
			 }
			if ( nrc>0 )
			 {	constraint_force(dvector,nvar,nc,
				lf->fconstraint.invlmatrix,lf->fconstraint.cvector);
			 }
		 }

		/* we get a new chi square value */
		nchi=data_get_chisquare_separated_linear(fitpnt,fitdep,fitwght,fiteval,nrow,dvector,&ffd,0,nvar,idx,nseparated,dvector);
		/* note: the above call, if nseparated>0 does likely change some elements in dvector[]. */

		/* just save it if it is minimal (useful as the fitted values) and   */
		/* to normalize the chisquare function (i.e. multiply the errors by  */
		/* some values if desired...)					     */
		if ( nchi<mchi )
		 {	mchi=nchi;
			for ( i=0 ; i<nvar ; i++ )
			 {	mvector[i]=dvector[i];		}
		 }
	
		vx=exp(-0.5*(nchi-chi2));
		alpha=(vx<1.0?vx:1.0);

		/* it is accepted, good: */
		if ( alpha>0.0 && (u=drand48()) <= alpha )
		 {	for ( i=0 ; i<nvar ; i++ )
			 {	bvector[i]=dvector[i];
				dcvector[i]=bvector[i]-ddvector[i];
			 }
			chi2=nchi;

			for ( i=0 ; i<nvar ; i++ )
			 {	double	d;
				d=dcvector[i];
				if ( ccnt<window )
				 {	corr_ftmp[i][ccnt]=d;		}
				corr_prev[i][corr_c]=d;
				for ( j=0 ; j<window && j<=ccnt ; j++ )
				 {	corr_carr[i][j]+=d*corr_prev[i][(corr_c+window-j)%window];	}
				corr_sum1[i]+=d;
				corr_sum2[i]+=d*d;
			 }
			corr_c=(corr_c+1)%window;

			/* increase the acceptance counter: */
			ccnt++;

			/* write a posteriori distr. element: */
			if ( fw != NULL )
			 {	if ( iiter<niter )
					fprintf(fw,"#");
				else
					fprintf(fw," ");

			 	fwrite_vars_chi2(fw,bvector,vars,nvar,chi2);
				fflush(fw);
			 }

			for ( i=0 ; i<nvar ; i++ )
			 {	for ( j=0 ; j<nvar ; j++ )
				 {	dcmatrix[i][j]+=dcvector[i]*dcvector[j];	}
				dsvector[i]+=dcvector[i];
			 }
		 }
		/* not accepted, set is_not_accepted to true is a symbolic thing...  */
		else
			is_not_accepted=1;	
	 }

	for ( i=0 ; i<nvar ; i++ )
	 {	double	s1,s2,sig2;
		int	j,n;
	
	 	if ( ccnt>0 )
		 {	corr_sum1[i]/=(double)ccnt;
			corr_sum2[i]/=(double)ccnt;
		 }
		s1=corr_sum1[i];
		s2=corr_sum2[i];

		sig2=s2-s1*s1;
		if ( sig2<0.0 )	sig2=0.0;
		
		for ( n=1 ; n<window ; n++ )
		 {	for ( j=0 ; j<n ; j++ )
			 {	corr_carr[i][n]+=corr_ftmp[i][j]*corr_prev[i][(corr_c+window-n)%window];	}
		 }
		for ( n=0 ; n<window ; n++ )
		 {	if ( ccnt>0 && sig2>0.0 )
				corr_carr[i][n]=(corr_carr[i][n]-s1*s1*ccnt)/(ccnt*sig2);
			else
				corr_carr[i][n]=0.0;
		 }
	 }

	for ( i=0 ; i<nvar && ccnt>0 ; i++ )
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	dcmatrix[i][j]/=(double)ccnt;		}
		dsvector[i]/=(double)ccnt;
	 }	
	for ( i=0 ; i<nvar ; i++ )
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	dcmatrix[i][j] -= dsvector[i]*dsvector[j];	}
	 }

	if ( fw != NULL )
	 {	char	*fstr;
		double	alpha0;

		if ( n>0 )
		 {	alpha=(double)ccnt/(double)n;
			dalpha=sqrt((double)ccnt)/(double)n;
		 }
		else
		 {	alpha=0.0;
			dalpha=0.0;
		 }
		fprintf(fw,"#\n");
		fprintf(fw,"# Accepted transitions / total iterations: %d/%d\n",ccnt,n);
		fprintf(fw,"# Total acceptance ratio : %.5f +/- %.5f\n",alpha,dalpha);
		alpha0=xmmc_theoretical_acceptance(nvar-nc-nseparated);	
		fprintf(fw,"# Theoretical probability: %.5f [independent:%d=%d-%d-%d (total-constrained-linear)]\n",alpha0,nvar-nc-nseparated,nvar,nc,nseparated);
		fprintf(fw,"#\n");

		fprintf(fw,"# Correlation lengths: \n");
		fprintf(fw,"#   ");
		for ( i=0 ; i<nvar ; i++ )
		 {	double	clen;
			clen=get_correlation_length(corr_carr[i],window);
			fprintf(fw," %7.2f",clen);
		 }
		fprintf(fw,"\n");
		fprintf(fw,"#\n");

		fprintf(fw,"# chi^2 values:\n");
		mchi=data_get_chisquare_separated_linear(fitpnt,fitdep,fitwght,fiteval,nrow,mvector,&ffd,1,nvar,idx,nseparated,mvector);
		fprintf(fw,"#\tminimal: %g\n",mchi);
		fprintf(fw,"# Appropriate values for this chi^2:\n# ");
		fwrite_vars(fw,mvector,vars,nvar);
		fprintf(fw,"\n#\n");

		if ( lf->fconstraint.nc>0 )	fstr="projected Fisher";
		else				fstr="Fisher";	
		fprintf(fw,"# Errors and correlations (%s):\n",fstr);
	
		cmatrix=data_get_covariance(fitpnt,fitdep,fitwght,nrow,nvar,mvector,&ffd);
		for ( i=0 ; i<nvar ; i++ )
		 {	tvector[i]=sqrt(cmatrix[i][i]);		}
		fprintf(fw,"# ");fwrite_vars(fw,tvector,vars,nvar);fprintf(fw,"\n#\n");
		for ( i=0 ; i<nvar ; i++ )
		 {	fprintf(fw,"# ");
			for ( j=0 ; j<nvar ; j++ )
			 {	if ( tvector[i]>0.0 && tvector[j]>0.0 )
					u=cmatrix[i][j]/(tvector[i]*tvector[j]);
				else
					u=(i==j?1.0:0.0);
				fprintf(fw," ");
				fprintf(fw,lf->corrfm,u);
			 }
			fprintf(fw,"\n");
		 }
		matrix_free(cmatrix);
		fprintf(fw,"#\n");

		fprintf(fw,"# Errors and correlations (statistical, around the best fit):\n");

		for ( i=0 ; i<nvar ; i++ )
		 {	tvector[i]=sqrt(dcmatrix[i][i]);		}
		fprintf(fw,"# ");fwrite_vars(fw,tvector,vars,nvar);fprintf(fw,"\n#\n");
		for ( i=0 ; i<nvar ; i++ )
		 {	fprintf(fw,"# ");
			for ( j=0 ; j<nvar ; j++ )
			 {	if ( tvector[i]>0.0 && tvector[j]>0.0 )
					u=dcmatrix[i][j]/(tvector[i]*tvector[j]);
				else
					u=(i==j?1.0:0.0);
				fprintf(fw," ");
				fprintf(fw,lf->corrfm,u);
			 }
			fprintf(fw,"\n");
		 }
		fprintf(fw,"#\n");
	 }

	/* this is the main point of the iterative XMMC: we update here      */
	/* the 'rcmatrix' which generetes the a priori transitions. In the   */
	/* first (zeroth) iteration, this matrix is derived from the Fisher  */
	/* infromation matrix, otherwise (later) it is simply the root of    */
	/* the covariance matrix (which is 'dcmatrix' here): 		     */
	if ( iiter<niter && rcmatrix != NULL )
	 {	for ( i=0 ; i<nvar ; i++ )
		 {	for ( j=0 ; j<nvar ; j++ )
			 {	rcmatrix[i][j]=dcmatrix[i][j];		}
		 }
		matrix_root(rcmatrix,nvar,NULL);
	 }

	/* just write a separator between the consecutive iterations: */
	if ( fw != NULL && iiter<niter )
	 {	char	buff[80];
		memset(buff,'#',79);
		buff[79]=0;
		fprintf(fw,"%s\n",buff);
	 }

  }

 if ( ! lf->parameters.tune.xmmc.is_adaptive )
	matrix_free(rcmatrix);

 if ( fiteval != NULL )		free(fiteval);
 if ( fitpnt  != NULL )		free(fitpnt);
 if ( fitinputrows != NULL )
  {	for ( i=0 ; i<nrow ; i++ )
	 {	free(fitinputrows[i].line);	}
	free(fitinputrows);
  }
 if ( fitwght  != NULL )	free(fitwght);
 if ( fitdep   != NULL )	free(fitdep);
 if ( fitdata  != NULL )	free(fitdata);

 vector_free(corr_sum2);
 vector_free(corr_sum1);
 matrix_free(corrlendd);

 if ( ifc.cmatrix  != NULL )	matrix_free(ifc.cmatrix);
 if ( ifc.cproject != NULL )	matrix_free(ifc.cproject);

 vector_free(dsvector);
 vector_free(ddvector);
 vector_free(dcvector);
 matrix_free(dcmatrix);

 vector_free(tvector);
 vector_free(gvector);
 vector_free(mvector);
 vector_free(dvector);
 vector_free(bvector);

 return(0);
}

/*****************************************************************************/

/* fit_fisher_matrix_analysis():
   This is not a real fit method but it provides a complete Fisher analysis: */
int fit_fisher_matrix_analysis(lfitdata *lf,fitout *off,
	variable *vars,int nvar,dvariable *dvars,int ndvar)
{
 int		i,j,n,nc,nrc;

 void		**fitpnt;
 double		*fitdata,*fitdep,*fitwght;
 double		*wvars,chi2,mchi;
 fitinputrow	*fitinputrows;
 int		nrow;
 double		*bvector,*dvector,*gvector;
 double		**rcmatrix;
 fitfunctdata	ffd;
 FILE		*fw;
 fitparam	*fp;

 if ( ndvar<=0 || dvars==NULL )
	ndvar=0,dvars=NULL;

 fp=&lf->parameters;

 nc=lf->fconstraint.nc;	/* number of constraints */
 nrc=nc;	/* number of "real" constraints, i.e. not explicit fixations */
 for ( i=0 ; i<nvar ; i++ )	
  {	if ( vars[i].flags & VAR_IS_CONSTANT )
		nrc--;
  }

 fitdata =NULL;
 fitdep  =NULL;
 fitwght =NULL;
 fitinputrows=NULL;
 nrow=0;

 bvector=vector_alloc(nvar+ndvar);	/* set of unknown variables	*/
 dvector=vector_alloc(nvar);		/* perturbed variable vector	*/
 gvector=vector_alloc(nvar);		/* gaussian random variables	*/
 
 read_all_fit_data(lf,nvar,&fitdata,&fitdep,&fitwght,&fitinputrows,&nrow);

 wvars=(double *)malloc(sizeof(double)*(nvar+lf->maxncol));

 fitpnt=(void **)malloc(nrow*sizeof(void *));
 for ( i=0 ; i<nrow ; i++ )
  {	fitinputrows[i].x=&fitdata[i*lf->maxncol];
	fitpnt[i]=(void *)(&fitinputrows[i]);
  }

 ffd.nvar=nvar;
 ffd.wvars=wvars;
 ffd.functs=lpg.pl_funct;
 ffd.lf=lf;

 for ( i=0 ; i<nvar ; i++ )
  {	bvector[i]=vars[i].init;			}
 for ( i=0 ; i<ndvar ; i++ )
  {	double	d;
	lfit_psn_double_calc(dvars[i].funct,NULL,lpg.pl_funct,&d,bvector);
	bvector[i+nvar]=d;
  }

 fw=off->f_write;

 mchi=data_get_chisquare(fitpnt,fitdep,fitwght,NULL,nrow,bvector,&ffd,1);
 data_get_scatters(fitpnt,fitdep,fitwght,NULL,nrow,bvector,&ffd,0);

 if ( fw != NULL )
  {	fprintf(fw,"#");
	fprint_variable_name_list(fw,vars,nvar);
	fprint_derived_variable_name_list(fw,dvars,ndvar);
	fprintf(fw,"\n");
  	fprintf(fw,"#");
	fprint_variable_name_index(fw,vars,nvar);
	fprint_derived_variable_name_index(fw,dvars,ndvar,nvar);
	fprintf(fw,"\n");
  }

 if ( fw != NULL && fp->tune.fima.write_original )
  {	fprintf(fw,"# Fisher analysis - original data:\n");
	if ( fp->tune.fima.do_montecarlo )
		fprintf(fw,"# ");
	else
		fprintf(fw,"  ");
	fwrite_vars(fw,bvector,vars,nvar);
	fwrite_dvars(fw,bvector+nvar,dvars,ndvar);
	fprintf(fw,"\n");
	fprintf(fw,"# Scatters (RMS, reduced \\chi^2): \n");
	for ( i=0 ; i<lf->ndatablock ; i++ )
	 {	fprintf(fw,"# [%s] %12g %12.8f\n",
			lf->datablocks[i].key,
			lf->datablocks[i].dbscatter,
			lf->datablocks[i].dbredchi2);
	 }
	fflush(fw);
  }

 rcmatrix=data_get_rootcovariance(fitpnt,fitdep,fitwght,nrow,nvar,bvector,&ffd,lf->fconstraint.nc,lf->fconstraint.cproject);

 for ( n=0 ; n<fp->mc_iterations && fp->tune.fima.do_montecarlo ; n++ )
  {	for ( j=0 ; j<nvar ; j++ )
	 {	gvector[j]=get_gaussian(0.0,1.0);		}
	for ( i=0 ; i<nvar ; i++ )
	 {	dvector[i]=bvector[i];
		for ( j=0 ; j<nvar ; j++ )
		 {	dvector[i]+=rcmatrix[i][j]*gvector[j];		}
	 }
	if ( nc > 0 )
	 {	for ( i=0 ; i<nvar ; i++ )
		 {	if ( vars[i].flags & VAR_IS_CONSTANT )
		 	 {	dvector[i]=vars[i].init;		}
		 }
		if ( nrc>0 )
		 {	constraint_force(dvector,nvar,nc,
			lf->fconstraint.invlmatrix,lf->fconstraint.cvector);
		 }
	 }
	if ( fw != NULL )
	 {	chi2=data_get_chisquare(fitpnt,fitdep,fitwght,NULL,nrow,dvector,&ffd,1);
		fwrite_vars_chi2(fw,dvector,vars,nvar,chi2);
	 }
  }

 matrix_free(rcmatrix);

 if ( fw != NULL )
  {	char	*fstr;
	double	**cmatrix,**ocmatrix,*tvector;
	double	u;

	if ( lf->fconstraint.nc>0 )	fstr="projected Fisher";
	else				fstr="Fisher";	

	if ( fp->tune.fima.write_uncert )
		fprintf(fw,"# Uncertainties (%s):\n",fstr);

	tvector=vector_alloc(nvar+ndvar);
	cmatrix=matrix_alloc(nvar+ndvar);
	ocmatrix=data_get_covariance(fitpnt,fitdep,fitwght,nrow,nvar,bvector,&ffd);

	if ( ocmatrix != NULL )
	 {	int	k,l,m,n;
		double	c,p1,p2,**diff;

		diff=matrix_alloc(nvar+ndvar);
		for ( k=0 ; k<ndvar ; k++ )
		 {	for ( m=0 ; m<nvar ; m++ )
			 {	lfit_psn_double_calc(dvars[k].diff[m],NULL,lpg.pl_funct,&p1,bvector);
				diff[k][m]=p1;
			 }
		 }

		for ( k=0 ; k<nvar+ndvar ; k++ )
		 {  for ( l=0 ; l<=k ; l++ )
		     {	c=0.0;
			for ( m=0 ; m<nvar ; m++ )
			 {  for ( n=0 ; n<nvar ; n++ )
			     {	if ( k<nvar && k==m )	p1=1.0;
				else if ( k<nvar )	p1=0.0;
				else			p1=diff[k-nvar][m];
				if ( l<nvar && l==n )	p2=1.0;
				else if ( l<nvar )	p2=0.0;
				else			p2=diff[l-nvar][n];
				c+=p1*ocmatrix[m][n]*p2;
			     }
			 }
			cmatrix[k][l]=c;
		     }
		 }
		for ( k=0 ; k<nvar+ndvar ; k++ )
		 {  for ( l=k+1 ; l<nvar+ndvar ; l++ )
		     {	cmatrix[k][l]=cmatrix[l][k];	}
		 }

		matrix_free(diff);
		matrix_free(ocmatrix);

	 }
	else
		cmatrix=NULL;

	if ( cmatrix==NULL )
	 {	for ( i=0 ; i<ndvar ; i++ )
		 {	tvector[i]=0.0;				}
		if ( fp->tune.fima.do_montecarlo )
			fprintf(fw,"# ");
		else
			fprintf(fw,"  ");
		if ( fp->tune.fima.write_uncert )
		 {	fwrite_vars(fw,tvector,vars,nvar);fprintf(fw,"\n");	}
	 }
	else
	 {	for ( i=0 ; i<nvar+ndvar ; i++ )
		 {	tvector[i]=sqrt(cmatrix[i][i]);		}

		if ( fp->tune.fima.write_uncert )
		 {	if ( fp->tune.fima.do_montecarlo )
				fprintf(fw,"# ");
			else
				fprintf(fw,"  ");
			fwrite_vars(fw,tvector,vars,nvar);
			fwrite_dvars(fw,tvector+nvar,dvars,ndvar);
			fprintf(fw,"\n");
		 }

		if ( fp->tune.fima.write_correl )
		 {	fprintf(fw,"# Correlations:\n");

			for ( i=0 ; i<nvar+ndvar ; i++ )
			 {	if ( fp->tune.fima.do_montecarlo )
					fprintf(fw,"# ");
				for ( j=0 ; j<nvar+ndvar ; j++ )
				 {	if ( tvector[i]>0.0 && tvector[j]>0.0 )
						u=cmatrix[i][j]/(tvector[i]*tvector[j]);
					else
						u=(i==j?1.0:0.0);
					fprintf(fw," ");
					fprintf(fw,lf->corrfm,u);
				 }
				fprintf(fw,"\n");
			 }
		 }

		matrix_free(cmatrix);
	 }
	vector_free(tvector);
  }

 if ( fitpnt  != NULL )		free(fitpnt);
 if ( fitinputrows != NULL )
  {	for ( i=0 ; i<nrow ; i++ )
	 {	free(fitinputrows[i].line);	}
	free(fitinputrows);
  }
 if ( fitwght  != NULL )	free(fitwght);
 if ( fitdep   != NULL )	free(fitdep);
 if ( fitdata  != NULL )	free(fitdata);

 vector_free(gvector);
 vector_free(dvector);
 vector_free(bvector);

 return(0);
}


/*****************************************************************************/

/* fit_downhill_simplex(): */
int fit_downhill_simplex(lfitdata *lf,fitout *off,variable *vars,int nvar)
{
 int		i,nc,nrc;

 void		**fitpnt;
 double		*fitdata,*fitdep,*fiteval,*fitwght;
 double		*wvars,chi2;
 fitinputrow	*fitinputrows;
 int		nrow;
 double		*bvector,*evector;
 double		**ematrix;
 fitfunctdata	ffd;
 FILE		*fw;
 fitparam	*fp;
 int		*idx,nseparated;

 fp=&lf->parameters;

 nc=lf->fconstraint.nc;	/* number of constraints */
 nrc=nc;	/* number of "real" constraints, i.e. not explicit fixations */
 for ( i=0 ; i<nvar ; i++ )	
  {	if ( vars[i].flags & VAR_IS_CONSTANT )
		nrc--;
  }

 fitdata =NULL;
 fitdep  =NULL;
 fitwght =NULL;
 fitinputrows=NULL;
 nrow=0;

 bvector=vector_alloc(nvar);	/* set of unknown variables	*/
 evector=vector_alloc(nvar);	/* temporary vector		*/
 ematrix=matrix_alloc(nvar);

 read_all_fit_data(lf,nvar,&fitdata,&fitdep,&fitwght,&fitinputrows,&nrow);

 if ( nrow < nvar-nc )	
  {	fprint_error("too few lines for fitting");
	exit(1);
  }

 wvars=(double *)malloc(sizeof(double)*(nvar+lf->maxncol));

 fitpnt=(void **)malloc(nrow*sizeof(void *));
 for ( i=0 ; i<nrow ; i++ )
  {	fitinputrows[i].x=&fitdata[i*lf->maxncol];
	fitpnt[i]=(void *)(&fitinputrows[i]);
  }

 ffd.nvar=nvar;
 ffd.wvars=wvars;
 ffd.functs=lpg.pl_funct;
 ffd.lf=lf;

 fiteval=(double *)malloc(sizeof(double)*nrow);

 for ( i=0 ; i<nvar ; i++ )
  {	bvector[i]=vars[i].init;
	evector[i]=vars[i].ierr;
  }

 ematrix=matrix_alloc(nvar);

 nseparated=get_separate_index(vars,nvar,&idx);

 if ( idx != NULL && 0<nseparated )
  {	fitconstraint	ifc;
	int		i,j;

	ifc.nc=nc+nseparated;
	ifc.cmatrix=matrix_alloc(nvar);
	/* explicit constraints are inherited: */
	for ( i=0 ; i<nc ; i++ )
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	ifc.cmatrix[i][j]=lf->fconstraint.cmatrix[i][j];		}
	 }
	/* linear terms imply other constraints which contribute to the projector: */
	for ( i=0 ; i<nseparated ; i++ )
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	ifc.cmatrix[nc+i][j]=0.0;		}
		ifc.cmatrix[nc+i][idx[i]]=1.0;
	 }
	ifc.cproject=matrix_alloc(nvar);
	constraint_initialize_proj_matrix(ifc.nc,nvar,ifc.cmatrix,ifc.cproject);
	/* well, that's it: */
 	get_error_matrix(ematrix,nvar,evector,NULL,ifc.nc,ifc.cproject);
	matrix_free(ifc.cproject);
	matrix_free(ifc.cmatrix);
  }
 else
 	get_error_matrix(ematrix,nvar,evector,NULL,lf->fconstraint.nc,lf->fconstraint.cproject);

 /* the simplex dimension in the hyperplane is nvar-nc-nseparated, i.e. */
 /* both the explicit constraints and the implied ones (by the linear 	*/
 /* separation) decrease the dimension of the embeddning hyperplane:	*/
 chi2=find_downhill_simplex_minimum(fitpnt,fitdep,fitwght,nrow,bvector,&ffd,ematrix,nvar,nvar-nc-nseparated,NULL,0);

 fw=off->f_write;

 if ( fw != NULL )
  {	fprintf(fw,"# Downhill simplex best fit values:\n");
	fwrite_vars_chi2(fw,bvector,vars,nvar,chi2);
  }

 if ( idx != NULL )	free(idx);

 matrix_free(ematrix);
 vector_free(evector);
 vector_free(bvector);

 return(0);
}

/*****************************************************************************/

int evaluate_and_dump(FILE *fr,FILE *fw,datablock *db,
	variable *vars,int nvar,int *colouts,int ncolout,
	dumpexpr *exprs,int nexpr)
{
 int		i,j,ln,len,is_cont,nt,mxcolout;
 double		*wvars,*lvars,*frets;
 char		*rbuff,**lvarstr;
 int		nrow,cnc;

 if ( db->funct != NULL )
	nt=db->funct->nseq;
 else
	nt=1;

 frets=(double *)malloc(sizeof(double)*nt);

 wvars=(double *)malloc(sizeof(double)*(nvar+db->ncol));
 lvars=(double *)malloc(sizeof(double)*(nvar+db->ncol));

 ln=0;nrow=0;

 mxcolout=0;
 if ( colouts != NULL && ncolout>0 )
  {	for ( i=0 ; i<ncolout ; i++ )
	 {	if ( colouts[i]>mxcolout )
			mxcolout=colouts[i];
	 }
	mxcolout++;
  }

 rbuff=NULL;lvarstr=NULL;
 while ( ! feof(fr) )
  {	ln++;
	if ( rbuff   != NULL )	{ free(rbuff);  rbuff=NULL;   }
	if ( lvarstr != NULL )	{ free(lvarstr);lvarstr=NULL; }

	rbuff=freadline(fr);
	if ( rbuff==NULL )	break;

	len=strlen(rbuff);
	for ( i=0 ; i<len ; i++ )
	 {	if ( rbuff[i]==',' || rbuff[i]==';' )	rbuff[i]=32;
		if ( rbuff[i]=='#' )			rbuff[i]=0,len=i;
	 }
	
	lvarstr=tokenize_spaces_dyn(rbuff);
	if ( lvarstr==NULL || lvarstr[0]==NULL )
	 	continue;
	for ( cnc=0 ; lvarstr[cnc] != NULL ; )	cnc++;

	if ( cnc<db->ncol )
	 {	if ( is_verbose>=0 )
			fprint_warning("missing column in line %d, skipped",ln);
		continue;
	 }
	is_cont=0;
	for ( i=0 ; i<db->ncol && (!is_cont) ; i++ )
	 {	if ( db->cols[i].name != NULL && db->cols[i].name[0] )
		 {	j=sscanf(lvarstr[i],"%lg",&lvars[i]);
			if ( j<1 || (! isfinite(lvars[i])) )
			 {	if ( is_verbose>=0 )
					fprint_warning("inappropriate field '%s' in line %d, skipped",lvarstr[i],ln);
				is_cont=1;
			 }
		 }
		else	lvars[i]=0.0;
	 }
	if ( is_cont )
		continue;

	if ( db->funct != NULL )
	 {	for ( i=0 ; i<nvar ; i++ )
		 {	wvars[i]=vars[i].init;		}
		for ( i=0 ; i<db->ncol ; i++ )
		 {	wvars[i+nvar]=lvars[i];		}
		lfit_psn_double_calc(db->funct,db->fchain,lpg.pl_funct,frets,wvars);
	 }
	else if ( db->lfunct != NULL )
	 {	for ( i=0 ; i<db->lfunct->nvar ; i++ )
		 {	wvars[i]=vars[db->map_var[i]].init;	}
		for ( i=0 ; i<db->lfunct->nidv ; i++ )
		 {	wvars[i+nvar]=lvars[db->map_idv[i]];	}
		db->lfunct->function(wvars,wvars+nvar,frets,NULL);
	 }

	if ( colouts != NULL && ncolout>0  )
	 {	for ( i=0 ; i<cnc || i<mxcolout ; i++ )
		 {	for ( j=0 ; j<ncolout ; j++ )
			 {	if ( colouts[j]==i )	break;		}
			if ( j<ncolout )
			 {	if ( exprs != NULL && j<nexpr )
				 {	fprintf(fw,exprs[j].format,frets[j]);
					fprintf(fw," ");
				 }
				else
					fprintf(fw,LFIT_DEFAULT_FORMAT " ",frets[j]);
			 }
			else if ( i<cnc )	fprintf(fw,"%s\t",lvarstr[i]);
			else			fprintf(fw,"-\t");
		 }
		fprintf(fw,"\n");
	 }
	else
	 {	for ( i=0 ; i<nt ; i++ )
		 {	if ( exprs != NULL && i<nexpr )
			 {	fprintf(fw,exprs[i].format,frets[i]);
				fprintf(fw," ");
			 }
			else
				fprintf(fw,LFIT_DEFAULT_FORMAT " ",frets[i]);
		 }
		fprintf(fw,"\n");
	 }
  }

 free(lvars);
 free(wvars);
 free(frets);

 return(0);
}

/*****************************************************************************/

int lfit_check_symbol_name(char *name)
{
 while ( *name )
  {	if ( ! ( isalnum(*name) || (*name)=='_' ) )
		return(1);
	name++;
  }
 return(0);
}

/*****************************************************************************/

int extract_variables(char *varstr,variable **rvars,int *rnvar)
{
 variable	*vars,*wv;
 int		nvar;
 char		*varstr0=varstr;
 double		init,ierr,imin,imax,istp,w;

 vars=NULL;
 nvar=0;

 while ( *varstr )
  {	vars=(variable *)realloc(vars,sizeof(variable)*(nvar+1));
	wv=&vars[nvar];
	wv->name=varstr;
	wv->flags=0;
	wv->diff=0.0;
	wv->is_linear=0;
	wv->is_separated=0;

	strcpy(wv->format,LFIT_DEFAULT_FORMAT);

	while ( *varstr && *varstr != '=' && *varstr != ',' )	varstr++;
 	if ( *varstr=='=' )
	 {	if ( varstr != varstr0 && *(varstr-1)==':' )
		 {	wv->flags |= VAR_IS_CONSTANT;
			*(varstr-1)=0;
		 }

		*varstr=0,varstr++;
		if ( sscanf(varstr,"[%lg:%lg:%lg]",&imin,&istp,&imax)==3 )
		 {	wv->init=0.0 ,wv->ierr=0.0;
			if ( istp<0.0  ) istp=-istp;
			if ( imin>imax ) w=imin,imin=imax,imax=w;
			wv->imin=imin,wv->imax=imax;
			wv->istp=istp;
			wv->flags |= (VAR_MIN|VAR_STEP|VAR_MAX);
		 }
		else if ( sscanf(varstr,"%lg:%lg[%lg:%lg]",&init,&ierr,&imin,&imax)==4 )
		 {	wv->init=init,wv->ierr=ierr;
			if ( imin>imax ) w=imin,imin=imax,imax=w;
			wv->imin=imin,wv->imax=imax;
			wv->flags |= (VAR_MIN|VAR_MAX|VAR_ERR);
		 }
		else if ( sscanf(varstr,"%lg[%lg:%lg]%lg",&init,&imin,&imax,&ierr)==4 )
		 {	wv->init=init,wv->ierr=ierr;
			if ( imin>imax ) w=imin,imin=imax,imax=w;
			wv->imin=imin,wv->imax=imax;
			wv->flags |= (VAR_MIN|VAR_MAX|VAR_ERR);
		 }
		else if ( sscanf(varstr,"%lg:%lg[:%lg]",&init,&ierr,&imax)==3 )
		 {	wv->init=init,wv->ierr=ierr;
			wv->imin=0.0 ,wv->imax=imax;
			wv->flags |= (VAR_MAX|VAR_ERR);
		 }
		else if ( sscanf(varstr,"%lg:%lg[%lg:]",&init,&ierr,&imin)==3 )
		 {	wv->init=init,wv->ierr=ierr;
			wv->imin=imin,wv->imax=0.0 ;
			wv->flags |= (VAR_MIN|VAR_ERR);
		 }
		else if ( sscanf(varstr,"%lg[:%lg]:%lg",&init,&imax,&ierr)==3 )
		 {	wv->init=init,wv->ierr=ierr;
			wv->imin=0.0 ,wv->imax=imax;
			wv->flags |= (VAR_MAX|VAR_ERR);
		 }
		else if ( sscanf(varstr,"%lg[%lg:]:%lg",&init,&imin,&ierr)==3 )
		 {	wv->init=init,wv->ierr=ierr;
			wv->imin=imin,wv->imax=0.0 ;
			wv->flags |= (VAR_MIN|VAR_ERR);
		 }
		else if ( sscanf(varstr,"%lg[%lg:%lg]",&init,&imin,&imax)==3 )
		 {	wv->init=init,wv->ierr=1.0;
			if ( imin>imax ) w=imin,imin=imax,imax=w;
			wv->imin=imin,wv->imax=imax;
			wv->flags |= (VAR_MIN|VAR_MAX);
		 }
		else if ( sscanf(varstr,"%lg[:%lg]",&init,&imax)==2 )
		 {	wv->init=init,wv->ierr=1.0;
			wv->imin=0.0 ,wv->imax=imax;
			wv->flags |= (VAR_MAX);
		 }
		else if ( sscanf(varstr,"%lg[%lg:]",&init,&imin)==2 )
		 {	wv->init=init,wv->ierr=1.0;
			wv->imin=imin,wv->imax=0.0;
			wv->flags |= (VAR_MIN);
		 }
		else if ( sscanf(varstr,"%lg:%lg",&init,&ierr)==2 )
		 {	wv->init=init,wv->ierr=ierr;
			wv->flags |= (VAR_ERR);
		 }
		else if ( sscanf(varstr,"%lg",&init)==1 )
			wv->init=init,wv->ierr=1.0;
		else
			wv->init=0.0, wv->ierr=1.0;

		/*
		r1=sscanf(varstr,"%lg[%lg]",&wv->init,&wv->ierr);
		r2=sscanf(varstr,"%lg:%lg" ,&wv->init,&wv->ierr);
		r=(r1>r2?r1:r2);
		if ( r<1 )
		 {	wv->init=0.0;
			wv->ierr=1.0;
		 }
		else if ( r<2 )
		 	wv->ierr=1.0;
		*/ /* obsolete (simpler) format parser */
	 }
	else
	 {	wv->init=0.0;
		wv->ierr=0.0;		/* unused */
		wv->imin=wv->imax=0.0;	/* unused */
		wv->flags=0;
	 }
	while ( *varstr && *varstr != ',' )	varstr++;
	if ( *varstr==',' )			*varstr=0,varstr++;
	nvar++;
  }

 if ( rvars != NULL )	*rvars=vars;
 if ( rnvar != NULL )	*rnvar=nvar;

 return(0);
}

int extract_derived_variables(char *dvrstr,dvariable **rdvars,int *rndvar)
{
 dvariable	*dvars,*wv;
 int		ndvar,pl;
 char		*dvrstr0;

 dvars=NULL;
 ndvar=0;

 while ( *dvrstr )
  {	dvars=(dvariable *)realloc(dvars,sizeof(dvariable)*(ndvar+1));
	wv=&dvars[ndvar];
	wv->name=dvrstr;

	strcpy(wv->format,LFIT_DEFAULT_FORMAT);

	while ( *dvrstr && *dvrstr != '=' && *dvrstr != ',' )	dvrstr++;

	if ( *dvrstr != '=' )
	 {	free(dvars);
		dvars=NULL;
		return(1);
	 }
	*dvrstr=0,dvrstr++;
	dvrstr0=dvrstr;
	pl=0;
	while ( *dvrstr )
	 {	if ( *dvrstr=='(' )			pl++;
		else if ( *dvrstr==')' )		pl--;
		else if ( *dvrstr==',' && pl<=0 )	break;
		dvrstr++;
	 }
	wv->expr=dvrstr0;
	if ( *dvrstr==',' )	*dvrstr=0,dvrstr++;

	/* fprintf(stderr,"[%s]=[%s]\n",wv->name,wv->expr); */

	ndvar++;
  }

 if ( rdvars != NULL )	*rdvars=dvars;
 if ( rndvar != NULL )	*rndvar=ndvar;

 return(0);
}

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

int extract_variable_formats(char *formatstr,
	variable *vars,int nvar,dvariable *dvars,int ndvar)
{
 char	*fstr,**formats,*args[4],*p;
 int	i,k0,k1,n,is_failed;

 if ( formatstr==NULL )
	return(0);

 fstr=strdup(formatstr);
 remove_spaces(fstr);

 is_failed=0;
 formats=tokenize_char_dyn(fstr,',');
 for ( i=0 ; formats != NULL && formats[i] != NULL && ! is_failed ; i++ )
  {	for ( p=formats[i] ; *p ; p++ )
	 {	if ( *p==':' )	*p='=';		}
	n=tokenize_char(formats[i],args,'=',2);
	if ( n<2 )
	 {	is_failed=1;
		break;
	 }
	for ( k0=0 ; vars  != NULL && k0<nvar  ; k0++ )
	 {	if ( strcmp(vars [k0].name,args[0])==0 ) break;	}
	for ( k1=0 ; dvars != NULL && k1<ndvar && k0>=nvar ; k1++ )
	 {	if ( strcmp(dvars[k1].name,args[0])==0 ) break;	}
	if ( k0>=nvar && k1>=ndvar )
	 {	is_failed=1;
		break;
	 }

	p=args[1];
	while ( *p=='%' )	p++;
	if ( strlen(p)>14 || format_check(p) )
	 {	is_failed=1;
		break;
	 }	
	
	if ( k0<nvar )
	 {	vars[k0].format[0]='%';
		strcpy(vars[k0].format+1,p);
	 }
	else 
	 {	dvars[k1].format[0]='%';
		strcpy(dvars[k1].format+1,p);
	 }
  }

 if ( formats != NULL )
	free(formats);

 free(fstr);

 return(is_failed);
}

int extract_variable_differences(char *vardiffstr,variable *vars,int nvar)
{
 char	*vstr,**diffs,*args[4],*p;
 int	i,k,n,is_failed;
 double	diff;

 if ( vardiffstr==NULL )
	return(0);

 vstr=strdup(vardiffstr);
 remove_spaces(vstr);

 is_failed=0;
 diffs=tokenize_char_dyn(vstr,',');
 for ( i=0 ; diffs != NULL && diffs[i] != NULL && ! is_failed ; i++ )
  {	for ( p=diffs[i] ; *p ; p++ )
	 {	if ( *p==':' )	*p='=';		}
	n=tokenize_char(diffs[i],args,'=',2);
	if ( n<2 )
	 {	is_failed=1;
		break;
	 }
	for ( k=0 ; k<nvar ; k++ )
	 {	if ( strcmp(vars[k].name,args[0])==0 )	break;	}
	if ( k>=nvar )
	 {	is_failed=1;
		break;
	 }
	if ( sscanf(args[1],"%lg",&diff)<1 )
	 {	is_failed=1;
		break;
	 }
	vars[k].diff=diff;;
  }

 if ( diffs != NULL )
	free(diffs);
 free(vstr);

 return(is_failed);
}

int extract_dumpexpr_formats(char *formatstr,dumpexpr **rexprs,int *rnexpr)
{
 dumpexpr	*exprs;
 int		nexpr;
 char		*fstr,**formats,*p;
 int		i,is_failed;

 if ( formatstr==NULL )
  {	*rexprs=NULL;
	*rnexpr=0;
	return(0);
  }

 exprs=NULL;
 nexpr=0;

 fstr=strdup(formatstr);
 remove_spaces(fstr);

 is_failed=0;
 formats=tokenize_char_dyn(fstr,',');
 for ( i=0 ; formats != NULL && formats[i] != NULL && ! is_failed ; i++ )
  {	p=formats[i];
	while ( *p=='%' )	p++;
	if ( strlen(p)>14 || format_check(p) )
	 {	is_failed=1;
		break;
	 }
	exprs=realloc(exprs,sizeof(dumpexpr)*(nexpr+1));
	exprs[nexpr].format[0]='%';
	strcpy(exprs[nexpr].format+1,p);
	nexpr++;
  }

 if ( formats != NULL )	free(formats);
 free(fstr);

 if ( is_failed )
  {	if ( exprs != NULL )	free(exprs);
	exprs=NULL;
	nexpr=0; 
  }

 *rexprs=exprs;
 *rnexpr=nexpr;

 return(is_failed);
}

int extract_columns(char *colstr,column **rcols,int *rncol)
{
 column		*cols,*wc;
 int		ncol,cc,pc; /* num of columns, current col., previous col. */
 int		i;
 char		*name;

 cols=NULL;
 ncol=0;

 if ( rcols != NULL )	*rcols=cols;
 if ( rncol != NULL )	*rncol=ncol;

 pc=-1;
 while ( *colstr )
  {	name=colstr;
	while ( *colstr && *colstr != ',' && *colstr != ':' )	colstr++;

	if ( *colstr==':' )
	 {	*colstr=0,colstr++;
		sscanf(colstr,"%d",&cc);
		cc--;
	 }
	else if ( pc<0 )
		cc=0;
	else
		cc=pc+1;

	/* invalid */
	if ( cc<0 )
	 {	if ( cols != NULL )	free(cols);
		return(1);
	 }

	/* new */
	if ( cc>=ncol )
	 {	int	pcol;
		pcol=ncol;
		ncol=cc+1;
		cols=(column *)realloc(cols,sizeof(cols)*ncol);
		for ( i=pcol ; i<ncol ; i++ )
		 {	cols[i].name=NULL; 		}
	 }

	/* already defined */
	else if ( cols[cc].name != NULL )
	 {	free(cols);
		return(1);
	 }

	wc=&cols[cc];
	wc->name=name;

	while ( *colstr && *colstr != ',' )	colstr++;
	if ( *colstr==',' )	*colstr=0,colstr++;

	pc=cc;
  };

 for ( i=0 ; i<ncol ; i++ )
  {	if ( cols[i].name==NULL )
		cols[i].name="";
  }

 if ( rcols != NULL )	*rcols=cols;
 if ( rncol != NULL )	*rncol=ncol;

 return(0);
}

static int is_var_in_psn(psn *w,int nvar)
{
 int	j;
 for ( j=0 ; j<w->nterm ; j++ )
  {	if ( w->terms[j].type==T_VAR && w->terms[j].major < nvar )
		return(1);
  }
 return(0);
}

int build_psn_replace_macros(psn **rfunct)
{
 psnmacro	*wmacro;
 psn		*funct;

 if ( lpg.pl_macro==NULL )
	return(0);

 for ( wmacro=lpg.pl_macro ; wmacro->name != NULL ; wmacro++ )
  {	funct=psn_replace(*rfunct,wmacro->major,wmacro->pmacro,NULL);
	psn_free(*rfunct);
	*rfunct=funct;
  }

 return(0);
}

static int build_psn_search_symbol(psnsym **mysyms,char *name)
{
 psnsym	*ws;
 for ( ; *mysyms != NULL ; mysyms++ )
  {	for ( ws=*mysyms ; ws->type && ws->name != NULL ; ws++ )
	 {	if ( ws->type==T_VAR && strcmp(ws->name,name)==0 )
		 {	return(ws->major);		}
	 }
  }
 return(-1);
}

int build_psn_base_sequences(datablock *db,char *fncstr,
	lfitfunction *lfuncts,int nlfunct,
	psnsym **mysyms,int nvar)
{
 if ( fncstr==NULL || db==NULL || mysyms==NULL )
	return(-1);

 else if ( fncstr[0]=='@' )
  {	char		*lstr,*p,**arr;
	int		i,l,n,m;
	lfitfunction	*ll;

	fncstr++;
	lstr=strdup(fncstr);

	remove_spaces(lstr);	
	l=strlen(lstr);
	p=strchr(lstr,'(');
	if ( p==NULL || l<1 || lstr[l-1] != ')' )
	 {	free(lstr);
		return(7);
	 }
	lstr[l-1]=0;
	*p=0;
	p++;

	for ( i=0,ll=NULL ; i<nlfunct ; i++ )
	 {	if ( strcmp(lfuncts[i].name,lstr)==0 )
		 {	ll=&lfuncts[i];
			break;
		 }
	 }

	if ( ll==NULL )
	 {	free(lstr);
		return(7);
	 }
	
	arr=tokenize_char_dyn(p,',');
	for ( n=0 ; arr != NULL && arr[n] != NULL ; )	n++;

	if ( n != (ll->nvar + ll->nidv) )
	 {	if ( arr != NULL )	free(arr);
		free(lstr);
		return(7);
	 }

	db->lfunct =ll;
	db->map_var=(int *)malloc(sizeof(int)*ll->nvar);
	db->map_idv=(int *)malloc(sizeof(int)*ll->nidv);

	for ( i=0 ; i<n ; i++ )
	 {	m=build_psn_search_symbol(mysyms,arr[i]);
		if ( m<0 || (i<ll->nvar && m>=nvar) )
		 {	free(arr);	
			free(lstr);	
			return(7);
		 }
		if ( i<ll->nvar )
			db->map_var[i]=m;
		else
			db->map_idv[i-ll->nvar]=m-nvar;
	 }

	if ( arr != NULL )	free(arr);
	free(lstr);

	db->funct=NULL;
	db->functorig=NULL;
	db->diff=NULL;
  }
 else
  {	psn	*w,*pfunct;

	w=psn_conv_string(fncstr,mysyms);
	if ( w==NULL )			return(7);
	if ( psn_init(w,lpg.pl_prop) )	return(8);
	pfunct=psn_conv(w,lpg.pl_prop);
	if ( pfunct==NULL )		return(8);
	if ( psn_test(pfunct) )	return(8);
	psn_free(w);

	build_psn_replace_macros(&pfunct);

	db->funct=pfunct;
	psn_optimize(db->funct);
	db->fchain=psn_argument_chain(db->funct);
	lfit_psn_chain_optimize(db->funct,db->fchain);

	db->functorig=pfunct;

	db->lfunct =NULL;
	db->map_var=NULL;
	db->map_idv=NULL;

  }

 db->diff =NULL;
 db->dep  =NULL;
 db->err  =NULL;
 
 db->is_linear=0;

 return(0);
}

int lfitfunction_is_differentiable(lfitfunction *ll)
{
 if ( ll->flags==LFITFUNCTION_DIFFERENTIABLE )
	return(1);
 else if ( ll->flags==LFITFUNCTION_LINEAR )
	return(1);
 else
	return(0);
}
int lfitfunction_is_linear(lfitfunction *ll)
{
 if ( ll->flags==LFITFUNCTION_LINEAR )
	return(1);
 else
	return(0);
}


int build_psn_fit_sequences(datablock *db,char *fncstr,char *depstr,
	lfitfunction *lfuncts,int nlfunct,
	char *errstr,psnsym **mysyms,
	variable *vars,int nvar,int force_nonlin,int calc_derivatives)
{
 psn	*w,*pfunctorig,*pfunct,**pdiff,*pdep,*perr,*psub;
 int	is_linear,i,j,is_fitvar_fnc,is_fitvar_dep;

 if ( fncstr==NULL || db==NULL || mysyms==NULL )
	return(-1);

 /* be optimistic: */
 for ( i=0 ; i<nvar ; i++ )
  {	vars[i].is_linear=1;		}	
 /* (okay, in practice a variable is excluded from the linear ones if it 
	is found to be non-linear by some reason) */

 if ( fncstr[0]=='@' )
  {	char		*lstr,*p,**arr;
	int		i,l,n,m;
	lfitfunction	*ll;

	fncstr++;
	lstr=strdup(fncstr);

	remove_spaces(lstr);	
	l=strlen(lstr);
	p=strchr(lstr,'(');
	if ( p==NULL || l<1 || lstr[l-1] != ')' )
	 {	free(lstr);
		return(7);
	 }
	lstr[l-1]=0;
	*p=0;
	p++;

	for ( i=0,ll=NULL ; i<nlfunct ; i++ )
	 {	if ( strcmp(lfuncts[i].name,lstr)==0 )
		 {	ll=&lfuncts[i];
			break;
		 }
	 }

	if ( ll==NULL )
	 {	free(lstr);
		return(7);
	 }

	if ( calc_derivatives && (!lfitfunction_is_differentiable(ll)) )
		return(16);
	
	arr=tokenize_char_dyn(p,',');
	for ( n=0 ; arr != NULL && arr[n] != NULL ; )	n++;

	if ( n != (ll->nvar + ll->nidv) )
	 {	if ( arr != NULL )	free(arr);
		free(lstr);
		return(7);
	 }

	db->lfunct =ll;
	db->map_var=(int *)malloc(sizeof(int)*ll->nvar);
	db->map_idv=(int *)malloc(sizeof(int)*ll->nidv);

	db->is_linear=lfitfunction_is_linear(ll);

	for ( i=0 ; i<n ; i++ )
	 {	m=build_psn_search_symbol(mysyms,arr[i]);
		if ( m<0 || (i<ll->nvar && m>=nvar) )
		 {	free(arr);	
			free(lstr);	
			return(7);
		 }
		if ( i<ll->nvar )
		 {	db->map_var[i]=m;
			if ( ! db->is_linear )
				vars[m].is_linear=0; /* excluded */
		 }
		else
			db->map_idv[i-ll->nvar]=m-nvar;
	 }

	if ( arr != NULL )	free(arr);
	free(lstr);

	w=psn_conv_string(depstr,mysyms);
	if ( w==NULL ) 		return(12);
	if ( psn_init(w,lpg.pl_prop) )	return(13);
	pdep=psn_conv(w,lpg.pl_prop);
	if ( pdep==NULL )		return(13);
	if ( psn_test(pdep) )		return(13);
	psn_free(w);

	db->dep   =pdep;

	db->funct=NULL;
	db->functorig=NULL;
	db->diff=NULL;
  }

 else
  {	w=psn_conv_string(fncstr,mysyms);
	if ( w==NULL )			return(7);
	if ( psn_init(w,lpg.pl_prop) )	return(8);
	pfunct=psn_conv(w,lpg.pl_prop);
	if ( pfunct==NULL )		return(8);
	if ( psn_test(pfunct) )	return(8);
	psn_free(w);

	build_psn_replace_macros(&pfunct);

	w=psn_conv_string(depstr,mysyms);
	if ( w==NULL ) 		return(12);
	if ( psn_init(w,lpg.pl_prop) )	return(13);
	pdep=psn_conv(w,lpg.pl_prop);
	if ( pdep==NULL )		return(13);
	if ( psn_test(pdep) )		return(13);
	psn_free(w);

	build_psn_replace_macros(&pdep);

	is_fitvar_fnc=is_var_in_psn(pfunct,nvar);

	is_fitvar_dep=is_var_in_psn(pdep,nvar);
	if ( ! is_fitvar_fnc && ! is_fitvar_dep )	return(23);
	else if ( is_fitvar_fnc && is_fitvar_dep )	return(24);
	else if ( is_fitvar_dep && ! is_fitvar_fnc )
		w=pfunct,pfunct=pdep,pdep=w;

	pfunctorig=psn_duplicate(pfunct);

	is_linear=1;
	pdiff=(psn **)malloc(sizeof(psn *)*nvar);

	for ( i=0 ; i<nvar && calc_derivatives ; i++ )
	 {	
		pdiff[i]=psn_diff_simplify(pfunct,lpg.pl_prop,lpg.pl_diff,
			lpg.pl_simp,lpg.pl_funct,i);
		if ( pdiff[i]==NULL )	return(16);
	
		/*
		pdiff[i]=psn_diff(pfunct,lpg.pl_prop,lpg.pl_diff,i);
		if ( pdiff[i]==NULL )		return(16);
		if ( psn_test(pdiff[i]) )	return(16);
		j=psn_simplify(pdiff[i],plt_simp,lpg.pl_prop,lpg.pl_funct);
		if ( j || psn_test(pdiff[i]) )	return(16);
		*/

		if ( is_var_in_psn(pdiff[i],nvar) )
		 {	is_linear=0;
			if ( vars[i].is_linear )
				vars[i].is_linear=-1;
		 }
	 }

	if ( ! force_nonlin && is_linear )
	 {	psub=psn_duplicate(pfunct);
		for ( j=0 ; j<psub->nterm ; j++ )
		 {	if ( psub->terms[j].type==T_VAR && psub->terms[j].major < nvar )
			 {	psub->terms[j].type=T_SCONST;
				psub->terms[j].major=0;
			 }
		 }
		j=psn_simplify(psub,lpg.pl_simp,lpg.pl_prop,lpg.pl_funct);
		if ( psub->nterm>1 )
			is_linear=0;
		else if ( ! ( psub->terms[0].type==T_SCONST && psub->terms[0].major==0 ) )
			is_linear=0;
		if ( ! is_linear )
		 {	if ( is_verbose>=0 )
				fprint_warning("fit function is nonlinear but affine, it is going to be re-ordered");
			psn_concate(pfunct,psub);
			psn_append_term(pfunct,T_OP,O_SUB,2,-1);
			pfunct->nseq--;
			is_linear=1;
			for ( i=0 ; i<nvar ; i++ )
			 {	if ( vars[i].is_linear<0 )
					vars[i].is_linear=1;
			 }
		 }
		else
		 {	psn_free(psub);
			psub=NULL;
		 }
	  }
	else	psub=NULL;

	for ( i=0 ; i<nvar ; i++ )
	 {	if ( vars[i].is_linear<0 )
			vars[i].is_linear=0;
	 }

	if ( ! is_linear && ! force_nonlin )	return(14);
	if ( force_nonlin )			is_linear=0;

	if ( psub != NULL )
	 {	psn_concate(pdep,psub);
		psn_append_term(pdep,T_OP,O_SUB,2,-1);
		pdep->nseq--;
	 }

	db->funct =pfunct;
	psn_optimize(db->funct);
	db->fchain=psn_argument_chain(db->funct);
	lfit_psn_chain_optimize(db->funct,db->fchain); 
	db->diff  =pdiff;
	db->dep   =pdep;
	db->functorig=pfunctorig;
	db->is_linear=is_linear;
  }

 if ( errstr != NULL )
  {	w=psn_conv_string(errstr,mysyms);
	if ( w==NULL )	 		return(12);
	if ( psn_init(w,lpg.pl_prop) )	return(12);
	perr=psn_conv(w,lpg.pl_prop);
	if ( perr==NULL )		return(12);
	if ( psn_test(perr) )		return(12);
	psn_free(w);
	build_psn_replace_macros(&perr);
  }
 else
	perr=NULL;

 db->err   =perr;

 return(0);
}

int build_psn_constraint_sequences(lfitdata *lf,char *cntstr,psnsym **mysyms,
	variable *vars,int nvar)
{
 psn		*w,*pcnt,*pcdf;
 int		i,j;
 fitconstraint	*fc;
 domconstraint	*dc;

 /* process the constraint argument after -t|--constraint <constraints>.     */
 /* this <constraints> is a comma-separated list of <expression>=<constant>  */
 /* terms where <expression> is an arbitrary linear expression where all     */
 /* symbolic variables are the name of variables to be fitted. 		     */

 fc=&lf->fconstraint;
 dc=&lf->dconstraint;

 if ( cntstr != NULL )
  {	char	*c,*cstr,**ctokens;
	int	nn,nc,l,ctype; 
	double	**cmatrix,*cvector,*wvector;

	cmatrix=matrix_alloc(nvar);
	cvector=vector_alloc(nvar);
	wvector=vector_alloc(nvar);
	for ( j=0 ; j<nvar ; j++ )
	 {	wvector[j]=0.0;		}
	nc=0;

	dc->ndcfunct=0;
	dc->dcfuncts=NULL;

	cstr=strdup(cntstr);
	for ( c=cstr,l=0 ; *c ; c++ )
	 {	if ( *c=='(' || *c=='[' || *c=='{' )		 l++;
		else if ( *c==')' || *c==']' || *c=='}' )	 l--;
		else if ( *c==',' && l <= 0 )	*c=';';
	 }
	ctokens=tokenize_char_dyn(cstr,';');

	for ( nn=0 ; ctokens[nn] != NULL ; nn++ )
	 {	char	*p1,*p2,*cstr;

		remove_spaces(ctokens[nn]);

		p1=ctokens[nn];
		p2=NULL;
		ctype=0;

		for ( c=p1 ; *c ; c++ )
		 {	if ( *c=='<' && *(c+1)=='=' )
				ctype=-1,p2=c+2;
			else if ( *c=='=' && *(c+1)=='<' )
				ctype=-1,p2=c+2;
			else if ( *c=='>' && *(c+1)=='=' )
				ctype=+1,p2=c+2;
			else if ( *c=='=' && *(c+1)=='>' )
				ctype=+1,p2=c+2;
			else if ( *c=='=' && *(c+1)=='=' )
				ctype=0,p2=c+2;
		 	else if ( *c=='<' )
				ctype=-2,p2=c+1;
			else if ( *c=='>' )
				ctype=+2,p2=c+1;
			else if ( *c=='=' )
				ctype=0,p2=c+1;
			if ( p2 != NULL )
				break;
		 }
		cstr=NULL;
		if ( p2 != NULL )
		 {	*c=0;
			if ( ctype<0 )
			 {	strappendf(&cstr,"(%s)-(%s)",p2,p1);
				ctype=-ctype;
			 }
			else
				strappendf(&cstr,"(%s)-(%s)",p1,p2);
		 }
		else
			strappendf(&cstr,"(%s)",p1);

		/* fprintf(stderr,"cstr='%s'\n",cstr); */
		w=psn_conv_string(cstr,mysyms);

		free(cstr);	

		if ( w==NULL )			return(19);
		if ( psn_init(w,lpg.pl_prop) )	return(19);
		pcnt=psn_conv(w,lpg.pl_prop);
		if ( pcnt==NULL )		return(19);
		if ( psn_test(pcnt) )		return(19);
		psn_free(w);
		build_psn_replace_macros(&pcnt);

		/* check if there's any non-fitting variable in the expression */
		w=pcnt;
		for ( j=0 ; j<w->nterm ; j++ )
		 {	if ( w->terms[j].type==T_VAR && w->terms[j].major>=nvar )
				break;
		 }
		/* if yes, break. */
		if ( w->nterm && j<w->nterm )	return(19);

		/* domain constraint: */
		if ( ctype )
		 {	dcfunct	*dcf;
			dc->dcfuncts=(dcfunct *)realloc(dc->dcfuncts,sizeof(dcfunct)*(dc->ndcfunct+1));
			dcf=&dc->dcfuncts[dc->ndcfunct];

			dcf->funct=pcnt;
			dcf->ctype=ctype;	

			dc->ndcfunct++;
		 }
		/* fit constraint: too many: */
		else if ( nc>=nvar )
			return(21);
		/* fit constraint: */
		else
		 {	int	i,j;
			double	cleft;

			cleft=0.0;

			for ( i=0 ; i<nvar ; i++ )
			 {	/* calculate the derivative of the constraint by the i'th fitting variable */
				pcdf=psn_diff(pcnt,lpg.pl_prop,lpg.pl_diff,i);
				if ( pcdf==NULL )	return(19);
				if ( psn_test(pcdf) )		return(9);
				j=psn_simplify(pcdf,lpg.pl_simp,lpg.pl_prop,lpg.pl_funct);
				if ( j || psn_test(pcdf) )	return(9);
				w=pcdf;
				/* is there any variable in the derivative? */
				for ( j=0 ; j<w->nterm ; j++ )
				 {	if ( w->terms[j].type==T_VAR )	break;	}
				/* if yes, then the constraint is non-linear, break */
				if ( w->nterm && j<w->nterm )	return(19);
				/* calculate the derivative */
				lfit_psn_double_calc(pcdf,NULL,lpg.pl_funct,&cmatrix[nc][i],wvector);
				psn_free(pcdf);
			 }
			lfit_psn_double_calc(pcnt,NULL,lpg.pl_funct,&cleft,wvector);
			cvector[nc]=-cleft;

			psn_free(pcnt);

			nc++;
		 }	
	 }

	free(ctokens);
	free(cstr);

	vector_free(wvector);

	fc->nc=nc;
	fc->cmatrix   =cmatrix;
	fc->cvector   =cvector;

  }
 else
  {	fc->nc=0;
	fc->cmatrix   =NULL;
	fc->cvector   =NULL;

	dc->ndcfunct=0;
	dc->dcfuncts=NULL;
  }

 /* add extra constraints, if the given variable (to be fitted) was declared */
 /* as <variable>:=<initial value> instead of <variable>=<initial value>.    */
 for ( i=0 ; i<nvar ; i++ )
  {	if ( ! (vars[i].flags & VAR_IS_CONSTANT) )
		continue;
	if ( fc->nc >= nvar )
		return(21);
	if ( fc->cmatrix==NULL )	fc->cmatrix=matrix_alloc(nvar);
	if ( fc->cvector==NULL )	fc->cvector=vector_alloc(nvar);
	for ( j=0 ; j<nvar ; j++ )
	 {	fc->cmatrix[fc->nc][j]=(i==j?1.0:0.0);		}
	fc->cvector[fc->nc]=vars[i].init;
	fc->nc++;
  }

 if ( fc->nc > 0 )
  {	double	**invlmatrix;
	double	**cproject;

	invlmatrix=matrix_alloc(nvar+fc->nc);
	cproject  =matrix_alloc(nvar);
	
	if ( constraint_initialize_lambda_matrix(fc->nc,nvar,fc->cmatrix,invlmatrix) )
		return(19);

	constraint_initialize_proj_matrix(fc->nc,nvar,fc->cmatrix,cproject);

	fc->cproject  =cproject;
	fc->invlmatrix=invlmatrix;
  }
 else
  {	double	**invlmatrix;
	double	**cproject;

	cproject  =matrix_alloc(nvar);
	invlmatrix=matrix_alloc(nvar);
	for ( i=0 ; i<nvar ; i++ )
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	cproject  [i][j]=0.0;
			invlmatrix[i][j]=0.0;
		 }
		cproject  [i][i]=1.0;
		invlmatrix[i][i]=1.0;
	 }

	fc->invlmatrix=invlmatrix;
	fc->cproject  =cproject;
  }

 return(0);
}

int build_psn_derived_variable_expression(dvariable *wv,psnsym **mysyms,int nvar,int is_calc_derivatives)
{ 
 psn	*w,*pfunct;
 int	i;

 if ( wv->expr==NULL )
	return(-1);

 w=psn_conv_string(wv->expr,mysyms);
 if ( w==NULL )			return(1);
 if ( psn_init(w,lpg.pl_prop) )	return(2);
 pfunct=psn_conv(w,lpg.pl_prop);
 if ( pfunct==NULL )		return(2);
 if ( psn_test(pfunct) )	return(2);
 psn_free(w);

 build_psn_replace_macros(&pfunct);

 wv->funct=pfunct;

 if ( is_calc_derivatives )
  {	wv->diff=(psn **)malloc(sizeof(psn *)*nvar);
	for ( i=0 ; i<nvar ; i++ )
	 {	wv->diff[i]=psn_diff_simplify(wv->funct,lpg.pl_prop,lpg.pl_diff,lpg.pl_simp,lpg.pl_funct,i);
		if ( wv->diff[i]==NULL )
			return(3);
	 }
  }
 else
	wv->diff=NULL;
 
 return(0);
}


/*****************************************************************************/

#define		F_RANDOM	124		/* random(): uniform [0:1]   */
#define		F_GAUSSIAN	125		/* gaussian(): normal dist   */

static int funct_f_random(double *s)
{
 *s=drand48();
 return(0);
}

static int funct_f_gaussian(double *s)
{
 *s=get_gaussian(0.0,1.0);
 return(0);
}

int lfit_register_random_functions(void)
{
 int	r;

 r=0;
 r|=lfit_register_function("random"  , F_RANDOM  ,0,funct_f_random  ,NULL,"random()"  );
 r|=lfit_register_function("gaussian", F_GAUSSIAN,0,funct_f_gaussian,NULL,"gaussian()");

 return(0);
}

int lfit_define_single_macro(lfitdata *lf,int major,char *macro,char **ereason)
{
 int		i,argnum,nmacro;
 psnsym		*argsyms;
 char		*name,*firstarg,**args;
 psnsym		*mysyms[8];
 psn		*w,*pmacro;
 psnmacro	*wmacro;

 remove_spaces_and_comments(macro);

 name=macro;
 while ( *macro && *macro != '(' )	macro++;
 if ( ! *macro )
  {	if ( ereason != NULL )
		*ereason="empty macro definition(?)";
	return(1);
  }
 *macro=0,macro++;

 if ( lfit_register_symbol_exists(&lpg,name) )
  {	if ( ereason != NULL )
		*ereason="macro name is used by a built-in function or by a previously defined macro";
	return(2);
  }
 else if ( lfit_check_symbol_name(name) )
  {	if ( ereason != NULL )
		*ereason="macro name contains unexpected character(s)";
	return(2);
  }

 firstarg=macro;
 while ( *macro && *macro != ')' )	macro++;
 if ( ! *macro )
  {	if ( ereason != NULL )
		*ereason="unexpected format/syntax for macro definition";
	return(3);
  }
 *macro=0,macro++;
 if ( *macro != '=' )
  {	if ( ereason != NULL )
		*ereason="unexpected format/syntax for macro definition";
	return(4);
  }
 macro++;
 args=tokenize_char_dyn(firstarg,',');
 for ( argnum=0 ; args[argnum] != NULL ; )	argnum++;
 if ( args==NULL || argnum<1 )
  {	if ( args != NULL )	free(args);
	if ( ereason != NULL )
		*ereason="macro should have at least one argument";
	return(5);
  }

 argsyms=(psnsym *)malloc(sizeof(psnsym)*(argnum+1));
 for ( i=0 ; i<argnum ; i++ )
  {	argsyms[i].type =T_VAR;
	argsyms[i].major=i;
	argsyms[i].name =args[i];
	argsyms[i].minor=0;
  }	
 argsyms[argnum].type =0;
 argsyms[argnum].major=0;
 argsyms[argnum].name =NULL;
 argsyms[argnum].minor=0;

 mysyms[0]=argsyms;
 mysyms[1]=lpg.pl_sym;
 mysyms[2]=NULL;

 w=psn_conv_string(macro,mysyms);

 free(argsyms);
 free(args);
 if ( w==NULL )
  {	if ( ereason != NULL )
		*ereason="symbolic error";
	return(7);
  }

 if ( psn_init(w,lpg.pl_prop) )
  {	psn_free(w);
	if ( ereason != NULL )
		*ereason="syntax error";
	return(8);
  }
 pmacro=psn_conv(w,lpg.pl_prop);

 psn_free(w);
 if ( pmacro==NULL || psn_test(pmacro) )
  {	if ( ereason != NULL )
		*ereason="syntax error";
	return(8);
  }

 build_psn_replace_macros(&pmacro);

 if ( lpg.pl_macro==NULL )
 	nmacro=0;
 else
  {	for ( nmacro=0 ; lpg.pl_macro[nmacro].name != NULL ; )
		nmacro++;
  }

 lpg.pl_macro=(psnmacro *)realloc(lpg.pl_macro,sizeof(psnmacro)*(nmacro+2));
 wmacro=&lpg.pl_macro[nmacro];
 wmacro->name=strdup(name);
 wmacro->pmacro=pmacro;
 wmacro->major=major;
 wmacro->argnum=argnum;

 /* fprintf(stderr,"lfit_register_function(): name=%s major=%d argnum=%d\n",wmacro->name,wmacro->major,argnum); */
 if ( lfit_register_function(wmacro->name,wmacro->major,argnum,NULL,NULL,NULL) )
  {	if ( ereason != NULL )
		*ereason="registration failed";
	return(9);
  }

 nmacro++; 
 wmacro=&lpg.pl_macro[nmacro];
 wmacro->name=NULL;
 wmacro->pmacro=NULL;
 wmacro->major=0;
 wmacro->argnum=0;

 return(0);
}

int lfit_define_macro(lfitdata *lf,char *imacro,int *omajor)
{
 char	*macro,*c,**macrolist,**m,*reason;
 int	plevel,ret;

 if ( imacro==NULL )
	return(0);
 macro=strdup(imacro);

 for ( c=macro,plevel=0 ; *c ; c++ )
  {	if ( *c == '(' )	plevel++;
	else if ( *c == ')' )	plevel--;
	else if ( *c == ',' && plevel <= 0 )	*c=';';
  }

 macrolist=tokenize_char_dyn(macro,';');
 ret=0;
 for ( m=macrolist ; macrolist != NULL && *m != NULL ; m++ )
  {	ret=lfit_define_single_macro(lf,*omajor,*m,&reason);
	if ( ret )
	 {	fprint_warning("unable to register macro '%s': %s",*m,reason);		}
	(*omajor)++;
  }
 if ( macrolist != NULL )
	free(macrolist);

 free(macro);

 return(ret);
}

/*****************************************************************************/

int lfit_is_method_requries_derivatives(int fit_method)
{
 int	ret;

 switch ( fit_method ) 
  {  case FIT_METHOD_CLLS:
     case FIT_METHOD_NLLM:
     case FIT_METHOD_XMMC:
     case FIT_METHOD_FIMA:
	ret=1;
	break;
     case FIT_METHOD_MCMC:
     case FIT_METHOD_MCHI:
     case FIT_METHOD_EMCE:
     case FIT_METHOD_DHSX:
     case FIT_METHOD_LMND:
	ret=0;
	break;
     case FIT_METHOD_NONE:
	ret=0;
	break;
     default:	/* unexpected */
	ret=0;
	break;
  }

 return(ret);
}

int lfit_is_method_random(int fit_method)
{
 int	ret;

 switch ( fit_method ) 
  {  case FIT_METHOD_MCMC:
     case FIT_METHOD_EMCE:
     case FIT_METHOD_XMMC:
     case FIT_METHOD_FIMA:
	ret=1;
	break;
     default:
	ret=0;
	break;
  }

 return(ret);
}

int lfit_set_separated_linears(fitconstraint *fc,variable *vars,int nvar,char *varseplist)
{
 char	*list,**cmd,*p;
 int	c,i,n,set;

 if ( varseplist==NULL || (*varseplist)==0 )
	return(0);

 list=strdup(varseplist);
 cmd=tokenize_char_dyn(varseplist,',');
 if ( cmd==NULL )
  {	free(list);
	return(0);
  }

 for ( n=0 ; cmd[n] != NULL ; n++ )
  {	p=cmd[n];
	if ( strcmp(p,"@")==0 )
	 {	for ( i=0 ; i<nvar ; i++ )
		 {	if ( vars[i].is_linear )
				vars[i].is_separated=1;
		 }
	 }
	else 
	 {	if ( *p=='!' )
		 {	p++;
			set=0;
		 }
		else
			set=1;
		for ( i=0 ; i<nvar ; i++ )
		 {	if ( strcmp(vars[i].name,p)==0 )
				break;
		 }
		if ( i<nvar && vars[i].is_linear )
			vars[i].is_separated=set;
		else if ( i<nvar && is_verbose>=0 )
			fprint_warning("(one of the) fit function(s) is nonlinear in the parametric variable '%s', separation ignored",vars[i].name);
		else if ( i>=nvar )
			return(1);
	 }

  }
 free(cmd);
 free(list);

/* 
 This loop below guarantees that the set of parametric variables affacted 
  by explicit constraints and the set of variables separated (because of 
  their lineary) are completely disjoint. That would be a great fun both
  for its theoretical aspects and for implementation, to see what's going
  on if these sets are not disjoint. However, in practice this occurs
  barely, eventually, coupled constraints (this is the main point) 
  are used when the model function is linear in these parametric variables.
*/

 for ( c=0 ; c<fc->nc && fc->cmatrix != NULL ; c++ )
  {	for ( i=0 ; i<nvar ; i++ )
	 {	if ( fc->cmatrix[c][i] != 0.0 && vars[i].is_separated )
		 {	/* tbd: print some warning here? */
			vars[i].is_separated=0;
		 }
	 }
  }

 return(0);
}

/*****************************************************************************/

int f_lfitfunction_generic(double *stack)
{
 return(0);
}

/*****************************************************************************/

int main(int argc,char *argv[])
{
 FILE		*fr,*fw;
 char		*ofname,*oaname,*obname,*outname,
		*inname,*oxname,*ovname;
 int		i,errdump,resdump,is_help,is_warn_obsolete,is_dump_delta,seed;
 int		is_fit,force_nonlin,
		is_weighted_sigma,*colouts,ncolout;
 fitout		of;

 char		*fncstr,*depstr,*errstr,*colstr,
		*varstr,*dvrstr,
		*cntstr,*clostr,
		*vardiffstr,*formatstr,*corrfmstr,*varseplist,
		**xdefinelist,*perturbstr,*submethodstr;
 char		**inp_keys,**inp_args,
		**col_keys,**col_args,
		**fnc_keys,**fnc_args,
		**dep_keys,**dep_args,
		**err_keys,**err_args,
		**wgt_keys,**wgt_args;
 char		**keys;
 int		nkey;
 int		errtype;	/* errtype=0: sigma, errtype=1: weight */

 variable	*vars;		/* primary variables to fit / to analyze */
 int		nvar;
 dvariable	*dvars;		/* derived variables of the analysis */
 int		ndvar;

 dumpexpr	*exprs;
 int		nexpr;

 char		**dliblist;	/* user argument		*/
 lfitdynlib	*ldyns;		/* dynamically loaded stuff	*/
 int		nldyn;
 lfitfunction	*lfuncts;	/* gathered 'lfitfunction' data	*/
 int		nlfunct;
 int		omajor;

 psnsym		*varsym,*mysyms[8];

 lfitdata	lf_static,*lf=&lf_static;

 char		*errorstr;

 progbasename=strrchr(argv[0],'/');
 if ( progbasename != NULL )	progbasename++;
 else				progbasename=argv[0];

 for ( i=0 ; psnlfit_list_builtin_normal_operators[i].name != NULL ; i++ )
  {	int	r;
	r=lfit_register_internal(&psnlfit_list_builtin_normal_operators[i]);
	if ( r )
	 {	fprint_error("unable to register a built-in operator, exiting");
		return(1);	
 	 }
	
  }

 for ( i=0 ; psnlfit_list_builtin_elementary_functions[i].name != NULL ; i++ )
  {	int	r;
	r=lfit_register_internal(&psnlfit_list_builtin_elementary_functions[i]);
	if ( r )
	 {	fprint_error("unable to register a built-in function, exiting");
		return(1);	
 	 }
  }

 for ( i=0 ; psnlfit_list_builtin_interpolators[i].name != NULL ; i++ )
  {	int	r;
	r=lfit_register_internal(&psnlfit_list_builtin_interpolators[i]);
	if ( r )
	 {	fprint_error("unable to register a built-in interpolator function, exiting");
		return(1);	
 	 }
  }

 for ( i=0 ; psnlfit_list_builtin_aa_functions[i].name != NULL ; i++ )
  {	int	r;
	r=lfit_register_internal(&psnlfit_list_builtin_aa_functions[i]);
	if ( r )
	 {	fprint_error("unable to register an additional built-in analytic function, exiting");
		return(1);	
 	 }
  }

#if defined _FI_SOURCE || defined _ASTRO_EXTEN
 for ( i=0 ; psnlfit_list_builtin_xfuncts[i].name != NULL ; i++ )
  {	int	r;
	r=lfit_register_internal(&psnlfit_list_builtin_xfuncts[i]);
	if ( r )
	 {	fprint_error("unable to register an extra built-in function, exiting");
		return(1);	
 	 }
	
  }
#endif

 is_verbose=is_dump_delta=is_help=0;
 varstr=colstr=fncstr=depstr=errstr=cntstr=NULL;
 formatstr=NULL;
 corrfmstr=NULL;
 vardiffstr=NULL;
 varseplist=NULL;
 xdefinelist=NULL;

 lf->parameters.fit_method=FIT_METHOD_CLLS;	/* default: CLLS	     */
 submethodstr=NULL;				/* -P|--parameters: finetune */

 /* common for many MC method (MCMC, XMMC, EMCE): */
 lf->parameters.mc_iterations=1000;

 /* NLLM, LMND: */
 lf->parameters.tune.nllm.lambda=0.001,		/* lambda initial	     */
 lf->parameters.tune.nllm.lambda_mpy=10.0,	/* lambda multiplicator	     */
 lf->parameters.tune.nllm.max_iter=10;		/* NLLM interations	     */
 lf->parameters.tune.nllm.numeric_derivs=0;	/* NLLM or LMND?	     */

 /* DHSX: */
 lf->parameters.tune.dhsx.use_fisher_sx=0;	/* default: no need for      */
						/* partial derivatives but   */
						/* some initial uncertainties*/

 /* MCMC: */
 lf->parameters.tune.mcmc.use_gibbs=0;		/* default: normal sampler   */
 lf->parameters.tune.mcmc.count_accepted=1;	/* default: count accepted   */

 /* XMMC: */
 lf->parameters.tune.xmmc.skip_initial_fit=0;	/* skip initial XMMC/DHSX fit*/
 lf->parameters.tune.xmmc.is_adaptive=0;	/* default: non-adaptive XMMC*/
 lf->parameters.tune.xmmc.count_accepted=1;	/* default: count accepted   */
 lf->parameters.tune.xmmc.window=16;		/* for a nice fit, it's fine */
 lf->parameters.tune.xmmc.niter=0;		/* no additional iteratiosn  */

 /* EMCE: */
 lf->parameters.tune.emce.skip_initial_fit=0;	/* skip initial EMCE fit */
 lf->parameters.tune.emce.sub_method=FIT_METHOD_DHSX;

 /* FIMA: */
 lf->parameters.tune.fima.do_montecarlo=0;
 lf->parameters.tune.fima.write_original=1;
 lf->parameters.tune.fima.write_uncert=1;
 lf->parameters.tune.fima.write_correl=1;
 
 lf->parameters.niter=0; 		/* default: no n-sigma rejection     */
 lf->parameters.sigma=3.0;		/* default sigma level: 3	     */

 ofname=obname=oaname=oxname=ovname=outname=inname=NULL;

 errdump=0;
 resdump=-1;
 
 /* dynamlic extensions: */
 dliblist=NULL;

 /* function, dependent expr, errors, columns, adjusted variables...:*/
 fncstr=depstr=errstr=colstr=varstr=dvrstr=NULL;
 inp_keys=inp_args=NULL;
 col_keys=col_args=NULL;
 fnc_keys=fnc_args=NULL;
 dep_keys=dep_args=NULL;
 err_keys=err_args=NULL;
 wgt_keys=wgt_args=NULL;
 
 clostr=NULL;
 perturbstr=NULL;

 is_fit=0;
 is_weighted_sigma=0;
 errtype=0;

 seed=0;	/* random seed for MCMC/EMCE/XMMC and the random functions */
 is_warn_obsolete=0;

 errorstr=NULL;

 i=scanarg(argc,argv,SCANARG_ALLOW_FLAGS,	/* free letters: [gm] */
	/* general informative options: */
	"-h|--help|--short-help|--help-short:%SN1f",&is_help,
	"--example|--examples:%SN2f",&is_help,
	"--version:%SN-1f",&is_help,
	"--function-list|--functions|--list|--function|--list-functions:%SN3f",&is_help,
	"--long-help|--help-long:%SN4f",&is_help,

	/* options declaring outputs: */
	"-o|--output:%s",&outname,
	"-u|--output-fitted:%s",&ofname,
	"-j|--output-rejected:%s",&obname,
	"-a|--output-all:%s",&oaname,
	"-p|--output-expression:%s",&oxname,
	"-l|--output-variables|--output-vars:%s",&ovname,
	"--delta:%SN1f",&is_dump_delta,
	"--delta-comment:%SN-1f",&is_dump_delta,
/* ! */	"-of:%s%f",&ofname,&is_warn_obsolete,
/* ! */	"-or:%s%f",&obname,&is_warn_obsolete,
/* ! */	"-oa:%s%f",&oaname,&is_warn_obsolete,
/* ! */	"-ox:%s%f",&oxname,&is_warn_obsolete,
/* ! */	"-ov:%s%f",&ovname,&is_warn_obsolete,
/* ! */	"-om:%s%f",&outname,&is_warn_obsolete,

	/* fitting methods: */
	"-L|--clls|--linear:"	SNf(FIT_METHOD_CLLS),&lf->parameters.fit_method,
	"-N|--nllm|--nonlinear:"SNf(FIT_METHOD_NLLM),&lf->parameters.fit_method,
	"-M|--mcmc:"		SNf(FIT_METHOD_MCMC),&lf->parameters.fit_method,
	"-K|--mchi|--chi2:"	SNf(FIT_METHOD_MCHI),&lf->parameters.fit_method,
	"-E|--emce:"		SNf(FIT_METHOD_EMCE),&lf->parameters.fit_method,
	"-D|--dhsx|--downhill:"	SNf(FIT_METHOD_DHSX),&lf->parameters.fit_method,
	"-X|--xmmc:"	  	SNf(FIT_METHOD_XMMC),&lf->parameters.fit_method,
	"-U|--lmnd:"		SNf(FIT_METHOD_LMND),&lf->parameters.fit_method,
	"-A|--fima:"		SNf(FIT_METHOD_FIMA),&lf->parameters.fit_method,

	/* Fine tune: */
/* ! */	"-S|--emce-fit-method:%s%f",&submethodstr,&is_warn_obsolete,
	"-P|--param|--params|--parameters:%s",&submethodstr,
	"-s|--seed:%d",&seed,
	"--perturbations:%s",&perturbstr,
	"-i|--mcmc-iterations|--emce-iterations|--xmmc-iterations|--fima-iterations:%d",&lf->parameters.mc_iterations,
	"-k|--separate:%Cr",&varseplist,

	/* variables to fit or analyze: */
	"-v|--var|--vars|--variable|--variables:%Cr",&varstr,
	"-g|--derived-variable|--derived-variables:%Cr",&dvrstr,
	"-F|--format:%Cr",&formatstr,
	"-C|--correlation-format:%s",&corrfmstr,
	"-q|--difference|--differences:%Cr",&vardiffstr,

#ifdef	LFIT_ENABLE_DYNAMIC_EXTENSIONS
	"-d|--dynamic:%Dt",&dliblist,
#endif

	/* input: input content format, specifications of the columns: */
	"-c|--col|--cols|--column|--columns:%m",&colstr,
	"-c[0-9a-zA-Z]*:%Dl%Dt",&col_keys,&col_args,
	"-z|--columns-output|--column-output|--col-out:%s",&clostr,

	/* input: function to be fitted: */
	"-f|--funct|--function:%Cr",&fncstr,
	"-f[0-9a-zA-Z]*:%Dl%Dt",&fnc_keys,&fnc_args,

	/* input: dependent variable: */
	"-y|--dep|--dependent:%m",&depstr,
	"-y[0-9a-zA-Z]*:%Dl%Dt",&dep_keys,&dep_args,

	/* input: errors or weights to be used: */
	"-e|--error:%SN0f%m",&errtype,&errstr,
	"-e[0-9a-zA-Z]*:%Dl%Dt",&err_keys,&err_args,
	"-w|--weight:%SN1f%m",&errtype,&errstr,
	"-w[0-9a-zA-Z]*:%Dl%Dt",&wgt_keys,&wgt_args,

	/* linear constraints between the fitted variables: */
	"-t|--constraint|--constraints:%Cr",&cntstr,

	/* additional, user defined macros/functions: */
	"-x|--define|--macro:%Dt",&xdefinelist,

	/* fit error outputs: */
	"--err|--errors|--error-rows:%SN1f",&errdump,
	"--err-lin|--error-line:%SN2f",&errdump,
	"--err-col|--error-columns:%SN3f",&errdump,
	"--res|--residual:" 	SNf(RESIDUAL_BIASED)		,&resdump,
	"--unbiased-residual:"	SNf(RESIDUAL_UNBIASED)		,&resdump,
	"--residual-unbiased:"	SNf(RESIDUAL_UNBIASED)		,&resdump,
	"--residual-chi2:"	SNf(RESIDUAL_CHI2)		,&resdump,
	"--unbiased-chi2:"	SNf(RESIDUAL_UNBIASED_CHI2)	,&resdump,

	/* parameters of the sigma rejection: */
	"-n|--iterations:%d",&lf->parameters.niter,
	"-r|--sigma|--rejection|--rejection-level:%g",&lf->parameters.sigma,
	"--weighted-sigma:%f",&is_weighted_sigma,

	/* other remaining options, error handling, input file, ...: */
	"--verbose:%i",&is_verbose,"(V):%i",&is_verbose,
	"--quiet:%SN-1f",&is_verbose,
	"-i[0-9a-zA-Z]*:%Dl%Dt",&inp_keys,&inp_args,
	"-:%w",&inname,
	"+*|-*:%e",
	"*:%w",&inname,
	NULL);

 if ( i )
  {	fprint_error("invalid command line argument near '%s'",argv[i]);
	return(1);
  }
 else if ( is_help==1 )
  {	fprint_lfit_usage(stderr);
	return(0);
  }
 else if ( is_help==2 )
  {	fprint_lfit_examples(stdout);
	return(0);
  }
 else if ( is_help==3 )
  {	fprint_lfit_function_list(stdout,
		(psnlfit *)psnlfit_list_builtin_normal_operators,
		(psnlfit *)psnlfit_list_builtin_elementary_functions,
		(psnlfit *)psnlfit_list_builtin_interpolators,
		(psnlfit *)psnlfit_list_builtin_aa_functions,
#if defined _FI_SOURCE || defined _ASTRO_EXTEN
		(psnlfit *)psnlfit_list_builtin_xfuncts,
#endif
		NULL);
	return(0);
  }
 else if ( is_help==4 )
  {	fprint_lfit_long_help(stdout);
	return(0);
  }
 else if ( is_help<0 )
  {	fprintf(stdout,"lfit %s (%s)\n",LFIT_VERSION,LFIT_LASTCHANGE);
	fprintf(stdout,"Copyright (C) 1996, 2002, 2004-2008, 2009; Pal, Andras <apal@szofi.net>\n");
	return(0);
  }

 if ( is_warn_obsolete )
  {	fprintf(stderr,"Warning: you have specified some of the following command line arguments:\n");
	fprintf(stderr,"\t-of\t-ox\t-ov\t-or\t-oa\t-om\t-S\n");
	fprintf(stderr,"Note that these arguments are obsolete and will have different effects in the\n");
 	fprintf(stderr,"further releases of `lfit`. Please use these respective arguments instead:\n");
	fprintf(stderr,"\t-u\t-p\t-l\t-j\t-a\t-o\t-P\n");
  }

 if ( lf->parameters.fit_method==FIT_METHOD_EMCE && submethodstr != NULL )
  {	int	is_default;
	is_default=0;	/* de facto dummy variable */
	i=scanpar(submethodstr,SCANPAR_DEFAULT,
		"default|defaults:%f",&is_default,
		"clls|linear:"		SNf(FIT_METHOD_CLLS), &lf->parameters.tune.emce.sub_method,
		"nllm|nonlinear:"	SNf(FIT_METHOD_NLLM), &lf->parameters.tune.emce.sub_method,
		"lmnd:"			SNf(FIT_METHOD_LMND), &lf->parameters.tune.emce.sub_method,
		"dhsx|downhill:"	SNf(FIT_METHOD_DHSX), &lf->parameters.tune.emce.sub_method,
		"mcmc|mc|montecarlo:"	SNf(FIT_METHOD_MCMC), &lf->parameters.tune.emce.sub_method,
		"fisher:%f",		&lf->parameters.tune.dhsx.use_fisher_sx,
		"skip:%f",&lf->parameters.tune.emce.skip_initial_fit,
		NULL);
	if ( i )
	 {	fprint_error("invalid parameter argument in '%s'",submethodstr);
		return(1);
	 }
  }
 else if ( (lf->parameters.fit_method==FIT_METHOD_NLLM || 
	    lf->parameters.fit_method==FIT_METHOD_LMND) && submethodstr != NULL )
  {	int	is_default;
	is_default=0;	/* de facto dummy variable */
 	i=scanpar(submethodstr,SCANPAR_DEFAULT,
		"default|defaults:%f",&is_default,
		"lambda:%g"	,&lf->parameters.tune.nllm.lambda,
		"multiply:%g"	,&lf->parameters.tune.nllm.lambda_mpy,
		"iterations:%d"	,&lf->parameters.tune.nllm.max_iter,
		NULL);
	if ( lf->parameters.fit_method==FIT_METHOD_LMND )
		lf->parameters.tune.nllm.numeric_derivs=1;
	else
		lf->parameters.tune.nllm.numeric_derivs=0;
	if ( i )
	 {	fprint_error("invalid parameter argument in '%s'",submethodstr);
		return(1);
	 }
  }
 else if ( lf->parameters.fit_method==FIT_METHOD_MCMC && submethodstr != NULL )
  {	int	is_default;
	is_default=0;	/* de facto dummy variable */
	i=scanpar(submethodstr,SCANPAR_DEFAULT,
		"default|defaults:%f",&is_default,
		"accepted:%f",	&lf->parameters.tune.mcmc.count_accepted,
		"nonaccepted:%SN0f",&lf->parameters.tune.mcmc.count_accepted,
		"gibbs:%f",&lf->parameters.tune.mcmc.use_gibbs,
		NULL);
	if ( i )
	 {	fprint_error("invalid parameter argument in '%s'",submethodstr);
		return(1);
	 }
  }
 else if ( lf->parameters.fit_method==FIT_METHOD_DHSX && submethodstr != NULL )
  {	int	is_default;
	is_default=0;	/* de facto dummy variable */
	i=scanpar(submethodstr,SCANPAR_DEFAULT,
		"default|defaults:%f",&is_default,
		"fisher:%f",&lf->parameters.tune.dhsx.use_fisher_sx,
		NULL);
	if ( i )
	 {	fprint_error("invalid parameter argument in '%s'",submethodstr);
		return(1);
	 }
  }
 else if ( lf->parameters.fit_method==FIT_METHOD_XMMC && submethodstr != NULL )
  {	int	is_default;
	is_default=0;	/* de facto dummy variable */
	i=scanpar(submethodstr,SCANPAR_DEFAULT,
		"default|defaults:%f",&is_default,
		"accepted:%f",	&lf->parameters.tune.xmmc.count_accepted,
		"nonaccepted:%SN0f",&lf->parameters.tune.xmmc.count_accepted,
		"skip:%f",	&lf->parameters.tune.xmmc.skip_initial_fit,
		"adaptive:%f",	&lf->parameters.tune.xmmc.is_adaptive,
		"window:%d",	&lf->parameters.tune.xmmc.window,
		"iterations:%d",&lf->parameters.tune.xmmc.niter,
		NULL);
	if ( i )
	 {	fprint_error("invalid parameter argument in '%s'",submethodstr);
		return(1);
	 }
	if ( lf->parameters.tune.xmmc.is_adaptive && lf->parameters.tune.xmmc.niter>0 )
	 {	fprint_error("sorry, simultaneous adaptive and iterative XMMC is incompatible");
		return(1);
	 }
  }
 else if ( lf->parameters.fit_method==FIT_METHOD_FIMA && submethodstr != NULL )
  {	int	is_default,wo,wu,wc;
	is_default=0;	/* de facto dummy variable */
	wo=wu=wc=-1;
	i=scanpar(submethodstr,SCANPAR_DEFAULT,
		"default|defaults:%f",&is_default,
		"mc|montecarlo:%f",&lf->parameters.tune.fima.do_montecarlo,
		"orig|original:%SN1f",&wo,
		"noorig|nooriginal:%SN0f",&wo,
		"error|errors|uncert|uncertainty|uncertainties:%SN1f",&wu,
		"noerror|noerrors|nouncert|nouncertainty|nouncertainties:%SN0f",&wu,
		"corr|correl|correlation|correlations:%SN1f",&wc,
		"nocorr|nocorrel|nocorrelation|nocorrelations:%SN0f",&wc,
		NULL);
	if ( i )
	 {	fprint_error("invalid parameter argument in '%s'",submethodstr);
		return(1);
	 }
	if ( ! wo )	 lf->parameters.tune.fima.write_original=0;
	else if ( wo>0 ) lf->parameters.tune.fima.write_original=1;
	if ( ! wu )	 lf->parameters.tune.fima.write_uncert=0;
	else if ( wu>0 ) lf->parameters.tune.fima.write_uncert=1;
	if ( ! wc )	 lf->parameters.tune.fima.write_correl=0;
	else if ( wc>0 ) lf->parameters.tune.fima.write_correl=1;
  }


 if ( ! lf->parameters.fit_method )
	force_nonlin=0;
 else if ( lf->parameters.fit_method==FIT_METHOD_CLLS )
	force_nonlin=0;
 else if ( lf->parameters.fit_method==FIT_METHOD_EMCE && lf->parameters.tune.emce.sub_method==FIT_METHOD_CLLS )
	force_nonlin=0; 
 else
	force_nonlin=1;

 ldyns=NULL;
 nldyn=0;
 lfuncts=NULL;
 nlfunct=0;
#ifdef	LFIT_ENABLE_DYNAMIC_EXTENSIONS
 for ( i=0 ; dliblist != NULL && dliblist[i] != NULL ; i++ )
  {	void		*handle;
	char		*dlibspec,*p,*n;
	lfitfunction	*lff;
	int		j;

	dlibspec=strdup(dliblist[i]);

	p=strchr(dlibspec,':');
	if ( p==NULL )
	 {	fprint_warning("no function array specifications after '%s', skipped",dlibspec);
		free(dlibspec);
		continue;
	 }
	*p=0;
	p++;

	handle=dlopen(dlibspec,RTLD_NOW);
	if ( handle==NULL )
	 {	fprint_warning("unable to load shared library '%s' (%s), skipped",dlibspec,dlerror());
		free(dlibspec);
		continue;
	 }
	ldyns=(lfitdynlib *)realloc(ldyns,sizeof(lfitdynlib)*(nldyn+1));
	ldyns[nldyn].file=dliblist[i];
	ldyns[nldyn].handle=handle;
	nldyn++;

	while ( p != NULL && (*p) )
	 {	while ( isspace((int)(*p)) )	p++;

		n=strchr(p,',');
		if ( n != NULL )
		 {	*n=0;
			n++;
		 }
		
		lff=(lfitfunction *)dlsym(handle,p);
		if ( lff==NULL )
		 {	fprint_warning("function array '%s' cannot be found in '%s', skipped",p,dlibspec);
			j=0;
		 }
		else
		 {	for ( j=0 ; lff[j].name != NULL && lff[j].function != NULL ; )
				j++;
		 }
		if ( j>0 )
		 {	lfuncts=(lfitfunction *)realloc(lfuncts,sizeof(lfitfunction)*(nlfunct+j));
			memcpy(lfuncts+nlfunct,lff,sizeof(lfitfunction)*j);
			nlfunct+=j;
		 }
		p=n;
	 };

	free(dlibspec);
  }
#endif
 lpg.lffregs=NULL;
 lpg.nlffreg=0;

 omajor=LFIT_F_MAX_BUILTIN;

 if ( lfuncts != NULL && 0<nlfunct )
  {
 	for ( i=0 ; i<nlfunct ; i++ )
	 {	lfitfunction	*lff;
		lffreg		*lr;
		short		*diffrule;
		lff=&lfuncts[i];
		if ( lfit_check_symbol_name(lff->name) )
		 {	fprint_error("function name '%s' (from a dynamically loaded library) contains unexpected character(s)",lff->name);
			return(1);
		 }
		else if ( lfit_register_symbol_exists(&lpg,lff->name) )
		 {	fprint_error("function name '%s' (from a dynamically loaded library) has already been in use",lff->name);
			return(1);
		 }
		if ( lff->flags & LFITFUNCTION_DIFFERENTIABLE )
		 {	int	i,j,p,n;
			p=0;
			n=lff->nvar+lff->nidv;
			diffrule=(short *)malloc(sizeof(short)*lff->nvar*(n+4));
			for ( i=0 ; i<lff->nvar ; i++ )
			 {	for ( j=0 ; j<n ; j++ )
				 {	diffrule[p]=-SS_REF-n+j;
					p++;
				 }
				diffrule[p+0]=omajor+i+1;
				diffrule[p+1]=-n+i;
				diffrule[p+2]=O_MUL;
				p+=3;
				if ( i>0 )
					diffrule[p]=O_ADD,p++;
			 }
			diffrule[p]=0;
		 }
		else
			diffrule=NULL;

		lfit_register_function(lff->name,omajor,lff->nvar+lff->nidv,f_lfitfunction_generic,diffrule,NULL);
		lpg.lffregs=(lffreg *)realloc(lpg.lffregs,sizeof(lffreg)*(lpg.nlffreg+1));
		lr=&lpg.lffregs[lpg.nlffreg];

		lr->lff=lff;
		lr->fmajor=omajor;
		
		omajor++;
		if ( lff->flags & LFITFUNCTION_DIFFERENTIABLE )
		 {	int	i;
			for ( i=0 ; i<lff->nvar ; i++ )
			 {	lfit_register_function(NULL,omajor+i,lff->nvar+lff->nidv,f_lfitfunction_generic,NULL,NULL);	}
			omajor += lff->nvar;
		 }

		lr->nmajor=omajor;
		lr->diff=NULL;

		lpg.nlffreg++;	
	 }
  }

 if ( xdefinelist != NULL )
  {	char	**xdef;
	for ( xdef=xdefinelist ; *xdef ; xdef++ )
	 {	lfit_define_macro(lf,*xdef,&omajor);		}
  }

 if ( varstr != NULL )
  {	remove_spaces_and_comments(varstr);
	extract_variables(varstr,&vars,&nvar);
  }
 else
  {	vars=NULL,
	nvar=0;
  }

 for ( i=0 ; vars != NULL && i<nvar ; i++ )
  {	int	j;
	if ( lfit_register_symbol_exists(&lpg,vars[i].name) )
	 {	fprint_error("variable name '%s' has already been reserved (for a built-in function or macro)",vars[i].name);
		return(1);
	 }
	if ( lfit_check_symbol_name(vars[i].name) )
	 {	fprint_error("variable name '%s' contains unexpected character(s)",vars[i].name);
		return(1);
	 }
	for ( j=0 ; j<i ; j++ )
	 {	if ( strcmp(vars[i].name,vars[j].name)==0 )
		 {	fprint_error("variable '%s' has been declared more than once",vars[i].name);
			return(1);
		 }
	 }
  }
 for ( i=0 ; vars != NULL && i<nvar ; i++ )
  {	lfit_register_symbol_add(&lpg,vars[i].name);		}

 if ( dvrstr != NULL )
  {	remove_spaces_and_comments(dvrstr);
	extract_derived_variables(dvrstr,&dvars,&ndvar);
  }
 else
  {	dvars=NULL;
	ndvar=0;
  }
 for ( i=0 ; dvars != NULL && i<ndvar ; i++ )
  {	int	j;
	if ( lfit_register_symbol_exists(&lpg,dvars[i].name) )
	 {	fprint_error("name of the derived variable '%s' has already been reserved (for a built-in function or macro) or used as a primary fit variable",dvars[i].name);
		return(1);
	 }
	if ( lfit_check_symbol_name(dvars[i].name) )
	 {	fprint_error("derived variable name '%s' contains unexpected character(s)",dvars[i].name);
		return(1);
	 }
	for ( j=0 ; j<i ; j++ )
	 {	if ( strcmp(dvars[i].name,dvars[j].name)==0 )
		 {	fprint_error("derived variable '%s' has been declared more than once",dvars[i].name);
			return(1);
		 }
	 }
  }
 for ( i=0 ; dvars != NULL && i<ndvar ; i++ )
  {	lfit_register_symbol_add(&lpg,dvars[i].name);		}

 keys=NULL;
 nkey=0;
 if ( inp_keys != NULL )	key_add_list(&keys,&nkey,inp_keys,2);
 if ( col_keys != NULL )	key_add_list(&keys,&nkey,col_keys,2);
 if ( fnc_keys != NULL )	key_add_list(&keys,&nkey,fnc_keys,2);
 if ( dep_keys != NULL )	key_add_list(&keys,&nkey,dep_keys,2);
 if ( err_keys != NULL )	key_add_list(&keys,&nkey,err_keys,2);
 if ( wgt_keys != NULL )	key_add_list(&keys,&nkey,wgt_keys,2);
 if ( nkey>0 && ( colstr != NULL || fncstr != NULL || depstr != NULL || errstr != NULL ) )
  {	fprint_error("invalid combination of input fit data");
	return(1);
  }

 if ( nkey>0 )
  {	int		i,k,l,r,is_l_fit;
	char		*key,*arg,*arg2,*c;
	datablock	*db;

	lf->ndatablock=nkey;
	lf->datablocks=(datablock *)malloc(sizeof(datablock)*lf->ndatablock);

	is_fit=0;
	for ( i=0 ; i<nkey ; i++ )
	 {	key=keys[i];
		db=&lf->datablocks[i];

		db->key=key;

		arg=key_search_argument(inp_keys,2,key,inp_args);

		/* no input argument found: */
		if ( arg==NULL )
		 {	fprint_error("input file name is missing, use -i<key> 'filename' to define");
			return(1);
		 }
		db->inparg=arg;

		/* no column definition found: */
		arg=key_search_argument(col_keys,2,key,col_args);
		if ( arg==NULL )
		 {	fprint_error("column definitions missing, use -c<key> 'cols' to define");
			return(1);
		 }
		db->colarg=strdup(arg);
		if ( extract_columns(db->colarg,&db->cols,&db->ncol) )
		 {	fprint_error("invalid column specification in '%s' (note that each column can only be defined once)",arg);
			return(1);
		 }
		for ( k=0,r=0 ; k<db->ncol && (! r) ; k++ )
		 {	if ( db->cols[k].name==NULL || strlen(db->cols[k].name)<=0 )
				continue;
			if ( lfit_check_symbol_name(db->cols[k].name) )
			 {	r=db->ncol+k+1;
				break;
			 }
			else if ( lfit_register_symbol_exists(&lpg,db->cols[k].name) )
			 {	r=k+1;
				break;
			 }
			for ( l=0 ; l<k && (! r); l++ )
			 {	if ( db->cols[l].name==NULL || strlen(db->cols[l].name)<=0 )
					continue;
				if ( strcmp(db->cols[l].name,db->cols[k].name)==0 )
				 {	r=-k-1;
					break;
				 }
			 }
		 }
		if ( r>db->ncol )
		 {	fprint_error("column name '%s' contains unexpected character(s)",db->cols[r-db->ncol-1].name);
			return(1);
		 }
		else if ( r>0 )
		 {	fprint_error("column name '%s' has already been used as a function, macro or variable name",db->cols[r-1].name);
			return(1);
		 }
		else if ( r<0 )
		 {	fprint_error("column name '%s' has already been used as a name of another column",db->cols[(-r)-1].name);
			return(1);
		 }
	
		arg=key_search_argument(fnc_keys,2,key,fnc_args);
		if ( arg==NULL )
		 {	fprint_error("fit function is missing, use -f<key> 'func' to define one");
			return(1);
		 }
		db->fncarg=strdup(arg);
		arg=key_search_argument(dep_keys,2,key,dep_args);
		if ( arg != NULL )
			db->deparg=strdup(arg);
		else
			db->deparg=NULL;

		for ( c=db->fncarg ; *c ; c++ )
		 {	if ( *c=='=' )	break;		}
		if ( *c=='=' && db->deparg != NULL )
		 {	fprint_error("dependent value missing or ambigous, use -y<key> 'expr' or -f<key> '...=expr' to define");
			return(1);
		 }
		else if ( *c=='=' )
		 {	*c=0,c++;
			db->deparg=strdup(c);
		 }
		if ( db->deparg != NULL )	is_l_fit=1,is_fit++;
		else				is_l_fit=0;
		if ( ! nvar && is_l_fit )
		 {	fprint_error("definition of fit variables missing, use -v 'vars' to define");
			return(1);
		 }

		arg =key_search_argument(err_keys,2,key,err_args);
		arg2=key_search_argument(wgt_keys,2,key,wgt_args);
		if ( arg != NULL && arg2 != NULL )
		 {	fprint_error("both error and weight have been defined, use only one of them");
			return(1);
		 }
		else if ( arg != NULL )
			db->errarg=strdup(arg),
			db->errtype=0;
		else if ( arg2 != NULL )
			db->errarg=strdup(arg2),
			db->errtype=!0;
		else
			db->errarg=NULL,
			db->errtype=0;
	 }
	if ( is_fit && is_fit < nkey )
	 {	fprint_error("input file name is missing, use -i<key> 'filename' to define");
		return(1);
	 }
  }
 else
  {	datablock	*db;
	char		*c;
	int		k,l,r;

	lf->ndatablock=1;
	lf->datablocks=(datablock *)malloc(sizeof(datablock));

	db=&lf->datablocks[0];

	if ( colstr==NULL )
	 {	fprint_error("column definitions missing, use -c 'cols' to define");
		return(1);
	 }
	if ( extract_columns(colstr,&db->cols,&db->ncol) )
	 {	fprint_error("invalid column specification in '%s' (note that each column can only be defined once)",colstr);
		return(1);
	 }

	for ( k=0,r=0 ; k<db->ncol && (! r) ; k++ )
	 {	if ( db->cols[k].name==NULL || strlen(db->cols[k].name)<=0 )
			continue;
		if ( lfit_check_symbol_name(db->cols[k].name) )
		 {	r=db->ncol+k+1;
			break;
		 }
		else if ( lfit_register_symbol_exists(&lpg,db->cols[k].name) )
		 {	r=k+1;
			break;
		 }
		for ( l=0 ; l<k && (! r); l++ )
		 {	if ( db->cols[l].name==NULL || strlen(db->cols[l].name)<=0 )
				continue;
			if ( strcmp(db->cols[l].name,db->cols[k].name)==0 )
			 {	r=-k-1;
				break;
			 }
		 }
	 }
	if ( r>db->ncol )
	 {	fprint_error("column name '%s' contains unexpected character(s)",db->cols[r-db->ncol-1].name);
		return(1);
	 }
	if ( r>0 )
	 {	fprint_error("column name '%s' has already been used as a function, macro or variable name",db->cols[r-1].name);
		return(1);
	 }
	else if ( r<0 )
	 {	fprint_error("column name '%s' has already been used as a name of another column",db->cols[(-r)-1].name);
		return(1);
	 }

	db->colarg=colstr;

	if ( fncstr==NULL )
	 {	fprint_error("fit function is missing, use -f 'func' to define one");
		return(1);
	 }

	for ( c=fncstr ; *c ; c++ )
	 {	if ( *c=='=' )	break;		}
	if ( *c=='=' && depstr != NULL )
	 {	fprint_error("dependent value missing or ambigous, use -y 'expr' or -f '...=expr' to define");
		return(1);
	 }
	else if ( *c=='=' )
	 {	*c=0,c++;
		depstr=strdup(c);
	 }

	if ( depstr != NULL )	is_fit=1;
	else			is_fit=0;

	if ( ! nvar && is_fit )
	 {	fprint_error("definition of fit variables missing, use -v 'vars' to define");
		return(1);
	 }

	db->fncarg=fncstr;
	db->deparg=depstr;
	db->errtype=errtype;
	db->errarg=errstr;

	db->inparg=inname;

	db->key="default";
  }


 /* random() and gaussian() functions are only available in evaluation mode: */
 if ( ! is_fit )
  {	int	r;
	r=lfit_register_random_functions();
	if ( r )
	 {	fprint_error("unable to register a random number generator function, exiting");
		return(1);	
 	 }
  }
 
 if ( formatstr != NULL && is_fit )
  {	i=extract_variable_formats(formatstr,vars,nvar,dvars,ndvar);
	if ( i )	
	 {	fprint_error("invalid variable format string");
		return(1);
	 }
  }
 if ( corrfmstr != NULL )
  {	if ( corrfmstr[0]=='%' )	corrfmstr++;
	if ( strlen(corrfmstr)>14 || format_check(corrfmstr) )
	 {	fprint_error("invalid variable format string");
		return(1);
	 }
	lf->corrfm[0]='%';
	strncpy(lf->corrfm+1,corrfmstr,15);
  }
 else
	strcpy(lf->corrfm,LFIT_DEFAULT_CORR_FORMAT);

 if ( vardiffstr != NULL )
  {	i=extract_variable_differences(vardiffstr,vars,nvar);
	if ( i )	
	 {	fprint_error("invalid parameter for variable differences (required by numeric derivative estimations0");
		return(1);
	 }
	for ( i=0 ; i<nvar ; i++ )
	 {	if ( vars[i].diff <= 0.0 )
		 {	fprint_error("negative, zero difference or some of the differences has not been defined");
			return(1);
		 }
	 }
  }

 if ( formatstr != NULL && ( ! is_fit ) )
  {	i=extract_dumpexpr_formats(formatstr,&exprs,&nexpr);
	if ( i )	
	 {	fprint_error("invalid variable format string");
		return(1);
	 }
  } 
 else
  {	exprs=NULL;
	nexpr=0;
  }

 lf->maxncol=0;
 for ( i=0 ; i<lf->ndatablock ; i++ )
  {	if ( lf->maxncol < lf->datablocks[i].ncol )
		lf->maxncol=lf->datablocks[i].ncol;
  }


 if ( outname != NULL )
  {	fw=fopenwrite(outname);
	if ( fw==NULL )	
	 {	fprint_error("unable to create output file '%s'",outname);
		return(1);
	 }
  }
 else			fw=NULL;

 varsym=(psnsym *)malloc(sizeof(psnsym)*(nvar+1));
 for ( i=0 ; i<nvar ; i++ )
  {	varsym[i].type=T_VAR;
	varsym[i].major=i;
	varsym[i].name=vars[i].name;
  }
 varsym[nvar].type=0;
 varsym[nvar].major=0;
 varsym[nvar].name=NULL;
 lf->varsym=varsym;

 for ( i=0 ; dvars != NULL && i<ndvar ; i++ )
  {	dvariable	*wv;
	int		r;

	mysyms[0]=varsym,
	mysyms[1]=lpg.pl_sym,
	mysyms[2]=NULL;

	wv=&dvars[i];

	r=build_psn_derived_variable_expression(wv,mysyms,nvar,1);
	if ( r==1 )
	 {	fprint_error("symbolic error in the definition function of the derived variable '%s'",wv->name);
		return(1);
	 }
	else if ( r==2 )
	 {	fprint_error("syntax error in the definition function of the derived variable '%s'",wv->name);
		return(1);
	 }
	else if ( r>=3 )
	 {	fprint_error("non-differentiable expression in the definition function of the derived variable '%s'",wv->name);
		return(1);
	 }
	else if ( r )
 	 {	fprint_error("unknown error in the definition function of the derived variable '%s'",wv->name);
		return(1);
	 }
  }

 for ( i=0 ; i<lf->ndatablock ; i++ )
  {	datablock	*db;
	int		j;
	psnsym		*colsym;

	db=&lf->datablocks[i];
	colsym=(psnsym *)malloc(sizeof(psnsym)*(db->ncol+1));
	for ( j=0 ; j<db->ncol ; j++ )
	 {	colsym[j].type=T_VAR;
		colsym[j].major=j+nvar;
		colsym[j].name=db->cols[j].name;
	 }
	colsym[db->ncol].type=0;
	colsym[db->ncol].major=0;
	colsym[db->ncol].name=NULL;
	db->colsym=colsym;
  }
 
 if ( is_fit )
  {	mysyms[0]=varsym,
	mysyms[1]=lpg.pl_sym,
	mysyms[2]=NULL;
	i=build_psn_constraint_sequences(lf,cntstr,mysyms,vars,nvar);
	if ( i )
		exit_with_usage(i);
	if ( lf->fconstraint.nc > nvar )
	 {	fprint_error("too many constraints");
		return(1);
	 }
  }

 for ( i=0 ; i<lf->ndatablock ; i++ )
	lf->datablocks[i].xnoise=0.0;	/* no extra (white)noise by default  */

 if ( perturbstr != NULL )
  {	
	int		i,n;
	double		xnoise;

	if ( nkey>0 )
	 {	char		*tmp,**cmd,*arg[4];
		datablock	*db;

		tmp=strdup(perturbstr);
		cmd=tokenize_char_dyn(tmp,',');

		for ( i=0 ; cmd != NULL && cmd[i] != NULL ; i++ )
		 {	n=tokenize_char(cmd[i],arg,':',2);
			if ( n<2 )
			 {	fprint_error("invalid syntax in the perturbation parameters");
				return(1);
			 }
			db=datablock_search_by_key(lf,arg[0]);
			if ( db==NULL )
			 {	fprint_error("invalid datablock key '%s'",arg[0]);
				return(1);
			 }
			if ( sscanf(arg[1],"%lg",&xnoise)<1 || xnoise<0.0 )
			 {	fprint_error("invalid noise level '%s' near datablock key '%s'",arg[1],arg[0]);
				return(1);
			 }
			db->xnoise=xnoise;
		 }
		if ( cmd != NULL )	free(cmd);
		free(tmp);
	 }
	else if ( lf->ndatablock <= 0 )
	 {	fprint_error("no datablocks are defined");
		return(1);
	 }
	else
	 {	if ( sscanf(perturbstr,"%lg",&xnoise)<1 || xnoise<0.0 )
		 {	fprint_error("invalid noise level '%s'",perturbstr);
			return(1);
		 }
		lf->datablocks[0].xnoise=xnoise;
	 }
		
  }

 mysyms[0]=varsym,		/* symbols for fit variables		     */
 mysyms[1]=NULL,		/* gonna be replaced by column symbol table  */
 mysyms[2]=lpg.pl_sym,		/* general operators			     */
 mysyms[3]=NULL;

 if ( seed<0 )
  {	struct timeval	tv;
	int		fd;

	if ( (fd=open("/dev/urandom",O_RDONLY)) >= 0 )
	 {	(void)read(fd,&seed,sizeof(int));
		close(fd);
	 }
	else
 	 {	gettimeofday(&tv,NULL);
		seed=tv.tv_usec | ((tv.tv_sec & 0xFFF) << 20);
	 }
  }
 srand48(seed);

 if ( ! is_fit ) 
  {	datablock	*db;
	int		j;

	for ( i=0 ; i<lf->ndatablock ; i++ )
	 {	db=&lf->datablocks[i];
		mysyms[1]=db->colsym;
		j=build_psn_base_sequences(db,db->fncarg,lfuncts,nlfunct,mysyms,nvar);
		if ( j )	exit_with_usage(j);
	 }


	if ( fw==NULL )	of.f_write=stdout;
	else		of.f_write=fw;

	if ( clostr != NULL /* && (! lf->parameters.map_chi2) */ )
 	 {	ncolout=0;
		colouts=NULL;
		while ( *clostr )
		 {	if ( sscanf(clostr,"%d",&i)<1 || i<=0 )
				break;
			colouts=(int *)realloc(colouts,sizeof(int)*(ncolout+1));
			colouts[ncolout]=i-1;
			ncolout++;
			while ( isdigit(*clostr) )	clostr++;
			if ( *clostr==',' )		clostr++;
		 };
	 }
	else
	 {	ncolout=0;
		colouts=NULL;
	 }

	for ( i=0 ; i<lf->ndatablock ; i++ )
	 {	db=&lf->datablocks[i];

		if ( db->inparg==NULL || strcmp(db->inparg,"-")==0 )
			fr=stdin;
		else
			fr=fopen(db->inparg,"rb");
		if ( fr==NULL )
		 {	fprint_error("unable to open input file");
			return(1);
	 	 }

		evaluate_and_dump(fr,of.f_write,db,vars,nvar,colouts,ncolout,exprs,nexpr);

		if ( i<lf->ndatablock-1 )
			fprintf(of.f_write,"\n");

		fcloseread(fr);
	 }

  }
 else
  {	FILE		*fo,*fb,*fa,*fx,*fv;
	int		i,j,calc_derivatives;
	datablock	*db;

	calc_derivatives=lfit_is_method_requries_derivatives(lf->parameters.fit_method);
	if ( lf->parameters.fit_method==FIT_METHOD_EMCE )
		calc_derivatives|=lfit_is_method_requries_derivatives(lf->parameters.tune.emce.sub_method);

	lf->is_linear=1;
	for ( i=0 ; i<lf->ndatablock ; i++ )
	 {	db=&lf->datablocks[i];
		mysyms[1]=db->colsym;
		j=build_psn_fit_sequences(db,db->fncarg,db->deparg,
			lfuncts,nlfunct,
			db->errarg,
			mysyms,vars,nvar,force_nonlin,calc_derivatives);
		if ( j )		exit_with_usage(j);
		if ( ! db->is_linear )	lf->is_linear=0;
	 }
	if ( varseplist != NULL )
		lfit_set_separated_linears(&lf->fconstraint,vars,nvar,varseplist);

	if ( ofname != NULL )
	 {	fo=fopenwrite(ofname);
		if ( fo==NULL )		
		 {	fprint_error("unable to create output list file '%s' for fitted lines",ofname);
			return(1);	
		 }
	 }
	else	fo=NULL;

	if ( obname != NULL )
	 {	fb=fopenwrite(obname);
		if ( fb==NULL )		
		 {	fprint_error("unable to create output list file '%s' for rejected lines",obname);
			return(1);
		 }
	 }
	else	fb=NULL;

	if ( oaname != NULL )
	 {	fa=fopenwrite(oaname);
		if ( fa==NULL )	
		 {	fprint_error("unable to create output list file '%s' for used lines",oaname);
			return(1);
		 }
	 }
	else	fa=NULL;

	if ( oxname != NULL )
	 {	fx=fopenwrite(oxname);
		if ( fx==NULL )		
		 {	fprint_error("unable to create output expression file '%s'",oxname);
			return(1);
		 }
	 }
	else	fx=NULL;

	if ( ovname != NULL )
	 {	fv=fopenwrite(ovname);
		if ( fv==NULL )
		 {	fprint_error("unable to create variable list file");
			return(1);
		 }
	 }
	else	fv=NULL;

	if ( fo==NULL && fb==NULL && fx==NULL && fw==NULL && fa==NULL && fv==NULL )
		fw=stdout;

	if ( lf->parameters.niter<=0 )
 		lf->parameters.niter=0;
	if ( lf->parameters.sigma<=0.0 )
	 {	lf->parameters.niter=0,
		lf->parameters.sigma=0.0;
	 }

	of.f_write=fw;
	of.f_lused=fo;
	of.f_lrejd=fb;
	of.f_lall =fa;
	of.f_expr =fx;
	of.f_vval =fv;

	switch ( lf->parameters.fit_method )
 	 {   case FIT_METHOD_MCMC:
		fit_markov_chain_monte_carlo(lf,&of,vars,nvar);
		break;
	     case FIT_METHOD_MCHI:
		fit_map_chi2(lf,&of,vars,nvar);
		break;
	     case FIT_METHOD_EMCE:
		fit_error_monte_carlo_estimation(lf,&of,vars,nvar);
		break;
	     case FIT_METHOD_XMMC:
		fit_extended_markov_chain_mc(lf,&of,vars,nvar);
		break;	
	     case FIT_METHOD_DHSX:
		fit_downhill_simplex(lf,&of,vars,nvar);
		break;	
	     case FIT_METHOD_CLLS:
	     case FIT_METHOD_NLLM:
		fit_linear_or_nonlinear(lf,&of,vars,nvar,
		errtype,errdump,resdump,is_dump_delta,is_weighted_sigma,0);
		break;
	     case FIT_METHOD_LMND:
		fit_linear_or_nonlinear(lf,&of,vars,nvar,
		errtype,errdump,resdump,is_dump_delta,is_weighted_sigma,1);
		break;
	     case FIT_METHOD_FIMA:
		fit_fisher_matrix_analysis(lf,&of,vars,nvar,dvars,ndvar);
		break;
	     default:
		fprint_error("internal: unimplemented fit method (code: %d)",lf->parameters.fit_method);
		break;
	 }

	if ( fx != NULL )	fclosewrite(fx);
	if ( fb != NULL )	fclosewrite(fb);
	if ( fo != NULL )	fclosewrite(fo);
	if ( fa != NULL )	fclosewrite(fa);
	if ( fw != NULL )	fclosewrite(fw);
   }

#ifdef	LFIT_ENABLE_DYNAMIC_EXTENSIONS
 for ( i=nldyn-1 ; i>=0 && ldyns != NULL ; i-- )
  {	dlclose(ldyns[i].handle);		}
 if ( ldyns != NULL )
	free(ldyns);
 ldyns=NULL;
 nldyn=0;
#endif

 return(0);
}

/*****************************************************************************/
    
