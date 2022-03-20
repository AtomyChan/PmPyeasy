/*****************************************************************************/
/* scanarg.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Standalone library for parsing command-line arguments (or other lists).   */
/* (c) 2005, 2006, Pal, A. (apal@szofi.elte.hu)				     */
/* Version: 1.1pre3, last modified: 2006.05.23				     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in scanarg.h	     */
/*****************************************************************************/

#ifdef	HAVE_NO_FNMATCH_H
#define SCANARG_WITHOUT_FNMATCH
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#ifndef	SCANARG_WITHOUT_FNMATCH
#include <fnmatch.h>
#endif

#include "scanarg.h"

#define		MBUFF_LEN		128
#define		is_flag(f,x)		(f&(1<<(x-'A')))

/*****************************************************************************/

static int scanarg_get_format_element(char **rfp,int *rn,int *rf)
{
 int	is_neg,n,ret,flags;
 char	*fp;

 fp=*rfp;
 if ( *fp == '%' )	fp++;

 if ( *fp == '(' )
  {	fp++;
	*rfp=fp;
	return((int)('('));
  }
 flags=0;
 while ( 'A' <= *fp && *fp <= 'Z' )	flags=flags|(1<<(*fp-'A')),fp++;

 if ( *fp == '-' )	is_neg=1,fp++;
 else			is_neg=0;
 n=0;
 while ( '0' <= *fp && *fp <= '9' )	n=n*10+(*fp-'0'),fp++;
 if ( is_neg )		n=-n;
 if ( rn != NULL )	*rn=n;

 while ( 'A' <= *fp && *fp <= 'Z' )	flags=flags|(1<<(*fp-'A')),fp++;
 if ( rf != NULL )	*rf=flags;

 if ( ! isalpha((int)*fp) )	return(-1);

 ret=*fp;
 fp++;
 *rfp=fp;
 return(ret);
}
static int scanarg_is_sep(int t)
{
 if ( t==',' || t==';' || t==':' )
	return(1);
 else	
	return(0);
}
static int scanarg_format_argnum(int fc)
{
 switch ( fc )
  {	case 'e': case 'q': case 'c':
		return(0);
	case 'd': case 'g': case 's': case 'm': case 'r': case 'w':
	case 'a': case 't': case 'l': case 'f': case 'i': case '(':
		return(1);
	default:
		return(-1);
  }
}
static int matchcmp(char *pattern,char *buff)
{
 int	l; 
 if ( pattern==NULL || pattern[0]==0 )	return(1);
 l=strlen(pattern);
 if ( pattern[l-1]=='*' )
  {	while ( l>0 && pattern[l-1]=='*' )	l--;
	if ( memcmp(buff,pattern,l)==0 )	return(0);
	else					return(1);
  }
 else if ( strcmp(buff,pattern)==0 )
	return(0);
 else
	return(1);
}

static int scanarg_check_integer(int farg,int ival)
{
 int	is_ok;
 if ( ! ( is_flag(farg,'P') || is_flag(farg,'Z') || is_flag(farg,'N') ) )
	return(0);	/* no constraints have been set: ival can be anything*/
 is_ok=0;
 if ( is_flag(farg,'P') && ival> 0 )	is_ok=1;
 if ( is_flag(farg,'Z') && ival==0 )	is_ok=1;
 if ( is_flag(farg,'N') && ival< 0 )	is_ok=1;
 if ( is_ok )	return(0);
 else		return(1);
}

int scanarg(int argc,char **argv,int flags,...)
{
 char		*fp,*wp,*mp,*carg,*aptr;
 char		mbuff[MBUFF_LEN],flagmask[128];
 int		a,i,l,m,len,is_matched,argdone,morearg;
 int		is_first,is_scan_sep,is_dd_exp,is_only_flag,is_flag_allowed;
 va_list	ap;

 if ( ! ( flags & SCANARG_NO_SKIP_FIRST ) )	/* skip argv[0], if desired */
  {	argv++;
	argc--;
  }
 argdone=1;is_dd_exp=0;

 memset(flagmask,0,128);
 va_start(ap,flags);
 while ( 1 )
  {	aptr=va_arg(ap,char *);
	if ( aptr==NULL )	break;
	fp=aptr;
	while ( *fp && *fp != ':' )
	 {	if ( fp[0]=='(' && isalnum((int)fp[1]) && fp[2]==')' )
		 {	if ( ! (flags&SCANARG_ALLOW_FLAGS) )
				return(-10);
			else
			 {	i=fp[1];
				flagmask[i]=1;
				fp+=3;
			 }
		 }
		while ( *fp && *fp != ':' && *fp != '|' )	fp++;
		if ( *fp=='|' )	fp++;
	 };
	while ( *fp )
	 {	if ( *fp=='%' )
		 {	int	fc,fn;
			fc=scanarg_get_format_element(&fp,NULL,NULL);
			if ( fc<0 )	return(-2);
			fn=scanarg_format_argnum(fc);
			if ( fn<0 )	return(-2);	
			for ( ; fn>0 ; fn-- )	(void)va_arg(ap,void *);
		 }
		else	fp++;
	 };
  }
 va_end(ap);

 is_only_flag=0;
 while ( argc>0 )
  {	
	if ( ( flags & SCANARG_DASHDASH_EXP ) && strcmp(*argv,"--")==0 )
	 {	argc--,argdone++;
		argv++;
		is_dd_exp=1;
		continue;
	 }

	va_start(ap,flags);
	is_matched=0;

	if ( (*argv)[0]=='-' && (flags&SCANARG_ALLOW_FLAGS) && 
	is_only_flag==0 && (!is_dd_exp) )
	 {	for ( i=1,m=1 ; (*argv)[i] && m ; i++ )
		 {	l=(*argv)[i];
			if ( ! isalnum((int)l) )	m=0;
			else if ( ! flagmask[l] )	m=0;
		 }
		if ( m && i>1 ) is_only_flag=1;
	 }

	while ( ! is_matched )
	 {	aptr=va_arg(ap,char *);
		if ( aptr==NULL )	break;

		fp=aptr;
		wp=fp;len=0;
		while ( *wp && *wp != ':' )	wp++,len++;
		if ( *wp != ':' )
			return(-1);
		wp=fp;l=0;
		while ( l<len && ! is_matched )
		 {	mp=wp;i=0;
			while ( l<len && *wp != '|' )	wp++,l++,i++;
			if ( i>MBUFF_LEN-1 )	i=MBUFF_LEN-1;
			memcpy(mbuff,mp,i);mbuff[i]=0;

			if ( is_dd_exp ) /* if desired, after a `--' all args 
					 are going to be matched by only '*' */
			 {	if ( strcmp(mbuff,"*")==0 )
					is_matched=1;
			 }
			else if ( is_only_flag>0 )
			 {	if ( mbuff[0]=='(' && mbuff[2]==')' )
				 {	if ( (*argv)[is_only_flag]==mbuff[1] )
					 {	is_matched=1,
						is_only_flag++;
					 }
				 }
			 }
			else
			 {	if ( strcmp(mbuff,"*")==0 )
				 {  if ( flags & SCANARG_WILDCARD_DASH ) 
					is_matched=1;
				    else if ( (*argv)[0] != '-' )
					is_matched=1;
				 }	
				else
				 {	
#ifdef	SCANARG_WITHOUT_FNMATCH
				    if ( ! matchcmp(mbuff,*argv) )
					is_matched=1;
#else	
				    if ( flags & SCANARG_USE_STRCMP )
				     {	if ( ! matchcmp(mbuff,*argv) )
						is_matched=1;
				     }
				    else
				     {	if ( ! fnmatch(mbuff,*argv,0) )
						is_matched=1;
				     }
#endif
				 }
			 }
			if ( *wp=='|' && ! is_matched )		wp++,l++;
		 };
		fp=fp+len+1;
		if ( is_matched )	break;
		else
		 {	int	fc,fn;
			while ( *fp )
			 {	if ( *fp=='%' )
				 {	fc=scanarg_get_format_element(&fp,
						NULL,NULL);
					if ( fc<0 )	return(-2);
					fn=scanarg_format_argnum(fc);
					if ( fn<0 )	return(-2);	
					for ( ; fn>0 ; fn-- )
						(void)va_arg(ap,void *);
				 }
				else	fp++;
			 };
		 }
	 };
	if ( ! is_matched )	return(argdone);

	a=0;carg=argv[1];
	is_first=1;is_scan_sep=is_flag_allowed=0;
	while ( *fp )
	 {	while ( *fp==32 )	fp++;
		if ( *fp == '%' )
		 {	int	fc,fn,narg,farg,ival,is_brk,is_static,*ipnt;
			double	dval;
			char	*pval,**ppval,**lval,***plval,*constr;

			narg=farg=0;
			fc=scanarg_get_format_element(&fp,&narg,&farg);
			if ( fc<0 )	return(-2);
			fn=scanarg_format_argnum(fc);
			if ( fn<0 )	return(-2);

			switch ( fc )
			 {   case 'e':			/* error */
				return(argdone);
				break;
			     case 'q':			/* quit, stop scan */
				return(0);
				break;
			     case 'c':			/* continue */
				break;

			     case 'd':		/* -d <integer> 	*/
				if ( is_first )	a++,carg=argv[a],is_first=0;
				if ( argc<a+1 || sscanf(carg,"%d",&ival)<1 )
					return(argdone);
				if ( scanarg_check_integer(farg,ival) )
					return(argdone);
				*(va_arg(ap,int *))=ival;
				is_scan_sep=1;
				break;
			     case 'g':		/* -g <real>		*/
				if ( is_first )	a++,carg=argv[a],is_first=0;
				if ( argc<a+1 || sscanf(carg,"%lg",&dval)<1 )
					return(argdone);
				*(va_arg(ap,double *))=dval;
				is_scan_sep=1;
				break;
			     case 's':		/* -s <string>		*/
				if ( is_first )	a++,carg=argv[a],is_first=0;
				if ( argc<a+1 ) return(argdone);
				pval=carg;
				*(va_arg(ap,char **))=pval;
				is_scan_sep=0;
				break;
			     case 'm':		/* -m <dynamic>		*/
				if ( is_first )	a++,carg=argv[a],is_first=0;
				if ( argc<a+1 )	return(argdone);
				pval=(char *)malloc(strlen(carg)+1);
				strcpy(pval,carg);
				*(va_arg(ap,char **))=pval;
				is_scan_sep=0;
				break;
			     case 'r': case 'a':/* -r A -r B => A,B	*/
				if ( fc=='r' )
				 {	if ( is_first )
						a++,
						carg=argv[a],
						is_first=0;
					if ( argc<a+1 )	
						return(argdone);
				 }
				else
					carg=argv[0];

				ppval=va_arg(ap,char **);
				pval=*ppval;

				l=0;
				if ( is_flag(farg,'C') )	constr=",",l++;
				else if ( is_flag(farg,'S') )	constr=" ",l++;
				else if ( is_flag(farg,'M') )	constr=";",l++;
				else if ( is_flag(farg,'L') )	constr=":",l++;
				else				constr=NULL;
				if ( l>1 )	return(-8);

				if ( pval==NULL )
				 {	len=strlen(carg);
					pval=(char *)malloc(len+1);
					strcpy(pval,carg);
				 }
				else if ( constr != NULL )
				 {	len=strlen(pval)+strlen(carg);
					len+=strlen(constr);
					pval=(char *)realloc(pval,len+1);
					strcat(pval,constr);
					strcat(pval,carg);
				 }
				else
				 {	len=strlen(pval)+strlen(carg);
					pval=(char *)realloc(pval,len+1);
					strcat(pval,carg);
				 }
				*ppval=pval;
				is_scan_sep=0;
				break;
			     case 'w':		/* <switch>		*/
				*(va_arg(ap,char **))=argv[0];
				is_scan_sep=0;
				break;
			     case 'l':		/* <sw> <sw> ...	*/
				plval=va_arg(ap,char ***);
				lval=*plval;
				if ( is_flag(farg,'D') && is_flag(farg,'S') )
					return(-7);
				else if ( is_flag(farg,'D') )	is_static=0;
				else if ( is_flag(farg,'S') )	is_static=1;
				else	is_static=(flags&SCANARG_STATIC_LISTS);
				if ( is_static )
				 {	if ( lval==NULL )	return(-7);
					for ( l=0 ; lval[l] != NULL ; )	l++;
					lval[l]=argv[0];
					lval[l+1]=NULL;
				 }
				else if ( lval==NULL )
				 {	lval=(char **)malloc(2*sizeof(char *));
					lval[0]=argv[0];
					lval[1]=NULL;
				 }
				else
				 {	for ( l=0 ; lval[l] != NULL ; ) l++;
					if ( narg<=0 || ( narg>0 && l<narg ) )
					 {	m=(l+2)*sizeof(char *);
						lval=(char **)realloc(lval,m);
						lval[l]=argv[0];
						lval[l+1]=NULL;
					 }
				 }
				*plval=lval;
				is_scan_sep=0;
				break;
			     case 't':		/* %t:  -i A -i B ...	*/
						/* %Mt: -i A B C ...	*/
				plval=va_arg(ap,char ***);
				lval=*plval;
				if ( is_flag(farg,'D') && is_flag(farg,'S') )
					return(-7);
				else if ( is_flag(farg,'D') )	is_static=0;
				else if ( is_flag(farg,'S') )	is_static=1;
				else	is_static=(flags&SCANARG_STATIC_LISTS);
				if ( is_flag(farg,'M') )	morearg=1;
				else				morearg=0;
				l=0;is_brk=1;
				while ( a+1 < argc )
				 {  if ( strcmp(argv[a+1],"--")==0 && morearg )
				     {	if ( flags & SCANARG_DASHDASH_ARG )
					 {	is_brk=0;
						a++;
						continue;
					 }
				     }	
				    if ( argv[a+1][0]=='-' && is_brk && morearg )
					break;

				    if ( is_static )
				     {	if ( lval==NULL )	return(-7);
					for ( l=0 ; lval[l] != NULL ; )	l++;
					if ( narg<=0 || ( narg>0 && l<narg ) )
					 {	lval[l]=argv[a+1];
						lval[l+1]=NULL;
					 }
				     }
				    else if ( lval==NULL )
				     {	len=2*sizeof(char *);
					lval=(char **)malloc(len);
					lval[0]=argv[a+1];
					lval[1]=NULL;
					l++;
				     }
				    else 
				     {	if ( l==0 )
					 {	while ( lval[l]!=NULL )
							 l++;
					 }
					if ( narg<=0 || ( narg>0 && l<narg ) )
					 {	m=(l+2)*sizeof(char *);
						lval=(char **)realloc(lval,m);
						lval[l]=argv[a+1];
						lval[l+1]=NULL;
						l++;
					 }
				     }
				    a++;
			            if ( ! morearg )	break;
				 };
				*plval=lval;
				is_scan_sep=0;
				break;
			     case 'f':
				ipnt=va_arg(ap,int *);
				if ( is_flag(farg,'N') )	l=narg;
				else				l=1<<narg;
				if ( is_flag(farg,'S') )	(*ipnt) =l;
				else				(*ipnt)|=l;
				is_scan_sep=0;
				is_flag_allowed=1;
				break;
			     case 'i':
				ipnt=va_arg(ap,int *);
				if ( narg==0 )	(*ipnt)++;
				else		(*ipnt)+=narg;
				is_scan_sep=0;
				is_flag_allowed=1;
				break;
			     case '(':
				if ( is_first )	a++,carg=argv[a],is_first=0;
				if ( argc<a+1 )
					return(argdone);
				i=0;
				while ( *fp != ')' && *fp )
				 {	wp=fp;len=0;
					while ( *wp != ')'&& *wp &&*wp != ',' )
						wp++,len++;
					if ( *wp != ')' && *wp != ',' )
						return(-3);
					wp=carg;l=0;
					while ( *wp && *wp != ',' )
						wp++,l++;
					if ( l==len && memcmp(fp,carg,l)==0 )
						break;
					else
					 {	i++;
						fp+=len;
						if ( *fp==',' )	fp++;
					 }
				 };
				if ( *fp==')' )	return(argdone);
				*(va_arg(ap,int *))=i;
				is_scan_sep=1;
				while ( *fp && *fp != ')' )	fp++;
				if ( *fp != ')' )	return(-6);
				fp++;
				break;
			     default:
				return(-3);
				break;
			 }
		 }
		if ( (!is_flag_allowed) && is_only_flag>0 )
			return(-8);
		else if ( scanarg_is_sep(*fp) && is_scan_sep )
		 {	while ( *fp != *carg && *carg )	carg++;
			if ( *fp != *carg )	return(argdone);
			else			carg++;
			fp++;
		 }
		else if ( scanarg_is_sep(*fp) && ! is_scan_sep )
			return(-4);
		else if ( *fp==' ' )
		 {	is_first=1;
			fp++;
		 }
		else if ( *fp && *fp != '%' )	/* unexpected char */
			return(-5);
	 };

	if ( is_only_flag==0 || (*argv)[is_only_flag]==0 )
	 {	a++;
		argdone+=a;
		argv+=a;
		argc-=a;
		is_only_flag=0;
	 }
	va_end(ap);
  };

 return(0);
}
            
/*****************************************************************************/

#define	scanpar_get_format_element(fp,n,f) scanarg_get_format_element(fp,n,f)

static int scanpar_is_sep(int t)
{
 if ( t==',' )
	return(1);
 else	
	return(0);
}
static int scanpar_format_argnum(int fc)
{
 switch ( fc )
  {	case 'e': case 'q': case 'c':
		return(0);
	case 'd': case 'g': case 's': case 'm': case 'f':
		return(1);
	default:
		return(-1);
  }
}

int scanpar(char *par,int flags,...)
{
 va_list	ap;
 char		*aptr,*fp,*wp,*mp,*aparpnt;
 char		mbuff[MBUFF_LEN],wpar[MBUFF_LEN],apar[MBUFF_LEN];
 int		len,l,i,argdone,is_matched;

 argdone=1;
 while ( *par )
  {	
	wp=par;i=0;
	while ( *wp && *wp != '=' && (! scanpar_is_sep(*wp)) && i<MBUFF_LEN-1 )
	 {	wpar[i]=*wp,wp++,i++;			};
	wpar[i]=0;
	if ( *wp=='=' )
	 {	wp++;
		i=0;aparpnt=wp;
		while ( *wp && (! scanpar_is_sep(*wp)) && i<MBUFF_LEN-1 )
		 {	apar[i]=*wp,wp++,i++;			};
		apar[i]=0;
	 }
	else
	 {	apar[0]=0;aparpnt=NULL;			}


	while ( scanpar_is_sep(*wp) )	wp++;
	par=wp;	/* okay...;] */

	va_start(ap,flags);
	is_matched=0;

	while ( ! is_matched )
	 {	aptr=va_arg(ap,char *);
		if ( aptr==NULL )	break;

		fp=aptr;
		wp=fp;len=0;
		while ( *wp && *wp != ':' )	wp++,len++;
		if ( *wp != ':' )		return(-1);
		wp=fp;l=0;

		while ( l<len && ! is_matched )
		 {	mp=wp;i=0;
			while ( l<len && *wp != '|' )	wp++,l++,i++;
			if ( i>MBUFF_LEN-1 )	i=MBUFF_LEN-1;
			memcpy(mbuff,mp,i);mbuff[i]=0;

			if ( strcmp(mbuff,"*")==0 )
				is_matched=1;
			else
			 {	
#ifdef	SCANARG_WITHOUT_FNMATCH
			    if ( ! matchcmp(mbuff,wpar) )
				is_matched=1;
#else	
			    if ( flags & SCANPAR_USE_STRCMP )
			     {	if ( ! matchcmp(mbuff,wpar) )
					is_matched=1;
			     }
			    else
			     {	if ( ! fnmatch(mbuff,wpar,0) )
					is_matched=1;
			     }
#endif
			 }
			if ( *wp=='|' && ! is_matched )		wp++,l++;
		 };
		fp=fp+len+1;
		if ( is_matched )	break;
		else
		 {	int	fc,fn;
			while ( *fp )
			 {	if ( *fp=='%' )
				 {	fc=scanpar_get_format_element(&fp,
						NULL,NULL);
					if ( fc<0 )	return(-2);
					fn=scanpar_format_argnum(fc);
					if ( fn<0 )	return(-2);	
					for ( ; fn>0 ; fn-- )
						(void)va_arg(ap,void *);
				 }
				else	fp++;
			 };
		 }
	 };

	if ( ! is_matched )	return(argdone);

	while ( *fp )
	 {	while ( *fp==32 || *fp==9 )	fp++;
		if ( *fp == '%' )
		 {	int	fc,fn,narg,farg,ival,*ipnt;
			double	dval;
			char	*pval;

			narg=0;
			fc=scanpar_get_format_element(&fp,&narg,&farg);
			if ( fc<0 )	return(-2);
			fn=scanpar_format_argnum(fc);
			if ( fn<0 )	return(-2);

			switch ( fc )
			 {   case 'e':			/* error */
				return(argdone);
				break;
			     case 'q':			/* quit, stop scan */
				return(0);
				break;
			     case 'c':			/* continue */
				break;
			     case 'd':
				if ( sscanf(apar,"%d",&ival)<1 )
					return(argdone);
				if ( scanarg_check_integer(farg,ival) )
					return(argdone);
				*(va_arg(ap,int *))=ival;
				break;
			     case 'g':
				if ( sscanf(apar,"%lg",&dval)<1 )
					return(argdone);
				*(va_arg(ap,double *))=dval;
				break;
			     case 'm':
				if ( aparpnt==NULL )
					pval=NULL;
				else
				 {	pval=(char *)malloc(strlen(apar)+1);
					strcpy(pval,apar);
				 }
				*(va_arg(ap,char **))=pval;
				break;
			     case 's':
				*(va_arg(ap,char **))=aparpnt;
				break;
			     case 'f':
				ipnt=va_arg(ap,int *);
				if ( is_flag(farg,'N') )	l=narg;
				else				l=1<<narg;
				if ( is_flag(farg,'S') )	(*ipnt) =l;
				else				(*ipnt)|=l;
				break;
		 	 }
		 }
		else if ( *fp )
			return(-2);
	 };
	argdone++;
	va_end(ap);


  };

 return(0);

}

/*****************************************************************************/

int scanflag(char *p,int flags,...)
{
 char		*wp,*name;
 int		len,retflag,negate,reset,flag;
 va_list	ap;
 
 retflag=0;

 while ( *p )
  {	while ( isspace((int)*p) )	p++;

	negate=reset=0;
	if ( *p=='/' && (flags & SCANFLAG_ALLOW_NEGATE) )
		negate=1,p++;
	else if ( *p=='-' && (flags & SCANFLAG_ALLOW_RESET) )
		reset=1,p++;

	while ( isspace((int)*p) )	p++;

	len=0;wp=p;
	while ( isalnum((int)*wp) )	wp++,len++;

	va_start(ap,flags);
	flag=0;
	while ( 1 ) 
	 {	name=va_arg(ap,char *);
		if ( name==NULL )
			break;
		flag=va_arg(ap,int);
		if ( memcmp(p,name,len)==0 && strlen(name)==len )
			break;
	 };
	va_end(ap);
	if ( name==NULL )	return(-1);
	else
	 {	if ( reset )		retflag &= ~flag;
		else if ( negate )	retflag |= ~flag;
		else			retflag |=  flag;
	 }
	p=wp;
	while ( isspace((int)*p) )	p++;
	if ( *p==',' )	p++;
	else if ( *p )	return(-1);
  };

 return(retflag & (((unsigned)(-1))>>1) );
}

/*****************************************************************************/
        
