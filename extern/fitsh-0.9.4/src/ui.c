/*****************************************************************************/
/* ui.c									     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (Very) high-level user-interface stuff for the programs of the 'fi' pkg.  */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2004-2009; Pal, A. (apal@szofi.net)				     */
/*****************************************************************************/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>

#include "fitsh.h"
#include "longhelp.h"
#include "io/scanarg.h"

/*****************************************************************************/

char *get_progname(char *a0)
{
 char   *progname,*rbuff;
 progname=a0;
 rbuff=strrchr(progname,'/');
 if ( rbuff != NULL )   progname=rbuff+1;
 return(progname);
}

int fprint_generic_version(FILE *fw,char *arg0,char *name,char *pv,int type)
{
 if ( name==NULL && arg0==NULL )	name=NULL;
 else if ( name==NULL )			name=get_progname(arg0);
 if ( name==NULL )			name="(?)";
 
 switch ( type )
  {   case -1:
	fprintf(fw,"%s %s (%s)\n",name,FITSH_VERSION,FITSH_RELEASE);
	fprintf(fw,"Copyright (C) 1996, 2002, 2004-2008, 2010-2016, 2018-2020; %s <%s>\n",FITSH_MAINT_RNAME,FITSH_MAINT_EMAIL);
	break;
      case -2:
	fprintf(fw,"%s [fitsh-%s]\n",name,FITSH_VERSION);
	break;
  }
 return(0);
}

int fprint_generic_long_help(FILE *fw,int is_wiki,longhelp_entry *help,char *synopsis,char *description)
{
 if ( is_wiki )
  {	fprintf(fw,"=== Synopsis ===\n");
	fprintf(fw,": %s\n\n",synopsis);
	fprintf(fw,"=== Description ===\n");
	fprintf(fw,": %s\n\n",description);
  }
 else
  {	fprintf(fw,"Usage:\t%s\n",synopsis);
	fprintf(fw,"%s\n",description);
	fprintf(fw,"\n");
  }

 if ( is_wiki )
 	longhelp_fprint_mediawiki(fw,help);
 else
  {	longhelp_fprint(fw,help,0,-1);
	fprintf(fw,"\n");
	fprintf(fw,"Report bugs to <%s>, see also https://fitsh.net/.\n",FITSH_MAINT_EMAIL);
  }

 return(0);
}

/*****************************************************************************/

int logmsg(int flag,char *msg,...)
{
 va_list ap;
 FILE	 *fp;
 if ( flag )
  {	fp=stderr;
	va_start(ap,msg);
	vfprintf(fp,msg,ap);
	va_end(ap);
 	fflush(fp);
	return(0);
  }
 else
	return(1);
}

/*****************************************************************************/

int parse_mask_flags(char *par)
{
 int	retflag;

 retflag=scanflag(par,SCANFLAG_ALLOW_NEGATE|SCANFLAG_ALLOW_RESET,
	"none",		MASK_OK,
	"clear",	MASK_CLEAR,
	"fault", 	MASK_FAULT,
	"hot",		MASK_HOT,
	"cosmic",	MASK_COSMIC,
	"outer",	MASK_OUTER,
	"oversaturated",MASK_OVERSATURATED,
	"leaked",	MASK_LEAKED,
	"bloomed",	MASK_LEAKED,
	"saturated",	MASK_SATURATED,
	"interpolated",	MASK_INTERPOLATED,
	"all",		MASK_ALL,
	"bad",		MASK_BAD,
	NULL);
	
 if ( retflag<0 )	return(-1);
 else			return( retflag & MASK_ALL );
}

int parse_fits_data_param(char *par,fitsdataparam *fdp)
{
 int	i,b,type;

 if ( par==NULL )	return(0);

 fdp->nquantizebit=0;

 type=-1;
 i=scanpar(par,SCANPAR_DEFAULT,
	"bscale:%g%f",&fdp->bscale,&fdp->is_scale,
	"bzero:%g%f",&fdp->bzero,&fdp->is_scale,
	"bitpix:%d",&fdp->bitpix,
	"char:%SN0f",&type,
	"byte|unsignedchar|unsigned_char:%SN1f",&type,
	"int|short:%SN2f",&type,
	"word|unsigned|unsignedint|unsigned_int|unsignedshort|unsigned_short:%SN3f",&type,
	"long:%SN4f",&type,
	"dword|unsignedlong|unsigned_long:%SN5f",&type,
	"float:%SN6f",&type,
	"double:%SN7f",&type,
	"quantizebits:%d",&fdp->nquantizebit,
	NULL);

 if ( i )	return(1);

 switch ( type )
  {  case 0:
	fdp->bitpix=  8;fdp->bscale=1.0;fdp->bzero=0.0;
	break;
     case 1:
	fdp->bitpix=  8;fdp->bscale=1.0;fdp->bzero=-(double)(1<<( 8-1));
	break;
     case 2:
	fdp->bitpix= 16;fdp->bscale=1.0;fdp->bzero=0.0;
	break;
     case 3:
	fdp->bitpix= 16;fdp->bscale=1.0;fdp->bzero=-(double)(1<<(16-1));
	break;
     case 4:
	fdp->bitpix= 32;fdp->bscale=1.0;fdp->bzero=0.0;
	break;
     case 5:
	fdp->bitpix= 32;fdp->bscale=1.0;fdp->bzero=-(double)(1<<(32-1));
	break;
     case 6:
	fdp->bitpix=-32;fdp->bscale=1.0;fdp->bzero=0.0;
	break;
     case 7:
	fdp->bitpix=-64;fdp->bscale=1.0;fdp->bzero=0.0;
	break;
  }
 if ( type>=0 )
	fdp->is_scale=1;

 b=fdp->bitpix;
 if ( b && ! ( b==8 || b==16 || b==32 || b==-32 || b==-64 ) )	return(1);

 return(0);
}

/*****************************************************************************/

static char nasty_chars[]="()*?[]&|;";

int is_nasty_char(int t)
{
 if ( strchr(nasty_chars,t)==NULL )	return(0);
 else					return(1); 
}
int is_any_nasty_char(char *buff)
{
 while ( *buff )
  {	if ( is_nasty_char(*buff) )	return(1);
	buff++;
  }
 return(0);
}

/*****************************************************************************/

size_t	parse_max_memory_string(char *maxmemstr)
{
 size_t	maxmem;

 if ( maxmemstr != NULL )
  {	double	mxd;
	int	len,t;
	size_t	pagesize;

	if ( sscanf(maxmemstr,"%lg",&mxd)<1 )
		return(1);	/* unexpected */

	len=strlen(maxmemstr);
	if ( len>0 )	t=maxmemstr[len-1];
	else		t=-1;
	if ( t=='b' || t=='B' )
	 {	if ( len>1 )	t=maxmemstr[len-2];
		else		t=-1;
	 }

	if ( t=='k' || t=='K' )		mxd*=1024.0;
	else if ( t=='m' || t=='M' )	mxd*=1024.0*1024.0;
	else if ( t=='g' || t=='G' )	mxd*=1024.0*1024.0*1024.0;
	else if ( t=='t' || t=='T' )	mxd*=1024.0*1024.0*1024.0*1024.0; 
	else
		return(1);	/* unexpected */

	if ( sizeof(size_t)<=4 && mxd>=1920.0*1024.0*1024.0 )
		mxd=1920.0*1024.0*1024.0;		/* 32bit: 1.875 giga */
	else if ( mxd>=1024*1024.0*1024.0*1024.0*1024.0 )
		mxd=1024*1024.0*1024.0*1024.0*1024.0;	/* 64bit: 1.000 peta */

	maxmem=(size_t)mxd;

#ifdef	HOST_WIN32
	pagesize=4096;
#else
	pagesize=(size_t)getpagesize();
#endif

	if ( maxmem<=pagesize )	maxmem=pagesize;
	maxmem=(maxmem)&(~(pagesize-1));
  }
 else
	maxmem=0;

 return(maxmem);
}

/*****************************************************************************/
           
