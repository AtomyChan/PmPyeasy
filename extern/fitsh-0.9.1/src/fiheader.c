/*****************************************************************************/
/* fiheader.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line tool for exporting FITS headers to file/standard output.     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2004, 2006; Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/
#define	FI_HEADER_VERSION	"0.9c1"
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include <fnmatch.h>

#include <fits/fits.h>

#include "longhelp.h"
#include "fi.h"

#include "io/scanarg.h"
#include "io/iof.h"
#include "io/tokenize.h"

/*****************************************************************************/

#define		FORMAT_FILE		1
#define		FORMAT_KEYWORD		2
#define		FORMAT_VALUE		3
#define		FORMAT_COMMENT		4

#define		ALTER_UPDATE		1
#define		ALTER_APPEND		2
#define		ALTER_DELETE		3

/*****************************************************************************/

typedef struct
 {	char	keyword[16];
	int	occurrence;
	char	value[72];
	char	comment[72];
 } alter;

/*****************************************************************************/

char *reserved_keywords[]= 
 { "SIMPLE", "EXTEND", "XTENSION", "NAXIS", "+NAXIS", "BITPIX", "END", NULL };

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

int fprint_header_value(FILE *fw,fitsheader *hdr)
{
 if ( hdr->vtype==FITS_VINT )
	fprintf(fw,"%d",hdr->vint);
 else if ( hdr->vtype==FITS_VBOOLEAN )
  {	if ( hdr->vint )	fprintf(fw,"T");
	else			fprintf(fw,"F");
  }
 else if ( hdr->vtype==FITS_VDOUBLE )
	fprintf(fw,"%.15g",hdr->vdouble);
 else if ( hdr->vtype==FITS_VSTR )
	fprintf(fw,"%s",hdr->vstr);

 return(0);
}

/*****************************************************************************/

int fprint_formatted_keyword(FILE *fw,char *sfbuff,char *name,char *keyw,int *format,fitsheader *hdr)
{
 int	f;

 for ( f=0 ; format[f]>=0 ; f++ )
  {	if ( f>0 )
	 {	fprintf(fw," ");	}	

	switch ( format[f] )
	 {   case FORMAT_FILE:
		fprintf(fw,sfbuff,name);
		break;
	     case FORMAT_KEYWORD:
		fprintf(fw,"%-8s",keyw);
		break;
	     case FORMAT_VALUE:
		if ( hdr==NULL )
		 {	fprintf(fw,"NULL");			}
		else
		 {	fprint_header_value(fw,hdr);		}
		break;
	     case FORMAT_COMMENT:
		if ( hdr != NULL )
		 {	fprintf(fw,"# %s",hdr->comment);	}
		break;
  	 }
  }
 fprintf(fw,"\n");
 return(0);
}

int is_any_wildcard(char *buff)
{
 while ( *buff )
  {	if ( *buff=='[' || *buff==']' || *buff=='*' || *buff=='?' )
		return(1);
	else
		buff++;
  };
 return(0);
}

int fprint_keyword(FILE *fw,fitsheaderset *header,char *name,char *keyw,int *format,int namelen)
{
 char		sfbuff[16];
 int		n,j;
 fitsheader	*hdr;

 sprintf(sfbuff,"%%-%ds",namelen);

 if ( ! is_any_wildcard(keyw) )
  {	n=fits_headerset_get_count(header,keyw);
	if ( n==0 )
	 {	fprint_formatted_keyword(fw,sfbuff,name,keyw,format,NULL);	}
	for ( j=0 ; j<n ; j++ )
	 {	hdr=fits_headerset_get_header(header,keyw,j);
		fprint_formatted_keyword(fw,sfbuff,name,keyw,format,hdr);
	 }
  }
 else
  {	n=0;
	for ( j=0 ; j<header->nhdr ; j++ )
	 {	hdr=&header->hdrs[j];
		if ( ! fnmatch(keyw,hdr->name,0) )
		 {	fprint_formatted_keyword(fw,sfbuff,name,hdr->name,format,hdr);
			n++;
		 }
	 }
  }
 if ( n>0 )	return(0);
 else		return(1);

}

int read_header_values_from_file(FILE *fr,char *extension,char *name,char **getkwlist,FILE *fw,int *format,int namelen)
{ 
 fits		*img;
 fitsheaderset	*header;
 int		i,k,nkw,ret;
 char		**cmd,*kwl;

 if ( getkwlist==NULL )		return(0);

 if ( extension==NULL )
  {	img=fits_create();
	fits_read_header(fr,img);
	header=&img->header;
  }
 else
  {	int	i,x;
	img=fits_read(fr);
	if ( sscanf(extension,"%d",&x)<1 )
		x=-1;
	if ( 0<=x && x<img->nxtn )
	 {	header=&img->xtns[x].header;		}
	else
	 {	fitsheader	*fh;
		int		l,l0;
		header=NULL;
		l0=strlen(extension);
		for ( i=0 ; i<img->nxtn ; i++ )
		 {	fh=fits_headerset_get_header(&img->xtns[i].header,"EXTNAME",0);
			if ( fh==0 || fh->vtype != FITS_VSTR )
				continue;
			l=strlen(fh->vstr);
			while ( 0<l && fh->vstr[l-1]==' ' )	l--;
			if ( strncmp(fh->vstr,extension,l)==0 && l==l0 )
				header=&img->xtns[i].header;
		 }
	 }
  }

 ret=0;

 for ( i=0 ; getkwlist[i] != NULL ; i++ )
  {	kwl=strdup(getkwlist[i]);
	remove_spaces(kwl);
	for ( k=0,nkw=1 ; kwl[k] ; k++ )
	 {	if ( kwl[k]==',' )	nkw++;
		else			kwl[k]=toupper(kwl[k]);
	 }
	cmd=(char **)malloc(sizeof(char *)*(nkw+1));
	tokenize_char(kwl,cmd,',',nkw+1);
	
	for ( k=0 ; k<nkw ; k++ )
	 {	if ( fprint_keyword(fw,header,name,cmd[k],format,namelen) )
			ret=1;
	 }
  }

/*
 fprintf(stderr,"[x]\n");
 fits_free(img);
 fprintf(stderr,"[y]\n");
*/

 return(ret);
}

/*****************************************************************************/

int alter_is_reserved_keyword(char *kw)
{
 int	i;
 char	*rk;
 for ( i=0 ; reserved_keywords[i] != NULL ; i++ )
  {	rk=reserved_keywords[i];
	if ( strcmp(kw,rk)==0 )	
		return(1);
	else if ( rk[0]=='+' && memcmp(kw,rk+1,strlen(rk+1))==0 )
		return(1);
  }
 return(0);
}

int alter_get_action(char *act)
{
 if ( strcmp(act,"-s")==0 || strcmp(act,"--set")==0 )
	return(ALTER_UPDATE);
 else if ( strcmp(act,"-u")==0 || strcmp(act,"--update")==0 )
	return(ALTER_UPDATE);
 else if ( strcmp(act,"-a")==0 || strcmp(act,"--append")==0 )
	return(ALTER_APPEND);
 else if ( strcmp(act,"-d")==0 || strcmp(act,"--delete")==0 )
	return(ALTER_DELETE);
 else
	return(-1);
}

int alter_is_keyword_char(int t)
{
 if ( t && t != '=' && t != ',' && (!isspace(t)) && t != '[' && t != ']' )
	return(1);
 else	
	return(0);
}

int alter_get_alternations(char *list,alter **ralters,int *rnalter,int is_set)
{
 alter	*alters,*wa;
 int	nalter;
 char	*p,sk[16],sv[96],sc[96];
 int	lk,lv,lc,inq,occurrence,val_set;

 p=list;
 
 alters=NULL;
 nalter=0;

 while ( *p )
  {	while ( isspace(*p) )	p++;
	lk=0;
	while ( alter_is_keyword_char(*p) && lk<8 )
	 {	sk[lk]=toupper(*p),lk++,p++;		}
	sk[lk]=0;
	while ( isspace(*p) )		p++;
	if ( *p=='[' )
	 {	p++;
		if ( sscanf(p,"%d",&occurrence)<1 )	occurrence=0;
		while ( *p && *p != ']' )	p++;
		if ( *p==']' )	p++;
		while ( isspace(*p) )		p++;
	 }
	else
	 {	occurrence=0;		}

	if ( *p=='=' )		p++,val_set=1;	/* we're trying to set value */
	else if ( *p==',' )	p++,val_set=0;	/* no, go to next keyword    */
	else			val_set=0;	/* this should not happen... */

	if ( val_set )
	 {	while ( isspace(*p) )		p++;
		lv=0;inq=0;
		while ( *p && lv<80 )
		 {	if ( *p=='\'') 	inq=!inq,sv[lv]=*p,lv++;
			else if ( ! inq && ( *p==',' || *p==0 || *p=='/' ) )	break;
			else 		sv[lv]=*p,lv++;
			p++;
		 }	
		sv[lv]=0;
		while ( isspace(*p) )	p++;
		if ( *p=='/' )
		 {	p++;
			lc=0;
			while ( *p && lc<80 && *p != ',' )
			 {	sc[lc]=*p,lc++,p++;		}
			sc[lc]=0;
		 }
		else
		 {	sc[0]=0,lc=0;		}

		while ( *p && *p != ',' )	p++;
		if ( *p==',' )	p++;
	 }
	else
	 {	sv[0]=0,lv=0;
		sc[0]=0,lv=0;
	 }

	if ( lk<=0 || ( is_set && lv<=0 ) )
		continue;

	alters=(alter *)realloc(alters,sizeof(alter)*(nalter+1));
	wa=&alters[nalter];
	nalter++;

	strcpy(wa->keyword,sk);
	wa->occurrence=occurrence;
	if ( is_set )
	 {	strcpy(wa->value,sv);
		strcpy(wa->comment,sc);
	 }
	else
	 {	wa->value[0]=0;
		wa->comment[0]=0;
	 }
  };

 if ( ralters != NULL )	*ralters=alters;
 if ( rnalter != NULL )	*rnalter=nalter;

 return(0);
}

int alter_header_set(fitsheaderset *header,char *updatelist,int policy,int force)
{
 alter	*alters,*wa;
 int	nalter,i,j,k,lv;
 char	*scomment,*sv;
 double	d;

 alter_get_alternations(updatelist,&alters,&nalter,1);
 if ( alters==NULL || nalter<=0 )	return(0);

 for ( i=0 ; i<nalter ; i++ )
  {	wa=&alters[i];

	if ( wa->comment[0] )	scomment=wa->comment;
	else			scomment=NULL;

	if ( alter_is_reserved_keyword(wa->keyword) )
	 {	if ( ! force )
		 {	fprintf(stderr,"Warning: altering reserved keyword '%s' is not allowed, skipping.\n",wa->keyword);
			continue;
		 }
		else
		 {	fprintf(stderr,"Warning: altering reserved keyword '%s'.\n",wa->keyword);	}
	 }

	sv=wa->value;
	lv=strlen(sv);

	if ( strcasecmp(sv,"f")==0 || strcasecmp(sv,"false")==0 )
	 {	fits_headerset_set_boolean(header,wa->keyword,policy,0,scomment);	}
	else if ( strcasecmp(sv,"t")==0 || strcasecmp(sv,"true")==0 )
	 {	fits_headerset_set_boolean(header,wa->keyword,policy,1,scomment);	}

	else if ( sv[0]=='\'' && sv[lv-1]=='\'' )
	 {	memmove(sv,sv+1,lv);lv--;
		sv[lv-1]=0,lv--;
		for ( j=0 ; j<lv ; j++ )
		 {	if ( sv[j]=='\'' )
			 {	memmove(sv+j,sv+j+1,lv-j);
				lv--;
			 }
		 }
		fits_headerset_set_string(header,wa->keyword,policy,sv,scomment);
	 }
	else
	 {	k=1;
		for ( j=0 ; j<lv ; j++ )
		 {	if ( ! (isdigit(sv[j]) || sv[j]=='+' || sv[j]=='-') )
				k=0;
		 }
		if ( k )
		 {	sscanf(sv,"%d",&j);
			fits_headerset_set_integer(header,wa->keyword,policy,j,scomment);
		 }
		else if ( sscanf(sv,"%lg",&d)==1 )
		 {	fits_headerset_set_double(header,wa->keyword,policy,d,scomment);	}
		else 
		 {	fits_headerset_set_string(header,wa->keyword,policy,sv,scomment);	}
	 }
	/*fprintf(stderr,"\"%s\"=\"%s\" [%s]\n",sk,sv,sc);*/
  };

 return(0);
}

int alter_header_delete(fitsheaderset *header,char *deletelist,int force)
{
 alter	*alters,*wa;
 int	nalter,i,n;

 alter_get_alternations(deletelist,&alters,&nalter,0);
 if ( alters==NULL || nalter<=0 )	return(0);
 
 for ( i=0 ; i<nalter ; i++ )
  {	wa=&alters[i];

	if ( alter_is_reserved_keyword(wa->keyword) )
	 {	if ( ! force )
		 {	fprintf(stderr,"Warning: deleting reserved keyword '%s' is not allowed, skipping.\n",wa->keyword);
			continue;
		 }
		else
		 {	fprintf(stderr,"Warning: deleting reserved keyword '%s'.\n",wa->keyword);	}
	 }

	if ( wa->occurrence>0 )
	 {	fits_headerset_delete(header,wa->keyword,wa->occurrence-1);	}
	else if ( wa->occurrence<0 )
	 {	n=fits_headerset_get_count(header,wa->keyword);
		fits_headerset_delete(header,wa->keyword,n+wa->occurrence);
	 }
	else
	 {	fits_headerset_delete_all(header,wa->keyword);			}
  };

 return(0);
}

int alter_header_values(fitsheaderset *header,char **setkwlist,int force)
{
 int	i,action;
 char	*kw;

 if ( setkwlist==NULL )	return(0);

 action=0;
 for ( i=0 ; setkwlist[i] != NULL ; i++ )
  {	if ( setkwlist[i][0]=='-' )
	 {	action=alter_get_action(setkwlist[i]);
		if ( action<0 )	return(-1);	/* it's incosistency! */
	 }
	else
	 {	kw=setkwlist[i];
		switch ( action )
	 	 {   case ALTER_UPDATE:
			alter_header_set(header,kw,FITS_SH_FIRST,force);
			break;
	 	     case ALTER_APPEND:
			alter_header_set(header,kw,FITS_SH_ADD,force);
			break;
		     case ALTER_DELETE:
			alter_header_delete(header,kw,force);
			break;
		 };
	 }
  }

 return(0);
}

/*****************************************************************************/

int fprint_fiheader_usage(FILE *fw)
{
 fprintf(fw,
"Usage:\tfiheader [-h|--help|--long-help] [--version]\n"
"\t[in.fits [-o out.fits|-w|--rewrite]] | [in1.fits [in2.fits ...]]\n");
 fprintf(fw,
"Read header/keyword values:\n"
"\t[--get|--read <KEYWORD1>[,<KEYWORD2>,...]] [--get|--read ...]\n"
"\t[-F|--format file[name],keyword,value,comment]\n");
 fprintf(fw,
"Altering keyword values:\n"
"\t[--set|--update <KEYWORD1>=<value>[/comment][,<KEYWORD2>=...,...]]\n"
"\t[--append <KEYWORD>=value[/comment][,<KEYWORD2>=...,...]]\n"
"\t[--delete <KEYWORD>[<occurrence>][,<KEYWORD2>,...]]\n"
"\t[--force-alter-reserved-keywords]\n");
 fprintf(fw,
"Notes:\n"
" - If a non-existent keyword is specified for --get|--read, it value is\n"
"   going to be the string \"NULL\" (w/o quotation), and the program will\n"
"   return with an exit status (errorlevel) of 2. Even if other specified\n"
"   keywords are exist and their values have been written properly!\n"
" - Altering and reading keyword values at the same time is not permitted.\n");
 fprintf(fw,
" - If header values are altered and only one filename has been specified,\n"
"   the altered FITS file will be written to stdout or into file given by -o.\n"
"   This output file can be the same as the original one, and --rewirte is also\n"
"   allowed.\n"
" - If header values are altered and more than one filename are specified, the\n"
"   same updating schema is used for all files (which are replaced sequentially\n"
"   by the updated ones).\n")  ;
 return(0);
}

longhelp_entry fiheader_long_help[]=
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

 { "Read keyword values:", NULL },
 { "-g, --get, --read <keyword>,...", 
	"List of keywords to get from FITS header." },
 { "-F, --format <format>", 
	"Format of the output, space separated list of the arguments "
	"``filename'', ``keyword'', ``value'' or ``comment'', in the "
	"desired order." },
 { "-x, --extension <extension>",
	"Read keyword values from the given extension." },
 { "Note that shell wildcard character patterns can be used in the "
   "keyword specifications.", NULL }, { "",NULL },

 { "Altering FITS headers:", NULL },
 { "-s, --set, --update <keyword>=<value>[/<comment>],...", 
	"List of keywords with their new values (and optional comments) to "
	"be altered." },
 { "-a, --append <keyword>=<value>[/<comment>],...", 
	"List of keywords with their new values (and optional comments) to "
	"be appended to the existing list of keywords." },
 { "-d, --delete <keyword>,...", "List of keywords to be removed from the "
	"header." },
 { "--force-alter-reserved-keywords", 
	"Enable to alter reserved FITS keywords (use with caution)." },
 { "-w, --rewrite",
	"Re-write the input file(s) directly. Use this option exclusively "
	"from -o|--output. " },

 { NULL, NULL }
};
 

int fprint_fiheader_long_help(FILE *fw)
{
 fprintf(fw,
"Usage:\tfiheader [opions] [input] [-o|--output <output>]\n"
"This program allows the user to read, set, alter or remove a set of values \n"
"associated with the specified keywords from the primary header unit\n"
"(a.k.a. simply ``header'') of a FITS image file.\n\n");

 longhelp_fprint(fw,fiheader_long_help,0,-1);
 
 fprintf(fw,"\n");
 fprintf(fw,"Report bugs to <%s>\n",FI_MAINT_EMAIL);

 return(0);
}

/*****************************************************************************/

static	int	format_v  []= { FORMAT_VALUE,-1};
static	int	format_fv []= { FORMAT_FILE,FORMAT_VALUE,-1 };
static	int	format_kv []= { FORMAT_KEYWORD,FORMAT_VALUE,-1 };
static	int	format_fkv[]= { FORMAT_FILE,FORMAT_KEYWORD,FORMAT_VALUE,-1 };

int main(int argc,char *argv[])
{
 fits		*img;
 FILE		*fr,*fw;
 int		is_help,ret;
 int		i,j,is_single_line,is_rewrite,is_force,incnt,gkwcnt,mxflen;
 char		*outname,**innames,**getkwlist,**setkwlist,*argformat;
 char		*extension;
 int		method;
 int		*format;

 progbasename=strrchr(argv[0],'/');
 if ( progbasename != NULL )	progbasename++;
 else				progbasename=argv[0];

 outname=NULL;
 innames=NULL;
 method=0;is_single_line=0;
 getkwlist=NULL;
 setkwlist=NULL;

 is_help=0,method=0;
 is_rewrite=is_force=0;
 argformat=NULL;
 extension=NULL;

 i=scanarg(argc,argv,SCANARG_ALLOW_FLAGS,
	"--version:%NS-1f%q",&is_help,
	"--version-short|--short-version:%NS-2f%q",&is_help,
	"-h|--help:%f%q",&is_help,
	"--long-help|--help-long:%SN2f%q",&is_help,
	"-o|--output:%s",&outname,
	"-i|--input:%Dt",&innames,
	"-g|--get|-r|--read:%0f%Dt",&method,&getkwlist,
	"-s|--set|-u|--update:%1f%Dl%Dt",&method,&setkwlist,&setkwlist,
	"-a|--append:%1f%Dl%Dt",&method,&setkwlist,&setkwlist,
	"-d|--delete:%1f%Dl%Dt",&method,&setkwlist,&setkwlist,
	"-w|--rewrite:%f",&is_rewrite,
	"-F|--format:%s",&argformat,
	"-x|--extension:%s",&extension,
	"--force-alter-reserved-keywords:%f",&is_force,
	"-*:%e",
	"*:%Dl",&innames,
	NULL);

 if ( i )		
  {	fprint_error("invalid command line argument near '%s'",argv[i]);
	return(1);
  }
 else if ( is_help<0 )
  {	fprint_generic_version(stdout,argv[0],"fiheader",FI_HEADER_VERSION,is_help);
	return(0);
  }
 else if ( is_help>1 )
  {	fprint_fiheader_long_help(stdout);
	return(0);
  }
 else if ( is_help )
  {	fprint_fiheader_usage(stdout);
	return(0);
  }
	
 if ( argformat==NULL )
	format=NULL;
 else
  {	char	*af,*cmd[8];
	int	i,n;
	format=(int *)malloc(8*sizeof(int));
	af=strdup(argformat);
	remove_spaces(af);
	n=tokenize_char(af,cmd,',',7);
	for ( i=0 ; cmd[i] != NULL && i<n ; i++ )
	 {	if ( strcmp(cmd[i],"file")==0 || strcmp(cmd[i],"filename")==0 ||
		strcmp(cmd[i],"name")==0 )
		 	format[i]=FORMAT_FILE;
		else if ( strcmp(cmd[i],"keyword")==0 )
			format[i]=FORMAT_KEYWORD;
		else if ( strcmp(cmd[i],"value")==0 )
			format[i]=FORMAT_VALUE;
		else if ( strcmp(cmd[i],"comment")==0 )
			format[i]=FORMAT_COMMENT;
		else
		 {	fprint_error("invalid format specification '%s'",argformat);
			return(1);
		 }
	 }
	format[i]=-1;
	free(af);
  }
	
 
 if ( (method & 3) == 3  )	
  {	fprint_error("invalid combination of command line arguments");
	return(1);
  }

 if ( innames==NULL )
  {	incnt=0;
	mxflen=0;
  }
 else
  {	mxflen=0;
	for ( incnt=0 ; innames[incnt] != NULL ; ) 
	 {	j=strlen(innames[incnt]);
		if ( j>mxflen )	mxflen=j;
		incnt++;
	 }
  }

 if ( (method & 1) && getkwlist != NULL )	/* read some headers */
  {	gkwcnt=0;
	for ( i=0 ; getkwlist[i] != NULL ; i++ )
	 {	for ( j=0 ; getkwlist[i][j] ; j++ )
		 {	if ( getkwlist[i][j]==',' )	gkwcnt++;	}
		gkwcnt++;
	 }

	if ( format==NULL )
	 {	if ( gkwcnt<=1 )
		 {	if ( incnt<=1 )		format=format_v;
			else			format=format_fv;
		 }
		else
		 {	if ( incnt<=1 )		format=format_kv;
			else			format=format_fkv;
		 }
	 }

	if ( incnt==0 )
	 {	fr=stdin;
		ret=read_header_values_from_file(fr,extension,"-",getkwlist,stdout,format,1);
	 }
	else	
	 {	ret=0;
		for ( i=0 ; i<incnt ; i++ )
		 {	fr=fopenread(innames[i]);
			if ( fr==NULL )	
			 {	fprint_error("unable to open input file '%s'",innames[i]);
				return(1);
			 }
			ret|=read_header_values_from_file(fr,extension,innames[i],getkwlist,stdout,format,1);
			fcloseread(fr);
		 }
	 }
	if ( ret )	return(2);	/* exit with errorlevel 2 */
  }
 else if ( (method & 2) && setkwlist != NULL )	/* alter some headers */
  {	if ( incnt==0 )
	 {	fr=stdin;
		img=fits_read_raw(fr);
		alter_header_values(&img->header,setkwlist,is_force);
		if ( outname==NULL )	fw=stdout;
		else			fw=fopenwrite(outname);
		if ( fw==NULL )		
		 {	fprint_error("unable to create output file '%s'",outname);
			return(1);
		 }
		fits_write(fw,img);
		fclosewrite(fw);
	 }
	else
	 {	if ( incnt>1 && outname != NULL )
		 {	fprint_error("invalid combination of command line arguments (output file name must be omitted)");
			return(1);
		 }
		if ( is_rewrite )	outname=NULL;

		for ( i=0 ; i<incnt ; i++ )
		 {	fr=fopenread(innames[i]);
			if ( fr==NULL )	
			 {	fprint_error("unable to open input file '%s'",innames[i]);
				return(1);
			 }
			img=fits_read_raw(fr);
			alter_header_values(&img->header,setkwlist,is_force);
			if ( outname != NULL )	fw=fopenwrite(outname);
			else			fw=fopenwrite(innames[i]);
			if ( fw==NULL )		
			 {	fprint_error("unable to create output file");
				return(1);
			 }
			fits_write(fw,img);
			fclosewrite(fw);
		 }
	 }
  }

 return(0);

}
