/*****************************************************************************/
/* fits.c 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Another simple standalone library for manipulating FITS files:	     */
/* Some functions for nice user-interfaces...				     */
/* (c) 2004-06, Pal, A. (apal@szofi.elte.hu). 				     */
/* See reference(s) at the end of this source code.			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in fits.h	     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#define	_FITS_SOURCE

#include <fits/fits.h>
#include "fits-common.h"

/*****************************************************************************/

static char *local_basename=NULL;

char * fits_basename(char *name,int *rframeno)
{
 int	frameno,len,pos,cnt;
 char	*rname;

 if ( name==NULL )
  {	if ( rframeno != NULL )	*rframeno=-1;	
	return(NULL);
  }

 len=strlen(name);

 if ( local_basename==NULL )
  {	local_basename=strdup(name);	}
 else 
  {	if ( strlen(local_basename)<len )
		local_basename=realloc(local_basename,len+1);
	strcpy(local_basename,name);
  }

 rname=local_basename;

 if ( rname[len-1]==']' )
  {	pos=len-2;cnt=0;
	while ( 0<pos && isdigit((int)rname[pos]) )	pos--,cnt++;
	if ( ! ( rname[pos]=='[' && cnt>0 && pos>0 ) )	pos=0;
  }
 else	pos=0;

 if ( pos>0 )
  {	sscanf(rname+pos+1,"%d",&frameno);
	if ( frameno<=0 )	frameno=-1,rname=name;
	else			frameno--,rname[pos]=0;
  }
 else	frameno=-1,rname=name;

 if ( rframeno != NULL )	*rframeno=frameno;

 return(rname);
}

/*****************************************************************************/

static char fits_nasty_chars[]="()*?[]&|;{}#'";

static int fits_is_nasty_char(int t)
{
 if ( strchr(fits_nasty_chars,t)==NULL )
	return(0);
 else
	return(1); 
}
static int fits_is_any_nasty_char(char *buff)
{
 for ( ; *buff ; buff++ )
  {	if ( fits_is_nasty_char(*buff) )	return(1);	}
 return(0);
}

int fits_header_export_command_line(fits *img,
	char *hdr,char *first,char *prefix,int argc,char *argv[])
{
 char	*buff,wbuff[80];
 int	i,tk,l,p,mxl,plen;

 if ( prefix==NULL )
  {	mxl=65;
	plen=0;
  }
 else
  {	plen=strlen(prefix);
	if ( plen>32 )	plen=32;
	mxl=65-plen;
  }

 if ( first != NULL )
	fits_headerset_set_string(&img->header,hdr,FITS_SH_ADD,first,NULL);

 buff=NULL;p=0;
 for ( i=0 ; i<argc ; i++ )
  {	l=strlen(argv[i]);
	tk=fits_is_any_nasty_char(argv[i]);
	buff=(char *)realloc(buff,p+l+4);
	if ( i==0 )	buff[0]=0;
	else		{ strcat(buff," ");p++; }
	if ( tk )	{ strcat(buff,"\"");p++; }
	strcat(buff,argv[i]);p+=l;
	if ( tk )	{ strcat(buff,"\"");p++; }
  }
 for ( l=0 ; l<p ; )
  {	i=p-l;
	if ( i>mxl )	i=mxl;
	if ( prefix != NULL )
	 {	strncpy(wbuff,prefix,plen);
		memcpy(wbuff+plen,buff+l,i);
		wbuff[plen+i]=0;
	 }
	else
	 {	memcpy(wbuff,buff+l,i);
		wbuff[i]=0;
	 }
	if ( p-l>mxl )	strcat(wbuff,"\\");
	fits_headerset_set_string(&img->header,hdr,FITS_SH_ADD,wbuff,NULL);
	l+=i;
  }
 if ( buff != NULL ) free(buff);
 return(0);
}

/*****************************************************************************/

typedef struct
 {	int	code;
	char	*message;
 } errmsg;

errmsg fits_error_messages[]=
{
 { FITSERR_OPEN		, "unable to open input file"			},
 { FITSERR_ALLOC	, "unsuccessful memory allocation" 		},
 { FITSERR_FRAME_MISSING, "missing frame number specification"		},
 { FITSERR_FRAME_AMBIG	, "ambigous frame number"			},
 { FITSERR_FRAME_INVALID, "invalid frame number"			},
 { FITSERR_IMAGE	, "unable to treat input as a FITS image"	},
 { FITSERR_SCALE	, "unable to rescale inpt image"		},
 { -1,NULL }
};

char *fits_error(int errcode)
{
 errmsg	*e;
 for ( e=fits_error_messages ; e->code>=0 && e->message != NULL ; e++ )
  {	if ( e->code == errcode )	return(e->message);		}
 return(NULL);
}

/*****************************************************************************/

FILE * fits_file_open(char *name)
{
 FILE	*fr;

 if ( name==NULL )
	fr=NULL;
 else if ( strcmp(name,"-")==0 )	
	fr=stdin;
 else
	fr=fopen(name,"rb");

 return(fr);
}
FILE * fits_file_create(char *name)
{
 FILE	*fw;

 if ( name==NULL )
	fw=NULL;
 else if ( strcmp(name,"-")==0 )	
	fw=stdout;
 else
	fw=fopen(name,"wb");

 return(fw);
}

int fits_file_close(FILE *f)
{
 int	ret;
 if ( fileno(f) != fileno(stdin) && fileno(f) != fileno(stdout) )
	ret=fclose(f);
 else
	ret=0;
 return(ret);
}

/*****************************************************************************/
                                                          
              
