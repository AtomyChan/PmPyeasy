/*****************************************************************************/
/* longhelp.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Copyright (C) 2008; Pal, A. (apal@szofi.elte.hu)			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Functions for creating nice ``long helps''. Output is intented to be      */
/* compatible with the `help2man` utility.				     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*  This library is free software: you can redistribute it and/or modify     */
/*  it under the terms of the GNU General Public License as published by     */
/*  the Free Software Foundation, either version 3 of the License, or	     */
/*  (at your option) any later version.					     */
/*									     */
/*  This program is distributed in the hope that it will be useful,	     */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of	     */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	     */
/*  GNU General Public License for more details.			     */
/*									     */
/*  You should have received a copy of the GNU General Public License	     */
/*  along with the program.  If not, see <http://www.gnu.org/licenses/>.     */
/*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifndef	HOST_WIN32
#include <sys/ioctl.h>
#endif

#include "longhelp.h"

#ifdef HAVE_LIBINTL_H
#include <libintl.h>
#define _(string)	gettext(string)
#else
#define _(string)	string
#endif

/*****************************************************************************/

#define	LONGHELP_ABORT_ON_MALLOC
#ifdef	LONGHELP_ABORT_ON_MALLOC
#define	realloc_check(ptr,size) \
	do \
	 {	if ( ptr==NULL && size>0 ) \
		 {	fprintf(stderr,"tokenize.c: %s.\n",_("memory exhausted"));\
			abort(); \
		 } \
	 } while(0)
#else
#define	realloc_check(ptr,size) 
#endif
#define	malloc_check(ptr)	realloc_check(ptr,1)

/*****************************************************************************/

static void remove_quotes(char *buff)
{
 int k;
 while ( *buff )
  {	for ( k=0 ; buff[k]=='"' ; )	k++;
	if ( k )	memmove(buff,buff+k,strlen(buff)+1-k);
	else		buff++;
  }
	
}

static int char_is_space(int c)
{
 if ( c==32 || c==13 || c==10 || c==9 )	return(1);
 else					return(0);
}

static char **tokenize_spaces_dyn(char *buff)
{
 int	intoken,inquota,i,n,nm;
 char 	**rtokens;

 nm=16;
 rtokens=(char **)malloc(sizeof(char *)*nm);
 malloc_check(rtokens);
 if ( rtokens==NULL )
	return(NULL);

 intoken=0,inquota=0;n=0;
 while ( *buff )
  {	if ( ( ! char_is_space(*buff) ) && ! intoken )
	 {	rtokens[n]=buff;
		intoken=!0,inquota=0;n++;
		if ( *buff=='"' )	inquota=!0;
		buff++;
		if ( n>=nm-1 )
		 {	nm+=16;
			rtokens=(char **)realloc(rtokens,sizeof(char *)*nm);
			realloc_check(rtokens,sizeof(char *)*nm);
		 }
	 }
	else if ( intoken && ( (char_is_space(*buff) && inquota) || (!char_is_space(*buff)) ) )
	 {	if ( *buff=='"' )	inquota=!inquota;
		buff++;
	 }
	else if ( intoken && ! inquota && char_is_space(*buff) )
	 {	*buff=0,buff++;
		intoken=0;
	 }
	else	buff++;
  };

 rtokens[n]=NULL;

 for ( i=0 ; i<n ; i++ )
  {	remove_quotes(rtokens[i]);	}

 return(rtokens);
}

/*****************************************************************************/

static int longhelp_fprint_description(FILE *fw,int width,int w,int fpad,int pad,char *desc)
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
			for ( i=0 ; i<fpad ; i++ )
			 {	fprintf(fw," ");	}
			if ( w>0 )	l--;
			w=0;
			p=fpad;
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

int longhelp_fprint_entry(FILE *fw,longhelp_entry *entry,int flags,int width)
{
 int	w,fpad,pad;

 w=fprintf(fw," %s",entry->options);

 pad=(w+4+7)&(~7);
 fpad=16;

 if ( width>0 ) 
 	longhelp_fprint_description(fw,width,w,fpad,pad,entry->description);
 else
  {	for ( ; w<pad ; w++ )
	 {	fprintf(fw," ");		}
	fprintf(fw,"%s\n",entry->description);
  }

 return(0);
}

int longhelp_fprint(FILE *fw,longhelp_entry *entry,int flags,int width)
{
 int	lcnt;

 if ( width<0 && isatty(fileno(fw)) )
  {     
#ifdef TIOCGWINSZ
	struct  winsize ws;
        if ( ! ioctl(fileno(fw),TIOCGWINSZ,&ws) )
                width=ws.ws_col-1;
        else
#endif
                width=0;
  }     

 lcnt=0;
 while ( entry != NULL && entry->options != NULL )
  {	if ( entry->description == NULL )
	 {	if ( lcnt>0 )
		 {	fprintf(fw,"\n");		}
		fprintf(fw,"%s\n", entry->options);
		lcnt=0;
	 }
	else
	 { 	longhelp_fprint_entry(fw,entry,flags,width);
		lcnt++;
	 }

	entry++;
  }

 return(0);
}

static char *longhelp_valid_html(char *in)
{
 char	*out;
 int	len,alen;
 len=0;
 alen=256;
 out=malloc(alen);
 while ( *in )
  {	if ( alen <= len+16 )
	 {	alen+=256;
		out=realloc(out,alen);
	 }
	if ( *in=='<' )
	 {	memcpy(out+len,"&lt;",4);len+=4;	}
	else if ( *in=='>' )
	 {	memcpy(out+len,"&gt;",4);len+=4;	}
	else
	 {	out[len]=*in,len++;			}
	in++;
  }
 out[len]=0;
 return(out);
}

static int longhelp_fprint_mediawiki_string(FILE *fw,char *string)
{
 char	*line,**cmd;
 int	i;
 line=longhelp_valid_html(string);
 cmd=tokenize_spaces_dyn(line);
 for ( i=0 ; cmd[i] != NULL ; i++ )
  {	if ( cmd[i][0]=='-' )	fprintf(fw,"'''%s''' ",cmd[i]);
	else			fprintf(fw,"%s ",cmd[i]);
  }

 free(cmd);
 free(line);
 return(0);
}

int longhelp_fprint_mediawiki(FILE *fw,longhelp_entry *entry)
{
 while ( entry != NULL && entry->options != NULL )
  {	if ( entry->description == NULL )
	 {	if ( (entry+1)->options != NULL && strcmp((entry+1)->options,"")==0 && (entry+1)->description==NULL )
			fprintf(fw,": %s\n",entry->options);
		else if ( strcmp(entry->options,"") != 0 )
			fprintf(fw,"=== %s ===\n",entry->options);
		fprintf(fw,"\n");
	 }
	else
	 { 	fprintf(fw,": ");
		longhelp_fprint_mediawiki_string(fw,entry->options);
		fprintf(fw,"\n");
		fprintf(fw,":: ");
		longhelp_fprint_mediawiki_string(fw,entry->description);
		fprintf(fw,"\n");
		fprintf(fw,"\n");
	 }
	entry++;
  }

 return(0);
}

/*****************************************************************************/
                                                                
                
