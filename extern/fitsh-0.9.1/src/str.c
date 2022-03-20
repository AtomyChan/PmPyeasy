/*****************************************************************************/
/* str.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Some functions related to (dynamic) string handling.			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2006; Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/

#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>

#include "str.h"

/*****************************************************************************/

char * strkcpy(char *out,char *in,int size)
{
 strncpy(out,in,size);
 out[size-1]=0;
 return(out);
}

/*****************************************************************************/

int strappend(char **str,char *cat)
{
 int	l1,l2;

 if ( str==NULL )
	return(-1);
 else if ( *str==NULL )
  {	*str=strdup(cat);
	return(0);
  }
 else
  {	l1=strlen(*str);
	l2=strlen(cat);
	*str=realloc(*str,l1+l2+1);
	strcpy((*str)+l1,cat);
	return(0);
  }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int vstrappendf(char **str,char *fmt,va_list ap)
{
 int		n,l,size;

 if ( str==NULL )	return(0);

 if ( *str==NULL )	l=0;
 else			l=strlen(*str);
 size=128;

 *str=realloc(*str,l+size);
 if ( *str==NULL )	return(-1);
 while ( 1 )
  {	n=vsnprintf((*str)+l,size,fmt,ap);
	if ( n>-1 && n<size )
		return(0);
	else if ( n>-1 )
		size=n+1;	
	else
		size=size*2;	
	if ( (*str=realloc(*str,l+size))==NULL )
		return(-1);
  };
 return(0);	
}

int strappendf(char **str,char *fmt,...)
{
 int		ret;
 va_list	ap;
 
 va_start(ap,fmt);
 ret=vstrappendf(str,fmt,ap);
 va_end(ap);

 return(ret);
}

/*****************************************************************************/
                                                                    
