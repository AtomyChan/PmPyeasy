/*****************************************************************************/
/* none.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* This library emulates the behaviour of the calling of the other parallel_ */
/* function(s) in a single-process environment.				     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2006; Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parallel.h"

int parallel_none(void *base,int nmemb,int size,int (*funct)(const void *param,void *dpnt),const void *param,int flags)
{
 char	*cbase;
 int	i,ret,cnt;
 void	*dpnt;
 
 cbase=(char *)base;
 
 cnt=0;
 for ( i=0 ; i<nmemb ; i++ )
  {	dpnt=(void *)(cbase+i*size);
	ret=funct(param,dpnt);
	if ( ret )
	 {	if ( flags & PARALLEL_BREAK )	return(1);
		else				cnt++;
	 }
  }
 return(cnt);
}

/*****************************************************************************/
               
