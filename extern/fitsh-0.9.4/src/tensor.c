/*****************************************************************************/
/* tensor.c 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Standalone library for allocating tensors with arbitrary rank and type.   */
/* (c) 2004, 2005; Pal, A. (apal@szofi.elte.hu).	 		     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in tensor.h	     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include "tensor.h"

/*****************************************************************************/

void *tensor_alloc_arr(int typesize,int rank,int *arr)
{
 size_t	tsize,psize,bsize,cd,cb,snext,pnext;
 size_t	i,j;
 void	*ret,**pret;

 if ( rank<=0 )
	return(NULL);

 psize=0;bsize=1;
 for ( i=rank-1 ; i>=1 ; i-- )
  {	bsize=bsize*arr[i];
	psize+=bsize;
  }
 psize=psize*sizeof(void *);
 tsize=psize+typesize*bsize*arr[0];
 pret=(void **)malloc(tsize);
 if ( pret==NULL )	return(NULL);

 psize=0;bsize=1;
 for ( i=rank-1 ; i>=1 ; i-- )
  {	cd=arr[i];
	if ( i>1 )	snext=sizeof(void *);
	else		snext=typesize;
	cb=cd*bsize;
	pnext=psize+cb;
	for ( j=0 ; j<cb ; j++ )
		pret[psize+j]=(void *)( ((char *)(&pret[pnext]))+j*snext*arr[i-1] );
	psize=pnext;
	bsize=cb;
  }
 ret=(void *)pret;
 return(ret);
}
void *tensor_alloc(int typesize,int rank,...)
{
 int		ip_static[16],*ip_dynamic,*ip,i;
 void		*ret;
 va_list	ap;

 if ( rank<=16 )
  {	ip=ip_static;
	ip_dynamic=NULL;
  }
 else
  {	ip_dynamic=(int *)malloc(sizeof(int)*rank);
	ip=ip_dynamic;
  }
 va_start(ap,rank);
 for ( i=0 ; i<rank ; i++ )
  {	ip[i]=va_arg(ap,int);		}
 va_end(ap);
 ret=tensor_alloc_arr(typesize,rank,ip);
 if ( ip_dynamic != NULL )	free(ip_dynamic);
 return(ret);
}

int tensor_free(void *tensor)
{
 free(tensor);
 return(0);
}

/*****************************************************************************/
                                                   
