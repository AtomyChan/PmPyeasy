/*****************************************************************************/
/* parallel.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Common functions for parallelization.				     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2006; Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/

#include <stdio.h>

#include "parallel.h"

/*****************************************************************************/

int parallel_general(void *base,int nmemb,int size,
	int (*funct)(const void *param,void *dpnt),const void *param,
	int flags,parallel *pg)
{
 int	ret;

 if ( pg==NULL || funct==NULL || base==NULL )	return(-1);

 switch ( pg->type )
  {  case PARALLEL_TYPE_NONE:
	ret=parallel_none(base,nmemb,size,funct,param,flags);
	break;
     case PARALLEL_TYPE_IPC:
	ret=parallel_ipc(base,nmemb,size,funct,param,flags,pg->p.ipc.nthread);
	break;
     default:
	ret=-1;
	break;
  }

 return(ret);
}
                                 
