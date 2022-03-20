/*****************************************************************************/
/* linkblock.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Functions related to linking/connecting rectangular-shaped blocks	     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2006; Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "linkblock.h"

static int linkblock_sort(const void *p1,const void *p2)
{
 if ( ((linkblock *)p1)->x1 < ((linkblock *)p2)->x1 )	return(-1);
 else							return(1);
}

int linkblock_connect(linkblock *lblocks,int nlblock)
{
 int		i,j,y1,y2,x2;
 linkblock	*lc,*ln,*wc,*wn,*ww;

 if ( lblocks==NULL || nlblock<=0 )	return(1);

 for ( i=0 ; i<nlblock ; i++ )
  {	lblocks[i].prev=NULL;
	lblocks[i].next=NULL;
	lblocks[i].id=i;
  }
 qsort(lblocks,nlblock,sizeof(linkblock),linkblock_sort);
 for ( i=0 ; i<nlblock ; i++ )
  {	lc=&lblocks[i];
	x2=lc->x2;
	y1=lc->y1;
	y2=lc->y2;
	for ( j=i+1 ; lblocks[j].x1<x2 && j<nlblock ; j++ )
	 {	ln=&lblocks[j];
		if ( ln->y2<y1 || y2<ln->y1 )	continue;
		for ( wn=ln ; wn->prev != NULL ; )	wn=wn->prev;
		for ( ww=wn ; ww != NULL ; ww=ww->next )
		 {	if ( ww==lc )	break;		}
		if ( ww != NULL )	continue;
		for ( wc=lc ; wc->next != NULL ; )	wc=wc->next;
		wc->next=wn;
		wn->prev=wc;
	 }
  }
 return(0);
}

/*****************************************************************************/
           
