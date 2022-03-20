/*****************************************************************************/
/* statistics.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Standalone library for calculating median and mode of a data set.	     */
/* (c) 2001, 2004, Pal, A. (apal@szofi.elte.hu). 			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in statistics.h    */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "statistics.h"

/*****************************************************************************/

static int median_compare(const void *px,const void *py)
{
 if ( *(double *)px < *(double *)py ) 	return(-1);
 else					return(1);
}
double median(double *data,int n)
{
 double m;

 if ( data==NULL || n<=0 )
	return(0.0);

 qsort((void *)data,n,sizeof(double),median_compare);
 if ( n%2 )	m=data[n/2];
 else		m=0.5*(data[n/2-1]+data[n/2]);

 return(m);
}

/*****************************************************************************/
          
