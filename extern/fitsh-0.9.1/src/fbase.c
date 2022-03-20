/*****************************************************************************/
/* fbase.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2008; Pal, A.							     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include "math/spline/spline.h"

#include "fbase.h"

/*****************************************************************************/

int fbase_spline(double **fbase,int order,int n)
{
 double	*y1,*y2,x;
 int	i,o;

 y1=(double *)malloc(sizeof(double)*(order+1));
 y2=(double *)malloc(sizeof(double)*(order+1));

 for ( o=0 ; o<=order ; o++ )
  {	for ( i=0 ; i<=order ; i++ )
	 {	y1[i]=(i==o?1.0:0.0);	}
	natspline_coeff(y1,order+1,y2);
	for ( i=0 ; i<n ; i++ )
	 {	x=(double)order*(0.5+(double)i)/(double)n;
		fbase[o][i]=natspline_inter(y1,y2,order+1,x);
	 }
  }

 free(y2);
 free(y1);
 
 return(0);
}

int fbase_polynomial(double **fbase,int order,int n)
{
 int	i,o;
 double	x,p0,p1,p2;

 for ( i=0 ; i<n ; i++ )
  {	x=2.0*((0.5+(double)i)/(double)n)-1.0;
	p0=1.0;
	p1=x;
	for ( o=0 ;  o<=order ; o++ )
	 {	fbase[o][i]=p0;
		p2=((double)(2*o+3)*x*p1-(double)(o+1)*p0)/(double)(o+2);
		p0=p1;
		p1=p2;
	 }
  }

 return(0);
}

/*****************************************************************************/
                                                   
                    
