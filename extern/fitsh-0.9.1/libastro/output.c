/*****************************************************************************/
/* output.c 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Prints sexagesimal-formatted angle into a string.			     */
/* This part of the library (output.[ch]) is standalone!		     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2001, 2006; Pal, A. (apal@szofi.elte.hu).			     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>

#include <astro/output.h>

/*****************************************************************************/

char *strdeg(char *buff,double a,int dl,int prec,int type)
{
 int	i,l;
 char	wb[8];
 double	x;

 if ( (type&STRDEG_SIGN) && a>=0.0 )	*buff='+',buff++;
 if ( a<0.0 )				*buff='-',buff++,a=-a;
 for ( l=0 ; l<dl ; l++ )
  {	i=(int)a;
	if ( l==0 )	sprintf(buff,"%d",i);
	else		sprintf(buff,"%.2d",i);
	buff+=strlen(buff);
	if ( type & STRDEG_COLON )	*buff=':',buff++;
	else				*buff=' ',buff++;
	a=(a-(double)i)*60.0;
  };

 x=0.5;
 for ( i=0 ; i<prec ; i++ )	x/=10.0;
 a+=x;

 if ( dl>0 && a<10.0 )	sprintf(wb,"0%%.%df",prec);
 else			sprintf(wb,"%%.%df",prec);
 sprintf(buff,wb,a);

 buff+=strlen(buff);
 return(buff);
}

/*****************************************************************************/
                                                             
                                  
