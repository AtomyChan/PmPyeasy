/*****************************************************************************/
/* basis.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Some functions and constans related to the definition of coordinate 	     */	
/* bases in star lists... to be compatible with programs using non-native    */
/* coordinating systems (e.g. IRAF).					     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "basis.h"

char	*basis_comment_string="##CoordBasisShift:";

int basis_set_shift_coords(int basistype,double *rbshx,double *rbshy)
{
 if ( basistype>0 )		*rbshx=*rbshy=+0.5;
 else if ( basistype<0 )	*rbshx=*rbshy=-0.5;
 else				*rbshx=*rbshy=0.0;
 return(0);
}
int basis_write_comment(FILE *fw,double bshx,double bshy)
{
 if ( bshx != 0.0 || bshy != 0.0 )
  {	fprintf(fw,"%s %g, %g",basis_comment_string,bshx,bshy);
	if ( bshx==0.5 && bshy==0.5 )
		fprintf(fw,"\t# Convetion of IRAF/daofind,daophot; ds9; ...");
	fprintf(fw,"\n");
  }
 return(0);
}
int basis_read_comment(char *buff,double *rbshx,double *rbshy)
{
 double	x,y;
 int	len;

 len=strlen(basis_comment_string);

 if ( buff==NULL )
	return(0);
 else if ( memcmp(buff,basis_comment_string,len)==0 )
  {	if ( sscanf(buff+len,"%lg,%lg",&x,&y)<2 )
		return(0);
	else if ( ! isfinite(x) || ! isfinite(y) )
		return(0);
	else
	 {	*rbshx=x,
		*rbshy=y;
		return(1);
	 }
  }
 else
	return(0);
  
}
                                                                
