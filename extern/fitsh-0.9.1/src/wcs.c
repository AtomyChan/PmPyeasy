/*****************************************************************************/
/* wcs.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Some functions related to sky coordinate conversions & projections	     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "wcs.h"

/*****************************************************************************/

int wcs_get_projection_matrix(double ra0,double de0,matrix mproj)
{
 double	sr0,sd0,cr0,cd0;

 ra0=M_D2R*ra0;
 de0=M_D2R*de0;

 sr0=sin(ra0);
 cr0=cos(ra0);
 sd0=sin(de0);
 cd0=cos(de0);

 mproj[0][0]=+sr0    ,mproj[0][1]=-cr0    ,mproj[0][2]=0.0;
 mproj[1][0]=-sd0*cr0,mproj[1][1]=-sd0*sr0,mproj[1][2]=+cd0;
 mproj[2][0]=-cd0*cr0,mproj[2][1]=-cd0*sr0,mproj[2][2]=-sd0;

 return(0);
}

int wcs_get_projected_coords_matrix(matrix mproj,double ra,double de,double *rx,double *ry)
{
 double	sr,sd,cr,cd,x,y,z;

 ra=M_D2R*ra;
 de=M_D2R*de;

 sr=sin(ra);
 cr=cos(ra);
 sd=sin(de);
 cd=cos(de);

 x=cd*cr;
 y=cd*sr;
 z=sd;

 *rx=mproj[0][0]*x+mproj[0][1]*y+mproj[0][2]*z;
 *ry=mproj[1][0]*x+mproj[1][1]*y+mproj[1][2]*z;
 
 return(0);
}

int wcs_invert_projected_coords_matrix(matrix mproj,double x,double y,double *rra,double *rde)
{
 double	z,px,py,pz;

 z=1.0-x*x-y*y;
 if ( z<0.0 )	return(1);
 else		z=-sqrt(z);

 px=mproj[0][0]*x+mproj[1][0]*y+mproj[2][0]*z;
 py=mproj[0][1]*x+mproj[1][1]*y+mproj[2][1]*z;
 pz=mproj[0][2]*x+mproj[1][2]*y+mproj[2][2]*z;

 *rde=asin(pz)*M_R2D;
 *rra=atan2(py,px)*M_R2D;
 if ( *rra<0.0 )	*rra+=360.0;
 
 return(0);
}

int wcs_project_distort(int type,double *rx,double *ry)
{
 double	d,m;
 switch ( type )
  {   case WCS_TAN:	/* gnomonic projection				*/
 	m=1.0/sqrt(1-(*rx)*(*rx)-(*ry)*(*ry));
	*rx=(*rx)*m;
	*ry=(*ry)*m;
	break;
     case WCS_ARC:	/* arc projection				*/
	d=sqrt((*rx)*(*rx)+(*ry)*(*ry));
	if ( d>0.0 )	m=asin(d)/d;
	else		m=1.0;
	*rx=(*rx)*m;
	*ry=(*ry)*m;
	break;
     case WCS_SIN:	/* orthographic: do nothing, successfully	*/
	break;
     default:		/* unknown projection 				*/
	return(1);
	break;
 }
 return(0);
}

int wcs_invert_project_distort(int type,double *rx,double *ry)
{
 double	d,m;
 switch ( type )
  {   case WCS_TAN:	/* gnomonic projection				*/
 	m=1.0/sqrt(1+(*rx)*(*rx)+(*ry)*(*ry));
	*rx=(*rx)*m;
	*ry=(*ry)*m;
	break;
     case WCS_ARC:	/* arc projection				*/
	d=sqrt((*rx)*(*rx)+(*ry)*(*ry));
	if ( d>0.0 )	m=sin(d)/d;
	else		m=1.0;
	*rx=(*rx)*m;
	*ry=(*ry)*m;
	break;
     case WCS_SIN:	/* orthographic: do nothing, successfully	*/
	break;
     default:		/* unknown projection 				*/
	return(1);
	break;
 }
 return(0);
}

int wcs_get_projected_coords(double ra,double de,double ra0,double de0,double *rx,double *ry)
{
 double	sinda,cosda,sind0,sind,cosd0,cosd;

 ra =M_D2R*(ra-ra0);
 de0=M_D2R*de0;
 de =M_D2R*de;

 sinda=sin(ra);
 cosda=cos(ra);
 sind0=sin(de0);
 cosd0=cos(de0);
 sind =sin(de);
 cosd =cos(de);

 *rx=+cosd*sinda;
 *ry=-sind0*cosd*cosda+cosd0*sind;

 return(0);
}

/*****************************************************************************/
                     
