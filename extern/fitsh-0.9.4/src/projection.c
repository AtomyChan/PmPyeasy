/*****************************************************************************/
/* wcs.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Some functions related to sky coordinate conversions & projections	     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "projection.h"

/*****************************************************************************/

int projection_get_matrix_quat(double qr,double qi,double qj,double qk,matrix mproj)
{
 double	inv_s,s;
 double	rxx,rxy,rxz,ryx,ryy,ryz,rzx,rzy,rzz;

 if ( (inv_s=qr*qr+qi*qi+qj*qj+qk*qk)<=0.0 )
	return(-1);

 s=1.0/inv_s;

 rxx=1-2*s*(qj*qj+qk*qk);
 rxy=2*s*(qi*qj-qk*qr);
 rxz=2*s*(qi*qk+qj*qr);
 ryx=2*s*(qi*qj+qk*qr);
 ryy=1-2*s*(qi*qi+qk*qk);
 ryz=2*s*(qj*qk-qi*qr);
 rzx=2*s*(qi*qk-qj*qr);
 rzy=2*s*(qj*qk+qi*qr);
 rzz=1-2*s*(qi*qi+qj*qj);

 mproj[0][0]=+rxy,mproj[0][1]=+ryy,mproj[0][2]=+rzy;
 mproj[1][0]=+rxx,mproj[1][1]=+ryx,mproj[1][2]=+rzx;
 mproj[2][0]=-rxz,mproj[2][1]=-ryz,mproj[2][2]=-rzz;

 return(0); 
}

int projection_get_matrix_rdr(double ra0,double de0,double roll0,matrix mproj)
{
 double	sa0,sd0,sr0,ca0,cd0,cr0;

 ra0=M_D2R*ra0;
 de0=M_D2R*de0;
 roll0=M_D2R*roll0;

 sa0=sin(ra0);
 ca0=cos(ra0);
 sd0=sin(de0);
 cd0=cos(de0);
 sr0=sin(roll0);
 cr0=cos(roll0);

/*
 mproj[0][0]=+sa0    ,mproj[0][1]=-ca0    ,mproj[0][2]=0.0;
 mproj[1][0]=-sd0*ca0,mproj[1][1]=-sd0*sa0,mproj[1][2]=+cd0;
 mproj[2][0]=-cd0*ca0,mproj[2][1]=-cd0*sa0,mproj[2][2]=-sd0;
*/

 mproj[0][0]=+cr0*sa0+sr0*sd0*ca0,mproj[0][1]=-cr0*ca0+sr0*sd0*sa0,mproj[0][2]=-sr0*cd0;
 mproj[1][0]=+sr0*sa0-cr0*sd0*ca0,mproj[1][1]=-sr0*ca0-cr0*sd0*sa0,mproj[1][2]=+cr0*cd0;
 mproj[2][0]=-cd0*ca0,            mproj[2][1]=-cd0*sa0,            mproj[2][2]=-sd0;

 return(0);
}

int projection_do_matrix_coord(matrix mproj,double ra,double de,double *rx,double *ry,double *rz)
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
 *rz=mproj[2][0]*x+mproj[2][1]*y+mproj[2][2]*z;
 
 return(0);
}

int projection_do_inverse_matrix_coord(matrix mproj,double x,double y,double *rra,double *rde)
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

int projection_do_distortion(int type,projection_distort *dist,double *rx,double *ry,double *rz)
{
 double	d,m;
 switch ( type )
  {   case PROJECTION_TAN:	/* gnomonic projection */
 	m=1.0/sqrt(1-(*rx)*(*rx)-(*ry)*(*ry));
	*rx=(*rx)*m;
	*ry=(*ry)*m;
	break;
     case PROJECTION_ARC:	/* arc projection */
	d=sqrt((*rx)*(*rx)+(*ry)*(*ry));
	if ( d>0.0 && *rz < 0.0 )	m=asin(d)/d;
	else if ( d>0.0 )		m=(M_PI-asin(d))/d;
	else				m=1.0;
	*rx=(*rx)*m;
	*ry=(*ry)*m;
	break;
     case PROJECTION_SIN:	/* orthographic: do nothing, successfully */
	break;
     default:		/* unknown projection */
	return(1);
	break;
 }
 if ( dist != NULL && ( 0<dist->brown_conrady_radial_ncoeff || 0<dist->brown_conrady_tangential_ncoeff ) )
  {	int	i;
	double	r2,rc,w,ux,uy,x,y;
	x=(*rx);
	y=(*ry);
	rc=r2=x*x+y*y;
	w=0.0;
	for ( i=0 ; i<dist->brown_conrady_radial_ncoeff ; i++ )
	 {	w += dist->brown_conrady_radial_coeffs[i]*rc;
		rc = rc * r2;
	 }
	ux=uy=0.0;
	if ( 0<dist->brown_conrady_tangential_ncoeff )
	 {	double	p1,p2,p3,p4;
		p1=dist->brown_conrady_tangential_coeffs[0];
		p2=dist->brown_conrady_tangential_coeffs[1];
		ux=p1*(r2+2*x*x)+p2*(2*x*y);
		uy=p1*(2*x*y)+p2*(r2+2*y*y);
		p3=dist->brown_conrady_tangential_coeffs[2];
		p4=dist->brown_conrady_tangential_coeffs[3];
		if ( 2<dist->brown_conrady_tangential_ncoeff )
		 {	ux=ux*(1+p3*r2+p4*r2*r2);
			uy=uy*(1+p3*r2+p4*r2*r2);
		 }
	 }
	*rx=x*(1+w)+ux;
	*ry=y*(1+w)+uy;
  }
 return(0);
}

#define		MAX_ITER_INVERSE_DISTORT	12
#define		PRECISION_INVERSE_DISTORT	(1e-12)

int projection_do_inverse_distortion(int type,projection_distort *dist,double *rx,double *ry,double *rz)
{
 double	d,m;

 if ( dist != NULL && 0<dist->brown_conrady_radial_ncoeff )
  {	int	i;
	double	r2,w,q1,q2,q3,q4;
	double	w1,w2,w4,w6,w8,dw;

	r2=(*rx)*(*rx)+(*ry)*(*ry);

	q1=dist->brown_conrady_radial_coeffs[0]*r2;
	q2=(1<dist->brown_conrady_radial_ncoeff?dist->brown_conrady_radial_coeffs[1]*r2*r2:0.0);
	q3=(2<dist->brown_conrady_radial_ncoeff?dist->brown_conrady_radial_coeffs[2]*r2*r2*r2:0.0);
	q4=(3<dist->brown_conrady_radial_ncoeff?dist->brown_conrady_radial_coeffs[3]*r2*r2*r2*r2:0.0);
	w=-q1+3*q1*q1-q2;

	if ( dist->brown_conrady_radial_ncoeff==1 )
	 {	for ( i=0 ; i<MAX_ITER_INVERSE_DISTORT ; i++ )
		 {      w1=(1+w);
			w2=w1*w1;
			dw=(w+w1*(q1*w2))/(1+3*q1*w2);
			w -= dw;
			if ( fabs(dw)<=PRECISION_INVERSE_DISTORT )	break;
		 }
         }
	else if ( dist->brown_conrady_radial_ncoeff==2 )
	 {	for ( i=0 ; i<MAX_ITER_INVERSE_DISTORT ; i++ )
		 {      w1=(1+w);
			w2=w1*w1;
			w4=w2*w2;
			dw=(w+w1*(q1*w2+q2*w4))/(1+3*q1*w2+5*q2*w4);
			w -= dw;
			if ( fabs(dw)<=PRECISION_INVERSE_DISTORT )	break;
		 }
         }
	else if ( dist->brown_conrady_radial_ncoeff==3 )
	 {	for ( i=0 ; i<MAX_ITER_INVERSE_DISTORT ; i++ )
		 {      w1=(1+w);
			w2=w1*w1;
			w4=w2*w2;
			w6=w2*w4;
			dw=(w+w1*(q1*w2+q2*w4+q3*w6))/(1+3*q1*w2+5*q2*w4+7*q3*w6);
			w -= dw;
			if ( fabs(dw)<=PRECISION_INVERSE_DISTORT )	break;
		 }
         }
	else if ( dist->brown_conrady_radial_ncoeff==4 )
	 {	for ( i=0 ; i<MAX_ITER_INVERSE_DISTORT ; i++ )
		 {      w1=(1+w);
			w2=w1*w1;
			w4=w2*w2;
			w6=w2*w4;
			w8=w4*w4;
			dw=(w+w1*(q1*w2+q2*w4+q3*w6+q4*w8))/(1+3*q1*w2+5*q2*w4+7*q3*w6+9*q4*w8);
			w -= dw;
			if ( fabs(dw)<=PRECISION_INVERSE_DISTORT )	break;
		 }
         }
	*rx=(*rx)*(1+w);
	*ry=(*ry)*(1+w);
  }

 switch ( type )
  {   case PROJECTION_TAN:	/* gnomonic projection	*/
 	m=1.0/sqrt(1+(*rx)*(*rx)+(*ry)*(*ry));
	*rx=(*rx)*m;
	*ry=(*ry)*m;
	break;
     case PROJECTION_ARC:	/* arc projection	*/
	d=sqrt((*rx)*(*rx)+(*ry)*(*ry));
	if ( d>0.0 )	m=sin(d)/d;
	else		m=1.0;
	*rx=(*rx)*m;
	*ry=(*ry)*m;
	break;
     case PROJECTION_SIN:	/* orthographic: do nothing, successfully */
	break;
     default:		/* unknown projection */
	return(1);
	break;
 }
 return(0);
}

int projection_do_coord(double ra,double de,double ra0,double de0,double roll0,double *rx,double *ry,double *rz)
{
 double	sinda,cosda,sind0,sind,cosd0,cosd;
 double	cr,sr,x,y;

 ra =M_D2R*(ra-ra0);
 de0=M_D2R*de0;
 de =M_D2R*de;

 sinda=sin(ra);
 cosda=cos(ra);
 sind0=sin(de0);
 cosd0=cos(de0);
 sind =sin(de);
 cosd =cos(de);

 cr=cos(roll0*M_D2R);
 sr=sin(roll0*M_D2R);

 x=+cosd*sinda;
 y=-sind0*cosd*cosda+cosd0*sind;

 *rx=x*cr-y*sr;
 *ry=x*sr+y*cr;
 *rz=+cosd0*cosd*cosda+sind0*sind;

 return(0);
}

/*****************************************************************************/
                                                                   
