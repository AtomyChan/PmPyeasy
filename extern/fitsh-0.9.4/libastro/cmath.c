/*****************************************************************************/
/* cmath.c 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Common basic functions used by other modules.			     */
/* This part of the library (cmath.[ch]) is standalone.		             */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2003-2006, Pal A. (apal@szofi.elte.hu).				     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>

#include <astro/cmath.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841968
#endif

/*****************************************************************************/

double sina(double x)	{	return(sin(x*M_PI/180.0));	}
double cosa(double x) 	{	return(cos(x*M_PI/180.0));	}
double asina(double x)	{	return(asin(x)*180.0/M_PI); 	}
double acosa(double x)	{	return(acos(x)*180.0/M_PI); 	}

/*****************************************************************************/

void rotate(double *x,double *y,double phi)
{
 double wx,wy;
 phi=phi*M_PI/180.0;
 wx=(*x)*cos(phi)-(*y)*sin(phi);
 wy=(*x)*sin(phi)+(*y)*cos(phi);
 *x=wx,*y=wy;
}
double getpcoords(double x,double y)
{
 if ( x==0 && y==0 )	return(0.0);
 x=atan2(y,x)*180.0/M_PI;
 if ( x<0.0 )	x+=360.0;
 return(x);
}
void get3dpcoords(double x,double y,double z,double *longi,double *latit)
{
 double	r;

 r=sqrt(x*x+y*y+z*z);
 if ( r==0.0 )	*latit=*longi=0;
 else 		*latit=asina(z/r),*longi=getpcoords(x,y);
}
double	get_angular_distance(double lon1,double lat1,double lon2,double lat2)
{
 double d,r;
 d=sina(lat1)*sina(lat2)+cosa(lat1)*cosa(lat2)*cosa(lon1-lon2);
 if ( d>1.0 )		r=0.0;
 else if ( d<-1.0 )	r=180.0;
 else			r=acosa(d);
 return(r);
}

/*****************************************************************************/

double	hypot3d(double x,double y,double z)
{
 return( sqrt(x*x+y*y+z*z) );
}

/*
double solve_kepler_equ(double m,double ex)
{
 double	e,e0;
 int	n;

 e=m+ex*sin(m)/(1-sin(m+ex)+sin(m));

 for ( n=15,e0=0.0 ; e != e0 && n>0 ; n-- )
  {	e0=e;
	e=e+(m+ex*sin(e)-e)/(1-ex*cos(e));
  };
 return(e);
}
*/

double solve_kepler_equ(double m,double e)
{
 double s,d,d0;
 int    n;

 s=sin(m);

 if ( s==0.0 )
        return(m);

 else if ( e>=0.8 && s*s<0.1 && cos(m)>0.0 )
  {     if ( s>0 )      d=+pow(+6*s,1.0/3.0);
        else            d=-pow(-6*s,1.0/3.0);
  }
 else
  {     d=e*s/(1-sin(m+e)+s);
	if ( d<-e )	d=-e;
	else if ( d>e )	d=+e;
  }

 for ( n=15,d0=0.0 ; d != d0 && n>0 ; n-- )
  {     d0=d;
        d=d-(d-e*sin(m+d))/(1-e*cos(m+d));
  };

 return(m+d);
}


/*****************************************************************************/

void matrix_rotation(matrix m,int ax,double fi)
{
 int	i,j,vi,vj;
 double	cfi,sfi;

 cfi=cosa(fi),sfi=sina(fi);
 ax=ax%3;
 
 for ( i=0 ; i<3 ; i++ )
  {	for ( j=0 ; j<3 ; j++ )
	 {	vi=(i-ax+3)%3,vj=(j-ax+3)%3;
		if ( vi==0 && vj==0 )	m[i][j]=1.0;
		else if ( vi==0 || vj==0 )	m[i][j]=0.0;
		else if ( vi==vj )	m[i][j]=cfi;
		else if ( vi==1 )	m[i][j]=-sfi;
		else			m[i][j]=+sfi;
	 }
  }
}

void matrix_mul(matrix r,matrix m1,matrix m2)	/* r = m1 * m2 ... */
{
 int	i,j;
 matrix	w;

 for ( i=0 ; i<3 ; i++ )
  {	for ( j=0 ; j<3 ; j++ )
		w[i][j]=m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j];
  } 
 for ( i=0 ; i<3 ; i++ )
  {	for ( j=0 ; j<3 ; j++ )
		r[i][j]=w[i][j];
  }
}

/*****************************************************************************/
    
