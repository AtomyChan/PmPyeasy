/*****************************************************************************/
/* biquad.c 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Standalone library for biquadratic spline interpolation in 2 dimension.   */
/* (c) 2004, Pal, A. (apal@szofi.elte.hu).	 			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in biquad.h	     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "biquad.h"

/*****************************************************************************/

static int intspline_coeff_2d(double **y,int px,int py,int dx,int dy,int n,double **l,int lx,int ly,double *tmp)
{
 double	bet,yprev,ycurr=0.0,lfoll;
 int	j;

 if ( n<=0 )	return(1);

 else if ( n==1 )
  {	l[ly][lx]=l[ly+dy][lx+dx]=y[py][px];
	return(0);
  }
 yprev=y[py][px];
 l[ly][lx]=2.0*yprev;

 bet=1.0;
 for ( j=1 ; j<n ; j++ )
  {	px+=dx,
	py+=dy;

	ycurr=y[py][px];
	tmp[j]=1.0/bet;
	bet=4.0-tmp[j];
	l[ly+dy][lx+dx]=(3.0*(ycurr+yprev)-l[ly][lx])/bet;

	lx+=dx,
	ly+=dy;
	yprev=ycurr;
  }

 tmp[n]=1.0/bet;
 bet=1.0-tmp[n];
 lfoll=l[ly+dy][lx+dx]=(2.0*ycurr-l[ly][lx])/bet;
 for ( j=n-1 ; j>=0 ; j-- )
  {	l[ly][lx] -= tmp[j+1]*lfoll;
	lfoll=l[ly][lx];
	lx-=dx,
	ly-=dy;
  }

 return(0);
}

int biquad_smooth(double **c,int sx,int sy,char **mask)
{
 int	i,j,tmpsize;
 double	*tmp;

 tmpsize=sx;
 if ( tmpsize < sy ) tmpsize=sy;
 tmpsize++;
 tmp=(double *)malloc(sizeof(double)*tmpsize);

 if ( mask==NULL )
  {	for ( i=1 ; i<2*sx+1 ; i+=2 )
	 {	intspline_coeff_2d(c,i,1,0,2,sy,c,i,0,tmp);	}
	for ( i=0 ; i<2*sy+1 ; i++  )
	 {	intspline_coeff_2d(c,1,i,2,0,sx,c,0,i,tmp);	}
  }
 else
  {	int	i0,j0;
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; )
		 {	while ( j<sx &&   mask[i][j] )	j++;
			j0=j;
			while ( j<sx && ! mask[i][j] )	j++;
			if ( j>j0 )
			 {	intspline_coeff_2d(c,2*j0+1,2*i+1,2,0,j-j0,c,2*j0,2*i+1,tmp);	}
		 }
	 }
	for ( j=0 ; j<sx ; j++ )
	 {	for ( i=0 ; i<sy ; )
		 {	while ( i<sy &&   mask[i][j] )	i++;
			i0=i;
			while ( i<sy && ! mask[i][j] )	i++;
			if ( i>i0 )
			 {	intspline_coeff_2d(c,2*j+1,2*i0+1,0,2,i-i0,c,2*j+1,2*i0,tmp);	}
		 }
	 }
	for ( j=0 ; j<sx ; )
	 {	while ( j<sx &&   mask[0][j] )	j++;
		j0=j;
		while ( j<sx && ! mask[0][j] )	j++;
		if ( j>j0 )
		 {	intspline_coeff_2d(c,2*j0+1,0,2,0,j-j0,c,2*j0,0,tmp);	}
	 }
	for ( i=1 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; )
		 {	while ( j<sx &&   ( mask[i-1][j] && mask[i][j] ) )	j++;
			j0=j;
			while ( j<sx && ! ( mask[i-1][j] && mask[i][j] ) )	j++;
			if ( j>j0 )
			 {	intspline_coeff_2d(c,2*j0+1,2*i,2,0,j-j0,c,2*j0,2*i,tmp);	}
		 }
	 }
	for ( j=0 ; j<sx ; )
	 {	while ( j<sx &&   mask[sy-1][j] )	j++;
		j0=j;
		while ( j<sx && ! mask[sy-1][j] )	j++;
		if ( j>j0 )
		 {	intspline_coeff_2d(c,2*j0+1,2*sy,2,0,j-j0,c,2*j0,2*sy,tmp);	}
	 }

  }
 free(tmp);
 return(0);
}

int biquad_coeff(double **y,int sx,int sy,double **c,char **mask)
{
 int	i,j;

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	c[2*i+1][2*j+1]=y[i][j];	}
  }

 biquad_smooth(c,sx,sy,mask);
 return(0);
}

/*****************************************************************************/

static double quad_eval(double *cf,double x)
{
 double	a,b,c,r;
 a= 3*cf[0]-6*cf[1]+3*cf[2];
 b=-4*cf[0]+6*cf[1]-2*cf[2];
 c=   cf[0];
 r=(a*x+b)*x+c;
 return(r);
}

double biquad_eval(double **c,double x,double y)
{
 int	ix,iy;
 double	ycf[3],r;

 ix=(int)x,x-=(double)ix;
 iy=(int)y,y-=(double)iy;

 ix=2*ix,iy=2*iy;

 ycf[0]=quad_eval(&c[iy+0][ix],x);
 ycf[1]=quad_eval(&c[iy+1][ix],x);
 ycf[2]=quad_eval(&c[iy+2][ix],x);
 r=quad_eval(ycf,y);
 return(r);
}

/*****************************************************************************/

static int biquad_quad_fit_coeff_naive(double **c,int x,int y,double *coeff)
{
 x=x*2,y=y*2;
 coeff[0]=0.25*c[y][x]+c[y][x+1]-0.25*c[y][x+2]+c[y+1][x]-c[y+1][x+1]-0.25*c[y+2][x]+0.25*c[y+2][x+2];
 coeff[1]=-0.5*c[y][x]+0.5*c[y][x+2]-4*c[y+1][x]+6*c[y+1][x+1]-2*c[y+1][x+2]+0.5*c[y+2][x]-0.5*c[y+2][x+2];
 coeff[2]=-0.5*c[y][x]-4*c[y][x+1]+0.5*c[y][x+2]+6*c[y+1][x+1]+0.5*c[y+2][x]-2*c[y+2][x+1]-0.5*c[y+2][x+2];
 coeff[3]=3*(c[y+1][x]-2*c[y+1][x+1]+c[y+1][x+2]);
 coeff[4]=c[y][x]-(c[y][x+2]+c[y+2][x])+c[y+2][x+2];
 coeff[5]=3*(c[y][x+1]-2*c[y+1][x+1]+c[y+2][x+1]);
 return(0);
}
int biquad_quad_fit_coeff(double **c,int x,int y,double *coeff)
{
 biquad_quad_fit_coeff_naive(c,x,y,coeff);
 coeff[3]=coeff[3]*2.0,
 coeff[5]=coeff[5]*2.0;
 return(0);
}
double biquad_quad_fit_eval(double **c,double x,double y)
{
 int	ix,iy;
 double	coeff[6],r;
 ix=(int)x,x-=(double)ix;
 iy=(int)y,y-=(double)iy;
 biquad_quad_fit_coeff_naive(c,ix,iy,coeff);
 r=coeff[0]+x*(coeff[1]+coeff[3]*x+coeff[4]*y)+y*(coeff[2]+coeff[5]*y);
 return(r);
}

int biquad_quad_fit_minmax(double **c,int x,int y,double *mmr,double *rx,double *ry)
{
 double	coeff[6],det,tr,mx,my;
 int	ret;
 biquad_quad_fit_coeff(c,x,y,coeff);
 det=coeff[3]*coeff[5]-coeff[4]*coeff[4];
 if ( det<=0.0 )	return(0);
 tr=coeff[3]+coeff[5];
 if ( tr>0.0 )		ret=-1;
 else if ( tr<0.0 )	ret=+1;
 else			return(0);
 mx=-(+coeff[5]*coeff[1]-coeff[4]*coeff[2])/det;
 my=-(-coeff[4]*coeff[1]+coeff[3]*coeff[2])/det;
 if ( mx<0 || my<0 || mx>1.0 || my>1.0 )	return(0);
 if ( mmr != NULL )
	*mmr=coeff[0]-(coeff[1]*coeff[1]*coeff[5]-2*coeff[1]*coeff[2]*coeff[4]+coeff[2]*coeff[2]*coeff[3])/(2.0*det);
 if ( rx != NULL )
	*rx=mx;
 if ( ry != NULL )
	*ry=my;
 return(ret);
}

/*****************************************************************************/

static int biquad_poly_coefficients(double **cf,double *pcf)
{
 double	wa1,wa2,wb1,wb2;

 pcf[0]=cf[0][0];
 pcf[1]=-4.0*cf[0][0]+6.0*cf[0][1]-2.0*cf[0][2];
 wa1   =-4.0*cf[1][0]+6.0*cf[1][1]-2.0*cf[1][2];
 wa2   =-4.0*cf[2][0]+6.0*cf[2][1]-2.0*cf[2][2];
 pcf[2]=3.0*(cf[0][0]-2.0*cf[0][1]+cf[0][2]);
 wb1   =3.0*(cf[1][0]-2.0*cf[1][1]+cf[1][2]);
 wb2   =3.0*(cf[2][0]-2.0*cf[2][1]+cf[2][2]);
 pcf[3]=-4.0*cf[0][0]+6.0*cf[1][0]-2.0*cf[2][0];
 pcf[4]=-4.0*pcf[1]+6.0*wa1-2.0*wa2;
 pcf[5]=-4.0*pcf[2]+6.0*wb1-2.0*wb2;
 pcf[6]=3.0*(cf[0][0]-2.0*cf[1][0]+cf[2][0]);
 pcf[7]=3.0*(pcf[1]-2.0*wa1+wa2);
 pcf[8]=3.0*(pcf[2]-2.0*wb1+wb2);

 return(0);
}

double biquad_int_pixel(double **c,int x,int y)
{
 double	r;
 r=c[2*y+1][2*x+1];
 return(r);
}
double biquad_int_square_pixel(double **c,int x,int y)
{
 double	*cf[3],p[9],s2;

 x=2*x,y=2*y;
 cf[0]=&c[y+0][x];
 cf[1]=&c[y+1][x];
 cf[2]=&c[y+2][x];
 biquad_poly_coefficients(cf,p);
 
 s2=	(p[0]*p[0])+
	(p[0]*(p[1]+p[3]))+
	(p[1]*p[1]+p[3]*p[3]+2.0*p[0]*(p[2]+p[6]))/3.0+
	(p[1]*p[3]+p[0]*p[4]+p[1]*p[2]+p[3]*p[6])/2.0+
	(p[0]*p[7]+p[3]*p[4]+p[1]*p[6]+p[0]*p[5]+p[1]*p[4]+p[3]*p[2])/3.0+
	(p[2]*p[2]+p[6]*p[6])/5.0+
	(p[3]*p[7]+p[4]*p[6]+p[5]*p[1]+p[4]*p[2])/4.0+
	(p[4]*p[4]+2.0*(p[0]*p[8]+p[1]*p[7]+p[3]*p[5]+p[2]*p[6]))/9.0+
	(p[8]*p[3]+p[7]*p[4]+p[5]*p[6]+p[8]*p[1]+p[5]*p[4]+p[7]*p[2])/6.0+
	(p[2]*p[5]+p[6]*p[7])/5.0+
	(p[7]*p[7]+p[5]*p[5]+2.0*(p[6]*p[8]+p[2]*p[8]))/15.0+
	(p[5]*p[7]+p[4]*p[8])/8.0+
	(p[8]*p[7]+p[8]*p[5])/10.0+
	(p[8]*p[8])/25.0;

 return(s2);
}

double biquad_scatter(double **c,int x,int y)
{
 double	s1,s2,s;
 s1=biquad_int_pixel(c,x,y); 
 s2=biquad_int_square_pixel(c,x,y);
 s=s2-s1*s1;
 if ( s<0.0 )	return(0.0);
 else		return(sqrt(s));
}

/*****************************************************************************/

int biquad_diff(double **c,int sx,int sy,double **d,char **mask,int var)
{
 int	i,j;
 double	l,t,r;

 if ( mask==NULL )
  {	if ( ! var )	/* var is zero/false: differentiation by 'x'	*/
	 {	for ( i=0 ; i<=2*sy ; i++ )
		 {	l=c[i][0],t=c[i][1],r=c[i][2];
			for ( j=0 ; j<2*sx ; j+=2 )
			 {	l=c[i][j],t=c[i][j+1],r=c[i][j+2];
				d[i][j+0]=-4*l+6*t-2*r;
				d[i][j+1]=-l+r;
			 }
			d[i][j]=+2*l-6*t+4*r;
		 }
	 }
	else		/* var is nonzero/true: differentiation by 'x'	*/
	 {	for ( j=0 ; j<=2*sx ; j++ )
		 {	l=c[0][j],t=c[1][j],r=c[2][j];
			for ( i=0 ; i<2*sy ; i+=2 )
			 {	l=c[i][j],t=c[i+1][j],r=c[i+2][j];
				d[i+0][j]=-4*l+6*t-2*r;
				d[i+1][j]=-l+r;
			 }
			d[i][j]=+2*l-6*t+4*r;
		 }
	 }
  }
 else	/* **TBD** now it is the same code which is okay but slower	*/
  {	if ( ! var )	/* var is zero/false: differentiation by 'x'	*/
	 {	for ( i=0 ; i<=2*sy ; i++ )
		 {	l=c[i][0],t=c[i][1],r=c[i][2];
			for ( j=0 ; j<2*sx ; j+=2 )
			 {	l=c[i][j],t=c[i][j+1],r=c[i][j+2];
				d[i][j+0]=-4*l+6*t-2*r;
				d[i][j+1]=-l+r;
			 }
			d[i][j]=+2*l-6*t+4*r;
		 }
	 }
	else		/* var is nonzero/true: differentiation by 'x'	*/
	 {	for ( j=0 ; j<=2*sx ; j++ )
		 {	l=c[0][j],t=c[1][j],r=c[2][j];
			for ( i=0 ; i<2*sy ; i+=2 )
			 {	l=c[i][j],t=c[i+1][j],r=c[i+2][j];
				d[i+0][j]=-4*l+6*t-2*r;
				d[i+1][j]=-l+r;
			 }
			d[i][j]=+2*l-6*t+4*r;
		 }
	 }
  }

  return(0);
}

/*****************************************************************************/
                    
