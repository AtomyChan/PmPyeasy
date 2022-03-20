/*****************************************************************************/
/* polygon.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2015; Pal, A. <apal@szofi.net>					     */
/*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "polygon.h"

/*****************************************************************************/

#define		ix(i)		(poly[2*((i)%(n))+0])
#define		iy(i)		(poly[2*((i)%(n))+1])
#define		inhalf(i)	((ix(i)-px)*nx+(iy(i)-py)*ny)
#define		flag(i)		flags[(i)%(n)]

/*****************************************************************************/

double polygon_area_signed(double *poly,int n)
{
 double	mx,my,a;
 int	i;

 if ( n <= 2 )
	return(0.0);

 mx=my=0;
 for ( i=0 ; i<n ; i++ )
  {	mx+=ix(i);
	my+=iy(i);
  }
 mx /= (double)n;
 my /= (double)n;

 a=0.0;
 for ( i=0 ; i<n ; i++ )
  {	a += (ix(i)-mx)*(iy(i+1)-my)-(ix(i+1)-mx)*(iy(i)-my);	
  }

 return(a/2.0);
}

double polygon_area(double *poly,int n)
{
 double	a;
 a=polygon_area_signed(poly,n);
 return(0.0<=a?a:-a);
}

int polygon_integrate_linear_monoms(double *poly,int n,double *fcoeff)
{
 double	mx,my;
 int	i;

 fcoeff[0]=0.0; 
 fcoeff[1]=0.0;
 fcoeff[2]=0.0;

 if ( n <= 2 )
	return(0);

 mx=my=0;

 for ( i=0 ; i<n ; i++ )
  {	fcoeff[0] += (ix(i)-mx)*(iy(i+1)-my)-(ix(i+1)-mx)*(iy(i)-my)/2.0;
	fcoeff[1] += +(iy(i+1)-iy(i))*((ix(i)-mx)*(ix(i)-mx)+(ix(i)-mx)*(ix(i+1)-mx)+(ix(i+1)-mx)*(ix(i+1)-mx))/3.0;
	fcoeff[2] += -(ix(i+1)-ix(i))*((iy(i)-my)*(iy(i)-my)+(iy(i)-my)*(iy(i+1)-my)+(iy(i+1)-my)*(iy(i+1)-my))/3.0;

  }

 return(0);
}

int polygon_convexity(double *poly,int n)
{
 int	i,np,nn,c;
 double	dx1,dx2,dy1,dy2;

 np=nn=0;

 for ( i=0 ; i<n ; i++ )
  {	dx2=ix(i+2)-ix(i+1);
	dy2=iy(i+2)-iy(i+1);
	dx1=ix(i+1)-ix(i);
	dy1=iy(i+1)-iy(i);
	if ( dx1*dy2 - dx2*dy1 >= 0.0 )	np++;
	else				nn++;
  }
 c=(np<nn?np:nn);

 return(c);
}

static int polygon_line_intersection(double x0,double y0,double x1,double y1,double px,double py,double nx,double ny,double *rx,double *ry)
{
 double	a;

 a=((px-x1)*nx+(py-y1)*ny)/((x0-x1)*nx+(y0-y1)*ny);
 *rx=x0*a+x1*(1-a);
 *ry=y0*a+y1*(1-a);

 return(0);
}

int polygon_intersection_halfplane(double *poly,int n,double px,double py,double nx,double ny)
{
 int		i,j,c;
 unsigned char	*flags;

 c=0;
 for ( i=0 ; i<n ; i++ )
  {	if ( 0.0 <= inhalf(i) )
		c++;
  }

 if ( c==0 )
	return(0);
 else if ( c==n )
	return(n);

 flags=(unsigned char *)malloc(2*n);
 for ( i=0 ; i<n ; i++ )
  {	if ( 0.0 <= inhalf(i) )
		flags[i]=1;
	else
		flags[i]=0;
  }

 for ( i=0 ; i<n ; )
  {	
	if ( flag(i) && ! flag(i+1) )
	 {	int	p0,p1,d,l;
		double	xn,yn,xp,yp;
		p0=i;
		for ( j=i+2 ; j<i+2+n ; j++ )
		 {	if ( flag(j) )	break;		}
		p1=j;

		polygon_line_intersection(ix(p0),iy(p0),ix(p0+1),iy(p0+1),px,py,nx,ny,&xn,&yn);
		polygon_line_intersection(ix(p1+n-1),iy(p1+n-1),ix(p1),iy(p1),px,py,nx,ny,&xp,&yp);

		if ( p1<n )
		 {	if ( p1-p0 < 3 )	/* equivalent to p1-p0 == 2, add a point: */
			 {	memmove(poly+2*p1+2,poly+2*p1,sizeof(double)*2*(n-p1));
				memmove(flags+p1+1,flags+p1,n-p1);
				p1++;
				n++;
			 }
			else if ( 3 < p1-p0 )	/* remove p1-p0-3 points: */
			 {	d=p1-p0-3;
				memmove(poly+2*(p0+3),poly+2*p1,sizeof(double)*2*(n-p1));
				memmove(flags+p0+3,flags+p1,n-p1);
				p1 -= d;
				n  -= d;
			 }
		 }
		else
		 {	if ( p1-p0 < 3 )	/* equivalent to p1-p0 == 2, add a point: */
				n++;
			else if ( 3 < p1-p0 )	/* remove p1-p0-3 points: */
			 {	d=p1-p0-3;
				l=n-p0-1;
				if ( d <= l )
					n -= d;
				else
				 {	d -= l;
					n -= l;
					memmove(poly,poly+2*d,sizeof(double)*2*(n-d));
					memmove(flags,flags+d,n-d);
					n -= d;
					p0 -= d;
					p1 -= d;
				 }
			 }
		 }

		ix(p0+1)=xn;
		iy(p0+1)=yn;
		flag(p0+1)=1;	/* new points are always inside */
		ix(p0+2)=xp;
		iy(p0+2)=yp;
		flag(p0+2)=1;

		i=0;	/* start over */
	 }
	else
		i++;	/* proceed */
  }

 free(flags);

 return(n);
}

int polygon_intersection_square(double *poly,int n,double x0,double y0,double dx,double dy)
{
 if ( 0<n )	n=polygon_intersection_halfplane(poly,n,x0+dx/2,y0,0,1);
 if ( 0<n )	n=polygon_intersection_halfplane(poly,n,x0+dx,y0+dy/2,-1,0);
 if ( 0<n )	n=polygon_intersection_halfplane(poly,n,x0+dx/2,y0+dy,0,-1);
 if ( 0<n )	n=polygon_intersection_halfplane(poly,n,x0,y0+dy/2,1,0);
 return(n);
}

/*
int main(int argc,char *argv[])
{
 double	*poly,*poly2,a1,a2;
 int	i,j,n,s;

 poly=(double *)malloc(sizeof(double)*100);
 poly2=(double *)malloc(sizeof(double)*100);
 n=6;
 ix(0)=0.77323540;iy(0)=6.82546684;
 ix(1)=2.08656504;iy(1)=3.51133654;
 ix(2)=6.84659813;iy(2)=1.27056651;
 ix(3)=8.63656006;iy(3)=7.24079234;
 ix(4)=5.05650154;iy(4)=4.64808897;
 ix(5)=3.43235789;iy(5)=8.70133549; 

 a1=polygon_area(poly,n);
 a2=0.0;

 for ( i=0 ; i<10 ; i+=2 )
  {	for ( j=0 ; j<10 ; j+=2 )
	 {	memmove(poly2,poly,sizeof(double)*2*n);
		s=polygon_intersection_square(poly2,n,j,i,2,2);
		a2 += polygon_area(poly2,s);
	 }
  }
 fprintf(stdout,"%12g %12g %12g\n",a1,a2,a1-a2);
 
 return(0);
}
*/

/*****************************************************************************/
                               
