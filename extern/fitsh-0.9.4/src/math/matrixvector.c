/*****************************************************************************/
/* matrixvector.c							     */
/*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

#include "matrixvector.h"

int matrix_mul(matrix a,matrix b,matrix c)
{
 matrix r;
 int	i,j,k;
 for ( i=0 ; i<3 ; i++ )
  {	for ( j=0 ; j<3 ; j++ )
	 {	r[i][j]=0.0;
		for ( k=0 ; k<3 ; k++ )
		 {	r[i][j] += b[i][k]*c[k][j];	}
	 }
  }
 for ( i=0 ; i<3 ; i++ )
  {	for ( j=0 ; j<3 ; j++ )
	 {	a[i][j]=r[i][j];		}
  }
 return(0);
}

int vector_mul(vector a,matrix b,vector c)
{
 vector r;
 int	i,k;
 for ( i=0 ; i<3 ; i++ )
  {	r[i]=0.0;
	for ( k=0 ; k<3 ; k++ )
	 {	r[i] += b[i][k]*c[k];	}
  }
 for ( i=0 ; i<3 ; i++ )
  {	a[i]=r[i];	}
 return(0);
}

double vector_matrix_vector_mul(vector a,matrix b,vector c)
{
 double	r;
 int	i,k;
 r=0;
 for ( i=0 ; i<3 ; i++ )
  {	for ( k=0 ; k<3 ; k++ )
	 {	r += a[i]*b[i][k]*c[k];		}
  }
 return(r);
}

double vector_vector_mul(vector a,vector c)
{
 double	r;
 int	i;
 r=0;
 for ( i=0 ; i<3 ; i++ )
  {	r += a[i]*c[i];		}
 return(r);
}

int matrix_add(matrix a,matrix b,matrix c)
{
 int	i,j;
 for ( i=0 ; i<3 ; i++ )
  {	for ( j=0 ; j<3 ; j++ )
	 {	a[i][j]=b[i][j]+c[i][j];		}
  }
 return(0);
}

int matrix_diad(matrix a,vector b,vector c)
{
 int	i,j;
 for ( i=0 ; i<3 ; i++ )
  {	for ( j=0 ; j<3 ; j++ )
	 {	a[i][j]=b[i]*c[j];		}
  }
 return(0);
}

int matrix_sub(matrix a,matrix b,matrix c)
{
 int	i,j;
 for ( i=0 ; i<3 ; i++ )
  {	for ( j=0 ; j<3 ; j++ )
	 {	a[i][j]=b[i][j]-c[i][j];		}
  }
 return(0);
}

int matrix_copy(matrix a,matrix b)
{
 int	i,j;
 for ( i=0 ; i<3 ; i++ )
  {	for ( j=0 ; j<3 ; j++ )
	 {	a[i][j]=b[i][j];		}
  }
 return(0);
}

int matrix_real_mul(matrix a,matrix b,double d)
{
 int	i,j;
 for ( i=0 ; i<3 ; i++ )
  {	for ( j=0 ; j<3 ; j++ )
	 {	a[i][j]=b[i][j]*d;		}
  }
 return(0);
}

int matrix_unity(matrix a)
{
 int	i,j;
 for ( i=0 ; i<3 ; i++ )
  {	for ( j=0 ; j<3 ; j++ )
	 {	a[i][j]=(i==j?1.0:0.0);		}
  }
 return(0);
}

int matrix_zero(matrix a)
{
 int	i,j;
 for ( i=0 ; i<3 ; i++ )
  {	for ( j=0 ; j<3 ; j++ )
	 {	a[i][j]=0.0;		}
  }
 return(0);
}

int matrix_antisym(matrix m,double a,double b,double c)
{
 m[0][0]=m[1][1]=m[2][2]=0.0;
 m[0][1]=-c;
 m[0][2]=+b;
 m[1][0]=+c;
 m[1][2]=-a;
 m[2][0]=-b;
 m[2][1]=+a;
 return(0);
}

int matrix_is_diff(matrix a,matrix b)
{
 int	i,j;
 for ( i=0 ; i<3 ; i++ )
  {	for ( j=0 ; j<3 ; j++ )
	 {	if ( a[i][j] != b[i][j] )
			return(1);
	 }
  }
 return(0);
}

int matrix_exp(matrix e,matrix a)
{
 int	t;
 matrix	w,s;

 matrix_unity(e);
 matrix_copy(w,e);

 for ( t=1 ; t<=100 ; t++ )
  {	matrix_mul(w,w,a);
	matrix_real_mul(w,w,1.0/(double)t);
	matrix_add(s,e,w);
	if ( ! matrix_is_diff(s,e) )
		break;
	matrix_copy(e,s);
  }
/* fprintf(stderr,"t=%d\n",t); */
 return(0);
}

int matrix_rot(matrix o,double a,double b,double c)
{
 double	n,a0,b0,c0;
 vector	v0;
 matrix	vv,m1;
 n=sqrt(a*a+b*b+c*c);
 v0[0]=a0=a/n;
 v0[1]=b0=b/n;
 v0[2]=c0=c/n;

 matrix_diad(vv,v0,v0);
 matrix_unity(o);
 matrix_sub(o,o,vv);
 matrix_real_mul(o,o,cos(n));
 matrix_antisym(m1,a0,b0,c0);
 matrix_real_mul(m1,m1,sin(n));
 matrix_add(o,o,m1);
 matrix_add(o,o,vv);

 return(0);
}

static double hypot3d(double x,double y,double z)
{
 return(sqrt(x*x+y*y+z*z));
}

double vector_length(vector v)
{
 double	r;
 r=hypot3d(v[0],v[1],v[2]);
 return(r);
}

int vector_normalize(vector v)
{
 double r;
 r=vector_length(v);
 if ( 0.0<r )
  {	v[0]/=r;
	v[1]/=r;
	v[2]/=r;
  }
 return(0);
}

/*****************************************************************************/
                                                          

