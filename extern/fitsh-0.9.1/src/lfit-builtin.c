/*****************************************************************************/
/* lfit-builtin.c							     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* LFIT - Symbolic fitting & arithmetic evaluating utility. 		     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* This file defines the built-in functions available in LFIT by default.    */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 1996, 2002, 2004-2005, 2006, 2007-2008; Pal, A. (apal@szofi.elte.hu)  */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#if defined     _CHDET_SOURCE
	#include "psn.h"
#elif defined   _FI_SOURCE
	#include <psn/psn.h>
#else
	#include <psn/psn.h>
#endif

#if defined   _FI_SOURCE
        #include "math/elliptic/elliptic.h"
        #include "math/elliptic/ntiq.h"
	#include "math/spline/spline.h"
	#include "xfunct.h"
#elif defined	_ASTRO_EXTEN
        #include "elliptic/elliptic.h"
        #include "elliptic/ntiq.h"
	#include "spline/spline.h"
	#include "xfunct.h"
#endif

#include <lfit/lfit.h>

#include "lfit-builtin.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841968
#endif

#ifndef	M_LN10
#define	M_LN10 2.3025850929940456840179914546843642076011
#endif

#ifndef	M_LOG10E
#define	M_LOG10E 0.4342944819032518276511289189166050822943
#endif

/*****************************************************************************/

static int bltf_o_add(double *s) { *(s-2)=(*(s-2))+(*(s-1));return(0);	}
static int bltf_o_sub(double *s) { *(s-2)=(*(s-2))-(*(s-1));return(0);	}
static int bltf_o_mul(double *s) { *(s-2)=(*(s-2))*(*(s-1));return(0);	}
static int bltf_o_div(double *s) { if ( *(s-1)==0.0 ) return(1);
				   *(s-2)=(*(s-2))/(*(s-1));return(0);	}
static int bltf_o_pow(double *s) { *(s-2)=pow(*(s-2),*(s-1));
				   return(0);				}

static int bltf_o_chs(double *s) { *(s-1)=-(*(s-1));return(0);		}

static int bltf_o_psq(double *s) { s--,(*s)*=*s;return(0);		}
static int bltf_o_rcp(double *s) { s--;if ( *s==0.0 ) return(1);
				   *s=1.0/(*s);return(0);		}

static int bltf_o_sqr(double *s) { s--;if ( *s<0.0 ) return(1);
				   *s=sqrt(*s);return(0);		}
static int bltf_o_abs(double *s) { s--;*s=fabs(*s);return(0);		}
static int bltf_o_sgn(double *s) { s--;	if ( *s>0 ) *s=1.0;
					else if ( *s<0 ) *s=-1.0;
					else *s=0.0;
				   return(0);				}
static int bltf_o_theta(double *s) { s--;if ( *s>0 ) *s=1.0;
					else if ( *s<0 ) *s=0.0;
					else *s=0.5;
				   return(0);				}
static int bltf_o_fmod(double *s) { double m,d;int k;
				   s-=2;m=*s,d=*(s+1);
				   if ( d<0.0 ) d=-d; else if ( d==0.0 ) return(1);
				   if ( m<0.0 ) k=(int)((-m)/d),m+=d*(double)(k+2);
				   k=(int)(m/d),m-=(double)k*d;
				   *s=m;return(0);			}
static int bltf_o_fint(double *s) { s--;*s=floor(*s);return(0);		}

static int bltf_o_fdiv(double *s) { double m,d;
				   s-=2;m=*s,d=*(s+1);
				   if ( d<0.0 ) d=-d; else if ( d==0.0 ) return(1);
				   m=floor(m/d);
				   *s=m;return(0);			}

static int bltf_o_vpi(double *s) { *s=M_PI;return(0); }

static int bltf_f_arg(double *s)   { s-=2; *s=atan2(*(s+1),*(s));return(0); }
static int bltf_f_atan2(double *s) { s-=2; *s=atan2(*(s),*(s+1));return(0); }
static int bltf_f_hypot(double *s) { s-=2; *s=sqrt((*s)*(*s)+(*(s+1))*(*(s+1)));return(0); }

static int bltf_f_exp(double *s) { s--;*s=exp(*s);return(0);		}
static int bltf_f_log(double *s) { s--;if ( *s<=0.0 ) return(1);
				   *s=log(*s);return(0);		}
static int bltf_f_exp10(double *s) { s--;*s=exp(*s*(M_LN10));return(0);	}
static int bltf_f_log10(double *s) { s--;if ( *s<=0.0 ) return(1);
				   *s=log(*s)*M_LOG10E;return(0);	}

static int bltf_f_sin(double *s) { *(s-1)=sin(*(s-1));return(0);	}
static int bltf_f_cos(double *s) { *(s-1)=cos(*(s-1));return(0);	}
static int bltf_f_tan(double *s) { *(s-1)=tan(*(s-1));return(0);	}
static int bltf_f_ctn(double *s) { *(s-1)=1.0/tan(*(s-1));return(0);	}

static int bltf_f_sina(double *s) { *(s-1)=sin(*(s-1)*M_PI/180.0);return(0);	 }
static int bltf_f_cosa(double *s) { *(s-1)=cos(*(s-1)*M_PI/180.0);return(0);	 }
static int bltf_f_tana(double *s) { *(s-1)=tan(*(s-1)*M_PI/180.0);return(0);	 }
static int bltf_f_ctna(double *s) { *(s-1)=1.0/tan(*(s-1)*M_PI/180.0);return(0); }

static int bltf_f_asin(double *s) { *(s-1)=asin(*(s-1));return(0); }
static int bltf_f_acos(double *s) { *(s-1)=acos(*(s-1));return(0); }
static int bltf_f_atan(double *s) { *(s-1)=atan(*(s-1));return(0); }
static int bltf_f_actn(double *s) { *(s-1)=atan(1.0/(*(s-1)));return(0); }

static int bltf_f_asina(double *s) { *(s-1)=asin(*(s-1))*180.0/M_PI;return(0); }
static int bltf_f_acosa(double *s) { *(s-1)=acos(*(s-1))*180.0/M_PI;return(0); }
static int bltf_f_atana(double *s) { *(s-1)=atan(*(s-1))*180.0/M_PI;return(0); }
static int bltf_f_actna(double *s) { *(s-1)=atan(1.0/(*(s-1)))*180.0/M_PI;return(0); }

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static int bltf_f_ilinear(double *s)
{
 int	n,b;
 double	x,r;

 s-=3;
 n=(int)s[0];
 b=(int)s[1];
 x=s[2];
 if ( x <= b-1 )
	r=0.0;
 else if ( x <= b )
	r=1.0-(b-x);
 else if ( x < b+1 )
	r=1.0-(x-b);
 else	
	r=0.0;
 *s=r;
 return(0);
}

static int bltf_f_ilinear_dx(double *s)
{
 int	n,b;
 double	x,r;

 s-=3;
 n=(int)s[0];
 b=(int)s[1];
 x=s[2];
 if ( x <= b-1 )
	r=0.0;
 else if ( x <= b )
	r=+1.0;
 else if ( x < b+1 )
	r=-1.0;
 else	
	r=0.0;
 *s=r;
 return(0);
}

static int bltf_f_ispline(double *s)
{
 static	int	t_save=0;
 static double	*arr=NULL;
 int	n,t,b;
 double	x,r;

 s-=3;
 n=(int)s[0];
 b=(int)s[1];
 x=s[2];
 if ( n<0 || ! ( 0<=b && b<=n ) )
	return(-1);
 if ( x<0 )		x=0.0;
 else if ( x>n )	x=(double)n;

 t=n+1;
 if ( t_save != t )
  {	int	i,j;
	t_save=t;
	if ( arr != NULL )
	 {	free(arr);	
		arr=NULL;
	 }
	arr=(double *)malloc(sizeof(double)*(2*t*t));
	for ( i=0 ; i<t ; i++ )
	 {	for ( j=0 ; j<t ; j++ )
		 {	arr[i*t+j]=(i==j?1.0:0.0);		}
		natspline_coeff(&arr[i*t],t,&arr[(i+t)*t]);
	 }
  }
 r=natspline_inter(&arr[b*t],&arr[(b+t)*t],t,x);

 s[0]=r;

 return(0);
}

static int bltf_f_ispline_dx(double *s)
{
 static	int	t_save=0;
 static double	*arr=NULL;
 int	n,t,b;
 double	x,r;

 s-=3;
 n=(int)s[0];
 b=(int)s[1];
 x=s[2];
 if ( n<0 || ! ( 0<=b && b<=n ) )
	return(-1);
 if ( x<0 )		x=0.0;
 else if ( x>n )	x=(double)n;

 t=n+1;
 if ( t_save != t )
  {	int	i,j;
	t_save=t;
	if ( arr != NULL )
	 {	free(arr);	
		arr=NULL;
	 }
	arr=(double *)malloc(sizeof(double)*(2*t*t));
	for ( i=0 ; i<t ; i++ )
	 {	for ( j=0 ; j<t ; j++ )
		 {	arr[i*t+j]=(i==j?1.0:0.0);		}
		natspline_coeff(&arr[i*t],t,&arr[(i+t)*t]);
	 }
  }
 r=natspline_inter_der(&arr[b*t],&arr[(b+t)*t],t,x);

 s[0]=r;

 return(0);
}

static int bltf_f_icyspline(double *s)
{
 static	int	n_save=0;
 static double	*arr=NULL;
 int	n;
 double	x,r;

 s-=2;
 n=(int)s[0];
 x=s[1];
 if ( n<0 )
	return(-1);

 if ( n_save != n )
  {	int	i;
	n_save=n;
	if ( arr != NULL )
	 {	free(arr);	
		arr=NULL;
	 }
	arr=(double *)malloc(sizeof(double)*(2*n));
	for ( i=0 ; i<n ; i++ )
	 {	arr[i]=0.0;		}
	arr[0]=1.0;
	cyspline_coeff(arr,n,arr+n);
  }
 r=cyspline_inter(arr,arr+n,n,x);

 s[0]=r;

 return(0);
}

static int bltf_f_icyspline_dx(double *s)
{
 static	int	n_save=0;
 static double	*arr=NULL;
 int	n;
 double	x,r;

 s-=2;
 n=(int)s[0];
 x=s[1];
 if ( n<0 )
	return(-1);

 if ( n_save != n )
  {	int	i;
	n_save=n;
	if ( arr != NULL )
	 {	free(arr);	
		arr=NULL;
	 }
	arr=(double *)malloc(sizeof(double)*(2*n));
	for ( i=0 ; i<n ; i++ )
	 {	arr[i]=0.0;		}
	arr[0]=1.0;
	cyspline_coeff(arr,n,arr+n);
  }
 r=cyspline_inter_der(arr,arr+n,n,x);

 s[0]=r;

 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static int bltf_f_jbessel(double *s)
{
 double	x,y;
 int	n;

 s-=2;
 n=(int)floor(s[0]);
 x=s[1];

 if ( n==0 )		y=j0(x);
 else if ( n==1 )	y=j1(x);
 else if ( n==-1 )	y=-j1(x);
 else			y=jn(n,x);

 s[0]=y;

 return(0);
}

static int bltf_f_ybessel(double *s)
{
 double	x,y;
 int	n;

 s-=2;
 n=(int)floor(s[0]);
 x=s[1];

 if ( x<=0.0 )		return(1);

 if ( n==0 )		y=y0(x);
 else if ( n==1 )	y=y1(x);
 else if ( n==-1 )	y=-y1(x);
 else			y=yn(n,x);

 s[0]=y;

 return(0);
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* Definitions of additional operators: relations and boolean operators: */

static int bltf_o_equ(double *s) { if ( *(s-1)==*(s-2) ) *(s-2)=1.0; else *(s-2)=0.0; return(0); }
static int bltf_o_neq(double *s) { if ( *(s-1)!=*(s-2) ) *(s-2)=1.0; else *(s-2)=0.0; return(0); }
static int bltf_o_ls(double *s)  { if ( *(s-2)< *(s-1) ) *(s-2)=1.0; else *(s-2)=0.0; return(0); }
static int bltf_o_gr(double *s)  { if ( *(s-2)> *(s-1) ) *(s-2)=1.0; else *(s-2)=0.0; return(0); }
static int bltf_o_le(double *s)  { if ( *(s-2)<=*(s-1) ) *(s-2)=1.0; else *(s-2)=0.0; return(0); }
static int bltf_o_ge(double *s)  { if ( *(s-2)>=*(s-1) ) *(s-2)=1.0; else *(s-2)=0.0; return(0); }
static int bltf_o_and(double *s) { if ( *(s-1)!=0 && *(s-2)!=0 ) *(s-2)=1.0; else *(s-2)=0.0; return(0); }
static int bltf_o_or(double *s)  { if ( *(s-1)!=0 || *(s-2)!=0 ) *(s-2)=1.0; else *(s-2)=0.0; return(0); }

/*****************************************************************************/

/* Derivative rules */

static	short	diffrule_o_add[]={ SS_D2,SS_D1,O_ADD,0 };
static	short	diffrule_o_sub[]={ SS_D2,SS_D1,O_SUB,0 };
static	short	diffrule_o_mul[]={ SS_D2,SS_N1,O_MUL,SS_D1,SS_N2,O_MUL,O_ADD,0 };
static	short	diffrule_o_div[]={ SS_D2,SS_N1,O_MUL,SS_D1,SS_N2,O_MUL,O_SUB,SS_N1,O_PSQ,O_DIV,0 };
static	short	diffrule_o_pow[]={ SS_N2,SS_N1,O_POW,SS_D1,O_MUL,SS_N2,F_LOG,O_MUL,SS_N2,SS_N1,CON_1,O_SUB,O_POW,SS_N1,O_MUL,SS_D2,O_MUL,O_ADD,0 };

static	short	diffrule_o_chs[]={ SS_D1,O_CHS,0 };
static	short	diffrule_o_psq[]={ SS_N1,CON_2,O_MUL,SS_D1,O_MUL,0 };
static	short	diffrule_o_rcp[]={ SS_N1,O_PSQ,O_RCP,O_CHS,SS_D1,O_MUL,0 };

static	short	diffrule_f_sqr[]={ SS_N1,F_SQR,CON_2,O_MUL,O_RCP,SS_D1,O_MUL,0 };
static	short	diffrule_f_abs[]={ SS_N1,F_SGN,SS_D1,O_MUL,0 };
static	short	diffrule_f_sgn[]  ={ CON_0,0 };
static	short	diffrule_f_theta[]={ CON_0,0 };
static	short	diffrule_f_fmod[]={ SS_D2,SS_N2,SS_N1,O_DIV,F_FINT,SS_D1,O_MUL,O_SUB,0 };
static	short	diffrule_f_fint[]={ CON_0,0 };
static	short	diffrule_f_fdiv[]={ CON_0,0 };
static	short	diffrule_f_vpi[]={ CON_0,0 };

static	short	diffrule_f_arg  []={ SS_N2,SS_D1,O_MUL,SS_N1,SS_D2,O_MUL,O_SUB,SS_N1,O_PSQ,SS_N2,O_PSQ,O_ADD,O_DIV,0 };
static	short	diffrule_f_atan2[]={ SS_N1,SS_D2,O_MUL,SS_N2,SS_D1,O_MUL,O_SUB,SS_N1,O_PSQ,SS_N2,O_PSQ,O_ADD,O_DIV,0 };
static	short	diffrule_f_hypot[]={ SS_N1,SS_D1,O_MUL,SS_N2,SS_D2,O_MUL,O_ADD,SS_N1,SS_N2,F_HYPOT,O_DIV,0 };

static	short 	diffrule_f_sin[]={ SS_N1,F_COS,SS_D1,O_MUL,0 };
static	short	diffrule_f_cos[]={ SS_N1,F_SIN,O_CHS,SS_D1,O_MUL,0 };
static	short	diffrule_f_tan[]={ SS_N1,F_COS,O_PSQ,O_RCP,SS_D1,O_MUL,0 };
static	short	diffrule_f_ctn[]={ SS_N1,F_SIN,O_PSQ,O_RCP,O_CHS,SS_D1,O_MUL,0 };

static	short 	diffrule_f_asin[]={ SS_D1,CON_1,SS_N1,O_PSQ,O_SUB,F_SQR,O_DIV,0 };
static	short	diffrule_f_acos[]={ SS_D1,CON_1,SS_N1,O_PSQ,O_SUB,F_SQR,O_DIV,O_CHS,0 };
static	short	diffrule_f_atan[]={ SS_D1,CON_1,SS_N1,O_PSQ,O_ADD,O_DIV,0 };
static	short	diffrule_f_actn[]={ SS_D1,CON_1,SS_N1,O_PSQ,O_ADD,O_DIV,O_CHS };

static	short	diffrule_f_exp[]={ SS_N1,F_EXP,SS_D1,O_MUL,0 };
static	short	diffrule_f_log[]={ SS_N1,O_RCP,SS_D1,O_MUL,0 };

static	short	diffrule_f_ilinear[]   = { SS_N3,SS_N2,SS_N1,F_ILINEAR_DX,SS_D1,O_MUL,0 };
static	short	diffrule_f_ispline[]   = { SS_N3,SS_N2,SS_N1,F_ISPLINE_DX,SS_D1,O_MUL,0 };
static	short	diffrule_f_icyspline[] = { SS_N2,SS_N1,F_ICYSPLINE_DX,SS_D1,O_MUL,0 };

static	short	diffrule_f_jbessel[] = { SS_N2,CON_1,O_SUB,SS_N1,F_JBESSEL,SS_N2,CON_1,O_ADD,SS_N1,F_JBESSEL,O_SUB,CON_2,O_DIV,SS_D1,O_MUL,0 };
static	short	diffrule_f_ybessel[] = { SS_N2,CON_1,O_SUB,SS_N1,F_YBESSEL,SS_N2,CON_1,O_ADD,SS_N1,F_YBESSEL,O_SUB,CON_2,O_DIV,SS_D1,O_MUL,0 };

/* Simplification rules: */
									/* Type */
		
static	short sim_o_add1[]={ S_EXP1,CON_0,O_ADD,0,S_EXP1,0 };		/* expr+0=expr		*/
static	short sim_o_add2[]={ CON_0,S_EXP1,O_ADD,0,S_EXP1,0 };		/* 0+expr=expr		*/

static	short sim_o_sub1[]={ S_EXP1,CON_0,O_SUB,0,S_EXP1,0 };		/* expr-0=expr		*/
static	short sim_o_sub2[]={ CON_0,S_EXP1,O_SUB,0,S_EXP1,O_CHS,0 };	/* 0-expr=-expr 	*/

static	short sim_o_mul1[]={ S_EXP1,CON_1 ,O_MUL,0,S_EXP1,0 };		/* expr*1=expr		*/
static	short sim_o_mul2[]={ CON_1 ,S_EXP1,O_MUL,0,S_EXP1,0 };		/* 1*expr=expr		*/
static	short sim_o_mul3[]={ S_EXP1,CON_0 ,O_MUL,0,CON_0,0  };		/* expr*0=0		*/
static	short sim_o_mul4[]={ CON_0 ,S_EXP1,O_MUL,0,CON_0,0  };		/* 0*expr=0		*/
static	short sim_o_mul5[]={ S_EXP1,S_EXP1,O_MUL,0,S_EXP1,O_PSQ,0 };	/* ex*ex=ex^(2)		*/

static	short sim_o_div1[]={ S_EXP1,CON_1 ,O_DIV,0,S_EXP1,0 };		/* expr/1=expr      	*/
static	short sim_o_div2[]={ CON_1 ,S_EXP1,O_DIV,0,S_EXP1,O_RCP,0 };	/* 1/expr=expr^(-1) 	*/
static	short sim_o_div3[]={ S_EXP1,CON_0 ,O_DIV,0,0 };			/* expr/0=UNDEF!     	*/
static	short sim_o_div4[]={ CON_0 ,S_EXP1,O_DIV,0,CON_0 ,0 };		/* 0/expr=0       	*/

static	short sim_o_pow1[]={ S_EXP1,CON_0,O_POW,0,CON_1 ,0 };		/* ex^0=1 		*/
static	short sim_o_pow2[]={ S_EXP1,CON_1,O_POW,0,S_EXP1,0 };		/* ex^1=ex		*/
static	short sim_o_pow3[]={ S_EXP1,CON_2,O_POW,0,S_EXP1,O_PSQ,0 };	/* ex^2=ex PSQ		*/

static	short sim_addmul_dis1[]={ S_EXP1,S_EXP2,O_MUL,S_EXP1,S_EXP3,O_MUL,O_ADD,0,S_EXP1,S_EXP2,S_EXP3,O_ADD,O_MUL,0 };	/* a*b+a*c=a*(b+c) */
static	short sim_addmul_dis2[]={ S_EXP2,S_EXP1,O_MUL,S_EXP1,S_EXP3,O_MUL,O_ADD,0,S_EXP1,S_EXP2,S_EXP3,O_ADD,O_MUL,0 }; /* a*b+c*a=a*(b+c) */
static	short sim_addmul_dis3[]={ S_EXP1,S_EXP2,O_MUL,S_EXP3,S_EXP1,O_MUL,O_ADD,0,S_EXP1,S_EXP2,S_EXP3,O_ADD,O_MUL,0 }; /* b*a+a*c=a*(b+c) */
static	short sim_addmul_dis4[]={ S_EXP2,S_EXP1,O_MUL,S_EXP3,S_EXP1,O_MUL,O_ADD,0,S_EXP1,S_EXP2,S_EXP3,O_ADD,O_MUL,0 }; /* b*a+c*a=a*(b+c) */
static	short sim_submul_dis1[]={ S_EXP1,S_EXP2,O_MUL,S_EXP1,S_EXP3,O_MUL,O_SUB,0,S_EXP1,S_EXP2,S_EXP3,O_SUB,O_MUL,0 }; /* a*b-a*c=a*(b-c) */
static	short sim_submul_dis2[]={ S_EXP2,S_EXP1,O_MUL,S_EXP1,S_EXP3,O_MUL,O_SUB,0,S_EXP1,S_EXP2,S_EXP3,O_SUB,O_MUL,0 }; /* a*b-c*a=a*(b-c) */
static	short sim_submul_dis3[]={ S_EXP1,S_EXP2,O_MUL,S_EXP3,S_EXP1,O_MUL,O_SUB,0,S_EXP1,S_EXP2,S_EXP3,O_SUB,O_MUL,0 }; /* b*a-a*c=a*(b-c) */
static	short sim_submul_dis4[]={ S_EXP2,S_EXP1,O_MUL,S_EXP3,S_EXP1,O_MUL,O_SUB,0,S_EXP1,S_EXP2,S_EXP3,O_SUB,O_MUL,0 }; /* b*a-c*a=a*(b-c) */
static	short sim_adddiv_dis []={ S_EXP2,S_EXP1,O_DIV,S_EXP3,S_EXP1,O_DIV,O_ADD,0,S_EXP2,S_EXP3,O_ADD,S_EXP1,O_DIV,0 }; /* b/a+c/a=(b+c)/a */
static	short sim_subdiv_dis []={ S_EXP2,S_EXP1,O_DIV,S_EXP3,S_EXP1,O_DIV,O_SUB,0,S_EXP2,S_EXP3,O_SUB,S_EXP1,O_DIV,0 }; /* b/a-c/a=(b-c)/a */
static	short sim_mulpow_dis []={ S_EXP2,S_EXP1,O_POW,S_EXP3,S_EXP1,O_POW,O_MUL,0,S_EXP2,S_EXP3,O_MUL,S_EXP1,O_POW,0 }; /* b^a*c^a=(b*c)^a */
static	short sim_divpow_dis []={ S_EXP2,S_EXP1,O_POW,S_EXP3,S_EXP1,O_POW,O_DIV,0,S_EXP2,S_EXP3,O_DIV,S_EXP1,O_POW,0 }; /* b^a/c^a=(b/c)^a */

static	short sim_tau1[]= { S_EXP1,S_EXP2,O_CHS,O_SUB,0,S_EXP1,S_EXP2,O_ADD,0 };	/* a-(-b)=a+b */

static  short sim_gon1[]= { S_EXP1,F_SIN,O_PSQ,S_EXP1,F_COS,O_PSQ,O_ADD,0,CON_1,0 };	/* sin^2(x)+cos^2(x)=1 */


short	*psn_lfit_simp[]=
 {	sim_o_add1,
	sim_o_add2,
	sim_o_sub1,
	sim_o_sub2,
	sim_o_mul1,
	sim_o_mul2,
	sim_o_mul3,
	sim_o_mul4,
	sim_o_mul5,
	sim_o_div1,
	sim_o_div2,
	sim_o_div3,
	sim_o_div4,
	sim_o_pow1,
	sim_o_pow2,
	sim_o_pow3,
	sim_addmul_dis1,
	sim_addmul_dis2,
	sim_addmul_dis3,
	sim_addmul_dis4,
	sim_submul_dis1,
	sim_submul_dis2,
	sim_submul_dis3,
	sim_submul_dis4,
	sim_adddiv_dis,
	sim_subdiv_dis,
	sim_mulpow_dis,
	sim_divpow_dis,
	sim_tau1,
	sim_gon1,
	NULL
 };

static psnfunctinfo	pfi_add		= { "addition" 		     		};
static psnfunctinfo	pfi_sub		= { "subtraction" 		     	};
static psnfunctinfo	pfi_chs		= { "negation" 		     		};
static psnfunctinfo	pfi_mul		= { "multiplication" 	     		};
static psnfunctinfo	pfi_div		= { "division"	 	     		};
static psnfunctinfo	pfi_pow		= { "power (right-associative)" 	};

static psnfunctinfo	pfi_sin		= { "sine function (argument in radians)" };
static psnfunctinfo	pfi_cos		= { "cosine function (argument in radians)" };
static psnfunctinfo	pfi_tan		= { "tangent function (argument in radians)" };
static psnfunctinfo	pfi_cot		= { "cotangent function (argument in radians)" };
static psnfunctinfo	pfi_asin	= { "inverse sine function (results radians)" };
static psnfunctinfo	pfi_acos	= { "inverse cosine function (results radians)" };
static psnfunctinfo	pfi_atan	= { "inverse tangent function (results radians)" };
static psnfunctinfo	pfi_acot	= { "inverse cotangent function (results radians)" };

static psnfunctinfo	pfi_exp 	= { "exponential function (natural, e-based)" };
static psnfunctinfo	pfi_log		= { "natural logarithm" };
static psnfunctinfo	pfi_exp10	= { "exponential function (base 10)" };
static psnfunctinfo	pfi_log10	= { "logarithm to the base 10" };

static psnfunctinfo	pfi_sgn		= { "sign function" };
static psnfunctinfo	pfi_theta	= { "Heaviside step function" };
static psnfunctinfo	pfi_abs		= { "absolute value function" };
static psnfunctinfo	pfi_sqrt	= { "square root function" };
static psnfunctinfo	pfi_fint	= { "integer part function (i.e. results the largest integer which is not larger than the argument)" };
static psnfunctinfo	pfi_fdiv	= { "integer division, div(x,y) is equivalent to int(x/y)" };
static psnfunctinfo	pfi_fmod	= { "real fractional remainder function" };
static psnfunctinfo	pfi_arg		= { "two dimensional argument function (results radians)" };
static psnfunctinfo	pfi_atan2	= { "two dimensional argument function, as defined in some programming languages, i.e. atan2(y,x)=arg(x,y) (results radians)" };
static psnfunctinfo	pfi_hypot	= { "hypotenuse function, it results the hypotenuse of a right triangle of which catheti are the two arguments of the function" };
static psnfunctinfo	pfi_vpi		= { "results the value of \\pi (a function with no arguments)" };

static psnfunctinfo	pfi_ilinear	= { "base function of linear interpolation" };
static psnfunctinfo	pfi_ispline	= { "base function of cubic spline interpolation" };
static psnfunctinfo	pfi_icyspline	= { "base function of cyclic cubic spline interpolation" };

static psnfunctinfo	pfi_jbessel	= { "Bessel function of the first kind" };
static psnfunctinfo	pfi_ybessel	= { "Bessel function of the second kind" };

/*****************************************************************************/

/*| sym     type  major   minor    | #|  P| assoc      | (*funct)()| diffrule      |  symstring   S  affix    | */
psnlfit	psnlfit_list_builtin_normal_operators[]=										
{ 
  { "+"    , T_OP, O_ADD  , TO_INFIX , 2, 20, ASSOC_LEFT , bltf_o_add  , diffrule_o_add  , "#(1)+#(2)", 0, TO_INFIX , &pfi_add	},
  { "+"    , T_OP, 0      , TO_PREFIX, 0,  0, 0          , NULL        , NULL            , NULL       , 0, 0        , NULL	},
  { "-"    , T_OP, O_SUB  , TO_INFIX , 2, 20, ASSOC_LEFT , bltf_o_sub  , diffrule_o_sub  , "#(1)-#[2]", 0, TO_INFIX , &pfi_sub	},
  { "-"    , T_OP, O_CHS  , TO_PREFIX, 1, 22, 0          , bltf_o_chs  , diffrule_o_chs  , "(-#(1))"  , 1, 0        , &pfi_chs	},
  { "*"    , T_OP, O_MUL  , TO_INFIX , 2, 21, ASSOC_LEFT , bltf_o_mul  , diffrule_o_mul  , "#(1)*#(2)", 0, TO_INFIX , &pfi_mul	},
  { "/"    , T_OP, O_DIV  , TO_INFIX , 2, 21, ASSOC_LEFT , bltf_o_div  , diffrule_o_div  , "#(1)/#[2]", 0, TO_INFIX , &pfi_div	},
  { "^"    , T_OP, O_POW  , TO_INFIX , 2, 23, ASSOC_RIGHT, bltf_o_pow  , diffrule_o_pow  , "#(1)^#(2)", 1, 0        , &pfi_pow	},
  { "RCP"  , T_OP, O_RCP  , 0        , 1, 23, 0          , bltf_o_rcp  , diffrule_o_rcp  , "1.0/#[1]" , 0, TO_INFIX , NULL	},
  { "PSQ"  , T_OP, O_PSQ  , 0        , 1, 23, 0          , bltf_o_psq  , diffrule_o_psq  , "#(1)^2"   , 1, 0        , NULL	}, 
  { NULL   , 0   , 0      , 0        , 0,  0, 0          , NULL        , NULL            , NULL       , 0, 0        , NULL	} 
};															
															
psnlfit	psnlfit_list_builtin_relational_operators[]=									
{															
  { "=="   , T_OP, O_EQU  , TO_INFIX , 2, 10, ASSOC_LEFT , bltf_o_equ  , NULL            , NULL          , 0, 0     , NULL	},
  { "!="   , T_OP, O_NEQ  , TO_INFIX , 2, 10, ASSOC_LEFT , bltf_o_neq  , NULL            , NULL          , 0, 0     , NULL	},
  { "<"    , T_OP, O_LS   , TO_INFIX , 2, 10, ASSOC_LEFT , bltf_o_ls   , NULL            , NULL          , 0, 0     , NULL	},
  { ">"    , T_OP, O_GR   , TO_INFIX , 2, 10, ASSOC_LEFT , bltf_o_gr   , NULL            , NULL          , 0, 0     , NULL	},
  { "<="   , T_OP, O_LE   , TO_INFIX , 2, 10, ASSOC_LEFT , bltf_o_le   , NULL            , NULL          , 0, 0     , NULL	},
  { ">="   , T_OP, O_GE   , TO_INFIX , 2, 10, ASSOC_LEFT , bltf_o_ge   , NULL            , NULL          , 0, 0     , NULL	},
  { "&&"   , T_OP, O_AND  , TO_INFIX , 2,  5, ASSOC_LEFT , bltf_o_and  , NULL            , NULL          , 0, 0     , NULL	},
  { "||"   , T_OP, O_OR   , TO_INFIX , 2,  4, ASSOC_LEFT , bltf_o_or   , NULL            , NULL          , 0, 0     , NULL	},
  { NULL   , 0   , 0      , 0        , 0,  0, 0          , NULL        , NULL            , NULL          , 0, 0     , NULL	},
};															
															
psnlfit psnlfit_list_builtin_elementary_functions[]=								
{															
  { "sin"  , T_FN, F_SIN  , 0        , 1,  0, 0          , bltf_f_sin  , diffrule_f_sin  , "sin(#1)"     , 0, 0     , &pfi_sin	},
  { "cos"  , T_FN, F_COS  , 0        , 1,  0, 0          , bltf_f_cos  , diffrule_f_cos  , "cos(#1)"     , 0, 0     , &pfi_cos	},
  { "tan"  , T_FN, F_TAN  , 0        , 1,  0, 0          , bltf_f_tan  , diffrule_f_tan  , "tan(#1)"     , 0, 0     , &pfi_tan	},
  { "cot"  , T_FN, F_CTN  , 0        , 1,  0, 0          , bltf_f_ctn  , diffrule_f_ctn  , "cot(#1)"     , 0, 0     , &pfi_cot	},
  { "asin" , T_FN, F_ASIN , 0        , 1,  0, 0          , bltf_f_asin , diffrule_f_asin , "asin(#1)"    , 0, 0     , &pfi_asin	},
  { "acos" , T_FN, F_ACOS , 0        , 1,  0, 0          , bltf_f_acos , diffrule_f_acos , "acos(#1)"    , 0, 0     , &pfi_acos	},
  { "atan" , T_FN, F_ATAN , 0        , 1,  0, 0          , bltf_f_atan , diffrule_f_atan , "atan(#1)"    , 0, 0     , &pfi_atan	},
  { "acot" , T_FN, F_ACTN , 0        , 1,  0, 0          , bltf_f_actn , diffrule_f_actn , "acot(#1)"    , 0, 0     , &pfi_acot	},
  { "sina" , T_FN, F_SINA , 0        , 1,  0, 0          , bltf_f_sina , NULL            , "sina(#1)"    , 0, 0     , NULL	},
  { "cosa" , T_FN, F_COSA , 0        , 1,  0, 0          , bltf_f_cosa , NULL            , "cosa(#1)"    , 0, 0     , NULL	},
  { "tana" , T_FN, F_TANA , 0        , 1,  0, 0          , bltf_f_tana , NULL            , "tana(#1)"    , 0, 0     , NULL	},
  { "ctga" , T_FN, F_CTNA , 0        , 1,  0, 0          , bltf_f_ctna , NULL            , "ctna(#1)"    , 0, 0     , NULL	},
  { "asina", T_FN, F_ASINA, 0        , 1,  0, 0          , bltf_f_asina, NULL            , "asina(#1)"   , 0, 0     , NULL	},
  { "acosa", T_FN, F_ACOSA, 0        , 1,  0, 0          , bltf_f_acosa, NULL            , "acosa(#1)"   , 0, 0     , NULL	},
  { "atana", T_FN, F_ATANA, 0        , 1,  0, 0          , bltf_f_atana, NULL            , "atana(#1)"   , 0, 0     , NULL	},
  { "actga", T_FN, F_ACTNA, 0        , 1,  0, 0          , bltf_f_actna, NULL            , "actna(#1)"   , 0, 0     , NULL	},
															
  { "exp"  , T_FN, F_EXP  , 0        , 1,  0, 0          , bltf_f_exp  , diffrule_f_exp  , "exp(#1)"     , 0, 0     , &pfi_exp	},
  { "log"  , T_FN, F_LOG  , 0        , 1,  0, 0          , bltf_f_log  , diffrule_f_log  , "log(#1)"     , 0, 0     , &pfi_log	},
  { "exp10", T_FN, F_EXP10, 0        , 1,  0, 0          , bltf_f_exp10, NULL            , "exp10(#1)"   , 0, 0     , &pfi_exp10},
  { "log10", T_FN, F_LOG10, 0        , 1,  0, 0          , bltf_f_log10, NULL            , "log10(#1)"   , 0, 0     , &pfi_log10},
  { "sqrt" , T_FN, F_SQR  , 0        , 1,  0, 0          , bltf_o_sqr  , diffrule_f_sqr  , "sqrt(#1)"    , 0, 0     , &pfi_sqrt	},
  { "abs"  , T_FN, F_ABS  , 0        , 1,  0, 0          , bltf_o_abs  , diffrule_f_abs  , "abs(#1)"     , 0, 0     , &pfi_abs	},
  { "sign" , T_FN, F_SGN  , 0        , 1,  0, 0          , bltf_o_sgn  , diffrule_f_sgn  , "sign(#1)"    , 0, 0     , &pfi_sgn	},
  { "theta", T_FN, F_THETA, 0        , 1,  0, 0          , bltf_o_theta, diffrule_f_theta, "theta(#1)"   , 0, 0     , &pfi_theta},
  { "int"  , T_FN, F_FINT , 0        , 1,  0, 0          , bltf_o_fint , diffrule_f_fint , "int(#1)"     , 0, 0     , &pfi_fint	},
  { "div"  , T_FN, F_FDIV , 0        , 2,  0, 0          , bltf_o_fdiv , diffrule_f_fdiv , "div(#1,#2)"  , 0, 0     , &pfi_fdiv	},
  { "mod"  , T_FN, F_FMOD , 0        , 2,  0, 0          , bltf_o_fmod , diffrule_f_fmod , "mod(#1,#2)"  , 0, 0     , &pfi_fmod	},
  { "arg"  , T_FN, F_ARG  , 0        , 2,  0, 0          , bltf_f_arg  , diffrule_f_arg  , "arg(#1,#2)"  , 0, 0     , &pfi_arg	},
  { "atan2", T_FN, F_ATAN2, 0        , 2,  0, 0          , bltf_f_atan2, diffrule_f_atan2, "atan2(#1,#2)", 0, 0     , &pfi_atan2},
  { "hypot", T_FN, F_HYPOT, 0        , 2,  0, 0          , bltf_f_hypot, diffrule_f_hypot, "hypot(#1,#2)", 0, 0     , &pfi_hypot},
  { "pi"   , T_FN, F_VPI  , 0        , 0,  0, 0          , bltf_o_vpi  , diffrule_f_vpi  , "pi()"        , 0, 0     , &pfi_vpi	},
  { NULL   , 0   , 0      , 0        , 0,  0, 0          , NULL        , NULL            , NULL          , 0, 0     , NULL	},
};																	

psnlfit psnlfit_list_builtin_interpolators[]=								
{															
  { "ilinear"     , T_FN, F_ILINEAR     , 0        , 3,  0, 0          , bltf_f_ilinear     , diffrule_f_ilinear  , "ilinear(#1,#2,#3)"     , 0, 0  , &pfi_ilinear	},
  { "ilinear_dx"  , T_FN, F_ILINEAR_DX  , 0        , 3,  0, 0          , bltf_f_ilinear_dx  , NULL                , "ilinear_dx(#1,#2,#3)"  , 0, 0  , NULL 		},
  { "ispline"     , T_FN, F_ISPLINE     , 0        , 3,  0, 0          , bltf_f_ispline     , diffrule_f_ispline  , "ispline(#1,#2,#3)"     , 0, 0  , &pfi_ispline	},
  { "ispline_dx"  , T_FN, F_ISPLINE_DX  , 0        , 3,  0, 0          , bltf_f_ispline_dx  , NULL                , "ispline_dx(#1,#2,#3)"  , 0, 0  , NULL 		},
  { "icyspline"   , T_FN, F_ICYSPLINE   , 0        , 2,  0, 0          , bltf_f_icyspline   , diffrule_f_icyspline, "icyspline(#1,#2)"      , 0, 0  , &pfi_icyspline	},
  { "icyspline_dx", T_FN, F_ICYSPLINE_DX, 0        , 2,  0, 0          , bltf_f_icyspline_dx, NULL                , "icyspline_dx(#1,#2)"   , 0, 0  , NULL 		},
  { NULL          , 0   , 0             , 0        , 0,  0, 0          , NULL               , NULL                , NULL                    , 0, 0  , NULL		},
};

psnlfit psnlfit_list_builtin_aa_functions[]=
{
  { "jbessel",	T_FN, F_JBESSEL, 0, 2, 0, 0, bltf_f_jbessel, diffrule_f_jbessel	, "jbessel(#1,#2)", 0, 0  , &pfi_jbessel	},
  { "ybessel",	T_FN, F_YBESSEL, 0, 2, 0, 0, bltf_f_ybessel, diffrule_f_ybessel	, "ybessel(#1,#2)", 0, 0  , &pfi_ybessel	},
  { NULL     ,  0   , 0        , 0, 0, 0, 0, NULL	   , NULL		, NULL		  , 0, 0  , NULL		}
};


/*****************************************************************************/
/* If 'lfit' is the part of the 'fi' package, we register some extra	     */
/* functions, declared and defined in xfunct.[ch]. 			     */
/*****************************************************************************/

#if defined _FI_SOURCE || defined _ASTRO_EXTEN

#define		F_EOQ		96		/* eoq(), eccentric offset   */
#define		F_EOP		97		/* eop(), eccentric offset   */
#define		F_ETC		98		/* etc(), eccentric trig.    */
#define		F_ETS		99		/* ets(), eccentric trig.    */

#define		F_ELL_K		100		/* E() , ell. integral, 1st  */
#define		F_ELL_E		101		/* K() , ell. integral, 2nd  */
#define		F_ELL_PI	102		/* Pi(), ell. integral, 3rd  */
#define		F_NTIU		104		/* ntiu(), norm. trans. int. */
#define		F_NTIQ		105		/* ntiq(), norm. trans. int. */
#define		F_NTIQ_DP	106		/* ntiq_diff(): [0]	     */
#define		F_NTIQ_DZ	107		/* ntiq_diff(): [1]	     */
#define		F_NTIQ_DG1	108		/* ntiq_diff(): [2]	     */
#define		F_NTIQ_DG2	109		/* ntiq_diff(): [3]	     */

#define		F_HJD		110		/* hjd(), helocentric JD     */
#define		F_BJD		111		/* bjd(), barycentric JD     */

static int funct_f_eoq(double *s)
 {	s-=3;
	*s=eccentric_offset_q(*s,*(s+1),*(s+2));
	return(0);
 }

static short    diffrule_f_eoq[]=
 {	SS_N3,SS_N2,SS_N1,F_EOP,SS_N3,O_ADD,F_COS,SS_N2,O_SUB,SS_D2,O_MUL,
	SS_N3,SS_N2,SS_N1,F_EOP,SS_N3,O_ADD,F_SIN,SS_N1,O_SUB,SS_D1,O_MUL,O_ADD,
	SS_N3,SS_N2,SS_N1,F_EOP,SS_D3,O_MUL,O_SUB,
	CON_1,SS_N3,SS_N2,SS_N1,F_EOQ,O_SUB,O_DIV,0 
 };

static int funct_f_eop(double *s) 
 {	s-=3;
	*s=eccentric_offset_p(*s,*(s+1),*(s+2));
	return(0);
 }

static short    diffrule_f_eop[]=
 {	SS_N3,SS_N2,SS_N1,F_EOP,SS_N3,O_ADD,F_SIN,SS_D2,O_MUL,
	SS_N3,SS_N2,SS_N1,F_EOP,SS_N3,O_ADD,F_COS,SS_D1,O_MUL,O_SUB,
	SS_N3,SS_N2,SS_N1,F_EOQ,SS_D3,O_MUL,O_ADD,
	CON_1,SS_N3,SS_N2,SS_N1,F_EOQ,O_SUB,O_DIV,0 
 };

static int funct_f_etc(double *s)
 {	s-=3;
	*s=cos((*s)+eccentric_offset_p(*s,*(s+1),*(s+2)));
	return(0);
 }

static short    diffrule_f_etc[]=
 {	SS_N3,SS_N2,SS_N1,F_ETS,O_CHS,SS_D3,O_MUL,
	SS_N3,SS_N2,SS_N1,F_ETS,O_PSQ,SS_D2,O_MUL,O_SUB,
	SS_N3,SS_N2,SS_N1,F_ETS,SS_N3,SS_N2,SS_N1,F_ETC,O_MUL,SS_D1,O_MUL,O_ADD,
	CON_1,SS_N3,SS_N2,SS_N1,F_EOQ,O_SUB,O_DIV,0 
 };

static int funct_f_ets(double *s)
 {	s-=3;
	*s=sin((*s)+eccentric_offset_p(*s,*(s+1),*(s+2)));
	return(0);
 }

static short    diffrule_f_ets[]=
 {	SS_N3,SS_N2,SS_N1,F_ETC,SS_D3,O_MUL,
	SS_N3,SS_N2,SS_N1,F_ETS,SS_N3,SS_N2,SS_N1,F_ETC,O_MUL,SS_D2,O_MUL,O_ADD,
	SS_N3,SS_N2,SS_N1,F_ETC,O_PSQ,SS_D1,O_MUL,O_SUB,
	CON_1,SS_N3,SS_N2,SS_N1,F_EOQ,O_SUB,O_DIV,0 
 };


static int funct_f_ntiu(double *s)
 {	s-=2;
	*s=1.0-ntiq(*s,*(s+1),0.0,0.0);
  	return(0);
 }


static short	diffrule_f_ntiq[]=
 {	SS_N4,SS_N3,SS_N2,SS_N1,F_NTIQ_DP ,SS_D4,O_MUL,
	SS_N4,SS_N3,SS_N2,SS_N1,F_NTIQ_DZ ,SS_D3,O_MUL,O_ADD,
	SS_N4,SS_N3,SS_N2,SS_N1,F_NTIQ_DG1,SS_D2,O_MUL,O_ADD,
	SS_N4,SS_N3,SS_N2,SS_N1,F_NTIQ_DG2,SS_D1,O_MUL,O_ADD,0
 };

static int funct_f_ntiq_diff(double *s,double *out)
{	static	double	st_args[4]={1,5,0,0},st_diff[5]={0,0,0,0,0};
	if ( memcmp(s,st_args,sizeof(double)*4) != 0 )
	 {	memcpy(st_args,s,sizeof(double)*4);
		st_diff[4]=ntiq_diff(s[0],s[1],s[2],s[3],st_diff);
	 }
	memcpy(out,st_diff,sizeof(double)*5);
	return(0);
}

static int funct_f_ntiq(double *s)
 {	double	st_diff[5];
	s-=4;
	funct_f_ntiq_diff(s,st_diff);
	*s=(1.0-st_diff[4]);
  	return(0);
 }

static int funct_f_ntiq_dp(double *s)
 {	double	st_diff[5];
	s-=4;
	funct_f_ntiq_diff(s,st_diff);
	*s=(-st_diff[0]);
	return(0);
 }
static int funct_f_ntiq_dz(double *s)
 {	double	st_diff[5];
	s-=4;
	funct_f_ntiq_diff(s,st_diff);
	*s=(-st_diff[1]);
	return(0);
 }
static int funct_f_ntiq_dg1(double *s)
 {	double	st_diff[5];
	s-=4;
	funct_f_ntiq_diff(s,st_diff);
	*s=(-st_diff[2]);
	return(0);
 }
static int funct_f_ntiq_dg2(double *s)
 {	double	st_diff[5];
	s-=4;
	funct_f_ntiq_diff(s,st_diff);
	*s=(-st_diff[3]);
	return(0);
 }

static int funct_f_hjd(double *s)
 {	s-=3;
	*s=get_hjd(*s,*(s+1),*(s+2));
	return(0);
 }

static int funct_f_bjd(double *s)
 {	s-=3;
	*s=get_bjd(*s,*(s+1),*(s+2));
	return(0);
 }

static int funct_f_elliptic_k(double *s)
{
	s-=1;
	if ( fabs(*s) < 1.0 )
	 {	*s=elliptic_complete_first(*s);
		return(0);
	 }
	else
		return(1);
}

static int funct_f_elliptic_e(double *s)
{
	s-=1;
	if ( fabs(*s) <= 1.0 )
	 {	*s=elliptic_complete_second(*s);
		return(0);
	 }
	else
		return(1);
}
static int funct_f_elliptic_pi(double *s)
{
	s-=2;
	if ( *s < 1.0 && fabs(*(s+1)) < 1.0 )
	 {	*s=elliptic_complete_third(*s,*(s+1));
		return(0);
	 }
	else
		return(1);
}

static psnfunctinfo	pfi_eoq		= { "eccentric offset function, `q` component (arguments: mean longitude in radians, k and h Lagrangian elements)" };
static psnfunctinfo	pfi_eop		= { "eccentric offset function, `p` component (arguments: mean longitude in radians, k and h Lagrangian elements)" };
static psnfunctinfo	pfi_etc		= { "eccentric trigonometric function, cosine component (arguments: mean longitude in radians, k and h Lagrangian elements)" };
static psnfunctinfo	pfi_ets		= { "eccentric trigonometric function, sine component (arguments: mean longitude in radians, k and h Lagrangian elements)" };
static psnfunctinfo	pfi_ntiu	= { "normalized transit intensity - uniform flux distribution (arguments: fractional radius, normalized distance)" };
static psnfunctinfo	pfi_ntiq	= { "normalized transit intensity - quadratic limb darkening assumption (arguments: fracional radius, normalized distance, the two limb darkening coefficients)" };
static psnfunctinfo	pfi_hjd		= { "heliocentric Julian date (arguments: julian day, RA and DEC, in degrees)" };
static psnfunctinfo	pfi_bjd		= { "barycentric Julian date (arguments: julian day, RA and DEC, in degrees)" };
static psnfunctinfo	pfi_ell_k	= { "complete elliptic integral of the first kind" };
static psnfunctinfo	pfi_ell_e	= { "complete elliptic integral of the second kind" };
static psnfunctinfo	pfi_ell_pi	= { "complete elliptic integral of the third kind" };

psnlfit psnlfit_list_builtin_xfuncts[]=
{                                                                               
  { "eoq"       , T_FN, F_EOQ     , 0        , 3,  0, 0          , funct_f_eoq        , diffrule_f_eoq  , "eoq(#1,#2,#3)"        , 0, 0   , &pfi_eoq  	},
  { "eop"       , T_FN, F_EOP     , 0        , 3,  0, 0          , funct_f_eop        , diffrule_f_eop  , "eop(#1,#2,#3)"        , 0, 0   , &pfi_eop  	},
  { "etc"       , T_FN, F_ETC     , 0        , 3,  0, 0          , funct_f_etc        , diffrule_f_etc  , "etc(#1,#2,#3)"        , 0, 0   , &pfi_etc  	},
  { "ets"       , T_FN, F_ETS     , 0        , 3,  0, 0          , funct_f_ets        , diffrule_f_ets  , "ets(#1,#2,#3)"        , 0, 0   , &pfi_ets  	},
  { "ntiu"      , T_FN, F_NTIU    , 0        , 2,  0, 0          , funct_f_ntiu       , NULL            , "ntiu(#1,#2)"          , 0, 0   , &pfi_ntiu  	},
  { "ntiq"      , T_FN, F_NTIQ    , 0        , 4,  0, 0          , funct_f_ntiq       , diffrule_f_ntiq , "ntiq(#1,#2,#3,#4)"    , 0, 0   , &pfi_ntiq  	},
  { "ntiq_dp"   , T_FN, F_NTIQ_DP , 0        , 4,  0, 0          , funct_f_ntiq_dp    , NULL            , "ntiq_dp(#1,#2,#3,#4)" , 0, 0   , NULL  	},
  { "ntiq_dz"   , T_FN, F_NTIQ_DZ , 0        , 4,  0, 0          , funct_f_ntiq_dz    , NULL            , "ntiq_dz(#1,#2,#3,#4)" , 0, 0   , NULL  	},
  { "ntiq_dg1"  , T_FN, F_NTIQ_DG1, 0        , 4,  0, 0          , funct_f_ntiq_dg1   , NULL            , "ntiq_dg1(#1,#2,#3,#4)", 0, 0   , NULL  	},
  { "ntiq_dg2"  , T_FN, F_NTIQ_DG2, 0        , 4,  0, 0          , funct_f_ntiq_dg2   , NULL            , "ntiq_dg2(#1,#2,#3,#4)", 0, 0   , NULL  	},
  { "hjd"       , T_FN, F_HJD     , 0        , 3,  0, 0          , funct_f_hjd        , NULL            , "hjd(#1,#2,#3)"        , 0, 0   , &pfi_hjd	},
  { "bjd"       , T_FN, F_BJD     , 0        , 3,  0, 0          , funct_f_bjd        , NULL            , "bjd(#1,#2,#3)"        , 0, 0   , &pfi_bjd	},
  { "ellipticK" , T_FN, F_ELL_K   , 0        , 1,  0, 0          , funct_f_elliptic_k , NULL            , "ellipticK(#1)"        , 0, 0   , &pfi_ell_k	},
  { "ellipticE" , T_FN, F_ELL_E   , 0        , 1,  0, 0          , funct_f_elliptic_e , NULL            , "ellipticE(#1)"        , 0, 0   , &pfi_ell_e	},
  { "ellipticPi", T_FN, F_ELL_PI  , 0        , 2,  0, 0          , funct_f_elliptic_pi, NULL            , "ellipticPi(#1)"       , 0, 0   , &pfi_ell_pi	},
  { NULL        , 0   , 0         , 0        , 0,  0, 0          , NULL               , NULL            , NULL                   , 0, 0   , NULL  	}
};

#endif

/*****************************************************************************/
                                                                
