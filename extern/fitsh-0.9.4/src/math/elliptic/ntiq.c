/*****************************************************************************/
/* ntiq.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Functions related to the calculations of obscured stellar flux taking     */
/* into account quadratic limb darkening. These functions also provide the   */
/* parametric partial derivatives of the obscuring function. A bit redundant */
/* coding, however, straight and hopefully clean.			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* If you use these code in your research, please cite astro-ph:0805.2157,   */
/* and if you are using the flux decrease itself (i.e. not only the partial  */
/* derivatives), please cite 2002, ApJ, 580, 171 too.			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2007; Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "elliptic.h"
#include "ntiq.h"

/*****************************************************************************/

#define		SMALL		(1e-13)

/*****************************************************************************/

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841968
#endif

/*****************************************************************************/

double ntiq(double p,double z,double g1,double g2)
{
 double	w0,w1,w2;
 double	f0,f1,fk,fe,fp,f2,n,k;
 double	df;

 w0=(6-6*g1-12*g2)/(6-2*g1-g2);
 w1=( 6*g1+12*g2 )/(6-2*g1-g2);
 w2=(	6*g2	 )/(6-2*g1-g2);

 z=fabs(z);

 if ( z<=0.0 && p<=1.0 )				/* A	*/ /* M&A 10 */
  {	f0=p*p;
	f1=2.0/3.0*(1.0-(1.0-f0)*sqrt(1.0-f0));
	fk=fe=fp=0.0;
	f2=0.5*f0*f0;
	k=n=0.0;
  }
 else if ( z<=p-1.0 )					/* A_G	*/ /* M&A 11 */
  {	f0=1.0;
	fk=fe=fp=0.0;
	f1=2.0/3.0;
	f2=0.5;
	k=n=0.0;
  }	
 else if ( z<p && z<1.0-p-SMALL )			/* B	*/ /* M&A  9 */
  {	double	a,b,ci,cik,cie,cip;
	a=(p-z)*(p-z);
	b=(p+z)*(p+z);
	
	f0=p*p;
	f1=2.0/3.0;
	ci=2.0/(9.0*M_PI*sqrt(1.0-a));
	cik=(1.0-5*z*z+f0+a*b);
	cie=(z*z+7*f0-4)*(1.0-a);
	cip=-3.0*(p+z)/(p-z);
	fk=ci*cik;
	fe=ci*cie;
	fp=ci*cip;
	f2=0.5*f0*(f0+2*z*z);
	k=sqrt(4.0*z*p/(1.0-a));
	n=(a-b)/a;
  }
 else if ( z<p && fabs(z-1.0+p)<=SMALL )		/* B_T	*/ /* M&A  - */
  {	f0=p*p;
	f1=2.0/(3.0*M_PI)*acos(1-2*p)-4.0/(9.0*M_PI)*(3+2*p-8*f0)*sqrt(p*(1-p));
	fk=0.0;
	fe=0.0;
	fp=0.0;
	f2=0.5*f0*(f0+2*z*z);
	k=n=0.0;
  }
 else if ( z<p )					/* B_G	*/ /* M&A  8 */
  {	double	a,b,k0,k1,cg,cgk,cge,cgp;
	a=(p-z)*(p-z);
	b=(p+z)*(p+z);

	k0=acos((p*p+z*z-1)/(2*p*z));
	k1=acos((1+z*z-p*p)/(2*z));
	f0=(p*p*k0+k1-sqrt(z*z-0.25*(1+z*z-p*p)*(1+z*z-p*p)))/M_PI;
	f1=2.0/3.0;		/* difference between 'B_G' and 'F' cases! */
	cg=1.0/(9.0*M_PI*sqrt(p*z));
	cgk=((1-b)*(2*b+a-3)-3*(p+z)*(p-z)*(b-2));
	cge=4*p*z*(z*z+7*p*p-4);
	cgp=-3*(p+z)/(p-z);
	fk=cg*cgk;
	fe=cg*cge;
	fp=cg*cgp;
	f2=(k1+p*p*(p*p+2*z*z)*k0-0.25*(1+5*p*p+z*z)*sqrt((1-a)*(b-1)))/(2.0*M_PI);
	k=sqrt((1-a)/(4*p*z));
	n=(a-1)/a;
  }
 else if ( fabs(z-p)<=SMALL && z<1.0-p-SMALL )		/* C	*/ /* M&A  5 */
  {	double	t;

	f0=p*p;
	f1=1.0/3.0;
	
	t=2.0/(9.0*M_PI);
	fk=t*(1.0-4.0*f0);
	fe=t*4*(2.0*f0-1.0);
	fp=0.0;
	f2=1.5*f0*f0;
	k=2*p;
	n=0.0;
  }
 else if ( fabs(z-p)<=SMALL && fabs(z-1.0+p)<=SMALL )	/* C_T	*/ /* M&A  6 */
  {	f0=0.25;
	f1=1.0/3.0-4.0/(9.0*M_PI);
	f2=3.0/32.0;
	fk=fe=fp=0.0;
	n=k=0.0;
  }
 else if ( fabs(z-p)<=SMALL )				/* C_G	*/ /* M&A  7 */
  {	double	a,b,k0,k1;
	a=(p-z)*(p-z);
	b=(p+z)*(p+z);

	k0=acos((p*p+z*z-1)/(2*p*z));
	k1=acos((1+z*z-p*p)/(2*z));
	f0=(p*p*k0+k1-sqrt(z*z-0.25*(1+z*z-p*p)*(1+z*z-p*p)))/M_PI;
	f1=1.0/3.0;
	fk=-(1-4*p*p)*(3-8*p*p)/(9*M_PI*p);
	fe=16*p*(2*p*p-1)/(9*M_PI);
	fp=0.0;
	f2=(k1+p*p*(p*p+2*z*z)*k0-0.25*(1+5*p*p+z*z)*sqrt((1-a)*(b-1)))/(2.0*M_PI);
	k=1.0/(2.0*p);
	n=0.0;	
  }
 else if ( z<1-p-SMALL )				/* D	*/ /* M&A  3 */
  {	double	a,b,ci,cik,cie,cip;
	a=(p-z)*(p-z);
	b=(p+z)*(p+z);
	
	f0=p*p;
	f1=0.0;
	ci=2.0/(9.0*M_PI*sqrt(1.0-a));
	cik=(1.0-5*z*z+f0+a*b);
	cie=(z*z+7*f0-4)*(1.0-a);
	cip=-3.0*(p+z)/(p-z);
	fk=ci*cik;
	fe=ci*cie;
	fp=ci*cip;
	f2=0.5*f0*(f0+2*z*z);
	k=sqrt(4.0*z*p/(1.0-a));
	n=(a-b)/a;
  }
 else if ( fabs(z-1.0+p)<=SMALL )			/* E	*/ /* M&A  4 */
  {	f0=p*p;
	f1=2.0/(3.0*M_PI)*acos(1-2*p)-4.0/(9.0*M_PI)*(3+2*p-8*f0)*sqrt(p*(1-p));
	fk=0.0;
	fe=0.0;
	fp=0.0;
	f2=0.5*f0*(f0+2*z*z);
	k=n=0.0;
  }
 else if ( z<1+p-SMALL )				/* F	*/ /* M&A  2 */
  {	double	a,b,k0,k1,cg,cgk,cge,cgp;
	a=(p-z)*(p-z);
	b=(p+z)*(p+z);

	k0=acos((p*p+z*z-1)/(2*p*z));
	k1=acos((1+z*z-p*p)/(2*z));
	f0=(p*p*k0+k1-sqrt(z*z-0.25*(1+z*z-p*p)*(1+z*z-p*p)))/M_PI;
	f1=0.0;
	cg=1.0/(9.0*M_PI*sqrt(p*z));
	cgk=((1-b)*(2*b+a-3)-3*(p+z)*(p-z)*(b-2));
	cge=4*p*z*(z*z+7*p*p-4);
	cgp=-3*(p+z)/(p-z);
	fk=cg*cgk;
	fe=cg*cge;
	fp=cg*cgp;
	f2=(k1+p*p*(p*p+2*z*z)*k0-0.25*(1+5*p*p+z*z)*sqrt((1-a)*(b-1)))/(2.0*M_PI);
	k=sqrt((1-a)/(4*p*z));
	n=(a-1)/a;
  }
 else							/* G	*/ /* M&A  1 */
  {	f0=0.0;
	f1=0.0;
	fk=0.0;
	fe=0.0;
	fp=0.0;
	f2=0.0;
	k=0.0;
	n=0.0;
  }

 df=w0*f0+w1*f1+w2*f2;

 if ( fk != 0.0 || fe != 0.0 )
  {	double	q,ek,ee,ep;

	q=(1-k)*(1+k);

	ek=carlson_elliptic_rf(0.0,q,1.0);
 	df+=w1*fk*ek;
	ee=ek-k*k*carlson_elliptic_rd(0.0,q,1.0)/3.0;
	df+=w1*fe*ee;
	if ( fp != 0.0 )
	 {	ep=ek+n*carlson_elliptic_rj(0.0,q,1.0,1.0-n)/3.0;
		df+=w1*fp*ep;
	 }
  }	
 
 return(df);
}

double ntiq_diff(double p,double z,double g1,double g2,double *diff)
{
 double	signz,w,w0,w1,w2,w0g1,w0g2,w1g1,w1g2,w2g1,w2g2;
 double	f0,f1,fk,fe,fp,f2,n,k;
 double	f0pp,f0zz,f1pp,f1zz,fkpp,fkzz,fepp,fezz,fppp,fpzz,f2pp,f2zz,
	npp,nzz,kpp,kzz;

 double	df;
 double	ek,ee,ep;
 int	calc_eke,calc_epi;

 if ( diff==NULL )
  {	df=ntiq(p,z,g1,g2);
	return(df);
  }

 w=(6-2*g1-g2);
 w0=(6-6*g1-12*g2)/w;
 w0g1=(2*w0-6)/w;
 w0g2=(w0-12)/w;
 w1=( 6*g1+12*g2 )/w;
 w1g1=(6+2*w1)/w;
 w1g2=(12+w1)/w; 
 w2=(	6*g2	 )/w;
 w2g1=2*w2/w;
 w2g2=(6+w2)/w;

 if ( z<0.0 )	signz=-1.0,z=-z;
 else		signz=+1.0;

 if ( z<=0.0 && p<=1.0 )				/* A	*/
  {	f0=p*p;
	f0pp=2*p;
	f0zz=0.0;

	f1=2.0/3.0*(1.0-(1.0-f0)*sqrt(1.0-f0));
	f1pp=2*p*sqrt(1.0-f0);
	f1zz=0.0;

	fk=fe=fp=0.0;
	fkpp=fepp=fppp=0.0;
	fkzz=fezz=fpzz=0.0;

	f2=0.5*f0*f0;
	f2pp=2*f0*p;
	f2zz=0.0;

	k=0.0;
	n=0.0;
	kpp=kzz=0.0;
	npp=nzz=0.0;
	calc_eke=0;
	calc_epi=0;
  }
 else if ( z<=p-1.0 )					/* A_G	*/
  {	f0=1.0;
	f0pp=0.0;
	f0zz=0.0;
	fk=fe=fp=0.0;
	fkpp=fepp=fppp=0.0;
	fkzz=fezz=fpzz=0.0;
	f1=2.0/3.0;
	f1pp=0.0;
	f1zz=0.0;
	f2=0.5;
	f2pp=0.0;
	f2zz=0.0;
	k=n=0.0;
	kpp=kzz=0.0;
	npp=nzz=0.0;
	calc_eke=0;
	calc_epi=0;
  }	
 else if ( z<p && z<1.0-p-SMALL )			/* B	*/
  {	double	a,b,ci,cik,cie,cip;
	a=(p-z)*(p-z);
	b=(p+z)*(p+z);
	
	f0=p*p;
	f0pp=2*p;
	f0zz=0.0;

	f1=2.0/3.0;
	f1pp=0.0;
	f1zz=0.0;

	ci=2.0/(9.0*M_PI*sqrt(1.0-a));
	cik=(1.0-5*z*z+f0+a*b);
	cie=(z*z+7*f0-4)*(1.0-a);
	cip=-3.0*(p+z)/(p-z);
	fk=ci*cik;
	fkpp=+ci*2*p*(1+2*(p*p-z*z))+ci*cik*(p-z)/(1-a);
	fkzz=-ci*2*z*(5+2*(p*p-z*z))-ci*cik*(p-z)/(1-a);
	fe=ci*cie;
	fepp=+ci*(14*p*(1-a)-2*(p-z)*(z*z+7*f0-4))+ci/(1-a)*(p-z)*cie;
	fezz=+ci*( 2*z*(1-a)+2*(p-z)*(z*z+7*f0-4))-ci/(1-a)*(p-z)*cie;
	fp=ci*cip;
	fppp=+ci*3.0*(+(2*z)/a-(p+z)/(1-a));
	fpzz=+ci*3.0*(-(2*p)/a+(p+z)/(1-a));

	f2=0.5*f0*(f0+2*z*z);
	f2pp=p*(f0+2*z*z)+f0*p;
	f2zz=2*f0*z;
	
	k=sqrt(4.0*z*p/(1.0-a));
	kpp=2*z*(1+p*p-z*z)/((1-a)*(1-a)*k);
	kzz=2*p*(1-p*p+z*z)/((1-a)*(1-a)*k);
	n=-4.0*z*p/a;
	npp=+4*z*(p+z)/((p-z)*a);
	nzz=-4*p*(p+z)/((p-z)*a);

	calc_eke=1;
	calc_epi=1;
  }
 else if ( z<p && fabs(z-1.0+p)<=SMALL )		/* B_T	*/
  {	f0=p*p;
	f0pp=2*p;
	f0zz=0.0;

	f1=2.0/(3.0*M_PI)*acos(1-2*p)-4.0/(9.0*M_PI)*(3+2*p-8*f0)*sqrt(p*(1-p));
	f1pp=(8.0/M_PI)*p*sqrt(p*(1-p));
	f1zz=-f1pp/3.0;

	fk=fe=fp=0.0;
	fkpp=fepp=fppp=0.0;
	fkzz=fezz=fpzz=0.0;

	f2=0.5*f0*(f0+2*z*z);
	f2pp=p*(f0+2*z*z)+f0*p;
	f2zz=2*f0*z;

	k=0.0;
	kpp=0.0;
	kzz=0.0;
	n=0.0;
	npp=0.0;
	nzz=0.0;

	calc_eke=0;
	calc_epi=0;
  }
 else if ( z<p )					/* B_G	*/
  {	double	a,b,k0,k1,cg,cgk,cge,cgp,t;
	a=(p-z)*(p-z);
	b=(p+z)*(p+z);

	k0=acos((p*p+z*z-1)/(2*p*z));
	t=(1+z*z-p*p)/(2*z);
	k1=acos(t);

	t=sqrt(z*z-0.25*(1+z*z-p*p)*(1+z*z-p*p));
	f0=(p*p*k0+k1-t)/M_PI;
	f0pp=(2*p*k0)/M_PI;
	f0zz=(-2*p*sin(k0))/M_PI;
	
	f1=2.0/3.0;
	f1pp=0.0;
	f1zz=0.0;

	cg=1.0/(9.0*M_PI*sqrt(p*z));
	cgk=((1-b)*(2*b+a-3)-3*(p+z)*(p-z)*(b-2));
	cge=4*p*z*(z*z+7*p*p-4);
	cgp=-3*(p+z)/(p-z);
	fk=cg*cgk;
	fkpp=-fk/(2*p)-cg*2*(p*p*(12*p+21*z)+z*(z*z-4)+2*p*(5*z*z-6));
	fkzz=-fk/(2*z)-cg*2*p*(-4+7*p*p+10*p*z+3*z*z);
	fe=cg*cge;
	fepp=-fe/(2*p)+cg*4*z*(-4+21*p*p+  z*z);
	fezz=-fe/(2*z)+cg*4*p*(-4+ 7*p*p+3*z*z);
	fp=cg*cgp;
	fppp=-fp/(2*p)+cg*6*z/a;
	fpzz=-fp/(2*z)-cg*6*p/a;

	t=sqrt((1-a)*(b-1));
	f2=(k1+p*p*(p*p+2*z*z)*k0-0.25*(1+5*p*p+z*z)*t)/(2.0*M_PI);
	f2pp=(2*p*(p*p+z*z)*k0-4*z*p*p*sin(k0))/M_PI;
	f2zz=(-2*p*(p*p+z*z)*sin(k0)+z*p*p*(2*k0+sin(2*k0)))/M_PI;

	k=sqrt((1-a)/(4*p*z));
	kpp=-(1+p*p-z*z)/(8*k*p*p*z);
	kzz=-(1-p*p+z*z)/(8*k*p*z*z);
	n=(a-1)/a;
	npp=2/(a*(p-z));
	nzz=-npp;

	calc_eke=1;
	calc_epi=1;
  }

 else if ( fabs(z-p)<=SMALL && z<1.0-p-SMALL )		/* C	*/
  {	double	t;

	f0=p*p;
	f0pp=2*p;
	f0zz=0.0;

	f1=1.0/3.0;
	f1pp=0.0;
	f1zz=0.0;
	
	t=2.0/(9.0*M_PI);

	fk=t*(1.0-4.0*f0);
	fkpp=+2.0/(9.0*M_PI)* 2*p-1.0/(3*M_PI*p);
	fkzz=-2.0/(9.0*M_PI)*10*p+1.0/(3*M_PI*p);

	fe=t*4*(2.0*f0-1.0);
	fepp=+t*(14*p)+1.0/(3*M_PI*p);
	fezz=+t*( 2*p)-1.0/(3*M_PI*p);

	fp=0.0;
	fppp=0.0;
	fpzz=0.0;

	f2=1.5*f0*f0;
	f2pp=4*f0*p;
	f2zz=2*f0*p;

	k=2*p;
	kpp=1.0;
	kzz=1.0;
	n=0.0;
	npp=0.0;
	nzz=0.0;

	calc_eke=1;
	calc_epi=0;
  }

 else if ( fabs(z-p)<=SMALL && fabs(z-1.0+p)<=SMALL )	/* C_T	*/
  {	f0=0.25;
	f0pp=1.0;
	f0zz=0.0;

	f1=1.0/3.0-4.0/(9*M_PI);
	f1pp=+2.0/M_PI;
	f1zz=-2.0/(3*M_PI);
	
	fk=0.0;
	fkpp=0.0;
	fkzz=0.0;

	fe=0.0;
	fepp=0.0;
	fezz=0.0;

	fp=0.0;
	fppp=0.0;
	fpzz=0.0;

	f2=3.0/32.0;
	f2pp=0.5;
	f2zz=0.25;

	k=0.0;
	kpp=1.0;
	kzz=1.0;
	n=0.0;
	npp=0.0;
	nzz=0.0;

	calc_eke=0;
	calc_epi=0;
  }


 else if ( fabs(z-p)<=SMALL )				/* C_G	*/
  {	double	k0,k1,cg,cgk,cge,t;

	k0=acos((p*p+z*z-1)/(2*p*z));
	t=(1+z*z-p*p)/(2*z);
	k1=acos(t);
	t=sqrt(z*z-0.25*(1+z*z-p*p)*(1+z*z-p*p));
	f0=(p*p*k0+k1-t)/M_PI;
	f0pp=(2*p*k0)/M_PI;
	f0zz=(-2*p*sin(k0))/M_PI;

	fk=-(1-4*p*p)*(3-8*p*p)/(9*M_PI*p);
	fe=16*p*(2*p*p-1)/(9*M_PI);
	fp=0.0;

	f1=1.0/3.0;
	f1pp=0.0;
	f1zz=0.0;

	cg=1.0/(9.0*M_PI*p);
	cgk=(1-4*p*p)*(8*p*p-3);
	cge=16*p*p*(2*p*p-1);
	fk=cg*cgk;
	fkpp=(3+16*p*p*(2-9*p*p))/(18*M_PI*p*p);
	fkzz=(3+ 8*p*p*(1-6*p*p))/(18*M_PI*p*p);
	fe=cg*cge;
	fepp=( -2+72*p*p)/(9*M_PI);
	fezz=(-14+24*p*p)/(9*M_PI);
	fp=0.0;
	fppp=0.0;
	fpzz=0.0;

	t=sqrt(4*p*p-1);
	f2=(k1+p*p*(3*p*p)*k0-0.25*(1+6*p*p)*t)/(2.0*M_PI);
	f2pp=(2*p*(p*p+z*z)*k0-4*z*p*p*sin(k0))/M_PI;
	f2zz=(-2*p*(p*p+z*z)*sin(k0)+z*p*p*(2*k0+sin(2*k0)))/M_PI;

	k=1.0/(2*p);
	kpp=-(1.0)/(4*p*p);
	kzz=-(1.0)/(4*p*p);
	n=0.0;
	npp=0.0;
	nzz=0.0;

	calc_eke=1;
	calc_epi=0;
  }
 else if ( z<1-p-SMALL )				/* D	*/
  {	double	a,b,ci,cik,cie,cip;
	a=(p-z)*(p-z);
	b=(p+z)*(p+z);
	
	f0=p*p;
	f0pp=2*p;
	f0zz=0.0;

	f1=0.0;
	f1pp=0.0;
	f1zz=0.0;

	ci=2.0/(9.0*M_PI*sqrt(1.0-a));
	cik=(1.0-5*z*z+f0+a*b);
	cie=(z*z+7*f0-4)*(1.0-a);
	cip=-3.0*(p+z)/(p-z);
	fk=ci*cik;
	fkpp=+ci*2*p*(1+2*(p*p-z*z))+ci*cik*(p-z)/(1-a);
	fkzz=-ci*2*z*(5+2*(p*p-z*z))-ci*cik*(p-z)/(1-a);
	fe=ci*cie;
	fepp=+ci*(14*p*(1-a)-2*(p-z)*(z*z+7*f0-4))+ci/(1-a)*(p-z)*cie;
	fezz=+ci*( 2*z*(1-a)+2*(p-z)*(z*z+7*f0-4))-ci/(1-a)*(p-z)*cie;
	fp=ci*cip;
	fppp=+ci*3.0*(+(2*z)/a-(p+z)/(1-a));
	fpzz=+ci*3.0*(-(2*p)/a+(p+z)/(1-a));

	f2=0.5*f0*(f0+2*z*z);
	f2pp=p*(f0+2*z*z)+f0*p;
	f2zz=2*f0*z;
	
	k=sqrt(4.0*z*p/(1.0-a));
	kpp=2*z*(1+p*p-z*z)/((1-a)*(1-a)*k);
	kzz=2*p*(1-p*p+z*z)/((1-a)*(1-a)*k);
	n=-4.0*z*p/a;
	npp=+4*z*(p+z)/((p-z)*a);
	nzz=-4*p*(p+z)/((p-z)*a);

	calc_eke=1;
	calc_epi=1;
  }

 else if ( fabs(z-1.0+p)<=SMALL )			/* E	*/
  {	f0=p*p;
	f0pp=2*p;
	f0zz=0.0;

	f1=2.0/(3.0*M_PI)*acos(1-2*p)-4.0/(9.0*M_PI)*(3+2*p-8*f0)*sqrt(p*(1-p));
	f1pp=(8.0/M_PI)*p*sqrt(p*(1-p));
	f1zz=-f1pp/3.0;

	fk=fe=fp=0.0;
	fkpp=fepp=fppp=0.0;
	fkzz=fezz=fpzz=0.0;

	f2=0.5*f0*(f0+2*z*z);
	f2pp=p*(f0+2*z*z)+f0*p;
	f2zz=2*f0*z;

	k=0.0;
	kpp=0.0;
	kzz=0.0;
	n=0.0;
	npp=0.0;
	nzz=0.0;

	calc_eke=0;
	calc_epi=0;
  }
 else if ( z<1+p-SMALL )				/* F	*/
  {	double	a,b,k0,k1,cg,cgk,cge,cgp,t;
	a=(p-z)*(p-z);
	b=(p+z)*(p+z);

	k0=acos((p*p+z*z-1)/(2*p*z));
	t=(1+z*z-p*p)/(2*z);
	k1=acos(t);

	t=sqrt(z*z-0.25*(1+z*z-p*p)*(1+z*z-p*p));
	f0=(p*p*k0+k1-t)/M_PI;
	f0pp=(2*p*k0)/M_PI;
	f0zz=(-2*p*sin(k0))/M_PI;
	
	f1=0.0;
	f1pp=0.0;
	f1zz=0.0;

	cg=1.0/(9.0*M_PI*sqrt(p*z));
	cgk=((1-b)*(2*b+a-3)-3*(p+z)*(p-z)*(b-2));
	cge=4*p*z*(z*z+7*p*p-4);
	cgp=-3*(p+z)/(p-z);
	fk=cg*cgk;
	fkpp=-fk/(2*p)-cg*2*(p*p*(12*p+21*z)+z*(z*z-4)+2*p*(5*z*z-6));
	fkzz=-fk/(2*z)-cg*2*p*(-4+7*p*p+10*p*z+3*z*z);
	fe=cg*cge;
	fepp=-fe/(2*p)+cg*4*z*(-4+21*p*p+  z*z);
	fezz=-fe/(2*z)+cg*4*p*(-4+ 7*p*p+3*z*z);
	fp=cg*cgp;
	fppp=-fp/(2*p)+cg*6*z/a;
	fpzz=-fp/(2*z)-cg*6*p/a;

	t=sqrt((1-a)*(b-1));
	f2=(k1+p*p*(p*p+2*z*z)*k0-0.25*(1+5*p*p+z*z)*t)/(2.0*M_PI);
	f2pp=(2*p*(p*p+z*z)*k0-4*z*p*p*sin(k0))/M_PI;
	f2zz=(-2*p*(p*p+z*z)*sin(k0)+z*p*p*(2*k0+sin(2*k0)))/M_PI;

	k=sqrt((1-a)/(4*p*z));
	kpp=-(1+p*p-z*z)/(8*k*p*p*z);
	kzz=-(1-p*p+z*z)/(8*k*p*z*z);
	n=(a-1)/a;
	npp=2/(a*(p-z));
	nzz=-npp;

	calc_eke=1;
	calc_epi=1;
  }
 else				/* G	*/
  {	f0=0.0;
	f0pp=0.0;
	f0zz=0.0;

	f1=0.0;
	f1pp=0.0;
	f1zz=0.0;

	fk=fe=fp=0.0;
	fkpp=fepp=fppp=0.0;
	fkzz=fezz=fpzz=0.0;

	f2=0.0;
	f2pp=0.0;
	f2zz=0.0;

	k=n=0.0;
	kpp=kzz=0.0;
	npp=nzz=0.0;

	calc_eke=0;
	calc_epi=0;
  }
 
 df=w0*f0+w1*f1+w2*f2;
 diff[0]=w0*f0pp+w1*f1pp+w2*f2pp;
 diff[1]=w0*f0zz+w1*f1zz+w2*f2zz;
 diff[2]=w0g1*f0+w1g1*f1+w2g1*f2;
 diff[3]=w0g2*f0+w1g2*f1+w2g2*f2;

 if ( calc_eke )
  {	double	q;

	q=(1-k)*(1+k);

	/*
	ek=elliptic_complete_first (k);
	*/
	ek=carlson_elliptic_rf(0.0,q,1.0);
	df+=w1*fk*ek;
	diff[0]+=w1*ek*(fkpp-(fk+fe)*kpp/k+(fp!=0.0?fp/(2*n*(n-1))*npp:0));
	diff[1]+=w1*ek*(fkzz-(fk+fe)*kzz/k+(fp!=0.0?fp/(2*n*(n-1))*nzz:0));
	diff[2]+=w1g1*fk*ek;
	diff[3]+=w1g2*fk*ek;
	/*
 	ee=elliptic_complete_second(k);
	*/
	ee=ek-k*k*carlson_elliptic_rd(0.0,q,1.0)/3.0;
	df+=w1*fe*ee;
	diff[0]+=w1*ee*(fepp+(fk/(k*(1-k*k))+fe/k+(fp!=0.0?fp*k/((n-k*k)*(k*k-1)):0))*kpp+(fp!=0.0?fp/(2*(k*k-n)*(n-1))*npp:0));
	diff[1]+=w1*ee*(fezz+(fk/(k*(1-k*k))+fe/k+(fp!=0.0?fp*k/((n-k*k)*(k*k-1)):0))*kzz+(fp!=0.0?fp/(2*(k*k-n)*(n-1))*nzz:0));
	diff[2]+=w1g1*fe*ee;
	diff[3]+=w1g2*fe*ee;

	if ( calc_epi )
	 {	/*
		ep=elliptic_complete_third (n,k);
		*/
		ep=ek+n*carlson_elliptic_rj(0.0,q,1.0,1.0-n)/3.0;
		df+=w1*fp*ep;
		diff[2]+=w1g1*fp*ep;
		diff[3]+=w1g2*fp*ep;
	 }
  }

 diff[1]=diff[1]*signz;

 return(df);
}

/*****************************************************************************/
    
