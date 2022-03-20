/*****************************************************************************/
/* dtc.c 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Calculations between different types of date formats			     */
/* This part of the library (dtc.[ch]) is standalone!			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2003, 2006; Pal, A. (apal@szofi.elte.hu).			     */
/*****************************************************************************/

#include <math.h>
#include <time.h>

#include <astro/dtc.h>

double read_julian_date(void)
{
 time_t	ctime;
 double	jd;

/* tzset(); */
 time(&ctime);

 jd=((double)ctime)/86400.0 + 2440587.5 ;

 return(jd);
}

double get_julian_date(int ye,int mo,int da)
{
 int	lya;
 double jd;

 if ( mo<=2 )	ye--,mo+=12;

 if ( ye>1582 || ( ye==1582 && mo>10 ) || ( ye==1582 && mo==10 && da>14 ) )
	lya=(2)-(ye/100)+(ye/400);
 else
	lya=0;

 jd=	(double)lya+
	1720994.5+
	(double)(((long)ye*1461L)/4L)+
	(double)((153*(mo+1))/5)+
	(double)da;

 return(jd);
}

void get_calendar_date(double jd,int *ye,int *mo,double *da)
{
 long	w,a,b,c,d,e,z;
 
 jd+=0.5;
 z=(long)jd;

 if ( z>=2299161L )
	w=(4L*z-7468865L)/146097L,a=1L+z+w-(w/4L);
 else	
	a=z;

 b=a+1524L;
 c=(20L*b-2442L)/7305L;
 d=(1461L*c)/4L;
 e=((b-d)*10000L)/306001L;

 *da=(double)(b-d-z-((e*153L)/5L))+jd;
 if ( e<=13L )	*mo=(int)(e-1L);
 else		*mo=(int)(e-13L);

 *ye=(int)(c-4716L);
 if ( *mo<=2 )	(*ye)++;
}

int get_easter_day(int year)
{
 int	x,a,b,c,d,e,f,g,h,i,k,l,m,p;
 x=year;

 a=x%19;
 b=x/100,c=x%100;
 d=b/4,e=b%4;
 f=(b+8)/25;
 g=(b-f+1)/3;
 h=(19*a+b-d-g+15)%30;
 i=c/4,k=c%4;
 l=(32+2*e+2*i-h-k)%7;
 m=(a+11*h+22*l)/451;
 p=h+l-7*m+22;
 return(p);
}
                                                 
