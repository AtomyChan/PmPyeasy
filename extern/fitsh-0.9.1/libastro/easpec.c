/*****************************************************************************/
/* easpec.c 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Calculates Earth-specific quantities (inclination of the rotation axis,   */
/* sideral time, geocentric position of an observer, ...)		     */
/* This part of the library is not standalone (depends on earth and cmath).  */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2003, 2006; Pal, A. (apal@szofi.elte.hu).			     */ 
/*****************************************************************************/

#include <stdio.h>
#include <math.h>

#include <astro/astro.h>
#include <astro/earth.h>
#include <astro/easpec.h>

/*****************************************************************************/

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841968
#endif

/*****************************************************************************/

double get_earth_axis_angle(double jdate)
{
 double	jc;
 static double cjdate=0.0,ce=0.0;

 if ( cjdate==jdate )	return(ce);

 jc=(jdate-2451545.0)/36525.0;
 ce=23.439291+(-0.0130042+(-0.00000016+0.000000504*jc)*jc)*jc;
 cjdate=jdate;
 return(ce);
}

double get_sideral_time(double jdate)
{
 double	jc,jday;
 static double cjdate=0.0,csid=0.0;

 if ( cjdate==jdate )	return(csid);

 jday=(double)((long)(jdate-0.5))+0.5;
 jc=(jday-2415020.0)/36525.0;
 csid=0.276919398+(100.0021359+0.000001075*jc)*jc;

 csid+=1.002737908*(jdate-jday);
 if ( csid<0 )	csid+=(double)((long)(-csid))+10.0;
 csid=csid-(double)((long)csid);
 csid=csid*360.0;

 cjdate=jdate;

 return(csid);
}

void get_observer_coords(double jdate,double lon,double lat,double *x,double *y,double *z)
{
 double	wx,wy,n,d,lst;

 lst=get_sideral_time(jdate)+lon;

 wx=r_earth_equ*cosa(lat);
 wy=r_earth_pol*sina(lat);
 n=sqrt(wx*wx+wy*wy);
 d =r_earth_equ*wx/n;
 *z=r_earth_pol*wy/n;
 *x=d*cosa(lst);
 *y=d*sina(lst);
}

double get_heliocentric_julian_date(double jd,double ra,double dec)
{
 double rx,ry,rz;
 double nx,ny,nz;
 double	eps,dis,hjd;

 get_earth_coords(jd,&rx,&ry,&rz);
 eps=get_earth_axis_angle(jd);

 nx=cosa(dec)*cosa(ra);
 ny=cosa(dec)*sina(ra);
 nz=sina(dec);

 rotate(&ny,&nz,-eps);

/* jd0=get_julian_date(2000,1,1); */
/* rotate(&nx,&ny,(jd-jd0)/(365.25*25800.0)*360.0); */

 dis=nx*rx+ny*ry+nz*rz;

 hjd=jd + dis * 0.0057755178;

 return(hjd);
}

double get_heliocentric_julian_date_epoch(double jd,double ra,double dec,double y0)
{
 double rx,ry,rz;
 double nx,ny,nz;
 double	eps,dis,hjd,jd0;

 get_earth_coords(jd,&rx,&ry,&rz);
 eps=get_earth_axis_angle(jd);

 nx=cosa(dec)*cosa(ra);
 ny=cosa(dec)*sina(ra);
 nz=sina(dec);

 rotate(&ny,&nz,-eps);

 jd0=2451545.0+(y0-2000.0)*365.25;
 rotate(&nx,&ny,(jd-jd0)/(365.25*25800.0)*360.0);

 dis=nx*rx+ny*ry+nz*rz;

 hjd=jd + dis * 0.0057755178;

 return(hjd);
}

double get_barycentric_julian_date(double jd,double ra,double dec)
{
 double rx,ry,rz;
 double nx,ny,nz;
 double	bx,by,bz,tmass;
 double	jupx,jupy,jupz,
	satx,saty,satz,
	urax,uray,uraz,
	nepx,nepy,nepz;
 double	eps,dis,hjd;

 get_planet_coords(jd,5,&jupx,&jupy,&jupz);
 get_planet_coords(jd,6,&satx,&saty,&satz);
 get_planet_coords(jd,7,&urax,&uray,&uraz);
 get_planet_coords(jd,8,&nepx,&nepy,&nepz);
 tmass=1+m_jupiter+m_saturn+m_uranus+m_neptune;
 bx=(jupx*m_jupiter+satx*m_saturn+urax*m_uranus+nepx*m_neptune)/tmass;
 by=(jupy*m_jupiter+saty*m_saturn+uray*m_uranus+nepy*m_neptune)/tmass;
 bz=(jupz*m_jupiter+satz*m_saturn+uraz*m_uranus+nepz*m_neptune)/tmass;

 get_earth_coords(jd,&rx,&ry,&rz);

 rx-=bx;
 ry-=by;
 rz-=bz; 

 eps=get_earth_axis_angle(jd);

 nx=cosa(dec)*cosa(ra);
 ny=cosa(dec)*sina(ra);
 nz=sina(dec);

 rotate(&ny,&nz,-eps);

/* jd0=get_julian_date(2000,1,1); */
/* rotate(&nx,&ny,(jd-jd0)/(365.25*25800.0)*360.0); */

 dis=nx*rx+ny*ry+nz*rz;

 hjd=jd + dis * 0.0057755178;

 return(hjd);
}

int get_annular_parallax_correction(double jd,double ra,double dec,double dis,double *rra,double *rdec)
{
 double	x0,y0,z0,dx,dy,dz,eps;
 
 x0=dis*cosa(dec)*cosa(ra);
 y0=dis*cosa(dec)*sina(ra);
 z0=dis*sina(dec);
 
 get_earth_coords(jd,&dx,&dy,&dz);
 eps=get_earth_axis_angle(jd);
 rotate(&dy,&dz,+eps);
 x0+=dx*M_PI/(180.0*3600.0);
 y0+=dy*M_PI/(180.0*3600.0);
 z0+=dz*M_PI/(180.0*3600.0);

 *rra=getpcoords(x0,y0);
 dis=sqrt(x0*x0+y0*y0+z0*z0);
 *rdec=asina(z0/dis);
 
 return(0);
}

/*****************************************************************************/
                                  
