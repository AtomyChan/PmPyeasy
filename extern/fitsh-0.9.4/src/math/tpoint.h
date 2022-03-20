/*****************************************************************************/
/* tpoint.h								     */
/*****************************************************************************/

/* point.h : operation on 2D points and arrays (domsa@konkoly.hu (2002)) */
/* point.h,v 5.5 2003/05/13 16:09:43 domsa Exp */

#ifndef __TPOINT_H_INCLUDED
#define __TPOINT_H_INCLUDED	1

typedef struct
{	int	id;
	double	xcoord,ycoord;
} tpoint;

typedef struct
{	int	length;
	tpoint	*points;
} tpointarr;

tpointarr	*tpoint_createarr(void);
tpointarr	*tpoint_buildarr(int length,tpoint *points);
void		tpoint_destroyarr(tpointarr *pa);

double		tpoint_eucdist(tpoint *p1,tpoint *p2);
double		tpoint_eucdist2(tpoint *p1,tpoint *p2);
double		tpoint_minx(int nelem,tpoint *points);
double		tpoint_miny(int nelem,tpoint *points);
double		tpoint_maxx(int nelem,tpoint *points);
double		tpoint_maxy(int nelem,tpoint *points);
void		tpoint_minmaxx(int nelem,tpoint *points,double *minx,double *maxx);
void		tpoint_minmaxy(int nelem,tpoint *points,double *miny,double *maxy);
double		tpoint_lengthx(int nelem,tpoint *points);
double		tpoint_lengthy(int nelem,tpoint *points);

void		tpoint_flipx(int nelem,tpoint *points);
void		tpoint_flipy(int nelem,tpoint *points);
void		tpoint_magx(int nelem,tpoint *points,double xmag);
void		tpoint_magy(int nelem,tpoint *points,double ymag);
void		tpoint_zoom(int nelem,tpoint *points,double zoom);
void		tpoint_shiftx(int nelem,tpoint *points,double xshift);
void		tpoint_shifty(int nelem,tpoint *points,double yshift);
void		tpoint_rotate(int nelem,tpoint *points,double rotang);
void		tpoint_rrotate(int nelem,tpoint *points,double rotangdeg);

int		tpoint_sortx(const void *, const void *);
int		tpoint_sorty(const void *, const void *);
int		tpoint_sortid(const void *, const void *);

tpoint	*	tpoint_copy(int nelem,tpoint *points);
void		tpoint_remove(int idx, int nelem,tpoint *points);

/*****************************************************************************/

#endif

/*****************************************************************************/
              
