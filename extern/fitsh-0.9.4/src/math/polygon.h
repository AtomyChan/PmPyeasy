/*****************************************************************************/
/* polygon.h								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2015; Pal, A. <apal@szofi.net>					     */
/*****************************************************************************/

#ifndef	__POLYGON_H_INCLUDED
#define	__POLYGON_H_INCLUDED	1

/*****************************************************************************/

double	polygon_area(double *poly,int n);
double	polygon_area_signed(double *poly,int n);
int	polygon_integrate_linear_monoms(double *poly,int n,double *fcoeff);
int	polygon_convexity(double *poly,int n);
int	polygon_intersection_halfplane(double *poly,int n,double px,double py,double nx,double ny);
int	polygon_intersection_square(double *poly,int n,double x0,double y0,double dx,double dy);

/*****************************************************************************/

#endif

/*****************************************************************************/
                                               
                           
