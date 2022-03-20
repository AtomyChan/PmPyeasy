/*****************************************************************************/
/* point.h								     */
/*****************************************************************************/

#ifndef __POINT_H_INCLUDED
#define	__POINT_H_INCLUDED	1

/* The 'point'  sturcture:
    - x,y: the 2D coordinates of the point (two real numbers)
    - value: an arbitrary quantity of the point
    - weight: the weight of a point (used by the fitting routines)
*/
typedef struct
 {	double	x,y;
	double	value;
	double	weight;
 } point;

#endif
                                                                     
