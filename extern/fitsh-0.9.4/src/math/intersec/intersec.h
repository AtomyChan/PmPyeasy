/*****************************************************************************/
/* intersec.h 								     */
/*****************************************************************************/

#ifndef __INTERSEC_H_INCLUDED
#define	__INTERSEC_H_INCLUDED	1

typedef struct
 {	double	x,y;
	double	sx,sy;
 } drectangle;

typedef struct
 {	double	x,y;
	double	radius;
 } dcircle;

/* intersec_rectangles():
   Calculates the instersection of the two rectangles 'r1' and 'r2'. If
   they are disjoint, the function returns a non-zero value, otherwise
   it returns 0 and the intersection (obviously, it's a recangle) is 
   stored in 'ret'.							     */
int	intersec_rectangles(drectangle *r1,drectangle *r2,drectangle *ret);

/* area_rectangle(), area_circle():
   Calculates the area of the rectangle 'r' and the circle 'c'.	Obviously,
   these functions uses only the structure fields 'sx', 'sy' and 'radius'.   */
double	area_rectangle(drectangle *r);
double	area_circle(dcircle *c);

/* area_intersec_of_rect_circ():
   Calculates the exact area of the intersection of the rectangle 'r' and the
   circle 'c'. If they are disjoint, the function returns 0, otherwise it
   returns a positive value.						     */
double	area_intersec_of_rect_circ(drectangle *r,dcircle *c);

#endif
                     
