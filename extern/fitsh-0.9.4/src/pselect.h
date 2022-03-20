/*****************************************************************************/
/* pselect.h								     */
/*****************************************************************************/

#ifndef	__PSELECT_H_INCLUDED
#define	__PSELECT_H_INCLUDED	1

#include "math/point.h"

/* point_select():
   Selects points from the set 'points' (which contains 'npoint' points). The 
   function tries to select 'nselect' points from the set as homogeneously, as 
   it is declared by 'level'. It is not guaranteed that, even if 
   'nselect' <= 'npoint', the desired number of points will be selected
   (it might happen when the initial set is very inhomogeneous). Only
   the point in the rectangle [x0:x1,y0:y1] will be selected. The array 'ret'
   is filled by this function with false (zero) or true (nonzero) values:
   if 'ret[i]' is true, the point 'points[i]' should be treated as a selected
   one. The number of selected points are returned (which is always less or
   equal to 'nselect'). The parameter 'level' indicates the desired level,
   whether homogeneousity or the weightness (indicaded for each point by
   the field 'weight' in the structure 'point') is more important during the
   selection. If 'level' is 0, only the 'nselect' points with the highest
   weights will be selected. The more higher 'level' is, the more homogeneous
   set will be returned. If 'level' is negative, the function does not take
   the weights into account.						     */
int point_select(point *points,int npoint,char *ret,int nselect,
		 double x0,double y0,double x1,double y1,int level);

#endif
         
