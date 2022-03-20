/*****************************************************************************/
/* delaunay.h								     */
/*****************************************************************************/

/* triangulation.c : triangulate data points (domsa@konkoly.hu (2002)) */
/* written for the HAT Project */
/* triangulation.h,v 5.5 2003/05/13 16:09:44 domsa Exp */

#ifndef	__DELAUNAY_H_INCLUDED
#define __DELAUNAY_H_INCLUDED	1

#include "tpoint.h"

typedef struct 
 {	tpoint		*vertices[3];
	double		trix,triy;	/* triangulate space coordinates     */
 } triangle;				/* (originally, ratios b/a and c/a)  */

typedef struct
 {	tpoint		***neigs;
	tpoint		**neigpoints;
 } triinfo;

/* delaunay_triangulation():
   Creates the Delaunay triangule list (using divide and conquer method) of
   the point list 'ps' (with 'np' points, if it is not less than 3). 
   The triangulation is returned in 'rtris' and the number of triangles
   is stored in 'rntri'. If the triangulation is successful, the function
   returns zero, otherwise a non-zero value. If 'ti' is not NULL, the
   function stores some auxiliary information in the fields of 'ti'
   (e.g. list of neighbouring points, see `typedef ... triinfo` above.	     */
int	delaunay_triangulation(tpoint *ps,int np,
		triangle **rtris,int *rntri,triinfo *ti);

/* free_triangulation_info():
   Releases the dynamically allocated information in the triinfo structure
   'ti' which should be previously set by delaunay_triangulation().	     */
int	free_triangulation_info(triinfo *ti);

/* expand_triangulation():
   Expands the Delaunay triangulation stored in 'ti' to the level 'level'. 
   The expanded set of triangles are going to be stored in 'rtris' 
   (with 'rntri' triangles) and the list 'points' and 'npoint' should be 
   the same what was passed to delaunay_triangulation().		     */
int	expand_triangulation(triinfo *ti,triangle **rtris,int *rntri,
		tpoint *points,int npoint,int level);

/* vertice_sort():
   Auxiliary function, used only by the trimatch library.		     */
void	vertice_sort(double *sides,tpoint **vertices);

#endif
           
