/*****************************************************************************/
/* trimatch.h								     */
/*****************************************************************************/

/* trimatch.h : triangle match of 2D pointsets (domsa@konkoly.hu (2002)) */
/* written for the HAT Project */
/* trimatch.h,v 5.5 2003/05/13 16:09:44 domsa Exp */

#ifndef __TRIMATCH_H_INCLUDED
#define __TRIMATCH_H_INCLUDED	1

#include "tpoint.h"		/* required by 'tpoint'		             */
#include "cpmatch.h"		/* required by 'cphit'		  	     */
#include "poly.h"		/* useful for handling the results '[xy]fit' */

typedef struct
 {	int	level;		/* triangulation level: 		     */
				/*  - <0: full triangulation, 		     */
				/*  - 0 : Delaunay-triangulation,	     */
				/*  - 0<: expanded Delaunay triangulation    */

	double	unitarity;	/* automatic triangulation expansion until   */
				/* the unitarity of the linear part of the   */
				/* transformation reaches this value, setting*/
				/* it to positive overrides the triangulation*/
				/* level specified by 'level'	  	     */

	double	maxdist;	/* maximal distance of matched point-pairs   */

	int	parity;		/* parameter for trinagle-space coordinate   */
				/* derivation. Possible values can be:	     */
				/*  - positive: force the same chirality     */
				/*  - negative: force different chirality    */
				/*  - zero: mixed coordinates (left-handed   */
				/*    and right-handed triangles appear      */
				/*    to be the same in the triangle-space)  */

 } trimatchpar;

typedef struct
 { 	double	unitarity;	/* final unitarity of the transformation     */
	int	level_init,	/* initial level (see trimatchpar: level)    */
		level_used;	/* final level when points have been matched */
	double	residual;	/* residual of the initial transformation fit*/
	double	time_used[8];	/* net cpu usage (in seconds): [0]: delaunay */
				/* [<0]: expansion of the triangulation	     */
 } trimatchlog;

/* trimatch():
   Matches the two point sets p1[] and p2[] (with n1 and n2 points, 
   respectively) using triangulate match algorithm. The resulting polynomial 
   fit (up to order specified by 'order') is stored in xfit[] and yfit[]
   while the residual is stored in 'sigma'. The resulting set of matched pairs
   is stored in 'hits' (which is dynamically allocated, can be released
   by free). If 'hits' is NULL, the set of matched pairs is not stored. 
   The number of matched pairs is stored in 'nhit'. Some fine-tune parameters
   can be passed in the optional parameter 'tmp', if it is not NULL. 
   If 'tmp' is NULL, the default parameters (raw Delaunay-triangulation with
   no maximal distance constraint) are used. The function returns 0 if
   the matching was successful, otherwise it returns a non-zero value.
   Some log about the triangulation and the matching are written into
   'tml' if it is desired (so, 'tml' is not NULL).		             */
int trimatch(tpoint *p1,int n1,tpoint *p2,int n2,int order,trimatchpar *tmp,
	cphit **hits,int *nhit,double *xfit,double *yfit,trimatchlog *tml);

#endif
                           
