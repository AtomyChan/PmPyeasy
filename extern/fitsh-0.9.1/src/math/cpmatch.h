/*****************************************************************************/
/* cpmatch.h								     */
/*****************************************************************************/

/* cpmatch.h : match two 2D point array (domsa@konkoly.hu (2003)) */
/* cpmatch.h,v 5.5 2003/05/13 16:09:43 domsa Exp */

#ifndef	__CPMATCH_H_INCLUDED
#define	__CPMATCH_H_INCLUDED	1

#include "tpoint.h"

typedef struct
 {	int	idx[2];		/* array indexes of the pairs		*/
	double	distance;	/* squared distance of the pairs	*/
 } cphit;

typedef enum 
 {	CPMTYPE_GUESS,
	CPMTYPE_ACC,
	CPMTYPE_DIST
 } cpmtype;

/*****************************************************************************/

/* functions for internal usage, but can be useful for other applications... */

/* cpmatch_find_xidx():
   Finds the point in the point array 'base' which is close to 'probe' along
   the x-axis. The function returns a negative number on error, otherwise it
   returns the array index of the matched point (which is going to be between 
   0 and base->length-1). Note that this point (the returned one) is _not_ 
   the nearest point, it is just a hint for further usage, e.g. 
   cpmatch_find_nearest(). Note also that the array 'base' should be sorted 
   (along the x-axis, in ascending order) before call. Do something like
	qsort((void *)base->points,base->length,sizeof(tpoint),tpoint_sortx),
   using the library function tpoint_sortx() of the tpoint.c module.	     */
int	cpmatch_find_xidx(tpoint *probe,tpointarr *base);
/* cpmatch_find_nearest():
   Finds the point in the point array 'base' which is nearest to 'probe'. 
   The function returns a negative number on error, otherwise it returns the 
   array index of the matched point (which is going to be between 0 and 
   base->length-1). The 'neighbour' is an index of a point which is probably
   close to 'probe'. Generally, it can be derived using cpmatch_find_xidx().
   The squared distance of the 'probe' point and the matched one is 
   returned in 'retdistance' if it was not NULL. Note that the array 'base'
   should be sorted (see cpmatch_find_xidx() above).			     */
int	cpmatch_find_nearest(tpoint *probe,tpointarr *base,int neighbour,
	double *retdistance,int is_exclude_self);

/* cpmatch_find_neighbour():
   Finds the neighbour point of the 'indx'th point in the point array 'base'.
   If the 'indx' has a bad value of 'base' contains less than 2 points,
   the function returns a negative value, otherwise it returns the desired
   index. Like the other cpmatch_find_*() functions, the array 'base' 
   should be sorted (see cpmatch_find_xidx() above).			     */ 
int	cpmatch_find_neighbour(tpointarr *base,int indx);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* cpmatch_symmetric():
   Matches two 2 dimensional point sets 'cat1' and 'cat2'. A P1 point from
   'cat1' and P2 from 'cat2' are treated as matched pairs if the nearest
   point from 'cat2' to P1 is P2 and the nearest point from 'cat1' to P2 
   is P1. The function returns an array of 'cphit' object (with 'rlen'
   elements) which is a dynamically allocated array, can be released by free().
   The parameter 'rev_sig' can be used for some post-matching rejection. 
   If it is greater than 0, the median of the distances between the matched
   pairs are calculated and pairs with points farer than 'rev_sig' times
   the median are rejected and removed from the returned array.		     */
cphit	*cpmatch_symmetric(tpointarr *cat1,tpointarr *cat2,
	int *rlen,double rev_sig,double rev_maxdist);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* cpmatch():
   The original matching function written by Istvan. It is not used, obsoleted
   by the functions (esp. cmatch_symmetric()) above. Note that cpmatch()
   returns a preliminary array of pairs which should be revised (God knows
   how) by cphit_revise(). Ask Istvan for details... The other functions
   (cphit_sort_dist and cphit_sort_id) may also have to be used somewhere.   */
cphit	*cpmatch(tpointarr *reference,tpointarr *input,int *hlength);
int	cphit_sort_dist(const void *cphit1, const void *cphit2);
int	cphit_sort_id(const void *cphit1, const void *cphit2);
int	cphit_revise(int hlen,cphit **hits, int maxarrlen,
	cpmtype matchtype,double *maxdist);

/*****************************************************************************/

#endif

/*****************************************************************************/
             
