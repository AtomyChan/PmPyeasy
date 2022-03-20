/*****************************************************************************/
/* linkpoint.h								     */
/*****************************************************************************/

#ifndef	__LINKPOINT_H_INCLUDED
#define	__LINKPOINT_H_INCLUDED	1

/*****************************************************************************/

typedef struct
 {	short	nx,ny;		/* indices of the next point in the link:    */
				/* - differs from the array index if normal  */
				/* - the same as the array index if extreme  */
				/* - negative if invalid		     */
	short	mx,my;		/* index of the extreme point in the link    */
 } linkpoint;

typedef struct
 {	int	mcnt;		/* number of extreme link references	     */
	int	identifier;	/* an identifier for the equivalence class   */
	short	ncnt;		/* number of next point references	     */
	short	nngh;		/* number of valid neighbors		     */
	short	nsgn;		/* number of neighbors within the same group */
				/* (latter two include the point itself)     */
	short	aux;		/* auxiliary/temporary variable		     */ 
 } linkreference;

typedef struct
 {	int	xmin,xmax;
	int	ymin,ymax;
 } linkrange;

/*****************************************************************************/

/* linkpoint_local_extreme():
   This function searches the local extreme value of the distribution 'arr'
   around the pixel (x,y). The array has the size of sx times sy. If the
   mask 'mask' is not NULL, the points where mask is non-zero are also
   excluded. The flag is positive if the extreme value is a maximum, otherwise
   the function searches for a minimum. The function returns the total number
   of points/pixels which were analyzed. This value is something between 0 and
   9. If the exclude_mask parameter is non-zero the pixels (respecing to the
   set bits of exclude_mask) around (x,y) are excluded from the search. 
   The coordinates of the pixel with the extreme value are returned in rnx 
   and rny, if at least one point is investigated (i.e. the return value is 
   positive), otherwise these values are unchanged. If an error occurs, the 
   function returns a negative value.
   The mask values are related to the neighbors as follows:
 	0x001: ( x-1 , y-1 )
	0x002: ( x   , y-1 )
	0x004: ( x+1 , y-1 )
	0x008: ( x-1 ,  y  )
	0x010: (  x  ,  y  )
	0x020: ( x+1 ,  y  )
	0x040: ( x-1 , y+1 )
	0x080: (  x  , y+1 )
	0x100: ( x+1 , y+1 )						     */
int	linkpoint_local_extreme(double **arr,int ox,int oy,int sx,int sy,char **mask,
	int x,int y,int flag,int exclude_mask,int *rnx,int *rny);

/* linkpoint_mask_same_group():
   This function creates and returns an exclusion mask based on the linkpoint
   array 'lparr' (with the size of sx times sy) for the given pixel (x,y).
   The points with invalid (negative) values of nx and ny and points which are
   in the same group as the pixel (x,y) are included in the mask. 	     */
int	linkpoint_mask_same_group(linkpoint **lparr,int sx,int sy,int x,int y);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* linkpoint_create():
   Creates a two dimensional (dynamically allocated) array of points with the 
   type of 'linkpoint' holding the information about the convergence to 
   local maximums ('flag' is positive) or local minimums ('flag' is negative). 
   The elements nx and ny in the structure 'linkpoint' are set, these are the
    array indices of the largest/smallest neighbour.
   The optional parameter 'lr' can be used to specify a region inside
   the rectangle [0,sx-1] x [0,sy-1] where this maximum/minimum search 
   should be performed. Thee optional argument 'mask' may contain a _bad_
   mask of points which should not be taken into account (bad means 
   the array element mask[i][j] should be non-zero if this point wants to
   be rejected). Both 'lr' and 'mask' can be NULL. If the linking was 
   successful, the function returns an array of 'sx' times 'sy' of linkpoints
   with the desired information, otherwise (when an error occurs) it returns
   NULL (even when 'flag' is zero: it must be any arbitrary positive or 
   negative integer, depending on the request, but not zero).		     */
linkpoint **	linkpoint_create(double **arr,int sx,int sy,linkrange *lr,
		char **mask,int flag);

/* linkpoint_get_link_end():
   Gets the endpoint of the pixel (nx,ny) based on the link array 'lparr'.
   If there is an endpoint, the function returns 0, otherwise it returns -1. */
int		linkpoint_get_link_end(linkpoint **lparr,
		int nx,int ny,int *rmx,int *rmy);

/* linkpoint_is_same_endpoint():
   This function returns a positive value if the endpoints of the links 
   defined by (x1,y1) and (x2,y2) are the same, otherwise, if the linkpoints
   are valid, the function returns 0; in the other cases (i.e. when any of 
   the points are invalid), the function returns a negative value. The 
   important role of this function is to check this property without the
   usage of the mx and my fields (i.e. if these are set, two points defines
   the same link if and only if their mx and my values are the same).	     */
int		linkpoint_is_same_endpoint(linkpoint **lparr,
		int x1,int y1,int x2,int y2);

/* linkpoint_reconnect():
   This function sets the linkpoint structure elements and mx and my, which 
   are the array indices of the local maximum/minimum can be reached via the 
   links of nx/ny points. The array 'lparr' should be (somehow) preinitialized,
   i.e. the nx and ny values should be proper ones. Generally this 
   initialization is done by linkpoint_create(). Note that this and the 
   previous functions are disjoint in the manner that the former sets only
   the nx and ny elements and sets mx and my to (-1) while the latter does
   not affect the values nx and ny, only sets mx and my. 		     */
int		linkpoint_reconnect(linkpoint **lparr,int sx,int sy);

/* linkpoint_free():
   Releases the array 'arr' which was dynamically allocated and returned by a 
   previous call of linkpoint_create(). This function always return zero.    */
int		linkpoint_free(linkpoint **arr);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* linkreference_create():
   This function creates a dynamically allocated array of reference counts
   of the linkpoint array 'lparr' (which has a size of sx times sy). The 
   linkpoint array should be a reconnected one, i.e. the mx and my values
   should be initialized by linkpoint_reconnect(). All values (ncnt, mcnt,
   nngh and nsgn) are set for a valid (non-masked) point. For an invalid one
   (where nx or ny is negative), all values are set to a negative value.     */
linkreference **linkreference_create(linkpoint **lparr,int sx,int sy);

/* linkreference_free():
   Releases the array 'arr' which was dynamically allocated and returned by a 
   previous call of linkreference_create(). This function always return zero.*/
int		linkreference_free(linkreference **arr);

/*****************************************************************************/

#endif

/*****************************************************************************/
          
